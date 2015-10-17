#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <string.h>
#include <fftw3.h>
#include <math.h>

#define BUFFER_LEN 1024

#define M_2PI 6.2831853071795
#ifndef M_PI
#define M_PI M_2PI/2.0
#endif
#define MAX_FRAME_LENGTH 4096

#define ROUND(a) ((float)((int)(a+0.5)))

enum{
  false = 0,
  true = 1
};

const int fftFrameSize = 2048;
const int overlap = 4;


void phaseVocAnalysis(fftw_complex *block, float *gLastPhase, double freqPerBin, double expct, float *gAnaMagn, float *gAnaFreq){

    double real, imag, magn, phase, tmp;
    long qpd, k;

    for (k = 0; k <= fftFrameSize/2; k++) {

        /* de-interlace FFT buffer */
        real = block[k][0];
        imag = block[k][1];

        /* compute magnitude and phase */
        magn = 2.*sqrt(real*real + imag*imag);
        phase = atan2(imag,real);

        /* compute phase difference */
        tmp = phase - gLastPhase[k];
        gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*expct;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        /* get deviation from bin frequency from the +/- Pi interval */
        tmp = overlap*tmp/(M_2PI);

        /* compute the k-th partials' true frequency */
        tmp = (double)k*freqPerBin + tmp*freqPerBin;

        /* store magnitude and true frequency in analysis arrays */
        gAnaMagn[k] = magn;
        gAnaFreq[k] = tmp;
    }
}

void phaseVocSynthesis(fftw_complex *block, float *gSumPhase, float *gSynMagn, float *gSynFreq, double freqPerBin, double expct){

    int k;
    double magn, tmp, phase;

    for (k = 0; k <= fftFrameSize/2; k++) {
        /* get magnitude and true frequency from synthesis arrays */
        magn = gSynMagn[k];
        tmp = gSynFreq[k];

        /* subtract bin mid frequency */
        tmp -= (double)k*freqPerBin;

        /* get bin deviation from freq deviation */
        tmp /= freqPerBin;

        /* take overlap into acnframes */
        tmp = M_2PI*tmp/overlap;

        /* add the overlap phase advance back in */
        tmp += (double)k*expct;

        /* accumulate delta phase to get bin phase */
        gSumPhase[k] += tmp;
        phase = gSumPhase[k];

        /* get real and imag part and re-interleave */
        block[k][0] = magn*cos(phase);
        block[k][1] = magn*sin(phase);
    }
}


int main(int argc, char **argv){ 
  static float data [BUFFER_LEN];
  SNDFILE      *infile, *outfile;
  SF_INFO      sfinfo;
  int          nframes;


  float   sPitchFactor = 1, sOutputGain = 1.0;
  float   cEffect = 0;

  float   *gInFIFO, *gOutFIFO, *gOutputAccum;
  float   *window;

  double  *fftTmpR;
  fftw_complex *fftTmpC;
  fftw_complex *fftOldC;
  fftw_plan fftPlan;

  char in_file[80];
  char out_file[80];

  if(argc!=5){
      printf("usage: %s <in_file> <out_file> <pitch_percent> <robot>\n", argv[0]);
      return -1;
  }
  strcpy(in_file, argv[1]);
  strcpy(out_file, argv[2]);
  sPitchFactor = (float)(atoi(argv[3]))/100.0;
  printf("pitch factor: %f\n", sPitchFactor);
  if(argv[4][0]=='1'){
      cEffect=1.0;
  }


  window=(float*)malloc(fftFrameSize*sizeof(float));
  for(int kk=0;kk<fftFrameSize;kk++){
      window[kk] = -.5*cos(M_2PI*(float)kk/(float)fftFrameSize)+.5;
  }

  gInFIFO=(float*)calloc(fftFrameSize, sizeof(float));
  gOutFIFO=(float*)calloc(fftFrameSize, sizeof(float));
  gOutputAccum=(float*)calloc(2*fftFrameSize, sizeof(float));

  // FFTW stuff
  fftTmpR=(double*)fftw_malloc(sizeof(double)*fftFrameSize);
  fftTmpC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);
  fftOldC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);

  static float gLastPhase[MAX_FRAME_LENGTH/2+1];
  static float gSumPhase[MAX_FRAME_LENGTH/2+1];
  static float gAnaFreq[MAX_FRAME_LENGTH];
  static float gAnaMagn[MAX_FRAME_LENGTH];
  static float gSynFreq[MAX_FRAME_LENGTH];
  static float gSynMagn[MAX_FRAME_LENGTH];

  static long gRover = false, gInit = false;

  double freqPerBin, expct;
  long i,k, index, inFifoLatency, stepSize, fftFrameSize2;

  float *fPointer, *fPointer2;
  double *dPointer;


  if(!(infile = sf_open(in_file, SFM_READ, &sfinfo))){
    printf ("Not able to open input file %s.\n", in_file);
    sf_perror (NULL);
    return  1;
  } 

  if(sfinfo.channels > 1){
    printf ("Not able to process more than 1 channel\n");
    return  1;
  }

  if(!(outfile = sf_open(out_file, SFM_WRITE, &sfinfo))){
    printf ("Not able to open output file %s.\n", out_file);
    sf_perror (NULL);
    return  1;
  }



  /* set up some handy variables */
  fftFrameSize2 = fftFrameSize/2;

  stepSize = fftFrameSize/overlap;
  freqPerBin = (double)sfinfo.samplerate/(double)fftFrameSize;
  expct = M_2PI*(double)stepSize/(double)fftFrameSize;

  inFifoLatency = fftFrameSize-stepSize;
  if (gRover == false) gRover = inFifoLatency;

  /* initialize our static arrays */
  if (gInit == false) {
    memset(gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
    memset(gSumPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
    memset(gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(float));
    memset(gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(float));
    gInit = true;
  }



  while((nframes = sf_read_float (infile, data, BUFFER_LEN))){
    


    /* main processing loop */
    for (i = 0; i < nframes; i++){

      // As long as we have not yet collected enough data just read in
      gInFIFO[gRover] = data[i];

      data[i] = gOutFIFO[gRover-inFifoLatency];
      gRover++;

      // now we have enough data for processing
      if (gRover >= fftFrameSize) {
        gRover = inFifoLatency;

        dPointer=fftTmpR;
        fPointer=gInFIFO;
        fPointer2=window;
        float tmp_power=0;
        for (k = 0; k < fftFrameSize;k++) {
          *dPointer=*(fPointer++) * *(fPointer2++);
          dPointer++;
          tmp_power += *dPointer * *dPointer;
        }
        tmp_power/=(float)fftFrameSize;

        // do transform
        fftPlan=fftw_plan_dft_r2c_1d(fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
        fftw_execute(fftPlan);
        fftw_destroy_plan(fftPlan);

        memcpy(fftOldC, fftTmpC, fftFrameSize*sizeof(fftw_complex));

        // pitch shifting with phase vocoder
        phaseVocAnalysis(fftTmpC, gLastPhase, freqPerBin, expct, gAnaMagn, gAnaFreq);
        memset(gSynMagn, 0, fftFrameSize*sizeof(float));
        memset(gSynFreq, 0, fftFrameSize*sizeof(float));
        for (k = 0; k <= fftFrameSize2; k++) {
          index = k*sPitchFactor;
          if (index <= fftFrameSize2) {
            gSynMagn[index] += gAnaMagn[k];
            gSynFreq[index] = gAnaFreq[k] * sPitchFactor;
            if(cEffect)
              gSynFreq[index] = 0;
          }
        }
        phaseVocSynthesis(fftTmpC, gSumPhase, gSynMagn, gSynFreq, freqPerBin, expct);


        // do inverse transform
        fftPlan=fftw_plan_dft_c2r_1d(fftFrameSize, fftTmpC, fftTmpR, FFTW_ESTIMATE);
        fftw_execute(fftPlan);
        fftw_destroy_plan(fftPlan);

        fPointer=gOutputAccum; dPointer=fftTmpR; fPointer2=window;
        for(k=0; k < fftFrameSize; k++) {
          *fPointer += 0.7 * *(dPointer++) / (fftFrameSize2*overlap) * sOutputGain * *(fPointer2++);
          fPointer++;
        }

        memcpy(gOutFIFO, gOutputAccum, stepSize*sizeof(float));

        // shift accumulator
        memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(float));

        // move input FIFO
        for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k+stepSize];

      }
    }
    sf_write_float(outfile, data, nframes);
  }

  sf_close (infile);
  sf_close (outfile);
  return 0;
}


