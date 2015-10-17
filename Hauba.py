"""
Hauba kostur V1
"""

##############
import os, random
import time
import pyaudio
import wave
import sys
import RPi.GPIO as GPIO  
GPIO.setmode(GPIO.BCM)  
##############

LED_R = 9
LED_G = 10
LED_B = 11
BTN_1 = 23
BTN_2 = 24
BTN_3 = 25

# GPIO 23 24 25 set up as input. It is pulled up to stop false signals  
GPIO.setup(BTN_1, GPIO.IN, pull_up_down=GPIO.PUD_UP)  
GPIO.setup(BTN_2, GPIO.IN, pull_up_down=GPIO.PUD_UP)
GPIO.setup(BTN_3, GPIO.IN, pull_up_down=GPIO.PUD_UP)
# GPIO 9 10 11 set up as output
GPIO.setup(LED_R, GPIO.OUT)
GPIO.setup(LED_G, GPIO.OUT)
GPIO.setup(LED_B, GPIO.OUT)

# setting output value RGB
GPIO.output(LED_R, GPIO.LOW)
GPIO.output(LED_G, GPIO.LOW)
GPIO.output(LED_B, GPIO.LOW)

# recording settings
CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 1
RATE = 44100
RECORD_SECONDS = 5
WAVE_OUTPUT_FILENAME = "secrets/secret_"

# recording file while button is pressed file will be placed in secrets folder
def record_secret():
    p = pyaudio.PyAudio()
    stream = p.open(format=FORMAT,
                    channels=CHANNELS,
                    rate=RATE,
                    input=True,
                    frames_per_buffer=CHUNK)
    print("* recording")
    frames = []
    btn2_value = GPIO.input(BTN_2)
    while(not btn2_value):
        data = stream.read(RATE/CHUNK)
        if data != '':
            frames.append(data)
        btn2_value = GPIO.input(BTN_2)
    print("* done recording")
    stream.stop_stream()
    stream.close()
    p.terminate()
    wf = wave.open(WAVE_OUTPUT_FILENAME + time.strftime("%d%m%y%H%M%S", time.localtime()) + ".wav", 'wb')
    wf.setnchannels(CHANNELS)
    wf.setsampwidth(p.get_sample_size(FORMAT))
    wf.setframerate(RATE)
    wf.writeframes(b''.join(frames))
    wf.close()
    time.sleep(1)

# Play file, once or in loop
def playSound(file,loop):
    print ("Playing file: " + file)
    print loop
    # open the file for reading.
    wf = wave.open(file, 'rb')
    # create an audio object
    p = pyaudio.PyAudio()
    # open stream based on the wave object which has been input.
    stream = p.open(format =
                    p.get_format_from_width(wf.getsampwidth()),
                    channels = wf.getnchannels(),
                    rate = wf.getframerate(),
                    output = True)
    # read data (based on the chunk size)
    data = wf.readframes(CHUNK)
    # PLAYBACK LOOP
    btn1_value = GPIO.input(BTN_1)
    btn3_value = GPIO.input(BTN_3)
    donePlaying = False
    while  (not donePlaying):
        data = wf.readframes(CHUNK)
        while (data != ''):
            if (loop):
                btn1_value = GPIO.input(BTN_1)
                if (btn1_value):
                    data = ''
            else:            
                btn3_value = GPIO.input(BTN_3)
                if (btn3_value):
                    data = ''
            stream.write(data)
            data = wf.readframes(CHUNK)
        if (btn1_value):
            donePlaying = True
        if (loop): # If file is over and loop is on then rewind.
            wf.rewind()
            print ("Loop on")
        else:
            donePlaying = True
    # cleanup stuff.
    stream.close()    
    p.terminate()
    print "Sound finished"
    time.sleep(1)

while True:
    try:
        btn1_value = GPIO.input(BTN_1)
        btn2_value = GPIO.input(BTN_2)
        btn3_value = GPIO.input(BTN_3)
        if (not btn1_value):
            playSound("music/music.wav",True)
        if (not btn2_value):
            record_secret()
        if (not btn3_value):
            playSound("secrets/" + random.choice(os.listdir("secrets/")),False);
    except KeyboardInterrupt:
        print "Thank you!!!"
        GPIO.cleanup()
        raise
GPIO.cleanup()
