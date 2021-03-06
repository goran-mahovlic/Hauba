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

# led: pin

PIN_LED_RED = 9
PIN_LED_GREEN = 10
PIN_LED_BLUE = 11

# button: pin, inverted

PIN_BUTTON_HEAD = [23, True]
PIN_BUTTON_RECORD = [24, True]
PIN_BUTTON_PLAY = [25, True]

# GPIO 23 24 25 set up as input. It is pulled up to stop false signals  
GPIO.setup(PIN_BUTTON_HEAD[0], GPIO.IN, pull_up_down=GPIO.PUD_UP)  
GPIO.setup(PIN_BUTTON_RECORD[0], GPIO.IN, pull_up_down=GPIO.PUD_UP)
GPIO.setup(PIN_BUTTON_PLAY[0], GPIO.IN, pull_up_down=GPIO.PUD_UP)
# GPIO 9 10 11 set up as output
GPIO.setup(PIN_LED_RED, GPIO.OUT)
GPIO.setup(PIN_LED_GREEN, GPIO.OUT)
GPIO.setup(PIN_LED_BLUE, GPIO.OUT)

# setting output value RGB
GPIO.output(PIN_LED_RED, GPIO.LOW)
GPIO.output(PIN_LED_GREEN, GPIO.LOW)
GPIO.output(PIN_LED_BLUE, GPIO.LOW)

# recording settings
CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 1
RATE = 44100
RECORD_SECONDS = 5
WAVE_OUTPUT_FILENAME = "secret_"

# some buttons are inverted, this returns True when buttons are pressed
def read_button(button):
    if button[1]:
        return not GPIO.input(button[0])
    return GPIO.input(button[0])

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
    while read_button(PIN_BUTTON_RECORD):
        data = stream.read(RATE/CHUNK)
        if data != '':
            frames.append(data)
    print("* done recording")
    stream.stop_stream()
    stream.close()
    p.terminate()
    file_name = "/home/pi/Hauba/tmp/" +  WAVE_OUTPUT_FILENAME + time.strftime("%d%m%y%H%M%S", time.localtime()) + ".wav"
    wf = wave.open(file_name, 'wb')
    wf.setnchannels(CHANNELS)
    wf.setsampwidth(p.get_sample_size(FORMAT))
    wf.setframerate(RATE)
    wf.writeframes(b''.join(frames))
    wf.close()
    time.sleep(1)
    if os.path.getsize(file_name) < 180 * 1024:
        os.remove(file_name)
        print("Removing: ", file_name)

# Play guide
def play_guide(file):
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
    # play stream (looping from beginning of file to the end)
    while data != '':
    # writing to the stream is what *actually* plays the sound.
        if read_button(PIN_BUTTON_RECORD) or read_button(PIN_BUTTON_PLAY):
            break
        stream.write(data)
        data = wf.readframes(CHUNK)

    # cleanup stuff.
    stream.close()    
    p.terminate()

# Play file, once or in loop
def play_sound(file,loop):
    print ("Playing file: ", file)
    print ("Loop", loop)
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
    # PLAYBACK LOOP
    while (True):
        # read chunk
        data = wf.readframes(CHUNK)
        if loop:
            # head button pressed while looping: play guide
            if read_button(PIN_BUTTON_HEAD):
                play_guide("guide/guide.wav")
                break
        else:
            # head button depressed and not looping: stop
            if not read_button(PIN_BUTTON_PLAY):
                break
        # have data: play
        if (data != ''):
            stream.write(data)
        # end of loop: rewind
        elif (loop):
            wf.rewind()
        # end of song: stop
        else:
            break
    # cleanup stuff.
    stream.close()    
    p.terminate()
    print "Sound finished"
    time.sleep(1)

while True:
    try:
        if (not read_button(PIN_BUTTON_HEAD)):
            play_sound("music/EssereDonna.wav", True)
        if (read_button(PIN_BUTTON_RECORD)):
            record_secret()
        if (read_button(PIN_BUTTON_PLAY)):
            play_sound("secrets/" + random.choice(os.listdir("secrets/")),False);
        #print("HEAD", read_button(PIN_BUTTON_HEAD)),
        #print("REC", read_button(PIN_BUTTON_RECORD)),
        #print("PLAY", read_button(PIN_BUTTON_PLAY))
    except KeyboardInterrupt:
        print("Thank you!!!")
        GPIO.cleanup()
        raise

GPIO.cleanup()
