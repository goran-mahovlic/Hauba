import os
from subprocess import call

for file in os.listdir("/home/pi/Hauba/tmp/"):
    if file.endswith(".wav"):
        print(file)
        call(["/home/pi/Hauba/fx/fx", "/home/pi/Hauba/tmp/" + file , "/home/pi/Hauba/secrets/" + file, "80", "1"])
        call(["rm","/home/pi/Hauba/tmp/" + file])
