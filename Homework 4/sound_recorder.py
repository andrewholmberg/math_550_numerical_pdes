import sounddevice as sd
import numpy as np
import scipy.io.wavfile as wav

# Set recording parameters
fs = 44100  # Sample rate
seconds = 20  # Duration of recording

# Record audio
myrecording = sd.rec(int(seconds * fs), samplerate=fs, channels=2)
sd.wait()  # Wait until recording is finished

# Save as WAV file 
wav.write("output.wav", fs, myrecording) 