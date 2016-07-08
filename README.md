# pixelsound

Cnvert an image to sound and visualize the produced sound with a spectrogram, where the original image is seen.

**Encoding**:
The conversion is done by encoding each pixel in the image to a sine wave, where frequency corresponds to the row of the pixel and amplitude corresponds to the intensity of the pixel.

**Decoding**:
The audio is decoded by Fourier transform, which decomposes the signal into the frequencies that make it up. The original image is restored with a spectrogram, which shows the spectrum of frequencies as a function of time. Amplitudes of frequencies are represented by the intensity or color in the spectrogram.

## Try it out

Visualize [this audio clip](data/im.wav) with this online spectrum analyzer: https://academo.org/demos/spectrum-analyzer/.
