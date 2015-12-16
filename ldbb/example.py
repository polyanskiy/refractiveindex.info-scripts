#!/usr/bin/python
# example of how to use ldbb.py

from ldbb import LDBB
from numpy import *
import matplotlib.pyplot as plt

# wavelength range of interest
lambda0 = linspace(200e-9,2000e-9,100)
epsD  = LDBB('Au','D',lambda0)
epsLD = LDBB('Au','LD',lambda0)
epsBB = LDBB('Au','BB',lambda0)

# scale wavelength for plot
lambda0 = lambda0*1e6

# LDBB returns permittivity, but we can easily convert it to refractive index

plt.subplot(2,1,1)
plt.title("Refractive Index of Au")
plt.ylabel('n')
plt.xlabel('wavelength (microns)')
plt.plot(lambda0,real(sqrt(epsD)), label="D")
plt.plot(lambda0,real(sqrt(epsBB)), label="LD")
plt.plot(lambda0,real(sqrt(epsLD)), label="BB")
plt.legend()

plt.subplot(2,1,2)
plt.ylabel('k')
plt.xlabel('wavelength (microns)')
plt.plot(lambda0,imag(sqrt(epsD)),label="D")
plt.plot(lambda0,imag(sqrt(epsBB)),label="LD")
plt.plot(lambda0,imag(sqrt(epsLD)),label="BB")
plt.legend()

plt.show()
