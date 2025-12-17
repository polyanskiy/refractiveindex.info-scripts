# Author: Mikhail Polyanskiy
# Last modified: 2025-11-11
# Original data: Chen et al. 2009, https://doi.org/10.1364/JOSAB.26.000A58

import numpy as np
import matplotlib.pyplot as plt

# model parameters
A  = 122.3e-27  # cm^6
B  = -22.88e-18 # cm^4
C  = 3.879e-9   # cm^2
S3 = 5.76
Γ3 = 2.8
ωL = 245.5
ωT = 237


def Epsilon(ω):
    ε = A*ω**6 + B*ω**4 + C*ω**2 + S3 \
        + (ωL**2-ωT**2)*S3 / (ωT**2 - ω**2 - 1j*Γ3*ω)
    return ε


ω_min = 6.67128  # cm^-1 (0.2 THz)
ω_max = 3335.64  # cm^-1 (100 THz)
npoints = 500
ω = np.logspace(np.log10(ω_min), np.log10(ω_max), npoints)
μm = 10000/ω
ε = Epsilon(ω)
n = (ε**.5).real
k = (ε**.5).imag


#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(npoints-1, -1, -1):
    file.write('\n        {:.4e} {:.4e} {:.4e}'.format(μm[i],n[i],k[i]))
file.close()
    

#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

#plot ε vs eV
plt.figure(1)
plt.plot(ω, ε.real, label="ε1")
plt.plot(ω, ε.imag, label="ε2")
plt.xlabel('Wave number (1/cm)')
plt.ylabel('ε')
plt.xscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs ω
plt.figure(2)
plt.plot(ω, n, label="n")
plt.plot(ω, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.xscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs μm
plt.figure(3)
plt.plot(μm, n, label="n")
plt.plot(μm, k, label="k")
plt.xlabel('Wave number (1/cm)')
plt.ylabel('n, k')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)