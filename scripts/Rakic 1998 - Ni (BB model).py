# -*- coding: utf-8 -*-
# Original data: https://doi.org/10.1364/AO.37.005271
# Author: Mikhail Polyanskiy
# Last modified: 2017-01-16

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz as w
π = np.pi

# Brendel-Bormann (BB) model parameters
ωp = 15.92 #eV
f0 = 0.083
Γ0 = 0.022 #eV

f1 = 0.357
Γ1 = 2.820 #eV
ω1 = 0.317 #eV
σ1 = 0.606 #eV

f2 = 0.039
Γ2 = 0.120 #eV
ω2 = 1.059 #eV
σ2 = 1.454 #eV

f3 = 0.127
Γ3 = 1.822 #eV
ω3 = 4.583 #eV
σ3 = 0.379 #eV

f4 = 0.654
Γ4 = 6.637 #eV
ω4 = 8.825 #eV
σ4 = 0.510 #eV

Ωp = f0**.5 * ωp  #eV

def BB(ω):  #ω: eV
    ε = 1-Ωp**2/(ω*(ω+1j*Γ0))

    α = (ω**2+1j*ω*Γ1)**.5
    za = (α-ω1)/(2**.5*σ1)
    zb = (α+ω1)/(2**.5*σ1)
    ε += 1.j*π**.5*f1*ωp**2 / (2**1.5*α*σ1) * (w(za)+w(zb)) #χ1
    
    α = (ω**2+1j*ω*Γ2)**.5
    za = (α-ω2)/(2**.5*σ2)
    zb = (α+ω2)/(2**.5*σ2)
    ε += 1j*π**.5*f2*ωp**2 / (2**1.5*α*σ2) * (w(za)+w(zb)) #χ2
    
    α = (ω**2+1j*ω*Γ3)**.5
    za = (α-ω3)/(2**.5*σ3)
    zb = (α+ω3)/(2**.5*σ3)
    ε += 1j*π**.5*f3*ωp**2 / (2**1.5*α*σ3) * (w(za)+w(zb)) #χ3
    
    α = (ω**2+1j*ω*Γ4)**.5
    za = (α-ω4)/(2**.5*σ4)
    zb = (α+ω4)/(2**.5*σ4)
    ε += 1j*π**.5*f4*ωp**2 / (2**1.5*α*σ4) * (w(za)+w(zb)) #χ4
    
    return ε
  
ev_min=0.1
ev_max=6
npoints=200
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε = BB(eV)
n = (ε**.5).real
k = (ε**.5).imag


#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(npoints-1, -1, -1):
    file.write('\n        {:.4e} {:.4e} {:.4e}'.format(μm[i],n[i],k[i]))
file.close()
    
    
#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

plt.figure(1)
plt.plot(eV, -ε.real, label="ε1")
plt.plot(eV, ε.imag, label="ε2")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.2,5])
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs eV
plt.figure(2)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs μm
plt.figure(3)
plt.plot(μm, n, label="n")
plt.plot(μm, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n, k')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)