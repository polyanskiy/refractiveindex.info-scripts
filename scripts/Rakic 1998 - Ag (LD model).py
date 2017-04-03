# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-02
# Original data: Rakić et al. 1998, https://doi.org/10.1364/AO.37.005271

import numpy as np
import matplotlib.pyplot as plt

# Lorentz-Drude (LD) model parameters
ωp = 9.01  #eV
f0 = 0.845
Γ0 = 0.048 #eV

f1 = 0.065
Γ1 = 3.886 #eV
ω1 = 0.816 #eV

f2 = 0.124
Γ2 = 0.452 #eV
ω2 = 4.481 #eV

f3 = 0.011
Γ3 = 0.065 #eV
ω3 = 8.185 #eV

f4 = 0.840
Γ4 = 0.916 #eV
ω4 = 9.083 #eV

f5 = 5.646
Γ5 = 2.419 #eV
ω5 = 20.29 #eV

Ωp = f0**.5 * ωp  #eV

def LD(ω):  #ω: eV
    ε = 1-Ωp**2/(ω*(ω+1j*Γ0))
    ε += f1*ωp**2 / ((ω1**2-ω**2)-1j*ω*Γ1)
    ε += f2*ωp**2 / ((ω2**2-ω**2)-1j*ω*Γ2)
    ε += f3*ωp**2 / ((ω3**2-ω**2)-1j*ω*Γ3)
    ε += f4*ωp**2 / ((ω4**2-ω**2)-1j*ω*Γ4)
    ε += f5*ωp**2 / ((ω5**2-ω**2)-1j*ω*Γ5)
    return ε
  
ev_min=0.1
ev_max=5
npoints=200
eV = np.logspace(np.log10(ev_min), np.log10(ev_max), npoints)
μm = 4.13566733e-1*2.99792458/eV
ε = LD(eV)
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
plt.plot(eV, -ε.real, label="-ε1")
plt.plot(eV, ε.imag, label="ε2")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.1,5])
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