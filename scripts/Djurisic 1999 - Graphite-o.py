# -*- coding: utf-8 -*-
# Original data: https://doi.org/10.1063/1.369370
# Author: Mikhail Polyanskiy
# Last modified: 2017-01-15

import numpy as np
import matplotlib.pyplot as plt

# LD model parameters - Normal polarization (ordinary)
ωp =   27
εinf = 1.070

f0 =   0.014
Γ0 =   6.365
ω0 =   0
α0 =   0

f1 =   0.073
Γ1 =   4.102
ω1 =   0.275
α1 =   0.505

f2 =   0.056
Γ2 =   7.328
ω2 =   3.508
α2 =   7.079

f3 =   0.069
Γ3 =   1.414
ω3 =   4.451
α3 =   0.362

f4 =   0.005
Γ4 =   0.46 # 0.046 in the original paper!
ω4 =   13.591
α4 =   7.426

f5 =   0.262
Γ5 =   1.862
ω5 =   14.226
α5 =   3.82e-4

f6 =   0.460
Γ6 =   11.922
ω6 =   15.550
α6 =   1.387

f7 =   0.200
Γ7 =   39.091
ω7 =   32.011
α7 =   28.963


def LD(ω):
    ε = εinf;
    Γ = Γ0*np.exp(-α0*((ω-ω0)/Γ0)**2)
    ε -= f0*ωp**2 / ((ω**2-ω0**2)+1j*ω*Γ)
    Γ = Γ1*np.exp(-α1*((ω-ω1)/Γ1)**2)
    ε -= f1*ωp**2 / ((ω**2-ω1**2)+1j*ω*Γ)
    Γ = Γ2*np.exp(-α2*((ω-ω2)/Γ2)**2)
    ε -= f2*ωp**2 / ((ω**2-ω2**2)+1j*ω*Γ)
    Γ = Γ3*np.exp(-α3*((ω-ω3)/Γ3)**2)
    ε -= f3*ωp**2 / ((ω**2-ω3**2)+1j*ω*Γ)
    Γ = Γ4*np.exp(-α4*((ω-ω4)/Γ4)**2)
    ε -= f4*ωp**2 / ((ω**2-ω4**2)+1j*ω*Γ)
    Γ = Γ5*np.exp(-α5*((ω-ω5)/Γ5)**2)
    ε -= f5*ωp**2 / ((ω**2-ω5**2)+1j*ω*Γ)
    Γ = Γ6*np.exp(-α6*((ω-ω6)/Γ6)**2)
    ε -= f6*ωp**2 / ((ω**2-ω6**2)+1j*ω*Γ)
    Γ = Γ7*np.exp(-α7*((ω-ω7)/Γ7)**2)
    ε -= f7*ωp**2 / ((ω**2-ω7**2)+1j*ω*Γ)
    return ε


ev_min=0.12
ev_max=40
npoints=1000
eV = np.linspace(ev_min, ev_max, npoints)
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

#plot ε vs eV
plt.figure(1)
plt.plot(eV, ε.real, label="ε1")
plt.plot(eV, ε.imag, label="ε2")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs eV
plt.figure(2)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
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