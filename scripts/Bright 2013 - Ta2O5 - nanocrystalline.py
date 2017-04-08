# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-06
# Original data: Bright et al. 2013, https://doi.org/10.1063/1.4819325

import numpy as np
import matplotlib.pyplot as plt

# model parameters
ωj =  [0,     91,  214, 324, 530,  842]
ωpj = [6490,  260, 844, 391, 1019, 372]
γj =  [6.5e5, 74,  61,  73,  142,  114]
A = 2.10
B = 0.04


def Epsilon(ω):
    ε = 0
    for j in range (0, len(ωj)):
        ε += ωpj[j]**2 / (ωj[j]**2 - ω**2 - 1j*γj[j]*ω)
    return ε


ω_min = 10   # cm^-1
ω_max = 20000 # cm^-1
npoints = 500
ω = np.logspace(np.log10(ω_min), np.log10(ω_max), npoints)
μm = 10000/ω
ε = (A+B/μm**2)**2 + Epsilon(ω)
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