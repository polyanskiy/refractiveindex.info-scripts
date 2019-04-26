# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-02
# Original data: Kitamura et al. 2007, https://doi.org/10.1364/AO.46.008118

import numpy as np
from scipy.special import dawsn
import matplotlib.pyplot as plt
π = np.pi

# Model parameters
εinf = 2.1232
α =  [3.7998, .46089, 1.2520, 7.8147, 1.0313,  5.3757, 6.3305, 1.2948]
η0 = [1089.7, 1187.7, 797.78, 1058.2, 446.13,  443.00, 465.80, 1026.7] #1/cm
σ =  [31.454, 100.46, 91.601, 63.153, 275.111, 45.220, 22.680, 232.14] #1/cm


D = lambda x: dawsn(x)

def M(η):  
     ε = εinf
     for i in range(0,8):
         gc = α[i]*np.exp(-4*np.log(2)*((η-η0[i])/σ[i])**2) - α[i]*np.exp(-4*np.log(2)*((η+η0[i])/σ[i])**2)
         gckkg = 2*α[i]/np.sqrt(π) * (D(2*np.log(2)**.5*(η+η0[i])/σ[i]) - D(2*np.log(2)**.5*(η-η0[i])/σ[i])) 
         ε += gckkg + 1j*gc
     return ε
  
η_min=10000/50 #μm -> 1/cm
η_max=10000/7  #μm -> 1/cm
npoints=200
η = np.linspace(η_min, η_max, npoints)
λ = 10000/η #1/cm -> μm
ε = M(η)
n = (ε**.5).real
k = (ε**.5).imag


#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(npoints-1, -1, -1):
    file.write('\n        {:.4e} {:.4e} {:.4e}'.format(λ[i],n[i],k[i]))
file.close()
    
    
#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

plt.figure(1)
plt.plot(η, ε.real, label="ε1")
plt.plot(η, ε.imag, label="ε2")
plt.xlabel('Wavenumber (1/cm)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n vs wavelength (Fig. 3 top in the paper)
plt.figure(2)
plt.plot(λ, n, label="n")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n')
plt.xlim([1,15])
plt.ylim([0,4])

#plot k vs wavelength  (Fig. 3 bottom in the paper)
plt.figure(3)
plt.plot(λ, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('k')
plt.yscale('log')
plt.xlim([1,15])
plt.ylim([1e-8 ,10])

#plot n vs wavelength (Fig. 4 top in the paper)
plt.figure(4)
plt.plot(λ, n, label="n")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n')
plt.xlim([15,100])
plt.ylim([0,3])


#plot k vs wavelength (Fig. 4 bottom in the paper)
plt.figure(5)
plt.plot(λ, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('k')
plt.yscale('log')
plt.xlim([15,100])
plt.ylim([1e-2 ,10])
