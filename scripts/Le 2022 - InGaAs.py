# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy (2017)
# modified by S.Vassant for InGaAs data
# Last modified: 2022-10-23
# Original data: Le et al. 2022, https://doi.org/10.1364/OME.455445

import numpy as np
from scipy.special import dawsn
import matplotlib.pyplot as plt
π = np.pi

# Model parameters
εinf = 11.31
α =  [18.65, 21.23, 10.13]
η0 = [225.31, 245.23, 254.14] #1/cm
σ =  [12.11,23.08,7.16] #1/cm


D = lambda x: dawsn(x)

def M(η):  
     ε = εinf
     for i in range(0,3):
         gc = α[i]*np.exp(-4*np.log(2)*((η-η0[i])/σ[i])**2) - α[i]*np.exp(-4*np.log(2)*((η+η0[i])/σ[i])**2)
         gckkg = 2*α[i]/np.sqrt(π) * (D(2*np.log(2)**.5*(η+η0[i])/σ[i]) - D(2*np.log(2)**.5*(η-η0[i])/σ[i]))
         ε += gckkg + 1j*gc
     return ε
  
η_min=150 #μm -> 1/cm
η_max=350  #μm -> 1/cm
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

# plot Re(ε), Im(ε), Fig. 4 of the paper
plt.figure(1)
plt.plot(η, ε.real, label="ε1")
plt.plot(η, ε.imag, label="ε2")
plt.xlabel('Wavenumber (1/cm)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
plt.show()

# plot n,k (λ)
plt.figure(2)
plt.plot(λ, n, label="n")
plt.plot(λ, k, label="k")
plt.xlabel('Wavelength μm')
plt.ylabel('n,k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
plt.show()
