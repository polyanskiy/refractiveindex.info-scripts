# -*- coding: utf-8 -*-

# Original data: Bright et al. 2012, https://doi.org/10.1016/j.tsf.2012.07.037

# Version history
# 2022-12-18 first version (Mikhail Polyanskiy)
#

import numpy as np
import matplotlib.pyplot as plt

# Model parameters
A = 1.956
C = 6.73e-10 #cm^2
ν  = [187.3, 254.9, 336.9, 402.9, 506.0, 594.8] #1/cm
νp = [247.1, 373.1, 683.3, 537.8, 371.1, 118.4] #1/cm
γ  = [215.9, 45.1,  62.4,  56.6,  54.1,  26.2 ] #1/m

def M(η):  
     ε = A**2 + C*η**2 + 0*1j
     for i in range(6):
         ε += νp[i]**2 / (ν[i]**2 - η**2 - 1j*γ[i]*η)
     return ε
  
η_min=20    #1/cm
η_max=26000 #1/cm
npoints=1000
η = np.logspace(np.log10(η_min), np.log10(η_max), npoints)
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

# plot n,k - Fig. 9 of the paper
plt.figure(1)
plt.xscale('log')
plt.ylim([0,5])
plt.grid()
plt.plot(η, n, label="n")
plt.plot(η, k, label="k")
plt.xlabel('Wavenumber (1/cm)')
plt.ylabel('Optical constants')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
plt.show()


