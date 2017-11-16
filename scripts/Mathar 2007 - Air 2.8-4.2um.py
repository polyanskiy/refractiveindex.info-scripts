# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-15
# Original data: Mathar 2007, https://doi.org/10.1088/1464-4258/9/5/008

############################### 2.8-4.2 μm ####################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# adjustable parameters
T = 273.15+15   # Temperature: K
p = 101325      # Pressure: Pa
H = 0           # Humidity: 0-100%

# model parameters
cref = [ 0.200049e-3,  0.145221e-9,   0.250951e-12, -0.745834e-15, -0.161432e-17,  0.352780e-20] # cm^j
cT   = [ 0.588432e-1, -0.825182e-7,   0.137982e-9,   0.352420e-13, -0.730651e-15, -0.167911e-18] # cm^j · K
cTT  = [-3.13579,      0.694124e-3,  -0.500604e-6,  -0.116668e-8,   0.209644e-11,  0.591037e-14] # cm^j · K^2
cH   = [-0.108142e-7,  0.230102e-11, -0.154652e-14, -0.323014e-17,  0.630616e-20,  0.173880e-22] # cm^j · %^-1
cHH  = [ 0.586812e-12, 0.312198e-16, -0.197792e-19, -0.461945e-22,  0.788398e-25,  0.245580e-27] # cm^j · %^-2
cp   = [ 0.266900e-8,  0.168162e-14,  0.353075e-17, -0.963455e-20, -0.223079e-22,  0.453166e-25] # cm^j · Pa^-1
cpp  = [ 0.608860e-17, 0.461560e-22,  0.184282e-24, -0.524471e-27, -0.121299e-29,  0.246512e-32] # cm^j · Pa^-2
cTH  = [ 0.517962e-4, -0.112149e-7,   0.776507e-11,  0.172569e-13, -0.320582e-16, -0.899435e-19] # cm^j · K · %^-1
cTp  = [ 0.778638e-6,  0.446396e-12,  0.784600e-15, -0.195151e-17, -0.542083e-20,  0.103530e-22] # cm^j · K · Pa^-1
cHp  = [-0.217243e-15, 0.104747e-20, -0.523689e-23,  0.817386e-26,  0.309913e-28, -0.363491e-31] # cm^j · %^-1 · Pa^-1
σref = 1e4/3.4     # cm^−1
Tref = 273.15+17.5 # K
pref = 75000       # Pa
Href = 10          #%

# model
def n(λ):
    σ = 1e4/λ # cm^-1
    n = 1
    for j in range(0, 6):
        n += ( cref[j] + cT[j]*(1/T-1/Tref) + cTT[j]*(1/T-1/Tref)**2
            + cH[j]*(H-Href) + cHH[j]*(H-Href)**2
            + cp[j]*(p-pref) + cpp[j]*(p-pref)**2
            + cTH[j]*(1/T-1/Tref)*(H-Href)
            + cTp[j]*(1/T-1/Tref)*(p-pref)
            + cHp[j]*(H-Href)*(p-pref) ) * (σ-σref)**j   
    return n 


# output - modify code below the line to match your needs
###############################################################################
    
λ = np.arange(2.8, 4.2, 0.01)
n = n(λ)

# write data file
file = open('out.txt', 'w')
for i in range(0, len(λ)):
    file.write('\n        {:.2f} {:.12f}'.format(λ[i],n[i]))
file.close()
    
#plot n vs μm
plt.rc('font', family='Arial', size='14')
plt.figure()
plt.plot(λ, n-1)
plt.xlabel('Wavelength (μm)')
plt.ylabel('n-1')
