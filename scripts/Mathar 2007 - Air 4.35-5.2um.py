# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-17
# Original data: Mathar 2007, https://doi.org/10.1088/1464-4258/9/5/008

############################### 4.35-5.2 μm ###################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# adjustable parameters
T = 273.15+15   # Temperature: K
p = 101325      # Pressure: Pa
H = 0           # Humidity: 0-100%

# model parameters
cref = [ 0.200020e-3,  0.275346e-9,   0.325702e-12, -0.693603e-14,  0.285610e-17,  0.338758e-18] # cm^j
cT   = [ 0.590035e-1, -0.375764e-6,   0.134585e-9,   0.124316e-11,  0.508510e-13, -0.189245e-15] # cm^j · K
cTT  = [-4.09830,      0.250037e-2,   0.275187e-6,  -0.653398e-8,  -0.310589e-9,   0.127747e-11] # cm^j · K^2
cH   = [-0.140463e-7,  0.839350e-11, -0.190929e-14, -0.121399e-16, -0.898863e-18,  0.364662e-20] # cm^j · %^-1
cHH  = [ 0.543605e-12, 0.112802e-15, -0.229979e-19, -0.191450e-21, -0.120352e-22,  0.500955e-25] # cm^j · %^-2
cp   = [ 0.266898e-8,  0.273629e-14,  0.463466e-17, -0.916894e-23,  0.136685e-21,  0.413687e-23] # cm^j · Pa^-1
cpp  = [ 0.610706e-17, 0.116620e-21,  0.244736e-24, -0.497682e-26,  0.742024e-29,  0.224625e-30] # cm^j · Pa^-2
cTH  = [ 0.674488e-4, -0.406775e-7,   0.289063e-11,  0.819898e-13,  0.468386e-14, -0.191182e-16] # cm^j · K · %^-1
cTp  = [ 0.778627e-6,  0.593296e-12,  0.145042e-14,  0.489815e-17,  0.327941e-19,  0.128020e-21] # cm^j · K · Pa^-1
cHp  = [-0.211676e-15, 0.487921e-20, -0.682545e-23,  0.942802e-25, -0.946422e-27, -0.153682e-29] # cm^j · %^-1 · Pa^-1
σref = 1e4/4.8     # cm^−1
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
    
λ = np.arange(4.35, 5.201, 0.01)
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
