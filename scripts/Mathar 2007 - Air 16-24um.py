# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-18
# Original data: Mathar 2007, https://doi.org/10.1088/1464-4258/9/5/008

################################ 16-24 μm #####################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# adjustable parameters
T = 273.15+17.5   # Temperature: K
p = 101325      # Pressure: Pa
H = 0           # Humidity: 0-100%

# model parameters
cref = [ 0.199436e-3,  0.299123e-8,  -0.214862e-10,  0.143338e-12,  0.122398e-14, -0.114628e-16] # cm^j
# something seems to be wrong with cT...
cT   = [ 0.621723e-1, -0.177074e-4,   0.152213e-6,  -0.954584-9,   -0.996706e-11,  0.921476e-13] # cm^j · K
cTT  = [-23.2409,      0.108557,     -0.102439e-2,   0.634072e-5,   0.762517e-7,  -0.675587e-9 ] # cm^j · K^2
cH   = [-0.772707e-7,  0.347237e-9,  -0.272675e-11,  0.170858e-13,  0.156889e-15, -0.150004e-17] # cm^j · %^-1
cHH  = [-0.326604e-12, 0.463606e-14, -0.364272e-16,  0.228756e-18,  0.209502e-20, -0.200547e-22] # cm^j · %^-2
cp   = [ 0.266827e-8,  0.120788e-14,  0.522646e-17,  0.783027e-19,  0.753235e-21, -0.228819e-24] # cm^j · Pa^-1
cpp  = [ 0.613675e-17, 0.585494e-22,  0.286055e-24,  0.425193e-26,  0.413455e-28, -0.812941e-32] # cm^j · Pa^-2
cTH  = [ 0.375974e-3, -0.171849e-5,   0.146704e-7,  -0.917231e-10, -0.955922e-12,  0.880502e-14] # cm^j · K · %^-1
cTp  = [ 0.778436e-6,  0.461840e-12,  0.306229e-14, -0.623183e-16, -0.161119e-18,  0.800756e-20] # cm^j · K · Pa^-1
cHp  = [-0.272614e-15, 0.304662e-18, -0.239590e-20,  0.149285e-22,  0.136086e-24, -0.130999e-26] # cm^j · %^-1 · Pa^-1
σref = 1e4/20      # cm^−1
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
    
λ = np.arange(16, 24.01, 0.05)
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
