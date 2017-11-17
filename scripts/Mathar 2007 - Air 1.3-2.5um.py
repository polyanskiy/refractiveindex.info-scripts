# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-16
# Original data: Mathar 2007, https://doi.org/10.1088/1464-4258/9/5/008

############################### 1.3-2.5 μm ####################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# adjustable parameters
T = 273.15+15   # Temperature: K
p = 101325      # Pressure: Pa
H = 0           # Humidity: 0-100%

# model parameters
cref = [ 0.200192e-3,  0.113474e-9,  -0.424595e-14,  0.100957e-16, -0.293315e-20,  0.307228e-24] # cm^j
cT   = [ 0.588625e-1, -0.385766e-7,   0.888019e-10, -0.567650e-13,  0.166615e-16, -0.174845e-20] # cm^j · K
cTT  = [-3.01513,      0.406167e-3,  -0.514544e-6,   0.343161e-9,  -0.101189e-12,  0.106749e-16] # cm^j · K^2
cH   = [-0.103945e-7,  0.136858e-11, -0.171039e-14,  0.112908e-17, -0.329925e-21,  0.344747e-25] # cm^j · %^-1
cHH  = [ 0.573256e-12, 0.186367e-16, -0.228150e-19,  0.150947e-22, -0.441214e-26,  0.461209e-30] # cm^j · %^-2
cp   = [ 0.267085e-8,  0.135941e-14,  0.135295e-18,  0.818218e-23, -0.222957e-26,  0.249964e-30] # cm^j · Pa^-1
cpp  = [ 0.609186e-17, 0.519024e-23, -0.419477e-27,  0.434120e-30, -0.122445e-33,  0.134816e-37] # cm^j · Pa^-2
cTH  = [ 0.497859e-4, -0.661752e-8,   0.832034e-11, -0.551793e-14,  0.161899e-17, -0.169901e-21] # cm^j · K · %^-1
cTp  = [ 0.779176e-6,  0.396499e-12,  0.395114e-16,  0.233587e-20, -0.636441e-24,  0.716868e-28] # cm^j · K · Pa^-1
cHp  = [-0.206567e-15, 0.106141e-20, -0.149982e-23,  0.984046e-27, -0.288266e-30,  0.299105e-34] # cm^j · %^-1 · Pa^-1
σref = 1e4/2.25    # cm^−1
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
    
λ = np.arange(1.30, 2.501, 0.01)
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
