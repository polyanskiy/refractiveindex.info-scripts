# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-15
# Original data: Mathar 2007, https://doi.org/10.1088/1464-4258/9/5/008

############################### 7.5-14.1 μm ###################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# adjustable parameters
T = 273.15+15   # Temperature: K
p = 101325      # Pressure: Pa
H = 0           # Humidity: 0-100%

# model parameters
cref = [ 0.199885e-3,  0.344739e-9,  -0.273714e-12,  0.393383e-15, -0.569488e-17,  0.164556e-19] # cm^j
cT   = [ 0.593900e-1, -0.172226e-5,   0.237654e-8,  -0.381812e-11,  0.305050e-14, -0.157464e-16] # cm^j · K
cTT  = [-6.50355,      0.103830e-1,  -0.139464e-4,   0.220077e-7,  -0.272412e-10,  0.126364e-12] # cm^j · K^2
cH   = [-0.221938e-7,  0.347377e-10, -0.465991e-13,  0.735848e-16, -0.897119e-19,  0.380817e-21] # cm^j · %^-1
cHH  = [ 0.393524e-12, 0.464083e-15, -0.621764e-18,  0.981126e-21, -0.121384e-23,  0.515111e-26] # cm^j · %^-2
cp   = [ 0.266809e-8,  0.695247e-15,  0.159070e-17, -0.303451e-20, -0.661489e-22,  0.178226e-24] # cm^j · Pa^-1
cpp  = [ 0.610508e-17, 0.227694e-22,  0.786323e-25, -0.174448e-27, -0.359791e-29,  0.978307e-32] # cm^j · Pa^-2
cTH  = [ 0.106776e-3, -0.168516e-6,   0.226201e-9,  -0.356457e-12,  0.437980e-15, -0.194545e-17] # cm^j · K · %^-1
cTp  = [ 0.77368e-6,   0.216404e-12,  0.581805e-15, -0.189618e-17, -0.198869e-19,  0.589381e-22] # cm^j · K · Pa^-1
cHp  = [-0.206365e-15, 0.300234e-19, -0.426519e-22,  0.684306e-25, -0.467320e-29,  0.126117e-30] # cm^j · %^-1 · Pa^-1
σref = 1e4/10.1    # cm^−1
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

λ = np.arange(7.5, 14.101, 0.05)
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
