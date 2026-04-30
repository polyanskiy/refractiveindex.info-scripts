# Original data: Ferrini 2002, https://doi.org/10.1140/epjb/e2002-00177-x

# CHANGELOG
# 2026-04-29 [Misha Polyanskiy] original version

###############################################################################

import numpy as np
import matplotlib.pyplot as plt


def n_func(λ):
    # complex n (n+ik)

    ν = 10000/λ       # 1/cm
    
    ε_inf = 9.4
    
    ν_TO1 = 329.1     # 1/cm
    S1    = 526.2     # 1/cm
    Γ1    = 9.92      # 1/cm
    
    ν_TO2 = 371       # 1/cm
    S2    = 113.5     # 1/cm
    Γ2    = 8.1       # 1/cm
    
    ν_p   = 235.2     # 1/cm
    γ     = 20        # 1/cm
    
    ε = (
        ε_inf
        + S1**2 / (ν_TO1**2 - ν**2 - 1j*Γ1*ν) # correct last sign for positive k convention
        + S2**2 / (ν_TO2**2 - ν**2 - 1j*Γ2*ν) # correct last sign for positive k convention
        - ν_p**2 / (ν**2 + 1j*γ*ν)
    )

    return np.sqrt(ε)


#============================   CALCULATIONS   ================================
eV_min = 0.01
eV_max = 0.08
npoints=256
eV = np.linspace(eV_max, eV_min, npoints)

h = 4.135667662e-15 # eV·s
c = 2.99792458e14   # μm/s

μm = h * c / eV

complex_n = n_func(μm)
n = np.real(complex_n)
k = np.imag(complex_n)

#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(npoints):
    file.write('\n        {:.6f} {:.6f} {:.4e}'.format(μm[i],n[i],k[i]))
file.close()
    
    
#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

#plot n vs eV
plt.figure(1)
plt.plot(eV, n, label='n')
plt.xlabel('Photon energy (eV)')
plt.ylabel('n')
plt.show()

#plot k vs eV
plt.figure(2)
plt.plot(eV, k, label='k')
plt.xlabel('Photon energy (eV)')
plt.ylabel('k')
plt.show()

#plot n vs μm
plt.figure(3)
plt.plot(μm, n, label='n')
plt.xlabel('Wavelength (μm)')
plt.ylabel('n')
plt.show()

#plot k vs μm
plt.figure(4)
plt.plot(μm, k, label='k')
plt.xlabel('Wavelength (μm)')
plt.ylabel('k')
plt.show()
