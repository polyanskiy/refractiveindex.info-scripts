# -*- coding: utf-8 -*-

# Original data: Chernova et al. 2017, https://doi.org/10.1364/OME.7.003844

# Kramers-Kroning Integration https://doi.org/10.1366/0003702884430380

# Version history
# 2025-02-19 first version (Pavel Dmitriev)
# 2025-02-27 simplify output (Misha Polyanskiy)
#

import numpy as np
import matplotlib.pyplot as plt

# Model parameters
E = [6.3, 6.96, 8.62, 11.47]
Amplitude = [2.95, 3.26, 2.13, 5.96]
FWHM = [0.78, 0.91, 1.78, 4.65]

UV_E = 12.
UV_Amplitude = 148

eps_inf = 1

# Simulate range
num_points = 1000
eV = np.linspace(0.1, 15.0, num_points, True)
dEv = eV[1]-eV[0]

# Model range
fit_points = 100
fit_eV = np.linspace(0.74, 8.8, fit_points, True)

epsilon_1_inf = eps_inf * np.ones(eV.shape)
epsilon_1_UV = [UV_Amplitude/(UV_E**2 - e**2) for e in eV]

eps_1 = np.add(eps_inf, epsilon_1_UV)
eps_2 = np.zeros(eV.shape)

for i in range(len(E)):
    Br = FWHM[i] / (2*np.sqrt(np.log(2)))
    eps_2_osc = [
        Amplitude[i]*np.exp(
            -((e-E[i])/Br)**2
        ) - \
        Amplitude[i]*np.exp(
            -((e+E[i])/Br)**2
        )
        for e in eV
    ]

    eps_2 = np.add(eps_2, eps_2_osc)

    eps_1_osc = np.zeros(eV.shape)
    for i, e in enumerate(eV):
        prefactor = (2/np.pi)*2*dEv

        maclaurin_sum = 0
        if i % 2 == 0:
            js = [z for z in range(eV.shape[-1])[1::2]]
        else:
            js = [z for z in range(eV.shape[-1])[0::2]]

        for j in js:
            maclaurin_sum += (1./2.)*(
                    eps_2_osc[j]/(eV[j]-e) +
                    eps_2_osc[j]/(eV[j]+e)
            )
        eps_1_osc[i] = prefactor*maclaurin_sum
    eps_1 = np.add(eps_1, eps_1_osc)


epsilon = eps_1 + 1j * eps_2
n = (epsilon**.5).real
k = (epsilon**.5).imag

#
# Interpolate to data range
#
n_interp = np.interp(fit_eV, eV, n)
k_interp = np.interp(fit_eV, eV, k)

wl_um = np.divide(1.23984193, fit_eV)

#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(fit_points-1, -1, -1):
    file.write('\n        {:.4e} {:.4e} {:.4e}'.format(wl_um[i],n_interp[i],k_interp[i]))
file.close()
    
    
#===============================   PLOT   =====================================
#plot n vs eV
plt.figure(1)
plt.plot(fit_eV, n_interp)
plt.xlabel('Photon energy (eV)')
plt.ylabel('n')

#plot n vs eV
plt.figure(2)
plt.plot(fit_eV, k_interp)
plt.xlabel('Photon energy (eV)')
plt.ylabel('k')

#plot n,k vs μm
plt.figure(3)
plt.plot(wl_um, n_interp, label="n")
plt.plot(wl_um, k_interp, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n, k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
