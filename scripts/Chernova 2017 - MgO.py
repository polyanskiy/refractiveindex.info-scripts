# -*- coding: utf-8 -*-

# Original data: Chernova et al. 2017, https://doi.org/10.1364/OME.7.003844

# Kramers-Kroning Integration https://doi.org/10.1366/0003702884430380

# Tauc-Lorentz oscillators https://doi.org/10.1063/1.118064

# Version history
# 2025-02-27 first version (Pavel Dmitriev)
# 2025-02-27 simplify output (Misha Polyanskiy)
#

import numpy as np
import matplotlib.pyplot as plt

UV_E = 13.7
UV_Amplitude = 47

# Tauc-Lorentz
E_TL = [7.07, 7.64]
A_TL = [472.4, 818.3]
C_TL = [5.99, 0.19]
Eg = [8.01, 7.45]

# Lorentz
Lorentz_Amplitude = 1.62
Lorentz_FWHM = 0.1
Lorentz_E = 7.6

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

# Tauc-Lorentz
for i in range(len(E_TL)):
    eps_2_osc = np.zeros(eV.shape)
    for j, e in enumerate(eV):
        if e > Eg[i]:
            eps_2_osc[j] = (1/e)*(
                    (A_TL[i]*E_TL[i]*C_TL[i]*(e-Eg[i])**2) / ((e**2 - E_TL[i]**2)**2 + C_TL[i]**2*e**2)
            )

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

#
# Lorentz-Drude
#
eps_lor = np.zeros(epsilon.shape)
for i, e in enumerate(eV):
    eps_lor[i] = 1 + Lorentz_Amplitude**2 / (Lorentz_E**2 - e**2 - 1j * e * Lorentz_FWHM)

epsilon = np.add(epsilon, eps_lor)    

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

