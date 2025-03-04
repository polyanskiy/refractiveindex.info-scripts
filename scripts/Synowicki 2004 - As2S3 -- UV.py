# -*- coding: utf-8 -*-

# Original data: Synowicki et al. 2004, https://doi.org/10.1016/j.tsf.2004.02.028
# NB! Lorentz equation 1a seems to have a typo -- denominator parenthesis should be squared (?)

# Kramers-Kroning Integration https://doi.org/10.1366/0003702884430380

# Version history
# 2025-03-03 first version (Pavel Dmitriev)
#

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("TkAgg")


def generate_epsilon():
    #
    # Model parameters
    #

    # # Gaussian cm-1 -- only IR
    # Gaussian_Amplitude = [0.72443, 4.0012, 1.0615, 8.8976]
    # Gaussian_E0 = [387.09, 330.31, 324.31, 308.96]
    # Gaussian_Br = [25.814, 68.616, 33.397, 23.55]


    # Gaussian eV -- only UV
    Gaussian_Amplitude = [0.50989, 4.92*10**-5, 6.15*10**-6, 0.69671]
    Gaussian_E0 = [6.045, 2.1442, 4.4274, 4.5874]
    Gaussian_Br = [1.7501, 0.21914, 4.2673, 1.188]


    #Tauc Lorentz eV -- only UV
    TL_A = [191.03, 178.58]
    TL_E0 = [8.4056, 3.5588]
    TL_C = [3.3764, 3.674]
    TL_Eg = [7.3762, 2.4054]

    #Poles

    eps_inf = 1

    # Simulate range
    num_points = 2000
    eV = np.linspace(0.01, 15.0, num_points, True)
    dEv = eV[1] - eV[0]

    # Epsilon infinity

    epsilon_1_inf = eps_inf * np.ones(eV.shape)

    #
    # Oscillators
    #
    eps_2 = np.zeros(eV.shape)

    #
    # Gaussian oscillators -- only for UV
    #
    for i in range(len(Gaussian_E0)):
        f = (0.5 / np.sqrt(np.log(2)))
        eps_2_osc = [
            Gaussian_Amplitude[i] * np.exp(
                -((e - Gaussian_E0[i]) / (f * Gaussian_Br[i])) ** 2
            ) - \
            Gaussian_Amplitude[i] * np.exp(
                -((e + Gaussian_E0[i]) / (f * Gaussian_Br[i])) ** 2
            )
            for e in eV
        ]

        eps_2 = np.add(eps_2, eps_2_osc)

    #
    # Tauc-Lorentz oscillators -- only for UV
    #
    for i in range(len(TL_E0)):
        eps_2_osc = np.zeros(eV.shape)
        for j, e in enumerate(eV):
            if e > TL_Eg[i]:
                eps_2_osc[j] = (1 / e) * (
                        (TL_A[i] * TL_E0[i] * TL_C[i] * (e - TL_Eg[i]) ** 2) / (
                        (e ** 2 - TL_E0[i] ** 2) ** 2 + TL_C[i] ** 2 * e ** 2)
                )

        eps_2 = np.add(eps_2, eps_2_osc)

    #
    # KK integral
    #
    eps_1_osc = np.zeros(eV.shape)
    for i, e in enumerate(eV):
        prefactor = (2 / np.pi) * 2 * dEv

        maclaurin_sum = 0
        if i % 2 == 0:
            js = [z for z in range(eV.shape[-1])[1::2]]
        else:
            js = [z for z in range(eV.shape[-1])[0::2]]

        for j in js:
            maclaurin_sum += (1. / 2.) * (
                    eps_2[j] / (eV[j] - e) +
                    eps_2[j] / (eV[j] + e)
            )
        eps_1_osc[i] = prefactor * maclaurin_sum

    eps_1 = np.add(epsilon_1_inf, eps_1_osc)

    epsilon = np.asarray(eps_1 + 1j * eps_2, dtype=np.complex128)

    return eV, epsilon


if __name__ == "__main__":
    eV, epsilon = generate_epsilon()

    # Model range
    fit_points = 200
    fit_eV = np.linspace(0.073, 9.53, fit_points, True)

    n = (epsilon ** .5).real
    k = (epsilon ** .5).imag

    #
    # Interpolate to data range
    #
    n_interp = np.interp(fit_eV, eV, n)
    k_interp = np.interp(fit_eV, eV, k)

    wl_um = np.divide(1.23984193, fit_eV)

    # ============================   DATA OUTPUT   =================================
    file = open('out.txt', 'w')
    for i in range(fit_points - 1, -1, -1):
        file.write('\n        {:.4e} {:.4e} {:.4e}'.format(wl_um[i], n_interp[i], k_interp[i]))
    file.close()

    # ===============================   PLOT   =====================================
    # plot eps1 vs eV
    plt.figure(1)
    plt.plot(eV, epsilon.real)
    plt.xlabel('Photon energy (eV)')
    plt.xlim([0.073, 9.53])
    plt.ylim([0., 10.])
    plt.ylabel('epsilon1')

    # plot n vs eV
    plt.figure(2)
    plt.plot(eV, epsilon.imag)
    plt.xlabel('Photon energy (eV)')
    plt.xlim([0.073, 9.53])
    plt.ylim([0., 8.])
    plt.ylabel('epsilon2')

    # plot n,k vs μm
    plt.figure(3)
    plt.plot(eV, epsilon.real, label="eps_1")
    plt.plot(eV, epsilon.imag, label="eps_2")
    plt.legend()
    plt.xlabel('Energy (eV)')
    plt.xlim([0.073, 9.53])
    plt.ylim([0., 10.])

    # plot n,k vs μm
    plt.figure(4)
    plt.plot(wl_um, n_interp, label="n")
    plt.plot(wl_um, k_interp, label="k")
    plt.legend()
    plt.xlabel('Wavelength (μm)')

    plt.show()
