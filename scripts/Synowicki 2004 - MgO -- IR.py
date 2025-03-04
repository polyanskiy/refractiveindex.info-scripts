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

    # Gaussian eV -- only UV
    Gaussian_Amplitude = [0.55291, 4.2344, 6.71*10**-4, 6.18*10**-5, 3.56*10**-5]
    Gaussian_E0 = [8.6734, 7.6359, 7.2652, 5.0889, 4.3675]
    Gaussian_Br = [1.0501, 0.18511, 2.1479, 0.80321, 0.60936]

    # Lorentz cm-1 -- only IR
    Lorentz_Amplitude = [1702.3, 0.26993]
    Lorentz_Eg = [396.73, 643.56]
    Lorentz_Br = [1.5498, 106.86]

    #Tauc Lorentz eV -- only UV
    TL_A = [297.39, 2274.6]
    TL_C = [4.4633, 0.17911]
    TL_E0 = [9.623, 7.7091]
    TL_Eg = [7.9862, 7.6602]

    #Poles

    eps_inf = 1

    # Simulate range
    num_points = 5601
    # eV = np.linspace(0.001, 0.9, num_points, True)
    # dEv = eV[1] - eV[0]

    waveNumber = np.linspace(300, 5900, num_points, True) #cm-1
    dcm = waveNumber[1] - waveNumber[0]


    # Epsilon infinity

    epsilon_1_inf = eps_inf * np.ones(waveNumber.shape)

    #
    # Oscillators
    #
    eps_2 = np.zeros(waveNumber.shape)

    #
    # Lorentz oscillators -- Only for IR
    #
    eps_1_lor = np.zeros(shape=waveNumber.shape)
    for i in range(len(Lorentz_Eg)):
        eps_2_osc = [
            Lorentz_Amplitude[i] * Lorentz_Br[i] * Lorentz_Eg[i] * e /
            ( (Lorentz_Eg[i]**2 - e**2)**2 + Lorentz_Br[i]**2*e**2)
            for e in waveNumber
        ]

        eps_2 = np.add(eps_2, eps_2_osc)

        eps_1_osc = [
            Lorentz_Amplitude[i] * Lorentz_Br[i] * Lorentz_Eg[i] * (Lorentz_Eg[i]**2 - e**2) /
            ( (Lorentz_Eg[i]**2 - e**2)**2 + Lorentz_Br[i]**2*e**2)
            for e in waveNumber
        ]

        eps_1_lor = np.add(eps_1_lor, eps_1_osc)

    eps_1 = np.add(epsilon_1_inf, eps_1_lor)

    epsilon = np.asarray(eps_1 + 1j * eps_2, dtype=np.complex128)

    return waveNumber, epsilon


if __name__ == "__main__":
    waveNumber, epsilon = generate_epsilon()

    # Model range
    fit_points = 200
    wl_um = np.linspace(1.7, 33, fit_points, True)
    fit_wavenumber = np.divide(1e4, wl_um)


    n = (epsilon ** .5).real
    k = (epsilon ** .5).imag

    #
    # Interpolate to data range
    #
    n_interp = np.interp(fit_wavenumber, waveNumber, n)
    k_interp = np.interp(fit_wavenumber, waveNumber, k)

    e2_interp = np.interp(fit_wavenumber, waveNumber, epsilon.imag)
    e1_interp = np.interp(fit_wavenumber, waveNumber, epsilon.real)

    # ============================   DATA OUTPUT   =================================
    file = open('out.txt', 'w')
    for i in range(fit_points - 1, -1, -1):
        file.write('\n        {:.4e} {:.4e} {:.4e}'.format(wl_um[i], n_interp[i], k_interp[i]))
    file.close()

    # ===============================   PLOT   =====================================
    # plot eps1 vs eV
    plt.figure(1)
    plt.plot(wl_um, e1_interp)
    plt.xlabel('Wavelength (um)')
    plt.ylabel('epsilon1')

    # plot eps2 vs eV
    plt.figure(2)
    plt.plot(wl_um, e2_interp)
    plt.xlabel('Wavelength (um)')
    plt.ylabel('epsilon2')

    # plot n,k vs μm
    plt.figure(3)
    plt.plot(wl_um, n_interp, label="n")
    plt.plot(wl_um, k_interp, label="k")
    plt.legend()
    plt.xlabel('Wavelength (μm)')

    plt.show()
