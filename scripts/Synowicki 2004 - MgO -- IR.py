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
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    # # Gaussian eV -- only UV
    # Gaussian_Amplitude = [0.55291, 4.2344, 6.71*10**-4, 6.18*10**-5, 3.56*10**-5]
    # Gaussian_E0 = [8.6734, 7.6359, 7.2652, 5.0889, 4.3675]
    # Gaussian_Br = [1.0501, 0.18511, 2.1479, 0.80321, 0.60936]

    # Lorentz cm-1 -- only IR
    Lorentz_Amplitude = [1702.3, 0.26993]
    Lorentz_Eg = [396.73, 643.56]
    Lorentz_Br = [1.5498, 106.86]

    # #Tauc Lorentz eV -- only UV
    # TL_A = [297.39, 2274.6]
    # TL_C = [4.4633, 0.17911]
    # TL_E0 = [9.623, 7.7091]
    # TL_Eg = [7.9862, 7.6602]

    #Poles

    eps_inf = 1

    # Simulate range
    num_points = 6000

    waveNumber = np.linspace(1, 6000, num_points, True) #cm-1

    #
    # Oscillators
    #
    eps_2 = np.zeros(waveNumber.shape)
    eps_1 = np.zeros(waveNumber.shape)

    #
    # Lorentz oscillators
    #
    for i in range(len(Lorentz_Amplitude)):
        eps_1_osc, eps_2_osc = auxfuncs.lorentz(waveNumber, Lorentz_Amplitude[i], Lorentz_Br[i], Lorentz_Eg[i])
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc


    eps_1 += eps_inf

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
