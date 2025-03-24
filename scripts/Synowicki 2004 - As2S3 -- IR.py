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


def generate_epsilon(num_points=6000, wavenum_min=300., wavenum_max=5900., kk="ML"):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    # Gaussian cm-1 -- only IR
    Gaussian_Amplitude = [0.72443, 4.0012, 1.0615, 8.8976]
    Gaussian_E0 = [387.09, 330.31, 324.31, 308.96]
    Gaussian_Br = [25.814, 68.616, 33.397, 23.55]


    # # Gaussian eV -- only UV
    # Gaussian_Amplitude = [0.50989, 4.92*10**-5, 6.15*10**-6, 0.69671]
    # Gaussian_E0 = [6.045, 2.1442, 4.4274, 4.5874]
    # Gaussian_Br = [1.7501, 0.21914, 4.2673, 1.188]
    #
    #
    # #Tauc Lorentz eV -- only UV
    # TL_A = [191.03, 178.58]
    # TL_E0 = [8.4056, 3.5588]
    # TL_C = [3.3764, 3.674]
    # TL_Eg = [7.3762, 2.4054]

    #Poles

    eps_inf = 1.

    # Simulate range
    waveNumber = np.linspace(wavenum_min, wavenum_max, num_points, True)  # cm-1

    #
    # Oscillators
    #
    eps_2 = np.zeros(waveNumber.shape)
    eps_1 = np.zeros(waveNumber.shape)

    #
    # Gaussian oscillators -- only for IR
    #
    for i in range(len(Gaussian_E0)):
        eps_1_osc, eps_2_osc = auxfuncs.gaussian(waveNumber, Gaussian_E0[i], Gaussian_Amplitude[i], Gaussian_Br[i], kk)
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc


    eps_1 += eps_inf

    epsilon = np.asarray(eps_1 + 1j * eps_2, dtype=np.complex128)

    return waveNumber, epsilon


if __name__ == "__main__":
    waveNumber, epsilon = generate_epsilon()
    model_um = np.divide(waveNumber, 1e4)

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
    #plt.ylim([0, 12])
    plt.ylabel('epsilon1')

    # plot eps2 vs eV
    plt.figure(2)
    plt.plot(wl_um, e2_interp)
    plt.xlabel('Wavelength (um)')
    plt.ylim([0, 15])
    plt.ylabel('epsilon2')

    # plot n,k vs μm
    plt.figure(3)
    plt.plot(wl_um, n_interp, label="n")
    plt.plot(wl_um, k_interp, label="k")
    plt.legend()
    plt.xlabel('Wavelength (μm)')

    plt.show()
