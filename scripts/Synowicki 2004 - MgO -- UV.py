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


def generate_epsilon(num_points=10000, min_eV=0.01, max_eV=30.0, kk="ML"):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    # Gaussian eV -- only UV
    Gaussian_Amplitude = [0.55291, 4.2344, 6.71*10**-4, 6.18*10**-5, 3.56*10**-5]
    Gaussian_E0 = [8.6734, 7.6359, 7.2652, 5.0889, 4.3675]
    Gaussian_Br = [1.0501, 0.18511, 2.1479, 0.80321, 0.60936]

    # # Lorentz cm-1 -- only IR
    # Lorentz_Amplitude = [1702.3, 0.26993]
    # Lorentz_Eg = [396.73, 643.56]
    # Lorentz_Br = [1.5498, 106.86]

    #Tauc Lorentz eV -- only UV
    TL_A = [297.39, 2274.6]
    TL_C = [4.4633, 0.17911]
    TL_E0 = [9.623, 7.7091]
    TL_Eg = [7.9862, 7.6602]

    #Poles

    eps_inf = 1.

    # Simulate range
    eV = np.linspace(min_eV, max_eV, num_points, True)


    #
    # Oscillators
    #
    eps_2 = np.zeros(eV.shape)
    eps_1 = np.zeros(eV.shape)

    #
    # Gaussian oscillators -- only for UV
    #
    for i in range(len(Gaussian_E0)):
        eps_1_osc, eps_2_osc = auxfuncs.gaussian(eV, Gaussian_E0[i], Gaussian_Amplitude[i], Gaussian_Br[i], kk)
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    #
    # Tauc-Lorentz oscillators -- only for UV
    #
    for i in range(len(TL_E0)):
        eps_1_osc, eps_2_osc = auxfuncs.taucLorentz_KK(eV, TL_E0[i], TL_A[i], TL_C[i], TL_Eg[i], kk)
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    eps_1 += eps_inf

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
    plt.ylim([2., 9.])
    plt.ylabel('epsilon1')

    # plot n vs eV
    plt.figure(2)
    plt.plot(eV, epsilon.imag)
    plt.xlabel('Photon energy (eV)')
    plt.xlim([0.073, 9.53])
    plt.ylim([0., 5.])
    plt.ylabel('epsilon2')

    # plot n,k vs μm
    plt.figure(3)
    plt.plot(eV, epsilon.real, label="eps_1")
    plt.plot(eV, epsilon.imag, label="eps_2")
    plt.legend()
    plt.xlabel('Energy (eV)')
    plt.xlim([0.073, 9.53])
    plt.ylim([0., 9.])

    # plot n,k vs μm
    plt.figure(4)
    plt.plot(wl_um, n_interp, label="n")
    plt.plot(wl_um, k_interp, label="k")
    plt.legend()
    plt.xlabel('Wavelength (μm)')

    plt.show()
