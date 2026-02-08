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


def generate_ir_oscillators(num_points=6000, wavenum_min=300., wavenum_max=5900.):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    # Lorentz cm-1 -- only IR
    Lorentz_Amplitude = [1702.3, 0.26993]
    Lorentz_Eg = [396.73, 643.56]
    Lorentz_Br = [1.5498, 106.86]

    # Simulate range
    waveNumber = np.linspace(wavenum_min, wavenum_max, num_points, True)  # cm-1

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

    return waveNumber, eps_1, eps_2


def generate_uv_oscillators(num_points=10000, min_eV=0.01, max_eV=30.0):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    # Gaussian eV -- only UV
    Gaussian_Amplitude = [0.55291, 4.2344, 6.71*10**-4, 6.18*10**-5, 3.56*10**-5]
    Gaussian_E0 = [8.6734, 7.6359, 7.2652, 5.0889, 4.3675]
    Gaussian_Br = [1.0501, 0.18511, 2.1479, 0.80321, 0.60936]

    #Tauc Lorentz eV -- only UV
    TL_A = [297.39, 2274.6]
    TL_C = [4.4633, 0.17911]
    TL_E0 = [9.623, 7.7091]
    TL_Eg = [7.9862, 7.6602]

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
        eps_1_osc, eps_2_osc = auxfuncs.gaussian(eV, Gaussian_E0[i], Gaussian_Amplitude[i], Gaussian_Br[i])
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    #
    # Tauc-Lorentz oscillators -- only for UV
    #
    for i in range(len(TL_E0)):
        eps_1_osc, eps_2_osc = auxfuncs.tauc_lorentz(eV, TL_E0[i], TL_A[i], TL_C[i], TL_Eg[i])
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    epsilon = np.asarray(eps_1 + 1j * eps_2, dtype=np.complex128)

    return eV, eps_1, eps_2


def generate_epsilon(fit_points=1000, min_um=1.7, max_um=33, gen_points=10000, min_ev=0.01, max_ev=30.0, min_wavenum=300., max_wavenum=5900.):
    waveNumber, osc_ir_1, osc_ir_2 = generate_ir_oscillators(num_points=gen_points, wavenum_min=min_wavenum, wavenum_max=max_wavenum)
    eV, osc_uv_1, osc_uv_2 = generate_uv_oscillators(num_points=gen_points, min_eV=min_ev, max_eV=max_ev)

    # Model range
    wl_um = np.linspace(max_um, min_um, fit_points, True)
    fit_wavenumber = np.divide(1e4, wl_um)
    fit_ev = np.divide(1.23984193, wl_um)

    ir_osc2_interp = np.interp(fit_wavenumber, waveNumber, osc_ir_2)
    ir_osc1_interp = np.interp(fit_wavenumber, waveNumber, osc_ir_1)
    uv_osc2_interp = np.interp(fit_ev, eV, osc_uv_2)
    uv_osc1_interp = np.interp(fit_ev, eV, osc_uv_1)

    eps_inf = 1.

    epsilon = eps_inf + ir_osc1_interp + uv_osc1_interp + 1j*ir_osc2_interp + 1j*uv_osc2_interp
    return wl_um, epsilon


if __name__ == "__main__":
    fit_points = 500

    #
    # IR
    #
    min_um = 1.7
    max_um = 33.
    wl_um, epsilon = generate_epsilon(fit_points=fit_points, min_um=min_um, max_um=max_um)

    # #
    # # UV
    # #
    # min_ev = 0.73
    # max_ev = 9.53
    # min_um = np.divide(1.23984193, min_ev)
    # max_um = np.divide(1.23984193, max_ev)
    # wl_um, epsilon = generate_epsilon(fit_points=fit_points, min_um=min_um, max_um=max_um)

    n = (epsilon ** .5).real
    k = (epsilon ** .5).imag


    # ============================   DATA OUTPUT   =================================
    file = open('out.txt', 'w')
    for i in range(fit_points - 1, -1, -1):
        file.write('\n        {:.4e} {:.4e} {:.4e}'.format(wl_um[i], n[i], k[i]))
    file.close()

    # ===============================   PLOT   =====================================
    # plot eps1 vs eV
    plt.figure(1)
    plt.plot(wl_um, epsilon.real)
    plt.xlabel('Wavelength (um)')
    plt.ylim([-20, 20])
    plt.ylabel('epsilon1')

    # plot eps2 vs eV
    plt.figure(2)
    plt.plot(wl_um, epsilon.imag)
    plt.xlabel('Wavelength (um)')
    plt.ylim([0, 50])
    plt.ylabel('epsilon2')

    # plot n,k vs μm
    plt.figure(3)
    plt.plot(wl_um, n, label="n")
    plt.plot(wl_um, k, label="k")
    plt.legend()
    plt.xlabel('Wavelength (μm)')

    plt.show()