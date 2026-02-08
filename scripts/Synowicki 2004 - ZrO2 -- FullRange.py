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


def generate_uv_oscillators(num_points=10000, min_eV=0.01, max_eV=30.0):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")
    #
    # Model parameters
    #

    #Tauc Lorentz eV -- only UV
    TL_A = [98.541, 145.09, 255.06]
    TL_E0 = [6.1261, 7.2693, 8.7795]
    TL_C = [1.2042, 2.2904, 2.8711]
    TL_Eg = [5.124, 5.3382, 6.8865]

    # Simulate range
    eV = np.linspace(min_eV, max_eV, num_points, True)

    #
    # Oscillators
    #
    eps_2 = np.zeros(eV.shape)
    eps_1 = np.zeros(eV.shape)

    #
    # Tauc-Lorentz oscillators -- only for UV
    #
    for i in range(len(TL_E0)):
        eps_1_osc, eps_2_osc = auxfuncs.tauc_lorentz(eV, TL_E0[i], TL_A[i], TL_C[i], TL_Eg[i])
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    return eV, eps_1, eps_2


def generate_ir_oscillators(num_points=6000, wavenum_min=300., wavenum_max=5900.):
    auxfuncs = __import__("Synowicki 2004 - Aux funcs")

    #
    # Model parameters
    #

    # Gaussian cm-1 -- IR
    Gaussian_Amplitude = np.asarray([24.263,    7.32,   32.304, 0.2135, 1.1005, 1.6355, 0.4633])
    Gaussian_E0 = np.asarray(       [69.044,    331.92, 339.2,  350.77, 476.18, 539.94, 709.72])
    Gaussian_Br = np.asarray(       [322.13,    43.363, 134.76, 907.53, 68.217, 122.4,  99.048])

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
        eps_1_osc, eps_2_osc = auxfuncs.gaussian(waveNumber, Gaussian_E0[i], Gaussian_Amplitude[i], Gaussian_Br[i])
        eps_2 += eps_2_osc
        eps_1 += eps_1_osc

    return waveNumber, eps_1, eps_2


def generate_epsilon(fit_points=1000, lin_wavelength=True, min_um=1.7, max_um=33, gen_points=10000, min_ev=0.01, max_ev=30.0, min_wavenum=300., max_wavenum=5900., kk="ML"):
    waveNumber, osc_ir_1, osc_ir_2 = generate_ir_oscillators(num_points=gen_points, wavenum_min=min_wavenum, wavenum_max=max_wavenum)
    eV, osc_uv_1, osc_uv_2 = generate_uv_oscillators(num_points=gen_points, min_eV=min_ev, max_eV=max_ev)

    # Model range
    if lin_wavelength:
        wl_um = np.linspace(max_um, min_um, fit_points, True)
        fit_wavenumber = np.divide(1e4, wl_um)
        fit_ev = np.divide(1.23984193, wl_um)
    else:
        fit_wavenumber = np.linspace(1e4/max_um, 1e4/min_um, fit_points, True)
        wl_um = np.divide(1e4, fit_wavenumber)
        fit_ev = np.divide(1.23984193, wl_um)

    ir_osc2_interp = np.interp(fit_wavenumber, waveNumber, osc_ir_2)
    ir_osc1_interp = np.interp(fit_wavenumber, waveNumber, osc_ir_1)
    uv_osc2_interp = np.interp(fit_ev, eV, osc_uv_2)
    uv_osc1_interp = np.interp(fit_ev, eV, osc_uv_1)

    eps_inf = 1.

    epsilon = eps_inf + ir_osc1_interp + uv_osc1_interp + 1j*ir_osc2_interp + 1j*uv_osc2_interp
    return wl_um, epsilon

if __name__ == "__main__":
    fit_points = 200

    #
    # IR
    #
    min_um = 1.7
    max_um = 33.
    wl_um, epsilon = generate_epsilon(fit_points=fit_points, lin_wavelength=True, min_um=min_um, max_um=max_um)

    # #
    # # UV
    # #
    # min_ev = 0.73
    # max_ev = 9.53
    # min_um = np.divide(1.23984193, max_ev)
    # max_um = np.divide(1.23984193, min_ev)
    #
    # wl_um, epsilon = generate_epsilon(fit_points=fit_points, lin_wavelength=False, min_um=min_um, max_um=max_um)

    n = (epsilon ** .5).real
    k = (epsilon ** .5).imag

    for i in range(len(k)):
        if k[i] < 1e-10:
            k[i] = 0.

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
