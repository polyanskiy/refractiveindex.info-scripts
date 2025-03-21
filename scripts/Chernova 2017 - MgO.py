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

import matplotlib
matplotlib.use("TkAgg")

def generate_epsilon():
    auxfuncs = __import__("Chernova 2017 - Aux Funcs")

    UV_E = 13.7
    UV_Amplitude = 47

    # Tauc-Lorentz
    E_TL = [7.07, 7.64]
    A_TL = [472.4, 818.3]
    C_TL = [5.99, 0.19]
    Eg_TL = [8.01, 7.45]

    # Lorentz
    Lorentz_Amplitude = [1.62]
    Lorentz_FWHM = [0.1]
    Lorentz_Eg = [7.6]

    eps_inf = 1

    # Simulate range
    num_points = 2000
    eV = np.linspace(0.1, 100.0, num_points, True)  # long TL tail requires model to go so high
    eps1 = np.zeros(eV.shape)
    eps2 = np.zeros(eV.shape)


    plt.figure(0)
    # Tauc-Lorentz
    for i in range(len(E_TL)):
        eps_1_TL, eps_2_TL = auxfuncs.taucLorentz_KK(eV, E_TL[i], A_TL[i], C_TL[i], Eg_TL[i])
        eps1 += eps_1_TL
        eps2 += eps_2_TL

        plt.plot(eV, eps_2_TL)
    #
    # Lorentz oscillators
    #
    for i in range(len(Lorentz_Eg)):
        eps_1_lor, eps_2_lor = auxfuncs.lorentz(eV, Lorentz_Amplitude[i], Lorentz_FWHM[i], Lorentz_Eg[i])
        eps1 += eps_1_lor
        eps2 += eps_2_lor

        plt.plot(eV, eps_2_lor)

    #
    # Poles
    #
    epsilon_1_UV = np.asarray([UV_Amplitude / (UV_E ** 2 - e ** 2) for e in eV])
    eps1 += epsilon_1_UV
    eps1 += eps_inf

    epsilon = eps1 + 1j * eps2

    plt.show()

    return eV, epsilon


if __name__ == "__main__":

    eV, epsilon = generate_epsilon()

    n = (epsilon ** .5).real
    k = (epsilon ** .5).imag

    #
    # Interpolate to data range
    #

    # Model range
    fit_points = 100
    fit_eV = np.linspace(0.74, 8.8, fit_points, True)


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

    # #plot n vs eV
    # plt.figure(2)
    # plt.plot(eV, epsilon.imag)
    # plt.xlabel('Photon energy (eV)')
    # plt.ylabel('eps2')


    #plot n,k vs μm
    plt.figure(3)
    plt.plot(wl_um, n_interp, label="n")
    plt.plot(wl_um, k_interp, label="k")
    plt.xlabel('Wavelength (μm)')
    plt.ylabel('n, k')
    plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

    #plot n vs eV
    plt.figure(4)
    plt.plot(fit_eV, k_interp)
    plt.xlabel('Photon energy (eV)')
    plt.ylabel('k')


    plt.show()