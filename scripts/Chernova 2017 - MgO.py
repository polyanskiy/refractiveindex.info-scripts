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
#
# def tauc_lorentz_analytic(eV, Eg, A, E0, C):
#
#     #
#     # Epsilon 2
#     #
#
#     eps2 = np.zeros(eV.shape)
#     for i, e in enumerate(eV):
#         if e > Eg:
#             eps2[i] = (1 / e) * ( (A * E0 * C * (e - Eg)**2)
#                              / ((e**2 - E0**2)**2 + C**2 * e**2))
#
#     #
#     # Epsilon 1 closed form
#     #
#     eps1 = np.zeros(eV.shape)
#     for i, e in enumerate(eV):
#
#         alpha_ln = (Eg**2 - E0**2)*e*2 + Eg**2*C**2 - E0**2 * (E0**2 + 3*Eg**2)
#         alpha_atan = (e**2 - E0**2) * (E0**2 +Eg**2) + Eg**2 * C**2
#
#         alpha = np.sqrt(4 * E0**2 - C**2)
#         gamma = np.sqrt(E0**2 - C**2/2)
#
#         xi4 = (e**2 - gamma**2)**2 + alpha**2 * C**2 / 4
#
#         eps1[i] = (
#             + (1/2)*(A/np.pi)*(C/xi4)*(alpha_ln/(alpha*E0))*np.log2(
#                 (E0**2 + Eg**2 + alpha*Eg) / (E0**2 + Eg**2 - alpha*Eg)
#             )
#             - (A / (np.pi * xi4)) * (alpha_atan/E0) * (np.pi
#                                                        - np.arctan((2*Eg+alpha)/C)
#                                                        + np.arctan((-2*Eg+alpha)/C))
#             + 2 * (A*E0*C)/(np.pi*xi4) * (Eg*(e**2-gamma**2)*(np.pi + 2*np.arctan((gamma**2 - Eg**2)/(alpha*C))))
#             - 2 * ((A*E0*C)/(np.pi * xi4))*((e**2 + Eg**2)/e)*np.log2(np.abs(e-Eg)/(e+Eg))
#             + 2 * (A*E0*C)/(np.pi*xi4)*Eg*np.log2((np.abs(e-Eg)*(e+Eg))/np.sqrt((E0**2 - Eg**2)**2+Eg**2*C**2))
#         )
#
#     return eps1, eps2
#


def TaucLorentz_KK(eV, E0, A, C, Eg):
    eps_2 = np.zeros(eV.shape)

    for j, e in enumerate(eV):
        if e > Eg:
            eps_2[j] = (1. / e) * (
                    (A * E0 * C * (e - Eg) ** 2) /
                    (
                        (e ** 2 - E0 ** 2) ** 2 + C ** 2 * e ** 2
                    )
            )

    eps_1 = kk_integral_maclaurin(eV, eps_2)

    return eps_1, eps_2


def Lorentz(eV, A, FWHM, Eg):
    Br = FWHM #/ (2 * np.sqrt(np.log(2)))
    # epsilon = [A*e**2 / (Eg**2-e**2-1j*e*Br) for e in eV]
    # return np.real(epsilon), np.imag(epsilon)

    eps_2 = [A * Br**2 * Eg * e / ((Eg**2 - e**2)**2 + Br**2*e**2) for e in eV]
    eps_1 = [A * Br * Eg * (Eg**2 - e**2) / ((Eg**2 - e**2)**2 + Br**2*e**2) for e in eV]
    return np.asarray(eps_1), np.asarray(eps_2)

def generate_epsilon():

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
        eps_1_TL, eps_2_TL = TaucLorentz_KK(eV, E_TL[i], A_TL[i], C_TL[i], Eg_TL[i])
        eps1 += eps_1_TL
        eps2 += eps_2_TL

        plt.plot(eV, eps_2_TL)
    #
    # Lorentz oscillators
    #
    for i in range(len(Lorentz_Eg)):
        eps_1_lor, eps_2_lor = Lorentz(eV, Lorentz_Amplitude[i], Lorentz_FWHM[i], Lorentz_Eg[i])
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


def kk_integral_maclaurin(f, eps2):
    #
    # KK integral
    #
    df = f[1] - f[0]
    eps1 = np.zeros(f.shape)
    for i, e in enumerate(f):
        prefactor = (2 / np.pi) * 2 * df

        maclaurin_sum = 0
        if i % 2 == 0:
            js = [z for z in range(f.shape[-1])[1::2]]
        else:
            js = [z for z in range(f.shape[-1])[0::2]]

        for j in js:
            maclaurin_sum += (1. / 2.) * (
                    eps2[j] / (f[j] - e) +
                    eps2[j] / (f[j] + e)
            )
        eps1[i] = prefactor * maclaurin_sum

    return eps1


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