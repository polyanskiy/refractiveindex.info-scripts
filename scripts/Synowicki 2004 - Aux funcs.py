# -*- coding: utf-8 -*-

# Kramers-Kroning Integration https://doi.org/10.1366/0003702884430380

# Tauc-Lorentz oscillators https://doi.org/10.1063/1.118064
#  eps1 correction erratum: https://doi.org/10.1063/1.118155

# Gaussian closed form
# https://doi.org/10.1016/J.JNONCRYSOL.2006.02.004

# Version history
# 2025-03-07 first version (Pavel Dmitriev)
#

import numpy as np
import matplotlib.pyplot as plt
import scipy.special

import matplotlib
matplotlib.use("TkAgg")

import elli

def tauc_lorentz(eV, E0, A, C, Eg):

    #
    # Epsilon 2
    #

    eps2 = np.zeros(eV.shape)

    for j, e in enumerate(eV):
        if e > Eg:
            eps2[j] = (1. / e) * (
                    (A * E0 * C * (e - Eg) ** 2) /
                    (
                        (e ** 2 - E0 ** 2) ** 2 + C ** 2 * e ** 2
                    )
            )

    #
    # Epsilon 1 closed form
    #
    eps1 = np.zeros(eV.shape)
    for i, e in enumerate(eV):

        alpha_ln = (Eg**2 - E0**2) * e**2 + Eg**2 * C**2 - E0**2 * (E0**2 + 3*Eg**2)
        alpha_atan = (e**2 - E0**2) * (E0**2 + Eg**2) + Eg**2 * C**2

        alpha = np.sqrt(4 * E0**2 - C**2)
        gamma = np.sqrt(E0**2 - C**2/2)

        xi4 = (e**2 - gamma**2)**2 + alpha**2 * C**2 / 4

        eps1[i] = (
            + ((A*C)/(np.pi*xi4))*(alpha_ln/(2*alpha*E0))*np.log(
                (E0**2 + Eg**2 + alpha*Eg) / (E0**2 + Eg**2 - alpha*Eg)
            )

            - (A / (np.pi * xi4)) * (alpha_atan/E0) * (np.pi
                                                       - np.arctan((2*Eg+alpha)/C)
                                                       + np.arctan((-2*Eg+alpha)/C))

            + 2 * (A*E0)/(np.pi*xi4*alpha) * (Eg*(e**2-gamma**2)*(np.pi + 2*np.arctan(2*(gamma**2 - Eg**2)/(alpha*C))))

            - ((A*E0*C)/(np.pi * xi4))*((e**2 + Eg**2)/e)*np.log(np.abs(e-Eg)/(e+Eg))

            + 2 * (A*E0*C)/(np.pi*xi4)*Eg*np.log((np.abs(e-Eg)*(e+Eg))/np.sqrt((E0**2 - Eg**2)**2+Eg**2*C**2))
        )

        # gamma2 = sqrt(Ei**2 - Ci**2 / 2) ** 2
        # alpha = sqrt(4 * Ei**2 - Ci**2)
        # aL = (Eg**2 - Ei**2) * E**2 + Eg**2 * Ci**2 - Ei**2 * (Ei**2 + 3 * Eg**2)
        # aA = (E**2 - Ei**2) * (Ei**2 + Eg**2) + Eg**2 * Ci**2
        # zeta4 = (E**2 - gamma2) ** 2 + alpha**2 * Ci**2 / 4
        #
        # # fmt: off
        # return (
        #     Ai*Ci*aL/2.0/np.pi/zeta4/alpha/Ei*np.log((Ei**2 + Eg**2 + alpha*Eg)/(Ei**2 + Eg**2 - alpha*Eg)) - \
        #     Ai*aA/np.pi/zeta4/Ei*(np.pi - np.arctan((2.0*Eg + alpha)/Ci) + np.arctan((alpha - 2.0*Eg)/Ci)) + \
        #     2.0*Ai*Ei*Eg/np.pi/zeta4/alpha*(E**2 - gamma2)*(np.pi + 2.0*np.arctan(2.0/alpha/Ci*(gamma2 - Eg**2))) - \
        #     Ai*Ei*Ci*(E**2 + Eg**2)/np.pi/zeta4/E*np.log(abs(E - Eg)/(E + Eg)) + \
        #     2.0*Ai*Ei*Ci*Eg/np.pi/zeta4 * \
        #     np.log(abs(E - Eg) * (E + Eg) / sqrt((Ei**2 - Eg**2)**2 + Eg**2 * Ci**2))
        # )


    return eps1, eps2


def tauc_lorentz_kk(eV, E0, A, C, Eg, kk="ML"):
    eps_2 = np.zeros(eV.shape)

    for j, e in enumerate(eV):
        if e > Eg:
            eps_2[j] = (1. / e) * (
                    (A * E0 * C * (e - Eg) ** 2) /
                    (
                        (e ** 2 - E0 ** 2) ** 2 + C ** 2 * e ** 2
                    )
            )

    if kk == "ML":
        eps_1 = kk_integral_maclaurin(eV, eps_2)
    elif kk == "FFT":
        eps_1 = kk_integral_fft(eV, eps_2)

    return eps_1, eps_2


def lorentz(eV, A, FWHM, Eg):
    Br = FWHM #/ (2 * np.sqrt(np.log(2)))
    # epsilon = [A*e**2 / (Eg**2-e**2-1j*e*Br) for e in eV]
    # return np.real(epsilon), np.imag(epsilon)

    eps_2 = [A * Br**2 * Eg * e / ((Eg**2 - e**2)**2 + Br**2*e**2) for e in eV]
    eps_1 = [A * Br * Eg * (Eg**2 - e**2) / ((Eg**2 - e**2)**2 + Br**2*e**2) for e in eV]
    return np.asarray(eps_1), np.asarray(eps_2)


def gaussian(eV, E0, Amplitude, Br):
    f = (0.5 / np.sqrt(np.log(2)))
    eps1= np.asarray([
        (2*Amplitude/np.sqrt(np.pi)) * (
            scipy.special.dawsn((e + E0)/(f*Br))
            - scipy.special.dawsn((e - E0)/(f*Br))
        )
        for e in eV
    ])

    eps2 = np.asarray([
        Amplitude * np.exp(
            -((e - E0)/(f*Br))**2
        )
        - Amplitude * np.exp(
            -((e + E0)/(f*Br))**2
        )
        for e in eV
    ])


    return eps1, eps2


def gaussian_kk(eV, E0, Amplitude, Br, kk="ML"):
    f = (0.5 / np.sqrt(np.log(2)))
    eps_2_osc = np.asarray([
        Amplitude*np.exp(
            -((e-E0)/(f*Br))**2
        ) +
        Amplitude*np.exp(
            -((e+E0)/(f*Br))**2
        )
        for e in eV
    ])

    if kk == "ML":
        eps_1_osc = kk_integral_maclaurin(eV, eps_2_osc)
    elif kk == "FFT":
        eps_1_osc = kk_integral_fft(eV, eps_2_osc)

    return eps_1_osc, np.asarray(eps_2_osc)


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


def kk_integral_fft_simple(f, eps2, num_pts = None):
    if not num_pts:
        num_pts = 2
        while num_pts < eps2.shape[0]:
            num_pts *= 2

    eps2_padded = np.zeros(num_pts)
    eps2_padded[0:eps2.shape[0]] = eps2

    # First FFT
    step1 = np.fft.ifft(eps2_padded)
    step1[eps2.shape[0]:] = np.zeros(num_pts-eps2.shape[0])

    # Second FFT
    step2 = np.fft.fft(step1)

    eps1 = 2 * np.imag(step2[0:eps2.shape[0]])

    return eps1


def kk_integral_fft(f, eps2, num_pts = None):

    df = f[1] - f[0]

    # Reverse to match paper
    r_eps2 = eps2[::-1]

    # Pad zero
    padding_f = []
    if f[0] > 0:
        padding_f = np.arange(-f[0], f[0], df)

    padding_zeros = np.zeros(padding_f.shape)

    full_r_eps2 = np.append(r_eps2, padding_zeros)

    # Pad negative
    neg_eps2 = np.multiply(eps2, -1)

    full_r_eps2 = np.append(full_r_eps2, neg_eps2)

    # Pad
    if not num_pts:
        num_pts = 2
        while num_pts < 4 * full_r_eps2.shape[0]:
            num_pts *= 2

    eps2_padded = np.zeros(num_pts)
    eps2_padded[0:full_r_eps2.shape[0]] = full_r_eps2

    # First FFT
    step1 = np.fft.ifft(eps2_padded)
    step1[full_r_eps2.shape[0]:] = np.zeros(num_pts-full_r_eps2.shape[0])

    # Second FFT
    step2 = np.fft.fft(step1)

    eps1 = -2 * np.imag(step2[0:eps2.shape[0]])

    return eps1[::-1]

