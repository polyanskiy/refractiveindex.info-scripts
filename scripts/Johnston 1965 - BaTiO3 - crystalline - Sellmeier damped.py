# -*- coding: utf-8 -*-

# Original data: Johnston et al. 1977, https://doi.org/10.1063/1.1660761

# Version history
# 2026-01-16 first version (Pavel Dmitriev)
#

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("TkAgg")


def generate_n(num_points=100, minL_um=0.4, maxL_um=1.0):
    #
    # Model parameters
    #
    S0 = 21.52  # um^-2
    lambda0 = 0.2799  # um
    G0 = 0.1

    S1 = 137.7  # um^-2
    lambda1 = 0.1377  # um
    G1 = 0.45

    # Simulate range
    lambdas = np.linspace(minL_um, maxL_um, num_points, True)

    #
    # 2 oscillator model
    #
    # n^2 - 1 = Seps * leps^2 / (1 - (leps^2 / lambda^2) + 2 i Geps leps/lambda)
    #             + Sgam * lgam^2 / (1 - (lgam^2 / lambda^2) + 2 i Geps lgam/lambda))
    #

    n = np.sqrt(
        1 + (S0 * lambda0**2) / (1 - (lambda0**2 / (lambdas)**2) + 2j * G0 * lambda0 / (lambdas))
        + (S1 * lambda1 ** 2) / (1 - (lambda1 ** 2 / (lambdas) ** 2) + 2j * G1 * lambda1 / (lambdas)), dtype=np.complex128
    )
    k = -np.imag(n)
    n = np.real(n)

    return n, k, lambdas


if __name__ == "__main__":
    n, k, lambdas = generate_n()

    # ============================   DATA OUTPUT   =================================
    file = open('out.txt', 'w')
    for i in range(len(n)):
        file.write('\n        {:.4e} {:.4e} {:.4e}'.format(lambdas[i], n[i], k[i]))
    file.close()

    # ===============================   PLOT   =====================================
    # plot eps1 vs eV
    plt.figure(1)
    plt.plot(lambdas, n)
    plt.plot(lambdas, k)
    plt.xlabel('Lambda, um')
    #plt.xlim([0.25, 1.05])
    #plt.ylim([1.9, 2.25])
    plt.ylabel('n')

    plt.show()
