# -*- coding: utf-8 -*-

# Original data: Johnston et al. 1971, https://doi.org/10.1063/1.1660761

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

    S1 = 137.7  # um^-2
    lambda1 = 0.1377  # um


    # Simulate range
    lambdas = np.linspace(minL_um, maxL_um, num_points, True)

    #
    # 2 oscillator model
    #
    # n^2 - 1 = Seps * leps^2 / (1 - (leps^2 / lambda^2)) + Sgam * lgam^2 / (1 - (lgam^2 / lambda^2))
    #

    n = np.sqrt(
        1 + (S0 * lambda0**2) / (1 - (lambda0**2 / (lambdas)**2))
        + (S1 * lambda1 ** 2) / (1 - (lambda1 ** 2 / (lambdas) ** 2))
    )


    return n, lambdas


if __name__ == "__main__":
    n, lambdas = generate_n()

    # ============================   DATA OUTPUT   =================================
    file = open('out.txt', 'w')
    for i in range(len(n) - 1):
        file.write('\n        {:.4e} {:.4e} {:.4e}'.format(lambdas[i], n[i], 0.0))
    file.close()

    # ===============================   PLOT   =====================================
    # plot eps1 vs eV
    plt.figure(1)
    plt.plot(lambdas, n)
    plt.xlabel('Lambda, um')
    #plt.xlim([0.25, 1.05])
    #plt.ylim([1.9, 2.25])
    plt.ylabel('n')

    plt.show()
