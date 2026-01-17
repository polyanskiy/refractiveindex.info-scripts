# -*- coding: utf-8 -*-

# Original data: Wohlecke et al. 1977, https://doi.org/10.1063/1.323822

# Version history
# 2026-01-16 first version (Pavel Dmitriev)
#

import numpy as np
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("TkAgg")


def generate_n(num_points=100, minL_um=0.3, maxL_um=1.0):
    #
    # Model parameters
    #
    E0 = 7.1 # eV  E0 = h*c/e*lambda0
    S0 = 0.96 * 10**14 # m^-2

    c = 2.99792458e8  # m/s
    h = 6.62607015e-34  # m^2 kg / s
    e = 1.602176634e-19  # C

    lambda0 = (h * c) / (e * E0)  # m


    # Simulate range
    lambdas = np.linspace(minL_um, maxL_um, num_points, True)

    #
    # DiDomenico and Wemple oscillator model
    #
    # n^2 - 1 = S0 * l0^2 / (1 - (l0^2 / lambda^2))
    #

    n = np.sqrt(
        1 + (S0 * lambda0**2) / (1 - (lambda0**2 / (lambdas * 1e-6)**2))
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
