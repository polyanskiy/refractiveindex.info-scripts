# -*- coding: utf-8 -*-
# Author: Andrea Gerini
# Last modified: 2025-07-31
# Original data: Gehrsitz et al. 2000, https://doi.org/10.1063/1.373462

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, e as q

def coth(x):
    return np.cosh(x) / np.sinh(x)


def nAlGaAs_Gehrsitz(x, λ, Temperature=20):
    """Al(x)Ga(1-x)As refractive index model according to [1].

    Parameters
    ----------
    x : float
        Al molar fraction (0 <= x <= 1)
    λ : float
        Wavelength in [m] (λ > 0)
    Temperature : float, optional
                  Temperature in [°C], default to 20°C

    Returns
    -------
    n + ik : complex
        Complex refractive index

    References
    ----------
    .. [1] S. Gehrsitz, F. K. Reinhart, C. Gourgon, N. Herres, A. Vonlanthen,
       and H. Sigg, "The refractive index of Al$_{x}$Ga$_{1-x}$As below the
       band gap: Accurate determination and empirical modeling," Journal of
       Applied Physics 87, 7825–7837 (2000).

    """
    T = Temperature + 273.15
    µm = 1.239856          # [eV µm], energy conversion to µm⁻¹
    kB = 0.0861708 * 1e-3  # [eV/K], Boltzmann constant
    V_T = 2 * kB * T       # [eV]

    E = h * c / q / λ / µm  # [µm⁻¹]

    # Energy bandgap (in Γ) for GaAs [eq. 11]
    EΓ0 = 1.5192  # [eV]
    EDeb = 15.9 * 1e-3  # [eV]
    ETO = 33.6 * 1e-3  # [eV]
    S = 1.8
    STO = 1.1
    EΓ_GaAs = EΓ0 + S * EDeb * (1 - coth(EDeb / V_T)) + STO * ETO * (1 - coth(ETO / V_T))  # [eV]

    # Wavelength-independent contribution to Sellmeier equation [Tab. II, GaAs FIT 2]
    A0 = 5.9613 + 7.7178e-4 * T - 0.953e-6 * T**2
    A = A0 - 16.159 * x + 43.511 * x**2 - 71.317 * x**3 + 57.535 * x**4 - 17.451 * x**5
    E1_2 = 4.7171 - 3.237e-4 * T - 1.358e-6 * T**2  # [µm⁻²]
    E2 = 0.724e-3  # [µm⁻²], GaAs E₂²
    C2 = 1.55e-3  # [µm⁻²], GaAs C₂
    E3 = 1.331e-3  # [µm⁻²], AlAs E₂²
    C3 = 2.61e-3  # [µm⁻²], AlAs C₂

    # Wavelength-dependent contributions to Sellmeier equation [Tab. IV]
    E0 = (EΓ_GaAs / µm + 1.1308 * x + 0.1436 * x**2)**2  # [µm⁻²]
    C0 = 1 / (50.535 - 150.7 * x - 62.209 * x**2 + 797.16 * x**3 - 1125 * x**4 + 503.79 * x**5)  # [µm⁻²]
    E1 = E1_2 + 11.006 * x - 3.08 * x**2  # [µm⁻²]
    C1 = 21.5647 + 113.74 * x - 122.5 * x**2 + 108.401 * x**3 - 47.318 * x**4  # [µm⁻²]

    # Sellmeier equation [eq. 12, 13]
    R = (1 - x) * C2 / (E2 - E**2) + x * C3 / (E3 - E**2)
    n2 = A + C0 / (E0 - E**2) + C1 / (E1 - E**2) + R
    return np.sqrt(n2, dtype=complex)

     
λ_min = 440e-9   # [m]
λ_max = 3100e-9  # [m]
npoints = 500
temperature = 23 # [°C]
Al_fraction = [0, 0.176, 0.334, 0.410, 0.427, 0.615, 0.753, 0.865, 1]

wavelengths = np.linspace(λ_min, λ_max, npoints)
eV = h * c / q / wavelengths  # [eV]


plt.rc('font', family='Arial', size='14')

fig_eV, ax_eV = plt.subplots()
fig_µm, ax_µm = plt.subplots()

for x in Al_fraction:
    n = nAlGaAs_Gehrsitz(x, wavelengths, Temperature=temperature)

    ax_eV.plot(eV, n.real, label=f'{x = :.3f}')
    ax_eV.set_xlabel('Photon energy (eV)')
    ax_eV.set_ylabel('Refractive index n')
    ax_eV.legend()
    ax_eV.set_ylim(2.85, 3.67)

    ax_µm.plot(wavelengths * 1e6, n.real, label=f'{x = :.3f}')
    ax_µm.set_xlabel('Wavelength (µm)')
    ax_µm.set_ylabel('Refractive index n')
    ax_µm.legend()
    ax_µm.set_ylim(2.85, 3.67)
    
plt.show()
