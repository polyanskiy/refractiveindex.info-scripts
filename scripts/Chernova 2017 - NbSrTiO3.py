# -*- coding: utf-8 -*-

# Original data: Chernova et al. 2017, https://doi.org/10.1364/OME.7.003844

# Kramers-Kroning Integration https://doi.org/10.1366/0003702884430380

# Version history
# 2025-02-27 first version (Pavel Dmitriev)
#

from ruamel.yaml.scalarstring import PreservedScalarString
from ruamel.yaml import YAML
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("TkAgg")

# Model parameters
E = [3.5, 3.85, 4.16, 4.7, 4.79, 5.3, 6.32, 6.43, 9.27]
Amplitude = [0.1, 0.27, 3.44, 4.1, 0.56, 4.17, 2.43, 1.32, 4.32]
FWHM = [0.24, 0.2, 0.56, 0.94, 0.39, 1.27, 2.45, 0.55, 2.55]

UV_E = 14.8
UV_Amplitude = 228

# Drude
Lorentz_Amplitude = 98.94
Lorentz_FWHM = 0.0045
Lorentz_E = 0

eps_inf = 1

# Simulate range
num_points = 1000
eV = np.linspace(0.1, 15.0, num_points, True)
dEv = eV[1]-eV[0]

# Model range
fit_points = 100
fit_eV = np.linspace(0.74, 8.8, fit_points, True)

epsilon_1_inf = eps_inf * np.ones(eV.shape)
epsilon_1_UV = [UV_Amplitude/(UV_E**2 - e**2) for e in eV]

eps_1 = np.add(eps_inf, epsilon_1_UV)
eps_2 = np.zeros(eV.shape)

for i in range(len(E)):
    Br = FWHM[i] / (2*np.sqrt(np.log(2)))
    eps_2_osc = [
        Amplitude[i]*np.exp(
            -((e-E[i])/Br)**2
        ) - \
        Amplitude[i]*np.exp(
            -((e+E[i])/Br)**2
        )
        for e in eV
    ]

    eps_2 = np.add(eps_2, eps_2_osc)

    eps_1_osc = np.zeros(eV.shape)
    for i, e in enumerate(eV):
        prefactor = (2/np.pi)*2*dEv

        maclaurin_sum = 0
        if i % 2 == 0:
            js = [z for z in range(eV.shape[-1])[1::2]]
        else:
            js = [z for z in range(eV.shape[-1])[0::2]]

        for j in js:
            maclaurin_sum += (1./2.)*(
                    eps_2_osc[j]/(eV[j]-e) +
                    eps_2_osc[j]/(eV[j]+e)
            )
        eps_1_osc[i] = prefactor*maclaurin_sum
    eps_1 = np.add(eps_1, eps_1_osc)


epsilon = eps_1 + 1j * eps_2

#
# Lorentz-Drude
#
eps_lor = np.zeros(epsilon.shape)
for i, e in enumerate(eV):
    eps_lor[i] = 1 + Lorentz_Amplitude**2 / (Lorentz_E**2 - e**2 - 1j * e * Lorentz_FWHM)

n = (epsilon**.5).real
k = (epsilon**.5).imag

#
# Interpolate to data range
#
n_interp = np.interp(fit_eV, eV, n)
k_interp = np.interp(fit_eV, eV, k)

wl_nm = np.divide(1239.84193, fit_eV)
wl_um = np.divide(wl_nm, 1000.)

datastring = ""

datalist = []

for l, n_, k_ in zip(wl_um, n_interp, k_interp):
    datastring = datastring + "{:.4e} {:.4e} {:.4e}\n".format(l, n_, k_)
    datalist.append([l, n_, k_])

data = {
    "REFERENCES": PreservedScalarString("""E. Chernova, C. Brooks, D. Chvostova, Z. Bryknar, A. Dejneka, and M. Tyunina.
    Optical NIR-VIS-VUV constants of advanced substrates for thin-film devices
    <a href="https://doi.org/10.1364/OME.7.003844"><i>Opt. Mater. Express</i>, Vol. 7, Issue 11, pp. 3844-3862 (2017)</a>"""),

    "COMMENTS": PreservedScalarString("""Single crystal Nb-doped (0.7 wt) Strontium titanate"""),

    "DATA": {
        "type": "tabulated nk",
        "data": PreservedScalarString(datastring)
    }
}

yml = YAML()
with open('Chernova_NbSrTiO3.yml', 'w') as outfile:
    yml.dump(data, outfile)


fig, axs = plt.subplots(2)
fig.suptitle('Nb:SrTiO3')

axs[0].plot(eV, n, label="n")
axs[1].plot(eV, k, label="k")

axs[0].set_ylim([1, 3.5])
axs[1].set_ylim([0, 1.7])

axs[0].set_xlim([0, 9])
axs[1].set_xlim([0, 9])

#plt.plot(eV, eps_2, label="k", linestyle="dotted")
axs[1].set_xlabel('Energy, eV')
axs[0].set_ylabel('Optical constants')
axs[1].set_ylabel('Optical constants')

plt.legend()
plt.show()
