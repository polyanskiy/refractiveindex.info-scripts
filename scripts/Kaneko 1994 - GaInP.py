# Author: Pedro Victor Sansoldo, University of Adelaide, Australia
# Data: https://doi.org/10.1063/1.357699

# Plot experimental points (two colors) + fitted curve

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Experimental data ---

# d = 851 nm
wvl_851 = np.array([645,670,695,725,765,810,865,920,995,1085,1195,1345])
n_851   = np.array([3.5800,3.5200,3.4550,3.4000,3.3675,3.3350,3.2950,3.2550,3.2275,3.2000,3.1675,3.1550])

# d = 580 nm
wvl_580 = np.array([650,680,725,785,850,940,1065,1225])
n_580   = np.array([3.5925,3.4825,3.4225,3.3550,3.2925,3.2475,3.2100,3.1650])

# Combine for fitting
wvl_all_nm = np.concatenate([wvl_851, wvl_580])
n_all = np.concatenate([n_851, n_580])

wvl_all_um = wvl_all_nm / 1000.0

# --- 2-term Sellmeier model ---
def sellmeier_1(lam_um, B1, C1, B2, C2):
    return np.sqrt(1 + (B1 * lam_um**2) / (lam_um**2 - C1) + (B2 * lam_um**2) / (lam_um**2 - C2))

# Fit
p0 = [8.0, 0.12, 1.0, 0.01]
c1_upper = (wvl_all_um.min()**2) * 0.95
c2_upper = (wvl_all_um.min()**2) * 0.95
bounds = ([0.0, 1e-8, 0.0, 1e-8], [200.0, c1_upper, 200.0, c2_upper])

params, _ = curve_fit(sellmeier_1, wvl_all_um, n_all, p0=p0, bounds=bounds, maxfev=500000)
B1_final, C1_final, B2_final, C2_final = params

# --- Plot ---
lam_plot_um = np.linspace(wvl_all_um.min(), wvl_all_um.max(), 600)

plt.figure()

# Experimental points (two different colors as requested)
plt.scatter(wvl_851/1000, n_851, color='blue', label='Experimental d = 851 nm')
plt.scatter(wvl_580/1000, n_580, color='red', label='Experimental d = 580 nm')

# Fitted curve
plt.plot(lam_plot_um, sellmeier_1(lam_plot_um, B1_final, C1_final, B2_final, C2_final), 
         color='black', label='Sellmeier fit (combined)')

plt.xlabel("Wavelength (µm)")
plt.ylabel("Refractive Index n")
plt.title("GaInP Sellmeier Fit (Combined Thickness Data)")
plt.legend()
plt.show()

results = {
    "Final_B1": float(B1_final),
    "Final_C1_um2": float(C1_final),
    "Final_B2": float(B2_final),
    "Final_C3_um2": float(C2_final),
    "Resonance_lambda0_1_um": float(np.sqrt(C1_final)),
    "Resonance_lambda0_2_um": float(np.sqrt(C2_final))
}

print(results)