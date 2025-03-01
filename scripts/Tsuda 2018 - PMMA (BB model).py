# -*- coding: utf-8 -*-
# Author: Braden Czapla (2019)
# Last modified: 2019-04-30
# Original data: Tsuda et al. 2018, https://doi.org/10.1364/OE.26.006899

# Version history
# 2019-05-17 first version (Braden Czapla)
# 2025-02-28 simplify output (Misha Polyanskiy)

from __future__ import absolute_import, division, print_function
import numpy as np
from scipy.special import wofz
import matplotlib.pyplot as plt

###############################################################################

# Determine wavelengths to sample
def w(w_max, w_min, step):
    linspace_lower = (np.floor_divide(w_min, step)+1)*step
    N = np.floor_divide(w_max-w_min, step)
    linspace_upper = linspace_lower + N*step
    w = np.linspace(linspace_lower, linspace_upper, int(N)+1)
    
    if not np.isclose(w[0], w_min, atol=step/5.):
        w = np.concatenate((np.array([w_min]), w))
        
    if not np.isclose(w[-1], w_max, atol=step/5.):
        w = np.concatenate((w,np.array([w_max])))
    
    return w, len(w)

# Compute dielectric function using Brendel-Bormann (aka Gaussian or Gaussian-convoluted Drude–Lorentz) model.
# Units of w and ResFreq must match and must be directly proportional to angular frequency. All other parameters are unitless.
def Gaussian(w, ResFreq, Strength, Damping_L, Damping_G, EpsInf): # Brendel-Bormann model
    # Model Source: https://doi.org/10.1063/1.350737; https://doi.org/10.1364/AO.37.005271
   
    Permittivity = EpsInf*np.ones(len(w), dtype=complex)
    for ii in range(len(ResFreq)):
        w_bar = w/ResFreq[ii]
        a = w_bar/np.sqrt(2.)*( np.sqrt( np.sqrt( 1. + (Damping_L[ii]/w_bar)**2 ) + 1. ) + 1j*np.sqrt( np.sqrt( 1. + (Damping_L[ii]/w_bar)**2 ) - 1. ) )
        
        coeff = 1j*np.sqrt(np.pi/2.)*Strength[ii]/(2.*a*Damping_G[ii])
        Permittivity += coeff*wofz( (a-1.)/(np.sqrt(2.)*Damping_G[ii]) )
        Permittivity += coeff*wofz( (a+1.)/(np.sqrt(2.)*Damping_G[ii]) )
        
    return Permittivity

###############################################################################

## Wavelengths to sample ##
w_um_max = 10000./550. # [um]
w_um_min = 10000./4000. # [um]
step_um = 0.01 # [um]

w_um, N_freq = w(w_um_max, w_um_min, step_um)
#w_um = np.linspace(10000./4000., 10000./550., 10000)
w_invcm = 10000./w_um
## ##

## Model Parameters ##
# See Table 2
ResFreq = np.array([752.35, 807.94, 825.16, 843.27, 916.18, 967.56, 991.66, 1067.00, 1149.39, 1190.59, 1240.53, 1269.65, 1366.54, 1388.26, 1434.65, 1450.70, 1482.33, 1730.54, 2840.70, 2928.24, 2951.16, 2997.47, 3440.08]) # [cm^-1]
Strength = np.array([3.14E-03, 7.52E-04, 7.92E-05, 3.11E-03, 2.14E-03, 3.57E-03, 1.96E-03, 1.07E-03, 3.00E-02, 9.82E-03, 4.80E-03, 6.43E-03, 2.14E-03, 6.03E-04, 8.89E-04, 4.87E-03, 1.51E-03, 1.23E-02, 4.70E-05, 1.23E-03, 3.60E-04, 8.83E-04, 3.95E-05])
Damping_G = np.array([5.67, 0.87, 5.04, 0.41, 26.67, 35.26, 18.17, 21.11, 0.30, 19.38, 25.11, 0.52, 1.35, 17.00, 12.66, 4.69, 19.33, 16.23, 17.94, 4.06, 18.31, 24.82, 20.73])/ResFreq/(2.*np.sqrt(2.*np.log(2.)))
Damping_L = np.array([12.26, 16.78, 0.09, 24.88, 23.95, 0.02, 0.01, 0.02, 32.76, 13.60, 0.11, 25.53, 65.85, 0.13, 0.02, 27.42, 3.54, 7.31, 0.25, 66.46, 0.39, 25.65, 25.81])/ResFreq
EpsInf = 2.162
## ##

## Generate and Save Data ##
eps = Gaussian(w_invcm, ResFreq, Strength, Damping_L, Damping_G, EpsInf)
RefInd = np.sqrt(eps)

export = np.column_stack((w_um, np.real(RefInd), np.imag(RefInd)))
np.savetxt('out.txt', export, fmt='        %4.3f %#.6g %#.3e')

## Plotting ##
plt.figure('Figure 3a - Real(ϵ)')
plt.plot(w_invcm, np.real(eps), label='PMMA')
plt.legend(loc=1)
plt.xlim(500,4000)
plt.ylim(1,3)
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Real(ϵ)')


plt.figure('Figure 3b - Imag(ϵ)')
plt.plot(w_invcm, np.imag(eps), label='PMMA')

plt.legend(loc=1)
plt.xlim(500,4000)
plt.ylim(0,2)
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('Imag(ϵ)')
## ##