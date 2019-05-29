# -*- coding: utf-8 -*-
# Author: Braden Czapla (2019)
# Last modified: 2019-04-30
# Original data: Tsuda et al. 2018, https://doi.org/10.1364/OE.26.006899

from __future__ import absolute_import, division, print_function
import numpy as np
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

# Compute dielectric function using Lorentzian model.
# Units of w and ResFreq must match and must be directly proportional to angular frequency. All other parameters are unitless.
def Lorentzian(w, ResFreq, Strength, Damping, Eps_Inf):
    Permittivity = Eps_Inf*np.ones(len(w), dtype=np.complex)
    for ii in range(len(ResFreq)):
        Permittivity += Strength[ii]/( 1. - (w/ResFreq[ii])**2 - 1j*Damping[ii]*(w/ResFreq[ii]) )   
    return Permittivity

# Save w, n, k to YML file
def SaveYML(w_um, RefInd, filename, references='', comments=''):
    
    header = np.empty(9, dtype=object)
    header[0] = '# this file is part of refractiveindex.info database'
    header[1] = '# refractiveindex.info database is in the public domain'
    header[2] = '# copyright and related rights waived via CC0 1.0'
    header[3] = ''
    header[4] = 'REFERENCES:' + references
    header[5] = 'COMMENTS:' + comments
    header[6] = 'DATA:'
    header[7] = '  - type: tabulated nk'
    header[8] = '    data: |'
    
    export = np.column_stack((w_um, np.real(RefInd), np.imag(RefInd)))
    np.savetxt(filename, export, fmt='%4.2f %#.4g %#.3e', delimiter=' ', header='\n'.join(header), comments='',newline='\n        ')
    return

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
# See Table 1
ResFreq  = np.array([752.25, 808.09, 825.19, 843.16, 913.82, 965.31, 989.60, 1066.27, 1149.37, 1190.32, 1241.23, 1269.59, 1361.50, 1387.61, 1434.77, 1450.59, 1481.89, 1730.18, 2840.98, 2920.93, 2950.55, 2997.71, 3440.07]) # [cm^-1]
Strength = np.array([3.18E-03, 6.94E-04, 1.13E-04, 2.86E-03, 1.68E-03, 3.94E-03, 2.79E-03, 1.10E-03, 2.92E-02, 1.04E-02, 6.64E-03, 5.49E-03, 1.09E-03, 1.07E-03, 1.34E-03, 4.11E-03, 2.12E-03, 1.56E-02, 6.66E-05, 8.42E-04, 6.60E-04, 9.53E-04, 4.15E-05])
Damping  = np.array([13.65, 15.51, 4.21, 23.09, 32.50, 26.70, 14.68, 13.56, 31.12, 22.12, 21.38, 24.66, 41.97, 15.80, 10.56, 25.15, 19.17, 9.40, 15.32, 60.94, 18.80, 36.68, 33.89])/ResFreq
EpsInf = 2.162
## ##

## Generate and Save Data ##
eps = Lorentzian(w_invcm, ResFreq, Strength, Damping, EpsInf)
RefInd = np.sqrt(eps)

references = ' "S. Tsuda, S. Yamaguchi, Y. Kanamori, and H. Yugami. Spectral and angular shaping of infrared radiation in a polymer resonator with molecular vibrational modes, <a href=\"https://doi.org/10.1364/OE.26.006899\"><i>Opt. Express</i> <b>26</b>, 6899-6915 (2018)</a>"'
comments = ' "MicroChem PMMA resist with a molecular weight of 950,000; Baked at 100°C for 10 min on a hot plate; Lorentz model parameters provided in Table 1 of manuscript."'
SaveYML(w_um, RefInd, 'Tsuda-PMMA (Lorentz Model).yml', references, comments)
## ##

## Plotting ## Note: This is not an exact match because Fig. 3 plots the results of the Brendel-Bormann model. The Lorentz model is very similar, but with more extreme peaks.
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