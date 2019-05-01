# -*- coding: utf-8 -*-
# Author: Braden Czapla (2019)
# Last modified: 2019-04-29
# Original data: Kaiser et al. 1962, https://doi.org/10.1103/PhysRev.127.1950

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
    np.savetxt(filename, export, fmt='%4.2f %#.4g %#.4g', delimiter=' ', header='\n'.join(header), comments='',newline='\n        ')
    return

###############################################################################

## Wavelengths to sample ##
w_um_max = 80. # [um]
w_um_min = 10. # [um]
step_um = 0.05 # [um]

w_um, N_freq = w(w_um_max, w_um_min, step_um)
w_invcm = 10000./w_um
## ##

## Model Parameters ##
# See Table I
ResFreq = np.array([257., 328.]) # [cm^-1]
Strength = np.array([4.20, 0.40])
Damping = np.array([0.018, 0.35])
Eps_Inf = 2.045
## ##

## Generate and Save Data ##
eps = Lorentzian(w_invcm, ResFreq, Strength, Damping, Eps_Inf)
RefInd = np.sqrt(eps)

references = ' "W. Kaiser, W. G. Spitzer, R. H. Kaiser, and L. E. Howarth. Infrared Properties of CaF2, SrF2, and BaF2, <a href=\"https://doi.org/10.1103/PhysRev.127.1950\"><i>Phys. Rev.</i> <b>127</b>, 1950 (1962)</a>"'
comments = ' "Single crystal; Room temperature; Lorentz oscillator model parameters provided."'
SaveYML(w_um, RefInd, 'Kaiser-CaF2.yml', references, comments)
## ##

## Plotting ##
plt.figure('Figure 7 - n')
plt.plot(w_um, np.real(RefInd), label='CaF$_{2}$')

plt.legend(loc=1)
plt.xlim(10,80)
plt.ylim(0,14)


plt.figure('Figure 8 - k')
plt.plot(w_um, np.imag(RefInd), label='CaF$_{2}$')

plt.legend(loc=1)
plt.xlim(10,80)
plt.ylim(0,14)
## ##