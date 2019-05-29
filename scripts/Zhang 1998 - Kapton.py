# -*- coding: utf-8 -*-
# Author: Braden Czapla (2019)
# Last modified: 2019-05-17
# Original data: Zhang 1998, https://doi.org/10.1023/A:1022655309574

from __future__ import absolute_import, division, print_function
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
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
    np.savetxt(filename, export, fmt='%4.3f %#.4g %#.3e', delimiter=' ', header='\n'.join(header), comments='',newline='\n        ')
    return

###############################################################################

## Wavelengths to sample ##
w_um_min = 10000./6000. # [um]
w_um_max = 10000./500. # [um]
step_um = 0.005 # [um]

w_um, N_freq = w(w_um_max, w_um_min, step_um)
w_invcm = 10000./w_um
## ##

## Model Parameters ##
ResFreq = np.array([529.5, 551., 568, 590., 636.5, 673., 698.5, 736.5, 752., 763.5, 795.5, 832., 865.5, 890., 1023., 1082.5, 1123., 1177.5, 1224., 1267., 1356., 1421., 1472., 1515., 1620., 1717., 1773., 2718., 3070., 3476., 3630.]) # [cm^-1]
Strength = np.array([0.004, 0.0065, 0.002, 0.0035, 0.001, 0.0018, 0.01, 0.0106, 0.00035, 0.0008, 0.0007, 0.02, 0.0016, 0.0035, 0.001, 0.016, 0.004, 0.00065, 0.0068, 0.0025, 0.06, 0.0032, 0.0004, 0.011, 0.002, 0.026, 0.0054, 0.00018, 0.00055, 0.00017, 0.00025])
Damping = np.array([8., 10., 10., 14., 10.5, 5.5, 12.5, 6.5, 3., 8., 5., 17., 18., 8.5, 9., 35., 14., 10., 19., 17., 25., 11., 15., 10., 19., 20., 10.5, 40., 80., 40., 80.])/ResFreq
EpsInf = 3.075
## ##

## Generate and Save Data ##
eps = Lorentzian(w_invcm, ResFreq, Strength, Damping, EpsInf)
RefInd = np.sqrt(eps)

references = ' "Z. M. Zhang, G. Lefever-Button, and F. R. Powell. Infrared Refractive Index and Extinction Coefficient of Polyimide Films, <a href=\"https://doi.org/10.1023/A:1022655309574\"><i>Int. J. Thermophys.</i> <b>19</b>, 905-916 (1998)</a>"'
comments = ' "Lorentz model parameters provided in Table II of manuscript."'
SaveYML(w_um, RefInd, 'Zhang-PI (Lorentz Model).yml', references, comments)
## ##

## Plotting ##
plt.figure('Figure 5a - n')
plt.plot(w_invcm, np.real(RefInd), label='PMMA - $n$')
plt.legend(loc=1)
plt.xlim(2000,500)
plt.ylim(0.8,2.4)
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('$n$')


plt.figure('Figure 5b - k')
plt.plot(w_invcm, np.imag(RefInd), label='PMMA - $k$')

plt.legend(loc=1)
plt.xlim(2000,500)
plt.ylim(0.0,0.9)
plt.xlabel('Wavenumber (cm$^{-1}$)')
plt.ylabel('$k$')
## ##