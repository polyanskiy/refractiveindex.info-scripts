# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-10
# Original data: Rakić and Majewski 1996, https://doi.org/10.1063/1.363586

import numpy as np
import matplotlib.pyplot as plt

# model parameters
E0   = 2.993   #eV
Δ0   = 3.201-E0#eV
E1   = 3.888   #eV
Δ1   = 4.087-E1#eV

εinf = 0.002

A    = 18.0    #eV**1.5
Γ0   = 0.046   #eV
α0   = 1.11

B1   = 3.50
B1s  = 0.11
B1x  = 0.69
B2x  = 0.82
Γ1   = 0.11    #eV
α1   = 0.01

f0pr = 2.78    #eV (not eV**2!!!)
Γ0pr = 0.35    #eV
E0pr = 4.64    #eV

f2x  = 2.88    #eV (not eV**2!!!)
Γ2x  = 0.53    #eV
E2x  = 4.73    #eV

f2Σ  = 7.04    #eV (not eV**2!!!)
Γ2Σ  = 0.65    #eV
E2Σ  = 4.89    #eV

α2   = 0.27


def Epsilon_I(ħω):
    Γ = Γ0*np.exp(-α0*((ħω-E0)/Γ0)**2)
    χ0 = (ħω + 1j*Γ) / E0 + 1e-100j #1e-100j: to avoid ambiguity in (1-χ)**0.5
    χ0s = (ħω + 1j*Γ) / (E0+Δ0) + 1e-100j #1e-100j: --"--
    fχ0 =  χ0**-2  * ( 2 - (1+χ0)**0.5  - (1-χ0)**0.5 )
    fχ0s = χ0s**-2 * ( 2 - (1+χ0s)**0.5 - (1-χ0s)**0.5 )
    return A*E0**-1.5 * (fχ0+0.5*(E0/(E0+Δ0))**1.5*fχ0s)

def Epsilon_II(ħω):
    Γ = Γ1*np.exp(-α1*((ħω-E1)/Γ1)**2)
    χ1 = (ħω + 1j*Γ) / E1
    χ1s = (ħω + 1j*Γ) / (E1+Δ1)
    return -B1*χ1**-2*np.log(1-χ1**2) - B1s*χ1s**-2*np.log(1-χ1s**2)

def Epsilon_III(ħω):
    Γ = Γ1*np.exp(-α1*((ħω-E1)/Γ1)**2)
    y=0
    for n in range(1,1000):
       y += 1/(2*n-1)**3 * ( B1x/(E1-ħω-1j*Γ) + B2x/(E1+Δ1-ħω-1j*Γ) )
    return y

def Epsilon_IV(ħω):
    Γ = Γ0pr*np.exp(-α2*((ħω-E0pr)/Γ0pr)**2)
    ε0pr = f0pr**2 / (E0pr**2 - ħω**2 - 1j*ħω*Γ)
    Γ = Γ2x*np.exp(-α2*((ħω-E2x)/Γ2x)**2)
    ε2x = f2x**2 / (E2x**2 - ħω**2 - 1j*ħω*Γ)  
    Γ = Γ2Σ*np.exp(-α2*((ħω-E2Σ)/Γ2Σ)**2)
    ε2Σ = f2Σ**2 / (E2Σ**2 - ħω**2 - 1j*ħω*Γ)
    return εinf + ε0pr + ε2x + ε2Σ


ev_min=0.5
ev_max=5.6
npoints=200
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε = Epsilon_I(eV) + Epsilon_II(eV) + Epsilon_III(eV) + Epsilon_IV(eV)
n = (ε**.5).real
k = (ε**.5).imag


#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(npoints-1, -1, -1):
    file.write('\n        {:.4e} {:.4e} {:.4e}'.format(μm[i],n[i],k[i]))
file.close()
    

#===============================   PLOT   =====================================
plt.rc('font', family='Arial', size='14')

#plot ε vs eV
plt.figure(1)
plt.plot(eV, ε.real, label="ε1")
plt.plot(eV, ε.imag, label="ε2")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs eV
plt.figure(2)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs μm
plt.figure(3)
plt.plot(μm, n, label="n")
plt.plot(μm, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n, k')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)