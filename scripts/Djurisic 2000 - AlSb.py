# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-02
# Original data: Djurišić et al. 2000, https://doi.org/10.1007/s003390050006

import numpy as np
import matplotlib.pyplot as plt

# model parameters from table I
E0   = 2.27    #eV
Δ0   = 2.93-E0 #eV
E1   = 2.82    #eV
Δ1   = 3.22-E1 #eV
A    = 0.009   #eV**1.5

εinf = 0.068
Γ0   = 0.006   #eV
α0   = 6.579

B1   = 6.025
B1s  = 0.328
Γ1   = 0.146   #eV
α1   = 0.062

F2   = 5.261   #eV
Γ2   = 1.209   #eV
α2   = 1.486
E2   = 3.175   #eV

F3   = 1.540   #eV
Γ3   = 0.426   #eV
α3   = 1.318
E3   = 5.254   #eV

F4   = 5.338   #eV
Γ4   = 0.491   #eV
α4   = 0.429
E4   = 4.016   #eV


def Epsilon_I(ħω):
    Γ0ω = Γ0*np.exp(-α0*((ħω-E0)/Γ0)**2)
    χ0 = (ħω + 1j*Γ0ω) / E0
    χ0s = (ħω + 1j*Γ0ω) / (E0+Δ0)
    fχ0 = χ0**-2 * ( 2  -(1+χ0)**0.5 - (1-χ0)**0.5 )
    fχ0s = χ0s**-2 * ( 2 - (1+χ0s)**0.5 - (1-χ0s)**0.5 )
    return A*E0**-1.5 * (fχ0+0.5*(E0/(E0+Δ0))**1.5*fχ0s)

def Epsilon_II(ħω):
    Γ1ω = Γ1*np.exp(-α1*((ħω-E1)/Γ1)**2)
    χ1 = (ħω + 1j*Γ1ω) / E1
    χ1s = (ħω + 1j*Γ1ω) / (E1+Δ1)
    return -B1*χ1**-2*np.log(1-χ1**2) - B1s*χ1s**-2*np.log(1-χ1s**2)

def Epsilon_III(ħω):
    Γ2ω = Γ2*np.exp(-α2*((ħω-E2)/Γ2)**2)
    ε2 = F2**2 / (E2**2 - ħω**2 - 1j*ħω*Γ2ω)    
    Γ3ω = Γ3*np.exp(-α3*((ħω-E3)/Γ3)**2)
    ε3 = F3**2 / (E3**2 - ħω**2 - 1j*ħω*Γ3ω)
    Γ4ω = Γ4*np.exp(-α4*((ħω-E4)/Γ4)**2)
    ε4 = F4**2 / (E4**2 - ħω**2 - 1j*ħω*Γ4ω)
    return ε2+ε3+ε4


ev_min=0.1
ev_max=6
npoints=200
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε = εinf + Epsilon_I(eV) + Epsilon_II(eV) + Epsilon_III(eV)
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