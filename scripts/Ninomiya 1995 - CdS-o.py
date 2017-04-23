# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-22
# Original data: Ninomiya and Adachi 1995, https://doi.org/10.1063/1.360355

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# model parameters
E0A = 2.482  #eV
E0B = 2.496  #eV
E0C = 2.555  #eV
G0  = 0.03   #eV
A0A = 13.0   #eV**1.5
A0B = 6.6    #eV**1.5
A0C = 6.4    #eV**1.5
A0xA= 0.045  #eV
A0xB= 0.023  #eV
A0xC= 0.022  #eV
Γ0A = 0.055  #eV
Γ0B = 0.055  #eV
Γ0C = 0.055  #eV

E1xA= 4.84   #eV
E1xB= 5.45   #eV
B1xA= 1.60   #eV
B1xB= 1.25   #eV
Γ1A = 0.42   #eV
Γ1B = 0.38   #eV

E0pr= 6.2    #eV
C   = 0.30
γ   = 0.12

εinf= 1.58

def Epsilon_0(ħω):
    χ0 = (ħω+1j*Γ0A) / E0A
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    ε = A0A*E0A**-1.5 * fχ0
    χ0 = (ħω+1j*Γ0B) / E0B
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    ε += A0B*E0B**-1.5 * fχ0
    χ0 = (ħω+1j*Γ0C) / E0C
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    ε += A0C*E0C**-1.5 * fχ0
    return ε

def Epsilon_0x(ħω):
    ε  = A0xA*(E0A-G0-ħω-1j*Γ0A)**-1
    ε += A0xB*(E0B-G0-ħω-1j*Γ0B)**-1
    ε += A0xC*(E0C-G0-ħω-1j*Γ0C)**-1
    return ε

def Epsilon_1(ħω):
    ε  = B1xA*(E1xA-ħω-1j*Γ1A)**-1
    ε += B1xB*(E1xB-ħω-1j*Γ1B)**-1
    return ε

def Epsilon_0pr(ħω):
    χ = ħω/E0pr
    return C / (1 - χ**2 - 1j*χ*γ)


ev_min=1.2
ev_max=5.7
npoints=500

eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε0   = Epsilon_0(eV)
ε0x  = Epsilon_0x(eV)
ε1   = Epsilon_1(eV)
ε0pr = Epsilon_0pr(eV)
ε = ε0 + ε0x + ε1 + ε0pr + εinf
n = (ε**.5).real
k = (ε**.5).imag
α = 4*π*k/μm*1e4 #1/cm


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

#plot intermediate data (for debugging)
plt.figure(2)
plt.plot(eV, ε0.real, label="Re(ε0)")
plt.plot(eV, ε0x.real, label="Re(ε0x)")
plt.plot(eV, ε1.real, label="Re(ε1)")
plt.plot(eV, ε0pr.real, label="Re(ε0')")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, ε0.imag, label="Im(ε0)")
plt.plot(eV, ε0x.imag, label="Im(ε0x)")
plt.plot(eV, ε1.imag, label="Im(ε1)")
plt.plot(eV, ε0pr.imag, label="Im(ε0')")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs eV
plt.figure(4)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs μm
plt.figure(5)
plt.plot(μm, n, label="n")
plt.plot(μm, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n, k')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot α vs eV
plt.figure(6)
plt.plot(eV,α)
plt.yscale('log')
plt.ylim([1e3,1e7])
plt.xlabel('Photon energy (eV)')
plt.ylabel('α (1/cm)')