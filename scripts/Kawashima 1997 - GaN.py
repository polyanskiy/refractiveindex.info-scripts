# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-23
# Original data: Kawashima et al. 1997, https://doi.org/10.1063/1.365671

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# parameters from table II
E0  = 3.38    #eV
G0  = 0.02    #eV
A0  = 27      #eV**1.5
A0x = 0.055   #eV**-1
Γ0  = 0.06    #eV

EminusG1A=6.8 #eV
EminusG1B=7.9 #eV
EminusG1C=9.0 #eV
B1xA = 6.2    #eV
B1xB = 0.6    #eV
B1xC = 3.0    #eV
Γ1A  = 0.78   #eV
Γ1B  = 0.35   #eV
Γ1C  = 1.0    #eV

ε1   = 2.20

def Epsilon_A(ħω):
    χ0 = (ħω+1j*Γ0) / E0
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    return A0*E0**-1.5 * fχ0

def Epsilon_Ax(ħω):
    y=0
    for n in range(1,1000):
       y += A0x/n**3 / (E0-G0/n**2-ħω-1j*Γ0) 
    return y

def Epsilon_B(ħω):
    ε  = B1xA / (EminusG1A-ħω-1j*Γ1A)
    ε += B1xB / (EminusG1B-ħω-1j*Γ1B)
    ε += B1xC / (EminusG1C-ħω-1j*Γ1C)
    return ε

ev_min=1.25
ev_max=10
npoints=500

eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
εA = Epsilon_A(eV)
εAx = Epsilon_Ax(eV)
εB = Epsilon_B(eV)    
ε = εA + εAx + εB + ε1
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
plt.plot(eV, εA.real, label="Re(εA)")
plt.plot(eV, εAx.real, label="Re(εAx)")
plt.plot(eV, εB.real, label="Re(εB)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, εA.imag, label="Im(εA)")
plt.plot(eV, εAx.imag, label="Im(εAx)")
plt.plot(eV, εB.imag, label="Im(εB)")
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