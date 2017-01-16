# -*- coding: utf-8 -*-
#Original data: https://doi.org/10.1063/1.345115
#Author: Mikhail Polyanskiy
#Last modified: 2017-01-16

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# parameters from table II
E0  = 2.69    #eV
Δ0  = 3.10-E0 #eV
G0  = 0.017   #eV
A   = 23.4    #eV**1.5
A0x = 0.03    #eV
Γ0  = 0.030   #eV

E1  = 4.75    #eV
Δ1  = 5.05-E1 #eV
B1x = 2.31    #eV
B2x = 1.16    #eV
Γ1  = 0.37    #eV

E2  = 6.7     #eV
C   = 1.6
γ   = 0.2

ε1  = 1.2

def Epsilon_A(ħω):
    χ0 = (ħω+1j*Γ0) / E0
    χs0 = (ħω+1j*Γ0) / (E0+Δ0)
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    fχs0 = χs0**-2 * ( 2-(1+χs0)**0.5-(1-χs0)**0.5 )
    return A*E0**-1.5 * (fχ0+0.5*(E0/(E0+Δ0))**1.5*fχs0)

def Epsilon_Ax(ħω):
    y=0
    for n in range(1,100):
       y += A0x/n**3 * ( 1/((E0-G0)/n**2-ħω-1j*Γ0) + 0.5/(((E0+Δ0)-G0)/n**2-ħω-1j*Γ0) )
    return y

def Epsilon_Bx(ħω):
    y=0
    for n in range(1,100):
       y += 1/(2*n-1)**3 * (B1x/(E1/(2*n-1)**2-ħω-1j*Γ1) + B2x/((E1+Δ1)/(2*n-1)**2-ħω-1j*Γ1) )
    return y

def Epsilon_C(ħω):
    χ2 = ħω/E2
    return C/((1-χ2**2)-1j*χ2*γ)

ev_min=1.5
ev_max=5.3
npoints=1000

eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
εA = Epsilon_A(eV)
εAx = Epsilon_Ax(eV)
εBx = Epsilon_Bx(eV)
εC = Epsilon_C(eV)    
ε = εA + εAx + εBx + εC + ε1
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
plt.plot(eV, εBx.real, label="Re(εBx)")
plt.plot(eV, εC.real, label="Re(εC)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, εA.imag, label="Im(εA)")
plt.plot(eV, εAx.imag, label="Im(εAx)")
plt.plot(eV, εBx.imag, label="Im(εBx)")
plt.plot(eV, εC.imag, label="Im(εC)")
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