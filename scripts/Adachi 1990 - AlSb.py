# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-02
# Original data: Adachi 1990, https://doi.org/10.1063/1.345115

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# parameters from table I
E0   = 2.27    #eV
Δ0   = 2.9-E0  #eV
A    = 29.0    #eV**1.5
Γ0   = 0.015   #eV

E1   = 2.84    #eV
Δ1   = 3.23-E1 #eV
B1   = 1.45
B2   = 1.17
B1x  = 1.98    #eV
B2x  = 1.60    #eV
Γ1   = 0.14    #eV

E0pr = 3.70    #eV
C0   = 1.83
γ0   = 0.18

E1pr = 5.25    #eV
C1   = 0.06
γ1   = 0.06

E2   = 4.05    #eV
C2   = 1.23
γ2   = 0.09

Eg   = 1.61    #eV
D    = 1.00
Γg   = 0.040   #eV


def Epsilon_A(ħω):
    χ0 = (ħω + 1j*Γ0) / E0
    χs0 = (ħω + 1j*Γ0) / (E0+Δ0)
    fχ0 = χ0**-2 * ( 2  -(1+χ0)**0.5 - (1-χ0)**0.5 )
    fχs0 = χs0**-2 * ( 2 - (1+χs0)**0.5 - (1-χs0)**0.5 )
    return A*E0**-1.5 * (fχ0+0.5*(E0/(E0+Δ0))**1.5*fχs0)

def Epsilon_B(ħω):
    χ1d = (ħω + 1j*Γ1) / E1
    χ1sd = (ħω + 1j*Γ1) / (E1+Δ1)
    return -B1*χ1d**-2*np.log(1-χ1d**2) - B2*χ1sd**-2*np.log(1-χ1sd**2)

def Epsilon_Bx(ħω):
    return B1x/(E1-ħω-1j*Γ1) +  B2x/(E1+Δ1-ħω-1j*Γ1)

def Epsilon_C0(ħω):
    χ2 = ħω/E0pr
    return C0  /((1-χ2**2)-1j*χ2*γ0)

def Epsilon_C1(ħω):
    χ2 = ħω/E1pr
    return C1 / ((1-χ2**2)-1j*χ2*γ1)
    
def Epsilon_C2(ħω):
    χ2 = ħω/E2
    return C2 / ((1-χ2**2)-1j*χ2*γ2)

def Epsilon_D(ħω):
    Ec = E1
    return 2*D/π * ( -Eg**2/(ħω+1j*Γg)**2*np.log(Ec/Eg)
    + 0.5*(1+Eg/(ħω+1j*Γg))**2*np.log((ħω+1j*Γg+Ec)/(ħω+1j*Γg+Eg))
    + 0.5*(1-Eg/(ħω+1j*Γg))**2*np.log((ħω+1j*Γg-Ec)/(ħω+1j*Γg-Eg)) )
    
    

ev_min=0.1
ev_max=6
npoints=200
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
εA  = Epsilon_A(eV)
εB  = Epsilon_B(eV)
εBx = Epsilon_Bx(eV)
εC0 = Epsilon_C0(eV)
εC1 = Epsilon_C1(eV)
εC2 = Epsilon_C2(eV)
εD  = Epsilon_D(eV)    
ε = εA + εB + εBx + εC0 + εC1 + εC2 + εD
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
plt.plot(eV, εB.real, label="Re(εB)")
plt.plot(eV, εBx.real, label="Re(εBx)")
plt.plot(eV, εC0.real, label="Re(εC0)")
plt.plot(eV, εC1.real, label="Re(εC1)")
plt.plot(eV, εC2.real, label="Re(εC2)")
plt.plot(eV, 10*εD.real, ls='--', label="10*Re(εD)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, εA.imag, label="Im(εA)")
plt.plot(eV, εB.imag, label="Im(εB)")
plt.plot(eV, εBx.imag, label="Im(εBx)")
plt.plot(eV, εC0.imag, label="Im(εC0)")
plt.plot(eV, εC1.imag, label="Im(εC1)")
plt.plot(eV, εC2.imag, label="Im(εC2)")
plt.plot(eV, 10*εD.imag, ls='--', label="10*Im(εD)")
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
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot α vs eV
plt.figure(6)
plt.plot(eV,α)
plt.yscale('log')
plt.ylim([1e3,1e7])
plt.xlabel('Photon energy (eV)')
plt.ylabel('α (1/cm)')