# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-09
# Original data: Ozaki and Adachi 1995, https://doi.org/10.1063/1.359966

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# model parameters from table I
E0   = 1.42    #eV
Δ0   = 1.75-E0 #eV
A    = 7.0     #eV**1.5
Γ0   = 0.03    #eV

E1   = 2.91    #eV
Δ1   = 3.14-E1 #eV
B1   = 3.50
B2   = 1.75
B1x  = 1.2     #eV
B2x  = 0.6     #eV
Γ1   = 0.12    #eV

E0pr = 4.45    #eV
C0pr = 0.80
Γ0pr = 0.67

E21  = 4.77    #eV
C21  = 1.35
Γ21  = 0.62

E22  = 5.00    #eV
C22  = 0.35
Γ22  = 0.4     #eV

E1pr = 6.6     #eV
C1pr = 0.7
Γ1pr = 0.6     #eV


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
    y=0
    for n in range(1,1000):
       y += 1/(2*n-1)**3 * ( B1x/(E1-ħω-1j*Γ1) + B2x/(E1+Δ1-ħω-1j*Γ1) )
    return y

def Epsilon_C0pr(ħω):
    return C0pr*E0pr**2/(E0pr**2-ħω**2-1j*ħω*Γ0pr)

def Epsilon_C1pr(ħω):
    return C1pr*E1pr**2/(E1pr**2-ħω**2-1j*ħω*Γ1pr)

def Epsilon_C21(ħω):
    return C21*E21**2/(E21**2-ħω**2-1j*ħω*Γ21)

def Epsilon_C22(ħω):
    return C22*E22**2/(E22**2-ħω**2-1j*ħω*Γ22)
     

ev_min = 1.2
ev_max = 5.6
npoints = 500
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
εA    = Epsilon_A(eV)
εB    = Epsilon_B(eV)
εBx   = Epsilon_Bx(eV)
εC0pr = Epsilon_C0pr(eV)
εC1pr = Epsilon_C1pr(eV)
εC21  = Epsilon_C21(eV)
εC22  = Epsilon_C22(eV)
ε = εA + εB + εBx + εC0pr + εC1pr + εC21 + εC22
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
plt.plot(eV, εC0pr.real, label="Re(εC0pr)")
plt.plot(eV, εC1pr.real, label="Re(εC1pr)")
plt.plot(eV, εC21.real, label="Re(εC21)")
plt.plot(eV, εC22.real, label="Re(εC22)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, εA.imag, label="Im(εA)")
plt.plot(eV, εB.imag, label="Im(εB)")
plt.plot(eV, εBx.imag, label="Im(εBx)")
plt.plot(eV, εC0pr.imag, label="Re(εC0pr)")
plt.plot(eV, εC1pr.imag, label="Re(εC1pr)")
plt.plot(eV, εC21.imag, label="Re(εC21)")
plt.plot(eV, εC22.imag, label="Re(εC22)")
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