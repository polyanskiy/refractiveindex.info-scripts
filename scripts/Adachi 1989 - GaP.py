# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-05-25
# Original data: Adachi 1989, https://doi.org/10.1063/1.343580

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# model parameters
E0   = 2.74    #eV
Δ0   = 2.84-E0 #eV
E1   = 3.70    #eV
Δ1   = 0       #eV
E2   = 5.0     #eV
Eg   = 2.26    #eV
A    = 13.76   #eV**1.5
B1   = 6.35
B11  = 9.49    #eV**-0.5
Γ    = 0.06    #eV
C    = 2.08
γ    = 0.132
D    = 4.6
εinf = 0.0

def H(x): #Heviside function
    return 0.5 * (np.sign(x) + 1)

def Epsilon_A(ħω): #E0
    χ0 = ħω/E0
    χso = ħω / (E0+Δ0)
    H0 = H(1-χ0)
    Hso = H(1-χso)
    fχ0 = χ0**-2 * ( 2 -(1+χ0)**0.5 - ((1-χ0)*H0)**0.5 )
    fχso = χso**-2 * ( 2 - (1+χso)**0.5 - ((1-χso)*Hso)**0.5 )
    H0 = H(χ0-1)
    Hso = H(χso-1)
    ε2 = A/(ħω)**2 * ( ((ħω-E0)*H0)**0.5 + 0.5*((ħω-E0-Δ0)*Hso)**0.5)
    ε1 = A*E0**-1.5 * (fχ0+0.5*(E0/(E0+Δ0))**1.5*fχso)
    return ε1 + 1j*ε2
    
def Epsilon_B(ħω): #E1
    # ignoring E1+Δ1 contribution - no data on B2 & B21 in the paper
    # result seems to reproduce graphical data from the paper
    χ1 = ħω/E1
    H1 = H(1-χ1)
    ε2 = π*χ1**-2*(B1-B11*((E1-ħω)*H1)**0.5)
    ε2 *= H(ε2) #undocumented trick: ignore negative ε2
    χ1 = (ħω+1j*Γ)/E1
    ε1 = -B1*χ1**-2*np.log(1-χ1**2)
    return ε1.real + 1j*ε2.real

def Epsilon_C(ħω): #E2
    χ2 = ħω/E2
    ε2 = C*χ2*γ / ((1-χ2**2)**2+(χ2*γ)**2)
    ε2 *= H(ħω-3.4) #undocumented trick: ignore low-energy shoulder
                    #(see E2 curve in Fig. 2a)
    ε1 = C*(1-χ2**2) / ((1-χ2**2)**2+(χ2*γ)**2)
    return ε1 + 1j*ε2

def Epsilon_D(ħω): #Eg
    # ignoring ħωq - no data in the paper
    # result seems to reproduce graphical data from the paper
    Ech = E1
    χg = Eg/ħω
    χch = ħω/Ech
    Hg = H(1-χg)
    Hch = H(1-χch)
    ε2 = D/ħω**2 * (ħω-Eg)**2 * Hg * Hch
    return 1j*ε2
    

ev_min=0.1
ev_max=6
npoints=200
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
εA  = Epsilon_A(eV)
εB  = Epsilon_B(eV)
εC  = Epsilon_C(eV)
εD  = Epsilon_D(eV)    
ε = εA + εB + εC + εD + εinf
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

#plot ε1 vs eV
plt.figure(1)
plt.plot(eV, ε.real, label="ε1")
plt.plot(eV, εA.real, label="Re(εA)")
plt.plot(eV, εB.real, label="Re(εB)")
plt.plot(eV, εC.real, label="Re(εC)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε1')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot ε2 vs eV
plt.figure(2)
plt.plot(eV, ε.imag, label="ε2")
plt.plot(eV, εA.imag, label="Im(εA)")
plt.plot(eV, εB.imag, label="Im(εB)")
plt.plot(eV, εC.imag, label="Im(εC)")
plt.plot(eV, εD.imag, label="Im(εD)")
plt.yscale('log')
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε2')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
plt.ylim([1e-2,1e2])

#plot n,k vs eV
plt.figure(3)
plt.plot(eV, n, label="n")
plt.plot(eV, k, label="k")
plt.xlabel('Photon energy (eV)')
plt.ylabel('n, k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot n,k vs μm
plt.figure(4)
plt.plot(μm, n, label="n")
plt.plot(μm, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n, k')
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot α vs eV
plt.figure(7)
plt.plot(eV,α)
plt.yscale('log')
plt.ylim([1e3,1e7])
plt.xlabel('Photon energy (eV)')
plt.ylabel('α (1/cm)')