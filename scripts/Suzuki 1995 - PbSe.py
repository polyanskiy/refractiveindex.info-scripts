# -*- coding: utf-8 -*-
# Author: Mikhail Polyanskiy
# Last modified: 2017-04-13
# Original data: Suzuki et al. 1995, https://doi.org/10.1063/1.358926

import numpy as np
import matplotlib.pyplot as plt
π = np.pi

# model parameters
E0  = 0.31   #eV
A   = 0.8    #eV**1.5
Γ0  = 0.01   #eV

E1  = 1.6    #eV
C1  = 3.5
γ1  = 0.36

E2  = 2.73   #eV
C2  = 15.4
γ2  = 0.48

E3  = 4.10   #eV
B   = 1.2
Γ3  = 0.30   #eV

εinf= 1.5

def Epsilon_0(ħω):
    χ0 = (ħω+1j*Γ0) / E0
    fχ0 = χ0**-2 * ( 2-(1+χ0)**0.5-(1-χ0)**0.5 )
    return A*E0**-1.5 * fχ0

def Epsilon_1(ħω):
    χ1 = ħω/E1
    return C1 / (1 - χ1**2 - 1j*χ1*γ1)

def Epsilon_2(ħω):
    χ1 = ħω/E2
    return C2 / (1 - χ1**2 - 1j*χ1*γ2)

def Epsilon_3(ħω):
    χ3 = (ħω+1j*Γ3) / E3
    return -B*χ3**-2*np.log(1-χ3**2)


ev_min=1.15
ev_max=5.4
npoints=500

eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε0 = Epsilon_0(eV)
ε1 = Epsilon_1(eV)
ε2 = Epsilon_2(eV)
ε3 = Epsilon_3(eV)
ε = ε0 + ε1 + ε2 + ε3 + εinf
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
plt.plot(eV, ε1.real, label="Re(ε1)")
plt.plot(eV, ε2.real, label="Re(ε2)")
plt.plot(eV, ε3.real, label="Re(ε3)")
plt.xlabel('Photon energy (eV)')
plt.ylabel('ε')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

plt.figure(3)
plt.plot(eV, ε0.imag, label="Im(ε0)")
plt.plot(eV, ε1.imag, label="Im(ε1)")
plt.plot(eV, ε2.imag, label="Im(ε2)")
plt.plot(eV, ε3.imag, label="Im(ε3)")
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