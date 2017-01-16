# -*- coding: utf-8 -*-
# Original data: https://doi.org/10.1063/1.369370
# Author: Mikhail Polyanskiy
# Last modified: 2017-01-15

import numpy as np
import matplotlib.pyplot as plt

# LD model parameters - Parallel polarization (extraordinary)
ωp =   19
εinf = 0.731

f1 =   0.034
Γ1 =   2.096
ω1 =   11.418
α1 =   0.138

f2 =   0.003
Γ2 =   3.492
ω2 =   4.095
α2 =   29.728

f3 =   0.078
Γ3 =   2.442
ω3 =   10.003
α3 =   0.516

f4 =   0.131
Γ4 =   2.529
ω4 =   14.991
α4 =   1.78e-6

f5 =   0.280
Γ5 =   6.829
ω5 =   17.516
α5 =   1.78e-6 # in the original: "1.01.78e-6"

f6 =   0.855
Γ6 =   14.541
ω6 =   30.712
α6 =   1.180

f7 =   0.972
Γ7 =   20.314
ω7 =   46.004
α7 =   9.388



def LD(ω):
    ε = εinf;
    Γ = Γ1*np.exp(-α1*((ω-ω1)/Γ1)**2)
    ε -= f1*ωp**2 / ((ω**2-ω1**2)+1j*ω*Γ)
    Γ = Γ2*np.exp(-α2*((ω-ω2)/Γ2)**2)
    ε -= f2*ωp**2 / ((ω**2-ω2**2)+1j*ω*Γ)
    Γ = Γ3*np.exp(-α3*((ω-ω3)/Γ3)**2)
    ε -= f3*ωp**2 / ((ω**2-ω3**2)+1j*ω*Γ)
    Γ = Γ4*np.exp(-α4*((ω-ω4)/Γ4)**2)
    ε -= f4*ωp**2 / ((ω**2-ω4**2)+1j*ω*Γ)
    Γ = Γ5*np.exp(-α5*((ω-ω5)/Γ5)**2)
    ε -= f5*ωp**2 / ((ω**2-ω5**2)+1j*ω*Γ)
    Γ = Γ6*np.exp(-α6*((ω-ω6)/Γ6)**2)
    ε -= f6*ωp**2 / ((ω**2-ω6**2)+1j*ω*Γ)
    Γ = Γ7*np.exp(-α7*((ω-ω7)/Γ7)**2)
    ε -= f7*ωp**2 / ((ω**2-ω7**2)+1j*ω*Γ)
    return ε


ev_min=2.1
ev_max=40
npoints=1000
eV = np.linspace(ev_min, ev_max, npoints)
μm = 4.13566733e-1*2.99792458/eV
ε = LD(eV)
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