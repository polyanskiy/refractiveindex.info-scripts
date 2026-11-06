# Original data: Galander et al. 2026, https://doi.org/10.1016/j.matlit.2026.100067

# CHANGELOG
# 2026-11-05 [Misha Polyanskiy] original version


import numpy as np
import matplotlib.pyplot as plt

# model parameters from tables 1 and 4
E0   = 0.72    #eV
Δ0   = 0.74    #eV
E1   = 2.17    #eV
Δ1   = 0.46    #eV
A    = 2.1957  #eV**1.5

Γ0   = 0.002   #eV
α0   = 0.0101

B1   = 12.7828
B1s  = 0.5606
Γ1   = 0.6042  #eV
α1   = 0.5594

F2   = 0.2265  #eV
Γ2   = 0.2695  #eV
α2   = 0.997
E2   = 0.4348  #eV

F3   = 0.142   #eV
Γ3   = 0.1222  #eV
α3   = 2.9659
E3   = 0.704   #eV


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
    return ε2+ε3

μm_min=1
μm_max=3.1
npoints=1051
μm = np.linspace(μm_min, μm_max, npoints)
eV = 4.13566733e-1*2.99792458/μm

ε = Epsilon_I(eV) + Epsilon_II(eV) + Epsilon_III(eV)
n = (ε**.5).real
k = (ε**.5).imag


#============================   DATA OUTPUT   =================================
file = open('out.txt', 'w')
for i in range(0, npoints):
    file.write('\n        {:.3f} {:.4f} {:.5f}'.format(μm[i],n[i],k[i]))
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

#plot n vs μm
plt.figure(3)
plt.plot(μm, n, label="n")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)

#plot k vs μm
plt.figure(4)
plt.plot(μm, k, label="k")
plt.xlabel('Wavelength (μm)')
plt.ylabel('k')
plt.legend(bbox_to_anchor=(0,1.02,1,0),loc=3,ncol=2,borderaxespad=0)
