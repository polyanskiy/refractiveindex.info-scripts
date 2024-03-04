# Original data: Birch and Downs 1994, https://doi.org/10.1088/0026-1394/31/4/006

# CHANGELOG
# 2024-02-21 [aspen138] original version
# 2024-03-03 [Misha Polyanskiy] minor refactoring

###############################################################################

import numpy as np
import matplotlib.pyplot as plt

π = np.pi

def n(λ, t, p, f):
    # λ: wavelength, 0.35 to 0.65 μm
    # t: temperature, °C
    # p: pressure, Pa
    # f: partial vapour pressure, Pa

    #esT= 6.112*1e-3*np.exp(17.62*t /(243.12+t))   # the saturation vapor pressure, using the Magnus-Tetens approximation
    #pv= esT*h # pressure of water vapour

    A = 8342.54
    B = 2406147
    C = 15998
    D = 96095.43
    E = 0.601
    F = 0.00972
    G = 0.003661

    ns = 1 + 1e-8 * (A + B / (130 - λ**(-2) ) + C / (38.9 - λ**(-2)))
    ntp=1+p*(ns-1)*(1+1e-8*(E-F*t)*p)/(1+G*t)/D
    n=ntp-f*(3.73345-0.0401/ (λ)**2 )*1e-10

    return n



# output - modify code below the line to match your needs
##############################################################################

# calculate n at particular conditions
print("n =",n(0.633,19.526,102094.8,1065))

# write n(λ) data to a file
λ = np.arange(0.35, 0.65, 0.01)
n0 = n(λ,15,101325,0) # Standard air: 15 °C, 101325 Pa, dry
file = open('out.txt', 'w')
for i in range(0, len(λ)):
    file.write('\n        {:.2f} {:.12f}'.format(λ[i],n0[i]))
file.close()


#plot n vs μm
λ = np.arange(0.35, 0.65, 0.01)
n1 = n(λ,15,101325,0) #dry air, 15 °C
n2 = n(λ,26.85,101325,0) #dry air, 300K
plt.rc('font', family='Arial', size='14')
plt.figure(1)
plt.plot(λ, n1-1, label="dry air, 15 °C, 101325 Pa")
plt.plot(λ, n2-1, label="300K (26.75 °C)")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n-1')
plt.legend()

t = np.arange(-40, 100.1, 1)
n3 = n(0.6328,t,101325,0) #dry air @ HeNe wavelength
plt.figure(2)
plt.plot(t, n3-1, label="dry air, 101325 Pa, 632.8 nm")
plt.xlabel('Temperature (°C)')
plt.ylabel('n-1')
plt.legend()

p = np.arange(80000, 120001, 250)
n4 = n(0.6328,15,p,0) #dry air, 15 ° @ HeNe wavelength
plt.figure(3)
plt.plot(p, n4-1, label="dry air, 15 °, 632.8 nm")
plt.xlabel('Pressure (Pa)')
plt.ylabel('n-1')
plt.legend()

f = np.arange(0, 1500, 10)
n5 = n(0.6328, 15, 101325, f) #15 °C, HeNe wavelength
plt.figure(4)
plt.plot(f, n5-1, label="15 °, 101325 Pa, 632.8 nm")
plt.xlabel('Vapour pressure (Pa)')
plt.ylabel('n-1')
plt.legend()

plt.show()
