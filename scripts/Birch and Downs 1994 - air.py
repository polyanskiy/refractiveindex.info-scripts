# -*- coding: utf-8 -*-
# Author: aspen138
# Last modified: 2024-02-21
# Original data: Birch and Downs 1994, https://doi.org/10.1088/0026-1394/31/4/006

###############################################################################

import numpy as np
import matplotlib.pyplot as plt

π = np.pi

def n(λ, t, p, h, xc):
    # λ: wavelength, 0.3 to 1.69 μm
    # t: temperature, -40 to +100 °C
    # p: pressure, 80000 to 120000 Pa
    # h: fractional humidity, 0 to 1
    # xc: CO2 concentration, 0 to 2000 ppm

    esT= 6.112*1e-3*np.exp(17.62*t /(243.12+t))   # the saturation vapor pressure, using the Magnus-Tetens approximation
    pv= esT*h # pressure of water vapour


    A = 8342.54
    B = 2406147
    C = 15998
    D = 96095.43
    E = 0.601
    F = 0.00972
    G = 0.003661

    ns = 1 + 1e-8 * (A + B / (130 - λ**(-2) ) + C / (38.9 - λ**(-2)))
    ntp=1+p*(ns-1)*(1+1e-8*(E-F*t)*p)/(1+G*t)/D
    n=ntp-pv*(3.73345-0.0401/ (λ)**2 )*1e-10

    return n


# output - modify code below the line to match your needs
##############################################################################

# use this to calculate n at particular conditions
print("n =",n(0.633,20,101325,0,450))


#plot n vs μm
λ = np.arange(0.3, 1.691, 0.01)
n1 = n(λ,15,101325,0,450) #dry air, 15 °C, 450 ppm
n2 = n(λ,15,101325,0.5,450) #50% humidity, 15 °C, 450 ppm
n3 = n(λ,15,101325,0,370) #dry air, 15 °C, 370 ppm
n4 = n(λ,26.85,101325,0,450) #dry air, 300K, 450 ppm
plt.rc('font', family='Arial', size='14')
plt.figure(1)
# plt.plot(λ, n1-1, label="dry air, 15 °C, 101325 Pa, 450 ppm CO2")
plt.plot(λ, n2-1, label="50% humidity")
plt.plot(λ, n3-1, label="370 ppm CO2")
plt.plot(λ, n4-1, label="300K (26.75 °C)")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n-1')
plt.legend()

t = np.arange(-40, 100.1, 1)
n5 = [None] * len(t)
for i in range(0, len(t)):
    n5[i] = n(0.6328,t[i],101325,0,450)-1 #dry air, 450 ppm @ HeNe wavelength
plt.figure(2)
plt.plot(t, n5, label="dry air, 101325 Pa, 450 ppm CO2, 632.8 nm")
plt.xlabel('Temperature (°C)')
plt.ylabel('n-1')
plt.legend()

p = np.arange(80000, 120001, 250)
n6 = n(0.6328,15,p,0,450) #dry air, 15 °, 450 ppm @ HeNe wavelength
plt.figure(3)
plt.plot(p, n6-1, label="dry air, 15 °, 450 ppm CO2, 632.8 nm")
plt.xlabel('Pressure (Pa)')
plt.ylabel('n-1')
plt.legend()

h = np.arange(0, 1.001, 0.01)
n7 = n(0.6328, 15, 101325, h, 450) #dry air, 15 °, 450 ppm @ HeNe wavelength
plt.figure(4)
plt.plot(h*100, n7-1, label="15 °, 101325 Pa, 450 ppm CO2, 632.8 nm")
plt.xlabel('Humidity (%)')
plt.ylabel('n-1')
plt.legend()

# --- the formula doesn't depend on xc: CO2 concentration, so I comment the code
# xc = np.arange(0, 2001, 100)
# n8 = n(0.6328, 15, 101325, 0, xc) #dry air, 15 °, 450 ppm @ HeNe wavelength
# plt.figure(5)
# plt.plot(xc, n8-1, label="dry air, 15 °, 101325 Pa, 632.8 nm")
# plt.xlabel('CO2 concentration (ppm)')
# plt.ylabel('n-1')
# plt.legend()

plt.show()