# Original data: Ciddor 1996, https://doi.org/10.1364/AO.35.001566

# CHANGELOG
# 2017-11-23 [Misha Polyanskiy] original version
# 2024-03-03 [Misha Polyanskiy] minor refactoring

###############################################################################

import numpy as np
import matplotlib.pyplot as plt
π = np.pi


def Z(T,p,xw): #compressibility
    t=T-273.15
    a0 = 1.58123e-6   #K·Pa^-1
    a1 = -2.9331e-8   #Pa^-1
    a2 = 1.1043e-10   #K^-1·Pa^-1
    b0 = 5.707e-6     #K·Pa^-1
    b1 = -2.051e-8    #Pa^-1
    c0 = 1.9898e-4    #K·Pa^-1
    c1 = -2.376e-6    #Pa^-1
    d  = 1.83e-11     #K^2·Pa^-2
    e  = -0.765e-8    #K^2·Pa^-2
    return 1-(p/T)*(a0+a1*t+a2*t**2+(b0+b1*t)*xw+(c0+c1*t)*xw**2) + (p/T)**2*(d+e*xw**2)


def n(λ,t,p,h,xc):
    # λ: wavelength, 0.3 to 1.69 μm 
    # t: temperature, -40 to +100 °C
    # p: pressure, 80000 to 120000 Pa
    # h: fractional humidity, 0 to 1
    # xc: CO2 concentration, 0 to 2000 ppm

    σ = 1/λ           #μm^-1
    
    T= t + 273.15     #Temperature °C -> K
    
    R = 8.314510      #gas constant, J/(mol·K)
    
    k0 = 238.0185     #μm^-2
    k1 = 5792105      #μm^-2
    k2 = 57.362       #μm^-2
    k3 = 167917       #μm^-2
 
    w0 = 295.235      #μm^-2
    w1 = 2.6422       #μm^-2
    w2 = -0.032380    #μm^-4
    w3 = 0.004028     #μm^-6
    
    A = 1.2378847e-5  #K^-2
    B = -1.9121316e-2 #K^-1
    C = 33.93711047
    D = -6.3431645e3  #K
    
    α = 1.00062
    β = 3.14e-8       #Pa^-1,
    γ = 5.6e-7        #°C^-2

    #saturation vapor pressure of water vapor in air at temperature T (Pa)
    svp = np.where(t>=0,
        np.exp(A*T**2 + B*T + C + D/T), # if t>=0
        10**(-2663.5/T+12.537))         # if t<0
    
    #enhancement factor of water vapor in air
    f = α + β*p + γ*t**2
    
    #molar fraction of water vapor in moist air
    xw = f*h*svp/p
    
    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, 450 ppm CO2
    nas = 1 + (k1/(k0-σ**2)+k3/(k2-σ**2))*1e-8
    
    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, xc ppm CO2
    naxs = 1 + (nas-1) * (1+0.534e-6*(xc-450))
    
    #refractive index of water vapor at standard conditions (20 °C, 1333 Pa)
    nws = 1 + 1.022*(w0+w1*σ**2+w2*σ**4+w3*σ**6)*1e-8
    
    Ma = 1e-3*(28.9635 + 12.011e-6*(xc-400)) #molar mass of dry air, kg/mol
    Mw = 0.018015                            #molar mass of water vapor, kg/mol
    
    Za = Z(288.15, 101325, 0)                #compressibility of dry air
    Zw = Z(293.15, 1333, 1)                  #compressibility of pure water vapor
    
    #Eq.4 with (T,P,xw) = (288.15, 101325, 0)
    ρaxs = 101325*Ma/(Za*R*288.15)           #density of standard air
    
    #Eq 4 with (T,P,xw) = (293.15, 1333, 1)
    ρws  = 1333*Mw/(Zw*R*293.15)             #density of standard water vapor
    
    # two parts of Eq.4: ρ=ρa+ρw
    ρa   = p*Ma/(Z(T,p,xw)*R*T)*(1-xw)       #density of the dry component of the moist air    
    ρw   = p*Mw/(Z(T,p,xw)*R*T)*xw           #density of the water vapor component
    
    nprop = 1 + (ρa/ρaxs)*(naxs-1) + (ρw/ρws)*(nws-1)
    
    return nprop



# output - modify code below the line to match your needs
##############################################################################
    
#use this to calculate n at particular conditions
print("n =",n(0.6328,15,101325,0,450))

#plot n vs μm
λ = np.arange(0.3, 1.691, 0.01)
n1 = n(λ,15,101325,0,450) #dry air, 15 °C, 450 ppm
n2 = n(λ,15,101325,0.5,450) #50% humidity, 15 °C, 450 ppm
n3 = n(λ,15,101325,0,370) #dry air, 15 °C, 370 ppm
n4 = n(λ,26.85,101325,0,450) #dry air, 300K, 450 ppm
plt.rc('font', family='Arial', size='14')
plt.figure(1)
plt.plot(λ, n1-1, label="dry air, 15 °C, 101325 Pa, 450 ppm CO2")
plt.plot(λ, n2-1, label="50% humidity")
plt.plot(λ, n3-1, label="370 ppm CO2")
plt.plot(λ, n4-1, label="300K (26.75 °C)")
plt.xlabel('Wavelength (μm)')
plt.ylabel('n-1')
plt.legend()

t = np.arange(-40, 100.1, 1)
n5 = n(0.6328,t,101325,0,450) #dry air, 450 ppm @ HeNe wavelength
plt.figure(2)
plt.plot(t, n5-1, label="dry air, 101325 Pa, 450 ppm CO2, 632.8 nm")
plt.xlabel('Temperature (°C)')
plt.ylabel('n-1')
plt.legend()

p = np.arange(80000, 120001, 250)
n6 = n(0.6328,15,p,0,450) #dry air, 15 °, 450 ppm @ HeNe wavelength
plt.figure(3)
plt.plot(p, n6-1, label="dry air, 15 °C, 450 ppm CO2, 632.8 nm")
plt.xlabel('Pressure (Pa)')
plt.ylabel('n-1')
plt.legend()

h = np.arange(0, 1.001, 0.01)
n7 = n(0.6328, 15, 101325, h, 450) #dry air, 15 °, 450 ppm @ HeNe wavelength
plt.figure(4)
plt.plot(h*100, n7-1, label="15 °C, 101325 Pa, 450 ppm CO2, 632.8 nm")
plt.xlabel('Humidity (%)')
plt.ylabel('n-1')
plt.legend()

xc = np.arange(0, 2001, 100)
n8 = n(0.6328, 15, 101325, 0, xc) #dry air, 15 °, 450 ppm @ HeNe wavelength
plt.figure(5)
plt.plot(xc, n8-1, label="dry air, 15 °C, 101325 Pa, 632.8 nm")
plt.xlabel('CO2 concentration (ppm)')
plt.ylabel('n-1')
plt.legend()

plt.show()
