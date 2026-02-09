# -*- coding: utf-8 -*-
# Fresnel reflection calculator
# (see https://en.wikipedia.org/wiki/Fresnel_equations for equations)
# Author: Mikhail Polyanskiy
# Last modified: 2017-11-24

import numpy as np

########################## input parameters ###################################
n1 = 1.0                 #complex ior of first medium (1 for vacuum)
n2 = 3.4209 + 7.7479e-7j #complex ior of second medium
θi = 0                   #incidence angle (degrees)
###############################################################################

θi = np.deg2rad(θi)              # incidence angle (radians)
θt = np.arcsin(n1/n2*np.sin(θi)) # refraction angle (radians)

rs = (n1*np.cos(θi)-n2*np.cos(θt)) / (n1*np.cos(θi)+n2*np.cos(θt))
rp = (n2*np.cos(θi)-n1*np.cos(θt)) / (n1*np.cos(θt)+n2*np.cos(θi))

Rs = np.abs(rs)**2
Rp = np.abs(rp)**2

Φs = np.rad2deg(np.angle(rs))
Φp = np.rad2deg(np.angle(rp))

print('Rs = {:f}\nRp = {:f}\nΦs = {:f}\nΦp = {:f}'.format(Rs, Rp, Φs, Φp))
