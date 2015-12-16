### BEGIN EXAMPLE ###
#!/usr/bin/python
# example of how to use ldbb.py - text version

from ldbb import LDBB
from numpy import *

# wavelength range of interest (in meters) and number of points
lambda0 = logspace(log10(0.2066e-6),log10(12.40e-6),200)
#epsD = LDBB('Au','D',lambda0)
#epsLD = LDBB('Au','LD',lambda0)
epsBB = LDBB('Ti','BB',lambda0)


# LDBB returns permittivity, but we can easily convert it to refractive index
for l,e in zip(lambda0,epsBB):
    print(l*1e6, real(sqrt(e)), imag(sqrt(e)))

### END EXAMPLE ###
