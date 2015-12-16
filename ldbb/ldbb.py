#!/usr/bin/python
# ldbb.py 
# Copyright (C) 2011 Aaron Webster
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Returns the complex permittivity for different metals using the Drude,
# Lorentz-Drude, or Brendel-Bormann models.
#
# The implementaion of the Drude and Lorentz-Drude models is based on
# 'LD.m", courtesy of Bora Ung of Ecole Polytechnique de Montreal.
#
# All coefficents and formalism are from "Aleksandar D. Rakic, Aleksandra B.
# Djurisic, Jovan M. Elazar, and Marian L. Majewski. Optical properties of
# metallic films for vertical-cavity optoelectronic devices.  Appl. Opt.,
# 37(22):5271-5283, Aug 1998"

from scipy import constants
from scipy.special import wofz
from numpy import *
import sys, os

def LDBB(metal, model, lambda0):

	# supported metals
	metals = ['Ag', 'Au', 'Cu', 'Al', 'Be', 'Cr', 'Ni', 'Pd', 'Pt', 'Ti', 'W']
	# supported models:
	# D - Drude
	# LD - Lorentz-Drude
	# BB - Brebhal
	models = ['D', 'LD', 'BB' ]

	# material data for the LD model, from rakic et al
	LDdata = array([
	[9.010,9.030,10.83,14.98,18.51,10.75,15.92,9.72,9.59,7.29,13.22],
	[0.845,0.760,0.575,0.523,0.084,0.168,0.096,0.330,0.333,0.148,0.206],
	[0.048,0.053,0.030,0.047,0.035,0.047,0.048,0.008,0.080,0.082,0.064],
	[0.065,0.024,0.061,0.227,0.031,0.151,0.100,0.649,0.191,0.899,0.054],
	[3.886,0.241,0.378,0.333,1.664,3.175,4.511,2.950,0.517,2.276,0.530],
	[0.816,0.415,0.291,0.162,0.100,0.121,0.174,0.336,0.780,0.777,1.004],
	[0.124,0.010,0.104,0.050,0.140,0.150,0.135,0.121,0.659,0.393,0.166],
	[0.452,0.345,1.056,0.312,3.395,1.305,1.334,0.555,1.838,2.518,1.281],
	[4.481,0.830,2.957,1.544,1.032,0.543,0.582,0.501,1.314,1.545,1.917],
	[0.011,0.071,0.723,0.166,0.530,1.149,0.106,0.638,0.547,0.187,0.706],
	[0.065,0.870,3.213,1.351,4.454,2.676,2.178,4.621,3.668,1.663,3.332],
	[8.185,2.969,5.300,1.808,3.183,1.970,1.597,1.659,3.141,2.509,3.580],
	[0.840,0.601,0.638,0.030,0.130,0.825,0.729,0.453,3.576,0.001,2.590],
	[0.916,2.494,4.305,3.382,1.802,1.335,6.292,3.236,8.517,1.762,5.836],
	[9.083,4.304,11.18,3.473,4.604,8.775,6.089,5.715,9.249,19.43,7.498],
	[5.646,4.384,0,0,0,0,0,0,0,0,0],
	[2.419,2.214,0,0,0,0,0,0,0,0,0],
	[20.29,13.32,0,0,0,0,0,0,0,0,0]
	])

	# experimental data for the BB model, from rakic et al
	BBdata = array([
	[9.01,9.03,10.83,14.98,18.51,10.75,15.92,9.72,9.59,7.29,13.22],
	[0.821,0.770,0.562,0.526,0.081,0.154,0.083,0.330,0.333,0.126,0.197],
	[0.049,0.050,0.030,0.047,0.035,0.048,0.022,0.009,0.080,0.067,0.057],
	[0.050,0.054,0.076,0.213,0.066,0.338,0.357,0.769,0.186,0.427,0.006],
	[0.189,0.074,0.056,0.312,2.956,4.256,2.820,2.343,0.498,1.877,3.689],
	[2.025,0.218,0.416,0.163,0.131,0.281,0.317,0.066,0.782,1.459,0.481],
	[1.894,0.742,0.562,0.013,0.277,0.115,0.606,0.694,0.031,0.463,3.754],
	[0.133,0.050,0.081,0.060,0.067,0.261,0.039,0.093,0.665,0.218,0.022],
	[0.067,0.035,0.047,0.315,3.962,3.957,0.120,0.497,1.851,0.100,0.277],
	[5.185,2.885,2.849,1.561,0.469,0.584,1.059,0.502,1.317,2.661,0.985],
	[0.665,0.349,0.469,0.042,3.167,0.252,1.454,0.027,0.096,0.506,0.059],
	[0.051,0.312,0.324,0.182,0.346,0.817,0.127,0.309,0.551,0.513,0.136],
	[0.019,0.083,0.113,1.587,2.398,2.218,1.822,2.022,2.604,0.615,1.433],
	[4.343,4.069,4.819,1.827,2.827,1.919,4.583,2.432,3.189,0.805,1.962],
	[0.189,0.830,1.131,0.256,1.446,0.225,0.379,1.167,0.766,0.799,0.273],
	[0.467,0.719,0.726,0.014,0.311,0.105,0.654,0.409,2.214,0.0002,2.648],
	[0.117,0.125,0.172,2.145,3.904,6.983,6.637,0.119,2.891,4.109,4.555],
	[9.809,6.137,8.136,4.495,4.318,6.997,8.825,5.987,8.236,19.86,5.442],
	[1.170,1.246,1.719,1.735,0.893,4.903,0.510,1.331,1.146,2.854,1.912],
	[4.000,1.648,0.0,0,0,0,0,0,0,0,0],
	[0.052,0.179,0.0,0,0,0,0,0,0,0,0],
	[18.56,27.97,0.0,0,0,0,0,0,0,0,0],
	[0.516,1.795,0.0,0,0,0,0,0,0,0,0]
	])

	# select the right material, make sure it exists
	try:
		idx = metals.index(metal)
	except ValueError:
		print("Error:", metal, "is not a supported metal.")
		print("Supported metals:")
		print(metals)
		print("")
		raise

	# make sure the selected model exists
	try:
		test = models.index(model)
	except ValueError:
		print("Error:", model, "is not a supported model.")
		print("Supported Models:")
		print(models)
		print("")
		raise

	# chose appropriate model
	if model == "LD" or model == "D":
		materialdata = LDdata[:,idx]
		wp = materialdata[0]
		f0 = materialdata[1]
		gamma0 = materialdata[2]

		if idx < 2:
			fj = array(materialdata[3:18:3])
			gammaj = array(materialdata[4:18:3])
			wj = array(materialdata[5:18:3])
		# the rest have 4
		else:
			fj = array(materialdata[3:15:3])
			gammaj = array(materialdata[4:15:3])
			wj = array(materialdata[5:15:3])

		wp *= constants.e/constants.hbar
		wj *= (constants.e/constants.hbar)
		gamma0 *= (constants.e/constants.hbar)
		gammaj *= (constants.e/constants.hbar)

		# drude model contributions
		def geteps(lambda0):
			omega=2*pi*constants.c/lambda0
			epsilon_D = 1 - (f0*wp**2)/(omega*(omega + 1.0J*gamma0))
			if model == "D":
				return epsilon_D
			else:
				epsilon_L = sum((fj*wp**2)/((wj**2-omega**2)-1.0J*omega*gammaj))
				return epsilon_D + epsilon_L
		try:
			epsilon = [ geteps(x) for x in lambda0 ]
		except TypeError:
			epsilon = geteps(lambda0)

		return  epsilon

	# BB model
	elif model == "BB":
		materialdata = BBdata[:,idx]

		# numerical constants from article
		wp = materialdata[0]
		f0 = materialdata[1]
		gamma0 = materialdata[2]

		# Ag and Au have 5 parameters
		if idx < 2:
			fj = array(materialdata[3:23:4])
			gammaj = array(materialdata[4:23:4])
			wj = array(materialdata[5:23:4])
			sigmaj = array(materialdata[6:23:4])
		# the rest have 4
		else:
			fj = array(materialdata[3:19:4])
			gammaj = array(materialdata[4:19:4])
			wj = array(materialdata[5:19:4])
			sigmaj = array(materialdata[6:19:4])

		wp *= constants.e/constants.hbar
		gamma0 *= (constants.e/constants.hbar)
		gammaj *= (constants.e/constants.hbar)
		wj *= (constants.e/constants.hbar)
		sigmaj *=(constants.e/constants.hbar)
		omegap = sqrt(f0)*wp

		def geteps(lambda0):
			omega=2*pi*constants.c/lambda0
			epsf = 1 - omegap**2/(omega*(omega+gamma0*1.0J))
			aj = sqrt(omega**2 + 1.0J*omega*gammaj)
			zplus = (aj+wj)/(sqrt(2)*sigmaj)
			zminus = (aj-wj)/(sqrt(2)*sigmaj)
			epsb = (1.0J*sqrt(pi)*fj*wp**2/(2*sqrt(2)*aj*sigmaj))*(wofz(zplus)+wofz(zminus))
			eps = epsf + sum(epsb)
			return eps

		try:
			epsilon = [ geteps(x) for x in lambda0 ]
		except TypeError:
			epsilon = geteps(lambda0)

		return epsilon
