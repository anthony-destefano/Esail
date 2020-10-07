#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 


def force(x, a, b, c, d):
	return np.log10( a / (((10**x)/c)**b + (10**x)**d) )

def force2(x, a, b, c,d):
	return np.log10( a * (10**x)**(b-1) *np.exp(-(10**x / c)**d))

filename = "ion_E_electron_Te.txt"

Ep, ETe, F = np.loadtxt(filename, unpack=True, usecols=(1,2,8))

Ep  = Ep.reshape((60, 50))
ETe = ETe.reshape((60, 50))
F   = F.reshape((60, 50))

for i in range(0,60):
	
	param, param_cov = curve_fit(force2, np.log10(Ep[0,:]), np.log10(F.T[:,i]), p0=(33., 1.5, 3000., 1.) )
	plt.loglog(Ep[0,:], F.T[:,i])
	plt.loglog(Ep[0,:], 10**force2(np.log10(Ep[0,:]), param[0], param[1], param[2], param[3]))
	plt.show()
	print(param)

plt.loglog(Ep[0,:], F.T)
plt.show()