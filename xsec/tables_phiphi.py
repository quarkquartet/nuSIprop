"""This code numerically computes the absorption integrals for double scalar production

When interpolating this, it's better to interpolate gamma as a function of log10(sbar_minus)
This will increase the accuracy of spline interpolation, as the function changes logarithmically
"""

import numpy as np
import scipy.integrate as integ
import warnings
import ctypes
import os

lib_C = ctypes.cdll.LoadLibrary(os.path.abspath('funcs.so'))

dsigma_over_tauphi = lib_C.dsigma_phiphi_over_tauphi
dsigma_over_tauphi.restype = ctypes.c_double
dsigma_over_tauphi.argtypes = (ctypes.c_int, ctypes.c_double)

warnings.filterwarnings('error') # Convert warnings into exceptions

tbar_plus = np.geomspace(-1e4, -4, 5000)[::-1]
log10_delta = np.linspace(0.005, 0.05, 100)

with open("alphatilde_phiphi.dat", 'w') as ofile:
    ofile.write("#|tbar_plus|  log10(delta)  \int_{tbar_plus}^{tbar_minus} dtbar \int_{-tbar}^{-tbar_plus} dsbar \int_{tauphibar_min}^tbar dtauphibar/(-tauphibar) dsigma/dtauphibar\n")
    for t in tbar_plus:
        for d in log10_delta:
            delta = 10**d
            tplus = t
            tmin = tplus/delta

            integral = integ.dblquad(dsigma_over_tauphi,
                                     tplus, tmin,
                                     lambda t: max(-t, 4, -t**2/(1+t)), lambda t: -tplus)[0]

            ofile.write("%.7e %.7e %.7e\n" % (abs(t), d, integral))
            ofile.flush()

sbar_plus = np.geomspace(4, 1e4, 1000)
log_sbar_minus_over_tbar_minus_over_log_delta = np.arange(1, 1001)
log10_delta = np.linspace(0.005, 0.05, 100)

with open("alpha_phiphi.dat", 'w') as ofile:
    ofile.write("#sbar_plus    log10(delta) log(sbar_minus/-tbar_minus)/log(delta)   \int_{tbar_plus}^{tbar_minus} dtbar \int_{sbar_minus}^{sbar_plus} dsbar \int_{tauphibar_min}^tbar dtauphibar/(-tauphibar) dsigma/dtauphibar\n")
    for s in sbar_plus:
        for n in log_sbar_minus_over_tbar_minus_over_log_delta:
            for d in log10_delta:
                delta = 10**d
                splus = s
                smin = splus/delta
                tmin = - smin/(delta**n)
                tplus = tmin*delta
                integral = integ.dblquad(dsigma_over_tauphi,
                                         tplus, tmin,
                                         lambda t: max(smin, 4), lambda t: splus)[0]
                if integral < 1e-37:
                    integral = 0
                    
                ofile.write("%.7e %4u %.7e %.7e\n" % (s, n, d, integral))
