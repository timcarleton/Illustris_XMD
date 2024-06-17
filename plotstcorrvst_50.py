import numpy as np
import matplotlib.pyplot as plt
from astropy import cosmology
astrocosmo = cosmology.FlatLambdaCDM(H0=67.74,Om0=0.307)

plt.style.use('../python/stylesheet.txt')

dat0=np.loadtxt('corr_all_50.txt')

dat1=np.loadtxt('corr_1.txt')
dat2=np.loadtxt('corr_2.txt')
dat3=np.loadtxt('corr_3.txt')
dat4=np.loadtxt('corr_4.txt')

zs=np.array([20.05,14.99,11.98,10.976,10.0,9.389,9.002,8.449,8.012,7.595,7.236,7.005,6.491,6.011,5.847,5.530,5.228,4.996,4.665,4.428,4.177,4.008,3.709,3.491,3.283,3.008,2.896,2.733,2.578,2.444,2.316,2.208,2.103,2.002,1.904,1.823,1.743,1.667,1.604,1.531,1.496,1.414,1.358,1.302,1.248,1.206,1.155,1.114,1.074,1.036,0.997,0.951,0.923,0.887,0.851,0.817,0.791,0.757,0.733,0.700,0.676,0.644,0.621,0.598,0.576,0.546,0.525,0.503,0.482,0.461,0.440,0.420,0.40,0.380,0.361,0.348,0.329,0.310,0.298,0.273,0.261,0.244,0.226,0.214,0.1973,0.1804,0.1693,0.1527,0.1419,0.1258,0.1099,0.0994,0.0839,0.0737,0.0585,0.0485,0.0337,0.0240,0.0095,0])
age=astrocosmo.age(zs)

plt.clf()
plt.plot(age[-1]-age,dat4,'C3D')
plt.plot(age[-1]-age,dat3,'C2s')
plt.plot(age[-1]-age,dat2,'C1x')
plt.plot(age[-1]-age,dat1,'C0*')
plt.xlabel(r'$M_*$ formed in last $t$ Gyr')
plt.ylabel(r'$M_*-Z_{\rm SFR}$ Correlation Coefficient')
plt.savefig('corretimescale_50.png')

plt.clf()
plt.plot(age[-1]-age,dat0,'C0o')
plt.xlabel(r'$t$ Gyr')
plt.ylabel(r'Correlation Coeff.: ${M_*{\rm ~formed~in~last}~t~{\rm Gyr}}-Z_{\rm SFR}$')
plt.savefig('corretimescaleall_50.png')
