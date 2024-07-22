#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset

plt.figure(0)
iky = None
fname = sys.argv[1]

data = Dataset("%s"%fname, mode='r')
t = data.groups['Grids'].variables['time'][:]
ky = data.groups['Grids'].variables['ky'][:]
tNorm = 10

if iky is None:
    for i in np.arange(0, len(ky)):
      y = data.groups['Diagnostics'].variables['Phi2_kyt'][:,i]
      if i==0:
         fmt = '--'
      else:
         fmt = '-'
      y_norm = np.interp( tNorm, t, y )
      normalised_y = y/y_norm # normalise at a fixed time
      plt.plot(t, normalised_y, fmt, label='ky = %.3f' % ky[i])
else:
    y = data.groups['Diagnostics'].variables['Phi2_kyt'][:,iky]
    plt.plot(t, y, '-', label='ky = %.3f' % ky[iky])

analytic_data = np.loadtxt( sys.argv[2] )
t_analytic = analytic_data[:,0]/np.sqrt(2)
phi2_analytic = analytic_data[:,1]

analytic_norm = np.interp( tNorm, t_analytic, phi2_analytic )

plt.plot( t_analytic, 1.1*phi2_analytic/analytic_norm, '-', label='Analytic Data')

plt.yscale('log')
plt.legend()

#plt.figure(1)
#avg = np.mean(data.groups['Spectra'].variables['Pkyst'][:,0,:], axis=0)
#plt.plot(ky, avg, 'o-')
#plt.xscale('log')
#plt.yscale('log')
plt.show()
