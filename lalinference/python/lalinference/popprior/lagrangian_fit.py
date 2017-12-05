# Copyright (C) 2016  Heather Fong, Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# copied from calculations_updated.py, restructured to input different rho for select t_k
# sample mass distribution directly from template bank
# new P_tj calculation: using Voronoi diagram
# use Chebyshev nodes to compute values of rho
# Edited with Kipp, July 25

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
from scipy.misc import logsumexp
from time import time
import random
from scipy.optimize import curve_fit
from scipy.interpolate import lagrange

def open_hdf5(filename):
    f = h5py.File(filename)
    return f['rho'].value, f['ln_P_jk'].value, f['masses'].value

start_time = time()

rho, lnP, mass = open_hdf5("logP_vs_rho_cheb_p_j_only_20160727.hdf5")
rho_interp = np.linspace(0,60,50)
rho_interp = np.append(rho_interp, rho)
rho_interp.sort()
cutoff = -3
idx = np.where(rho[cutoff]==rho_interp)[0][0]
#fits = np.zeros((len(rho_interp), len(lnP[0])))
#for k in range(len(lnP[0])):
#    fits[:idx,k]=lagrange(rho, lnP[:,k])(rho_interp)[:idx]
#    fits[idx:,k] = np.full(len(fits[idx:,k]),lnP[cutoff,k])

for i in range(len(lnP[0])):
    #plt.plot(rho_interp, np.exp(fits[:,i]))
    plt.plot(rho, np.exp(lnP[:,i]), '-')

#plt.xlim([0,6])
plt.yscale('log')

print "Elapsed time:", time()-start_time
'''
idx = []
for i in range(len(lnP[0])):
    if np.exp(lnP[-3][i]) < 1e-90:
        idx.append(i)

c= ['blue', 'green', 'red', 'magenta', 'orange']
for i in range(len(idx)):
    plt.plot(rho, np.exp(lnP[:,idx[i]]), 'o-', color=c[i], label=mass[idx[i]])

plt.legend(loc="lower left")
plt.yscale('log')
plt.xlim([0,6])
plt.show()
'''
