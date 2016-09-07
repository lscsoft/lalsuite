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
from scipy import special
from scipy.misc import logsumexp
import time
import random
from scipy.optimize import curve_fit
import collections
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d
import types
from scipy.interpolate import lagrange
from glue.text_progress_bar import ProgressBar
#from shapely.ops import polygonize
#from shapely.geometry import LineString, MultiPolygon, MultiPoint, Point, Polygon

def find_ln_p_j_voronoi(m, f, popt):
    vor = Voronoi(m)
    all_vols = []
    all_p = []
    all_i = []
    vtot = 0
    new_points = []
    for i, p in enumerate(vor.points):
        region = vor.regions[vor.point_region[i]]
        if -1 not in region:
            pvol = vol(vor, region)
            if pvol > 0.05:
                continue
            all_vols.append(pvol)
            all_p.append(p)
            all_i.append(i)
            vtot += pvol
    print "Voronoi diagram, total volume = "+str(vtot)
    all_p = np.array(all_p)
    all_vols = np.array(all_vols)
    all_i = np.array(all_i)
    p_m1 = f(all_p[:,0],*popt)/np.trapz(f(all_p[:,0][np.argsort(all_p[:,0])],*popt),all_p[:,0][np.argsort(all_p[:,0])])
    p_m2 = f(all_p[:,1],*popt)/np.trapz(f(all_p[:,1][np.argsort(all_p[:,1])],*popt),all_p[:,1][np.argsort(all_p[:,1])])
    p_tj = p_m1*p_m2*all_vols
    p_tj = p_tj/sum(p_tj)
    return np.log(p_tj), all_p, all_i

def trivol((a, b, c)):
    # Calculates area of triangle
    return abs(np.cross((c-a),(c-b)))/2.

def tetravol((a, b, c, d)):
    # Calculates volume of tetrahedron, given vertices a, b, c, d (triplets)
    return abs(np.dot((a-d), np.cross(b-d),(c-d)))/6.

def vol(vor, region):
    dpoints = vor.vertices[region]
    return sum(trivol(dpoints[simplex]) for simplex in Delaunay(dpoints).simplices)

def ln_p_k(ovrlp, rho, t_k, acc=0.001):
    # Finds probability of template t_k fitting data, given that the signal is rho*t_j
    # ovrlp = vector of the overlaps of (t_j, t_k) for one t_j
    # rho = SNR (float)
    # t_k = template
    ln_num = ln_p_k_num(rho*ovrlp[t_k]) # compute the numerator of p_k
    # for full template bank, don't need ln_p_k_den. just need to do a logexpsum of ln_num.
    ln_den = ln_p_k_den(ovrlp, rho, acc=acc) # compute the denominator of p_k
    return ln_num-ln_den

def ln_p_k_num(x, sqrtpiover2 = np.sqrt(np.pi/2.), sqrt2 = np.sqrt(2.)):
    halfxsquared = 0.5*x**2.
    lny = halfxsquared + np.log( sqrtpiover2*(x**3.+3.*x)*(1.+special.erf(x/sqrt2))+np.exp(-halfxsquared)*(x**2.+2.) ) # N = 4
    #lny = halfxsquared + np.log( sqrtpiover2*(x**2.+1.)*(1.+special.erf(x/sqrt2))+np.exp(-halfxsquared)*x ) # N = 3
    return lny

def ln_p_k_den(tjtk, rho, acc=0.001):
    # Finds the denominator part of the probability P(t_k | t_j) for N=4 dimensions (m1, m2, chi_eff, rho)
    # P(t_k | t_j) is a template for the signal rho*t_j (tjtk is the overlap between t_k and t_j)
    # acc = accuracy we wish to obtain
    # Denominator is the sum of all numerator terms
    x = rho*tjtk
    if rho < 10: # must process full array
        return logsumexp(ln_p_k_num(x))
    x.sort()
    lny_list = []
    limit = np.log(acc/len(x))
    for w in x[::-1]:
        lny_list.append(ln_p_k_num(w))
        if lny_list[-1]-lny_list[0] < limit:
            break
    return logsumexp(lny_list)

def source_population(srcfile):
    data = np.loadtxt(srcfile,unpack=True)
    arr = np.zeros(np.shape(data))
    arr[0] = np.linspace(0.5,2.5,len(data[0])) # masses
    arr[1] = np.cumsum(data[1]) # probability density
    x = arr[0]
    y = arr[1]/np.trapz(arr[1],x)
    f = lambda x, *p: p[0]*np.exp(-p[1]*(x-p[2])**2)
    popt, pcov = curve_fit(f, x, y, [2.5,1.,1.3]) # fit p(m) to a Gaussian curve
    return f, popt

def load_overlaps(h5_file):
    return h5_file[h5_file.keys()[0]]['overlaps'].value

def load_bank(h5_file):
    # return masses and chis
    return np.vstack((h5_file[h5_file.keys()[0]]['mass1'].value, h5_file[h5_file.keys()[0]]['mass2'].value)).T
    #return np.vstack((h5_file[h5_file.keys()[0]]['mass1'].value, h5_file[h5_file.keys()[0]]['mass2'].value, h5_file[h5_file.keys()[0]]['spin1z'].value, h5_file[h5_file.keys()[0]]['spin2z'].value))
