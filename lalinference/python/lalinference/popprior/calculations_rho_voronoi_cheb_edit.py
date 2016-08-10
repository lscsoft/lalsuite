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

def source_population():
    data = np.loadtxt("1309_6635_DNS_Posterior.txt",unpack=True)
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
    
############################################
############################################

start_time = time.time()

plots=False 
datetime=time.strftime('%Y%m%d')
bank = h5py.File("overlap_first2years.hdf","r")
save_data = True

f, popt = source_population() # find probability density of masses for BNS source population
mass = load_bank(bank) # load mass
#mass = load_bank(bank) # load mass and chi
overlap = load_overlaps(bank)
print "Overlap data loaded. Total time elapsed:", time.time()-start_time

############################################
############################################
    
# CALCULATE PROBABILITIES

n = 10
k = np.arange(1,n+1)
rho = 30*np.cos((2*k-1)*np.pi/(2*n)) + 30 # min(rho)~0, max(rho)~80
rho.sort()
rho = np.insert(rho,0,0)

ln_p_j, p_j_templates, p_j_indices = find_ln_p_j_voronoi(mass, f, popt) #p_j_indices = indices of the templates for m1, m2 arrays (eg. p_j_indices[0] = 1, which means p_j_template[0] corresponds to the 1th template in m1, m2)
t_k = np.array(range(len(mass)))
ln_p_jk = np.log(np.zeros((len(rho), len(t_k)))) # P(signal t_j is recovered by template t_k)

#p_j_indices.tolist().sort(key=ln_p_j.__getitem__) # doing loop in order of p_j (smallest to largest) for numerical accuracy
order = np.argsort(ln_p_j)
ln_p_j, p_j_templates, p_j_indices = ln_p_j[order], p_j_templates[order], p_j_indices[order]

progress = ProgressBar(max=len(p_j_indices))
for i in range(len(p_j_indices)): # loop over all signal population
    progress.increment(text=str(p_j_templates[order][i]))
    ovrlp = overlap[p_j_indices[i]]
    for r in range(len(rho)): # loop over all rho
        ln_p_jk[r,:] = np.logaddexp(ln_p_jk[r,:], ln_p_j[order][i]+ln_p_k(ovrlp, rho[r], t_k))
        
print "ln(P_jk) computed for all templates. Time elapsed:", time.time()-start_time

# Save data to hdf5 file
if save_data:
    directory = ""
    filename = "logP_vs_rho_"+datetime
    f = h5py.File(directory+filename+".hdf5","w")
    f.create_dataset('rho', data=rho)
    f.create_dataset('ln_P_jk', data=ln_p_jk)
    f.create_dataset('masses', data=mass)
    f.close()
    print "Data saved. Time elapsed:", time.time()-start_time

