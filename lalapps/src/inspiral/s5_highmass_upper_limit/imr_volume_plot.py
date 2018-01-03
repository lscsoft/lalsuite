#!/usr/bin/python
import sys
from pylal import rate
from pylal import SimInspiralUtils
import scipy
import numpy
import matplotlib
matplotlib.use('Agg')
import pylab
from math import *
import sys
import glob
import copy
from glue.ligolw.utils import ligolw_add
from glue.ligolw import ligolw, utils

# FIXME, I apparently don't know how to read XML as cleverly as I thought, how do I find the table
# without specifying the childnode number?
def get_combined_array(tablename, childnode):
  # FIXME assumes that all the xml files have the same binned array tables
  # Figure out the shape of the arrays in the file, make an array with one more
  # dimension, the number of files from sys.argv[1:]
  xmldoc = utils.load_filename(sys.argv[1], verbose=True, gz = (sys.argv[1] or "stdin").endswith(".gz"))
  xmldoc = xmldoc.childNodes[0]
  A  = rate.binned_array_from_xml(xmldoc.childNodes[childnode], tablename) 
  bins = rate.bins_from_xml(xmldoc.childNodes[childnode])
  out = numpy.zeros((len(sys.argv[1:]),)+A.array.shape,dtype="float")
  # Read the data
  for i, f in enumerate(sys.argv[1:]):
    xmldoc = utils.load_filename(f, verbose=True, gz = (f or "stdin").endswith(".gz"))
    xmldoc = xmldoc.childNodes[0]
    out[i] = rate.binned_array_from_xml(xmldoc.childNodes[childnode], tablename).array
  A.array = numpy.zeros(A.array.shape)
  return bins, out, A 


def istrue(arg):
  return True


def posterior(VT, sigmasq, Lambda):
        '''
        This function implements the analytic marginalization in 
        Biswas, Creighton, Brady, Fairhurst (25)
        Cause hey, why not?
        This takes arrays of VT, sigma and Lambda to combine.
        '''

        length = 100000
        #FIXME does this need to be an integer???
        K = VT**2 / sigmasq

        #FIXME, drew said this was cool?
        mu = numpy.arange(length) * 100.0 / VT.sum() / length

        post = numpy.ones(len(mu), dtype="float")

        for vt, k, lam in zip(VT, K, Lambda):
                post *= vt / (1.0 + lam) * ( (1.0 + mu * vt / k)**(-k-1) + (mu * vt * lam * (1.0 + 1.0/k) /(1.0 + mu * vt / k)**(k+2)) )

        return mu, post

def integrate_posterior(mu, post, conf):
        cumpost = post.cumsum()/post.sum()
	#if you can do it, maybe you can't cause the volume is zero
        try: val = [idx for idx in range(len(cumpost)) if cumpost[idx] >= conf][0]
	except: val = 0
        return mu[val]

# test case

#VT = numpy.array([10.0**8])
#sigmasq = numpy.array([4.0*10**14])
#Lambda = numpy.array([1.0])

#mu, post = posterior(VT, sigmasq, Lambda)

#pylab.semilogy(mu, post.cumsum()/post.sum())
#pylab.show()
#print integrate_posterior(mu, post, 0.90)

# bins is the same for each call and ulA is an empty binnedArray that has the right shape
# that can hold the upperlimit when we get around to computing it later, so it is okay
# that bins, and ulA are overwritten in each call. vA, vA2 and dvA are the important ones
bins, vA, ulA = get_combined_array("2DsearchvolumeFirstMoment", 0)
#FIXME Hack to give volume that is zero a value = 0.01
vA[vA==0] = 0.01
for i, l in enumerate(vA):
  for j, m in enumerate(l):
    for k, n in enumerate(m):
      if n == 0: print i,j,k,n #vA[i][j][k] = 0.01
    #if k = 0 k = 0.01
print vA.shape
bins, vA2, ulA = get_combined_array("2DsearchvolumeSecondMoment", 1)
#bins, dvA, ulA = get_combined_array("2DsearchvolumeDerivative", 2)
#bins, vAD, ulA = get_combined_array("2DsearchvolumeDistance", 3)
 
#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(bins.lower()[0]) + [bins.upper()[0][-1]] )
Y = numpy.array( list(bins.lower()[1]) + [bins.upper()[1][-1]] )

#compute combined posterior over m1, m2
#for m1 in range(len(bins.lower()[0])):
#  for m2 in range(len(bins.lower()[1])): 
#    mu, post = posterior(vA[...,m1,m2], vA2[...,m1,m2], dvA[...,m1,m2])
#    ulA.array[m1][m2] = integrate_posterior(mu, post, 0.90)

log_vol = pylab.log10(vA[0])

#log_ul = pylab.log10(ulA.array)

vol_error = vA2[0]**0.5 / (vA[0] + 0.0001)

#der = dvA[0] #pylab.log10(eA.array)

fn = sys.argv[1]

pylab.figure(1)
#pylab.gray()
pylab.pcolor(X,Y, log_vol, vmin=0, vmax=10)
#pylab.hold(1)
#cn = pylab.contour(bins.centres()[0], bins.centres()[1], ul)
#pylab.clabel(cn)
#pylab.hold(0)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Log10[< Volume >] Mpc^3",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig(fn.split('-')[-1].replace('.xml','volume_time.png'))

#pylab.show()

pylab.figure(2)
pylab.pcolor(X,Y, vol_error, vmin=0, vmax=2)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Fractional Error on Volume [std/mean]",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig(fn.split('-')[-1].replace('.xml','fractional_error.png'))

if vA.shape[0] == 2:
  pylab.figure(3)
  #pylab.gray()
  vA0 = vA[0]
  print vA0.min()
  vA0[vA0 == vA0.min()] = 1.0
  vA1 = vA[1]
  vA1[vA1 == vA1.min()] = 1.0
  print vA1.min()
  pylab.pcolor(X,Y, pylab.log10(vA0 /  vA1), vmin=-1, vmax=1)
  #pylab.hold(1)
  #cn = pylab.contour(bins.centres()[0], bins.centres()[1], ul)
  #pylab.clabel(cn)
  #pylab.hold(0)
  pylab.colorbar()
  pylab.ylim([0, 51])
  pylab.xlim([11, 101])
  pylab.title("Log10[Volume \n" + sys.argv[1] + "/\n" + sys.argv[2] +']',fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(sys.argv[1].replace('/','_') + "_" + sys.argv[2].replace('/','_') + 'volume_time_ratio.png')

