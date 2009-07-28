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

VT = numpy.array([10.0**8])
sigmasq = numpy.array([4.0*10**14])
Lambda = numpy.array([1.0])

mu, post = posterior(VT, sigmasq, Lambda)

#pylab.semilogy(mu, post.cumsum()/post.sum())
#pylab.show()
print integrate_posterior(mu, post, 0.90)

# bins is the same for each call and ulA is an empty binnedArray that has the right shape
# that can hold the upperlimit when we get around to computing it later, so it is okay
# that bins, and ulA are overwritten in each call. vA, vA2 and dvA are the important ones
bins, vA, ulA = get_combined_array("2DsearchvolumeFirstMoment", 0)
bins, vA2, ulA = get_combined_array("2DsearchvolumeSecondMoment", 1)
bins, dvA, ulA = get_combined_array("2DsearchvolumeDerivative", 2)

 
#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(bins.lower()[0]) + [bins.upper()[0][-1]] )
Y = numpy.array( list(bins.lower()[1]) + [bins.upper()[1][-1]] )

#compute combined posterior over m1, m2
for m1 in range(len(bins.lower()[0])):
  for m2 in range(len(bins.lower()[1])): 
    mu, post = posterior(vA[...,m1,m2], vA2[...,m1,m2], dvA[...,m1,m2])
    ulA.array[m1][m2] = integrate_posterior(mu, post, 0.90)

log_vol = pylab.log10(vA[0])
log_ul = pylab.log10(ulA.array)


vol_error = vA2[0]**0.5 / (vA[0] + 0.0001)

der = dvA[0] #pylab.log10(eA.array)

fn = sys.argv[1]

pylab.figure(1)
masses = bins[15,15]
print masses
mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
pylab.loglog(mu,post/post.max())
pylab.hold(1)
masses = bins[50,50]
print masses
mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
pylab.loglog(mu,post/post.max())
masses = bins[1,99]
print masses
mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
pylab.loglog(mu,post/post.max())
masses = bins[1,24]
print masses
mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
pylab.loglog(mu,post/post.max())
pylab.hold(0)
pylab.title("Combined posteriors for a few mass bins",fontsize=14)
pylab.legend(["15,15", "50,50", "1,99", "1,24"])
pylab.ylabel("Prob (unnormalized)",fontsize=14)
pylab.xlabel("Rate",fontsize=14)
pylab.ylim([0.0001, 1])
pylab.grid()
if len(sys.argv) == 2: pylab.savefig(fn.split('-')[-1].replace('.xml','posterior.png'))
else: pylab.savefig("combinedposterior.png")

pylab.figure(2)
pylab.gray()
pylab.pcolor(X,Y, log_vol)
#pylab.hold(1)
#cn = pylab.contour(bins.centres()[0], bins.centres()[1], ul)
#pylab.clabel(cn)
#pylab.hold(0)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Log10[< Volume * Time>] in mergers/Mpc^3/yr",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig(fn.split('-')[-1].replace('.xml','volume_time.png'))

#pylab.show()

pylab.figure(3)
pylab.pcolor(X,Y, vol_error)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Fractional Error on Volume * Time [std/mean]",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig(fn.split('-')[-1].replace('.xml','fractional_error.png'))


pylab.figure(4)
pylab.pcolor(X,Y, der )
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Volume derivative, Lambda",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig(fn.split('-')[-1].replace('.xml','lambda.png'))


pylab.figure(5)
pylab.gray()
pylab.pcolor(X,Y, log_ul)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Log10[90% upper limit] in mergers/Mpc^3/yr",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
if len(sys.argv) == 2: pylab.savefig(fn.split('-')[-1].replace('.xml','upper_limit.png'))
else: pylab.savefig("combinedupper_limit.png")


#pylab.show()


