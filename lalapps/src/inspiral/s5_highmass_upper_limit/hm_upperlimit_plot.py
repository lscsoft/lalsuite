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

# bins are the same for each call and ulA is an empty binnedArray that has the right shape
# that can hold the upperlimit when we get around to computing it later, so it is okay
# that bins, and ulA are overwritten in each call. vA, vA2 and dvA are the important ones
bins, vA, ulA = get_combined_array("2DsearchvolumeFirstMoment", 0)
#FIXME Hack to give volume that is zero a value = 0.01

vA[vA==0] = 0.01

bins, vA2, ulA = get_combined_array("2DsearchvolumeSecondMoment", 1)
bins, dvA, ulA = get_combined_array("2DsearchvolumeDerivative", 2)
bins, vAD, ulA = get_combined_array("2DsearchvolumeDistance", 3)
 
#bin edges Number of bins + 1 for pcolor
X = numpy.array( list(bins.lower()[0]) + [bins.upper()[0][-1]] )
Y = numpy.array( list(bins.lower()[1]) + [bins.upper()[1][-1]] )

numfiles = len(sys.argv[1:])
f = pylab.figure(1)
###############################################################################
# loop over all the filenames and masses and compute the posteriors separately
###############################################################################
for i, f in enumerate(sys.argv[1:]):
  legend_str = []
  lines = []
  #FIXME we can't get ifos from filenames, we should put it in the xml :(
  ifos = "_".join(sys.argv[1+i].split('-')[-2:]).replace('.xml','')
  wiki = open(ifos+'range_summary.txt',"w")
  wiki.write("||Masses||Range||\n")
  # loop over mass bins
  for j, m1 in enumerate(bins.centres()[0]):
    for k, m2 in enumerate(bins.centres()[1]):
      masses = bins[m1,m2]
      if vAD[i,masses[0],masses[1]] == 0: continue
      legend_str.append("%.1f, %.1f" % (m1, m2))
      wiki.write("||%.2f,%.2f||%.2f||\n" % (m1, m2, vAD[i,masses[0],masses[1]] ) )
      mu,post = posterior(vA[i:i+1,masses[0],masses[1]], vA2[i:i+1,masses[0],masses[1]], dvA[i:i+1,masses[0],masses[1]])
      lines.append(pylab.loglog(mu,post/post.max()))
      ulA.array[j][k] = integrate_posterior(mu, post, 0.90)

  
  fudge = 0.01 * min (ulA.array[ulA.array !=0])
  print fudge
  log_vol = pylab.log10(vA[i])
  #HACKS FOR LOG PLOT :(
  log_ul = pylab.log10(ulA.array + fudge)
  vol_error = vA2[i]**0.5 / (vA[i] + 0.0001)
  der = dvA[i]

  ##
  # Make posterior plots
  ##

  #f = pylab.figure(i)
  pylab.title("%s posteriors for a few mass bins" % (ifos,),fontsize=14)
  leg = pylab.figlegend(lines, legend_str, 'lower right')
  leg.prop.set_size(8)
  pylab.ylabel("Prob (unnormalized)",fontsize=14)
  pylab.xlabel("Rate",fontsize=14)
  pylab.ylim([0.0001, 1])
  pylab.grid()
  #FIXME hardcoded rate limits are bad for advanced ligo
  pylab.xlim([1e-8, 1])
  pylab.savefig(ifos + 'posterior.png')
  pylab.clf()

  ##
  # Make log volume plot
  ##

  #FIXME make it gray for pub?
  #pylab.gray()
  pylab.pcolor(X,Y, log_vol)
  pylab.colorbar()
  pylab.ylim([0, 51])
  pylab.xlim([11, 101])
  pylab.title(ifos + " Log10[< Volume * Time>] in mergers/Mpc^3/yr",fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(ifos+'volume_time.png')
  pylab.clf()

  ## 
  # Make vol error plot
  ##
  pylab.pcolor(X,Y, vol_error)
  pylab.colorbar()
  pylab.ylim([0, 51])
  pylab.xlim([11, 101])
  pylab.title(ifos + " Fractional Error on Volume * Time [std/mean]",fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(ifos+'fractional_error.png')
  pylab.clf()

  ##
  # Make lambda plot
  ##

  pylab.pcolor(X,Y, der)
  pylab.colorbar()
  pylab.ylim([0, 51])
  pylab.xlim([11, 101])
  pylab.title(ifos + " Volume derivative, Lambda",fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(ifos + 'lambda.png')
  pylab.clf()

  ##
  # Make UL plot
  ##
  pylab.pcolor(X,Y, log_ul)
  pylab.colorbar()
  pylab.ylim([0, 51])
  pylab.xlim([11, 101])
  pylab.title(ifos + " Log10[90% upper limit] in mergers/Mpc^3/yr",fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(ifos + 'upper_limit.png')
  pylab.clf()


###############################################################################
# now write out the special combined case
###############################################################################
lines = []
legend_str = []
for j, m1 in enumerate(bins.centres()[0]):
  for k, m2 in enumerate(bins.centres()[1]):
    masses = bins[m1,m2]
    if vAD[i,masses[0],masses[1]] == 0: continue
    legend_str.append("%.1f, %.1f" % (m1, m2))
    mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
    lines.append(pylab.loglog(mu,post/post.max()))
    ulA.array[j][k] = integrate_posterior(mu, post, 0.90)

# Make posterior plots
pylab.title("Combined posteriors for a few mass bins",fontsize=14)
leg = pylab.figlegend(lines, legend_str, 'lower right')
leg.prop.set_size(8)
pylab.ylabel("Prob (unnormalized)",fontsize=14)
pylab.xlabel("Rate",fontsize=14)
pylab.ylim([0.0001, 1])
pylab.grid()
#FIXME hardcoded rate limits are bad for advanced ligo
pylab.xlim([1e-8, 1])
pylab.savefig('combinedposterior.png')
pylab.clf()


##
# Make UL plot
##
pylab.pcolor(X,Y, log_ul)
pylab.colorbar()
pylab.ylim([0, 51])
pylab.xlim([11, 101])
pylab.title("Combined Log10[90% upper limit] in mergers/Mpc^3/yr",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig('combinedupper_limit.png')
pylab.clf()


print >> sys.stderr, "ALL FINNISH!"
sys.exit(0)
