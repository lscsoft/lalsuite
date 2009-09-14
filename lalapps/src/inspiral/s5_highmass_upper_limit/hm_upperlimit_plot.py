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
# FIXME:  the plots should use LaTeX to generate the labels, but some don't
# work yet
matplotlib.rcParams.update({
	"text.usetex": False
})

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
		# FIXME maybe this should be smarter?
		# a total lack of injections can screw this up.  If vt = 0 move on. 
                if vt <= 0 or k <= 0: continue
                post *= vt / (1.0 + lam) * ( (1.0 + mu * vt / k)**(-k-1) + (mu * vt * lam * (1.0 + 1.0/k) /(1.0 + mu * vt / k)**(k+2)) )

        return mu, post

def integrate_posterior(mu, post, conf):
        cumpost = post.cumsum()/post.sum()
	#if you can do it, maybe you can't cause the volume is zero
        try: val = [idx for idx in range(len(cumpost)) if cumpost[idx] >= conf][0]
	except: val = 0
        return mu[val]

def get_mass_ranges(bins, mbA):
	mass1 = []
	mass2 = []
	for i,m1 in enumerate(bins.lower()[0]):
		for j,m2 in enumerate(bins.lower()[1]):
			if not mbA[i][j]:
				mass1.append(m1)
				mass1.append(bins.upper()[0][i])	
				mass2.append(m2)
				mass2.append(bins.upper()[1][j])
	return [min(mass1), max(mass1)], [min(mass2), max(mass2)]
				

# bins are the same for each call and ulA is an empty binnedArray that has the right shape
# that can hold the upperlimit when we get around to computing it later, so it is okay
# that bins, and ulA are overwritten in each call. vA, vA2 and dvA are the important ones
bins, vA, ulA = get_combined_array("2DsearchvolumeFirstMoment", 0)
#FIXME Hack to give volume that is zero a value = 0.01
vA[vA==0] = 0.01

bins, vA2, ulA = get_combined_array("2DsearchvolumeSecondMoment", 1)
bins, dvA, ulA = get_combined_array("2DsearchvolumeDerivative", 2)
bins, vAD, ulA = get_combined_array("2DsearchvolumeDistance", 3)
bins, ltA, tmp= get_combined_array("2DsearchvolumeLiveTime", 4)

# track livetime bins that are 0, they are outside of the search space
# this is a useful array for computing only meaningful things
mbA = numpy.zeros(ltA[0].shape)
mbA[ltA[0] == 0] = 1


m1range, m2range = get_mass_ranges(bins, mbA)

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
  #FIXME we shouldn't get ifos from filenames, we should put it in the xml :(
  ifos = sys.argv[1+i].replace('.xml','').replace("2Dsearchvolume-","")
  wiki = open(ifos+'range_summary.txt',"w")
  wiki.write("||Masses||Range(Mpc)||Time||UL @ 90%||Fractional error||\n")
  # loop over mass bins
  
  ulA.array *= 0

  for j, m1 in enumerate(bins.centres()[0]):
    for k, m2 in enumerate(bins.centres()[1]):
      masses = bins[m1,m2]
      if mbA[masses[0], masses[1]]: continue
      legend_str.append("%.1f, %.1f" % (m1, m2))
      mu,post = posterior(vA[i:i+1,masses[0],masses[1]], vA2[i:i+1,masses[0],masses[1]], dvA[i:i+1,masses[0],masses[1]])
      lines.append(pylab.loglog(mu,post/post.max()))
      ulA.array[j][k] = integrate_posterior(mu, post, 0.90)
      wiki.write("||%.2f,%.2f||%.2f|| %.2f || %.2e || %.2f ||\n" % (m1, m2, vAD[i,masses[0],masses[1]], ltA[i, masses[0], masses[1]], ulA.array[j][k], vA2[i, masses[0], masses[1]]**0.5 / (vA[i, masses[0], masses[1]])) )


  
  fudge = 0.01 * min (ulA.array[ulA.array !=0])
  log_vol = pylab.log10(vA[i])
  #HACKS FOR LOG PLOT :(
  log_ul = pylab.log10(ulA.array + fudge)
  vol_error = vA2[i]**0.5 / (vA[i] + 0.0001)
  der = dvA[i]
  log_der = pylab.log10(dvA[i] + 0.000000001)
  # set points outside mass space to 1 (log(1) = 0)
  log_ul[mbA == 1] = 0
  log_der[mbA == 1] = -3

  ##
  # Make posterior plots
  ##

  #f = pylab.figure(i)
  pylab.title("%s posteriors for a few mass bins" % (ifos,),fontsize=14)
  leg = pylab.figlegend(lines, legend_str, 'lower right')
  leg.prop.set_size(6)
  leg.pad = 0
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
  #FIXME setting these limits on the scale won't work in adv LIGO :)
  pylab.pcolor(X,Y, log_vol, vmin=0, vmax=10)
  pylab.colorbar()
  pylab.ylim(m1range)
  pylab.xlim(m2range)
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
  #FIXME we don't show errors greater than 100%
  pylab.pcolor(X,Y, vol_error, vmin=0, vmax=1)
  pylab.colorbar()
  pylab.ylim(m1range)
  pylab.xlim(m2range)
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

  #FIXME hard limits
  pylab.pcolor(X,Y, log_der, vmin=-3, vmax=3)
  pylab.colorbar()
  pylab.ylim(m1range)
  pylab.xlim(m2range)
  pylab.title(ifos + " Log10[fg/bg likelihood ratio(Lambda)]",fontsize=14)
  pylab.xlabel("Mass 2",fontsize=14)
  pylab.ylabel("Mass 1",fontsize=14)
  pylab.gca().set_aspect(1)
  pylab.grid()
  pylab.savefig(ifos + 'lambda.png')
  pylab.clf()

  ##
  # Make UL plot
  ##
  #FIXME hard limits at log(rate) of -7, 0
  pylab.pcolor(X,Y, log_ul, vmin=-7, vmax=0)
  pylab.colorbar()
  pylab.ylim(m1range)
  pylab.xlim(m2range)
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
ulA.array *= 0

for j, m1 in enumerate(bins.centres()[0]):
  for k, m2 in enumerate(bins.centres()[1]):
    masses = bins[m1,m2]
    if mbA[masses[0], masses[1]]: continue
    legend_str.append("%.1f, %.1f" % (m1, m2))
    mu,post = posterior(vA[...,masses[0],masses[1]], vA2[...,masses[0],masses[1]], dvA[...,masses[0],masses[1]])
    lines.append(pylab.loglog(mu,post/post.max()))
    ulA.array[j][k] = integrate_posterior(mu, post, 0.90)

#HACKS FOR LOG PLOT :(
# fudge is 1% effect
fudge = 0.01 * min (ulA.array[ulA.array !=0])
log_ul = pylab.log10(ulA.array + fudge)

# set points outside mass space to 1 (log(1) = 0)
log_ul[mbA == 1] = 0

# Make posterior plots
pylab.title("Combined posteriors for a few mass bins",fontsize=14)
leg = pylab.figlegend(lines, legend_str, 'lower right')
leg.prop.set_size(6)
leg.pad = 0
leg.ncol = 2
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
#FIXME hard limits at log(rate) of -7, 0
pylab.pcolor(X,Y, log_ul, vmin=-7.0, vmax=0)
pylab.colorbar()
pylab.ylim(m1range)
pylab.xlim(m2range)
pylab.title("Combined Log10[90% upper limit] in mergers/Mpc^3/yr",fontsize=14)
pylab.xlabel("Mass 2",fontsize=14)
pylab.ylabel("Mass 1",fontsize=14)
pylab.gca().set_aspect(1)
pylab.grid()
pylab.savefig('combinedupper_limit.png')
pylab.clf()


print >> sys.stderr, "ALL FINNISH!"
sys.exit(0)
