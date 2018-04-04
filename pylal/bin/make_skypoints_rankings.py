#!/usr/bin/python

from pylal import git_version

__author__ = "Larry Price <larry.price@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


import sys
import cPickle
import glob
from math import sqrt
from optparse import *
from pylal import skylocutils
import numpy as np


def compute_bw(data,npts):
  """
  compute the 'optimal' bandwidth for kernel density estimator
  currently uses Silverman's rule of thumb
  """

  return 1.0592 * min(np.std(data),skylocutils.iqr(data)/1.34)*pow(npts,-1./5.)

def compute_nbins(data,npts):
  """
  compute the number of bins appropriate for a histogram
  currently uses Scott's method
  """
  span = max(data)-min(data)
  denom = 3.5 * np.std(data)*pow(npts,-(1./3.))
  return int(np.ceil(span/denom))


usage = """
usage: %prog [options]

Create a pickle of rankings for use by run_skypoints.py.
See ligovirgocvs/cbc/protected/projects/s6/sky_localization for an example.

"""

def parse_command_line():
  """
  Parser function dedicated
  """
  parser = OptionParser(usage=usage)
  parser.add_option("-g","--glob",action="store",type="string",\
    default=None, metavar=" GLOB",help="GLOB of coire or coinc (xml) files to read" )
  parser.add_option("-z","--input-type",action="store",default="coire",\
    help="specify the type of input in the glob.  valid options are coinctable (DEFAULT) and coire")
  parser.add_option("-s","--snr-threshold",action="store_true",\
    default=False, help="use snr-dependent quantities for timing and effective distance" )
  parser.add_option("-f","--reference-frequency",action="store",type="float", default=0.0, metavar=" REFERENCE_FREQUENCY", \
    help="reference frequency for signal timing" )
  parser.add_option("-N","--Npts",action="store",type="float", default=10000, metavar=" NUM_POINTS", \
    help="number of points to use in pdf estimation" )
  parser.add_option("-p","--plot-pdfs",action="store_true",default=False, \
    help="show plots of the pdfs obtianed using a kernel density estimator against normalized histograms" )
  parser.add_option("-F","--output-file",action="store",type="string",\
    default='rankings.pkl', metavar=" OUTFILE",help="name of the file to write output to" )


  (options,args) = parser.parse_args()

  return options, sys.argv[1:]

opts, args = parse_command_line()

#deal with the glob 
files = []
if opts.glob is not None:
  for gl in opts.glob.split(" "):
    files.extend(glob.glob(gl))
  if len(files) < 1:
    print >>sys.stderr, "The glob for " + opts.glob + " returned no files" 
    sys.exit(1)
else:
  print >>sys.stderr, "Need to specify a glob"
  sys.exit(1)

#put the files into the coinc data structure
coincs = skylocutils.Coincidences(files,opts.input_type)

dt = []
dD = []

for coinc in coincs:
  inj_pt = (coinc.latitude_inj,coinc.longitude_inj)
  
  if len(coinc.ifo_list) < 3:
    continue
  else:
    if opts.snr_threshold:
      rhosquared = 0.0
      for ifo in coinc.ifo_list:
        rhosquared += coinc.snr[ifo]*coinc.snr[ifo]
      dtrss_inj = sqrt(rhosquared)*skylocutils.get_delta_t_rss(inj_pt,coinc,opts.reference_frequency)/10.0
    else:
      dtrss_inj = skylocutils.get_delta_t_rss(inj_pt,coinc,opts.reference_frequency)
    dDrss_inj = skylocutils.get_delta_D_rss(inj_pt,coinc)
    dt.append(dtrss_inj)
    dD.append(dDrss_inj)

#for the kdes
npts = opts.Npts

dtbw = compute_bw(dt,npts)
dtdat = np.array(dt,'float').reshape(len(dt),1)
pdtx = np.linspace(0.0,max(dt),npts)
pdty = [skylocutils.gaussian_kde(dtdat,xn,dtbw) for xn in pdtx]
pdtnorm = np.trapz(pdty,pdtx)
Ptx = pdtx
Pty = np.array([np.trapz(pdty[:i],Ptx[:i])/pdtnorm for i in np.arange(1,len(Ptx)+1)])

dDbw = compute_bw(dD,npts)
dDdat = np.array(dD,'float').reshape(len(dD),1)
pdDx = np.linspace(0.0,max(dD),npts)
pdDy = np.array([skylocutils.gaussian_kde(dDdat,xn,dDbw) for xn in pdDx])


rankings = {}
rankings['dt'] = skylocutils.Ranking(pdtx,pdty)
rankings['Pdt'] = skylocutils.Ranking(Ptx,Pty)
rankings['dD'] = skylocutils.Ranking(pdDx,pdDy)
rankings['snr_threshold'] = opts.snr_threshold
if opts.reference_frequency:
  rankings['ref_freq'] = opts.reference_frequency
else:
  rankings['ref_freq'] = None

Ls = [rankings['dt'].get_rank(t)*rankings['dD'].get_rank(D) for t,D in zip(dt,dD)]
Lbw = compute_bw(Ls,npts)
Ldat = np.array(Ls,'float').reshape(len(Ls),1)
#first construct the pdf for the values of L
#denote the pdf by f
fLx = np.linspace(0.0,max(Ls),npts)
fLy = np.array([skylocutils.gaussian_kde(Ldat,xn,Lbw) for xn in fLx])
#make sure the probabilities end at 1
fLnorm = np.trapz(fLy,fLx)
#now integrate to find the probability
Px = fLx
Py = np.array([np.trapz(fLy[:i],Px[:i])/fLnorm for i in np.arange(1,len(Px)+1)])

rankings['P'] = skylocutils.Ranking(Px,Py)

f = open(opts.output_file,'w')
cPickle.dump(rankings,f,protocol=2)
f.close()


if opts.plot_pdfs:
  
  #compute bin widths for the histograms
  dtbins = compute_nbins(dt,len(dt))
  dDbins = compute_nbins(dD,len(dD))
  Lbins = compute_nbins(Ls,len(Ls))
  
  from pylab import *

  figure()
  hist(dt,bins=dtbins,normed=1)
  plot(pdtx,pdty,'ro')
  title(r'p($\Delta t_{\rm rss}|\alpha,\delta)$', fontsize=20)
  savefig('pdt.png')

  figure()
  hist(dD,bins=dDbins,normed=1)
  plot(pdDx,pdDy,'ro')
  title(r'p($\Delta \tilde{D}_{\rm rss}|\alpha,\delta)$', fontsize=20)
  savefig('pdD.png')

  figure()
  hist(dt,bins=dtbins,normed=1,cumulative=1)
  plot(Ptx,Pty,'ro')
  title(r'Cumulative distribution of $p(\Delta t_{\rm rss}|\alpha,\delta)$', fontsize=20)
  savefig('cumpdt.png')

  figure()
  hist(Ls,bins=Lbins,normed=1,cumulative=1)
  plot(Px,Py,'ro')
  title(r'Cumulative distribution of $p(\Delta t_{\rm rss}|\alpha,\delta)p(\Delta \tilde{D}_{\rm rss}|\alpha\delta)$', fontsize=20)
  savefig('PdtdD.png')
