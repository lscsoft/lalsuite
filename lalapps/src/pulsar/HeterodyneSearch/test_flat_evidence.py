#!/usr/bin/env python

"""
A script to run lalapps_pulsar_parameter_estimation_nested with increasing amplitude prior
ranges. For a range of h0 ranges the nested sampling will be run and the odds ratio extracted.
It will calculate the value of the log(odds ratio)-log(prior) and check that it is roughly
flat as the prior range increases.
"""

from __future__ import division

import os
import glob
import sys
import numpy as np
import subprocess as sp
import scipy.stats as ss
import matplotlib.pyplot as pl
import h5py

lalapps_root = os.environ['LALAPPS_PREFIX'] # install location for lalapps
execu = lalapps_root+'/bin/lalapps_pulsar_parameter_estimation_nested' # executable

# create files needed to run the code

# par file
parfile="\
PSRJ J0000+0000\n\
RAJ 00:00:00.0\n\
DECJ 00:00:00.0\n\
F0 100\n\
PEPOCH 54000"

parf = 'test.par'
f = open(parf, 'w')
f.write(parfile)
f.close()

# data file
datafile = 'data.txt.gz'
dlen = 1440 # number of data points
dt = 60     # number of seconds between each data point
startgps = 900000000. # start GPS time
endgps = startgps + dt*(dlen-1)
ds = np.zeros((dlen,3))
ds[:,0] = np.linspace(startgps, endgps, dlen) # time stamps
noisesigma = 1e-24
ds[:,-2:] = noisesigma*np.random.randn(dlen,2)

ulest = 10.8*np.sqrt(noisesigma**2/dlen)
print("Estimated upper limit is %.4e" % (ulest))

# output data file
np.savetxt(datafile, ds, fmt='%.12e');

# range of upper limits on h0 in prior file
h0uls = np.logspace(np.log10(5.*ulest), np.log10(500.*ulest), 4)

# some default inputs
dets='H1'
Nlive=1024
Nmcmcinitial=0
outfile='test.hdf'

# test two different proposals - the default proposal (which is currently --ensembleWalk 3 --uniformprop 1)
# against just using the ensemble walk proposal
proposals = ['', '--ensembleWalk 1 --uniformprop 0']
labels = ['Default', 'Walk']
pcolor = ['b', 'r']

Ntests = 10 # number of times to run nested sampling for each h0 value to get average

fig, ax = pl.subplots(1, 1)

for i, prop in enumerate(proposals):
  odds_prior = []
  std_odds_prior = []
  logpriors = []

  for h0ul in h0uls:
    # prior file
    priorfile="\
H0 uniform 0 %e\n\
PHI0 uniform 0 %f\n\
COSIOTA uniform -1 1\n\
PSI uniform 0 %f" % (h0ul, np.pi, np.pi/2.)

    priorf = 'test.prior'
    f = open(priorf, 'w')
    f.write(priorfile)
    f.close()

    logprior = -np.log(h0ul*np.pi*2.*(np.pi/2.))
    logpriors.append(logprior)

    hodds = []
    # run Ntests times to get average
    for j in range(Ntests):
      # run code
      commandline="%s --detectors %s --par-file %s --input-files %s --outfile %s --prior-file %s --Nlive %d --Nmcmcinitial %d %s" \
% (execu, dets, parf, datafile, outfile, priorf, Nlive, Nmcmcinitial, prop)

      sp.check_call(commandline, shell=True)

      # get odds ratio
      f = h5py.File(outfile, 'r')
      a = f['lalinference']['lalinference_nest']
      hodds.append(a.attrs['log_bayes_factor'])
      f.close()

    odds_prior.append(np.mean(hodds)-logprior)
    std_odds_prior.append(np.std(hodds))

  # use reduced chi-squared value to test for "flatness"
  ns = np.array(odds_prior)
  p = np.sum((ns-np.mean(ns))**2/np.array(std_odds_prior)**2)/float(len(h0uls))
  stdchi = np.sqrt(2.*float(len(h0uls)))/float(len(h0uls)) # standard deviation of chi-squared distribution
  nsigma = np.abs(p-1.)/stdchi

  print "Reduced chi-squared test for linear relation = %f" % (p)

  if nsigma > 2.:
    print "This is potentially significantly (%f sigma) different from a flat line" % nsigma

  # plot figure
  ax.errorbar(-np.array(logpriors), odds_prior, yerr=std_odds_prior, fmt='o', label=labels[i], color=pcolor[i])
  ax.set_xlabel('log(prior volume)')
  ax.set_ylabel('log(odds ratio)-log(prior)')

pl.legend()
pl.show()

# clean up temporary files
os.remove(parf)
os.remove(priorf)
os.remove(datafile)
ofs = glob.glob(outfile+'*')
for fs in ofs:
  os.remove(fs)

sys.exit(0)
