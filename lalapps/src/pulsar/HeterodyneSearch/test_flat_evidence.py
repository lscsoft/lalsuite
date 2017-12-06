#!/usr/bin/env python

"""
A script to run lalapps_pulsar_parameter_estimation_nested with increasing amplitude prior
ranges. For a range of h0 ranges the nested sampling will be run and the odds ratio extracted.
It will calculate the value of the log(odds ratio)-log(prior) and check that it is roughly
flat as the prior range increases.
"""

import os
import glob
import sys
import numpy as np
import subprocess as sp
import scipy.stats as ss
import matplotlib.pyplot as pl

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
ds = np.zeros((1440,3))
ds[:,0] = np.linspace(900000000., 900000000.+86400.-60., 1440) # time stamps
ds[:,-2:] = 1.e-24*np.random.randn(1440,2)

# output data file
np.savetxt(datafile, ds, fmt='%.12e');

# range of upper limits on h0 in prior file
h0uls = [1e-22, 1e-21, 1e-20, 1e-19, 1e-18]

# some default inputs
dets='H1'
Nlive='1000'
Nmcmcinitial='200'
outfile='test.out'
#uniformprop='--uniformprop 0'
uniformprop=''

Ntests = 10 # number of times to run nested sampling for each h0 value to get average

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
    commandline="\
%s --detectors %s --par-file %s --input-files %s --outfile %s \
--prior-file %s --Nlive %s --Nmcmcinitial %s %s" \
% (execu, dets, parf, datafile, outfile, priorf, Nlive, Nmcmcinitial, uniformprop)

    sp.check_call(commandline, shell=True)

    # get odds ratio
    f = open(outfile+'_B.txt', 'r')
    hodds.append(float(f.readlines()[0].split()[0]))
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
fig, ax = pl.subplots(1, 1)
ax.errorbar(-np.array(logpriors), odds_prior, yerr=std_odds_prior, fmt='o')
ax.set_xlabel('log(prior volume)')
ax.set_ylabel('log(odds ratio)-log(prior)')
pl.show()

# clean up temporary files
os.remove(parf)
os.remove(priorf)
os.remove(datafile)
ofs = glob.glob(outfile+'*')
for fs in ofs:
  os.remove(fs)

sys.exit(0)
