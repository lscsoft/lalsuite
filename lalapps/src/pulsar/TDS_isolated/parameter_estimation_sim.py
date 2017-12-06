#!/usr/bin/python

"""
  script to run the pulsar parameter estimation code many times by generating
  fake noise, choosing a pulsar at random to use, generating condor_submit
  files and running the code in grid-based, 4-parameter MCMC, and
  multi-parameter MCMC modes
"""

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import exceptions
import random
import math

count = 1000
i = 0

# get a list of the pulsars in the par file directory
parlist = os.listdir('/home/matthew/analyses/S5/pulsars/fine')
parlist.sort()

covlist = os.listdir('/home/matthew/analyses/S5/pulsars/covariance')
covlist.sort()

# set default values
iterations = 200000
burnin = 200000
temperature = 0.001
det = 'H1'

begin = 815184013 # start time of S5
end = 877824013 # end time of S5
numsteps = 1440
dt = 60 # interval between point in the fake fine heterodyne file
stddev = 1 # standard deviation of the fake data

# set random seed
random.seed()

# create three condor submit files - one for grid-based mode, one for 4
# parameter MCMC mode and on for multi-parameter MCMC mode

# subfiles
try:
  f1 = open('/home/matthew/test/grid_based.sub', 'w')
  f2 = open('/home/matthew/test/mcmc_four.sub', 'w')
  f3 = open('/home/matthew/test/mcmc_multi.sub', 'w')
except Exception, e:
  print >> sys.stderr, "Can't open subfiles!"
  sys.exit(1)

univ = 'universe = standard\n'
f1.write(univ)
f2.write(univ)
f3.write(univ)

execu = 'executable = \
/home/matthew/lscsoft/lalapps/src/pulsar/TDS_isolated/\
lalapps_pulsar_parameter_estimation\n'
f1.write(execu)
f2.write(execu)
f3.write(execu)

noti = 'notification = never\n'
f1.write(noti)
f2.write(noti)
f3.write(noti)

log = 'log = /local/matthew/logs/param_sim.dag.log\n'
f1.write(log)
f2.write(log)
f3.write(log)

err = 'error = /local/matthew/logs/param_sim-$(cluster).err\n'
f1.write(err)
f2.write(err)
f3.write(err)

out = 'output = /local/matthew/logs/param_sim-$(cluster).out\n'
f1.write(out)
f2.write(out)
f3.write(out)

args1 = 'arguments = --pulsar $(macropulsar) --detectors %s --par-file \
$(macropar) --input-dir /home/matthew/test --output-dir $(macroout) --maxh0 0'\
% det
f1.write(args1)

args2 = args1 + ' --mcmc --iterations %d --burn-in %d --temperature %f \
--h0-width 0' % (iterations, burnin, temperature)
f2.write(args2)

args3 = args2 + ' --earth-ephem \
/home/matthew/lscsoft/lal/packages/pulsar/test/earth05-09-DE405.dat --sun-ephem\
 /home/matthew/lscsoft/lal/packages/pulsar/test/sun05-09-DE405.dat --covariance\
 $(macrocov)'
f3.write(args3)

f1.write('\n')
f2.write('\n')
f3.write('\n')

queue = 'queue 1\n'
f1.write(queue)
f2.write(queue)
f3.write(queue)

f1.close()
f2.close()
f3.close()

# open dag file
fdag = open('/home/matthew/test/parameter_estimation_sim.dag', 'w')

while i < count:
  # create a fake pulsar name (leave RA part as 0000 and increment through dec)
  psr = 'J0000+' + '%04d' % i

  # generate a start time for the fake data
  start = math.floor(begin + (end-begin)*random.random())

  finefile = '/home/matthew/test/dataH1/finehet_%s_%s' % (psr, det)
  try:
    # open file for output
    f = open(finefile, 'w')
  except Exception, e:
    print >> sys.stderr, "Can't open file %s!" % finefile
    sys.exit(1)

  # generate real and imaginary noise and output to file
  j = 0
  while j < numsteps:
    finedat = str(start+dt*j) + ' ' + str(random.gauss(0, stddev)) + ' ' +\
str(random.gauss(0, stddev)) + '\n'

    f.write(finedat)
    j += 1

  f.close()

  # pick a pulsar at random from parlist to use as parameters
  numpars = len(parlist)
  parname = parlist[int(math.floor(numpars*random.random()))]
  parfile = '/home/matthew/analyses/S5/pulsars/fine/' + parname

  # check whether theres a covariance file for that pulsar
  # get .mat files associated with .par file
  matfile = None
  tempmatfile = re.sub(r'\.par',r'',parname) + '.mat'
  k = 0
  while k < len(covlist):
    if tempmatfile in covlist[k]:
      matfile = '/home/matthew/analyses/S5/pulsars/covariance/' + tempmatfile
      break

    k += 1

  if matfile == None:
    continue

  # write jobs to dag
  outfile = '/home/matthew/test/output'
  job1 = 'JOB %d grid_based.sub\nRETRY %d 0\nVARS %d ' % ((3*i+1), (3*i+1),\
(3*i+1))
  vars1 = 'macropulsar=\"%s\" macropar=\"%s\" macroout=\"%s\"\n' % (psr,\
parfile, outfile)
  fdag.write(job1+vars1)
  
  job2 = 'JOB %d mcmc_four.sub\nRETRY %d 0\nVARS %d ' % ((3*i+2), (3*i+2),\
(3*i+2))
  fdag.write(job2+vars1)

  outfile = '/home/matthew/test/output_mcmc'
  job3 = 'JOB %d mcmc_multi.sub\nRETRY %d 0\nVARS %d ' % ((3*i+3), (3*i+3),\
(3*i+3))
  vars3 = 'macropulsar=\"%s\" macropar=\"%s\" macrocov=\"%s\" \
macroout=\"%s\"\n' % (psr, parfile, matfile, outfile)
  fdag.write(job3+vars3)

  i += 1

fdag.close()

sys.exit(0)
