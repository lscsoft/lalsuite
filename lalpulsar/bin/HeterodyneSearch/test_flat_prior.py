"""
A script to run lalpulsar_parameter_estimation_nested with increasing amplitude prior
ranges. Samples supposedly generated from the prior will be output and checked using a KS test
to see if they match the expected flat prior across the range.
"""

import os
import sys
import numpy as np
import subprocess as sp
import scipy.stats as ss
import h5py

if os.environ['LALINFERENCE_ENABLED'] == 'false':
    print('Skipping test: requires LALInference')
    sys.exit(77)

exit_code = 0

execu = './lalpulsar_parameter_estimation_nested' # executable

# lalpulsar_parameter_estimation_nested runs much slower with memory debugging
os.environ['LAL_DEBUG_LEVEL'] = os.environ['LAL_DEBUG_LEVEL'].replace('memdbg', '')
print("Modified LAL_DEBUG_LEVEL='%s'" % os.environ['LAL_DEBUG_LEVEL'])

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
h0uls = [1e-22, 1e-21, 1e-20, 1e-19, 1e-18, 1e-17]

# some default inputs
dets='H1'
Nlive='5000'
Nmcmcinitial='0'
outfile='test.hdf'
outfile_SNR='test_SNR'
outfile_Znoise='test_Znoise'
priorsamples=Nlive

for h, h0ul in enumerate(h0uls):
    print("--- h0=%i/%i ---" % (h+1, len(h0uls)), flush=True)

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

    # run code
    commandline="\
%s --detectors %s --par-file %s --input-files %s --outfile %s --prior-file %s --Nlive %s --Nmcmcinitial %s --sampleprior %s" \
% (execu, dets, parf, datafile, outfile, priorf, Nlive, Nmcmcinitial, priorsamples)

    sp.check_call(commandline, shell=True)

    # read in prior samples
    f = h5py.File(outfile, 'r')
    a = f['lalinference']
    h0samps = a['lalinference_nest']['nested_samples']['H0'][:]

    # get normed histogram of samples
    [n, nedges] = np.histogram(h0samps, bins=20, range=(0., h0ul), density=True)
    nc = np.cumsum(n)*(nedges[1]-nedges[0])

    stat, p = ss.kstest(nc, 'uniform')

    print("K-S test p-value for upper range of %e = %f" % (h0ul, p))

    if p < 0.005:
        print("There might be a problem for this prior distribution")
        try:
            import matplotlib as mpl
            mpl.use('Agg')
            import matplotlib.pyplot as pl
        except ModuleNotFoundError:
            print("matplotlib unavailable; skipping plot")
            exit_code = 1
            break
        fig, ax = pl.subplots(1, 1)
        ax.hist(h0samps, bins=20, density=True, cumulative=True, histtype='stepfilled', alpha=0.2)
        ax.plot([0., h0ul], [0., 1], 'k--')
        ax.set_xlim((0., h0ul))
        ax.set_ylim((0., 1.))
        ax.set_xlabel('h_0')
        ax.set_ylabel('Cumulative probability')
        fig.savefig('h0samps.png')
        print("Saved plot to 'h0samps.png'")
        exit_code = 1
        break

    # clean up per-run temporary files
    for fs in (outfile, outfile_SNR, outfile_Znoise):
        os.remove(fs)

# clean up temporary files
for fs in (priorf, parf, datafile):
    os.remove(fs)

sys.exit(exit_code)
