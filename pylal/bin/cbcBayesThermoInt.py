#!/usr/bin/env python

from optparse import OptionParser
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pp
import numpy as np
import scipy.integrate as si

matplotlib.rcParams['text.usetex'] = True

def extract_temp(filename):
    """Extracts the PTMCMC temperature from the header lines of the
    given file."""
    with open(filename, 'r') as inp:
        for line in inp:
            line=line.lstrip().split()
            try:
                i=line.index('Tchain')
                #####################################################
                # WARNING: hardcoded off-by-one adjustment to account
                # for 'null likelihood' header name that splits into
                # two list elements
                #####################################################
                return float(inp.next().split()[i-1])
            except ValueError:
                pass

        raise ValueError('extract_temp: did not find header line with \'Tchain\'')

def get_mean_logl(filename):
    """Returns the mean value of log(L) from the given filename,
    excluding the first 50% of samples as burn-in."""
    with open(filename, 'r') as inp:
        skip=0
        for line in inp:
            line=line.lstrip().split()
            skip+=1
            try:
                line.index('cycle')
                break
            except ValueError:
                pass

    col = line.index('logl')

    data=np.loadtxt(filename, skiprows=skip, usecols=(col,))
    N=data.shape[0]

    return np.mean(data[N/2:])

if __name__=='__main__':
    # Custom usage and help message
    usage = """%s [-h] [--plotfile FILE] [--evfile FILE] OUTPUT_FILE [OUTPUT_FILE ...]
    
Thermodynamic integration on PTMCMC samples.

positional arguments:
  OUTPUT_FILE      PTMCMC output files""" % (os.path.basename(sys.argv[0]))

    parser = OptionParser(usage=usage)
    parser.add_option('--plotfile', metavar='FILE', default='evidence-integrand.pdf', help='plot output file')
    parser.add_option('--evfile', metavar='FILE', default='evidence.dat', help='evidence output file')

    (options, args) = parser.parse_args()

    # Make positional arguments required
    if len(args)==0:
        parser.error('Positional filename argument(s) are required')

    betas = np.array([1.0/extract_temp(f) for f in args])
    logls = np.array([get_mean_logl(f) for f in args])

    inds = np.argsort(betas)[::-1]

    betas = betas[inds]
    logls = logls[inds]

    # Now extend to infinite temperature by copying the last <log(L)>.
    # This works as long as the chains have extended to high enough
    # temperature to sample the prior.
    ebetas = np.concatenate((betas, [0.0]))
    elogls = np.concatenate((logls, [logls[-1]]))

    ebetas2 = np.concatenate((betas[::2], [0.0]))
    elogls2 = np.concatenate((logls[::2], [logls[::2][-1]]))

    evidence = -si.trapz(elogls, ebetas)
    evidence2 = -si.trapz(elogls2, ebetas2)

    devidence = np.abs(evidence - evidence2)

    pp.plot(betas, betas*logls)
    pp.xscale('log')
    pp.xlabel(r'$\beta$')
    pp.ylabel(r'$\beta \left\langle \ln \mathcal{L} \right\rangle$')
    pp.savefig(options.plotfile)

    with open(options.evfile, 'w') as out:
        out.write('# ln-Z delta-ln-Z\n')
        out.write(str(evidence) + ' ' + str(devidence))
    
    
    
