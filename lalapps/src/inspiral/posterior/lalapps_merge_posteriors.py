#!/usr/bin/env python

# resample posteriors (C) John Veitch, 2015

from numpy import genfromtext,vstack,log
from optparse import OptionParser
import sys
import os
from pylal import bayespputils as bppu

from lalapps.nest2pos import draw_posterior, draw_N_posterior

usage='''%prog [-o output.dat] [-N NPOS] -p pos1.dat -w weight1 [-p pos2.dat -w weight2 ...]
%prog takes a list of posterior files and weights and draws samples from the combined,
reweighted distribution
'''

def load_data(filename,header=None):
    peparser=bppu.PEOutputParser('common')
    commonObj=peparser.parse(open(filename,'r'),info=[header,None])
    pos=bppu.Posterior(commonObj)
    return pos

if __name__=='__main__':
    parser=OptionParser(usage)
    parser.add_option('-o','--output',action='store',type='string',default=None,help='output file',metavar='output.dat')
    parser.add_option('-p','--posterior',action='append',type='string',default=[],metavar='pos.dat',help='posterior input file')
    parser.add_option('-w','--weight',action='append',type='float',default=[],metavar='NUM',help='weight of an input file')
    parser.add_option('-N','--npos',action='store',default=None,metavar='NPOS',help='Optional number of posterior samples to draw')
    parser.add_option('--verbose',action='store_true',default=False,help='Prints additional information')
    (opts,args) = parser.parse_args()

    if len(opts.posterior)==0:
        sys.stderr.write('No input files given\n')
        sys.exit(1)
    if len(opts.weights) != len(opts.posterior):
        sys.stderr.write('Error: must specify same number of weights and posteriors\n')
        sys.exit(1)

    # Create posterior samples for each input file
    datas=map(load_data,opts.posterior)
    weights = [ log(x) for _ in d for d in datas]
    bigdata=vstack(datas)
    
    # Call reweighting function
    if opts.npos is not None:
        merged=draw_N_posterior(bigdata,weights,opts.npos,verbose=opts.verbose)
    else:
        merged=draw_posterior(bigdata,weights,verbose=opts.verbose)

    outObj=bppu.Posterior((datas[0].names,merged))

    # Write output file
    if opts.output is not None:
        outObj.write_to_file(opts.output)

    return

