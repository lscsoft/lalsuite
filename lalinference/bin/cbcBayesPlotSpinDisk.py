# Ben Farr 2016

import os
from matplotlib import (use, rc)
use('agg')

rc('text', usetex=True)
rc('font', family='lmodern')

small, big = 10, 12
rc('axes', labelsize=big)
rc('font',size=big)
rc('legend', fontsize=small)
rc('xtick', labelsize=small)
rc('ytick', labelsize=small)

from lalinference.plot import make_disk_plot

USAGE='What this does is : '

if __name__=='__main__':
  import sys
  from optparse import OptionParser
  from lalinference import bayespputils as bppu

  parser=OptionParser(USAGE)
  parser.add_option("-o","--outpath", dest="outpath",default=None, help="make page and plots in DIR", metavar="DIR")
  parser.add_option("-d","--data",dest="data",help="Posteriors samples file (must be in common format)")
  (opts,args)=parser.parse_args()

  if opts.outpath is None:
    opts.outpath=os.getcwd()
  if not os.path.isfile(opts.data):
    print("Cannot find posterior file %s\n"%opts.data)
    sys.exit(1)
  else:
    peparser=bppu.PEOutputParser('common')
    commonResultsObj=peparser.parse(open(opts.data,'r'),info=[None,None])
    ps,samps = commonResultsObj
    pos = bppu.Posterior(commonResultsObj)
    make_disk_plot(pos,opts.outpath)
