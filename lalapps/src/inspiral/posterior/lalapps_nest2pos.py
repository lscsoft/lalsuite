from numpy import loadtxt,logaddexp,log
from optparse import OptionParser
import sys
import os
import gzip

from lalapps.nest2pos import draw_posterior_many, draw_N_posterior_many,compute_weights

usage='''%prog -N Nlive [-p posterior.dat] [-H header.txt] [--npos Npos] datafile1.dat [datafile2.dat ...]

%prog takes at least one nested sampling output file and outputs posterior samples.
\tIf more than one input file is specified, each file is converted, then posterior samples drawn
\taccording to the evidence of each. Will output to stdout if no -p option given.
\tIf the --npos option is used the algorithm will draw approximately that number of samples
\tfrom the posterior. This may give repeated samples in the output file. By default, the
\tnon-repeating algorithm is used, but that may not produce enough samples.
'''

def getBfile(Bfile):
    Bfile=Bfile
    if os.access(Bfile,os.R_OK):
        outstat = loadtxt(Bfile)
        return outstat
    else:
        return None

if __name__=='__main__':
  headerfilename=None
  parser=OptionParser(usage)
  parser.add_option('-N','--Nlive',action='store',type='int',dest='Nlive',default=None,help='Number of live points in each chain loaded',metavar='NUM')
  parser.add_option('-p','--pos',action='store',type='string',default=None,help='Output file for posterior samples',metavar='posterior.dat')
  parser.add_option('--npos',action='store',type='int',default=None,help='Draw a specific number of posteriors samples. May give repetitions. Disabled by default.',metavar='NUM')
  parser.add_option('-H','--headers',action='store',type='string',default=None,help='Header file explaining columns in data file',metavar='file.dat_params.txt')
  parser.add_option('-d','--prec',action='store',type='int',dest='prec',default=14,help='Number of decimal place required for output posterior samples. Default is 14.',metavar='NUM')
  parser.add_option('-z','--gzip',action="store_true",dest='gz',default=False,help='Gzip compress the output posterior samples (this will append .gz onto the posterior file). Default is no compression.')
  parser.add_option('-v','--verbose',action="store_true",default=False,help="Print some additional information")
  (opts,args)=parser.parse_args()

  datafiles=args
  if(len(datafiles)<1):
    print 'No input file specified, exiting'
    sys.exit(1)
  if opts.Nlive is None:
    print 'Must specify number of live points using the --Nlive option'
    sys.exit(1)

  if opts.headers is None:
    defaultparamfilename=args[0]+'_params.txt'
    if os.access(defaultparamfilename,os.R_OK):
        headerfilename=defaultparamfilename
  else:
    headerfilename=opts.headers
  if headerfilename is not None:
    if not os.access(headerfilename,os.R_OK):
      print 'Unable to open header file %s!'%(headerfilename)
      sys.exit(1)
    elif opts.verbose: print 'Reading headers from %s'%(headerfilename)
    headerfile=open(headerfilename,'r')
    headerstr=headerfile.readline()
    headerfile.close()
  else: # old inspnest defaults
    if opts.verbose:
      print 'No header file found, assuming inspnest default'
    headerstr='mchirp\t eta\t time\t phi0\t dist\t ra\t dec\t psi\t iota\t logL'
  headers=headerstr.lower().split()
  logLcol=headers.index('logl')

  if opts.npos is not None:
     sampler=lambda datas,Nlives,logLcol: draw_N_posterior_many(datas,Nlives,opts.npos,logLcol=logLcol)
  else:
     sampler=draw_posterior_many

  # set the output number format
  if opts.prec < 1 or opts.prec > 20: # precision is outside a reasonable range, so just set to default of 14
    prec = 14
  else:
    prec = opts.prec

  outformat = '%%.%de' % prec

  # Create posterior samples for each input file
  inarrays=map(loadtxt,datafiles)
  posterior=sampler(inarrays,[int(opts.Nlive) for d in datafiles],logLcol=logLcol)
  log_evs,log_wts=zip(*[compute_weights(data[:,logLcol],opts.Nlive) for data in inarrays])
  if opts.verbose:
    print 'Log evidences from input files: %s'%(str(log_evs))
  if opts.pos is not None:
    if opts.gz:
      posfile=gzip.open(opts.pos+'.gz','wb')
    else:
      posfile=open(opts.pos,'w')
  else:
    posfile=sys.stdout
  for h in headers:
    posfile.write('%s\t'%(h))
  posfile.write('\n')
  for row in posterior:
    for i in row:
      outval = outformat % i
      posfile.write('%s\t'%(outval))
    posfile.write('\n')
  if opts.pos is not None:
    posfile.close()

  # Write out an evidence file
  if opts.pos is not None:
    Bfile=opts.pos+'_B.txt'
    Bs=map(getBfile,[d.replace('.gz', '')+'_B.txt' for d in datafiles])
    if not any([b is None for b in Bs]):
      meanZ=reduce(logaddexp,log_evs)-log(len(log_evs))
      noiseZ=reduce(logaddexp,[b[2] for b in Bs]) - log(len(Bs))
      meanB=meanZ-noiseZ
      maxL=max([b[3] for b in Bs])
      outB=open(Bfile,'w')
      outB.write('%f %f %f %f\n'%(meanB,meanZ,noiseZ,maxL))
      outB.close()

