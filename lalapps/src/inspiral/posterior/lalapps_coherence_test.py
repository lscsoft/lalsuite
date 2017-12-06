'''
This program simply computes the bayes factor between the coherent and incoherent models
given as input the B files for the run.
'''

from optparse import OptionParser
import sys
import os
import matplotlib
matplotlib.use("Agg")
from pylab import loadtxt
from math import log,exp

usage='''%prog [options] coherent_B.txt incoherent1_B.txt incoherent2_B.txt ... incoherentN_B.txt
%prog takes a list of arguments which specify FIRST the coherent B.txt file, followed by all the
B.txt files from the individual runs on each IFO. The coherent file must be specified first.

The code will then perform either the test of coherent vs incoherent signal models, i.e.
B_{CI} = Z_co / Z_1 * Z_2 * Z_3

or the test of coherent vs ( incoherent signal or noise ) models, i.e.
B_{CIN} / Z_co / (Z_1 * Z_2 *Z_3 + Z_noise )

Output will be written to stdout, and to a text file specified by the --outfile option.

'''

parser=OptionParser(usage)
parser.add_option("-c","--coherent-incoherent",action="store_true",dest="ci",help="Perform coherence test of Coherent signal vs incoherent signal",default=False )
parser.add_option("-n","--coherent-incoherent-noise", dest="cin", action="store_true",  help="Perform coherence test of coherent signal vs incoherent signal or noise", default=False )
parser.add_option("-o","--outfile",default=None, help="Optional output file name", metavar="OUTFILE.txt",type="string",action="store")

(opts,args)=parser.parse_args()

# Sanity check input arguments
if len(args)==0:
	print 'No input files specified'
	sys.exit(1)

if((not opts.ci and not opts.cin) or (opts.ci and opts.cin)):
	print 'Please specify either -c or -n. See help for details.'
	sys.exit(1)

cofile=args[0]
incofiles=args[1:]

if len(incofiles)<2:
	print 'Cannot perform coherence test on less than 2 IFOs'
	sys.exit(0)

if len(incofiles)==0:
	print 'Error: You must specify incoherent input files'
	sys.exit(1)

def logadd(a,b):
    if(a>b): (a,b)=(b,a)
    return (b+log(1+exp(a-b)))

def getBfile(Bfile):
    Bfile=Bfile
    print 'Looking for '+Bfile
    if os.access(Bfile,os.R_OK):
        outstat = loadtxt(Bfile)
        return outstat
    else:
        return None

co=getBfile(cofile)
incos=map(getBfile, incofiles)

Zco=co[1]
Zinco=0
for inco in incos:
	Zinco+=inco[1]
Znoise=co[2]

if(opts.outfile is not None):
	out=open(opts.outfile,'w')

if(opts.ci):
	B=Zco-Zinco
if(opts.cin):
	B=Zco-logadd(Zinco,Znoise)

print '%.3f'%(B)

if(opts.outfile is not None ):
	print >>out,'%.3f'%(B)
	out.close()

