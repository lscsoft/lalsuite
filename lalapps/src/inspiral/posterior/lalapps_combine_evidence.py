import sys
import os

from optparse import OptionParser

from combine_evidence import combine_evidence

usage='''%prog [options] datafile1.dat [datafile2.dat ...]
%prog takes at least one nested sampling output file and outputs posterior samples,
nested samples, and a Bayes factor summary file. If more than one input files is
specified, the files are merged as if they have equal prior volume. The resulting
Bayes Factor is the mean of the individual bayes factors of the runs.
The output samples have their likelihoods reweighted so samples from different
input files can be treated equally.
'''

if __name__=='__main__':

    parser=OptionParser(usage)
    parser.add_option("-o","--outsamp",help="location of output file with nested samples", metavar="OUTFILE",default=None,type="string")
    parser.add_option("-p","--outpos",help="location of output file of posterior samples",metavar="POSFILE",default=None,type="string")
    parser.add_option("-B","--bayesfactor",help="Location of file to store global bayes factor",default=None,type="string")
    parser.add_option("-X","--x2iota", dest="xflag", action="store_true", help="If parameter x needs to be converted to iota")
    parser.add_option("-N","--Nlive", dest="Nlive", action="store",type="int",help="number of live points for each of the files")
    parser.add_option("-H","--headers",dest="headers",action="store",type="string",help="Text file with parameter names",default=None)
    (opts,args)=parser.parse_args()

    data=args

    if opts.headers is None:
        headerstr='mchirp\t eta\t time\t phi0\t dist\t ra\t dec\t psi\t iota\t logL'
    else:
        headerfile=open(opts.headers,'r')
        headerstr=headerfile.readline()
        headerfile.close()

    logLcolumn=headerstr.lower().split().index('logl')
    if not opts.outsamp and not opts.outpos and not opts.bayesfactor:
        print 'No output specified, nothing to do'
        sys.exit(0)

    pos,d_all,totalBayes,ZnoiseTotal=combine_evidence(data,opts.xflag,opts.Nlive,logLcolumn=logLcolumn)

    if opts.outpos is not None:
        posfilename=opts.outpos
        posfile=open(posfilename,'w')
        posfile.write(headerstr+'\n')
        for row in pos:
            for i in row:
                posfile.write('%f\t' %(i))
            posfile.write('\n')
        posfile.close()

    #write new output file(s)#
    if opts.outsamp is not None:
        outfilename=opts.outsamp
        outfilenew=open(outfilename,'w')
        for row in d_all:
            for i in row:
                outfilenew.write('%f\t' %(i))
            outfilenew.write('\n')
        outfilenew.close()
        outfileB=open(outfilename+'_B.txt','w')
        # Write out B Zsig Znoise maxL
        print >>outfileB,"%f %f %f %f\n"%(totalBayes, totalBayes+ZnoiseTotal, ZnoiseTotal, pos[-1,logLcolumn]-ZnoiseTotal)
        outfileB.close()

    if opts.bayesfactor is not None:
        bayesfile=open(opts.bayesfactor,'w')
        bayesfile.write('%f' %totalBayes)
        bayesfile.close()
