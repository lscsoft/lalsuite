import sys
import math
import os
from numpy import log,loadtxt,vstack,array

def logadd(a,b):
    if(a>b): (a,b)=(b,a)
    return (b+log(1+exp(a-b)))

def getBfile(datname):
    Bfile=datname+'_B.txt'
    print 'Looking for '+Bfile
    if os.access(Bfile,os.R_OK):
        outstat = loadtxt(Bfile)
        return outstat
    else:
        return None

def x2iota(chain):
    chainlen=size(chain)
    iota=[]
    for i in range(0,chainlen):
        if chain[i] < 0:
            new=chain[i]+math.pi
        else:
            new=chain[i]
        iota=iota + [new]

    return (iota)

def makelist(path):
    #datalist=[]
    #for i in range(0,N):
    #    datalistitem=path + '/outfile_' + str(i) + '.dat'
    #    datalist += [datalistitem]
    dirlist=os.listdir(path)
    datalist=[]
    for fname in dirlist:
        if 'outfile' in fname and 'B.txt' not in fname:
            fname=path + '/' + fname
            datalist.append(fname)
    return datalist

def loaddata(datalist):
    out = list(map(loadtxt,datalist))
    Bfiles = list(map(getBfile,datalist))
    return out,Bfiles

def prodnoise(B):
    """
    Calculates sum (logZnoise[i] for i!=j) for each j to get logZ to add to each of the input files
    """
    N=len(B)
    totalnoise=[]
    for i in range(0,N):
        tn=0
    for j in range(0,N):
        if i!=j:
            tn+=B[j,2]
        totalnoise.append(tn)
    return totalnoise

def weightsamp(d,Nlive):
    #creates an array of vectors containing weights for every sample in every segment
    total_weight=[]
    for outfile in d:
        N=len(outfile[:,0])
        N_weighted = N - Nlive
        segment_weight=[]
        for i in range(1,N_weighted + 1):
            logw = -(i)/Nlive
            segment_weight.append(logw)
        for i in range(N_weighted + 1,N + 1):
            logw = -N_weighted / Nlive
            segment_weight.append(logw)
        total_weight += segment_weight
    return total_weight

def nest2pos(samps,weights):
    randoms=rand(size(samps,0))
    wt=weights+samps[:,-1]
    maxwt=max(wt)
    posidx=find(wt>maxwt+log(randoms))
    pos=samps[posidx,:]
    return pos

def combine_evidence(data,x2iota,Nlive):

    nfiles=len(data)

    #load in seperate data files#
    #datalist = makelist(path)
    (d,Bfiles) = loaddata(data)
    Barray = reduce(lambda x,y: vstack([x,y]), Bfiles)

    #Calculate total Bayes factor#
    if len(data)>1:
        ZnoiseTotal=sum(Barray[:,2])-log(len(data))
        totalBayes= reduce(logadd,Barray[:,0])
        totalBayes= float(totalBayes) - log(len(data)) #divide by 60 because we used 60 priors
    else:
        totalBayes=Barray[0]
        ZnoiseTotal=Barray[2]
    print "Total Bayes Factor= %f" %totalBayes

    #Scale likelihoods for entire sample#
    #Make list of sum(noise evidence), excepting current noise evidence)
    if len(data)>1:
        totalnoise = prodnoise(Barray)
    else: totalnoise=array([0])
    # Add logZnoise for other files to likelihoods for each sample
    if not None in Bfiles:
        for (outfile,noise) in zip(d,totalnoise):
            outfile[:,-1]+=noise

    #Remapping Parameters#
    for outfile in d:
        outfile[:,0]=exp(outfile[:,0])
        outfile[:,4]=exp(outfile[:,4])
        if xflag:
            outfile[:,8]=x2iota(outfile[:,8])

    #Posterior Samples
    weights=weightsamp(d,Nlive)
    d_all = reduce(lambda x,y: vstack([x,y]), d)
    pos=nest2pos(d_all,weights)

    d_idx=argsort(d_all[:,-1])
    d_all=d_all[d_idx,:]

    return pos,d_all,totalBayes,ZnoiseTotal
