from numpy import log1p, log, logaddexp, array, digitize, loadtxt, zeros, exp, hstack, vstack, reshape
from numpy.random import uniform

def compute_weights(data, Nlive):
    """Returns log_ev, log_wts for the log-likelihood samples in data,
    assumed to be a result of nested sampling with Nlive live points."""

    start_data=data[:-Nlive]
    end_data=data[-Nlive:]

    log_wts=zeros(data.shape[0])

    log_vol_factor=log1p(-1.0/Nlive)
    log_dvol = -1.0/Nlive

    log_vol = 0.0
    log_ev = -float('inf')
    for i,log_like in enumerate(start_data):
        # Volume associated with this likelihood = Vol/Nlive:
        log_this_vol=log_vol+log_dvol
        log_wts[i] = log_like+log_this_vol
        log_ev = logaddexp(log_ev, log_wts[i])
        log_vol += log_vol_factor

    avg_log_like_end = -float('inf')
    for i,log_l in enumerate(end_data):
        avg_log_like_end = logaddexp(avg_log_like_end, log_l)
    avg_log_like_end-=log(Nlive)

    # Each remaining live point contributes (Vol/Nlive)*like to
    # integral, but have posterior weights Vol relative to the other samples
    log_wts[-Nlive:] = log_vol+end_data

    log_ev = logaddexp(log_ev, avg_log_like_end + log_vol)

    log_wts -= log_ev

    return log_ev, log_wts

def draw_posterior(data, log_wts):
    """Draw points from the given data (of shape (Nsamples, Ndim))
    with associated log(weight) (of shape (Nsamples,)). Draws uniquely so
    there are no repeated samples"""
    maxWt=max(log_wts)
    normalised_wts=log_wts-maxWt
    selection=[n > log(uniform()) for n in normalised_wts]
    idx=filter(lambda i: selection[i], range(len(selection)))
    return data[idx,:]
    
def draw_posterior_many(datas, Nlives,logLcol=-1):
    """Draw samples from the posteriors represented by the
    (Nruns, Nsamples, Nparams)-shaped array datas, each sampled with
    the corresponding Nlive number of live points. Will draw without repetition,
    and weight according to the evidence in each input run"""

    # list of log_evidences, log_weights
    log_evs,log_wts=zip(*[compute_weights(data[:,logLcol],Nlive) for data,Nlive in zip(datas, Nlives)])

    log_total_evidence=reduce(logaddexp, log_evs)
    log_max_evidence=max(log_evs)
    #print 'evidences: %s'%(str(log_evs))
    fracs=[exp(log_ev-log_max_evidence) for log_ev in log_evs]
    #print 'Relative weights of input files: %s'%(str(fracs))
    Ns=[fracs[i]/len(datas[i]) for i in range(len(fracs))]
    Ntot=max(Ns)
    fracs=[n/Ntot for n in Ns]
    #print 'Relative weights of input files taking into account their length: %s'%(str(fracs))

    bigpos=[]
    posts=[draw_posterior(data,logwt) for (data,logwt,logZ) in zip(datas,log_wts,log_evs)]
    #print 'Expected number of samples from each input file %s'%(str([int(f*len(p)) for f,p in zip(fracs,posts)]))
    for post,frac in zip(posts,fracs):
    	for samp in post:
	    if(uniform()<frac):
	    	bigpos.append(samp)
    return bigpos
    
def draw_N_posterior(data,log_wts, N):
    """
    Draw N samples from the input data, weighted by log_wt.
    For large N there may be repeated samples
    """
    if(N==0): return []
    log_cumsums=zeros(log_wts.shape[0]+1)
    log_cumsums[0]=-float('inf')
    for i,log_wt in enumerate(log_wts):
        log_cumsums[i+1]=logaddexp(log_cumsums[i], log_wt)

    us=log(uniform(size=N))
    idxs=digitize(us, log_cumsums)

    return data[idxs-1, :]

def draw_N_posterior_many(datas, Nlives, Npost, logLcol=-1):
    """
    Draw Npost samples from the posteriors represented by the
    (Nruns, Nsamples, Nparams)-shaped array datas, each sampled with
    the corresponding number of live points Nlive. The returned number
    of samples may not be exactly Npost due to rounding
    """
    # get log_evidences, log_weights.
    log_evs,log_wts=zip(*[compute_weights(data[:,logLcol],Nlive) for data,Nlive in zip(datas, Nlives)])
    
    log_total_evidence=reduce(logaddexp, log_evs)
    Ns=[int(round(Npost*exp(log_ev-log_total_evidence))) for log_ev in log_evs]
    posts=[draw_N_posterior(data,logwt,N) for (data,logwt,N) in zip(datas,log_wts,Ns)]
    bigpos=[]
    for post in posts:
        for samp in post:
            bigpos.append(samp)
    return bigpos
