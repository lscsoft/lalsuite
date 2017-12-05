#!/bin/python

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 22:20:43 2014

@author: Michalis Agathos, Tjonnie Li
"""

__author__ = "Michalis Agathos, Tjonnie Li"
__credits__ = ["Michalis Agathos", "Tjonnie Li"]
__maintainer__ = "Michalis Agathos, Tjonnie Li"
__email__ = "michalis.agathos@ligo.org, tjonnie.li@ligo.org"
__status__ = "Production"

###########################################
#
#  Define hardcoded
#
###########################################

#CUTS
Bcut = 3.0
NetSNRcut = 8.0
HNetSNRcut = 50.0
exclude = []


# PLOTTING STYLES
typemarkers = {} 
typecolors = {}

###########################################
#
#  Import modules
#
###########################################

import os
#import argparse
import optparse as op
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#from matplotlib import animation
from pylab import *

from scipy import stats
from scipy import integrate
from cPickle import dump, load
from ConfigParser import SafeConfigParser

#from sklearn.neighbors import KernelDensity
#from scipy.stats import gaussian_kde
#from statsmodels.nonparametric.kde import KDEUnivariate
#from statsmodels.nonparametric.kernel_density import KDEMultivariate


###########################################
#
#  Define preferences (plotting etc)
#
###########################################


plot_params = {'backend': 'png',
		'axes.labelsize': 16, 
		'axes.titlesize': 24, 
		'text.fontsize': 14, 
		'title.fontsize': 18, 
		'legend.fontsize': 12, 
		'xtick.labelsize': 14, 
		'ytick.labelsize': 14, 
		'axes.grid' : True,
		'text.usetex': True,
		'lines.markersize' : 5}
          
rcParams.update(plot_params)

###########################################
#
#  Class definitions
#
###########################################

class Parameter(object):
  '''Simple parameter class'''
  def __init__(self, name, prior_shape, min, max, latex_name="", latex_unit="", sigma=1., gaussmean=0, injected=None, nbins=1000, xfm=None):
    '''Generates the parameter.'''
    self.name = name
    self.nbins = nbins
    self.prior_shape = prior_shape
    self.gaussmean = gaussmean
    self.sigma = sigma
    self.inj = injected
    self.min = min
    self.max = max
    self.latex_name = latex_name
    self.latex_unit = latex_unit
    self.xfm = xfm
    self.setPrior()

  def setPrior(self):
    '''Generates prior distribution'''
    x = linspace(self.min, self.max, self.nbins+1)
    self.domain = (x[1:] + x[:-1])*0.5
    if self.prior_shape == "gaussian":
      prior = stats.norm.pdf(self.domain, self.gaussmean, self.sigma)
    elif self.prior_shape == "flat":
      prior = 1.0/(self.max-self.min)*ones(self.nbins)
    elif self.prior_shape == "cauchy": #FIXME: NORMALIZATION???
      prior = 1.0/(1 + pow(self.domain ,2.))
    self.prior = prior

  def setPosterior(self, data, weight=None):
    '''Read posterior from 1-d array'''
    self.postpoints = array(data)
    if weight is None:
      self.posterior = histogram(self.postpoints, self.nbins, (self.min, self.max))[0]
      postrestr = self.postpoints[ where( logical_and(self.min < self.postpoints, self.postpoints < self.max) )[0] ]
#      self.kde = stats.gaussian_kde(mirrorEdges(self.postpoints, min(self.postpoints), max(self.postpoints)))
      self.kde = stats.gaussian_kde(mirrorEdges(postrestr, self.min, self.max))
    else:
      print "Not ready to handle weights yet"
      exit(-1)
      self.posterior = histogram(self.postpoints, self.nbins, (self.min, self.max), weights=weight)[0]
      self.kde = stats.gaussian_kde(mirrorEdges(self.postpoints, self.min, self.max)) ###FIXME EEE
            
  def plotPriorPDF(self):
    return

class GlobalParam(Parameter): 
  '''A class for parameters that are common across sources. It inherits from Parameter'''
  def __init__(self, name, min, max, prior_shape='flat', sourcelist=None, latex_name="", latex_unit="", injected=None, color='k', nbins=1000, xfm=None):
    super(GlobalParam, self).__init__(name, prior_shape, min, max, latex_name, latex_unit, injected=injected, nbins=nbins, xfm=xfm)
    if sourcelist is None:
      self.sourcelist = []
    else:
      self.sourcelist = sourcelist
    self.postpoints = []
    self.posteriors = []
    self.kdes = []
    self.medians = []
    self.confs = []
    self.comb_post = []
    self.comb_kde = []
    self.color = color
    
  def addSource(self, source):
    if self.name in source.recparams.keys():
      self.sourcelist.append(source)
    else:
      print "This source does not have the variable " + self.name + " as recovered parameters."

  def setPosterior(self):
    '''Read posterior from list of sources'''
    self.posteriors = []
    self.kdes = []
    for s in self.sourcelist:
      index = s.recparams[self.name]
#      print s, index, len(s.posteriors[index]), min(s.posteriors[index]), max(s.posteriors[index]), self.min, self.max
      data = s.posteriors[index]
      postpoints = array(data)
      self.posteriors.append(histogram(postpoints, self.nbins, (self.min, self.max))[0])
      postrestr = postpoints[where( logical_and(postpoints > self.min, postpoints < self.max) )[0]]
#      self.kde = stats.gaussian_kde(mirrorEdges(self.postpoints, min(self.postpoints), max(self.postpoints)))
      self.kdes.append(stats.gaussian_kde(mirrorEdges(postrestr, self.min, self.max)))
#      self.kdes.append(stats.gaussian_kde(postpoints))
#      self.kdes.append(stats.gaussian_kde(mirrorEdges(postpoints, self.min, self.max)))

  def combinePosteriors(self):
    '''Combine posteriors from all the sources in the list'''
    print "Combining posteriors for " + self.name + " from " + str(len(self.sourcelist)) + " sources..."
    if len(self.posteriors) == 0:    
      self.setPosterior()
    #print len(self.posteriors)
    self.comb_post = [self.prior]
    self.comb_kde = [self.prior]
    self.comb_confs = []
    self.comb_kconfs = []

    cp = self.comb_post[0]
    ck = self.comb_kde[0]
    for p in self.posteriors:
      cp = cp*p/self.prior
      cp = cp/integrate.simps(cp, self.domain)
      self.comb_post.append(cp)
      self.comb_confs.append(getConfidenceIntervals(cp, self.domain, [0,68,95,99]))
    for k in self.kdes:
      ck = ck*k.evaluate(self.domain)
      ck = ck/integrate.simps(ck, self.domain)
      self.comb_kde.append(ck)
      self.comb_kconfs.append(getConfidenceIntervals(ck, self.domain, [0,68,95,99]))
  
  def replotSummary(self, sax):
    '''Plots the median and confidence interval for combined posteriors over sources'''
    sax.set_ylabel(self.latex_name + '$[$' + self.latex_unit + '$]$') #FIXME: transfer to parent
    sax.set_ylim(self.min, self.max)
    medians = []
    err95 = []
        
    '''Loop through sources'''
    for i in arange(len(self.sourcelist)):
      #s = self.sourcelist[i]
      #cc = self.comb_confs[i]
      cc = self.comb_kconfs[i]
      median = cc[0][1][0]
      medians.append(median)
      err95.append(abs(median - cc[2][1]))
      #ax_sum.annotate(str(s.snr)+"\n"+str("{0:0.1f}".format(s.bfac)), xy=(i+1,float(median)), xytext=(i+0.9,float(median)+0.4*((-1)**i)), size=5)

    err95 = array(err95).transpose()
    sax.errorbar(arange(len(medians))+1, medians, yerr=err95, fmt='o', mfc=self.color, ecolor=self.color, mec='k', label='95\% CI '+ self.sourcelist[0].batch.label)        
  
  def plotCombinedPDFs(self, outdir=None, sourceidx=None):
    '''Plots individual PDFs together with combined PDFs'''
    self.figs = {}
    x = self.domain

    if sourceidx is None:
      slist = self.sourcelist
    else:
      slist = list(sourceidx)
    
    '''Loop through sources'''
    for i in arange(len(slist)):

      '''Set up posterior plot for each source'''      
      #p = self.posteriors[i]
      k = self.kdes[i].evaluate(self.domain)
      k /= integrate.simps(k, self.domain)
      #cp = self.comb_post[i+1]
      ck = self.comb_kde[i+1]
      s = self.sourcelist[i]
      #pp = self.postpoints[i]
      pp = s.posteriors[s.recparams[self.name]]
      fig = plt.figure()
      ax = fig.add_subplot(111)
  
      '''Setup plot axes'''
      ax.set_title("Event "+str(s.id)+" of batch " + str(s.seed))
      ax.text(0.1,0.9,r" $\rho_{\rm net}="+str(s.SNR)+"$",fontsize=24, transform = ax.transAxes)
      ax.grid(color='grey', linestyle='--', linewidth=0.5,alpha=0.5)
      ax.set_xlim(self.min, self.max)
      ax.set_ylim(bottom=0.0)
      ax.set_xlabel(self.latex_name + '$[$' + self.latex_unit + '$]$')
      ax.set_ylabel("$p($" + self.latex_name + "$|d, I)$")
        
      '''Plot injected values as vertical lines'''
      if self.inj is not None:
        ax.axvline(x=self.inj, color="r", linestyle="--")
      elif self.name in s.injparams.keys():
        ax.axvline(x=s.injparams[self.name], color="r", linestyle="..")
        
      '''Plot things'''
      #ax.plot(x, p, color='g', label='single source hist')
      ax.plot(x, k, color='g', label='single source kde')
      ax.fill_between(x, zeros(size(x)), k, facecolor="g",alpha=0.3)
      #ax.plot(x, cp, color='k', label=str(i+1) + ' combined sources hist')
      ax.plot(x, ck, color='k', label=str(i+1) + ' combined sources kde')
      ax.fill_between(x, zeros(size(x)), ck, facecolor="k",alpha=0.3)
      ax.hist(pp, linspace(self.min, self.max, 30), normed=True, color='blue', histtype='stepfilled', alpha=0.3, label='single source hist')
      if eplot:
        ax.eventplot(pp[list(set(randint(0, len(pp), 100)))], colors=array([[1,0,0]]), linelengths = 0.01, lineoffsets=0.005)
      
      ax.legend()

      '''Save figure to file'''
      outname = os.path.join(outdir, s.batch.label, self.name + '-' + s.batch.label + '-' + str(i+1) + '-post.png')
      if outdir is None:
        print "Saving figure in memory"
        self.figs[outname] = fig
      else:
        print 'Saving figure to file ' + outname
        fig.savefig(outname, bbox_inches='tight')

  def animatePDFs(self, outdir=None):
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)
    return anim


'''
class Parameter2D(object):
class GlobalParam2D(Parameter2D):
ALMOST READY
'''
 

class Batch: 
  '''Generates a batch given at least a cluster name, a location, a seed, nlive and a label.'''
  def __init__(self, cluster, basedir, seed, nlive, label, xmlfile=None, description="", rundir=None, logdir=None, snrdir=None, postdir=None):
      self.cluster = cluster
      self.basedir = basedir
      self.seed = seed
      self.nlive = nlive
      self.label = label
      self.xml = xmlfile
      self.sources = []
      self.desc = description
      self.rundir = rundir
      self.logdir = logdir
      self.snrdir = snrdir
      if rundir is None:
          self.rundir = os.path.join(basedir, 'engine')
      if logdir is None:
          self.logdir = os.path.join(basedir, 'log')
      if snrdir is None:
          self.snrdir = os.path.join(basedir, 'SNR')
      if postdir is None:
          self.postdir = os.path.join(basedir, 'posterior_samples')

  def populateSources(self):
      self.sources = []
      paramfiles = matchFiles(self.rundir, ['lalinferencenest-', '_params.txt'])
      for f in paramfiles:
          chainfile = f.rpartition('_')[0] 
          (n, gps, ifos) = decChainFile(chainfile)
          self.sources.append(Source(self, n, gps, self.nlive, xmlfile=self.xml, chainfile=os.path.join(self.rundir, chainfile), IFOs=ifos))
      print "Added " + str(len(paramfiles)) + " sources for seed " + str(self.seed) + " from " + self.rundir + " ."

      return


class Source:
    def __init__(self, batch, eventno, gpstime, nlive, xmlfile=None, chainfile=None, posteriorfile=None, Bfile=None, paramFile=None, injFile=None, IFOs=['H1','L1','V1']):
      self.batch = batch
      self.seed = batch.seed
      self.id = eventno
      self.gpstime = gpstime
      self.xml = xmlfile
      self.nlive = nlive
      self.injparams = {}
      self.recparams = {}
      self.figs = {}
      self.posteriors = []
      self.postweights = {}
      self.IFOs = IFOs
      self.chainfile = chainfile
      if chainfile is None:
        chainfname = 'lalinferencenest-' + str(self.id + 1) + '-' + ''.join(self.IFOs) + '-' + str(self.gpstime) + '.0-' + str(self.id) + '.dat'
        self.chainfile = os.path.join(self.batch.rundir, chainfname)
      self.postfile = posteriorfile
      if posteriorfile is None:
        postfname = 'posterior_' + ''.join(self.IFOs) + '_' + str(self.gpstime) + '-' + str(self.id + 1) + '.dat'
        self.postfile = os.path.join(self.batch.postdir, postfname)
      self.isdone = False
      self.Bfac = None
      if Bfile is None:
        self.Bfile = self.chainfile + '_B.txt'
      self.paramfile = paramFile
      if paramFile is None:
        self.paramfile = self.chainfile + '_params.txt'
      self.injfile = injFile
      if injFile is None:
        self.injfile = self.chainfile + '.injection'
      self.snrfile = os.path.join(self.batch.snrdir, 'snr_' + ''.join(self.IFOs) + '_' + str(self.gpstime) + '.0.dat')

    def isDone(self):
      '''Check if run is finished for this source'''
      self.isdone = os.path.isfile(self.Bfile)
      return self.isdone
    
    def setBfac(self):
      if self.isDone():
        d = genfromtxt(self.Bfile)
        self.Bfac = d[0]
      else:
        self.Bfac = None

    def isFound(self, threshold=0.0):
      '''Check if run makes the Bayes factor threshold'''
      self.setBfac()
      if self.Bfac is None:
        return False
      else:
        return (self.Bfac > threshold)
      
    def SNRpass(self, threshold = 8.0):
      '''Check if run makes the SNR threshold'''
      if os.path.isfile(self.snrfile):
        self.SNR = genfromtxt(self.snrfile)[3,1]
        return (self.SNR >= threshold)
      else:
        return False

    def getParams(self):
      '''Populate recovery and injection parameters from files'''
      self.recparams = {}
      if os.path.isfile(self.paramfile):
        f = open(self.paramfile, 'r')
        l = f.readline().strip()
        f.close()
        parlist = l.split()
        self.recparams = dict(zip(parlist, arange(len(parlist))))
      self.injparams = {}
      if os.path.isfile(self.injfile):
        injf = open(self.injfile, 'r')
        nl = injf.readline().strip()
        vl = injf.readline().strip()
        injf.close()
        injnames = nl.split()
        injvalues = vl.split()
        self.injparams = dict(zip(injnames, injvalues))

    def getXfmParams(self, xfm):
      '''Add names and posteriors of transformed parameters'''
      if self.recparams == {}:
        self.getParams()
      if self.posteriors == []:
        self.setPosteriors()
      if set(xfm.varnames) < set(self.recparams.keys()):
        vardict = {}
        for vn in xfm.varnames:
          vardict[vn] = self.posteriors[self.recparams[vn]]
        postdict, weight = xfm.xfmFunc(vardict)
        for name in (set(xfm.names) - set(self.recparams.keys())):
          self.recparams[name] = len(self.recparams)
          self.posteriors.append(postdict[name])
          self.postweights[name] = weight
      else:
        print "Could not add params. Variables are not in the recovered params!"
      

    def setPosteriors(self, seed=1234):
      '''Generate posteriors or read from file'''
      if self.recparams == {}:
        self.getParams()
      if os.path.isfile(self.postfile):
        p = genfromtxt(self.postfile)[1:]
        np = len(p)
        print 'Found posterior file: ' + self.postfile
      elif os.path.isfile(self.chainfile):
        p, np, nc = nest2pos(self.chainfile, self.nlive, seed, self.recparams['logL'])
        print 'Did not find posterior file ' + self.postfile
        print 'Sampled ' + str(nc) + '-->' + str(np) +  ' posteriors from chain file ' + self.chainfile + ' using ' + str(self.nlive) + ' live points.'
      else:
        print "No posterior or chain file provided!"
        exit(-1)
      self.posteriors = list(p.transpose())
      self.npost = np
      return p, np

    def readChain(self):
      '''Read chain points from file'''
      if self.recparams == {}:
        self.getParams()
      if os.path.isfile(self.chainfile):
        self.chain = list(genfromtxt(self.chainfile, skip_header=1).transpose())
        self.nchain = len(self.chain[0])
      else:
        print "No chain file provided!"
        exit(-1)

    def getStats(self):
      '''Gets statistics of the recovered parameters'''
      for idx in arange(len(self.recparams)):
        data = self.chain[idx]
        stats.append([min(data), max(data), mean(data), stdev(data), getConfidenceIntervals(hist(data, 1000)[0], linspace(min(data),max(data),1000), [68.0,95.0,99.0])])

    def getOutName(self): #FIXME: get rid of this
      name = 'lalinference-' + str(self.eventno+1) + '-' + ''.join(self.IFOs) + '-' + str(self.gpstime) + '-' + str(self.eventno) + '.dat'
      return name
      
class Xfm:
  def __init__(self, names, varnames, xfmFunc):
    self.names = names
    self.varnames = varnames
    self.xfmFunc = xfmFunc

class Label:
  def __init__(self, name):
    self.name = name
    self.color = typecolors[self.name]
    self.marker = typemarkers[self.name]
    self.injvalue = None
    if name in injValues.keys():
      self.injvalue = injValues[self.name]
    self.batches = []
    self.globalpars = {}
    self.globalpars2d = {}
          
  def addGlobalParam(self, gp):
    self.globalpars[gp.name] = GlobalParam(gp.name, gp.min, gp.max, gp.prior_shape, latex_name=gp.latex_name, latex_unit=gp.latex_unit, color=self.color, injected=injValues[l][gp.name], xfm=gp.xfm)

  def addGlobalParam2D(self, gp):
    print "Adding " + str(gp.name)
    self.globalpars2d[gp.name[0]+'-'+gp.name[1]] = GlobalParam2D(gp.parameters, color=self.color, injectedx=injValues[l][gp.name[0]], injectedy=injValues[l][gp.name[1]])

  def addSummaryPlot(self, sax, parname):
    if parname in self.globalpars.keys():
      self.globalpars[parname].replotSummary(sax)
      if self.injvalue is not None:
        sax.axhline(self.injvalue[parname], 0, len(self.globalpars[parname].sourcelist)+1, color=self.color, linestyle='--', label='Injected value ' + self.name)



###########################################
#
#  Function definitions
#
###########################################

def mirrorEdges(data, left, right, margin=0.2):
  '''Mirror data left and right of the edges to reduce kde edge effects'''
  lmirr = left - (data - left)
  rmirr = right - (data - right)
  mwidth = (right - left)*margin
  ledge = left - mwidth
  redge = right + mwidth
  ldata = lmirr[where(lmirr > ledge)[0]]
  rdata = rmirr[where(rmirr < redge)[0]]
  extdata = hstack((ldata,data,rdata))
  return extdata

def mirrorEdges2D(data, left, right, high, low, hmargin=0.2, vmargin=0.2): #FIXME
  '''Mirror data around the 2D box edges to reduce kde edge effects'''
  lmirr = left - (data[:,0] - left) #FIXME
  rmirr = right - (data[:,0] - right) #FIXME
  lomirr = low - (data[:,1] - low) #FIXME
  himirr = high - (data[:,1] - high) #FIXME
  mwidth = (right - left)*hmargin
  mheight = (high - low)*vmargin
  ledge = left - mwidth
  redge = right + mwidth
  loedge = low - mheight
  hiedge = high + mheight
  lpatch = lmirr[where(lmirr[0] > ledge)[0]]
  rpatch = rmirr[where(rmirr[0] < redge)[0]]
  lopatch = lomirr[where(lomirr[:,1] > loedge)[0]]
  hipatch = himirr[where(himirr[:,1] < hiedge)[0]]
  # hirpatch = himirr[where(hirmirr[:,1] < hiedge and hirmirr[:,0] < redge)[0]]
  extdata = vstack((data, lpatch, rpatch, hipatch, lopatch))
  return extdata


def getConfidenceIntervals(post, x, pclist):
  '''Returns a list of value-interval pairs for a given list of confidence values'''
  res = []
  cumpost = integrate.cumtrapz(post, x)
  for pc in pclist:
    left = (1.0 - pc*0.01)*0.5
    right = 1.0 - left
    lowidx = where(cumpost/cumpost[-1] < left)[0] 
    if len(lowidx) < 2:
      lowidx = 0
    else:
      lowidx = lowidx.max()
    highidx = where(cumpost/cumpost[-1] > right)[0]
    if len(highidx) < 2:
      highidx = len(cumpost) - 1
    else:
      highidx = highidx.min()
    interval = array([x[lowidx], x[highidx]])
    res.append((pc, interval))
  return res


def nest2pos(datafile, Nlive, sd, loglidx):
    """
    READ INSPNEST DATAFILES TO OUTPUT POSTERIORS
    """
    seed(sd)
    samps = genfromtxt(datafile)[1:,:]
    #samps = genfromtxt(datafile, skip_header=1)
    length=len(samps)
    weight = -linspace(1,length/Nlive,length)
    weight = weight + samps[:, loglidx]
    maxwt = max(weight)
    randoms = rand(length)
    pos = zeros(size(samps,1))
    posidx = find(weight>maxwt+log(randoms))
    pos=samps[posidx,:]
    
    Nchains = len(samps)
    Nsamples = len(posidx) 
    
    return pos, Nsamples, Nchains

def ensure_dir(f):
    """
    CREATE FOLDER IF IT DOES NOT EXIST
    """
    if not os.path.exists(f):
        os.makedirs(f)


''' FILE RELATED FUNCTIONS '''  

def decChainFile(fname):
  '''Decomposes a chainfile name into its components'''
  parts = fname.rpartition('.')
  base = parts[0]
  #ext = parts[2]
  (pref, np1, ifos, gps, n) = base.split('-')
  gps = gps.partition('.')[0]
  return (int(n), int(gps), ifos)
  

def matchFiles(location, substrings):
    '''Returns list of filenames containing a list of substrings'''
    rawList=os.listdir(location)
    matchList=[]
    for f in rawList: 
        match = True
        for s in list(substrings):
            match = match and (s in f)
        if match:
            matchList.append(f)
    print "Found " + str(len(matchList)) + " events in " + location+" matching " + str(substrings)

    return matchList

def readBatches(fname):
    '''Reads list of batch info from file'''
    f = open(fname, 'r')
    l = []
    for line in f:
        l.append(line.strip())
    f.close()
    batches = []
    labels = []
    for line in l:
        (cluster, location, seed, nlive, injtype, xmlfile) = line.split()
        batches.append(Batch(cluster, location, int(seed), int(nlive), injtype, xmlfile))
        labels.append(injtype)
    return batches, list(set(labels))
    

def loadfrompickle(filename):
	"""
	LOAD FROM PICKLE
	"""
	print '... Loading file from pickle',filename
	fp = open(filename,'rb')
	return load(fp)

def savetopickle(self,dest):
	"""
	Pull data from remote or local locations
	"""

	# SAVING TO PICKLE
	f = open(dest,'wb')
	dump(self,f,2)

def loadconfigfile(configfile):
	"""
	Load configuration file
	"""

	# CONFIGURATION FILE: CHECKING FILE
	if access(configfile, R_OK): 
		config = SafeConfigParser()
		config.read(configfile)
	else:
		exit("Configuration file: checking file - file does not exist. Abort script\n")

	return config

''' TIGER RELATED FUNCTIONS''' 



if __name__ == "__main__":

###########################################
#
#  Parse arguments (argparser)
#
###########################################


  # ARGPARSER IS ONLY AVAILABLE IN v2.7 OR HIGHER WHICH CANNOT BE FOUND ON ALL CLUSTERS
  parser = op.OptionParser()
  parser.add_option("--individual", action="store_true", dest="individual", help="output individual posteriors",default=False)
  parser.add_option("-a", "--animate", action="store_true", dest="animate", help="compile posterior PDFs into an animated plot", default=False)
  parser.add_option("-b", "--batchfile", type='str', dest="batchfile",  help="File containing the list of batch info: cluster seed location xmlfile, nlive")
  parser.add_option("-B", "--BayesCut", type='float', dest="Bcut", default=0.0, help="Bayes factor cutoff (0.0)")
  parser.add_option("-c", "--config", type='str', dest="configfile", help="File containing configuration options", metavar="FILE",default="")
  parser.add_option("-e", "--eventplot", action="store_true", dest="eventplot", help="plot sample of posterior points (mpl v1.3)", default=False)
  parser.add_option("-k", "--kde", type='choice', dest="kde", help="Choose kde to use", choices=["GAUSS"], metavar="NAME",default=["GAUSS"])
  parser.add_option("-n", "--nMax", type='int', dest="nMax", help="number of sources to consider, 0 to plot all sources in range",default=0)
  parser.add_option("-o", "--output", type='str', dest="outputfolder", help="outputfolder", metavar="FOLDER",default=".")
  parser.add_option("-s", "--seed", type='int', dest="seed", help="seed number for RNG", metavar="INT", default=1234)
  parser.add_option("-S", "--SNRCut", type='float', dest="SNRcut", default=0.0, help="SNR cutoff (0.0)")

  (args, dumm) = parser.parse_args()

  '''
  parser = argparse.ArgumentParser(description="Post-process output files of tidal runs.")
  parser.add_argument("-b", "--batchfile", type=str, dest="batchfile", required=True, help="File containing the list of batch info: cluster location seed nLive label xmlfile")
  parser.add_argument("-o", "--output", type=str, dest="outputfolder", help="outputfolder", metavar="FOLDER",default=".")
  parser.add_argument("-s", "--seed", type=int, dest="seed", help="seed number for RNG", metavar="INT", default=1234)
  parser.add_argument("-n", "--nMax", type=int, dest="nMax", help="number of sources to consider, 0 to plot all sources in range",default=0)
  parser.add_argument("--individual", action="store_true", dest="indivdiual", help="output individual posteriors",default=False)
  parser.add_argument("-B", "--BayesCut", type=float, dest="Bcut", default=0.0, help="Bayes factor cutoff (0.0)")
  parser.add_argument("-S", "--SNRCut", type=float, dest="SNRcut", default=0.0, help="SNR cutoff (0.0)")
  parser.add_argument("-e", "--eventplot", action="store_true", dest="eventplot", help="plot sample of posterior points (mpl v1.3)", default=False)
  parser.add_argument("-a", "--animate", action="store_true", dest="animate", help="compile posterior PDFs into an animated plot", default=False)
  
  args = parser.parse_args()
  '''
  
  #infolder = os.path.normpath(args.inputlist) 
  batchfile = args.batchfile
  outputfolder = os.path.normpath(args.outputfolder)
  myseed = args.seed
  Bcut = args.Bcut
  NetSNRcut = args.SNRcut
  nMax = args.nMax
  plotIndividual = args.individual
  eplot = args.eventplot
  animate = args.animate
  configfile = args.configfile
  
  batches, labels = readBatches(batchfile)
  
  '''Create output folders'''
  ensure_dir(outputfolder)
  if plotIndividual:
    for l in labels:
      ensure_dir(os.path.join(outputfolder, 'combined_posteriors', l))



###########################################
#
#  Main body
#
###########################################

  ''' CREATE INJECTION VALUES FROM LABELS '''
  params=[]
  gparams=[]
  params2d=[]
  gparams2d=[]
  injValues=[]
  # params, gparams, params2d, gparams2d, injValues = setupTIGERparams(labels) # Set up parameters for TIGER (dchis) if posteriors are needed

  labeldict = {}

  ''' LOOP OVER LABELS AND CREATE DATA SILCES '''
  for l in labels:
    thisLabel = Label(l)
    for par in gparams:
      thisLabel.addGlobalParam(par)
    for par in gparams2d:
      thisLabel.addGlobalParam2D(par)

    print "Processing ", l
      
    LocationFile = open(os.path.join(outputfolder,'locations_'+l+'.txt'),'w')

    '''CREATE LISTS CONTAINING DATA TO BE PLOTTED'''

    iEventCount=0
    for b in batches:
      if b.label != l:
        continue
      thisLabel.batches.append(b)
      print "Processing batch with seed ", b.seed
      b.populateSources()
      
# HERE IS A GOOD PLACE TO LOOK FOR RESULTS FROM ALL HYPOTHESES AND POPULATE "COMBINED SOURCES"   
# A COMBINED SOURCE CAN BE A CLASS WITH A SOURCE DICTIONARY (WITH HYPOTHESES NAMES AS KEYS AND SOURCE OBJECTS AS VALUES),
# A COMBINED ODDS RATIO, AND POSSIBLY SNR, SEED, GPSTIME ETC
      
      # THIS PART NEEDS TO BE ADJUSTED FOR TIGER
      
      print "Populated " + str(len(b.sources)) + " sources"

      for s in b.sources:
        if iEventCount >= nMax and nMax!=0:
          break

        '''Apply Bayes factor cutoff'''
        if not s.isFound(Bcut):
          continue

        '''Apply SNR cutoff'''
        if not s.SNRpass(NetSNRcut):
          continue
                                
                                
        '''Derive the posterior points'''
        s.setPosteriors()
        
        '''Check for xfm'd parameters'''
        for par in params+gparams:
          #print par.name
          if par.xfm is not None and par.name not in s.recparams.keys():
            #print "Adding transformed parameters to source: " + str(par.xfm.names)
            s.getXfmParams(par.xfm)
            
        '''Add to list of sources to process'''
        for gp in thisLabel.globalpars.values():
          gp.addSource(s)
          

        iEventCount+=1
        LocationFile.write(str(iEventCount)+'\t'+b.basedir+'\t'+s.chainfile+'\t'+s.Bfile+'\t'+s.snrfile+'\t'+ s.postfile +'\n')
		
        '''Output progress'''
        print str(iEventCount)+": "+ os.path.split(s.chainfile)[1] + " event " + str(s.id) +" --> " + str(len(s.posteriors[0]))+" pts, SNR: "+str(s.SNR)+", Bfactor: "+str(s.Bfac)
               
    # CLOSE LOCATIONS FILE
    LocationFile.close()
    
    '''Combine posteriors accross sources for each parameter'''
    for par in gparams:
      thisLabel.globalpars[par.name].combinePosteriors()
    
    labeldict[l] = thisLabel
   
   
   
   
'''   
   # PLOT SUMMARY PLOTS
  for par in gparams:
    #Plot summary plot
    fig_sum = plt.figure()
    ax_sum = fig_sum.add_subplot(111)
    ax_sum.set_title("Combined posterior summary for " + par.latex_name)
    #Plot final combined posteriors and summary
    fig_multi = plt.figure()
    ax_msum = []
    ax_msum.append( fig_multi.add_subplot(2,2,1))
    nax = ceil(sqrt(len(gparams)))
    for j in arange(len(labels)):
      labeldict[labels[j]].addSummaryPlot(ax_sum, par.name)
      ax_msum = fig_multi.add_subplot(2,2,1+j)

    ax_sum.legend()
    fig_sum.savefig(os.path.join(outputfolder, par.name + '_summary.png'))
    
    
#    for l in labels:  
#      fig_comb = plt.figure()
#      ax_comb = fig_sum.add_subplot(2,2,j+1)

  for l in labels:
    for gp in labeldict[l].globalpars.values() + labeldict[l].globalpars2d.values():
      if plotIndividual:
        print "plotting label " + l
        gp.plotCombinedPDFs(os.path.join(outputfolder, 'combined_posteriors'))
      if animate:
        gp.animatePDFs(os.path.join(outputfolder, 'combined_posteriors'))
'''
