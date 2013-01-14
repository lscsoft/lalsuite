"""
  Classes needed for the time domain pulsar parameter estimation pipeline.

  This has been greatly inspired (well blatently hacked together from) the stochastic pipeline by
  Adam Mercer and the frequency domain binary pulsar search by Chris Messenger.
  
  Matt Pitkin 07/03/08
  
"""

__author__ = 'Matt Pitkin <matthew@astro.gla.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string
import exceptions
import os
from glue import pipeline

# And begin ...

class parameterJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_heterodyne_pulsar job to coarse heterodyne the data.
  """
  def __init__(self,cp):
    # cp = ConfigParser object from which options are read.
    self.__executable = cp.get('exec', 'parameter')#gets the executable from.ini
    self.__universe = cp.get('condor', 'universe')
    pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
    pipeline.AnalysisJob.__init__(self, cp)
    
    # set log files for job
    self.set_stdout_file('logs/pulsar_parameter_estimation-$(cluster).out')
    self.set_stderr_file('logs/pulsar_parameter_estimation-$(cluster).err')
    self.set_sub_file('pulsar_parameter_estimation.sub')

class parameterNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A heterodyneNode runs an instance of lalapps_heterodyne_pulsar in coarse heterodyne mode in
  a condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_CalculateSensitivity.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    
    # initilise job variables
    self.__detectors = None
    self.__param_file = None
    self.__output_dir = None
    self.__pulsar = None
    self.__input_dir = None
    self.__maxh0 = None
    self.__minh0 = None
    self.__h0steps = None
    self.__maxphi0 = None
    self.__minphi0 = None
    self.__phi0steps = None
    self.__maxci = None
    self.__minci = None
    self.__cisteps = None
    self.__maxpsi = None
    self.__minpsi = None
    self.__psisteps = None
    self.__iterations = None
    self.__burnin = None
    self.__temperature = None
    self.__outputrate = None
    self.__earthephem = None
    self.__sunephem = None
    self.__timeephem = None
    self.__covfile = None
    self.__ul = None
    self.__h0prior = None
    self.__phiprior = None
    self.__psiprior = None
    self.__iotaprior = None
    self.__h0mean = None
    self.__phimean = None
    self.__psimean = None
    self.__iotamean = None
    self.__h0sig = None
    self.__phisig = None
    self.__psisig = None
    self.__iotasig = None
    self.__h0width = None
    self.__phiwidth = None
    self.__psiwidth = None
    self.__ciwidth = None
    self.__priorfile = None

  def set_output_dir(self,output_dir):
    # set output directory
    self.add_var_opt('output-dir',output_dir)
    self.__output_dir = output_dir
  
  def set_priorfile(self,priorfile):
    # set a file containing a h0 vs cos(iota) distribution to be used as a prior
    self.add_var_opt('priorfile', priorfile)
    self.__priorfile = priorfile
  
  def set_detectors(self, detectors):
    # set detectors
    self.add_var_opt('detectors', detectors)
    self.__detectors = detectors
    
  def set_param_file(self, param_file):
    # set pulsar parameter file
    self.add_var_opt('par-file', param_file)
    self.__param_file = param_file
  
  def set_input_dir(self, input_dir):
    # set the input directory
    self.add_var_opt('input-dir', input_dir)
    self.__input_dir = input_dir

  def set_pulsar(self,pulsar):
    # set pulsar name
    self.add_var_opt('pulsar',pulsar)
    self.__pulsar = pulsar
    
  def set_ul(self,ul):
    # set to output the upper limit
    self.add_var_opt('dob-ul', ul) # no variable required
    self.__ul = ul

  def set_mcmc(self):
    # set to perform posterior calculation via an MCMC
    self.add_var_opt('mcmc', '')

  def set_h0width(self,h0width):
    # set width of h0 prior
    self.add_var_opt('h0-width',h0width)
    self.__h0width = h0width

  def set_phiwidth(self,phiwidth):
    # set width of phi0 prior
    self.add_var_opt('phi0-width',phiwidth)
    self.__phiwidth = phiwidth

  def set_psiwidth(self,psiwidth):
    # set width of psi prior
    self.add_var_opt('psi-width',psiwidth)
    self.__psiwidth = psiwidth

  def set_ciwidth(self,ciwidth):
    # set width of ci prior
    self.add_var_opt('ci-width',ciwidth)
    self.__ciwidth = ciwidth

  def set_usepriors(self):
    # set to use priors
    self.add_var_opt('use-priors', '')

  def set_h0prior(self,h0prior):
    # set the h0 prior
    self.add_var_opt('h0prior',h0prior)
    self.__h0prior = h0prior

  def set_phiprior(self,phiprior):
    # set the phi0 prior
    self.add_var_opt('phi0prior',phiprior)
    self.__phiprior = phiprior

  def set_psiprior(self,psiprior):
    # set the psi prior
    self.add_var_opt('psiprior',psiprior)
    self.__psiprior = psiprior

  def set_iotaprior(self,iotaprior):
    # set the iota prior
    self.add_var_opt('iotaprior',iotaprior)
    self.__iotaprior = iotaprior

  def set_h0mean(self,h0mean):
    # set mean of h0 prior
    self.add_var_opt('h0mean',h0mean)
    self.__h0mean = h0mean

  def set_phimean(self,phimean):
    # set mean of phi0 prior
    self.add_var_opt('phi0mean',phimean)
    self.__phimean = phimean

  def set_psimean(self,psimean):
    # set mean of psi prior
    self.add_var_opt('psimean',psimean)
    self.__psimean = psimean

  def set_iotamean(self,iotamean):
    # set mean of iota prior
    self.add_var_opt('iotamean',iotamean)
    self.__iotamean = iotamean

  def set_h0sig(self,h0sig):
    # set sigma of h0 prior
    self.add_var_opt('h0sig',h0sig)
    self.__h0sig = h0sig

  def set_phisig(self,phisig):
    # set sigma of phi0 prior
    self.add_var_opt('phi0sig',phisig)
    self.__phisig = phisig

  def set_psisig(self,psisig):
    # set sigma of psi prior
    self.add_var_opt('psisig',psisig)
    self.__psisig = psisig

  def set_iotasig(self,iotasig):
    # set sigma of iota prior
    self.add_var_opt('iotasig',iotasig)
    self.__iotasig = iotasig

  def set_covfile(self,covfile):
    # set the covariance matrix file
    self.add_var_opt('covariance',covfile)
    self.__covfile = covfile

  def set_usecov(self):
    # set to use a covariance matrix for prior/proposal
    self.add_var_opt('use-cov', '')
    
  def set_verbose(self):
    # set verbose flag
    self.add_var_opt('verbose', '') # no variable required
   
  def set_onlyjoint(self):
    # only perform the joint posterior calculation for MCMC
    self.add_var_opt('only-joint', '')

  def set_maxh0(self, maxh0):
    # set maximum h0 value in grid
    self.add_var_opt('maxh0',maxh0)
    self.__maxh0 = maxh0
    
  def set_minh0(self, minh0):
    # set minimum h0 value in grid
    self.add_var_opt('minh0',minh0)
    self.__minh0 = minh0
    
  def set_h0steps(self, h0steps):
    # set number of grid points in h0 grid
    self.add_var_opt('h0steps',h0steps)
    self.__h0steps = h0steps
  
  def set_maxphi0(self, maxphi0):
    # set maximum phi0 value in grid
    self.add_var_opt('maxphi0',maxphi0)
    self.__maxphi0 = maxphi0
    
  def set_minphi0(self, minphi0):
    # set minimum phi0 value in grid
    self.add_var_opt('minphi0',minphi0)
    self.__minphi0 = minphi0
  
  def set_phi0steps(self, phi0steps):
    # set number of grid points in phi0 grid
    self.add_var_opt('phi0steps', phi0steps)
    self.__phi0steps = phi0steps
  
  def set_maxci(self, maxci):
    # set maximum cos(iota) value in grid
    self.add_var_opt('maxci',maxci)
    self.__maxci = maxci
    
  def set_minci(self, minci):
    # set minimum cos(iota) value in grid
    self.add_var_opt('minci',minci)
    self.__minci = minci
  
  def set_cisteps(self, cisteps):
    # set number of grid points in cos(iota) grid
    self.add_var_opt('cisteps', cisteps)
    self.__cisteps = cisteps
  
  def set_maxpsi(self, maxpsi):
    # set maximum psi value in grid
    self.add_var_opt('maxpsi',maxpsi)
    self.__maxspi = maxspi
    
  def set_minpsi(self, minpsi):
    # set minimum psi value in grid
    self.add_var_opt('minpsi',minpsi)
    self.__minpsi = minpsi
  
  def set_psisteps(self, psisteps):
    # set number of grid points in psi grid
    self.add_var_opt('psisteps', psisteps)
    self.__psisteps = psisteps

  def set_iterations(self, iterations):
    # set the number of MCMC iterations
    self.add_var_opt('iterations', iterations)
    self.__iterations = iterations

  def set_burnin(self, burnin):
    # set the number of burn in iterations
    self.add_var_opt('burn-in', burnin)
    self.__burnin = burnin

  def set_temperature(self, temperature):
    # set the burn in temperature
    self.add_var_opt('temperature', temperature)
    self.__temperature = temperature

  def set_outputrate(self, outputrate):
    # set the MCMC output rate
    self.add_var_opt('output-rate', outputrate)
    self.__outputrate = outputrate

  def set_earth(self, earthephem):
    # set the earth ephemeris file
    self.add_var_opt('earth-ephem', earthephem)
    self.__earthephem = earthephem

  def set_sun(self, sunephem):
    # set the sun ephemeris file
    self.add_var_opt('sun-ephem', sunephem)
    self.__sunephem = sunephem
    
  def set_ephem_time(self, timeephem):
    # set the time correction ephemeris file
    self.add_var_opt('time-ephem',timeephem)
    self.__timeephem = timeephem
