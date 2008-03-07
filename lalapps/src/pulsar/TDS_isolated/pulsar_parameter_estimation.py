"""
  Classes needed for the time domain pulsar parameter estimation pipeline.

  This has been greatly inspired (well blatently hacked together from) the stochastic pipeline by
  Adam Mercer and the frequency domain binary pulsar search by Chris Messenger.
  
  Matt Pitkin 07/03/08
  
  $Id$
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
   
    self.add_condor_cmd('WantRemoteIO','false')
    self.add_condor_cmd('local_files','*')
    
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
     
  def set_output_dir(self,output_dir):
    # set output directory
    self.add_var_opt('output-dir',output_dir)
    self.__output_dir = output_dir
    
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
    
  def set_ul(self):
    # set calibration flag
    self.add_var_opt('dob-ul', '') # no variable required
    
  def set_verbose(self):
    # set verbose flag
    self.add_var_opt('verbose', '') # no variable required
   
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