"""
  Classes needed for the mismatch pipeline.
  
  Matt Pitkin 11/09/08
  
"""

__author__ = 'Matt Pitkin <matthew@astro.gla.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string
import exceptions
import os
from glue import pipeline

# And begin ...

class mismatchJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_heterodyne_pulsar job to coarse heterodyne the data.
  """
  def __init__(self,cp):
    # cp = ConfigParser object from which options are read.
    self.__executable = cp.get('exec', 'executable')#gets the executable from.ini
    self.__universe = cp.get('condor', 'universe')
    pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
    pipeline.AnalysisJob.__init__(self, cp)
   
    self.add_condor_cmd('WantRemoteIO','false')
    self.add_condor_cmd('local_files','*')
    
    # set log files for job
    self.set_stdout_file('logs/mismatch-$(cluster).out')
    self.set_stderr_file('logs/mismatch-$(cluster).err')
    self.set_sub_file('mismatch.sub')

class mismatchNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
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
    self.__detector = None
    self.__param_file = None
    self.__output_dir = None
    self.__pulsar = None
    self.__iterations = None
    self.__earthephem = None
    self.__sunephem = None
    self.__covfile = None
    self.__deltat = None
    self.__start = None
    self.__timespan = None

  def set_deltat(self,deltat):
    # set time steps
    self.add_var_opt('deltat',deltat)
    self.__deltat = deltat

  def set_start(self,start):
    # set start time
    self.add_var_opt('start',start)
    self.__start = start

  def set_timespan(self,timespan):
    # set time span
    self.add_var_opt('timespan',timespan)
    self.__timespan = timespan

  def set_output_dir(self,output_dir):
    # set output directory
    self.add_var_opt('output-dir',output_dir)
    self.__output_dir = output_dir
    
  def set_detector(self, detector):
    # set detectors
    self.add_var_opt('detector', detector)
    self.__detector = detector
    
  def set_param_file(self, param_file):
    # set pulsar parameter file
    self.add_var_opt('par-file', param_file)
    self.__param_file = param_file
  
  def set_pulsar(self,pulsar):
    # set pulsar name
    self.add_var_opt('pulsar',pulsar)
    self.__pulsar = pulsar
    
  def set_covfile(self,covfile):
    # set the covariance matrix file
    self.add_var_opt('covariance',covfile)
    self.__covfile = covfile

  def set_iterations(self, iterations):
    # set the number of MCMC iterations
    self.add_var_opt('iterations', iterations)
    self.__iterations = iterations

  def set_earth(self, earthephem):
    # set the earth ephemeris file
    self.add_var_opt('earth', earthephem)
    self.__earthephem = earthephem

  def set_sun(self, sunephem):
    # set the sun ephemeris file
    self.add_var_opt('sun', sunephem)
    self.__sunephem = sunephem
