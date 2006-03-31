"""
  Classes needed for the time domain pulsar heterodyne pipeline.

  This has been greatly inspired (well blatently hacked together from) the stochastic pipeline by
  Adam Mercer and the frequency domain binary pulsar search by Chris Messenger.
  
  Matt Pitkin 30/03/06
  
  $Id$
"""

__author__ = 'Matt Pitkin <matthew@astro.gla.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string
import exceptions
import os
from glue import pipeline

class binarypulsarError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args
    
# And begin ...

class heterodyneJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_heterodyne_pulsar job to coarse heterodyne the data.
  """
  def __init__(self,cp):
    # cp = ConfigParser object from which options are read.
    self.__executable = cp.get('exec', 'heterodyne') # gets the executable from .ini
    self.__universe = cp.get('condor', 'universe')
    pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
    pipeline.AnalysisJob.__init__(self, cp)
    
    self.add_condor_cmd('USE_NFS','True')
    self.add_condor_cmd('WantRemoteIO','false')
    self.add_condor_cmd('local_files','*')
    
    # set log files for job
    self.set_stdout_file('logs/heterodyne_pulsar-$(cluster).out')
    self.set_stderr_file('logs/heterodyne_pulsar-$(cluster).err')
    self.set_sub_file('heterodyne_pulsar.sub')

class heterodyneNode(pipeline.CondorDAGJob, pipeline.AnalysisJob):
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
    self.__ifo = None
    self.__param_file = None
    self.__filter_knee = None
    self.__sample_rate = None
    self.__resample_rate = None
    self.__data_file = None
    self.__channel = None
    self.__seg_list = None
    self.__data_file = None
    self.__output_dir = None
    self.__het_flag = None
    self.__pulsar = None
      
  def set_data_file(self,data_file):
    # set file containing data to be heterodyne (either list of frames or coarse het output)
    self.add_var_opt('data-file',data_file)
    self.__data_file = data_file
    
  def set_output_dir(self,output_dir)
    # set output directory
    self.add_var_opt('output-dir',output_dir)
    self.__output_dir = output_dir
    
  def set_ifo(self, ifo):
    # set detector
    self.add_var_opt('ifo', ifo)
    self.__ifo = ifo
    
  def set_param_file(self, param_file):
    # set pulsar parameter file
    self.add_var_opt('param-file', param_file)
    self.__param_file = param_file
  
  def set_param_file_update(self, param_file_update):
    # set file containing updated pulsar parameters
    self.add_var_opt('param-file-update',param_file_update)
    
  def set_ephem_earth_file(self, ephem_earth_file):
    # set the file containing the earth's ephemeris
    self.add_var_opt('ephem-earth-file', ephem_earth_file)
    
  def set_ephem_sun_file(self, ephem_sun_file):
    # set the file containing the sun's ephemeris
    self.add_var_opt('ephem-sun-file', ephem_sun_file)

  def set_pulsar(self,pulsar):
    # set pulsar name
    self.add_var_opt('pulsar',pulsar)
    self.__pulsar = pulsar
    
  def set_het_flag(self,het_flag):
    # set heterodyne flag
    self.add_var_opt('heterodyne-flag',het_flag)
    self.__het_flag = het_flag
    
  def set_filter_knee(self,filter_knee):
    # set filter knee frequency
    self.add_var_opt('filter-knee',filter_knee)
    self.__filter_knee = filter_knee
  
  def set_channel(self,channel):
    # set channel containing data from frames
    self.add_var_opt('channel',channel)
    self.__channel = channel
    
  def set_sample_rate(self,sample_rate):
    # set sample rate of input data
    self.add_var_opt('sample-rate',sample_rate)
    self.__sample_rate = sample_rate
    
  def set_resample_rate(self,resample_rate):
    # set resample rate for output data
    self.add_var_opt('resample_rate',resample_rate)
    self.__resample_rate = resample_rate
    
  def set_stddev_thresh(self,stddev_thresh):
    # set standard deviation threshold at which to remove outliers
    self.add_var_opt('stddev-thresh',stddev_thresh)
    
  def set_calibrate(self):
    # set calibration flag
    self.add_var_opt('calibrate')
    
  def set_response_function(self,response_function):
    # set reponse function file
    self.add_var_opt('response-file',response_function)
    
  def set_coefficient_file(self,coefficient_file):
    # set the file containing the calibration coefficients (e.g alpha and gammas)
    self.add_var_opt('coefficient-file',coefficient_file)
    
  def set_sensing_function(self,sensing_function):
    # set file containin the sensing function for calibration
    self.add_var_opt('sensing-function',sensing_function)
    
  def set_open_loop_gain(self,open_loop_gain):
    # set file containing the open loop gain for calibration
    self.add_var_opt('open-loop-gain',open_loop_gain)
   
