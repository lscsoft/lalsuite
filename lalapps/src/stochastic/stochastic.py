"""
$Id$

Classes needed for the stochastic analysis pipeline.

This script produced the necessary condor submit and dag files to run
the standalone stochastic code on LIGO data.

Based on inspiral.py by Duncan Brown
"""


__author__ = 'Adam Mercer <ram@star.sr.bham.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]


import string
import exceptions
import pipeline
import os


class StochasticError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class DataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A LALdataFind job used by the stochastic pipeline. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LALdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','datafind')
    self.__universe = 'scheduler'
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['datafind']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',
      """LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH);PYTHONPATH=$ENV(PYTHONPATH)""" )

    self.set_stderr_file('logs/datafind-$(macroinstrument)-$(macrostart)-$(macroend).err')
    self.set_stdout_file('cache/$(macroinstrument)-$(macrostart)-$(macroend).cache')
    self.set_sub_file('datafind.sub')


class StochasticJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_stochastic job used by the stochastic pipeline. The static 
	options are read from the section [stochastic] in the ini file.  The 
	stdout and stderr from the job are directed to the logs directory.
	The path to the executable and the universe is determined from the ini
  file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','stochastic')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['stochastic']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/stoch-$(macrogpsstarttime)-$(macrogpsendtime).out')
    self.set_stderr_file('logs/stoch-$(macrogpsstarttime)-$(macrogpsendtime).err')
    self.set_sub_file('stochastic.sub')


class DataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A DataFindNode runs an instance of datafind in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of LALdataFind.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__instrument = None
    self.__output = None
   
  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__instrument:
      self.__output = 'cache/' + self.__instrument + '-' + str(self.__start) 
      self.__output = self.__output + '-' + str(self.__end) + '.cache'

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    time = GPS start time of query.
    """
    self.add_var_opt('start', time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    time = GPS end time of query.
    """
    self.add_var_opt('end', time)
    self.__end = time
    self.__set_output()

  def set_ifo(self,ifo):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford 
    interferometers is stored in the same frame file, this takes the first 
    letter of the IFO (e.g. L or H) and passes it to the --instrument option
    of LALdataFind.
    ifo = IFO to obtain data for.
    """
    self.add_var_opt('instrument',ifo[0])
    self.__instrument = ifo[0]
    self.__set_output()

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    """
    return self.__output


class StochasticNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  An StochaticNode runs an instance of the stochastic code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_stochastic.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo_one = None
    self.__ifo_two = None
    
  def set_ifo_one(self, ifo):
    """
    Set the interferometer code to use as IFO One.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-one', ifo)
    self.__ifo_one = ifo

  def get_ifo_one(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_one

  def set_ifo_two(self, ifo):
    """
    Set the interferometer code to use as IFO Two.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-two', ifo)
    self.__ifo_two = ifo

  def get_ifo_two(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_two

  def set_cache_one(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache-one argument.
    file = calibration file to use.
    """
    self.add_var_opt('frame-cache-one', file)

  def set_cache_two(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache-two argument.
    file = calibration file to use.
    """
    self.add_var_opt('frame-cache-two', file)

  def set_calibration_one(self,ifo,start):
    """
    Set the path to the calibration cache file for the given IFO. During
    S2, the Hanford 2km IFO had two calibration epochs, so if the start
    time is during S2, we use the correct cache file.
    """
    cal_path = self.job().get_config('calibration','path')

    if (ifo == 'H2'):
      if start < int(self.job().get_config('calibration','H2-cal-epoch-boundary')):
        cal_file = self.job().get_config('calibration','H2-1')
      else:
        cal_file = self.job().get_config('calibration','H2-2')
    else:
      cal_file = self.job().get_config('calibration',ifo)

    cal = os.path.join(cal_path,cal_file)
    self.add_var_opt('calibration-cache-one',cal)

  def set_calibration_two(self,ifo,start):
    """
    Set the path to the calibration cache file for the given IFO. During
    S2, the Hanford 2km IFO had two calibration epochs, so if the start
    time is during S2, we use the correct cache file.
    """
    cal_path = self.job().get_config('calibration','path')

    if (ifo == 'H2'):
      if start < int(self.job().get_config('calibration','H2-cal-epoch-boundary')):
        cal_file = self.job().get_config('calibration','H2-1')
      else:
        cal_file = self.job().get_config('calibration','H2-2')
    else:
      cal_file = self.job().get_config('calibration',ifo)

    cal = os.path.join(cal_path,cal_file)
    self.add_var_opt('calibration-cache-two',cal)			

  def get_output(self):
    """
    Returns the file name of output from the stochastic code. This must be
		kept synchronized with the name of the output file in stochastic.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo_one() or not self.get_ifo_two():
      raise StochasticError, "Start time, end time, ifo one or ifo two has not been set"
    out = 'stoch-' + self.get_ifo_one() + self.get_ifo_two() + '-'
    out = out + str(self.get_start()) + '-' + str(self.get_stop())
    return out + '.dat'
