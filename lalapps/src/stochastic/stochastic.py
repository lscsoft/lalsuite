"""
$Id$

Classes needed for the stochastic analysis pipeline.

This script produced the necessary condor submit and dag files to run
the standalone stochastic code on LIGO data.
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
    Returns the IFO code of the first interferometer.
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
    Returns the IFO code of the second interferometer.
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

    if ((ifo == 'H2') and (start >= 729273613) and (start <= 734367613)):
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

    if ((ifo == 'H2') and (start >= 729273613) and (start <= 734367613)):
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
    out = self.get_ifo_one() + self.get_ifo_two() + '-'
    out = out + "stochastic" + '-' + str(self.get_start()) + '-' + str(self.get_stop())
    return out + '.xml'

# vim: et
