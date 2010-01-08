"""
Classes needed for the cosmic string analysis pipeline.
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import exceptions
from glue import pipeline


class StringError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args

class StringJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_StringSearch job used by the string pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','string')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['string']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/string-$(macrochannel)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/string-$(macrochannel)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('string.sub')

class StringNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def get_output(self):
    """
    Returns the file name of output from the ring code. This must be kept
    synchronized with the name of the output file in ring.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise StringError, "Start time, end time or ifo has not been set"

    basename = 'triggers/'+self.get_ifo() + '-STRINGSEARCH'

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'
    
class InjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_binj job used by the power pipeline. The static options
  are read from the sections [data] and [power] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','injection')
    self.__universe = 'local'
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['injection']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/inj-$(macrochannel)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inj-$(macrochannel)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('injection.sub')

class InjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A InjNode runs an instance of binj code
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_power.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    
  def get_output(self):
    """
    Returns the file name of output from the power code. This must be kept
    synchronized with the name of the output file in power.c.
    """
    if not self.get_start() or not self.get_end():
      raise InjError, "Start time or end time has not been set"

    basename = 'HL-INJECTIONS-'+ str(self.get_start()) + \
               '-' +str(self.get_end() - self.get_start()) + '.xml'

    return basename
