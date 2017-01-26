"""
Classes needed for the excess ring analysis pipeline.
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string
from glue import pipeline


class RingError(Exception):
  def __init__(self, args=None):
    self.args = args


class RingJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_ring job used by the ring pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','ring')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['ring']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/ring-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/ring-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('ring.sub')
    

class RingNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_ring.
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
      raise RingError, "Start time, end time or ifo has not been set"

    basename = self.get_ifo() + '-RING'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'
