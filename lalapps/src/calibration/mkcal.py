"""
Classes needed for the frame calibration generation programs.
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
import pipeline


class MkCalError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class MkCalFacJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  Make the calibration frames from SenseMonitor coefficients.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','mkcalfac')
    self.__universe = cp.get('condor','universe') 
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    self.set_stdout_file('logs/mkcalfac-$(macrorun)-$(macroifo)-$(macroversion)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/mkcalfac-$(macrorun)-$(macroifo)-$(macroversion)-$(cluster)-$(process).err')
    self.set_sub_file('mkcalfac.sub')


class MkCalFacNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A MkCalFacNode runs an instance of the calibration factor frame job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_tmpltbank.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo = None

  def set_ifo(self,ifo):
    """
    Set the ifo name.
    @param ifo: two letter ifo code (e.g. L1, H1 or H2).
    """
    self.__ifo = ifo
    self.add_var_opt('ifo', ifo)

  def get_ifo(self):
    """
    Returns the two letter IFO code for this node.
    """
    return self.__ifo


