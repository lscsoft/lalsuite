#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Classes needed for the cosmic string analysis pipeline.
"""


__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'[7:-2]
__version__ = '$Revision$'[11:-2]


from glue import pipeline
from glue import segments


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#


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
    self.__usertag = job.get_config('pipeline','user_tag')

  def get_output(self):
    """
    Returns the file name of output from the ring code. This must be kept
    synchronized with the name of the output file in ring.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise ValueError, "Start time, end time or ifo has not been set"

    return 'triggers/%s-STRINGSEARCH-%s-%s.xml' % (self.get_ifo(), str(self.get_start()), str(self.get_end() - self.get_start()))


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#


def clip_segment(seg, pad, short_segment_duration):
  # clip segment to the length required by lalapps_StringSearch.  if
  #
  #   duration = segment length - padding
  #
  # then
  #
  #   duration / short_segment_duration - 0.5
  #
  # must be an odd integer, therefore
  #
  #   2 * duration + short_segment_duration
  #
  # must be divisble by (4 * short_segment_duration)
  duration = float(abs(seg)) - 2 * pad
  extra = (2 * duration + short_segment_duration) % (4 * short_segment_duration)
  extra /= 2

  # clip segment
  seg = segments.segment(seg[0], seg[1] - extra)

  # bounds must be integers
  if int(seg[0]) != seg[0] or int(seg[1]) != seg[1]:
    raise ValueError, "segment %s does not have integer boundaries" % str(seg)
  seg = segments.segment(int(seg[0]), int(seg[1]))

  # done
  return seg


#
# =============================================================================
#
#                                DAG Job Types
#
# =============================================================================
#


stringjob = None


def init_job_types(config_parser, job_types = ("string",)):
  """
  Construct definitions of the submit files.
  """
  global stringjob

  # lalapps_StringSearch
  if "string" in job_types:
    stringjob = StringJob(config_parser)
