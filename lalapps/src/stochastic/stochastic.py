"""
Classes needed for the stochastic analysis pipeline.
This script produces the necessary condor submit and dag files to run
the standalone stochastic code on LIGO/GEO data.
"""


__author__ = 'Adam Mercer <ram@star.sr.bham.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
import os
from glue import pipeline


def version():
  """
  Return version
  """
  return __version__


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
    @param cp: ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','stochastic')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['stochastic']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/stochastic-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/stochastic-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('stochastic.sub')


class StochasticNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  An StochaticNode runs an instance of the stochastic code in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of lalapps_stochastic.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo_one = None
    self.__ifo_two = None
    self.__f_min = None
    self.__f_max = None
    self.__usertag = None
    self.__output_dir = None

  def set_ifo_one(self, ifo):
    """
    Set the interferometer code to use as IFO One.
    @param ifo: IFO code (e.g. L1, H1, H2 or G1).
    """
    self.add_var_opt('ifo-one', ifo)
    self.__ifo_one = ifo

  def set_ifo_two(self, ifo):
    """
    Set the interferometer code to use as IFO Two.
    @param ifo: IFO code (e.g. L1, H1, H2 or G1).
    """
    self.add_var_opt('ifo-two', ifo)
    self.__ifo_two = ifo

  def set_cache_one(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache-one argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('frame-cache-one', file)

  def set_cache_two(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache-two argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('frame-cache-two', file)

  def set_calibration_one(self,ifo,start):
    """
    Set the path to the calibration cache file for the given IFO. During
    S2, the Hanford 2km IFO had two calibration epochs, so if the start
    time is during S2, we use the correct cache file.
    @param ifo: IFO code (e.g. L1, H1, H2 or G1).
    @param start: GPS time for calibration.
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
    @param ifo: IFO code (e.g. L1, H1, H2 or G1).
    @param start: GPS time for calibration.
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

  def set_f_min(self,f_min):
    """
    Set the minimum frequency
    @param f_min: Minimum frequency.
    """
    self.add_var_opt('f-min',f_min)
    self.__f_min = f_min

  def set_f_max(self,f_max):
    """
    Set the maximum frequency
    @param f_max: Maximum frequency.
    """
    self.add_var_opt('f-max',f_max)
    self.__f_max = f_max

  def set_f_ref(self,f_ref):
    """
    Set the reference frequency
    @param f_ref: Reference frequency.
    """
    self.add_var_opt('f-ref',f_ref)

  def set_user_tag(self,usertag):
    """
    Set the user tag
    @param usertag: User tag to append to job.
    """
    self.add_var_opt('user-tag',usertag)
    self.__usertag = usertag

  def set_output_dir(self,dir):
    """
    Set the output directory
    @param dir: Directory for output files.
    """
    self.add_var_opt('output-dir',dir)
    self.__output_dir = dir

  def get_output(self):
    """
    Returns the file name of output from the stochastic code. This must be
    kept synchronized with the name of the output file in stochastic.c.
    """
    if not self.get_start() or not self.get_end() or not \
      self.__ifo_one or not self.__ifo_two:
        raise StochasticError, "Start time, end time, ifo one or ifo " \
          "two has not been set"

    basename = self.__ifo_one + self.__ifo_two + '-' + 'STOCHASTIC'

    if self.__usertag:
      basename += '_' + self.__usertag

    if self.__output_dir:
      filename = self.__output_dir + '/' + basename + '-' + \
                 str(self.get_start()) + '-' + \
                 str(self.get_end() - self.get_start()) + '.xml'
    else:
      filename = basename + '-' + str(self.get_start()) + '-' + \
                 str(self.get_end() - self.get_start()) + '.xml'

    return filename


class StoppJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_stopp job used by the stochastic pipeline. The static
  options are read from the section [stopp] in the ini file. The stdout
  and stderr from the job are directed to the logs directory. The path
  to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    @param cp: ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','stopp')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['stopp']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/stopp-$(cluster)-$(process).out')
    self.set_stderr_file('logs/stopp-$(cluster)-$(process).err')
    self.set_sub_file('stopp.sub')


class StoppNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An StoppNode runs an instance of the stochastic stopp code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDagNode that can run an instance of lalapps_stopp.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)


class LSCDataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  An LSCdataFind job used to locate data. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.

  Note: This class overrides the LSCDataFindJob class within
  glue.pipeline, it has support for doing runing datafind jobs for
  multiple frame types within the same DAG. This will be eventually be
  merged into the main glue.pipeline.
  """
  def __init__(self,cache_dir,log_dir,config_file):
    """
    @param cache_dir: the directory to write the output lal cache files to.
    @param log_dir: the directory to write the stderr file to.
    @param config_file: ConfigParser object containing the path to the LSCdataFind
    executable in the [condor] section and a [datafind] section from which
    the LSCdataFind options are read.
    """
    self.__executable = config_file.get('condor','datafind')
    self.__universe = 'scheduler'
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,config_file)
    self.__cache_dir = cache_dir

    # we need a lal cache for files on the localhost
    self.add_opt('match','localhost')
    self.add_opt('lal-cache','')
    self.add_opt('url-type','file')

    self.add_condor_cmd('getenv','True')

    self.set_stderr_file(log_dir + '/datafind-$(macroobservatory)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_stdout_file(self.__cache_dir + '/$(macroobservatory)-$(macrogpsstarttime)-$(macrogpsendtime).cache')
    self.set_sub_file('datafind.sub')

  def get_cache_dir(self):
    """
    returns the directroy that the cache files are written to.
    """
    return self.__cache_dir


class LSCDataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A DataFindNode runs an instance of LSCdataFind in a Condor DAG.

  Note: This class overrides the LSCDataFindNode class within
  glue.pipeline, it has support for doing runing datafind jobs for
  multiple frame types within the same DAG. This will be eventually be
  merged into the main glue.pipeline.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of LSCdataFind.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__observatory = None
    self.__type = None
    self.__output = None
    self.__job = job
   
  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__observatory:
      self.__output = self.__job.get_cache_dir() + '/' + self.__observatory + '-'  
      self.__output += str(self.__start) + '-' + str(self.__end) + '.cache'
      self.add_output_file(self.__output)

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    @param time: GPS start time of query.
    """
    self.add_var_opt('gps-start-time', time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    @param time: GPS end time of query.
    """
    self.add_var_opt('gps-end-time', time)
    self.__end = time
    self.__set_output()

  def set_observatory(self,obs):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford 
    interferometers is stored in the same frame file, this takes the first 
    letter of the IFO (e.g. L or H) and passes it to the --observatory option
    of LSCdataFind.
    @param obs: IFO to obtain data for.
    """
    self.add_var_opt('observatory',obs)
    self.__observatory = obs
    self.__set_output()

  def set_type(self,type):
    """
    Set the frame type to retrieve data for
    @param type: Frame type to obtain data for.
    """
    self.add_var_opt('type',type)
    self.__type = type
    self.__set_output()

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    """
    return self.__output

# vim: et
