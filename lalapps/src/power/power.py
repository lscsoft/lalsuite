"""
Classes needed for the excess power analysis pipeline.
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
from glue import pipeline


class PowerError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class MdcDataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A MdcdataFind job used by the power pipeline. The static options are
  read from the section [mdcdatafind] in the ini file. The stdout from
  MdcdataFind contains the paths to the frame files and is directed to a file
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

    for sec in ['mdcdatafind']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('getenv','True')

    self.set_stderr_file('logs/mdcdatafind-GHLTV-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_stdout_file('mdccache/GHLTV-$(macrogpsstarttime)-$(macrogpsendtime).cache')
    self.set_sub_file('mdcdatafind.sub')


class BurstInjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
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
    self.__executable = cp.get('condor','binjection')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['binjection']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('binjection.sub')


class InspInjJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
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
    self.__executable = cp.get('condor','iinjection')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['iinjection']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/iinj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/iinj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('iinjection.sub')
    


class PowerJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_power job used by the power pipeline. The static options
  are read from the sections [data] and [power] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','power')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['data','power']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/power-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('power.sub')

    
class BurcaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_burca job used by the power pipeline. The static options are
  read from the section [burca] in the ini file.  The stdout and stderr from
  the job are directed to the logs directory.  The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','burca')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['burca']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/burca-$(macrostarttime)-$(macrostoptime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/burca-$(macrostarttime)-$(macrostoptime)-$(cluster)-$(process).err')
    self.set_sub_file('burca.sub')

class VigilJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A vigilance job used to produce summary info
  once triggers have been produced
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','vigil')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['vigil']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/vigil-$(cluster)-$(process).out')
    self.set_stderr_file('logs/vigil-$(cluster)-$(process).err')
    self.set_sub_file('vigil.sub')

    
class MdcDataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A DataFindNode runs an instance of datafind in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of LSCdataFind.
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
    if self.__start and self.__end:
      self.__output = 'mdccache/' + 'GHLTV' + '-' + str(self.__start) 
      self.__output = self.__output + '-' + str(self.__end) + '.cache'

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    time = GPS start time of query.
    """
    self.add_var_opt('gps-start-time', time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    time = GPS end time of query.
    """
    self.add_var_opt('gps-end-time', time)
    self.__end = time
    self.__set_output()

  def set_ifo(self,ifo):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford 
    interferometers is stored in the same frame file, this takes the first 
    letter of the IFO (e.g. L or H) and passes it to the --instrument option
    of LSCdataFind.
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
  

class BurstInjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A InjNode runs an instance of binj code
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_power.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('binjection','user-tag')
    
  def get_output(self,seed):
    """
    Returns the file name of output from the power code. This must be kept
    synchronized with the name of the output file in power.c.
    """
    if not self.get_start() or not self.get_end():
      raise InjError, "Start time or end time has not been set"

    basename = 'HL' + '-INJECTIONS'

    if self.__usertag:
      basename += '_' + str(self.__usertag) 

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

class InspInjNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A InjNode runs an instance of binj code
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_power.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('iinjection','user-tag')
    
  def get_output(self,seed):
    """
    Returns the file name of output from the power code. This must be kept
    synchronized with the name of the output file in power.c.
    """
    if not self.get_start() or not self.get_end():
      raise InjError, "Start time or end time has not been set"

    basename = 'HL' + '-INJECTIONS' + '_' + str(seed)

    if self.__usertag:
      basename += '_' + str(self.__usertag) 

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'
  
  
class PowerNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A PowerNode runs an instance of the power code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_power.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def get_output(self):
    """
    Returns the file name of output from the power code. This must be kept
    synchronized with the name of the output file in power.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise PowerError, "Start time, end time or ifo has not been set"

    basename = self.get_ifo()
    basename += '-' + 'POWER'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

  def set_mdccache(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('mdc-cache', file)
    self.add_input_file(file)

  def set_burstinj(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('burstinjection-file', file)
    self.add_input_file(file)

  def set_inspinj(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('inspiralinjection-file', file)
    self.add_input_file(file)

  def set_siminj(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('siminjection-file', file)
    self.add_input_file(file)

class BurcaNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A BurcaNode runs an instance of the burca code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_burca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')
    self.__start = 0
    self.__end = 0
    self.__output = None
    self.__prefix = None
    
  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end :
      self.__output = 'H1'+ '-' + 'BURCA' + '_' + 'H1H2' + '_' + 'P' + '_' + '0' + '_' + 'H1H2' + '-' + str(self.__start) 
      self.__output = self.__output + '-' + str(self.__end - self.__start) + '.xml'
                     
  def set_ifoa(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('ifo-a', file)
    self.add_input_file(file)

  def set_ifob(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    @param file: calibration file to use.
    """
    self.add_var_opt('ifo-b', file)
    self.add_input_file(file)   

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    time = GPS start time of query.
    """
    self.add_var_opt('start-time',time)
    self.__start = time
    self.__set_output()

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    time = GPS end time of query.
    """
    self.add_var_opt('stop-time',time)
    self.__end = time
    self.__set_output()

  def set_prefix(self,prefix):

    self.__prefix = prefix
    
  def get_H1H2output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    """
    return self.__output

class VigilNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A InjNode runs an instance of binj code
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of vigil.sh.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')
    
