"""
Classes needed for the inspiral analysis pipeline.
This script produced the necessary condor submit and dag files to run
the standalone inspiral code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
import pipeline


class InspiralError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class DataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A LALdataFind job used by the inspiral pipeline. The static options are
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

    self.set_stderr_file('logs/datafind-$(macroinstrument)-$(macrostart)-$(macroend)-$(cluster)-$(process).err')
    self.set_stdout_file('cache/$(macroinstrument)-$(macrostart)-$(macroend).cache')
    self.set_sub_file('datafind.sub')


class TmpltBankJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_tmpltbank job used by the inspiral pipeline. The static options
  are read from the sections [data] and [tmpltbank] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specfied in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','tmpltbank')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['data','tmpltbank']:
      self.add_ini_opts(cp,sec)
  
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/tmpltbank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/tmpltbank-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('tmpltbank.sub')
    

class InspiralJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inspiral job used by the inspiral pipeline. The static options
  are read from the sections [data] and [inspiral] in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specfied in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','inspiral')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['data','inspiral']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/inspiral-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inspiral-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('inspiral.sub')
    

class TrigToTmpltJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_trigtotmplt job used by the inspiral pipeline. The static
  options are read from the section [trigtotmplt] in the ini file.  The
  stdout and stderr from the job are directed to the logs directory. The job
  always runs in the scheduler universe. The path to the executable is
  determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','trigtotmplt')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['trigtotmplt']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    
    self.set_stdout_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('trigtotmplt.sub')


class IncaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inca job used by the inspiral pipeline. The static options are
  read from the section [inca] in the ini file.  The stdout and stderr from
  the job are directed to the logs directory.  The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','inca')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    for sec in ['inca']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('inca.sub')

class SireJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_sire job used by the inspiral pipeline. The stdout and stderr from
  the job are directed to the logs directory. The path to the executable is 
  determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','sire')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/sire-$(macroifo)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/sire-$(macroifo)-$(cluster)-$(process).err')
    self.set_sub_file('sire.sub')

class Tama2LigoLwJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_tama2ligolw job used by the inspiral pipeline. The stdout and 
  stderr from the job are directed to the logs directory. The path to the 
  executable is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','tama2lw')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)
    
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/tama2ligolw-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/tama2ligolw-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('tama2ligolw.sub')

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


class TmpltBankNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A TmpltBankNode runs an instance of the template bank generation job in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_tmpltbank.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def get_output(self):
    """
    Returns the file name of output from the template bank code. This must
    be kept synchronized with the name of the output file in tmpltbank.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"
    if self.__usertag:
      bank = self.get_ifo() + '-TMPLTBANK_' + self.__usertag + '-' 
      bank = bank + str(self.get_start())
    else:
      bank = self.get_ifo() + '-TMPLTBANK-' + str(self.get_start())
    bank = bank + '-' + str(self.get_end() - self.get_start()) + '.xml'

    self.add_output_file(bank)

    return bank


class InspiralNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An InspiralNode runs an instance of the inspiral code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inspiral.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__usertag = job.get_config('pipeline','user-tag')

  def set_bank(self,bank):
    self.add_var_opt('bank-file', bank)
    self.add_input_file(bank)

  def get_output(self):
    """
    Returns the file name of output from the inspiral code. This must be kept
    synchronized with the name of the output file in inspiral.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo():
      raise InspiralError, "Start time, end time or ifo has not been set"

    basename = self.get_ifo() + '-INSPIRAL'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag

    filename = basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

    self.add_output_file(filename)

    return filename


class TrigToTmpltNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A TrigToTmpltNode runs an instance of the triggered bank generator in a
  Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of inca in trigtotmplt mode.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__output = None

  def make_trigbank(self,chunk,max_slide,source_ifo,dest_ifo,
    usertag=None,ifo_tag=None):
    """
    Sets the name of triggered template bank file.
    chunk = the analysis chunk that is being 
    max_slide = the maximum length of a time slide for background estimation
    source_ifo = the name of the ifo that the triggers come from
    dest_ifo = the name of the ifo that the templates will be used for
    usertag = usertag to tag the output filename with
    ifo_tag = string to tag source interferometers, overrides the source_ifo
    for naming files
    """
    if chunk.trig_start():
      self.set_start(chunk.trig_start() - max_slide)
    else:
      self.set_start(chunk.start() - max_slide)
    if chunk.trig_end():
      self.set_end(chunk.trig_end() + max_slide)
    else:
      self.set_end(chunk.end() + max_slide)

    self.add_var_opt('ifo-a',source_ifo)

    outfile = dest_ifo + '-TRIGBANK_'
    if ifo_tag:
      outfile += ifo_tag
    else:
      outfile += source_ifo
    if usertag:
      outfile += '_' + usertag 
    outfile += '-' + str(chunk.start()) + '-' + str(chunk.dur()) + '.xml'
    self.__output = outfile
    self.add_var_opt('triggered-bank',outfile)

  def get_output(self):
    """
    Returns the name of the triggered template bank file.
    """
    return self.__output
    


class IncaNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  An IncaNode runs an instance of the inspiral coincidence code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo_a = None
    self.__ifo_b = None
    self.__usertag = job.get_config('pipeline','user-tag')
    
  def set_ifo_a(self, ifo):
    """
    Set the interferometer code to use as IFO A.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-a', ifo)
    self.__ifo_a = ifo

  def get_ifo_a(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_a

  def set_ifo_b(self, ifo):
    """
    Set the interferometer code to use as IFO B.
    ifo = IFO code (e.g. L1, H1 or H2).
    """
    self.add_var_opt('ifo-b', ifo)
    self.__ifo_b = ifo

  def get_ifo_b(self):
    """
    Returns the IFO code of the primary interferometer.
    """
    return self.__ifo_b

  def get_output_a(self):
    """
    Returns the file name of output from inca for ifo a. This must be kept
    synchronized with the name of the output file in inca.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo_a():
      raise InspiralError, "Start time, end time or ifo a has not been set"

    basename = self.get_ifo_a() + '-INCA'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag 

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'

  def get_output_b(self):
    """
    Returns the file name of output from inca for ifo b. This must be kept
    synchronized with the name of the output file in inca.c.
    """
    if not self.get_start() or not self.get_end() or not self.get_ifo_b():
      raise InspiralError, "Start time, end time or ifo a has not been set"

    basename = self.get_ifo_b() + '-INCA'

    if self.get_ifo_tag():
      basename += '_' + self.get_ifo_tag()
    if self.__usertag:
      basename += '_' + self.__usertag 

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'


class SireNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A SireNode runs an instance of the single inspiral reader code in a Condor
  DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__ifo = None
    self.__usertag = job.get_config('pipeline','user-tag')

  def set_outputs(self,out_name,usertag=None,slide_time=None,tama_output=None):
    """
    Sets the name of the sire output file.
    out_name = name of sire output file
    usertag = usertag to tag the output filename with
    """
    outfile = out_name
    if usertag:
      outfile += '_' + usertag
    if slide_time:
      if slide_time < 0: outfile += '_SLIDEneg' + str(abs(slide_time))
      else: outfile += '_SLIDE' + str(slide_time)
    summ_file = outfile + '.txt' 
    if tama_output:
      tama_out = outfile + '.dat'
      self.add_var_opt('tama-output',tama_out)

    outfile +=  '.xml'
    self.__output = outfile
    self.add_var_opt('output',outfile)
    self.add_var_opt('summary',summ_file)

  def set_inj_outputs(self,out_name,usertag=None,tama_output=None,
			slide_time=None,cluster=None):
    """
    Sets the name of the sire output file.
    out_name = name of sire output file
    usertag = usertag to tag the output filename with
    """
    outfile = out_name
    if usertag:
      outfile += '_' + usertag
    if slide_time:
      if slide_time < 0: outfile += '_SLIDEneg' + str(abs(slide_time))
      else: outfile += '_SLIDE' + slide_time
    if cluster:
      outfile += '_CLUSTER'  
    summ_file = outfile + '_FOUND.txt' 
    missed_file = outfile + '_MISSED.xml'
    if tama_output:
      tama_out = outfile + '_FOUND.dat'
      self.add_var_opt('tama-output',tama_out)
    outfile +=  '_FOUND.xml'
    self.__output = outfile
    self.add_var_opt('output',outfile)
    self.add_var_opt('summary',summ_file)
    self.add_var_opt('missed-injections',missed_file)
    self.add_var_opt


  def get_output(self):
    """
    Returns the name of the sire output.
    """
    return self.__output


class Tama2LigoLwNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  A Tama2LigoLwNode runs an instance of the tama triggers to LIGO Lw XML 
  converter in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_inca.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__input = None
    self.__output = None
    self.__usertag = job.get_config('pipeline','user-tag')


##############################################################################
# some functions to make life easier later

def overlap_test(interval1, interval2, slide_sec=0):
  """
  Test whether the two intervals could possibly overlap with one of them being
  slid by a maximum time of slide_sec.  Perform three tests:
  1)  Does the start of interval 1 lie within interval 2's range (with the 
    start decremented by slide_sec and the end incremented by slide_sec)
  2)  Does the end of interval 1 lie within interval 2's range (with the start 
    decremented by slide_sec and the end incremented by slide_sec)
  3)  Does interval 1 completely cover (the extended) interval 2, 
    ie is interval_1 start before (interval 2 start - slide_sec) AND 
    interval 1 end after (interval 2 end + slide_sec)
  If any of the above conditions are satisfied then return 1, else 0.
  """
  if ( 
    interval1.start() >= interval2.start() - slide_sec
    and interval1.start() <= interval2.end() + slide_sec
    ) or (
    interval1.end() >= interval2.start() - slide_sec 
    and interval1.end() <= interval2.end() + slide_sec
    ) or (
    interval1.start() <= interval2.start() - slide_sec    
    and interval1.end() >= interval2.end() + slide_sec ):
    return 1
  else:
    return 0



