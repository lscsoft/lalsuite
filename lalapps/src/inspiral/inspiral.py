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
    
    self.set_stdout_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/trigtotmplt-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('trigtotmplt.sub')


class IncaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_inca job used by the inspiral pipeline. The static options are
  read from the section [inca] in the ini file.  The stdout and stderr from
  the job are directed to the logs directory. The job always runs in the
  scheduler universe. The path to the executable is determined from the ini
  file.
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

    self.set_stdout_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/inca-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('inca.sub')


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

    return basename + '-' + str(self.get_start()) + '-' + \
      str(self.get_end() - self.get_start()) + '.xml'


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

  def make_trigbank(self,chunk,source_ifo,dest_ifo,usertag=None,ifo_tag=None):
    """
    Sets the name of triggered template bank file.
    chunk = the analysis chunk that is being 
    source_ifo = the name of the ifo that the triggers come from
    dest_ifo = the name of the ifo that the templates will be used for
    usertag = usertag to tag the output filename with
    ifo_tag = string to tag source interferometers, overrides the source_ifo
    for naming files
    """
    if chunk.trig_start():
      self.set_start(chunk.trig_start())
    else:
      self.set_start(chunk.start())
    if chunk.trig_end():
      self.set_end(chunk.trig_end())
    else:
      self.set_end(chunk.end())

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


##############################################################################
# some functions to make life easier later
class AnalyzedIFOData:
  """
  Contains the information for the data that needs to be filtered.
  """
  def __init__(self,chunk,node):
    self.__analysis_chunk = chunk
    self.__dag_node = node

  def set_chunk(self,chunk):
    self.__analysis_chunk = chunk

  def get_chunk(self):
    return self.__analysis_chunk

  def set_dag_node(self,node):
    self.__dag_node = node

  def get_dag_node(self):
    return self.__dag_node

def chunk_in_segment(chunk,seg):
  if ( 
    chunk.start() >= seg.start() and chunk.start() <= seg.end() 
    ) or (
    chunk.end() >= seg.start() and chunk.end() <= seg.end()
    ) or (
    seg.start() >= chunk.start() and seg.end() <= chunk.end() ):
    return 1
  else:
    return 0

def chunks_overlap(chunk1,chunk2):
  if ( 
    chunk1.start() >= chunk2.start() and chunk1.start() <= chunk2.end()
    ) or (
    chunk1.end() >= chunk2.start() and chunk1.end() <= chunk2.end()
    ) or (
    chunk1.start() >= chunk2.start() and chunk1.end() <= chunk2.end() ):
    return 1
  else:
    return 0

def analyze_lho_ifo(l1_lho_data,l1_chunks_analyzed,ifo_data,ifo_name,
  insp_job,trig_job,df_job,snr,chisq,pad,prev_df,
  do_datafind,do_trigbank,do_inspiral,dag,usertag=None):
  """
  Analyze the data from a Hanford instrument. Since the way we treat H1 and
  H2 is symmetric, this function is the same for both instruments. Returns the
  last LALdataFind job that was executed and the chunks analyzed.
  l1_lho_data = the overlap of L1 data with the Hanford IFOs
  l1_chunks_analyzed = the L1 master chunks that have been analyzed
  lho_data = the master science segments the Hanford IFO
  ifo_name = the name of the Hanford IFO
  insp_job = the condor job that we should use to analyze Hanford data
  trig_job = the trigtomplt job used for the LHO inspiral search
  df_job = the condor job to find the data
  snr = the signal-to-noise threshold for this IFO
  chisq = the chi squared threshold for this IFO
  pad = data start/end padding
  prev_df = the previous LALdataFind job that was executed
  do_datafind = should we do the datafind?
  do_trigbank = should we do the triggered bank?
  do_inspiral = should we do the inspiral?
  dag = the DAG to attach the nodes to
  """
  chunks_analyzed = []
  # loop over the master science segments
  for seg in ifo_data:
    # make sure that we do not do a data find more than once per science segment
    done_df_for_seg = None
    # loop over the master analysis chunks in the science segment
    for chunk in seg:
      done_this_chunk = 0
      # now loop over all the L1 data that we need to filter
      for seg_to_do in l1_lho_data:
        # if the current chunk is in one of the segments we need to filter
        if chunk_in_segment(chunk,seg_to_do) and not done_this_chunk:
          # make sure we only filter the master chunk once
          done_this_chunk = 1
          # make sure we have done one and only one datafind for the segment
          if not done_df_for_seg:
            df = DataFindNode(df_job)
            df.set_ifo(ifo_name)
            df.set_start(seg.start() - pad)
            df.set_end(seg.end() + pad)
            if prev_df: df.add_parent(prev_df)
            if do_datafind: dag.add_node(df)
            prev_df = df
            done_df_for_seg = 1
          # make a trigbank for the master LHO chunk: to do this, we need
          # the master L1 chunks that overlap with this LHO chunk as input
          trigbank = TrigToTmpltNode(trig_job)
          trigbank.make_trigbank(chunk,'L1',ifo_name,usertag)
          for l1_done in l1_chunks_analyzed:
            if chunks_overlap(chunk,l1_done.get_chunk()):
              trigbank.add_var_arg(l1_done.get_dag_node().get_output())
              if do_inspiral: trigbank.add_parent(l1_done.get_dag_node())
          if do_trigbank: dag.add_node(trigbank)
          # analyze the LHO master chunk with this triggered template bank
          insp = InspiralNode(insp_job)
          insp.set_start(chunk.start())
          insp.set_end(chunk.end())
          insp.set_ifo(ifo_name)
          insp.set_ifo_tag('L1')
          insp.add_var_opt('snr-threshold',snr)
          insp.add_var_opt('chisq-threshold',chisq)
          insp.add_var_opt('trig-start-time',chunk.trig_start())
          insp.add_var_opt('trig-end-time',chunk.trig_end())
          insp.set_cache(df.get_output())
          insp.set_bank(trigbank.get_output())
          if do_datafind: insp.add_parent(df)
          if do_trigbank: insp.add_parent(trigbank)
          if do_inspiral: dag.add_node(insp)
          # store this chunk in the list of filtered L1 data
          chunks_analyzed.append(AnalyzedIFOData(chunk,insp))
  
  return tuple([prev_df,chunks_analyzed])
