"""
Classes needed for the cosmic string analysis pipeline.
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string,sys
import exceptions
from glue import pipeline


class StringError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class DataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A LSCdataFind job used by the string pipeline. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
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


class StrainJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
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
    self.__executable = cp.get('condor','strain')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    for sec in ['strain']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/strain-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/strain-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('strain.sub')
    
class NoiseJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_StringSearch job used by the string pipeline. The static options
  are read from the section in the ini file. The
  stdout and stderr from the job are directed to the logs directory. The job
  runs in the universe specified in the ini file. The path to the executable
  is determined from the ini file.
  """
  def __init__(self,cp,dax=False):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','noise')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,dax)

    for sec in ['noisecomp']:
      self.add_ini_opts(cp,sec)

    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/noisecomp-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/noisecomp-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err')
    self.set_sub_file('noisecomp.sub')

class EpochData(object):
  """ 
  Holds calibration data epochs
  """
  def __init__(self,cp,opts):
    self.file = open(cp.get('epochs','epochs_data'),'r')
    self.epoch_data = []
    for line in self.file:
      if line.strip()[0] != '#':
        self.epoch_data.append(line.split())  
    self.epoch_data.sort()

  def epoch_segs(self):
    tmpList = []
    for line in self.epoch_data:
      tmpList.append(tuple(map(int,(line[0:4]))))
    return tmpList
     
class DataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
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



class StrainNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)


class NoiseNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A RingNode runs an instance of the ring code in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_StringSearch.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

  def add_cache(self,file,frame):
    if isinstance( file, str ):
      # the name of a lal cache file created by a datafind node
      self.add_var_opt(frame, file)
      self.add_input_file(file)
    else:
      cacheFileName = 'cache/'+self.get_name()+'-'+frame+'-'+str(self.get_start())+'-'+str(self.get_end())+'.cache'
      cacheFile = open(cacheFileName,'w')
      self.add_var_opt(frame, cacheFileName)
      self.add_input_file(cacheFileName)
      # check we have an LFN list
      from glue import LDRdataFindClient
      if isinstance( file, LDRdataFindClient.lfnlist ):
        #self.add_var_opt('glob-frame-data',' ')
        # only add the LFNs that actually overlap with this job
        for lfn in file:
          a, b, c, d = lfn.split('.')[0].split('-')
          t_start = int(c)
          t_end = int(c) + int(d)
          if ( t_start <= self.get_end() and t_end >= self.get_start() ):
            self.add_input_file(lfn)
            cacheFile.write(a+' '+b+' '+c+' '+d+' '+lfn+'\n')
        # set the frame type based on the LFNs returned by datafind
        #self.add_var_opt('frame-type',b)
      else:
        raise CondorDAGNodeError, "Unknown LFN cache format"

# Convenience functions to cat together noise output files.
def open_noise_cat_file(dir):
  outfilename = dir + '/' + dir + '.cat'
  try:  outfile = open(outfilename,'w')
  except: 
    sys.stderr.write('Could not open '+outfilename+' for writing')
    sys.exit(1)
  return outfile

def cat_noise_jobs(file,node):
  out = node.get_output_files()
  tmpfile = open(out[0],'r')
  tmpstr = ''.join( tmpfile.readlines() ) 
  try:  file.write(tmpstr)
  except: pass
  return

