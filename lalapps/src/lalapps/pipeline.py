"""
This modules contains objects that make it simple for the user to 
create python scripts that build Condor DAGs to run code on the LSC
Data Grid.
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import os
import string, re
import exceptions
import time
import random
import md5


def s2play(t):
  """
  Return 1 if t is in the S2 playground, 0 otherwise
  t = GPS time to test if playground
  """
  if ((t - 729273613) % 6370) < 600:
    return 1
  else:
    return 0


class CondorError(exceptions.Exception):
  """Error thrown by Condor Jobs"""
  def __init__(self, args=None):
    self.args = args
class CondorJobError(CondorError):
  pass
class CondorSubmitError(CondorError):
  pass
class CondorDAGError(CondorError):
  pass
class CondorDAGNodeError(CondorError):
  pass
class SegmentError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args


class CondorJob:
  """
  Generic condor job class. Provides methods to set the options in the
  condor submit file for a particular executable
  """
  def __init__(self, universe, executable, queue):
    """
    universe = the condor universe to run the job in.
    executable = the executable to run.
    queue = number of jobs to queue.
    """
    self.__universe = universe
    self.__executable = executable
    self.__queue = queue

    # These are set by methods in the class
    self.__arguments = {}
    self.__notification = None
    self.__log_file = None
    self.__err_file = None
    self.__out_file = None
    self.__sub_file_path = None

  def add_arg(self, arg, value):
    """
    Add a command line argument to the executable.
    arg = command line argument to add.
    value = value to pass to the argument (None for no argument).
    """
    self.__arguments[arg] = value

  def add_ini_args(self, cp, section):
    """
    Parse command line arguments from a given section in an ini file and
    pass to the executable.
    cp = ConfigParser object pointing to the ini file.
    section = section of the ini file to add to the arguments.
    """
    for opt in cp.options(section):
      arg = string.strip(cp.get(section,opt))
      self.__arguments[opt] = arg

  def set_notifcation(self, value):
    """
    Set the email address to send notification to.
    value = email address or never for no notification.
    """
    self.__notification = value

  def set_log_file(self, path):
    """
    Set the Condor log file.
    path = path to log file.
    """
    self.__log_file = path
    
  def set_stderr_file(self, path):
    """
    Set the file to which Condor directs the stderr of the job.
    path = path to stderr file.
    """
    self.__err_file = path

  def get_stderr_file(self):
    """
    Get the file to which Condor directs the stderr of the job.
    """
    return self.__err_file

  def set_stdout_file(self, path):
    """
    Set the file to which Condor directs the stdout of the job.
    path = path to stdout file.
    """
    self.__out_file = path
  
  def get_stdout_file(self):
    """
    Get the file to which Condor directs the stdout of the job.
    """
    return self.__out_file

  def set_sub_file(self, path):
    """
    Set the name of the file to write the Condor submit file to when
    write_sub_file() is called.
    path = path to submit file.
    """
    self.__sub_file_path = path

  def get_sub_file(self):
    """
    Get the name of the file which the Condor submit file will be
    written to when write_sub_file() is called.
    path = path to submit file.
    """
    return self.__sub_file_path

  def write_sub_file(self):
    """
    Write a submit file for this Condor job.
    """
    if not self.__log_file:
      raise CondorSubmitError, "Log file not specified."
    if not self.__err_file:
      raise CondorSubmitError, "Error file not specified."
    if not self.__out_file:
      raise CondorSubmitError, "Output file not specified."
  
    if not self.__sub_file_path:
      raise CondorSubmitError, 'No path for submit file.'
    try:
      subfile = open(self.__sub_file_path, 'w')
    except:
      raise CondorSubmitError, "Cannot open file " + self.__sub_file_path

    subfile.write( 'universe = ' + self.__universe + '\n' )
    subfile.write( 'executable = ' + self.__executable + '\n' )

    if self.__arguments.keys():
      subfile.write( 'arguments =' )
    for c in self.__arguments.keys():
      if self.__arguments[c]:
        subfile.write( ' --' + c + ' ' + self.__arguments[c] )
      else:
        subfile.write( ' --' + c )
    subfile.write( '\n' )

    subfile.write( 'log = ' + self.__log_file + '\n' )
    subfile.write( 'error = ' + self.__err_file + '\n' )
    subfile.write( 'output = ' + self.__out_file + '\n' )
    if self.__notification:
      subfile.write( 'notification = ' + self.__notifcation + '\n' )
    subfile.write( 'queue ' + str(self.__queue) + '\n' )

    subfile.close()



class CondorDAGJob(CondorJob):
  """
  A Condor DAG job never notifies the user on completion and can have variable
  arguments that are set for a particular node in the DAG. Inherits methods
  from a CondorJob.
  """
  def __init__(self, universe, executable):
    """
    universe = the condor universe to run the job in.
    executable = the executable to run in the DAG.
    """
    CondorJob.__init__(self, universe, executable, 1)
    self.__notifcation = 'never'
    self.__var_args = []

  def add_var_arg(self, arg):
    """
    Add a variable (or macro) option to the condor job. The option is added 
    to the submit file and a different argument to the option can be set fot
    each node in the DAG.
    arg = name of option to add.
    """
    if arg not in self.__var_args:
      self.__var_args.append(arg)
      self.add_arg(arg,'$(' + arg + ')')



class CondorDAGNode:
  """
  A CondorDAGNode represents a node in the DAG. It corresponds to a particular
  condor job (and so a particular submit file). If the job has variable
  (macro) arguments, they can be set here so each nodes executes with the
  correct arguments.
  """
  def __init__(self, job):
    """
    job = the CondorJob that this node corresponds to.
    """
    if not isinstance(job, CondorDAGJob):
      raise CondorDAGNodeError, "A DAG node must correspond to a Condor DAG job"
    self.__name = None
    self.__job = job
    self.__vars = {}
    self.__retry = 0
    self.__parents = []
    self.set_name()

  def __repr__(self):
    return self.__name

  def job(self):
    """
    Return the CondorJob that this node is associated with.
    """
    return self.__job
  
  def set_name(self):
    """
    Generate a unique name for this node in the DAG.
    """
    t = str( long( time.time() * 1000 ) )
    r = str( long( random.random() * 100000000000000000L ) )
    a = str( self.__class__ )
    self.__name = md5.md5(t + r + a).hexdigest()
    
  def add_var(self,var,value):
    """
    Add the a variable (macro) arguments for this node. If the option
    specified does not exist in the CondorJob, it is added so the submit
    file will be correct when written.
    var = option name.
    value = value of the option for this node in the DAG.
    """
    self.__vars[var] = value
    self.__job.add_var_arg(var)

  def set_retry(self, retry):
    """
    Set the number of times that this node in the DAG should retry.
    retry = number of times to retry node.
    """
    self.__retry = retry

  def write_job(self,fh):
    """
    Write the DAG entry for this node's job to the DAG file descriptor.
    fh = descriptor of open DAG file.
    """
    fh.write( 'JOB ' + self.__name + ' ' + self.__job.get_sub_file() +  '\n' )
    fh.write( 'RETRY ' + self.__name + ' ' + str(self.__retry) + '\n' )

  def write_vars(self,fh):
    """
    Write the variable (macro) arguments to the DAG file descriptor.
    fh = descriptor of open DAG file.
    """
    if self.__vars.keys():
      fh.write( 'VARS ' + self.__name )
      for k in self.__vars.keys():
        fh.write( ' ' + str(k) + '="' + str(self.__vars[k]) + '"' )
      fh.write( '\n' )

  def write_parentss(self,fh):
    """
    Write the parent/child relations for this job to the DAG file descriptor.
    fh = descriptor of open DAG file.
    """
    for parent in self.__parents:
      fh.write( 'PARENT ' + parent + ' CHILD ' + str(self) + '\n' )

  def set_log_file(self,log):
    """
    Set the Condor log file to be used by this CondorJob.
    log = path of Condor log file.
    """
    self.__job.set_log_file(log)

  def add_parent(self,node):
    """
    Add a parent to this node. This node will not be executed until the
    parent node has run sucessfully.
    node = CondorDAGNode to add as a parent.
    """
    if not isinstance(node, CondorDAGNode):
      raise CondorDAGNodeError, "Parent must be a Condor DAG node"
    self.__parents.append( str(node) )
    


class CondorDAG:
  """
  A CondorDAG is a Condor Directed Acyclic Graph that describes a collection
  of Condor jobs and the order in which to run them. All Condor jobs in the
  DAG must write their Codor logs to the same file.
  NOTE: The log file must not be on an NFS mounted system as the Condor jobs
  must be able to get an exclusive file lock on the log file.
  """
  def __init__(self,log):
    """
    log = path to log file which must not be on an NFS mounted file system.
    """
    self.__log_file_path = log
    self.__dag_file_path = None
    self.__jobs = []
    self.__nodes = []

  def set_dag_file(self, path):
    """
    Set the name of the file into which the DAG is written.
    path = path to DAG file.
    """
    self.__dag_file_path = path

  def add_node(self,node):
    """
    Add a CondorDAGNode to this DAG. The CondorJob that the node uses is 
    also added to the list of Condor jobs in the DAG so that a list of the
    submit files needed by the DAG can be maintained. Each unique CondorJob
    will be added once to prevent duplicate submit files being written.
    node = CondorDAGNode to add to the CondorDAG.
    """
    if not isinstance(node, CondorDAGNode):
      raise CondorDAGError, "Nodes must be class CondorDAGNode or subclass"
    node.set_log_file(self.__log_file_path)
    self.__nodes.append(node)
    if node.job() not in self.__jobs:
      self.__jobs.append(node.job())

  def write_sub_files(self):
    """
    Write all the submit files used by the dag to disk. Each submit file is
    written to the file name set in the CondorJob.
    """
    for job in self.__jobs:
      job.write_sub_file()

  def write_dag(self):
    """
    Write all the nodes in the DAG to the DAG file.
    """
    if not self.__log_file_path:
      raise CondorDAGError, "No path for DAG file"
    try:
      dagfile = open( self.__dag_file_path, 'w' )
    except:
      raise CondorDAGError, "Cannot open file " + self.__dag_file_path
    for node in self.__nodes:
      node.write_job(dagfile)
    for node in self.__nodes:
      node.write_vars(dagfile)
    for node in self.__nodes:
      node.write_parents(dagfile)
    dagfile.close()



class AnalysisJob:
  """
  Describes a generic analysis job that filters LIGO data as configured by
  an ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object that contains the configuration for this job.
    """
    self.__cp = cp
    self.__channel = string.strip(self.__cp.get('input','channel'))
  def get_config(self,sec,opt):
    """
    Get the configration variable in a particular section of this jobs ini
    file.
    sec = ini file section.
    opt = option from section sec.
    """
    return string.strip(self.__cp.get(sec,opt))
  def channel(self):
    """
    Returns the name of the channel that this job is filtering. Note that 
    channel is defined to be IFO independent, so this may be LSC-AS_Q or
    IOO-MC_F. The IFO is set on a per node basis, not a per job basis.
    """
    return self.__channel
  def calibration(self,ifo):
    """
    Returns the name of the calibration file to use for the given IFO.
    ifo = name of interferomener (e.g. L1, H1 or H2).
    """
    cal_path = string.strip(self.__cp.get('calibration','path'))
    cal_file = string.strip(self.__cp.get('calibration',ifo))
    cal = os.path.join(cal_path,cal_file)
    return cal    



class AnalysisNode(CondorDAGNode):
  """
  Contains the methods that allow an object to be built to analyse LIGO
  data in a Condor DAG.
  """
  def __init__(self):
    self.__start = 0
    self.__end = 0
    self.__ifo = None
    self.__input = None
    self.__output = None

  def set_start(self,time):
    """
    Set the GPS start time of the analysis node by setting a --gps-start-time
    option to the node when it is executed.
    time = GPS start time of job.
    """
    self.add_var('gps-start-time',time)
    self.__start = time

  def get_start(self):
    """
    Get the GPS start time of the node.
    """
    return self.__start
    
  def set_end(self,time):
    """
    Set the GPS end time of the analysis node by setting a --gps-end-time
    option to the node when it is executed.
    time = GPS end time of job.
    """
    self.add_var('gps-end-time',time)
    self.__end = time

  def get_end(self):
    """
    Get the GPS end time of the node.
    """
    return self.__end

  def set_input(self,file):
    """
    Add an input to the node by adding a --input option.
    file = option argument to pass as input.
    """
    self.__input = file
    self.add_var('input', file)

  def get_input(self):
    """
    Get the file that will be passed as input.
    """
    return self.__input

  def set_output(self, file):
    """
    Add an output to the node by adding a --output option.
    file = option argument to pass as output.
    """
    self.__output = file
    self.add_var('output', file)

  def get_output(self):
    """
    Get the file that will be passed as output.
    """
    return self.__output

  def set_ifo(self,ifo):
    """
    Set the channel name to analyze and add a calibration file for that
    channel. The name of the ifo is prepended to the channel name obtained
    from the job configuration file and passed with a --channel-name option.
    A calibration file is obtained from the ini file and passed with a 
    --calibration-cache option.
    ifo = two letter ifo code (e.g. L1, H1 or H2).
    """
    self.__ifo = ifo
    self.add_var('channel-name', ifo + ':' + self.job().channel())
    self.add_var('calibration-cache', self.job().calibration(ifo))

  def get_ifo(self):
    """
    Returns the two letter IFO code for this node.
    """
    return self.__ifo

  def set_cache(self,file):
    """
    Set the LAL frame cache to to use. The frame cache is passed to the job
    with the --frame-cache argument.
    file = calibration file to use.
    """
    self.add_var('frame-cache', file)




class AnalysisChunk:
  """
  An AnalysisCunk is the unit of data that a node works with, usually some
  subset of a ScienceSegment.
  """
  def __init__(self, start, end):
    """
    start = GPS start time of the chunk.
    end = GPS end time of the chunk.
    """
    self.__start = start
    self.__end = end
    self.__length = end - start

  def __repr__(self):
    return '<AnalysisChunk: start %d, end %d>' % (self.start(), self.end())
    
  def start(self):
    """
    Returns the GPS start time of the chunk.
    """
    return self.__start

  def end(self):
    """
    Returns the GPS end time of the chunk.
    """
    return self.__end
    
  def dur(self):
    """
    Returns the length (duration) of the chunk in seconds.
    """
    return self.__length


class ScienceSegment:
  """
  A ScienceSegment is a period of time where the experimenters determine
  that the inteferometer is in a state where the data is suitable for 
  scientific analysis. A science segment can have a list of AnalysisChunks
  asscociated with it that break the segment up into (possibly overlapping)
  smaller time intervals for analysis.
  """
  def __init__(self,segment):
    """
    segemnt = a tuple containing the (segment id, gps start time, gps end
    time, duration) of the segment.
    """
    self.__segment = segment
    self.__chunks = []
    self.__unused = self.dur()
    self.__ifo = None

  def __getitem__(self,i):
    """
    Allows iteration over and direct access to the AnalysisChunks contained
    in this ScienceSegment.
    """
    if i < 0: raise IndexError, "list index out of range"
    return self.__chunks[i]
    
  def __len__(self):
    """
    Returns the number of AnalysisChunks contained in this ScienceSegment.
    """
    return len(self.__chunks)

  def __repr__(self):
    return '<ScienceSegment: id %d, start %d, end %d, duration %d>' % (
    self.id(),self.start(),self.end(),self.dur())

  def make_chunks(self,length=0,overlap=0,play=0):
    """
    Divides the science segment into chunks of length seconds overlapped by
    overlap seconds. If the play option is set, only chunks that contain S2
    playground data are generated. If the user has a more complicated way
    of generating chunks, this method should be overriden in a sub-class.
    Any data at the end of the ScienceSegment that is too short to contain a 
    chunk is ignored. The length of this unused data is stored and can be
    retrieved with the unused() method.
    length = length of chunk in seconds.
    overlap = overlap between chunks in seconds.
    play = only generate chunks that overlap with S2 playground data.
    """
    time_left = self.dur()
    start = self.start()
    increment = length - overlap
    while time_left >= length:
      middle = start + length / 2
      end = start + length
      if (not play) or ( play 
       and ( s2play(start) or s2play(middle) or s2play(end) ) ):
        self.__chunks.append(AnalysisChunk(start,end))
      start += increment
      time_left -= increment
    self.__unused = time_left

  def add_chunk(self,start,end):
    """
    Add an AnalysisChunk to the list associated with this ScienceSegment.
    start = GPS start time of chunk.
    end = GPS end time of chunk.
    """
    self.__analysis_segs.append(AnalysisChunk(start,end))

  def unused(self):
    """
    Returns the length of data in the science segment not used to make chunks.
    """
    return self.__unused

  def id(self):
    """
    Returns the ID of this ScienceSegment.
    """
    return self.__segment[0]
    
  def start(self):
    """
    Returns the GPS start time of this ScienceSegment.
    """
    return self.__segment[1]

  def end(self):
    """
    Returns the GPS end time of this ScienceSegment.
    """
    return self.__segment[2]

  def dur(self):
    """
    Returns the length (duration) in seconds of this ScienceSegment.
    """
    return self.__segment[3]


class ScienceData:
  """
  An object that can contain all the science data used in an analysis. Can
  contain multiple ScienceSegments and has a method to generate these from
  a text file produces by the LIGOtools segwizard program.
  """
  def __init__(self):
    self.__sci_segs = []
    self.__sci_times = []
    self.__file = None

  def __getitem__(self,i):
    """
    Allows direct access to or iteration over the ScienceSegments associated
    with the ScienceData.
    """
    return self.__sci_segs[i]

  def __repr__(self):
    return '<ScienceData: file %s>' % self.__file

  def __len__(self):
    """
    Returns the number of ScienceSegments associated with the ScienceData.
    """
    return len(self.__sci_segs)

  def read(self,file):
    """
    Parse the science segments from the segwizard output contained in file.
    file = input text file containing a list of science segments generated by
    segwizard.
    """
    self.__file = file
    octothorpe = re.compile(r'\A#')
    for line in open(file):
      if not octothorpe.match(line):
        self.__sci_times.append(tuple(map(int,line.split())))
        x = ScienceSegment(tuple(map(int,line.split())))
        self.__sci_segs.append(x)

  def make_chunks(self,length,overlap,play):
    """
    Divide each ScienceSegment contained in this object into AnalysisChunks.
    length = length of chunk in seconds.
    overlap = overlap between segments.
    play = if true, only generate chunks that overlap with S2 playground data.
    """
    for seg in self.__sci_segs:
      seg.make_chunks(length,overlap,play)
