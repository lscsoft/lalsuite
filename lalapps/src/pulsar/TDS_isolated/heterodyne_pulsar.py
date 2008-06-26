"""
  Classes needed for the time domain pulsar heterodyne pipeline.

  This has been greatly inspired (well blatently hacked together from) the stochastic pipeline by
  Adam Mercer and the frequency domain binary pulsar search by Chris Messenger.
  
  Matt Pitkin 30/03/06
  
  $Id$
"""

__author__ = 'Matt Pitkin <matthew@astro.gla.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

import string
import exceptions
import os
from glue import pipeline

class binarypulsarError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args
    
# And begin ...

class heterodyneJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_heterodyne_pulsar job to coarse heterodyne the data.
  """
  def __init__(self,cp):
    # cp = ConfigParser object from which options are read.
    self.__executable = cp.get('exec', 'heterodyne') # gets the executable from .ini
    self.__universe = cp.get('condor', 'universe')
    pipeline.CondorDAGJob.__init__(self, self.__universe, self.__executable)
    pipeline.AnalysisJob.__init__(self, cp)
    
    self.add_condor_cmd('USE_NFS','True')
    self.add_condor_cmd('WantRemoteIO','false')
    self.add_condor_cmd('local_files','*')
    
    # set log files for job
    self.set_stdout_file('logs/heterodyne_pulsar-$(cluster).out')
    self.set_stderr_file('logs/heterodyne_pulsar-$(cluster).err')
    self.set_sub_file('heterodyne_pulsar.sub')

class heterodyneNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A heterodyneNode runs an instance of lalapps_heterodyne_pulsar in coarse heterodyne mode in
  a condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_CalculateSensitivity.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    
    # initilise job variables
    self.__ifo = None
    self.__param_file = None
    self.__freq_factor = None
    self.__filter_knee = None
    self.__sample_rate = None
    self.__resample_rate = None
    self.__data_file = None
    self.__channel = None
    self.__seg_list = None
    self.__data_file = None
    self.__output_dir = None
    self.__het_flag = None
    self.__pulsar = None
    self.__high_pass = None
    self.__scale_fac = None
    self.__manual_epoch = None
      
  def set_data_file(self,data_file):
    # set file containing data to be heterodyne (either list of frames or coarse het output)
    self.add_var_opt('data-file',data_file)
    self.__data_file = data_file
    
  def set_output_dir(self,output_dir):
    # set output directory
    self.add_var_opt('output-dir',output_dir)
    self.__output_dir = output_dir
    
  def set_seg_list(self,seg_list):
    # set segment list to be used
    self.add_var_opt('seg-file',seg_list)
    self.__seg_list = seg_list
    
  def set_ifo(self, ifo):
    # set detector
    self.add_var_opt('ifo', ifo)
    self.__ifo = ifo
    
  def set_param_file(self, param_file):
    # set pulsar parameter file
    self.add_var_opt('param-file', param_file)
    self.__param_file = param_file
  
  def set_freq_factor(self, freq_factor):
    # set the factor by which to muliply the pulsar spin frequency (normally 2.0)
    self.add_var_opt('freq-factor', freq_factor)
    self.__freq_factor = freq_factor
    
  def set_param_file_update(self, param_file_update):
    # set file containing updated pulsar parameters
    self.add_var_opt('param-file-update',param_file_update)
    
  def set_manual_epoch(self, manual_epoch):
      # set manual pulsar epoch
      self.add_var_opt('manual-epoch',manual_epoch)
      self.__manual_epoch = manual_epoch
    
  def set_ephem_earth_file(self, ephem_earth_file):
    # set the file containing the earth's ephemeris
    self.add_var_opt('ephem-earth-file', ephem_earth_file)
    
  def set_ephem_sun_file(self, ephem_sun_file):
    # set the file containing the sun's ephemeris
    self.add_var_opt('ephem-sun-file', ephem_sun_file)

  def set_pulsar(self,pulsar):
    # set pulsar name
    self.add_var_opt('pulsar',pulsar)
    self.__pulsar = pulsar
    
  def set_het_flag(self,het_flag):
    # set heterodyne flag
    self.add_var_opt('heterodyne-flag',het_flag)
    self.__het_flag = het_flag
    
  def set_filter_knee(self,filter_knee):
    # set filter knee frequency
    self.add_var_opt('filter-knee',filter_knee)
    self.__filter_knee = filter_knee
  
  def set_channel(self,channel):
    # set channel containing data from frames
    self.add_var_opt('channel',channel)
    self.__channel = channel
    
  def set_sample_rate(self,sample_rate):
    # set sample rate of input data
    self.add_var_opt('sample-rate',sample_rate)
    self.__sample_rate = sample_rate
    
  def set_resample_rate(self,resample_rate):
    # set resample rate for output data
    self.add_var_opt('resample-rate',resample_rate)
    self.__resample_rate = resample_rate
    
  def set_stddev_thresh(self,stddev_thresh):
    # set standard deviation threshold at which to remove outliers
    self.add_var_opt('stddev-thresh',stddev_thresh)
    
  def set_calibrate(self):
    # set calibration flag
    self.add_var_opt('calibrate', '') # no variable required
    
  def set_verbose(self):
    # set verbose flag
    self.add_var_opt('verbose', '') # no variable required
    
  def set_bininput(self):
    # set binary input file flag
    self.add_var_opt('binary-input', '') # no variable required
    
  def set_binoutput(self):
    # set binary output file flag
    self.add_var_opt('binary-output', '') # no variable required
    
  def set_response_function(self,response_function):
    # set reponse function file
    self.add_var_opt('response-file',response_function)
    
  def set_coefficient_file(self,coefficient_file):
    # set the file containing the calibration coefficients (e.g alpha and gammas)
    self.add_var_opt('coefficient-file',coefficient_file)
    
  def set_sensing_function(self,sensing_function):
    # set file containing the sensing function for calibration
    self.add_var_opt('sensing-function',sensing_function)
    
  def set_open_loop_gain(self,open_loop_gain):
    # set file containing the open loop gain for calibration
    self.add_var_opt('open-loop-gain',open_loop_gain)

  def set_scale_fac(self,scale_fac):
    # set scale factor for calibrated data
    self.add_var_opt('scale-factor',scale_fac)

  def set_high_pass(self,high_pass):
    # set high-pass frequency for calibrated data
    self.add_var_opt('high-pass-freq',high_pass)

# taken from pipeline.py to add a log file to the .sub file  
class LSCDataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  An LSCdataFind job used to locate data. The static options are
  read from the section [datafind] in the ini file. The stdout from
  LSCdataFind contains the paths to the frame files and is directed to a file
  in the cache directory named by site and GPS start and end times. The stderr
  is directed to the logs directory. The job always runs in the scheduler
  universe. The path to the executable is determined from the ini file.
  """
  def __init__(self,cache_dir,log_dir,config_file,dax=0):
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
    self.__dax = dax
    self.__config_file = config_file

    # we have to do this manually for backwards compatibility with type
    for o in self.__config_file.options('datafind'):
      opt = string.strip(o)
      if opt is not "type":
        arg = string.strip(self.__config_file.get('datafind',opt))
        self.add_opt(opt,arg)

    if self.__dax:
      # only get the LFNs not the PFNs
      self.add_opt('names-only','')
    else:
      # we need a lal cache for file PFNs
      self.add_opt('lal-cache','')
      self.add_opt('url-type','file')

    self.add_condor_cmd('getenv','True')

    self.set_stderr_file(os.path.join(log_dir,\
'datafind-$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(\
process).err'))
    self.set_stdout_file(os.path.join(self.__cache_dir,\
'$(macroobservatory)-$(macrotype)-$(macrogpsstarttime)-$(macrogpsendtime).cache'))
    self.set_sub_file('datafind.sub')

  def get_cache_dir(self):
    """
    returns the directroy that the cache files are written to.
    """
    return self.__cache_dir

  def is_dax(self):      
    """          
    returns the dax flag         
    """          
    return self.__dax    

  def get_config_file(self):
    """
    return the configuration file object
    """
    return self.__config_file
  

class LSCDataFindNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  A DataFindNode runs an instance of LSCdataFind in a Condor DAG.
  """
  def __init__(self,job):
    """
    @param job: A CondorDAGJob that can run an instance of LALdataFind.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
    self.__start = 0
    self.__end = 0
    self.__observatory = None
    self.__output = None
    self.__job = job
    self.__dax = job.is_dax()    
    self.__lfn_list = None
     
    # try and get a type from the ini file and default to type None
    try:
      self.set_type(self.job().get_config_file().get('datafind','type'))
    except:
      self.__type = None
 
  def __set_output(self):
    """
    Private method to set the file to write the cache to. Automaticaly set
    once the ifo, start and end times have been set.
    """
    if self.__start and self.__end and self.__observatory and self.__type:
      self.__output = os.path.join(self.__job.get_cache_dir(), self.__observatory + '-' +
self.__type + '-' + str(self.__start) + '-' + str(self.__end) + '.cache')
      self.add_output_file(self.__output)

  def set_start(self,time):
    """
    Set the start time of the datafind query.
    @param time: GPS start time of query.
    """
    self.add_var_opt('gps-start-time', time)
    self.__start = time
    self.__set_output()

  def get_start(self):
    """
    Return the start time of the datafind query
    """
    return self.__start

  def set_end(self,time):
    """
    Set the end time of the datafind query.
    @param time: GPS end time of query.
    """
    self.add_var_opt('gps-end-time', time)
    self.__end = time
    self.__set_output()

  def get_end(self):
    """
    Return the start time of the datafind query
    """
    return self.__end

  def set_observatory(self,obs):
    """
    Set the IFO to retrieve data for. Since the data from both Hanford 
    interferometers is stored in the same frame file, this takes the first 
    letter of the IFO (e.g. L or H) and passes it to the --observatory option
    of LSCdataFind.
    @param obs: IFO to obtain data for.
    """
    self.add_var_opt('observatory',obs)
    self.__observatory = str(obs)
    self.__set_output()

  def get_observatory(self):
    """
    Return the start time of the datafind query
    """
    return self.__observatory

  def set_type(self,type):
    """
    sets the frame type that we are querying
    """
    self.add_var_opt('type',str(type))
    self.__type = str(type)
    self.__set_output()

  # add a set server function
  def set_server(self,server):
    # set the server manually rather than from the environment variable
    self.add_var_opt('server',server)
  
  def get_type(self):
    """
    gets the frame type that we are querying
    """
    return self.__type

  def get_output(self):
    """
    Return the output file, i.e. the file containing the frame cache data.
    or the files itself as tuple (for DAX)       
    """  
    if self.__dax:
      if not self.__lfn_list:
        # call the datafind client to get the LFNs
        from pyGlobus import security
        from glue import LDRdataFindClient
        from glue import gsiserverutils

        hostPortString = os.environ['LSC_DATAFIND_SERVER']
        print >>sys.stderr, ".",
        if hostPortString.find(':') < 0:
          # no port specified
          myClient = LDRdataFindClient.LSCdataFindClient(hostPortString)
        else:
          # server and port specified
          host, portString = hostPortString.split(':')
          port = int(portString)
          myClient = LDRdataFindClient.LSCdataFindClient(host,port)

        clientMethod = 'findFrameNames'
        clientMethodArgDict = {
          'observatory': self.get_observatory(),
          'end': str(self.get_end()),
          'start': str(self.get_start()),
          'type': self.get_type(),
          'filename': None,
          'urlType': None,
          'match': None,
          'limit': None,
          'offset': None,
          'strict' : None,
          'namesOnly' : True
          }

        print >>sys.stderr, ".",
        time.sleep( 1 )
        result = eval("myClient.%s(%s)" % (clientMethod, clientMethodArgDict))
        
        if not isinstance(result,LDRdataFindClient.lfnlist):
          msg = "datafind server did not return LFN list : " + str(result)
          raise SegmentError, msg
        if len(result) == 0:
          msg = "No LFNs returned for segment %s %s" % ( str(self.get_start()),
            str(self.get_end()) )
          raise SegmentError, msg
        self.__lfn_list = result
      return self.__lfn_list
    else:        
      return self.__output
