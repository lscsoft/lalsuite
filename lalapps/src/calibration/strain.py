"""
Classes needed for the cosmic string analysis pipeline.
"""

__author__ = 'Xavier Siemens<siemens@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string,sys
import exceptions
from glue import pipeline
from pylab import *
import operator

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

def plot_noise_jobs(filelist,cp,dir,epoch,dag,qjob):
  qfile = write_qscan_conf(dir,cp)
  web = open(dir+'/'+'index.html','w')
  filelist.sort()
  ifo = cp.get('pipeline','ifo')
  duration = cp.get('plot','duration')
  specList = []
  input = open(filelist[0],'r')
  start = input.readlines()[0].split()[0]
  input.close()
  time = []
  freq = []
  fignames = []
  timeFreqTuple = []
  freqfile = cp.get('noisecomp','freq-file')
  for line in open(freqfile,'r').readlines():
    freq.append(float(line.strip()))
  #STOP = 0;
  for file in filelist:
    try: input = open(file,'r')
    except: 
      print "WARNING: file " + file + " doesn't exist"
      continue
    #if STOP > 100: break
    #STOP+=1
    for line in input.readlines():
      tmp = line.split()
      if (float(tmp[0]) - float(start)) > float(duration):
        #fn,ft = plot_noise_spec(specList,time,freq,cp,dir,dag,qjob,qfile,timeFreqTuple)
        fn,ft = plot_noise_spec(specList,cp,dir,dag,qjob,qfile,timeFreqTuple)
        fignames.append([fn,ft])
        start = float(tmp[0])
        #time = []
        specList = []
        timeFreqTuple = []
      specCol = []
      #time.append(float(tmp[0]))
      #specCol.extend(float(tmp[ix]) for ix in range(3,len(tmp),3))
      #specList.append(specCol)
      for ix in range(0,len(freq)):
        #timeFreqTuple.append((tmp[0],freq[ix],specCol[ix]))
        timeFreqTuple.append((float(tmp[0]),float(freq[ix]),float(tmp[3+3*ix])))
        specCol.append(float(tmp[3+3*ix]))
      specList.append(specCol)
    input.close() 
  fignames.sort(reverse=True)
  web.write('<table>')
  web.write('<tr><td><b>IFO=</b>'+ifo+'</td></tr>\n')
  web.write('<tr><td><b>BAND=</b>'+cp.get('noisecomp','band')+'</td></tr>\n')
  web.write('<tr><td><b>TIME=</b>'+cp.get('noisecomp','time')+'</td></tr>\n')
  web.write('<tr><td><b>FCAL=</b>'+cp.get('noisecomp','fcal')+'</td></tr>\n')
  web.write('<tr><td><b>GAMMA-FUDGE-FACTOR=</b>'+cp.get('noisecomp','gamma-fudge-factor')+'</td></tr>\n')
  web.write('<tr><td><b>HOFT=</b>'+cp.get('noisecomp','hoft-channel')+'</td></tr>\n')
  web.write('<tr><td><b>ASQ=</b>'+cp.get('noisecomp','asq-channel')+'</td></tr>\n')
  web.write('<tr><td><b>EXC=</b>'+cp.get('noisecomp','exc-channel')+'</td></tr>\n')
  web.write('<tr><td><b>DARM=</b>'+cp.get('noisecomp','darm-channel')+'</td></tr>\n')
  web.write('<tr><td><b>DERR=</b>'+cp.get('noisecomp','derr-channel')+'</td></tr>\n')
  web.write('<tr><td><b>EPOCH=</b>'+epoch[0]+'</td></tr>\n')
  web.write('<tr><td><b>START=</b>'+epoch[1]+'</td></tr>\n')
  web.write('<tr><td><b>STOP=</b>'+epoch[2]+'</td></tr>\n')
  web.write('<tr><td><b>DURATION=</b>'+epoch[3]+'</td></tr>\n')
  web.write('<tr><td><b>OPEN LOOP GAIN (RE)=</b>'+epoch[4]+'</td></tr>\n')
  web.write('<tr><td><b>OPEN LOOP GAIN (IM)=</b>'+epoch[5]+'</td></tr>\n')
  web.write('<tr><td><b>SERVO (RE)=</b>'+epoch[6]+'</td></tr>\n')
  web.write('<tr><td><b>SERVO (IM)=</b>'+epoch[7]+'</td></tr>\n')
  web.write('<tr><td><b>WHITENER (RE)=</b>'+epoch[8]+'</td></tr>\n')
  web.write('<tr><td><b>WHITENER (IM)=</b>'+epoch[9]+'</td></tr>\n')
  web.write('<tr><td><b>OLG FILE=</b>'+epoch[11]+'</td></tr>\n')
  web.write('<tr><td><b>SENSING FILE=</b>'+epoch[12]+'</td></tr>\n')
  web.write('</table><br><br>')
  web.write('<h1>SPEC GRAMS, starting with the bad stuff between '+str(freq[0])+' and '+str(freq[-1])+' Hz</h1>')
  web.write('<table><tr>\n')
  cnter = 1
  for fig in fignames:
    if cnter == 6:
      web.write('</tr><tr>')
      cnter = 0
    web.write('<td><a href='+fig[0]+'><img src=thumb-'+fig[0]+'></a>\n')
    web.write('<br><a href='+str(fig[1])+'>'+str(fig[1])+'</a></td>\n')
    cnter += 1
  web.write('</tr></table>\n')
  web.close()

def plot_noise_spec(specList,cp,dir,dag,qjob,qfile,tftuple):
  fignames = []
  time = list(set(map(operator.itemgetter(0),tftuple)))
  time.sort()
  freq = list(set(map(operator.itemgetter(1),tftuple)))
  freq.sort()
  #print time, freq
  freq = array(freq,typecode='f')
  Time = array(time,typecode='d') - time[0]
  X,Y = meshgrid(Time,freq)
  start = str(time[0])
  end = str(time[-1])
  flat_specList = []
#  for item in specList:
#    flat_specList.extend(item)
#  MIN = min(flat_specList) 
#  MAX = max(flat_specList)
  tftuple.sort(key=operator.itemgetter(2))
  MIN = tftuple[0][2]
  MINTime = tftuple[0][0]
  tftuple.sort(key=operator.itemgetter(2),reverse=True)
  MAX = tftuple[0][2]
  MAXTime = tftuple[0][0]
  OUTLIER = [1-MIN, MAX-1]
  if (1-MIN) > (MAX-1):
    qscanTime = MINTime
    # flat_specList.index(MIN) / len(freq) * 60 + int(float(start)) +32
  else:
    qscanTime = MAXTime
    # flat_specList.index(MAX) / len(freq) * 60 + int(float(start)) +32
  dag.add_node(qscanNode(qjob,qscanTime,qfile,cp.get('pipeline','ifo'),dir,OUTLIER))
  figname = str(max(OUTLIER))+'-'+dir + '-' + start + '-' + end + '.png'
  A = array(specList,typecode='f')
  figure(1)
  pcolor(X,Y,A.transpose(),shading='flat',vmin=0.95,vmax=1.05)
  title('h(t) and h(f) power ratios per freq bin GPS '+start + '\n min = '+str(MIN) + ' max = '+str(MAX) )
  xlabel('Time')
  ylabel('Frequency')
  colorbar()
  savefig(dir + '/'+ figname)
  thumb = 'thumb-'+figname
  savefig(dir + '/'+ thumb,dpi=20)
  clf()
  close()
  return figname,qscanTime
  
def write_qscan_conf(epoch,cp):
  duration = cp.get('noisecomp','time')
  ifo = cp.get('pipeline','ifo')
  qfilename = str(ifo)+str(epoch)+duration+'.conf'
  qfile = open(str(ifo)+str(epoch)+duration+'.conf','w')
  qfile.write('[Context,Context]\n\n[Parameters,Parameter Estimation]\n\n[Notes,Notes]')
  qfile.write('[Gravitational,Gravitational wave data]\n{\n')
  qfile.write("channelName:\t'"+cp.get('noisecomp','hoft-channel')+"'\n")
  qfile.write("frameType:\t'"+cp.get('input','type-hoft')+"'\n")
  qfile.write("sampleFrequency:\t4096\n")
  qfile.write("searchTimeRange:\t64\n")
  qfile.write("searchFrequencyRange:\t[32 Inf]\n")
  qfile.write("searchQRange:\t[4 64]\n")
  qfile.write("searchMaximumEnergyLoss:\t0.2\n")
  qfile.write("whiteNoiseFalseRate:\t1e2\n")
  qfile.write("searchWindowDuration:\t0.5\n")
  qfile.write("plotTimeRanges:\t[1 4 "+duration+"]\n")
  qfile.write("plotFrequencyRange:\t[]\n")
  qfile.write("plotMaximumEnergyLoss:\t0.2\n")
  qfile.write("plotNormalizedEnergyRange:\t[0 25.5]\n}\n\n")
  qfile.write('{\n')
  qfile.write("channelName:\t'"+cp.get('noisecomp','derr-channel')+"'\n")
  qfile.write("frameType:\t'"+cp.get('input','type-derr')+"'\n")
  qfile.write("searchTimeRange:\t64\n")
  qfile.write("searchFrequencyRange:\t[32 Inf]\n")
  qfile.write("searchQRange:\t[4 64]\n")
  qfile.write("searchMaximumEnergyLoss:\t0.2\n")
  qfile.write("whiteNoiseFalseRate:\t1e2\n")
  qfile.write("searchWindowDuration:\t0.5\n")
  qfile.write("plotTimeRanges:\t[1 4 "+duration+"]\n")
  qfile.write("plotFrequencyRange:\t[]\n")
  qfile.write("plotMaximumEnergyLoss:\t0.2\n")
  qfile.write("plotNormalizedEnergyRange:\t[0 25.5]\n}\n")
  qfile.close()
  return qfilename

class qscanJob(pipeline.CondorDAGJob):
  """
  A qscan job
  """
  def __init__(self, cp, tag_base='QSCAN'):
    """
    """
    self.__executable = string.strip(cp.get('condor','qscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.add_condor_cmd('getenv','True')
    self.set_stdout_file('logs/qscan-$(macrochannelname)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/qscan-$(macrochannelname)-$(cluster)-$(process).err')
    self.set_sub_file('qscan.sub')



class qscanNode(pipeline.CondorDAGNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,qfile,ifo,dir,OUTLIER):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    self.id = ifo + '-qscan-' + repr(time)

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_arg(repr(time))
    if max(OUTLIER) < 0.10:
      # just look at darm and h(t) for puny outliers. 
      self.add_file_arg(qfile)
    else:
      print ".....found 10% outlier running full qscan\n"
      # run the standard qscan on outliers greater than 10%
      self.add_var_arg('@default')
    self.add_var_arg('@default')
    self.add_var_arg(dir)
    
