# DAG Class definitions for LALInference Pipeline
# (C) 2012 John Veitch, Kiersten Ruisard, Kan Wang

import glue
from glue import pipeline
import os

# We use the GLUE pipeline utilities to construct classes for each
# type of job. Each class has inputs and outputs, which are used to
# join together types of jobs into a DAG.

dummyCacheNames=['LALLIGO','LALVirgo','LALAdLIGO']

def chooseEngineNode(name):
  if name=='lalinferencenest':
    return LALInferenceNestNode
  if name=='lalinferencemcmc':
    return LALInferenceMCMCNode
  return EngineNode

class LALInferencePipelineDAG(pipeline.CondorDAG):
  def __init__(self,log,cp,dax=False):
    pipeline.CondorDAG.__init__(self,log,dax)
    self.subfiles=[]
    self.config=cp
    self.engine=cp.get_option('analysis','engine')
    self.EngineNode=chooseEngineNode(self.engine)
    if cp.has_option('paths','basedir'):
      self.basepath=cp.get_option('paths','basedir')
    else:
      self.basepath=os.getcwd()
      print 'No basepath specified, using current directory: %s'%(self.basepath)
    if cp.has_option('paths','cachedir'):
      self.cachepath=cp.get_option('paths','cachedir')
    else:
      self.cachepath=os.path.join(self.basepath,'caches')
    if cp.has_option('paths','logdir'):
      self.logpath=cp.get_option('paths','logdir')
    else:
      self.logpath=os.path.join(self.basepath,'log')
    if cp.has_option('analysis','ifos'):
      self.ifos=cp.get_option('analysis','ifos')
    else:
      self.ifos=['H1','L1','V1']
    self.segments={}
    for ifo in ifos:
      self.segments[ifo]=[]
    self.dq={}
    self.frtypes=cp.get_option('datafind','types')
    self.use_available_data=False
    # Set up necessary job files.
    self.datafind_job = pipeline.LSCDataFindJob(self.cachepath,self.logpath,self.config)
    self.datafind_job.add_opt('url-type','file')
    self.datafind_job.set_sub_file(os.path.join(self.basepath,'datafind.sub'))
    self.engine_job = EngineJob(self.config, os.path.join(self.basepath,'lalinference.sub'),self.logpath)
  
  def add_full_analysis_time(self,gpstime):
    """
    Analyse a given GPS time
    """
    datafindnode=self.add_datafind(gpstime)
    enginenode=self.add_engine(gpstime)
    self.add_resultspage(enginenode.posteriorfile)
  
  def get_science_segment(self,ifo,gpstime):
    # Check if time is in existing segment
    for seg in self.segments[ifo]:
      if gpstime in seg: return seg
    raise pipeline.CondorDAGError('Unable to find time in segments')
  
  def add_science_segments(self)
    # Query the segment database for science segments and
    # add them to the pool of segments
    for ifo in self.ifos:
      segFileName=inspiralutils.findSegmentsToAnalyze(self.config, ifo, self.veto_categories, generate_segments=True,\
	use_available_data=self.use_available_data , data_quality_vetoes=False)
      segfile=open(segFileName)
      (segs,self.dq[ifo])=segmentsUtils.fromsegwizard(segfile)
      segs.coalesce()
      segfile.close()
      for seg in segs:
	sciseg=pipeline.ScienceSegment((segs.index(seg),seg[0],seg[1],seg[1]-seg[0]))
	df_node=self.get_datafind_node(ifo,self.frtypes[ifo],int(sciseg.start()),int(sciseg.end()))
	sciseg.set_df_node(df_node)
	self.segments[ifo].append(sciseg)
  
  def get_datafind_node(self,ifo,frtype,gpsstart,gpsend):
    node=pipeline.LSCDataFindNode(self.datafind_job)
    node.set_observatory(ifo[0])
    node.set_type(frtype)
    node.set_start(gpsstart)
    node.set_end(gpsend)
    #self.add_node(node)
    return node
    
  def add_engine_node(self,end_time,extra_options=None):
    node=self.EngineNode(self.engine_job)
    for ifo in self.ifos:
      for seg in self.segments[ifo]:
	if end_time > seg.start and end_time < seg.end:
	  node.add_ifo_data(ifo,seg)
    if extra_options is not None:
      for opt in extra_options.keys():
	node.add_var_arg('--'+opt+' '+extra_options[opt])
    # Add the nodes it depends on
    for dfnode in node.__parents:
      if df_node not in self.__nodes:
	self.add_node(dfnode)
    self.add_node(node)
    

class EngineJob(pipeline.CondorDAGJob):
  def __init__(self,cp,submitFile,logdir):
    self.engine=cp.get('analysis','engine')
    exe=cp.get('condor',self.engine)
    pipeline.CondorDAGJob.__init__(self,"standard",exe)
    # Set the options which are always used
    self.set_sub_file(submitFile)
    self.add_ini_opts(cp,self.engine)
    self.set_stdout_file(os.path.join(logdir,'lalinference-$(cluster)-$(process).out'))
    self.set_stderr_file(os.path.join(logdir,'lalinference-$(cluster)-$(process).err'))

class LALInferenceNestNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencenest'

class LALInferenceMCMCNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencemcmc'
    
class EngineNode(pipeline.CondorDAGNode):
  def __init__(self,li_job):
    pipeline.CondorDAGNode.__init__(self,li_job)
  
  def set_seed(self,seed):
    self.add_var_opt('randomseed',seed)
  
  def set_dataseed(self,seed):
    self.add_var_opt('dataseed',seed)
    
  def get_ifos(self):
    return ''.join(map(str,self.__ifos))
      
  def set_trig_time(self,time):
    """
    Set the end time of the signal for the centre of the prior in time
    """
    self.__trigtime=float(time)
    self.add_var_opt('trigtime',str(time))
    
  def set_event_number(self,event):
    """
    Set the event number in the injection XML.
    """
    if event is not None:
      self.__event=int(event)
      self.add_var_opt('event',str(event))
    
  get_trig_time = lambda self: self.__trigtime
  
  def add_ifo_data(self,ifo,sciseg,timeslide=0):
    self.ifos.append(ifo)
    self.channels[ifo]=channel
    self.scisegs[ifo]=sciseg
    self.add_parent(sciseg.get_df_node())
    self.timeslides[ifo]=timeslide
  
  def finalize(self):
    self._finalize_ifo_data()
    pipeline.CondorDAGNode.finalize()
    
  def _finalize_ifo_data(self):
  """
  Add list of IFOs and data to analyse to command line arguments.
  """
  cp = self.job().get_cp()
  ifostring='['
  cachestring='['
  channelstring='['
  first=True
  for ifo in self.ifos:
	if first:
	  delim=''
	  first=False
	else: delim=','
	cache=self.scisegs[ifo].get_df_node().get_output_files()[0]
	self.add_parent(self.scisegs[ifo].get_df_node())
	ifostring=ifostring+delim+ifo
	cachestring=cachestring+delim+cache
	channelstring=channelstring+delim+self.job().get_cp().get('data',ifo.lower()+'-channel')
  ifostring=ifostring+']'
  cachestring=cachestring+']'
  channelstring=channelstring+']'
  self.add_var_arg('--IFO '+ifostring)
  self.add_var_arg('--channel '+channelstring)
  self.add_var_arg('--cache '+cachestring)
  # Start at earliest common time
  # NOTE: We perform this arithmetic for all ifos to ensure that a common data set is
  # Used when we are running the coherence test.
  # Otherwise the noise evidence will differ.
  starttime=max([int(self.scisegs[ifo].start()) for ifo in self.ifos])
  endtime=min([int(self.scisegs[ifo].end()) for ifo in self.ifos])
  self.__GPSstart=starttime
  self.__GPSend=endtime
  length=endtime-starttime
  
  # Now we need to adjust the start time and length to make sure the maximum data length
  # is not exceeded.
  trig_time=self.get_trig_time()
  maxLength=float(cp.get('analysis','analysis-chunk-length'))
  if(length > maxLength):
    while(self.__GPSstart+maxLength<trig_time and self.__GPSstart+maxLength<self.__GPSend):
	  self.__GPSstart+=maxLength/2.0
  # Override calculated start time if requested by user in ini file
  if self.job().get_cp().has_option(self.engine,'psdstart'):
    self.__GPSstart=self.job().get_cp().getfloat(self.engine,'psdstart')
    print 'Over-riding start time to user-specified value %f'%(self.__GPSstart)
    if self.__GPSstart<starttime or self.__GPSstart>endtime:
      print 'ERROR: Over-ridden time lies outside of science segment!'
      raise Exception('Bad psdstart specified')
  else: 
    self.add_var_opt('psdstart',str(self.__GPSstart))
  if self.job().get_cp().has_option(self.engine,'psdlength'):
    length=self.job().get_cp().getfloat(self.engine,'psdlength')
    print 'Over-riding PSD length to user-specified value %f'%(length)
  else:
    length=self.__GPSend-self.__GPSstart
    if(length>maxLength):
      length=maxLength
  self.add_var_opt('PSDlength',str(int(length)))
  self.add_var_opt('seglen',self.job().get_cp().get('analysis','psd-chunk-length'))

class ResultsPageJob(pipeline.CondorDAGJob):
  def __init__(self,cp,submitFile,logdir):
    exe=cp.get('condor','resultspage')
    pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
    self.set_sub_file(submitFile)
    self.set_stdout_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).out'))
    self.set_stderr_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).err'))
    self.add_condor_cmd('getenv','True')
    # self.add_opt('Nlive',cp.get('analysis','nlive'))
    
    if cp.has_option('results','skyres'):
    self.add_opt('skyres',cp.get('results','skyres'))
