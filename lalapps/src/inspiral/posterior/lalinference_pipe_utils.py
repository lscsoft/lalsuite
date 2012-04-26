# DAG Class definitions for LALInference Pipeline
# (C) 2012 John Veitch, Kiersten Ruisard, Kan Wang

import glue
from glue import pipeline,segmentsUtils
import os
from lalapps import inspiralutils
import uuid
import ast

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

def scan_timefile(self,timefile):
    import re
    p=re.compile('[\d.]+')
    times=[]
    timefilehandle=open(timefile,'r')
    for time in timefilehandle:
      if not p.match(time):
	continue
      if float(time) in times:
	print 'Skipping duplicate time %s'%(time)
	continue
      print 'Read time %s'%(time)
      times.append(float(time))
      timefilehandle.close()
      return times
  
class LALInferencePipelineDAG(pipeline.CondorDAG):
  def __init__(self,cp,dax=False):
    self.subfiles=[]
    self.config=cp
    self.engine=cp.get('analysis','engine')
    self.EngineNode=chooseEngineNode(self.engine)
    if cp.has_option('paths','basedir'):
      self.basepath=cp.get('paths','basedir')
    else:
      self.basepath=os.getcwd()
      print 'No basepath specified, using current directory: %s'%(self.basepath)
    daglogdir=cp.get('paths','daglogdir')
    self.daglogfile=os.path.join(daglogdir,'lalinference_pipeline-'+str(uuid.uuid1())+'.log')
    pipeline.CondorDAG.__init__(self,self.daglogfile,dax)
    if cp.has_option('paths','cachedir'):
      self.cachepath=cp.get('paths','cachedir')
    else:
      self.cachepath=os.path.join(self.basepath,'caches')
    if cp.has_option('paths','logdir'):
      self.logpath=cp.get('paths','logdir')
    else:
      self.logpath=os.path.join(self.basepath,'log')
    if cp.has_option('analysis','ifos'):
      self.ifos=ast.literal_eval(cp.get('analysis','ifos'))
    else:
      self.ifos=['H1','L1','V1']
    self.segments={}
    if cp.has_option('datafind','veto-categories'):
      self.veto_categories=cp.get('datafind','veto-categories')
    else: self.veto_categories=[]
    for ifo in self.ifos:
      self.segments[ifo]=[]
    self.dq={}
    self.frtypes=ast.literal_eval(cp.get('datafind','types'))
    self.use_available_data=False
    self.webdir=cp.get('paths','webdir')
    # Set up necessary job files.
    self.datafind_job = pipeline.LSCDataFindJob(self.cachepath,self.logpath,self.config)
    self.datafind_job.add_opt('url-type','file')
    self.datafind_job.set_sub_file(os.path.join(self.basepath,'datafind.sub'))
    self.engine_job = EngineJob(self.config, os.path.join(self.basepath,'lalinference.sub'),self.logpath)
    self.results_page_job = ResultsPageJob(self.config,os.path.join(self.basepath,'resultspage.sub'),self.logpath)


    # Process the input to build list of analyses to do
    self.times=[]
    if cp.has_option('input','gps-time-file'):
      times=scan_timefile(cp.get('input','gps-time-file'))
      for time in times:
	self.times.append(time)
    # SimInspiral Table
    if cp.has_option('input','injection-file'):
      from pylal import SimInspiralUtils
      injTable=SimInspiralUtils.ReadSimInspiralFromFiles([cp.get('input','injection-file')])
      map(self.times.append, [inj.get_end() for inj in injTable])
    # SnglInspiral Table
    if cp.has_option('input','sngl-inspiral-file'):
      from pylal import SnglInspiralUtils
      trigTable=SnglInspiralUtils.ReadSnglInspiralFromFiles([opts.single_triggers])
      map(self.times.append,[trig.get_end() for trig in trigTable])

    # CoincInspiralTable
    
    # Pipedown database

    # Sanity checking
    if len(self.times)==0:
      print 'No input times found, please check your config. Generating an empty DAG'
    
    # Set up the segments
    if not self.config.has_option('input','gps-start-time'):
      self.config.set('input','gps-start-time',str(min(self.times)))
    if not self.config.has_option('input','gps-end-time'):
      self.config.set('input','gps-end-time',str(max(self.times)))
    self.add_science_segments()
    
    # Save the final configuration that is being used
    conffilename=os.path.join(self.basepath,'config.ini')
    with open(conffilename,'wb') as conffile:
      self.cp.write(conffile)
    
    # Generate the DAG according to the config given
    
  
  def setup_from_times(self,times):
    """
    Generate a DAG from a list of times
    """
    for time in self.times:
      self.add_full_analysis_time(str(time))
  
  def add_full_analysis_time(self,gpstime):
    """
    Analyse a given GPS time
    """
    datafindnode=self.add_datafind(gpstime)
    enginenode=self.add_engine(gpstime)
    ifos=reduce(lambda a,b:a+b,enginenode.ifos)
    pagedir=os.path.join(ifos,str(gpstime)+'-'+str(id(enginenode)))
    self.add_results_page_node(outdir=pagedir,parent=enginenode)
  
  def get_science_segment(self,ifo,gpstime):
    # Check if time is in existing segment
    for seg in self.segments[ifo]:
      if gpstime in seg: return seg
    raise pipeline.CondorDAGError('Unable to find time in segments')
  
  def add_science_segments(self):
    # Query the segment database for science segments and
    # add them to the pool of segments
    for ifo in self.ifos:
      (segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(self.config, ifo, self.veto_categories, generate_segments=True,\
	use_available_data=self.use_available_data , data_quality_vetoes=False)
      self.dqVetoes=dqVetoes
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
    return node
    
  def add_results_page_node(self,outdir=None,parent=None,extra_options=None):
    node=ResultsPageNode(self.results_page_job)
    if parent is not None:
      node.add_parent(parent)
      infiles=parent.get_output_files()
      for infile in infiles:
	node.add_var_arg(infile)
    node.set_output_dir(outdir)
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

class LALInferenceNestNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencenest'
    
  def set_output_file(self,filename):
    self.add_file_opt(self.outfilearg,filename,file_is_output_file=True)
    self.paramsfile=filename+'_params.txt'
    self.Bfilename=filename+'_B.txt'

class LALInferenceMCMCNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencemcmc'
 
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

class ResultsPageNode(pipeline.CondorDAGNode):
    def __init__(self,results_page_job):
        pipeline.CondorDAGNode.__init__(self,results_page_job)
        self.webpath=self.job().get_cp().get('paths','webdir')
    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_arg('--eventnum '+str(event))
    def add_engine_parent(self,node):
      """
      Add a parent node which is one of the engine nodes
      And automatically set options accordingly
      """
      self.add_parent(node)
      for infile in node.get_output_files():
	    self.add_file_arg(infile)
      if isinstance(node, LALInferenceNestNode):
	    self.add_var_opt('ns','')
	
      if isinstance(node,LALInferenceMCMCNode):
	    self.add_var_opt('lalinfmcmc','')
    def set_output_dir(self,dir):
        self.add_var_opt('outpath',dir)
        inspiralutils.mkdir(dir)
         
