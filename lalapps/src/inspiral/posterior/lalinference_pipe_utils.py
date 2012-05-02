# DAG Class definitions for LALInference Pipeline
# (C) 2012 John Veitch, Kiersten Ruisard, Kan Wang

import glue
from glue import pipeline,segmentsUtils
import os
from lalapps import inspiralutils
import uuid
import ast
import pdb

# We use the GLUE pipeline utilities to construct classes for each
# type of job. Each class has inputs and outputs, which are used to
# join together types of jobs into a DAG.

dummyCacheNames=['LALLIGO','LALVirgo','LALAdLIGO']

def mkdirs(path):
  """
  Helper function. Make the given directory, creating intermediate
  dirs if necessary, and don't complain about it already existing.
  """
  if os.access(path,os.W_OK) and os.path.isdir(path): return
  else: os.makedirs(path)

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
    mkdirs(self.basepath)
    daglogdir=cp.get('paths','daglogdir')
    mkdirs(daglogdir)
    self.daglogfile=os.path.join(daglogdir,'lalinference_pipeline-'+str(uuid.uuid1())+'.log')
    pipeline.CondorDAG.__init__(self,self.daglogfile,dax)
    if cp.has_option('paths','cachedir'):
      self.cachepath=cp.get('paths','cachedir')
    else:
      self.cachepath=os.path.join(self.basepath,'caches')
    mkdirs(self.cachepath)
    if cp.has_option('paths','logdir'):
      self.logpath=cp.get('paths','logdir')
    else:
      self.logpath=os.path.join(self.basepath,'log')
    mkdirs(self.logpath)
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
    self.channels=ast.literal_eval(cp.get('data','channels'))
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
	self.times.append(float(time))
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
    (mintime,maxtime)=self.get_required_data(self.times)
    if not self.config.has_option('input','gps-start-time'):
      self.config.set('input','gps-start-time',str(mintime))
    if not self.config.has_option('input','gps-end-time'):
      self.config.set('input','gps-end-time',str(maxtime))
    self.add_science_segments()
    
    # Save the final configuration that is being used
    conffilename=os.path.join(self.basepath,'config.ini')
    with open(conffilename,'wb') as conffile:
      self.config.write(conffile)
    
    # Generate the DAG according to the config given
    self.setup_from_times(self.times)
    self.dagfilename="lalinference_%s-%s"%(self.config.get('input','gps-start-time'),self.config.get('input','gps-end-time'))
    self.set_dag_file(os.path.join(self.basepath,self.dagfilename))
    
  def get_required_data(self,times):
    """
    Calculate teh data that will be needed to process all events
    """
    psdlength=self.config.getint('input','max-psd-length')
    padding=self.config.getint('input','padding')
    # Assume that the PSD is estimated from the interval (end_time+ buffer , psdlength)
    # Also require padding before start time
    return (min(times)-padding,max(times)+padding+psdlength)

  def setup_from_times(self,times):
    """
    Generate a DAG from a list of times
    """
    for time in self.times:
      self.add_full_analysis_time(time)
  
  def add_full_analysis_time(self,gpstime):
    """
    Analyse a given GPS time
    """
    # datafindnode=self.get_datafind_node(gpstime)
    enginenode=self.add_engine_node(gpstime)
    ifos=reduce(lambda a,b:a+b,enginenode.ifos)
    pagedir=os.path.join(self.basepath,str(gpstime)+'-'+str(enginenode.uuid),ifos)
    mkdirs(pagedir)
    self.add_results_page_node(outdir=pagedir,parent=enginenode)
  
  def add_science_segments(self):
    # Query the segment database for science segments and
    # add them to the pool of segments
    for ifo in self.ifos:
      (segFileName,dqVetoes)=inspiralutils.findSegmentsToAnalyze(self.config, ifo, self.veto_categories, generate_segments=True,\
	use_available_data=self.use_available_data , data_quality_vetoes=False)
      self.dqVetoes=dqVetoes
      segfile=open(segFileName)
      segs=segmentsUtils.fromsegwizard(segfile)
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
    node.set_trig_time(end_time)
    for ifo in self.ifos:
      for seg in self.segments[ifo]:
        print 'Looking at segment %s for ifo %s to see if it contains end time %f...'%(str(seg),str(ifo),end_time)
        if end_time >= seg.start() and end_time < seg.end():
          print '  Adding segment'
          node.add_ifo_data(ifo,seg,self.channels[ifo])
    if extra_options is not None:
      for opt in extra_options.keys():
	    node.add_var_arg('--'+opt+' '+extra_options[opt])
    # Add the nodes it depends on
    for seg in node.scisegs.values():
      dfnode=seg.get_df_node()
      if dfnode not in self.get_nodes():
    	self.add_node(dfnode)
    self.add_node(node)
    # Add control options
    node.set_seglen(self.config.getint('lalinference','seglen'))
    if self.config.has_option('input','psd-length'):
      node.set_psdlength(self.config.getint('input','psd-length'))
    if self.config.has_option('input','psd-start-time'):
      node.set_psdstart(self.config.getint('input','psd-start-time'))
    node.set_max_psdlength(self.config.getint('input','max-psd-length'))
    out_dir=os.path.join(self.basepath,'engine')
    mkdirs(out_dir)
    node.set_output_file(os.path.join(out_dir,node.engine+'-'+node.get_ifos()+'-'+str(node.get_trig_time())+'-'+'%x'%(id(node))))
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
    self.ifos=[]
    self.scisegs={}
    self.channels={}
    self.timeslides={}
    self.seglen=None
    self.psdlength=None
    self.maxlength=None
    self.psdstart=None
    self.uuid=uuid.uuid1()

  def set_seglen(self,seglen):
    self.seglen=seglen

  def set_psdlength(self,psdlength):
    self.psdlength=psdlength

  def set_max_psdlength(self,psdlength):
    self.maxlength=psdlength

  def set_psd_start(self,psdstart):
    self.psdstart=psdstart

  def set_seed(self,seed):
    self.add_var_opt('randomseed',seed)
  
  def set_dataseed(self,seed):
    self.add_var_opt('dataseed',seed)

  def get_ifos(self):
    return ''.join(map(str,self.ifos))
      
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
  
  def add_ifo_data(self,ifo,sciseg,channelname,timeslide=0):
    self.ifos.append(ifo)
    self.scisegs[ifo]=sciseg
    self.add_parent(sciseg.get_df_node())
    self.timeslides[ifo]=timeslide
    self.channels[ifo]=channelname
  
  def finalize(self):
    self._finalize_ifo_data()
    pipeline.CondorDAGNode.finalize(self)
    
  def _finalize_ifo_data(self):
      """
      Add list of IFOs and data to analyse to command line arguments.
      """
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
        self.add_input_file(cache)
        self.add_parent(self.scisegs[ifo].get_df_node())
        ifostring=ifostring+delim+ifo
        cachestring=cachestring+delim+cache
        channelstring=channelstring+delim+self.channels[ifo]
      ifostring=ifostring+']'
      cachestring=cachestring+']'
      channelstring=channelstring+']'
      self.add_var_opt('IFO',ifostring)
      self.add_var_opt('channel',channelstring)
      self.add_var_opt('cache',cachestring)
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
      maxLength=self.maxlength
      if(length > maxLength):
        while(self.__GPSstart+maxLength<trig_time and self.__GPSstart+maxLength<self.__GPSend):
          self.__GPSstart+=maxLength/2.0
      # Override calculated start time if requested by user in ini file
      if self.psdstart is not None:
        self.__GPSstart=self.psdstart
        print 'Over-riding start time to user-specified value %f'%(self.__GPSstart)
        if self.__GPSstart<starttime or self.__GPSstart>endtime:
          print 'ERROR: Over-ridden time lies outside of science segment!'
          raise Exception('Bad psdstart specified')
      else: 
        self.add_var_opt('psdstart',str(self.__GPSstart))
      if self.psdlength is None:
        self.psdlength=self.__GPSend-self.__GPSstart
        if(self.psdlength>self.maxlength):
          self.psdlength=self.maxlength
      self.add_var_opt('psdlength',self.psdlength)
      self.add_var_opt('seglen',self.seglen)

class LALInferenceNestNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencenest'
    self.outfilearg='outfile'
    
  def set_output_file(self,filename):
    self.add_file_opt(self.outfilearg,filename+'.dat',file_is_output_file=True)
    self.paramsfile=filename+'_params.txt'
    self.Bfilename=filename+'_B.txt'

class LALInferenceMCMCNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencemcmc'
    self.outfilearg='outfile'

  def set_output_file(self,filename):
    self.add_file_opt(self.outfilearg,filename)
    self.add_output_file(filename+'.00')
 
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
    def __init__(self,results_page_job,outpath=None):
        pipeline.CondorDAGNode.__init__(self,results_page_job)
        if outpath is not None:
            self.set_output_path(path)
    def set_output_path(self,path):
        self.webpath=path
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
        mkdirs(dir)
         
