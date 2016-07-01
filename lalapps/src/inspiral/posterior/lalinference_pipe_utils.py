
#flow DAG Class definitions for LALInference Pipeline
# (C) 2012 John Veitch, Vivien Raymond, Kiersten Ruisard, Kan Wang

import itertools
import glue
from glue import pipeline,segmentsUtils,segments
import os
import socket
from lalapps import inspiralutils
import uuid
import ast
import pdb
import string
from math import floor,ceil,log,pow
import sys
import random
from itertools import permutations
import shutil
import numpy as np

# We use the GLUE pipeline utilities to construct classes for each
# type of job. Each class has inputs and outputs, which are used to
# join together types of jobs into a DAG.

class Event():
  """
  Represents a unique event to run on
  """
  new_id=itertools.count().next
  def __init__(self,trig_time=None,SimInspiral=None,SimBurst=None,SnglInspiral=None,CoincInspiral=None,event_id=None,timeslide_dict=None,GID=None,ifos=None, duration=None,srate=None,trigSNR=None,fhigh=None):
    self.trig_time=trig_time
    self.injection=SimInspiral
    self.burstinjection=SimBurst
    self.sngltrigger=SnglInspiral
    if timeslide_dict is None:
      self.timeslides={}
    else:
      self.timeslides=timeslide_dict
    self.GID=GID
    self.coinctrigger=CoincInspiral
    if ifos is None:
      self.ifos = []
    else:
      self.ifos = ifos
    self.duration = duration
    self.srate = srate
    self.trigSNR = trigSNR
    self.fhigh = fhigh
    if event_id is not None:
        self.event_id=event_id
    else:
        self.event_id=Event.new_id()
    if self.injection is not None:
        self.trig_time=self.injection.get_end()
        if event_id is None: self.event_id=int(str(self.injection.simulation_id).split(':')[2])
    if self.burstinjection is not None:
        self.trig_time=self.burstinjection.get_end()
        if event_id is None: self.event_id=int(str(self.burstinjection.simulation_id).split(':')[2])
    if self.sngltrigger is not None:
        self.trig_time=self.sngltrigger.get_end()
        self.event_id=int(str(self.sngltrigger.event_id).split(':')[2])
    if self.coinctrigger is not None:
        self.trig_time=self.coinctrigger.end_time + 1.0e-9 * self.coinctrigger.end_time_ns
    if self.GID is not None:
        self.event_id=int(''.join(i for i in self.GID if i.isdigit()))
    self.engine_opts={}
  def set_engine_option(self,opt,val):
    """
    Can set event-specific options for the engine nodes
    using this option, e.g. ev.set_engine_option('time-min','1083759273')
    """
    self.engine_opts[opt]=val

dummyCacheNames=['LALLIGO','LALVirgo','LALAdLIGO','LALAdVirgo']

def readLValert(SNRthreshold=0,gid=None,flow=40.0,gracedb="gracedb",basepath="./",downloadpsd=True):
  """
  Parse LV alert file, containing coinc, sngl, coinc_event_map.
  and create a list of Events as input for pipeline
  Based on Chris Pankow's script
  """
  output=[]
  from glue.ligolw import utils
  from glue.ligolw import lsctables
  from glue.ligolw import ligolw
  class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
    pass
  lsctables.use_in(LIGOLWContentHandler)
  from glue.ligolw import param
  from glue.ligolw import array
  from pylal import series as lalseries
  import subprocess
  from subprocess import Popen, PIPE
  cwd=os.getcwd()
  os.chdir(basepath)
  print "%s download %s coinc.xml"%(gracedb,gid)
  subprocess.call([gracedb,"download", gid ,"coinc.xml"])
  xmldoc=utils.load_filename("coinc.xml",contenthandler = LIGOLWContentHandler)
  coinctable = lsctables.CoincInspiralTable.get_table(xmldoc)
  coinc_events = [event for event in coinctable]
  sngltable = lsctables.SnglInspiralTable.get_table(xmldoc)
  sngl_events = [event for event in sngltable]
  #Issues to identify IFO with good data that did not produce a trigger
  #search_summary = lsctables.getTablesByType(xmldoc, lsctables.SearchSummaryTable)[0]
  #ifos = search_summary[0].ifos.split(",")
  #coinc_table = lsctables.getTablesByType(xmldoc, lsctables.CoincTable)[0]
  #ifos = coinc_table[0].instruments.split(",")
  trigSNR = 2.0*coinctable[0].snr #The factor of 2.0 is because detection pipelines recover SNR lower than PE can recover.
  # Parse PSD
  srate_psdfile=16384
  ifos=None
  if downloadpsd:
    print "gracedb download %s psd.xml.gz" % gid
    subprocess.call([gracedb,"download", gid ,"psd.xml.gz"])
    xmlpsd = lalseries.read_psd_xmldoc(utils.load_filename('psd.xml.gz',contenthandler = lalseries.LIGOLWContentHandler))
    # Note: This finds the active IFOs by looking for available PSDs
    # Is there another way of getting this info?
    ifos = xmlpsd.keys()
  psdasciidic=None
  fhigh=None
  if os.path.exists("psd.xml.gz"):
    psdasciidic=get_xml_psds(os.path.realpath("./psd.xml.gz"),ifos,os.path.realpath('./PSDs'),end_time=None)
    combine=np.loadtxt(psdasciidic[psdasciidic.keys()[0]])
    srate_psdfile = pow(2.0, ceil( log(float(combine[-1][0]), 2) ) ) * 2
  else:
    print "Failed to gracedb download %s psd.xml.gz. lalinference will estimate the psd itself." % gid
  # Logic for template duration and sample rate disabled
  coinc_map = lsctables.CoincMapTable.get_table(xmldoc)
  for coinc in coinc_events:
    these_sngls = [e for e in sngl_events if e.event_id in [c.event_id for c in coinc_map if c.coinc_event_id == coinc.coinc_event_id] ]
    dur=[]
    srate=[]
    for e in these_sngls:
      # Review: Replace this with a call to LALSimulation function at some point
      p=Popen(["lalapps_chirplen","--flow",str(flow),"-m1",str(e.mass1),"-m2",str(e.mass2)],stdout=PIPE, stderr=PIPE, stdin=PIPE)
      strlen = p.stdout.read()
      dur.append(pow(2.0, ceil( log(max(8.0,float(strlen.splitlines()[2].split()[5]) + 2.0), 2) ) ) )
      srate.append(pow(2.0, ceil( log(float(strlen.splitlines()[1].split()[5]), 2) ) ) * 2 )
    if max(srate)<srate_psdfile:
      srate = max(srate)
    else:
      srate = srate_psdfile
      fhigh = srate_psdfile/2.0 * 0.95 # Because of the drop-off near Nyquist of the PSD from gstlal
    ev=Event(CoincInspiral=coinc, GID=gid, ifos = ifos, duration = max(dur), srate = srate, trigSNR = trigSNR, fhigh = fhigh)
    if(coinc.snr>SNRthreshold): output.append(ev)

  print "Found %d coinc events in table." % len(coinc_events)
  os.chdir(cwd)
  return output

def open_pipedown_database(database_filename,tmp_space):
    """
    Open the connection to the pipedown database
    """
    if not os.access(database_filename,os.R_OK):
	raise Exception('Unable to open input file: %s'%(database_filename))
    from glue.ligolw import dbtables
    import sqlite3
    working_filename=dbtables.get_connection_filename(database_filename,tmp_path=tmp_space)
    connection = sqlite3.connect(working_filename)
    if tmp_space:
	dbtables.set_temp_store_directory(connection,tmp_space)
    #dbtables.DBTable_set_connection(connection)
    return (connection,working_filename)

def get_zerolag_lloid(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1, min_cfar=-1):
	"""
	Returns a list of Event objects
	from pipedown data base. Can dump some stats to dumpfile if given,
	and filter by gpsstart and gpsend to reduce the nunmber or specify
	max_cfar to select by combined FAR
	"""
	output={}
	if gpsstart is not None: gpsstart=float(gpsstart)
	if gpsend is not None: gpsend=float(gpsend)
	# Get coincs
	get_coincs = "SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
		FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name=='sngl_inspiral' and coinc_event_map.event_id ==\
		sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) \
		join coinc_inspiral on (coinc_event.coinc_event_id==coinc_inspiral.coinc_event_id) \
        WHERE coinc_event.time_slide_id=='time_slide:time_slide_id:1'\
		"
	if gpsstart is not None:
		get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 > %f'%(gpsstart)
	if gpsend is not None:
		get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 < %f'%(gpsend)
	if max_cfar !=-1:
		get_coincs=get_coincs+' and coinc_inspiral.combined_far < %f'%(max_cfar)
	if min_cfar != -1:
		get_coincs=get_coincs+' and coinc_inspiral.combined_far > %f'%(min_cfar)
	db_out=database_connection.cursor().execute(get_coincs)
    	extra={}
	for (sngl_time, ifo, coinc_id, snr, chisq, cfar) in db_out:
      		coinc_id=int(coinc_id.split(":")[-1])
	  	if not coinc_id in output.keys():
			output[coinc_id]=Event(trig_time=sngl_time,timeslide_dict={},event_id=int(coinc_id))
			extra[coinc_id]={}
		output[coinc_id].timeslides[ifo]=0
		output[coinc_id].ifos.append(ifo)
		extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
	if dumpfile is not None:
		fh=open(dumpfile,'w')
		for co in output.keys():
			for ifo in output[co].ifos:
				fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
		fh.close()
	return output.values()


def get_zerolag_pipedown(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1, min_cfar=-1):
	"""
	Returns a list of Event objects
	from pipedown data base. Can dump some stats to dumpfile if given,
	and filter by gpsstart and gpsend to reduce the nunmber or specify
	max_cfar to select by combined FAR
	"""
	output={}
	if gpsstart is not None: gpsstart=float(gpsstart)
	if gpsend is not None: gpsend=float(gpsend)
	# Get coincs
	get_coincs = "SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
		FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name=='sngl_inspiral' and coinc_event_map.event_id ==\
		sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) \
		join coinc_inspiral on (coinc_event.coinc_event_id==coinc_inspiral.coinc_event_id) \
		WHERE coinc_event.time_slide_id=='time_slide:time_slide_id:10049'\
		"
	if gpsstart is not None:
		get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 > %f'%(gpsstart)
	if gpsend is not None:
		get_coincs=get_coincs+' and coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1.0e-9 < %f'%(gpsend)
	if max_cfar !=-1:
		get_coincs=get_coincs+' and coinc_inspiral.combined_far < %f'%(max_cfar)
	if min_cfar != -1:
		get_coincs=get_coincs+' and coinc_inspiral.combined_far > %f'%(min_cfar)
	db_out=database_connection.cursor().execute(get_coincs)
    	extra={}
	for (sngl_time, ifo, coinc_id, snr, chisq, cfar) in db_out:
      		coinc_id=int(coinc_id.split(":")[-1])
	  	if not coinc_id in output.keys():
			output[coinc_id]=Event(trig_time=sngl_time,timeslide_dict={},event_id=int(coinc_id))
			extra[coinc_id]={}
		output[coinc_id].timeslides[ifo]=0
		output[coinc_id].ifos.append(ifo)
		extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
	if dumpfile is not None:
		fh=open(dumpfile,'w')
		for co in output.keys():
			for ifo in output[co].ifos:
				fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
		fh.close()
	return output.values()


def get_timeslides_pipedown(database_connection, dumpfile=None, gpsstart=None, gpsend=None, max_cfar=-1):
	"""
	Returns a list of Event objects
	with times and timeslide offsets
	"""
	output={}
	if gpsstart is not None: gpsstart=float(gpsstart)
	if gpsend is not None: gpsend=float(gpsend)
	db_segments=[]
	sql_seg_query="SELECT search_summary.out_start_time, search_summary.out_end_time from search_summary join process on process.process_id==search_summary.process_id where process.program=='thinca'"
	db_out = database_connection.cursor().execute(sql_seg_query)
	for d in db_out:
		if d not in db_segments:
			db_segments.append(d)
	seglist=segments.segmentlist([segments.segment(d[0],d[1]) for d in db_segments])
	db_out_saved=[]
	# Get coincidences
	get_coincs="SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,time_slide.offset,sngl_inspiral.ifo,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq,coinc_inspiral.combined_far \
		    FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name == 'sngl_inspiral' and coinc_event_map.event_id \
		    == sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id) join time_slide\
		    on (time_slide.time_slide_id == coinc_event.time_slide_id and time_slide.instrument==sngl_inspiral.ifo)\
		    join coinc_inspiral on (coinc_inspiral.coinc_event_id==coinc_event.coinc_event_id) where coinc_event.time_slide_id!='time_slide:time_slide_id:10049'"
	joinstr = ' and '
	if gpsstart is not None:
		get_coincs=get_coincs+ joinstr + ' coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1e-9 > %f'%(gpsstart)
	if gpsend is not None:
		get_coincs=get_coincs+ joinstr+' coinc_inspiral.end_time+coinc_inspiral.end_time_ns*1e-9 <%f'%(gpsend)
	if max_cfar!=-1:
		get_coincs=get_coincs+joinstr+' coinc_inspiral.combined_far < %f'%(max_cfar)
	db_out=database_connection.cursor().execute(get_coincs)
	from pylal import SnglInspiralUtils
	extra={}
	for (sngl_time, slide, ifo, coinc_id, snr, chisq, cfar) in db_out:
		coinc_id=int(coinc_id.split(":")[-1])
		seg=filter(lambda seg:sngl_time in seg,seglist)[0]
		slid_time = SnglInspiralUtils.slideTimeOnRing(sngl_time,slide,seg)
		if not coinc_id in output.keys():
			output[coinc_id]=Event(trig_time=slid_time,timeslide_dict={},event_id=int(coinc_id))
			extra[coinc_id]={}
		output[coinc_id].timeslides[ifo]=slid_time-sngl_time
		output[coinc_id].ifos.append(ifo)
		extra[coinc_id][ifo]={'snr':snr,'chisq':chisq,'cfar':cfar}
	if dumpfile is not None:
		fh=open(dumpfile,'w')
		for co in output.keys():
			for ifo in output[co].ifos:
				fh.write('%s %s %s %s %s %s %s\n'%(str(co),ifo,str(output[co].trig_time),str(output[co].timeslides[ifo]),str(extra[co][ifo]['snr']),str(extra[co][ifo]['chisq']),str(extra[co][ifo]['cfar'])))
		fh.close()
	return output.values()

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
  if name=='lalinferenceburst':
    return LALInferenceBurstNode
  if name=='lalinferencemcmc':
    return LALInferenceMCMCNode
  if name=='lalinferencebambi' or name=='lalinferencebambimpi':
    return LALInferenceBAMBINode
  if name=='lalinferencedatadump':
    return LALInferenceDataDumpNode
  return EngineNode

def get_engine_name(cp):
    name=cp.get('analysis','engine')
    if name=='random':
        engine_list=['lalinferencenest','lalinferencemcmc','lalinferencebambimpi']
        if cp.has_option('input','gid'):
            gid=cp.get('input','gid')
            engine_number=int(''.join(i for i in gid if i.isdigit())) % 2
        else:
            engine_number=random.randint(0,1)
        return engine_list[engine_number]
    else:
        return name


def scan_timefile(timefile):
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

def get_xml_psds(psdxml,ifos,outpath,end_time=None):
  """
  Get a psd.xml.gz file and:
  1) Reads it
  2) Converts PSD (10e-44) -> ASD ( ~10e-22)
  3) Checks the psd file contains all the IFO we want to analyze
  4) Writes down the ASDs into an ascii file for each IFO in psd.xml.gz. The name of the file contains the trigtime (if given) and the ifo name.
  Input:
    psdxml: psd.xml.gz file
    ifos: list of ifos used for the analysis
    outpath: path where the ascii ASD will be written to
    (end_time): trigtime for this event. Will be used a part of the ASD file name
  """
  lal=1
  from glue.ligolw import utils
  try:
    from pylal import series
    lal=0
  except ImportError:
    print "ERROR, cannot import pylal.series in bppu/get_xml_psds()\n"
    exit(1)

  out={}
  if not os.path.isdir(outpath):
    os.makedirs(outpath)
  if end_time is not None:
    time=repr(float(end_time))
  else:
    time=''
  #check we don't already have ALL the psd files #
  got_all=1
  for ifo in ifos:
    path_to_ascii_psd=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
    # Check we don't already have that ascii (e.g. because we are running parallel runs of the save event
    if os.path.isfile(path_to_ascii_psd):
      got_all*=1
    else:
      got_all*=0
  if got_all==1:
    #print "Already have PSD files. Nothing to do...\n"
    for ifo in ifos:
      out[ifo]=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
    return out

  # We need to convert the PSD for one or more IFOS. Open the file
  if not os.path.isfile(psdxml):
    print "ERROR: impossible to open the psd file %s. Exiting...\n"%psdxml
    sys.exit(1)
  xmlpsd =  series.read_psd_xmldoc(utils.load_filename(psdxml,contenthandler = series.LIGOLWContentHandler))
  # Check the psd file contains all the IFOs we want to analize
  for ifo in ifos:
    if not ifo in [i.encode('ascii') for i in xmlpsd.keys()]:
      print "ERROR. The PSD for the ifo %s does not seem to be contained in %s\n"%(ifo,psdxml)
      sys.exit(1)
  #loop over ifos in psd xml file
  for instrument in xmlpsd.keys():
    #name of the ascii file we are going to write the PSD into
    path_to_ascii_psd=os.path.join(outpath,instrument.encode('ascii')+'_psd_'+time+'.txt')
    # Check we don't already have that ascii (e.g. because we are running parallel runs of the save event
    if os.path.isfile(path_to_ascii_psd):
      continue
    # get data for the IFO
    ifodata=xmlpsd[instrument]
    #check data is not empty
    if ifodata is None:
      continue
    # we have data. Get psd array
    if lal==0:
      #pylal stores the series in ifodata.data
      data=ifodata
    else:
      # lal stores it in ifodata.data.data
      data=ifodata.data
    # Fill a two columns array of (freq, psd) and save it in the ascii file
    f0=ifodata.f0
    deltaF=ifodata.deltaF

    combine=[]

    for i in np.arange(len(data.data)) :
      combine.append([f0+i*deltaF,np.sqrt(data.data[i])])
    np.savetxt(path_to_ascii_psd,combine)
    ifo=instrument.encode('ascii')
    # set node.psds dictionary with the path to the ascii files
    out[ifo]=os.path.join(outpath,ifo+'_psd_'+time+'.txt')
  return out

def create_pfn_tuple(filename,protocol='file://',site='local'):
    return( (os.path.basename(filename),protocol+os.path.abspath(filename),site) )

class LALInferencePipelineDAG(pipeline.CondorDAG):
  def __init__(self,cp,dax=False,first_dag=True,previous_dag=None,site='local'):
    self.subfiles=[]
    self.config=cp
    self.engine=get_engine_name(cp)
    self.EngineNode=chooseEngineNode(self.engine)
    self.site=site
    if cp.has_option('paths','basedir'):
      self.basepath=cp.get('paths','basedir')
    else:
      self.basepath=os.getcwd()
      print 'No basepath specified, using current directory: %s'%(self.basepath)
    mkdirs(self.basepath)
    if dax:
        os.chdir(self.basepath)
    self.posteriorpath=os.path.join(self.basepath,'posterior_samples')
    mkdirs(self.posteriorpath)
    if first_dag:
      daglogdir=cp.get('paths','daglogdir')
      mkdirs(daglogdir)
      self.daglogfile=os.path.join(daglogdir,'lalinference_pipeline-'+str(uuid.uuid1())+'.log')
      pipeline.CondorDAG.__init__(self,self.daglogfile,dax=dax)
    elif not first_dag and previous_dag is not None:
      daglogdir=cp.get('paths','daglogdir')
      mkdirs(daglogdir)
      self.daglogfile=os.path.join(daglogdir,'lalinference_pipeline-'+str(uuid.uuid1())+'.log')
      pipeline.CondorDAG.__init__(self,self.daglogfile,dax=dax)
      for node in previous_dag.get_nodes():
        self.add_node(node)
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
    self.computeroqweightsnodes={}
    self.bayeslinenodes={}
    self.dq={}
    self.frtypes=ast.literal_eval(cp.get('datafind','types'))
    self.channels=ast.literal_eval(cp.get('data','channels'))
    self.use_available_data=False
    self.webdir=cp.get('paths','webdir')
    if cp.has_option('analysis','dataseed'):
      self.dataseed=cp.getint('analysis','dataseed')
    else:
      self.dataseed=None
    # Set up necessary job files.
    self.prenodes={}
    self.datafind_job = pipeline.LSCDataFindJob(self.cachepath,self.logpath,self.config,dax=self.is_dax())
    self.datafind_job.add_opt('url-type','file')
    if cp.has_option('analysis','accounting_group'):
      self.datafind_job.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    if cp.has_option('analysis','accounting_group_user'):
      self.datafind_job.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
    self.datafind_job.set_sub_file(os.path.abspath(os.path.join(self.basepath,'datafind.sub')))
    self.preengine_job = EngineJob(self.config, os.path.join(self.basepath,'prelalinference.sub'),self.logpath,engine='lalinferencedatadump',ispreengine=True,dax=self.is_dax())
    self.preengine_job.set_grid_site('local')
    self.preengine_job.set_universe('vanilla')
    if self.config.has_option('condor','computeroqweights'):
      self.computeroqweights_job = ROMJob(self.config,os.path.join(self.basepath,'computeroqweights.sub'),self.logpath,dax=self.is_dax())
      self.computeroqweights_job.set_grid_site('local')
    if self.config.has_option('condor','bayesline'):
      self.bayesline_job = BayesLineJob(self.config,os.path.join(self.basepath,'bayesline.sub'),self.logpath,dax=self.is_dax())
      self.bayesline_job.set_grid_site('local')
    # Need to create a job file for each IFO combination
    self.engine_jobs={}
    ifocombos=[]
    for N in range(1,len(self.ifos)+1):
        for a in permutations(self.ifos,N):
            ifocombos.append(a)
    for ifos in ifocombos:
        self.engine_jobs[ifos] = EngineJob(self.config, os.path.join(self.basepath,'engine_%s.sub'%(reduce(lambda x,y:x+y, map(str,ifos)))),self.logpath,engine=self.engine,dax=self.is_dax(), site=site)
    self.results_page_job = ResultsPageJob(self.config,os.path.join(self.basepath,'resultspage.sub'),self.logpath,dax=self.is_dax())
    self.results_page_job.set_grid_site('local')
    self.cotest_results_page_job = ResultsPageJob(self.config,os.path.join(self.basepath,'resultspagecoherent.sub'),self.logpath,dax=self.is_dax())
    self.cotest_results_page_job.set_grid_site('local')
    self.merge_job = MergeNSJob(self.config,os.path.join(self.basepath,'merge_runs.sub'),self.logpath,dax=self.is_dax())
    self.merge_job.set_grid_site('local')
    self.coherence_test_job = CoherenceTestJob(self.config,os.path.join(self.basepath,'coherence_test.sub'),self.logpath,dax=self.is_dax())
    self.coherence_test_job.set_grid_site('local')
    self.gracedbjob = GraceDBJob(self.config,os.path.join(self.basepath,'gracedb.sub'),self.logpath,dax=self.is_dax())
    self.gracedbjob.set_grid_site('local')
    # Process the input to build list of analyses to do
    self.events=self.setup_from_inputs()

    # Sanity checking
    if len(self.events)==0:
      print 'No input events found, please check your config if you expect some events'
    self.times=[e.trig_time for e in self.events]

    # Set up the segments
    if not (self.config.has_option('input','gps-start-time') and self.config.has_option('input','gps-end-time')) and len(self.times)>0:
      (mintime,maxtime)=self.get_required_data(self.times)
      if not self.config.has_option('input','gps-start-time'):
        self.config.set('input','gps-start-time',str(int(floor(mintime))))
      if not self.config.has_option('input','gps-end-time'):
        self.config.set('input','gps-end-time',str(int(ceil(maxtime))))
    self.add_science_segments()

    # Save the final configuration that is being used
    # first to the run dir
    conffilename=os.path.join(self.basepath,'config.ini')
    with open(conffilename,'wb') as conffile:
      self.config.write(conffile)
    if self.config.has_option('paths','webdir'):
      mkdirs(self.config.get('paths','webdir'))
      with open(os.path.join(self.config.get('paths','webdir'),'config.ini'),'wb') as conffile:
        self.config.write(conffile)

    # Generate the DAG according to the config given
    for event in self.events: self.add_full_analysis(event)
    self.add_skyarea_followup()
    if self.config.has_option('analysis','upload-to-gracedb'):
      if self.config.getboolean('analysis','upload-to-gracedb'):
        self.add_gracedb_FITSskymap_upload(self.events[0],engine=self.engine)
    self.dagfilename="lalinference_%s-%s"%(self.config.get('input','gps-start-time'),self.config.get('input','gps-end-time'))
    self.set_dag_file(self.dagfilename)
    if self.is_dax():
      self.set_dax_file(self.dagfilename)

  def add_skyarea_followup(self):
    # Add skyarea jobs if the executable is given
    # Do one for each results page for now
    if self.config.has_option('condor','skyarea'):
      self.skyareajob=SkyAreaJob(self.config,os.path.join(self.basepath,'skyarea.sub'),self.logpath,dax=self.is_dax())
      respagenodes=filter(lambda x: isinstance(x,ResultsPageNode) ,self.get_nodes())
      if self.engine=='lalinferenceburst':
          prefix='LIB_'
      else:
          prefix='LALInference_'
      for p in respagenodes:
          skyareanode=SkyAreaNode(self.skyareajob,prefix=prefix)
          skyareanode.add_resultspage_parent(p)
          skyareanode.set_ifos(p.ifos)
          self.add_node(skyareanode)

  def add_full_analysis(self,event):
    if self.engine=='lalinferencenest' or  self.engine=='lalinferenceburst':
      result=self.add_full_analysis_lalinferencenest(event)
    elif self.engine=='lalinferencemcmc':
      result=self.add_full_analysis_lalinferencemcmc(event)
    elif self.engine=='lalinferencebambi' or self.engine=='lalinferencebambimpi':
      result=self.add_full_analysis_lalinferencebambi(event)
    return result


  def create_frame_pfn_file(self):
    """
    Create a pegasus cache file name, uses inspiralutils
    """
    import inspiralutils as iu
    gpsstart=self.config.get('input','gps-start-time')
    gpsend=self.config.get('input','gps-end-time')
    pfnfile=iu.create_frame_pfn_file(self.frtypes,gpsstart,gpsend)
    return pfnfile

  def get_required_data(self,times):
    """
    Calculate the data that will be needed to process all events
    """
    #psdlength = self.config.getint('input','max-psd-length')
    padding=self.config.getint('input','padding')
    if self.config.has_option('engine','seglen') or self.config.has_option('lalinference','seglen'):
      if self.config.has_option('engine','seglen'):
        seglen = self.config.getint('engine','seglen')
      if self.config.has_option('lalinference','seglen'):
        seglen = self.config.getint('lalinference','seglen')

      if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')):
        psdlength = 0
      else:
        psdlength = 32*seglen
    else:
      seglen = max(e.duration for e in self.events)
      if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')):
        psdlength = 0
      else:
        psdlength = 32*seglen
    # Assume that the data interval is (end_time - seglen -padding , end_time + psdlength +padding )
    # -> change to (trig_time - seglen - padding - psdlength + 2 , trig_time + padding + 2) to estimate the psd before the trigger for online follow-up.
    # Also require padding before start time
    return (min(times)-padding-seglen-psdlength+2,max(times)+padding+2)

  def setup_from_times(self,times):
    """
    Generate a DAG from a list of times
    """
    for time in self.times:
      self.add_full_analysis(Event(trig_time=time))

  def select_events(self):
    """
    Read events from the config parser. Understands both ranges and comma separated events, or combinations
    eg. events=[0,1,5:10,21] adds to the analysis the events: 0,1,5,6,7,8,9,10 and 21
    """
    events=[]
    times=[]
    raw_events=self.config.get('input','events').replace('[','').replace(']','').split(',')
    for raw_event in raw_events:
        if ':' in raw_event:
            limits=raw_event.split(':')
            if len(limits) != 2:
                print "Error: in event config option; ':' must separate two numbers."
                exit(0)
            low=int(limits[0])
            high=int(limits[1])
            if low>high:
                events.extend(range(int(high),int(low)+1))
            elif high>low:
                events.extend(range(int(low),int(high)+1))
        else:
            events.append(int(raw_event))
    return events

  def setup_from_inputs(self):
    """
    Scan the list of inputs, i.e.
    gps-time-file, injection-file, sngl-inspiral-file, coinc-inspiral-file, pipedown-database
    in the [input] section of the ini file.
    And process the events found therein
    """
    events=[]
    gpsstart=None
    gpsend=None
    if self.config.has_option('input','gps-start-time'):
      gpsstart=self.config.getfloat('input','gps-start-time')
    if self.config.has_option('input','gps-end-time'):
      gpsend=self.config.getfloat('input','gps-end-time')
    inputnames=['gps-time-file','burst-injection-file','injection-file','sngl-inspiral-file','coinc-inspiral-file','pipedown-db','gid','gstlal-db']
    ReadInputFromList=sum([ 1 if self.config.has_option('input',name) else 0 for name in inputnames])
    # If no input events given, just return an empty list (e.g. for PP pipeline)
    if ReadInputFromList!=1 and (gpsstart is None or gpsend is None):
        return []
    # Review: Clean up this section
    if self.config.has_option('input','events'):
      selected_events=self.config.get('input','events')
      print 'Selected events %s'%(str(selected_events))

      if selected_events=='all':
          selected_events=None
      else:
          selected_events=self.select_events()
    else:
        selected_events=None

    # No input file given, analyse the entire time stretch between gpsstart and gpsend
    if self.config.has_option('input','analyse-all-time') and self.config.getboolean('input','analyse-all-time')==True:
        print 'Setting up for analysis of continuous time stretch %f - %f'%(gpsstart,gpsend)
        seglen=self.config.getfloat('engine','seglen')
        if(self.config.has_option('input','segment-overlap')):
          overlap=self.config.getfloat('input','segment-overlap')
        else:
          overlap=32.;
        if(overlap>seglen):
          print 'ERROR: segment-overlap is greater than seglen'
          sys.exit(1)
        # Now divide gpsstart - gpsend into jobs of seglen - overlap length
        t=gpsstart
        events=[]
        while(t<gpsend):
            ev=Event(trig_time=t+seglen-2)
            ev.set_engine_option('segment-start',str(t-overlap))
            ev.set_engine_option('time-min',str(t))
            tMax=t + seglen - overlap
            if tMax>=gpsend:
                tMax=gpsend
            ev.set_engine_option('time-max',str(tMax))
            events.append(ev)
            t=tMax
        return events

    # ASCII list of GPS times
    if self.config.has_option('input','gps-time-file'):
      times=scan_timefile(self.config.get('input','gps-time-file'))
      events=[Event(trig_time=time) for time in times]
    # Siminspiral Table
    if self.config.has_option('input','injection-file'):
      from pylal import SimInspiralUtils
      injTable=SimInspiralUtils.ReadSimInspiralFromFiles([self.config.get('input','injection-file')])
      events=[Event(SimInspiral=inj) for inj in injTable]
      self.add_pfn_cache([create_pfn_tuple(self.config.get('input','injection-file'))])
    # SimBurst Table
    if self.config.has_option('input','burst-injection-file'):
      #from pylal import SimBurstUtils
      from glue.ligolw import lsctables
      from glue.ligolw import utils
      from glue.ligolw import ligolw
      injfile=self.config.get('input','burst-injection-file')
      class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass
      lsctables.use_in(LIGOLWContentHandler)
      injTable=lsctables.SimBurstTable.get_table(utils.load_filename(injfile,contenthandler = LIGOLWContentHandler))
      events=[Event(SimBurst=inj) for inj in injTable]
      self.add_pfn_cache([create_pfn_tuple(self.config.get('input','burst-injection-file'))])
    # SnglInspiral Table
    if self.config.has_option('input','sngl-inspiral-file'):
      from pylal import SnglInspiralUtils
      trigTable=SnglInspiralUtils.ReadSnglInspiralFromFiles([self.config.get('input','sngl-inspiral-file')])
      events=[Event(SnglInspiral=trig) for trig in trigTable]
      self.add_pfn_cache([create_pfn_tuple(self.config.get('input','sngl-inspiral-file'))])
    if self.config.has_option('input','coinc-inspiral-file'):
      from pylal import CoincInspiralUtils
      coincTable = CoincInspiralUtils.readCoincInspiralFromFiles([self.config.get('input','coinc-inspiral-file')])
      events = [Event(CoincInspiral=coinc) for coinc in coincTable]
      self.add_pfn_cache([create_pfn_tuple(self.config.get('input','coinc-inspiral-file'))])
    # LVAlert CoincInspiral Table
    if self.config.has_option('input','gid'):
      gid=self.config.get('input','gid')
      flow=40.0
      if self.config.has_option('lalinference','flow'):
        flow=min(ast.literal_eval(self.config.get('lalinference','flow')).values())
      downloadgracedbpsd=True
      if self.config.has_option('input','ignore-gracedb-psd'):
        if self.config.getboolean('input','ignore-gracedb-psd'):
          downloadgracedbpsd=False
      events = readLValert(gid=gid,flow=flow,gracedb=self.config.get('condor','gracedb'),basepath=self.basepath,downloadpsd=downloadgracedbpsd)
    else: gid=None
    # pipedown-database
    if self.config.has_option('input','gstlal-db'):
            queryfunc=get_zerolag_lloid
            dbname=self.config.get('input','gstlal-db')
    elif self.config.has_option('input','pipedown-db'):
            queryfunc=get_zerolag_pipedown
            dbname=self.config.get('input','pipedown-db')
    else: dbname=None
    if dbname:
      db_connection = open_pipedown_database(dbname,None)[0]
      # Timeslides
      if self.config.has_option('input','time-slide-dump'):
        timeslidedump=self.config.get('input','time-slide-dump')
      else:
        timeslidedump=None
      if self.config.has_option('input','min-cfar'):
        mincfar=self.config.getfloat('input','min-cfar')
      else:
        mincfar=-1
      if self.config.has_option('input','max-cfar'):
	maxcfar=self.config.getfloat('input','max-cfar')
      else:
	maxcfar=-1
      if self.config.get('input','timeslides').lower()=='true':
	events=get_timeslides_pipedown(db_connection, gpsstart=gpsstart, gpsend=gpsend,dumpfile=timeslidedump,max_cfar=maxcfar)
      else:
	events=queryfunc(db_connection, gpsstart=gpsstart, gpsend=gpsend, dumpfile=timeslidedump,max_cfar=maxcfar,min_cfar=mincfar)
    if(selected_events is not None):
        used_events=[]
        for i in selected_events:
            e=events[i]
            e.event_id=i
            used_events.append(e)
        events=used_events
    if gpsstart is not None:
        events = filter(lambda e: not e.trig_time<gpsstart, events)
    if gpsend is not None:
        events = filter(lambda e: not e.trig_time>gpsend, events)
    return events

  def add_full_analysis_lalinferencenest(self,event):
    """
    Generate an end-to-end analysis of a given event (Event class)
    For LALinferenceNest code. Uses parallel runs if specified
    """
    evstring=str(event.event_id)
    if event.trig_time is not None:
        evstring=str(event.trig_time)+'-'+str(event.event_id)
    if self.config.has_option('analysis','nparallel'):
        Npar=self.config.getint('analysis','nparallel')
    else:
        Npar=4
    # Set up the parallel engine nodes
    enginenodes=[]
    for i in range(Npar):
      n=self.add_engine_node(event)
      if n is not None:
        enginenodes.append(n)
        if i>0:
          n.add_var_arg('--dont-dump-extras')
    if len(enginenodes)==0:
      return False
    myifos=enginenodes[0].get_ifos()
    # Merge the results together
    pagedir=os.path.join(self.webdir,evstring,myifos)
    #pagedir=os.path.join(self.basepath,evstring,myifos)
    mkdirs(pagedir)
    mergenode=MergeNSNode(self.merge_job,parents=enginenodes)
    mergenode.set_pos_output_file(os.path.join(self.posteriorpath,'posterior_%s_%s.dat'%(myifos,evstring)))
    self.add_node(mergenode)
    # Call finalize to build final list of available data
    enginenodes[0].finalize()
    enginenodes[0].set_psd_files()
    enginenodes[0].set_snr_file()
    if self.config.getboolean('analysis','coherence-test') and len(enginenodes[0].ifos)>1:
        if self.site!='local':
          zipfilename='postproc_'+evstring+'.tar.gz'
        else:
          zipfilename=None
        respagenode=self.add_results_page_node(resjob=self.cotest_results_page_job,outdir=pagedir,parent=mergenode,gzip_output=zipfilename,ifos=enginenodes[0].ifos)
        respagenode.set_psd_files(enginenodes[0].get_psd_files())
        respagenode.set_snr_file(enginenodes[0].get_snr_file())
        if os.path.exists(self.basepath+'/coinc.xml'):
          respagenode.set_coinc_file(self.basepath+'/coinc.xml')
        mkdirs(os.path.join(self.basepath,'coherence_test'))
        par_mergenodes=[]
        for ifo in enginenodes[0].ifos:
            cotest_nodes=[self.add_engine_node(event,ifos=[ifo]) for i in range(Npar)]
            for co in cotest_nodes:
              co.set_psdstart(enginenodes[0].GPSstart)
              co.set_psdlength(enginenodes[0].psdlength)
              if co!=cotest_nodes[0]:
                co.add_var_arg('--dont-dump-extras')
              else:
                co.set_psd_files()
                co.set_snr_file()
            pmergenode=MergeNSNode(self.merge_job,parents=cotest_nodes)
            pmergenode.set_pos_output_file(os.path.join(self.posteriorpath,'posterior_%s_%s.dat'%(ifo,evstring)))
            self.add_node(pmergenode)
            par_mergenodes.append(pmergenode)
            presultsdir=os.path.join(pagedir,ifo)
            mkdirs(presultsdir)
            if self.site!='local':
              pzipfilename='postproc_'+evstring+'_'+ifo+'.tar.gz'
            else:
              pzipfilename=None
            subresnode=self.add_results_page_node(outdir=presultsdir,parent=pmergenode, gzip_output=pzipfilename,ifos=ifo)
            subresnode.set_psd_files(cotest_nodes[0].get_psd_files())
            subresnode.set_snr_file(cotest_nodes[0].get_snr_file())
            if os.path.exists(self.basepath+'/coinc.xml'):
              subresnode.set_coinc_file(self.basepath+'/coinc.xml')
            subresnode.set_bayes_coherent_noise(pmergenode.get_B_file())
            if self.config.has_option('input','injection-file') and event.event_id is not None:
                subresnode.set_injection(self.config.get('input','injection-file'),event.event_id)
            elif self.config.has_option('input','burst-injection-file') and event.event_id is not None:
                subresnode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
        coherence_node=CoherenceTestNode(self.coherence_test_job,outfile=os.path.join(self.basepath,'coherence_test','coherence_test_%s_%s.dat'%(myifos,evstring)))
        coherence_node.add_coherent_parent(mergenode)
        map(coherence_node.add_incoherent_parent, par_mergenodes)
        self.add_node(coherence_node)
        respagenode.add_parent(coherence_node)
        respagenode.set_bayes_coherent_incoherent(coherence_node.get_output_files()[0])
    else:
        if self.site!='local':
          zipfilename='postproc_'+evstring+'.tar.gz'
        else:
          zipfilename=None
        # Note: Disabled gzip_output for now. Possibly need it for future Pegasus use
        respagenode=self.add_results_page_node(outdir=pagedir,parent=mergenode,gzip_output=None,ifos=enginenodes[0].ifos)
        respagenode.set_psd_files(enginenodes[0].get_psd_files())
        respagenode.set_snr_file(enginenodes[0].get_snr_file())
        if os.path.exists(self.basepath+'/coinc.xml'):
          respagenode.set_coinc_file(self.basepath+'/coinc.xml')
    respagenode.set_bayes_coherent_noise(mergenode.get_B_file())
    if self.config.has_option('input','injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','injection-file'),event.event_id)
    elif self.config.has_option('input','burst-injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)

    if self.config.has_option('analysis','upload-to-gracedb'):
      if self.config.getboolean('analysis','upload-to-gracedb') and event.GID is not None:
        self.add_gracedb_start_node(event.GID,'LALInference',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
        self.add_gracedb_log_node(respagenode,event.GID)
      elif self.config.has_option('analysis','ugid'):
        # LIB will want to upload info to gracedb but if we pass the gid in the usual way the pipeline
        # will try to pull inspiral-only XML tables from the gdb page, failing.
        # To avoid that, LIB will read the gracedDB id to upload info to as an ugid=ID option
        # in the analysis section.
        ugid=self.config.get('analysis','ugid')
        self.add_gracedb_start_node(ugid,'LIB',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
        self.add_gracedb_log_node(respagenode,ugid)
    return True

  def add_full_analysis_lalinferencemcmc(self,event):
    """
    Generate an end-to-end analysis of a given event
    For LALInferenceMCMC.
    """
    evstring=str(event.event_id)
    if event.trig_time is not None:
        evstring=str(event.trig_time)+'-'+str(event.event_id)
    if self.config.has_option('analysis','nparallel'):
        Npar=self.config.getint('analysis','nparallel')
    else:
        Npar=2
    enginenodes=[]
    for i in range(Npar):
        enginenodes.append(self.add_engine_node(event))
        if i>0:
            enginenodes[-1].add_var_arg('--dont-dump-extras')
    myifos=enginenodes[0].get_ifos()
    enginenodes[0].set_psd_files()
    enginenodes[0].set_snr_file()
    pagedir=os.path.join(self.webdir,evstring,myifos)
    mkdirs(pagedir)
    respagenode=self.add_results_page_node(outdir=pagedir,ifos=enginenodes[0].ifos)
    respagenode.set_psd_files(enginenodes[0].get_psd_files())
    respagenode.set_snr_file(enginenodes[0].get_snr_file())
    if os.path.exists(self.basepath+'/coinc.xml'):
        respagenode.set_coinc_file(self.basepath+'/coinc.xml')
    if self.config.has_option('input','injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','injection-file'),event.event_id)
    if self.config.has_option('input','burst-injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
    map(respagenode.add_engine_parent, enginenodes)
    if event.GID is not None:
      if self.config.has_option('analysis','upload-to-gracedb'):
        if self.config.getboolean('analysis','upload-to-gracedb'):
          self.add_gracedb_start_node(event.GID,'LALInference',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
          self.add_gracedb_log_node(respagenode,event.GID)

  def add_full_analysis_lalinferencebambi(self,event):
    """
    Generate an end-to-end analysis of a given event
    For LALInferenceBAMBI.
    """
    evstring=str(event.event_id)
    if event.trig_time is not None:
        evstring=str(event.trig_time)+'-'+str(event.event_id)
    enginenodes=[]
    enginenodes.append(self.add_engine_node(event))
    myifos=enginenodes[0].get_ifos()
    enginenodes[0].set_psd_files()
    enginenodes[0].set_snr_file()
    pagedir=os.path.join(self.webdir,evstring,myifos)
    mkdirs(pagedir)
    respagenode=self.add_results_page_node(outdir=pagedir,ifos=enginenodes[0].ifos)
    respagenode.set_psd_files(enginenodes[0].get_psd_files())
    respagenode.set_snr_file(enginenodes[0].get_snr_file())
    if os.path.exists(self.basepath+'/coinc.xml'):
      respagenode.set_coinc_file(self.basepath+'/coinc.xml')
    respagenode.set_bayes_coherent_noise(enginenodes[0].get_B_file())
    respagenode.set_header_file(enginenodes[0].get_header_file())
    if self.config.has_option('input','injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','injection-file'),event.event_id)
    elif self.config.has_option('input','burst-injection-file') and event.event_id is not None:
        respagenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
    map(respagenode.add_engine_parent, enginenodes)
    if event.GID is not None:
      if self.config.has_option('analysis','upload-to-gracedb'):
        if self.config.getboolean('analysis','upload-to-gracedb'):
          self.add_gracedb_start_node(event.GID,'LALInference',[sciseg.get_df_node() for sciseg in enginenodes[0].scisegs.values()])
          self.add_gracedb_log_node(respagenode,event.GID)

  def add_science_segments(self):
    # Query the segment database for science segments and
    # add them to the pool of segments
    if self.config.has_option('input','ignore-science-segments'):
      if self.config.getboolean('input','ignore-science-segments'):
        start=self.config.getfloat('input','gps-start-time')
        end=self.config.getfloat('input','gps-end-time')
        i=0
        for ifo in self.ifos:
          sciseg=pipeline.ScienceSegment((i,start,end,end-start))
          df_node=self.get_datafind_node(ifo,self.frtypes[ifo],int(sciseg.start()),int(sciseg.end()))
          sciseg.set_df_node(df_node)
          self.segments[ifo].append(sciseg)
          i+=1
        return
    # Look up science segments as required
    segmentdir=os.path.join(self.basepath,'segments')
    mkdirs(segmentdir)
    curdir=os.getcwd()
    os.chdir(segmentdir)
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
    os.chdir(curdir)

  def get_datafind_node(self,ifo,frtype,gpsstart,gpsend):
    node=pipeline.LSCDataFindNode(self.datafind_job)
    node.set_observatory(ifo[0])
    node.set_type(frtype)
    node.set_start(gpsstart)
    node.set_end(gpsend)
    #self.add_node(node)
    return node

  def add_engine_node(self,event,ifos=None,extra_options=None):
    """
    Add an engine node to the dag. Will find the appropriate cache files automatically.
    Will determine the data to be read and the output file.
    Will use all IFOs known to the DAG, unless otherwise specified as a list of strings
    """
    if ifos is None and len(event.ifos)>0:
      ifos=event.ifos
    if ifos is None:
      ifos=self.ifos
    end_time=event.trig_time
    seglen=self.config.getfloat('engine','seglen')
    segstart=end_time+2-seglen
    segend=segstart+seglen
    myifos=set([])
    for ifo in ifos:
      for seg in self.segments[ifo]:
        if segstart >= seg.start() and segend < seg.end():
          myifos.add(ifo)
    ifos=myifos
    if len(ifos)==0:
      print 'No data found for time %f - %f, skipping'%(segstart,segend)
      return

    computeroqweightsnode={}
    bayeslinenode={}
    prenode=LALInferenceDataDumpNode(self.preengine_job)
    node=self.EngineNode(self.engine_jobs[tuple(ifos)])
    roqeventpath=os.path.join(self.preengine_job.roqpath,str(event.event_id)+'/')
    if self.config.has_option('condor','bayesline') or self.config.has_option('condor','computeroqweights'):
      mkdirs(roqeventpath)
    node.set_trig_time(end_time)
    prenode.set_trig_time(end_time)
    randomseed=random.randint(1,2**31)
    node.set_seed(randomseed)
    prenode.set_seed(randomseed)
    srate=0
    if event.srate:
      srate=event.srate
    if self.config.has_option('lalinference','srate'):
      srate=ast.literal_eval(self.config.get('lalinference','srate'))
    elif self.config.has_option('engine','seglen'):
      srate=ast.literal_eval(self.config.get('engine','srate'))
    if srate is not 0:
      node.set_srate(srate)
      prenode.set_srate(srate)
    if event.trigSNR:
      node.set_trigSNR(event.trigSNR)
    if self.dataseed:
      node.set_dataseed(self.dataseed+event.event_id)
      prenode.set_dataseed(self.dataseed+event.event_id)
    gotdata=0
    for ifo in ifos:
      if event.timeslides.has_key(ifo):
        slide=event.timeslides[ifo]
      else:
        slide=0
      for seg in self.segments[ifo]:
        if segstart >= seg.start() and segend < seg.end():
            if not self.config.has_option('lalinference','fake-cache'):
              if self.config.has_option('condor','bayesline') or self.config.has_option('condor','computeroqweights'):
                prenode.add_ifo_data(ifo,seg,self.channels[ifo],timeslide=slide)
              gotdata+=node.add_ifo_data(ifo,seg,self.channels[ifo],timeslide=slide)
            else:
              fakecachefiles=ast.literal_eval(self.config.get('lalinference','fake-cache'))
              if self.config.has_option('condor','bayesline') or self.config.has_option('condor','computeroqweights'):
                prenode.add_fake_ifo_data(ifo,seg,fakecachefiles[ifo],self.channels[ifo],timeslide=slide)
              gotdata+=node.add_fake_ifo_data(ifo,seg,fakecachefiles[ifo],self.channels[ifo],timeslide=slide)
    if self.config.has_option('lalinference','psd-xmlfile'):
      psdpath=os.path.realpath(self.config.get('lalinference','psd-xmlfile'))
      node.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=end_time)
      prenode.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=end_time)
      if len(ifos)==0:
        node.ifos=node.cachefiles.keys()
        prenode.ifos=prenode.cachefiles.keys()
      else:
        node.ifos=ifos
        prenode.ifos=ifos
      gotdata=1
    if self.config.has_option('input','gid'):
      if os.path.isfile(os.path.join(self.basepath,'psd.xml.gz')):
        psdpath=os.path.join(self.basepath,'psd.xml.gz')
        node.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=None)
        prenode.psds=get_xml_psds(psdpath,ifos,os.path.join(self.basepath,'PSDs'),end_time=None)
    for ifo in ifos:
      prenode.flows[ifo]=str(40.0)
    if self.config.has_option('lalinference','flow'):
      node.flows=ast.literal_eval(self.config.get('lalinference','flow'))
      prenode.flows=ast.literal_eval(self.config.get('lalinference','flow'))
    if event.fhigh:
      for ifo in ifos:
        node.fhighs[ifo]=str(event.fhigh)
        prenode.fhighs[ifo]=str(event.fhigh)
    if self.config.has_option('lalinference','fhigh'):
      node.fhighs=ast.literal_eval(self.config.get('lalinference','fhigh'))
      prenode.fhighs=ast.literal_eval(self.config.get('lalinference','fhigh'))
    prenode.set_max_psdlength(self.config.getint('input','max-psd-length'))
    prenode.set_padding(self.config.getint('input','padding'))
    #prenode[ifo].set_output_file('/dev/null')
    prenode.add_var_arg('--outfile '+roqeventpath+'data-dump')
    prenode.add_var_arg('--data-dump')
    if self.config.has_option('lalinference','seglen'):
      prenode.set_seglen(self.config.getint('lalinference','seglen'))
    elif self.config.has_option('engine','seglen'):
      prenode.set_seglen(self.config.getint('engine','seglen'))
    else:
      prenode.set_seglen(event.duration)
    # Add the nodes it depends on
    for ifokey, seg in node.scisegs.items():
      dfnode=seg.get_df_node()

      if 1==1:
        #if dfnode is not None and dfnode not in self.get_nodes():
        if self.config.has_option('condor','bayesline') or self.config.has_option('condor','computeroqweights'):
          #if gotdata and seg.id() not in self.prenodes.keys():
          if gotdata and event.event_id not in self.prenodes.keys():
            if prenode not in self.get_nodes():
              self.add_node(prenode)
              for ifo in ifos:
                if self.config.has_option('condor','computeroqweights'):
                  computeroqweightsnode[ifo]=self.add_rom_weights_node(ifo,prenode)
                #self.add_node(computeroqweightsnode[ifo])
                if self.config.has_option('input','injection-file'):
                  freqDataFile=os.path.join(roqeventpath,'data-dump'+ifo+'-freqDataWithInjection.dat')
                else:
                  freqDataFile=os.path.join(roqeventpath,'data-dump'+ifo+'-freqData.dat')
                prenode.add_output_file(freqDataFile)
                prenode.add_output_file(os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                if self.config.has_option('condor','bayesline'):
                  bayeslinenode[ifo]=self.add_bayesline_node(ifo,prenode)
                  bayeslinenode[ifo].add_var_arg('-i '+freqDataFile)
                  bayeslinenode[ifo].add_input_file(freqDataFile)
                  bayeslinenode[ifo].add_var_arg('-o '+os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                  bayeslinenode[ifo].add_output_file(os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                if self.config.has_option('condor','computeroqweights'):
                  computeroqweightsnode[ifo].add_var_arg('-d '+freqDataFile)
                  computeroqweightsnode[ifo].add_input_file(freqDataFile)
                  if self.config.has_option('condor','bayesline'):
                    computeroqweightsnode[ifo].add_var_arg('-p '+os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                    computeroqweightsnode[ifo].add_input_file(os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat'))
                    computeroqweightsnode[ifo].add_parent(bayeslinenode[ifo])
                  else:
                    computeroqweightsnode[ifo].add_var_arg('-p '+os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                    computeroqweightsnode[ifo].add_input_file(os.path.join(roqeventpath,'data-dump'+ifo+'-PSD.dat'))
                  computeroqweightsnode[ifo].add_var_arg('-o '+roqeventpath)
                  computeroqweightsnode[ifo].add_output_file(os.path.join(roqeventpath,'weights_'+ifo+'.dat'))
              #self.prenodes[seg.id()]=(prenode,computeroqweightsnode)
              if self.config.has_option('condor','bayesline'):
                  self.prenodes[event.event_id]=(prenode,bayeslinenode)
              if self.config.has_option('condor','computeroqweights'):
                  self.prenodes[event.event_id]=(prenode,computeroqweightsnode)

      if self.config.has_option('condor','bayesline') or self.config.has_option('condor','computeroqweights'):
        #node.add_parent(self.prenodes[seg.id()][1][ifokey])
        node.add_parent(self.prenodes[event.event_id][1][ifokey])

      if dfnode is not None and dfnode not in self.get_nodes():
        if not self.config.has_option('lalinference','fake-cache'):
          self.add_node(dfnode)

    if gotdata:
      self.add_node(node)
    else:
      'Print no data found for time %f'%(end_time)
      return None
    if extra_options is not None:
      for opt in extra_options.keys():
	    node.add_var_arg('--'+opt+' '+extra_options[opt])
    # Add control options
    if self.config.has_option('input','injection-file'):
       node.set_injection(self.config.get('input','injection-file'),event.event_id)
       prenode.set_injection(self.config.get('input','injection-file'),event.event_id)
    if self.config.has_option('input','burst-injection-file'):
       node.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
       prenode.set_injection(self.config.get('input','burst-injection-file'),event.event_id)
    if self.config.has_option('lalinference','seglen'):
      node.set_seglen(self.config.getint('lalinference','seglen'))
    elif  self.config.has_option('engine','seglen'):
      node.set_seglen(self.config.getint('engine','seglen'))
    else:
      node.set_seglen(event.duration)
    if self.config.has_option('input','psd-length'):
      node.set_psdlength(self.config.getint('input','psd-length'))
    if self.config.has_option('input','psd-start-time'):
      node.set_psdstart(self.config.getfloat('input','psd-start-time'))
    node.set_max_psdlength(self.config.getint('input','max-psd-length'))
    node.set_padding(self.config.getint('input','padding'))
    out_dir=os.path.join(self.basepath,'engine')
    mkdirs(out_dir)
    node.set_output_file(os.path.join(out_dir,node.engine+'-'+str(event.event_id)+'-'+node.get_ifos()+'-'+str(node.get_trig_time())+'-'+str(node.id)))
    if self.config.has_option('condor','computeroqweights'):
      for ifo in ifos:
        node.add_var_arg('--'+ifo+'-roqweightsLinear '+os.path.join(roqeventpath,'weights_linear_'+ifo+'.dat'))
        node.add_input_file(os.path.join(roqeventpath,'weights_linear_'+ifo+'.dat'))
        node.add_var_arg('--'+ifo+'-roqweightsQuadratic '+os.path.join(roqeventpath,'weights_quadratic_'+ifo+'.dat'))
        node.add_input_file(os.path.join(roqeventpath,'weights_quadratic_'+ifo+'.dat'))
      node.add_var_arg('--roqtime_steps '+os.path.join(roqeventpath,'roq_sizes.dat'))
      node.add_input_file(os.path.join(roqeventpath,'roq_sizes.dat'))
      node.add_var_arg('--roq-times '+os.path.join(roqeventpath,'tcs.dat'))
      node.add_input_file(os.path.join(roqeventpath,'tcs.dat'))
      node.add_var_arg('--roqnodesLinear '+os.path.join(roqeventpath,'fnodes_linear.dat'))
      node.add_input_file(os.path.join(roqeventpath,'fnodes_linear.dat'))
      node.add_var_arg('--roqnodesQuadratic '+os.path.join(roqeventpath,'fnodes_quadratic.dat'))
      node.add_input_file(os.path.join(roqeventpath,'fnodes_quadratic.dat'))
    if self.config.has_option('condor','bayesline'):
      for ifo in ifos:
        node.psds[ifo]=os.path.join(roqeventpath,'BayesLine_PSD_'+ifo+'.dat')
    for (opt,arg) in event.engine_opts.items():
        node.add_var_opt(opt,arg)
    return node

  def add_results_page_node(self,resjob=None,outdir=None,parent=None,extra_options=None,gzip_output=None,ifos=None):
    if resjob is None:
        resjob=self.results_page_job
    node=ResultsPageNode(resjob)
    if parent is not None:
      node.add_parent(parent)
      infile=parent.get_pos_file()
      node.add_file_arg(infile)
    node.set_output_path(outdir)
    if gzip_output is not None:
        node.set_gzip_output(gzip_output)
    if ifos is not None:
      if isinstance(ifos,list):
        pass
      else:
        ifos=[ifos]
      node.set_ifos(ifos)

    self.add_node(node)
    return node

  def add_gracedb_start_node(self,gid,name='',parent=None):

    node=GraceDBNode(self.gracedbjob,parent=parent,gid=gid,command='log',tag='pe')
    node.set_message(name+' online parameter estimation started.')
    self.add_node(node)
    return node

  def add_gracedb_log_node(self,respagenode,gid):
    nodes=[]
    node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='log',tag='pe')
    resurl=respagenode.webpath.replace(self.gracedbjob.basepath,self.gracedbjob.baseurl)
    #node.set_message('online parameter estimation results:  '+resurl+'/posplots.html')
    node.set_message("LALInference online parameter estimation finished. <a href="+resurl+"/posplots.html>results</a>")
    self.add_node(node)
    nodes.append(node)

    tag='pe'
    if self.config.has_option('analysis','add-lvem-tag'):
      if self.config.getboolean('analysis','add-lvem-tag'):
        tag='pe,lvem'

    node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag=tag)
    node.set_filename(respagenode.webpath+'/corner/extrinsic.png')
    self.add_node(node)
    nodes.append(node)

    node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag='pe')
    node.set_filename(respagenode.webpath+'/corner/intrinsic.png')
    self.add_node(node)
    nodes.append(node)

    node=GraceDBNode(self.gracedbjob,parent=respagenode,gid=gid,command='upload',tag='pe')
    node.set_filename(respagenode.webpath+'/corner/sourceFrame.png')
    self.add_node(node)
    nodes.append(node)

    return nodes

  def add_gracedb_FITSskymap_upload(self,event,engine=None):
    gid=event.GID
    if gid is None:
      return
    if engine=='lalinferenceburst':
      prefix='LIB'
    elif engine is None:
      prefix=""
    else:
      prefix='LALInference'
    nodes=None
    if self.config.has_option('condor','skyarea'):
      if self.config.has_option('analysis','upload-to-gracedb'):
        if self.config.getboolean('analysis','upload-to-gracedb'):
          tag='sky_loc'
          if self.config.has_option('analysis','add-lvem-tag'):
            if self.config.getboolean('analysis','add-lvem-tag'):
              tag='sky_loc,lvem'
          skynodes=filter(lambda x: isinstance(x,SkyAreaNode) ,self.get_nodes())
          nodes=[]
          for sk in skynodes:
            if len(sk.ifos)>1:
              node=GraceDBNode(self.gracedbjob,parent=sk,gid=gid,tag=tag)
              #for p in sk.__parents:
              #  if isinstance(p,ResultPageNode):
              #    resultpagenode=p
              node.set_filename(sk.outdir+'/%s_skymap.fits.gz'%prefix)
              node.set_message('%s FITS sky map'%prefix)
              self.add_node(node)
              nodes.append(node)
    return nodes

  def add_rom_weights_node(self,ifo,parent=None):
    #try:
    #node=self.computeroqweightsnodes[ifo]
        #except KeyError:
    node=ROMNode(self.computeroqweights_job,ifo,parent.seglen,parent.flows[ifo])
    self.computeroqweightsnodes[ifo]=node
    if parent is not None:
      node.add_parent(parent)
    self.add_node(node)
    return node

  def add_bayesline_node(self,ifo,parent=None):
    node=BayesLineNode(self.bayesline_job)
    self.bayeslinenodes[ifo]=node
    if parent is not None:
      node.add_parent(parent)
    self.add_node(node)
    return node


class EngineJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
  def __init__(self,cp,submitFile,logdir,engine,ispreengine=False,dax=False,site=None):
    self.ispreengine=ispreengine
    self.engine=engine
    basepath=cp.get('paths','basedir')
    if ispreengine is True:
      roqpath=os.path.join(basepath,'ROQdata')
      self.roqpath=roqpath
      mkdirs(roqpath)
      exe=cp.get('condor',self.engine)
      universe='vanilla'
    else:
      if self.engine=='lalinferencemcmc':
        exe=cp.get('condor','mpiwrapper')
        universe="vanilla"
      elif self.engine=='lalinferencebambimpi':
        exe=cp.get('condor','mpiwrapper')
        universe="vanilla"
      elif self.engine=='lalinferencenest' or self.engine=='lalinferenceburst' or self.engine=='lalinferencebambi':
        exe=cp.get('condor',self.engine)
        if site is not None and site!='local':
          universe='vanilla'
        else:
          # Run in the vanilla universe when using resume
          if cp.has_option('engine','resume'):
            universe='vanilla'
          else:
            # This will go away, since condor_compile will not be supported in future
            universe='standard'
      else:
        print 'LALInferencePipe: Unknown engine node type %s!'%(self.engine)
        sys.exit(1)

    pipeline.CondorDAGJob.__init__(self,universe,exe)
    pipeline.AnalysisJob.__init__(self,cp,dax=dax)
    if cp.has_option('analysis','accounting_group'):
      self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    if cp.has_option('analysis','accounting_group_user'):
      self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
    try:
      hostname=socket.gethostbyaddr(socket.gethostname())[0]
    except:
      hostname='Unknown'
    if cp.has_option('condor','queue'):
      self.add_condor_cmd('+'+cp.get('condor','queue'),'True')
      self.add_condor_cmd('Requirements','(TARGET.'+cp.get('condor','queue')+' =?= True)')
    # Set grid site if needed
    if cp.has_option('engine','resume'):
        self.resume=True
    else:
        self.resume=False
    if site:
      self.set_grid_site(site)
      if site!='local':
        self.set_executable_installed(False)
    # Set the options which are always used
    self.set_sub_file(os.path.abspath(submitFile))
    self.add_condor_cmd('getenv','true')
    if self.engine=='lalinferencemcmc' or self.engine=='lalinferencebambimpi':
      self.binary=cp.get('condor',self.engine.replace('mpi',''))
      self.mpirun=cp.get('condor','mpirun')
      if cp.has_section('mpi'):
        if ispreengine is False:
          self.machine_count=cp.get('mpi','machine-count')
          self.machine_memory=cp.get('mpi','machine-memory')
          if cp.has_option('mpi','mpi_task_count'):
            self.mpi_task_count=cp.get('mpi','mpi_task_count')
          else:
            self.mpi_task_count=self.machine_count
        else:
          self.machine_count=str(1)
          self.machine_memory=cp.get('mpi','machine-memory')
          self.mpi_task_count=self.machine_count
      else:
        self.machine_count=str(1)
        self.machine_memory=str(1024) # default value if the user did not specify something
        self.mpi_task_count=self.machine_count
      #self.add_condor_cmd('machine_count',machine_count)
      #self.add_condor_cmd('environment','CONDOR_MPI_PATH=%s'%(openmpipath))
      if hostname=='pcdev1.phys.uwm.edu':
        self.add_condor_cmd('Requirements','CAN_RUN_MULTICORE')
        self.add_condor_cmd('+RequiresMultipleCores','True')
      self.add_condor_cmd('request_cpus',self.machine_count)
      self.add_condor_cmd('request_memory',str(float(self.machine_count)*float(self.machine_memory)))
    if cp.has_section(self.engine):
      if not ispreengine:
        self.add_ini_opts(cp,self.engine)
    if  cp.has_section('engine'):
      if not ispreengine:
        self.add_ini_opts(cp,'engine')
    self.set_stdout_file(os.path.join(logdir,'lalinference-$(cluster)-$(process)-$(node).out'))
    self.set_stderr_file(os.path.join(logdir,'lalinference-$(cluster)-$(process)-$(node).err'))
    # For LALInferenceNest demand only 1 thread (to be tuned later)
    if self.engine=='lalinferencenest':
            self.add_condor_cmd('environment','OMP_NUM_THREADS=1')

  def set_grid_site(self,site=None):
    """
    Over-load base class method to choose condor universe properly
    """
    if self.engine=='lalinferencenest' or self.engine=='lalinferenceburst' or self.engine=='lalinferencebambi':
      if site is not None and site!='local':
        self.set_universe('vanilla')
      else:
        if self.resume:
          self.set_universe('vanilla')
        else:
          self.set_universe('standard')
    else:
      self.set_universe('vanilla')
    pipeline.CondorDAGJob.set_grid_site(self,site)


class EngineNode(pipeline.CondorDAGNode):
  new_id = itertools.count().next
  def __init__(self,li_job):
    pipeline.CondorDAGNode.__init__(self,li_job)
    self.ifos=[]
    self.scisegs={}
    self.channels={}
    self.psds={}
    self.flows={}
    self.fhighs={}
    self.timeslides={}
    self.seglen=None
    self.psdlength=None
    self.padding=None
    self.maxlength=None
    self.psdstart=None
    self.snrfile=None
    self.psdfiles=None
    self.cachefiles={}
    if li_job.ispreengine is False:
      self.id=EngineNode.new_id()
    self.__finaldata=False
    self.fakedata=False
    self.lfns=[] # Local file names (for frame files and pegasus)

  def set_seglen(self,seglen):
    self.seglen=seglen

  def set_psdlength(self,psdlength):
    self.psdlength=psdlength

  def set_max_psdlength(self,psdlength):
    self.maxlength=psdlength

  def set_padding(self,padding):
    self.padding=padding

  def set_psdstart(self,psdstart):
    self.psdstart=psdstart

  def set_seed(self,seed):
    self.add_var_opt('randomseed',str(seed))

  def set_srate(self,srate):
    self.add_var_opt('srate',str(srate))

  def set_trigSNR(self,trigSNR):
    self.add_var_opt('trigger-snr',str(trigSNR))

  def set_dataseed(self,seed):
    self.add_var_opt('dataseed',str(seed))

  def get_ifos(self):
    return ''.join(map(str,self.ifos))

  def set_psd_files(self):
    #if node.psdspath is not None:
    try:
      #bambi store the outfile in fileroot
      pathroot=self.fileroot
    except:
      #mcmc and nest in postfile
      pathroot=self.posfile
      if pathroot[-3:]=='.00':
        pathroot=pathroot[:-3]
    st=""
    for i in self.ifos:
      tmpst="%s%s-PSD.dat,"%(pathroot,i)
      st+=tmpst
      self.add_output_file(tmpst[:-1])
    st=st[:-1]
    self.psdfiles=st

  def get_psd_files(self):
    return self.psdfiles

  def set_snr_file(self):
    try:
      #bambi store the outfile in fileroot
      pathroot=self.fileroot
    except:
      #mcmc and nest in postfile
      pathroot=self.posfile
      if pathroot[-3:]=='.00':
        pathroot=pathroot[:-3]
    st="%s_snr.txt"%pathroot
    self.add_output_file(st)
    self.snrfile=st

  def get_snr_file(self):
    return self.snrfile

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

  def set_injection(self,injfile,event):
    """
    Set a software injection to be performed.
    """
    self.add_file_opt('inj',injfile)
    self.set_event_number(event)

  def get_trig_time(self): return self.__trigtime

  def add_fake_ifo_data(self,ifo,sciseg,fake_cache_name,fake_channel_name,timeslide=0):
    """
    Dummy method to set up fake data without needing to run datafind
    """
    self.ifos.append(ifo)
    self.scisegs[ifo]=sciseg
    self.cachefiles[ifo]=fake_cache_name
    self.timeslides[ifo]=timeslide
    self.channels[ifo]=fake_channel_name
    self.fakedata=True
    return 1

  def add_ifo_data(self,ifo,sciseg,channelname,timeslide=0):
    self.ifos.append(ifo)
    self.scisegs[ifo]=sciseg
    parent=sciseg.get_df_node()
    if parent is not None:
        self.add_parent(parent)
        df_output=parent.get_output()
        self.set_cache(df_output,ifo)
        #self.cachefiles[ifo]=parent.get_output_files()[0]
        #self.add_input_file(self.cachefiles[ifo])
        self.timeslides[ifo]=timeslide
        self.channels[ifo]=channelname
        return 1
    else: return 0

  def set_cache(self,filename,ifo):
    """
    Add a cache file from LIGODataFind. Based on same method from pipeline.AnalysisNode
    """
    #print 'Adding cache files %s'%(str(filename))
    if isinstance(filename,str): # A normal lal cache file
        self.cachefiles[ifo]=filename
        self.add_input_file(filename)
    elif isinstance(filename,list): # A list of LFNs (for DAX mode)
        self.add_var_opt('glob-frame-data',' ')
        if len(filename) == 0:
          raise pipeline.CondorDAGNodeError, \
          "LDR did not return any LFNs for query: check ifo and frame type"
        for lfn in filename:
            self.lfns.append(lfn)

  def finalize(self):
    if not self.__finaldata:
      self._finalize_ifo_data()
    pipeline.CondorDAGNode.finalize(self)

  def _finalize_ifo_data(self):
      """
      Add final list of IFOs and data to analyse to command line arguments.
      """
      for ifo in self.ifos:
        self.add_var_arg('--ifo '+ifo)
        if self.fakedata:
            self.add_var_opt('%s-cache'%(ifo),self.cachefiles[ifo])
        elif not self.lfns:
            self.add_file_opt('%s-cache'%(ifo),self.cachefiles[ifo])
        self.add_var_opt('%s-channel'%(ifo),self.channels[ifo])
        if self.flows: self.add_var_opt('%s-flow'%(ifo),self.flows[ifo])
        if self.fhighs: self.add_var_opt('%s-fhigh'%(ifo),self.fhighs[ifo])
        if self.psds:
            self.add_var_opt('%s-psd'%(ifo),self.psds[ifo])
            #self.add_input_file(self.psds[ifo])
        if any(self.timeslides): self.add_var_opt('%s-timeslide'%(ifo),self.timeslides[ifo])

      # Start at earliest common time
      # NOTE: We perform this arithmetic for all ifos to ensure that a common data set is
      # Used when we are running the coherence test.
      # Otherwise the noise evidence will differ.
      #if self.scisegs!={}:
      starttime=max([int(self.scisegs[ifo].start()) for ifo in self.ifos])
      endtime=min([int(self.scisegs[ifo].end()) for ifo in self.ifos])
      #else:
      #  (starttime,endtime)=self.get_required_data(self.get_trig_time())
      #  starttime=floor(starttime)
      #  endtime=ceil(endtime)
        #starttime=self.get_trig_time()-self.padding-self.seglen-self.psdlength#-0.5*self.maxlength
        #endtime=starttime+self.padding#+self.maxlength
      self.GPSstart=starttime
      self.__GPSend=endtime
      length=endtime-starttime

      # Now we need to adjust the start time and length to make sure the maximum data length
      # is not exceeded.
      trig_time=self.get_trig_time()
      maxLength=self.maxlength
      if(length > maxLength):
        while(self.GPSstart+maxLength<trig_time and self.GPSstart+maxLength<self.__GPSend):
          self.GPSstart+=maxLength/2.0
      # Override calculated start time if requested by user in ini file
      if self.psdstart is not None:
        self.GPSstart=self.psdstart
        #print 'Over-riding start time to user-specified value %f'%(self.GPSstart)
        #if self.GPSstart<starttime or self.GPSstart>endtime:
        #  print 'ERROR: Over-ridden time lies outside of science segment!'
        #  raise Exception('Bad psdstart specified')
      self.add_var_opt('psdstart',str(self.GPSstart))
      if self.psdlength is None:
        self.psdlength=self.__GPSend-self.GPSstart-2*self.padding-self.seglen-1
        if(self.psdlength>self.maxlength):
          self.psdlength=self.maxlength
      self.add_var_opt('psdlength',self.psdlength)
      self.add_var_opt('seglen',self.seglen)
      for lfn in self.lfns:
        a, b, c, d = lfn.split('.')[0].split('-')
        t_start = int(c)
        t_end = int(c) + int(d)
        data_end=max(self.GPSstart+self.psdlength,trig_time+2)
        if( t_start <= data_end and t_end>self.GPSstart):
        #if (t_start <= (self.GPSstart+self.psdlength or t_start <=trig_time+2 or t_end >=) \
        #    and ( (t_end <= (self.GPSstart+self.psdlength )) or (t_end <= trig_time+2) ))  :
            self.add_input_file(lfn)
      self.__finaldata=True

class LALInferenceNestNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencenest'
    self.outfilearg='outfile'

  def set_output_file(self,filename):
    self.nsfile=filename+'.dat'
    self.posfile=self.nsfile
    self.add_file_opt(self.outfilearg,self.nsfile,file_is_output_file=True)
    self.Bfilename=self.nsfile+'_B.txt'
    self.add_output_file(self.Bfilename)
    self.headerfile=self.nsfile+'_params.txt'
    self.add_output_file(self.headerfile)
    if self.job().resume:
        self.add_checkpoint_file(self.nsfile+'_resume')

  def get_B_file(self):
    return self.Bfilename

  def get_ns_file(self):
    return self.nsfile

  def get_header_file(self):
    return self.headerfile

class LALInferenceBurstNode(EngineNode,LALInferenceNestNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferenceburst'
    self.outfilearg='outfile'
  def set_injection(self,injfile,event):
    """
    Set a software injection to be performed.
    """
    self.add_file_opt('binj',injfile)
    self.set_event_number(event)

class LALInferenceMCMCNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencemcmc'
    self.outfilearg='outfile'
    self.add_var_opt('mpirun',li_job.mpirun)
    self.add_var_opt('np',str(li_job.mpi_task_count))
    self.add_var_opt('executable',li_job.binary)

  def set_output_file(self,filename):
    self.posfile=filename+'.h5'
    # Should also take care of the higher temperature outpufiles with
    # self.add_output_file, getting the number of files from machine_count
    self.add_file_opt(self.outfilearg,self.posfile,file_is_output_file=True)

  def get_pos_file(self):
    return self.posfile

class LALInferenceDataDumpNode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencedatadump'
    self.outfilearg='outfile'
  def set_output_file(self,filename):
    pass

class LALInferenceBAMBINode(EngineNode):
  def __init__(self,li_job):
    EngineNode.__init__(self,li_job)
    self.engine='lalinferencebambi'
    self.outfilearg='outfile'
    if li_job.engine=='lalinferencebambimpi':
      self.add_var_opt('mpirun',li_job.mpirun)
      self.add_var_opt('np',str(li_job.machine_count))
      self.add_var_opt('executable',li_job.binary)

  def set_output_file(self,filename):
    self.fileroot=filename+'_'
    self.posfile=self.fileroot+'post_equal_weights.dat'
    self.paramsfile=self.fileroot+'params.txt'
    self.Bfilename=self.fileroot+'evidence.dat'
    self.headerfile=self.paramsfile
    self.add_file_opt(self.outfilearg,self.fileroot,file_is_output_file=True)
    self.add_output_file(self.Bfilename)
    self.add_output_file(self.posfile)
    self.add_output_file(self.headerfile)

  def get_B_file(self):
    return self.Bfilename

  def get_pos_file(self):
    return self.posfile

  def get_header_file(self):
    return self.headerfile

class ResultsPageJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
  def __init__(self,cp,submitFile,logdir,dax=False):
    exe=cp.get('condor','resultspage')
    pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
    pipeline.AnalysisJob.__init__(self,cp,dax=dax) # Job always runs locally
    if cp.has_option('analysis','accounting_group'):
      self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    if cp.has_option('analysis','accounting_group_user'):
      self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
    self.set_sub_file(os.path.abspath(submitFile))
    self.set_stdout_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).out'))
    self.set_stderr_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).err'))
    self.add_condor_cmd('getenv','True')
    self.add_condor_cmd('request_memory','2000')
    self.add_ini_opts(cp,'resultspage')
    # self.add_opt('Nlive',cp.get('analysis','nlive'))

    if cp.has_option('results','skyres'):
        self.add_opt('skyres',cp.get('results','skyres'))

class ResultsPageNode(pipeline.CondorDAGNode):
    def __init__(self,results_page_job,outpath=None):
        pipeline.CondorDAGNode.__init__(self,results_page_job)
        if outpath is not None:
            self.set_output_path(path)
        self.__event=0
        self.ifos=None
        self.injfile=None
    def set_gzip_output(self,path):
        self.add_file_opt('archive',path,file_is_output_file=True)
    def set_output_path(self,path):
        self.webpath=path
        #self.add_file_opt('outpath',path,file_is_output_file=True)
        self.add_var_opt('outpath',path)
        #self.add_file_opt('archive','results.tar.gz',file_is_output_file=True)
        mkdirs(path)
        self.posfile=os.path.join(path,'posterior_samples.dat')
        self.add_output_file(self.posfile)
    def get_output_path(self):
        return self.webpath
    def set_injection(self,injfile,eventnumber):
        self.injfile=injfile
        self.add_file_opt('inj',injfile)
        self.set_event_number(eventnumber)
    def get_injection(self):
        return self.injfile
    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_arg('--eventnum '+str(event))

    def get_event_number(self):
      return self.__event
    def set_psd_files(self,st):
      if st is None:
        return
      for i in st.split(','):
        self.add_input_file(i)
      self.add_var_arg('--psdfiles %s'%st)

    def set_snr_file(self,st):
      if st is None:
        return
      self.add_file_opt('snr',st)

    def set_coinc_file(self,coinc):
      if coinc is None:
        return
      self.add_var_arg('--trig '+coinc)

    def add_engine_parent(self,node):
      """
      Add a parent node which is one of the engine nodes
      And automatically set options accordingly
      """
      self.add_parent(node)
      self.add_file_arg(node.get_pos_file())

    def get_pos_file(self): return self.posfile
    def set_bayes_coherent_incoherent(self,bcifile):
        self.add_file_opt('bci',bcifile)
    def set_bayes_coherent_noise(self,bsnfile):
        self.add_file_opt('bsn',bsnfile)
    def set_header_file(self,headerfile):
        self.add_file_opt('header',headerfile)
    def set_ifos(self,ifos):
        self.ifos=ifos

class CoherenceTestJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class defining the coherence test job to be run as part of a pipeline.
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
      exe=cp.get('condor','coherencetest')
      pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
      pipeline.AnalysisJob.__init__(self,cp,dax=dax)
      if cp.has_option('analysis','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
      if cp.has_option('analysis','accounting_group_user'):
        self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
      self.add_opt('coherent-incoherent','')
      self.add_condor_cmd('getenv','True')
      self.set_stdout_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).out'))
      self.set_stderr_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).err'))
      self.set_sub_file(os.path.abspath(submitFile))

class CoherenceTestNode(pipeline.CondorDAGNode):
    """
    Class defining the node for the coherence test
    """
    def __init__(self,coherencetest_job,outfile=None):
      pipeline.CondorDAGNode.__init__(self,coherencetest_job)
      self.incoherent_parents=[]
      self.coherent_parent=None
      self.finalized=False
      if outfile is not None:
        self.add_file_opt('outfile',outfile,file_is_output_file=True)

    def add_coherent_parent(self,node):
      """
      Add a parent node which is an engine node, and process its outputfiles
      """
      self.coherent_parent=node
      self.add_parent(node)
    def add_incoherent_parent(self,node):
      """
      Add a parent node which provides one of the single-ifo evidence values
      """
      self.incoherent_parents.append(node)
      self.add_parent(node)
    def finalize(self):
      """
      Construct command line
      """
      if self.finalized==True: return
      self.finalized=True
      self.add_file_arg(self.coherent_parent.get_B_file())
      for inco in self.incoherent_parents:
        self.add_file_arg(inco.get_B_file())

class MergeNSJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class defining a job which merges several parallel nested sampling jobs into a single file
    Input arguments:
    cp        - A configparser object containing the setup of the analysis
    submitFile    - Path to store the submit file
    logdir        - A directory to hold the stderr, stdout files of the merge runs
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
      exe=cp.get('condor','mergescript')
      pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
      pipeline.AnalysisJob.__init__(self,cp,dax=dax)
      if cp.has_option('analysis','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
      if cp.has_option('analysis','accounting_group_user'):
        self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
      self.set_sub_file(os.path.abspath(submitFile))
      self.set_stdout_file(os.path.join(logdir,'merge-$(cluster)-$(process).out'))
      self.set_stderr_file(os.path.join(logdir,'merge-$(cluster)-$(process).err'))
      self.add_condor_cmd('getenv','True')
      if cp.has_option('engine','nlive'):
        self.add_opt('Nlive',cp.get('engine','nlive'))
      elif cp.has_option('engine','Nlive'):
        self.add_opt('Nlive',cp.get('engine','Nlive'))
      if cp.has_option('merge','npos'):
      	self.add_opt('npos',cp.get('merge','npos'))


class MergeNSNode(pipeline.CondorDAGNode):
    """
    Class defining the DAG node for a merge job
    Input arguments:
    merge_job = A MergeJob object
    parents = iterable of parent LALInferenceNest nodes (must have get_ns_file() method)
    """
    def __init__(self,merge_job,parents=None):
        pipeline.CondorDAGNode.__init__(self,merge_job)
        if parents is not None:
          for parent in parents:
            self.add_engine_parent(parent)

    def add_engine_parent(self,parent):
        self.add_parent(parent)
        self.add_file_arg(parent.get_ns_file())
        self.add_file_opt('headers',parent.get_header_file())
        self.add_input_file(parent.get_B_file())

    def set_pos_output_file(self,file):
        self.add_file_opt('pos',file,file_is_output_file=True)
        self.posfile=file
        self.Bfilename=self.posfile+'_B.txt'
        self.add_output_file(self.Bfilename)

    def get_pos_file(self): return self.posfile
    def get_B_file(self): return self.Bfilename

class GraceDBJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
    """
    Class for a gracedb job
    """
    def __init__(self,cp,submitFile,logdir,dax=False):
      exe=cp.get('condor','gracedb')
      #pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
      pipeline.CondorDAGJob.__init__(self,"scheduler",exe)
      pipeline.AnalysisJob.__init__(self,cp,dax=dax)
      if cp.has_option('analysis','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
      if cp.has_option('analysis','accounting_group_user'):
        self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
      self.set_sub_file(os.path.abspath(submitFile))
      self.set_stdout_file(os.path.join(logdir,'gracedb-$(cluster)-$(process).out'))
      self.set_stderr_file(os.path.join(logdir,'gracedb-$(cluster)-$(process).err'))
      self.add_condor_cmd('getenv','True')
      self.baseurl=cp.get('paths','baseurl')
      self.basepath=cp.get('paths','webdir')

class GraceDBNode(pipeline.CondorDAGNode):
    """
    Run the gracedb executable to report the results
    """
    def __init__(self,gracedb_job,gid=None,parent=None,message=None,upfile=None,command='upload',tag=None):
        # Message need to be a string
        # Upfile is the full path of the file to be uploaded
        pipeline.CondorDAGNode.__init__(self,gracedb_job)
        if gid: self.set_gid(gid)
        if parent:
          if isinstance(parent, list):
            for p in parent:
              self.add_parent(p)
          else:
            self.add_parent(parent)
        self.message=message
        self.filename=upfile
        self.command=command
        self.tag=tag
        self.__finalized=False

    def set_gid(self,gid):
        """
        Set the GraceDB ID to log to
        """
        self.gid=gid

    def set_message(self,message):
        self.message=message

    def set_filename(self,filename):
        self.filename=filename

    def finalize(self):
        if self.__finalized:
            return
        self.add_var_arg(self.command)
        if self.tag:
          self.add_var_arg('--tag-name='+self.tag)
        self.add_var_arg(str(self.gid))
        if self.filename:
          self.add_var_arg(self.filename+' ')
        if self.message:
          self.add_var_arg(self.message)
        self.__finalized=True

class ROMJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
  """
  Class for a ROM compute weights job
  """
  def __init__(self,cp,submitFile,logdir,dax=False):
    time_step=0.000172895418228
    #This ensures that, for SNR < 100, a signal with bandwidth 20-4096 Hz will
    #have a resolved time posterior assuming a chirp-like frequency evolution
    #and aLIGO_ZDHP PSD
    dt=0.1
    exe=cp.get('condor','computeroqweights')
    pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
    pipeline.AnalysisJob.__init__(self,cp,dax=dax)
    if cp.has_option('analysis','accounting_group'):
      self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    if cp.has_option('analysis','accounting_group_user'):
      self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
    self.set_sub_file(submitFile)
    self.set_stdout_file(os.path.join(logdir,'computeroqweights-$(cluster)-$(process).out'))
    self.set_stderr_file(os.path.join(logdir,'computeroqweights-$(cluster)-$(process).err'))
    self.add_condor_cmd('getenv','True')
    self.add_arg('-B '+str(cp.get('paths','roq_b_matrix_directory')))
    if cp.has_option('engine','dt'):
      dt=cp.get('engine','dt')
    self.add_arg('-t '+str(dt))
    if cp.has_option('engine','time_step'):
      time_step=cp.get('engine','time_step')
    self.add_arg('-T '+str(time_step))
    if cp.has_option('condor','computeroqweights_memory'):
      computeroqweights_memory=str(cp.get('condor','computeroqweights_memory'))
    else:
      params = np.genfromtxt(str(cp.get('paths','roq_b_matrix_directory')+'/params.dat'), names=True)
      computeroqweights_memory=str(int(
      os.path.getsize(str(cp.get('paths','roq_b_matrix_directory')+'/B_linear.npy'))/(1024*1024)
      + ((params['fhigh']-params['flow'])*params['seglen'])*(float(dt)*2/float(time_step))*2*8/(1024*1024)
      + os.path.getsize(str(cp.get('paths','roq_b_matrix_directory')+'/B_quadratic.npy'))/(1024*1024)
      ))
    self.add_condor_cmd('request_memory',computeroqweights_memory)

class ROMNode(pipeline.CondorDAGNode):
  """
  Run the ROM compute weights script
  """
  def __init__(self,computeroqweights_job,ifo,seglen,flow):
    pipeline.CondorDAGNode.__init__(self,computeroqweights_job)
    self.__finalized=False
    self.add_var_arg('-s '+str(seglen))
    self.add_var_arg('-f '+str(flow))
    self.add_var_arg('-i '+ifo)

  def finalize(self):
    if self.__finalized:
      return
    self.__finalized=True

class BayesLineJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
  """
  Class for a BayesLine job
  """
  def __init__(self,cp,submitFile,logdir,dax=False):
    exe=cp.get('condor','bayesline')
    pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
    pipeline.AnalysisJob.__init__(self,cp,dax=dax)
    if cp.has_option('analysis','accounting_group'):
      self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
    if cp.has_option('analysis','accounting_group_user'):
      self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
    self.set_sub_file(submitFile)
    self.set_stdout_file(os.path.join(logdir,'bayesline-$(cluster)-$(process).out'))
    self.set_stderr_file(os.path.join(logdir,'bayesline-$(cluster)-$(process).err'))
    self.add_condor_cmd('getenv','True')

class BayesLineNode(pipeline.CondorDAGNode):
  """
  Run the BayesLine code
  """
  def __init__(self,bayesline_job):
    pipeline.CondorDAGNode.__init__(self,bayesline_job)
    self.__finalized=False

  def finalize(self):
    if self.__finalized:
      return
    self.__finalized=True

class SkyAreaNode(pipeline.CondorDAGNode):
  """
  Node to run sky area code
  """
  def __init__(self,skyarea_job,posfile=None,parent=None,prefix=None):
      pipeline.CondorDAGNode.__init__(self,skyarea_job)
      if parent:
          self.add_parent(parent)
      if posfile:
          self.set_posterior_file(posfile)
      self.ifos=None
      self.outdir=None
      self.prefix=prefix

  def set_posterior_file(self,posfile):
      self.add_file_opt('samples',posfile,file_is_output_file=False)
      self.posfile=posfile
  def set_outdir(self,outdir):
      self.add_var_opt('outdir',outdir)
      self.outdir=outdir
  def get_outdir(self):
      return self.outdir
  def set_fits_name(self):
      name='skymap.fits.gz'
      if self.prefix is not None:
        name=self.prefix+name
      self.add_var_opt('fitsoutname',name)
  def set_injection(self,injfile,eventnum):
      if injfile is not None:
        self.add_file_opt('inj',injfile)
        self.add_var_opt('eventnum',str(eventnum))
  def set_objid(self,objid):
      self.add_var_opt('objid',objid)
  def add_resultspage_parent(self,resultspagenode):
      self.set_posterior_file(resultspagenode.get_pos_file())
      self.set_outdir(resultspagenode.get_output_path())
      self.add_parent(resultspagenode)
      self.set_fits_name()
      self.set_injection(resultspagenode.get_injection(),resultspagenode.get_event_number())
      self.set_objid(resultspagenode.get_event_number())
  def set_ifos(self,ifos=None):
      self.ifos=ifos

class SkyAreaJob(pipeline.CondorDAGJob,pipeline.AnalysisJob):
  """
  Class for Sky Area Jobs
  """
  def __init__(self,cp,submitFile,logdir,dax=False):
      exe=cp.get('condor','skyarea')
      pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
      pipeline.AnalysisJob.__init__(self,cp,dax=dax)
      if cp.has_option('analysis','accounting_group'):
        self.add_condor_cmd('accounting_group',cp.get('analysis','accounting_group'))
      if cp.has_option('analysis','accounting_group_user'):
        self.add_condor_cmd('accounting_group_user',cp.get('analysis','accounting_group_user'))
      self.set_sub_file(submitFile)
      self.set_stdout_file(os.path.join(logdir,'skyarea-$(cluster)-$(process).out'))
      self.set_stderr_file(os.path.join(logdir,'skyarea-$(cluster)-$(process).err'))
      self.add_condor_cmd('getenv','True')
      # Add user-specified options from ini file
      self.add_ini_opts(cp,'skyarea')
