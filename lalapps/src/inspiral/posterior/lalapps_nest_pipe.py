# DAG generation code for running inspnest pipeline
# (C) 2010 John Veitch

from lalapps import nest_utils as OddsPipeline
import ConfigParser
import os
import sys
from optparse import OptionParser,OptionValueError
import glue
from glue import pipeline
from glue import segmentsUtils
from glue import segments
from glue import segmentdb
from lalapps import inspiralutils
import scipy

#Presets
ifos=['H1','L1','V1']
fakeTypes=["LALVirgo","LALLIGO","LALAdLIGO"]
timeslides=None

def open_pipedown_database(database_filename,tmp_space):
    """
    Open the connection to the pipedown database
    """
    if not os.access(database_filename,os.R_OK):
	raise Exception('Unable to open input file: %s'%(database_filename))
    from glue.ligolw import dbtables
    try:
        import sqlite3
    except ImportError:
        # Pre-2.5
        from pysqlite2 import dbapi2 as sqlite3
    working_filename=dbtables.get_connection_filename(database_filename,tmp_path=tmp_space)
    connection = sqlite3.connect(working_filename)
    if tmp_space:
	dbtables.set_temp_store_directory(connection,tmp_space)
    dbtables.DBTable_set_connection(connection)
    return (connection,working_filename) 

def get_timeslides_pipedown(database_connection, dumpfile=None, gpsstart=None, gpsend=None):
	"""
	Returns a dict of dicts {trigtime: {'H1':0,'L1':5,'V1':10},...}
	of times and timeslide offsets
	"""
	output={}
	if gpsstart is not None: gpsstart=float(gpsstart)
	if gpsend is not None: gpsend=float(gpsend)
	if dumpfile is not None:
		outfile=open(dumpfile,'w')
	else:
		outfile=None
	db_segments=[]
	sql_seg_query="SELECT search_summary.out_start_time, search_summary.out_end_time from search_summary join process on process.process_id==search_summary.process_id where process.program=='thinca'"
	db_out = database_connection.cursor().execute(sql_seg_query)
	for d in db_out:
		if d not in db_segments:
			db_segments.append(d)
	seglist=segments.segmentlist([segments.segment(d[0],d[1]) for d in db_segments])
	db_out_saved=[]
	# Get coincidences
	get_coincs="SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9,coinc_event.coinc_event_id \
		    FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name == 'sngl_inspiral' and coinc_event_map.event_id \
		    == sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id==coinc_event_map.coinc_event_id)"
	if gpsstart is not None:
		get_coincs=get_coincs+ ' where sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9 > %f'%(gpsstart)
		joinstr=' and '
	else:
		joinstr=' where '
	if gpsend is not None:
		get_coincs=get_coincs+ joinstr+' sngl_inspiral.end_time+sngl_inspiral.end_time*1e-9 <%f'%(gpsend)
	db_out=database_connection.cursor().execute(get_coincs)
	co_ids=set([d[1] for d in db_out])
	print 'Found %i coinc events in pipedown database'%(len(co_ids))
	for co_id in co_ids:
		coinc_number=co_id.split(":")[-1]
		sql_query="SELECT time_slide.instrument,time_slide.offset from coinc_event join time_slide on \
			   (time_slide.time_slide_id==coinc_event.time_slide_id) where coinc_event.coinc_event_id == '%s' "%(str(co_id))
		db_out=database_connection.cursor().execute(sql_query)
		offsets={}
		for event in db_out:
			offsets[event[0]]=event[1]
		initialkey=offsets.keys()[0]
		sql_query='SELECT sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9 from sngl_inspiral join coinc_event_map on (coinc_event_map.table_name=="sngl_inspiral" and coinc_event_map.event_id == sngl_inspiral.event_id) where coinc_event_map.coinc_event_id == "%s"'%(str(co_id))
		db_out=database_connection.cursor().execute(sql_query)
		db_out_saved=[d for d in db_out]
		if(len(db_out_saved)==0):
			print ' Error! No events found for IFO %s for coinc %s'%(str(initialkey),str(co_id))
			continue
		true_time=db_out_saved[0][0]
		for ifo in offsets.keys():
			deltat=offsets[ifo]
			# The shifted time is also in the same segment as the original time
			segfilt=filter(lambda seg:true_time in seg, seglist)
			if(len(segfilt)!=1):
				print 'ERROR: time %f must appear in only one segment'%(shifted_time)
			else:
				seg=segfilt[0]
			# Original time is shifted_time-deltat
			while (true_time+deltat) not in seg:
				if deltat<0:
					deltat=deltat+(seg[1]-seg[0])
				else:
					deltat=deltat-(seg[1]-seg[0])
			offsets[ifo]=deltat
		realtime=true_time-offsets[initialkey]
		output[realtime]=offsets
		if outfile is not None:
			# Get the SNR and chi-squared to dump
			sql_query = 'SELECT sngl_inspiral.ifo,sngl_inspiral.snr,sngl_inspiral.chisq FROM sngl_inspiral join coinc_event_map on (coinc_event_map.table_name == "sngl_inspiral" and coinc_event_map.event_id==sngl_inspiral.event_id) where coinc_event_map.coinc_event_id == "%s"'%(str(co_id))
			db_out=database_connection.cursor().execute(sql_query)	
			for (ifo,snr,chisq) in db_out:
				print >>outfile,"%s %s %f %f %f %f"%(coinc_number,ifo,realtime-offsets[ifo],snr,chisq,offsets[ifo])
	if outfile is not None:
		outfile.close()
	print 'Found %i time slide coincs to run on'%(len(output))
	return output

def get_timeslides_pipedown_broken(database_connection, dumpfile=None, gpsstart=None, gpsend=None):
	"""
	returns a dict of dicts { trigtime: {'H1':0,'L1':5,'V1':10}) , ... }
	of times and timeslide offsets.
	Optionally dumps times, SNR, chi-squared statistics to a file
	"""
	if gpsstart is not None: gpsstart=float(gpsstart)
        if gpsend is not None: gpsend=float(gpsend)
        if dumpfile is not None:
		outfile=open(dumpfile,'w')
	else:
		outfile=None
	db_segments=[]
	sql_seg_query="SELECT search_summary.out_start_time, search_summary.out_end_time from search_summary join process on process.process_id==search_summary.process_id where process.program=='thinca'"
	db_out = database_connection.cursor().execute(sql_seg_query)
	for d in db_out:
		if d not in db_segments:
			db_segments.append(d)
	seglist=segments.segmentlist([segments.segment(d[0],d[1]) for d in db_segments])
	db_out_saved=[]
	sql_query="SELECT time_slide.time_slide_id,sngl_inspiral.end_time + sngl_inspiral.end_time_ns*1e-9,sngl_inspiral.ifo,time_slide.offset,coinc_event.coinc_event_id,sngl_inspiral.snr,sngl_inspiral.chisq\
		from sngl_inspiral join coinc_event_map on (coinc_event_map.table_name == 'sngl_inspiral' and coinc_event_map.event_id \
		== sngl_inspiral.event_id) join coinc_event on (coinc_event.coinc_event_id == coinc_event_map.coinc_event_id) join time_slide on \
		(time_slide.instrument == sngl_inspiral.ifo and time_slide.time_slide_id == coinc_event.time_slide_id)"
        if gpsstart is not None:
		sql_query=sql_query + ' where sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9 > %f'%(gpsstart)
		if gpsend is not None:
			sql_query=sql_query+' and sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9 < %f'%(gpsend)
	else:
		if gpsend is not None:
			sql_query=sql_query+' where sngl_inspiral.end_time+sngl_inspiral.end_time_ns*1e-9 < %f'%(gpsend)
	db_out = database_connection.cursor().execute(sql_query)
	for d in db_out:
		db_out_saved.append(d)
	output={}
	co_ids=set([tup[4] for tup in db_out_saved])
        print 'Found %i coinc events from %i events in pipedown database'%(len(co_ids),len(db_out_saved))
	for co in co_ids:
		events=[tup for tup in db_out_saved if tup[4]==co]
		d={}
		if outfile is not None:
			for ev in events:
				print >>outfile,'%s %f %f %f %f'%(ev[2],ev[1],ev[3],ev[5],ev[6])
			print >>outfile,'\n'
		# print ' Processing %i events for coinc ID %s'%(len(events),str(co))
		for e in events:
			segfilt=filter(lambda seg:e[1] in seg, seglist)
			if len(segfilt)!=1:
				print ' ERROR: time %f must appear in only one segment'%(e[3])
			else:
				seg=segfilt[0]
			deltat=e[3]
			while (e[1]+deltat) not in seg:
				#print ' adjusting timeslide for '+str(e[1])+' + '+str(deltat)
				if deltat<0:
					deltat=deltat+(seg[1]-seg[0])
				else:
					deltat=deltat-(seg[1]-seg[0])
			d[e[2]]=deltat
		t=events[0][1]+events[0][3]
		if t in output.keys():
			print 'Already found event at time %f, skipping'%(t)
		else:
			output[t]=d
	print 'Found a total of %i useable timeslides'%(len(output))
	if outfile is not None:
		outfile.close()
	return output

# Read command line options
usage=""" %prog [options]
Setup a Condor DAG file to run the nested sampling pipeline

The user must specify either an injection file to analyse, with the --inj option,
a list of SnglInspiralTable or CoincInspiralTable triggers with the --<x>-triggers options,
or a length of time between --gps-start-time and --gps-end-time,
or an ASCII list of GPS times with the --gps-time-file option.
Specifying start and end times in addition to triggers or injections will  restrict
the times analysed to those that fall in the indicated time period.

The user must also specify and ini file which will contain the main analysis setup.

The coherence test option will calculate Bayes factors for a coherent and incoherent
model where possible (>1 IFO available).

"""

parser=OptionParser(usage)
parser.add_option("-i","--ini-file",default=None,action="store",type="string",help="ini-file for nestGRB pipeline",metavar="CONFIG.ini")
parser.add_option("-r","--run-path",default="./",action="store",type="string",help="Directory to run pipeline in (default: %default)",metavar="RUNDIR")
parser.add_option("-p","--dag-log-path",default=None,action="store",type="string",help="Directory to store condor log files for the dag, defaults to RUNDIR/log/ SHOULD BE LOCAL TO SUBMIT NODE",metavar="DAGLOGDIR")
parser.add_option("-l","--jobs-log-path",default=None,action="store",type="string",help="Directory to store stderr and stdout files for the jobs, defaults to RUNDIR/log/ MUST BE VISIBLE FROM THE NODES",metavar="JOBLOGDIR")
parser.add_option("-s","--gps-start-time",action="store",type="string",default=None,help="Start time of analysis")
parser.add_option("-e","--gps-end-time",action="store",type="string",default=None,help="End time of analysis")
parser.add_option("-I","--inj",action="store",type="string",default=None,help="List of injections to perform and analyse",metavar="INJFILE.xml")
parser.add_option("-F","--disable-inject",action="store_true",default=False,help="Don't actually inject the signals - useful for hardware injection lists")
parser.add_option("-t","--single-triggers",action="store",type="string",default=None,help="SnglInspiralTable trigger list",metavar="SNGL_FILE.xml")
parser.add_option("-T","--coinc-triggers",action="store",type="string",default=None,help="CoincInspiralTable trigger list",metavar="COINC_FILE.xml")
parser.add_option("-g","--gps-time-file",action="store",type="string",default=None,help="Text file containing list of GPS times to analyse",metavar="TIMES.txt")
parser.add_option("-P","--pipedown-db",action="store",type="string",default=None,help="Pipedown database of triggers",metavar="PIPEDOWN-DB.sqlite")
parser.add_option("-C","--coherence-test",action="store_true",help="Perform the coherence test",default=False)
parser.add_option("--disable-pages",action="store_true",default=False,help="Disable the outpage page and postprocessing")
parser.add_option("--enable-pages-alltimes",action="store_true",default=False,help="Enable a results page for every time bin. May create a lot of dirs.")
parser.add_option("--condor-submit",action="store_true",default=False,help="Automatically submit the condor dag")
parser.add_option("--ignore-science-mode",action="store_true",default=False,help="Skip the segment database step and try and find data for all specified times.")
parser.add_option("--notify",action="store",default=None,help="Optional e-mail address to notify upon completion",metavar="you@yourdomain.edu")
parser.add_option("-X","--program",action="store",default="inspnest",help="Over-ride program to run, e.g. -X lalinferencenest",metavar="lalinferencenest")
parser.add_option("--timeslides",action="store_true",default=False,help="Analyse timeslides (requires --pipedown-db)")

(opts,args)=parser.parse_args()

cp=ConfigParser.ConfigParser()
cp.optionxform = str
cp.readfp(open(opts.ini_file))
#Get ifos from string list
ifos=eval(cp.get('analysis','ifos'))

if opts.program == "inspnest":
    nestclass=OddsPipeline.InspNestNode
if opts.program == "lalinferencenest":
    nestclass=OddsPipeline.LALInferenceNode

## This function is used to put the ifos.keys() in the same order than in the ifos section in the parser, whist passing them to the inspnest nodes
def sort_ifo(a,b):
    global ifos
    if ifos.index(a)<ifos.index(b):
        return -1
    elif ifos.index(a)>ifos.index(b):
        return 1
    else: 
        return 0

print 'Setting up pipeline with ifos '+' '.join(map(str,ifos))
IFOs=''.join(map(str,ifos))

def checkDir(dir):
        if not os.access(dir,os.F_OK):
                try:
                        print 'Creating %s\n'%(dir)
                        os.makedirs(dir)
                except:
                        print 'Error: Unable to create directory %s'%(dir)
                        sys.exit(1)
        if not os.access(dir,os.W_OK):
                print 'Unable to write to %s'%(dir)
                sys.exit(1)
        return True

# Check command line options
if opts.ini_file is None:
        parser.error("Must specify an ini file using --ini-file")
#if opts.trigger_time is None:
#       parser.error("Must specify the trigger time using --trigger-time")

checkDir(opts.run_path)
if opts.dag_log_path is None:
        opts.dag_log_path=os.path.join(opts.run_path,"log")
checkDir(opts.dag_log_path)
if opts.jobs_log_path is None:
        opts.jobs_log_path=os.path.join(opts.run_path,"log")
checkDir(opts.jobs_log_path)

# Get science segments

# Read the ini file using ConfigParser

cp = ConfigParser.ConfigParser()
cp.optionxform = str
cp.readfp(open(opts.ini_file))
dt=float(cp.get('inspnest','dt'))

# Simple function because lambda form doesn't work in python 2.4
def istoint(a):
    if a is not None:
        return 1
    else:
        return 0

if(sum(map(istoint, [getattr(opts,a) for a in ['single_triggers','coinc_triggers','inj','gps_time_file','pipedown_db']]))>1):
    parser.error('You can only specify one of --inj, --single-triggers, --coinc-triggers, --pipedown-db or --gps-time-file\n')

types_list=eval(cp.get('data','types'))
channels_list=eval(cp.get('data','channels'))


#Create types and channels dictionaries
types={}
channels={}
for ifo,type,channel in zip(ifos,types_list,channels_list):
    types[ifo]=type
    channels[ifo]=channel

for ifo,channel in channels.items():
	cp.set('data',ifo.lower()+'-channel',channel)

if cp.has_option('analysis','time-slide-dump'):
	timeslidefiledump= cp.get('analysis','time-slide-dump')
else:
	timeslidefiledump=None

time_event=None
if opts.inj and cp.has_option('analysis','events'):
    time_event={}
    events=[]
    times=[]
    raw_events=cp.get('analysis','events').replace('[','').replace(']','').split(',')

    from pylal import SimInspiralUtils
    injTable=SimInspiralUtils.ReadSimInspiralFromFiles([opts.inj])

    if 'all' is raw_events or 'all' in raw_events:
        events=range(len(injTable))

    else:
        for raw_event in raw_events:
            if ':' in raw_event:
                limits=raw_event.split(':')
                if len(limits) != 2:
                    print "Error: in event config option; ':' must separate two numbers."
                    exit(0)
                low=int(limits[0])
                high=int(limits[1])
                if low>high:
                    events.extend(range(int(high),int(low)))
                elif high>low:
                    events.extend(range(int(low),int(high)))
            else:
                events.append(int(raw_event))

    for event in events:
        times.append(injTable[event].get_end())
        time_event[times[-1]]=event
        

    starttime=min(times)-1.0
    endtime=max(times)+1.0

    cp.set('results','injXML',opts.inj)

elif opts.inj:
    from pylal import SimInspiralUtils
    injTable=SimInspiralUtils.ReadSimInspiralFromFiles([opts.inj])
    times=[inj.get_end() for inj in injTable]
    starttime=min(times)
    endtime=max(times)

if opts.gps_time_file:
    import re
    p=re.compile('[\d.]+')
    times=[]
    timefile=open(opts.gps_time_file,'r')
    for time in timefile:
        if not p.match(time):
            continue
        if float(time) in times:
            print 'Skipping duplicate time %s'%(time)
            continue
        print 'Read time %s'%(time)
        times.append(float(time))
    timefile.close()
    starttime=min(times)
    endtime=max(times)

if opts.single_triggers:
    from pylal import SnglInspiralUtils
    trigTable=SnglInspiralUtils.ReadSnglInspiralFromFiles([opts.single_triggers])
    times=[trig.get_end() for trig in trigTable]
    starttime=min(times)
    endtime=max(times)

if opts.coinc_triggers:
    from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
    from glue.ligolw import table
    from glue.ligolw import lsctables
    from glue.ligolw import utils

    CoincInspiralTriggers=None
    doc = utils.load_filename(opts.coinc_triggers, gz=(opts.coinc_triggers or "stdin").endswith(".gz"), verbose=False)
    # extract the sim inspiral table
    try: CoincInspiralTable = \
        table.get_table(doc, lsctables.CoincInspiralTable.tableName)
    except: CoincInspiralTable = None
    if CoincInspiralTriggers and CoincInspiralTable: 
        CoincInspiralTriggers.extend(CoincInspiralTable)
    elif not CoincInspiralTriggers:
        CoincInspiralTriggers = CoincInspiralTable

    times=[trig.get_end() for trig in CoincInspiralTriggers]
    starttime=min(times)
    endtime=max(times)

# Read pipedown database
if opts.pipedown_db:
	print 'establishing connection to pipedown database'
	(db_connection,db_working_filename)=open_pipedown_database(opts.pipedown_db,None)[0]
	timeslides=get_timeslides_pipedown(db_connection,dumpfile=timeslidefiledump,gpsstart=opts.gps_start_time,gpsend=opts.gps_end_time)
	times=timeslides.keys()
	starttime=min(times)
	endtime=max(times)

if opts.gps_start_time is not None:
    mystarttime=float(opts.gps_start_time)
    if opts.inj:
        starttime=max([starttime,mystarttime])
        times=[time for time in times if time>starttime]
    else:
        starttime=mystarttime
    cp.set('analysis','gps-start-time',str(starttime))


if opts.gps_end_time is not None:
    myendtime=float(opts.gps_end_time)
    if opts.inj:
        endtime=min([endtime,myendtime])
        times=[time for time in times if time<endtime]
    else:
        endtime=myendtime
    cp.set('analysis','gps-end-time',str(endtime))

if not opts.inj and not opts.coinc_triggers and not opts.gps_time_file and not opts.single_triggers and not opts.pipedown_db:
    if opts.gps_start_time is None and opts.gps_end_time is None:

        raise OptionValueError('Must specify both a start and end time if not using injection file or coinc file.\n')


    times=scipy.linspace(float(starttime)+dt/2.0,float(endtime)-dt/2.0,(float(endtime)-float(starttime))/dt)

padding=float(cp.get('analysis','padding'))
length=float(endtime)-float(starttime)+2*padding

datastart=starttime-padding
dataend=endtime+padding

# Check for manual over-ride of PSD start time
if cp.has_option('lalinference','psdstart'):
	manual_datastart=cp.getfloat('lalinference','psdstart')-padding
	if manual_datastart<datastart:
		datastart=manual_datastart

if cp.has_option('lalinference','psdlength'):
	manual_length=cp.getfloat('lalinference','psdlength')
	if dataend<datastart+manual_length+2.0*padding:
		dataend=datastart+manual_length+2.0*padding

cp.add_section('input')
cp.set('input','gps-start-time',str(int(datastart)))
cp.set('input','gps-end-time',str(int(dataend)))


print 'Setting up analysis for %i times between %.3f and %.3f\n'%(len(times),starttime,endtime)

cache_dir=os.path.join(opts.run_path,'data')
checkDir(cache_dir)


segdir=os.path.join(opts.run_path,'segments')
checkDir(segdir)
os.chdir(segdir)
science_segs={}
seg_files={}
segs={}

basename="nest_%.3f-%.3f"%(starttime,endtime)
# Create DAG #################
daglogfile=os.path.join(opts.dag_log_path,basename+'.log')
dagfile=os.path.join(opts.run_path,basename)

dag = pipeline.CondorDAG(daglogfile)
dag.set_dag_file(dagfile)

datafind_job = pipeline.LSCDataFindJob(cache_dir,opts.jobs_log_path,cp)
datafind_job.add_opt('url-type','file')
datafind_job.set_sub_file(os.path.join(opts.run_path,'datafind.sub'))

# Build list of science segments
# Covers entire time range. Not necessarily all used
for ifo in ifos:
    if not opts.ignore_science_mode and types[ifo] not in fakeTypes:
        seg_files[ifo]=inspiralutils.science_segments(ifo,cp)
        segfile=open(seg_files[ifo])
        #segs[ifo]=segmentsUtils.fromfilenames([seg_files[ifo]])
        segs[ifo]=segmentsUtils.fromsegwizard(segfile)
        segs[ifo].coalesce()
        segfile.close()
    else:   # If we skip the segdb step, just construct a large segment
        print 'Faking segment from %i to %i\n'%(datastart,dataend)
        segs[ifo]=segments.segmentlist([segments.segment(int(datastart),int(dataend))])


for ifo in ifos:
    science_segs[ifo]=[]
    if types[ifo] in fakeTypes:
        science_segs[ifo].append(None)
    else:
        # Setup find data jobs
        for seg in segs[ifo]:
            sciseg=pipeline.ScienceSegment((segs[ifo].index(seg),seg[0],seg[1],seg[1]-seg[0]))
            science_segs[ifo].append(sciseg)
            df_node=pipeline.LSCDataFindNode(datafind_job)
            df_node.set_start(int(sciseg.start()))
            df_node.set_end(int(sciseg.end()))
            df_node.set_observatory(ifo[0])
            df_node.set_type(types[ifo])
            sciseg.set_df_node(df_node)

os.chdir('../../')

# Now loop over times and add datafind nodes to the dag

filtered_time=filter(lambda t: reduce(lambda a,b:a or b, map(lambda ifo: t in segs[ifo],ifos)), times)
times=filtered_time
print 'Found segments for %i times\n'%(len(times))

df_nodes_by_time={}
ifos_by_time={}
segments_by_time={}
for time in times:
    ifos_by_time[time]=[]
    segments_by_time[time]={}
    for ifo in ifos:
        for (seg,sciseg) in zip(segs[ifo],science_segs[ifo]):
            if time in seg:
                segments_by_time[time][ifo]=(seg,sciseg)
            #else:
            #   segments_by_time[time][ifo]=None
        #if segments_by_time[time][ifo] is None:
        #   continue
                ifos_by_time[time].append(ifo)
                if segments_by_time[time][ifo][1] is not None:
                    node=segments_by_time[time][ifo][1].get_df_node()
                    if node not in dag.get_nodes():
                        dag.add_node(node)

# Filter out times with no data
#times=[time for time in times if len(segments_by_time[time].keys())>0 ]
print 'Found data for %i times\n'%(len(times))

if len(times)==0:
    print 'Unable to find any science data for times specified. Aborting!'
    sys.exit(1)

inspnest_subfile=os.path.join(opts.run_path,'inspnest.sub')
if opts.program == 'inspnest':
	inspnest_job=OddsPipeline.InspNestJob(cp,inspnest_subfile,opts.jobs_log_path)
else:
	if opts.program == 'lalinferencenest':
		inspnest_job=OddsPipeline.LALInferenceJob(cp,inspnest_subfile,opts.jobs_log_path)
	else:
		raise Exception('Unknown program %s'%(opts.program))

inspnest_job.add_opt('Nlive',cp.get('analysis','nlive'))
inspnest_job.add_opt('Nmcmc',cp.get('analysis','nmcmc'))

if opts.inj and not opts.disable_inject:
    inspnest_job.add_opt('inj',opts.inj)

nest_path=os.path.join(opts.run_path,'nest')
checkDir(nest_path)

combine_dir=os.path.join(opts.run_path,'combine')

# Make jobs write SNR files to the appropriate place
snr_dir=os.path.join(opts.run_path,'SNR')
checkDir(snr_dir)
inspnest_job.add_opt('snrpath',snr_dir)

time_combine_sub_file=os.path.join(opts.run_path,'combine_each_time.sub')
combine_each_time_job=OddsPipeline.CombineZJob(cp,time_combine_sub_file,opts.jobs_log_path)
combine_each_time_nodes={}

inspnest_nodes={}
inspnest_nodes_by_time={}
nest_files={}
# Now loop over times and set up cpmbine jobs for getting individual posteriors
nparallel=int(cp.get('analysis','nparallel'))
posfiles_by_time={}
bayesfiles_by_time={}
snrs_by_time={}
for time in times:
    posfile=os.path.join(combine_dir,'posterior_samples_%.3f'%(time))
    posfiles_by_time[time]=posfile
    combine_node=OddsPipeline.CombineZNode(combine_each_time_job)
    combine_node.add_file_opt('outpos',posfile,file_is_output_file=True)
    bayesfiles_by_time[time]=os.path.join(combine_dir,'bayesfactor_%.3f.txt'%(time))
    snrs_by_time[time]=os.path.join(snr_dir,'snr_%s_%9.1f.dat'%(IFOs,time))
    combine_node.add_file_opt('bayesfactor',bayesfiles_by_time[time])
    combine_each_time_nodes[time]=combine_node
# Set up inspnest jobs, using parallel if needed
if nparallel==1:
    for time in times:
        if time_event is not None and not opts.disable_inject:
                thisevent=time_event[time]
        else:
                thisevent=None
	if timeslides:
		if time in timeslides.keys():
			ts=timeslides[time]
	else:
		ts=None
        inspnest_node=OddsPipeline.setup_single_nest(cp,inspnest_job,time,segments_by_time[time],nest_path,ifos=sorted(segments_by_time[time].keys(),sort_ifo),event=thisevent,nodeclass=nestclass,timeslides=ts)
        inspnest_nodes_by_time[time]=inspnest_node
        dag.add_node(inspnest_node)
        combine_each_time_nodes[time].add_parent(inspnest_node)
        nest_files[time]=inspnest_node.get_output_files()[0]
else:
    merge_sub_file=os.path.join(opts.run_path,'merge.sub')
    merge_job=OddsPipeline.MergeJob(cp,merge_sub_file,opts.jobs_log_path)
    for time in times:
        if time_event is not None and not opts.disable_inject:
                thisevent=time_event[time]
        else:
                thisevent=None
	if timeslides:
		if time in timeslides.keys():
			ts=timeslides[time]
	else:
		ts=None
        (merge_node,inspnest_nodes)=OddsPipeline.setup_parallel_nest(cp,inspnest_job,merge_job,time,segments_by_time[time],nest_path,ifos=sorted(segments_by_time[time].keys(),sort_ifo),event=thisevent,nodeclass=nestclass,timeslides=ts)
        dag.add_node(merge_node)
        map(dag.add_node,inspnest_nodes)
        nest_files[time]=merge_node.get_output_files()[0]
        combine_each_time_nodes[time].add_parent(merge_node)
        inspnest_nodes_by_time[time]=merge_node

for time in times:
    dag.add_node(combine_each_time_nodes[time])
    combine_each_time_nodes[time].add_var_arg(nest_files[time])
    combine_each_time_nodes[time].add_file_opt('headers',nest_files[time]+'_params.txt')

# Set up coherence test jobs if required
if opts.coherence_test:
    coherence_test_sub_file=os.path.join(opts.run_path,'coherence_test.sub')
    coherence_test_job=OddsPipeline.CoherenceTestJob(cp,coherence_test_sub_file,opts.jobs_log_path)
    coherence_inspnest_nodes={}
    coherence_test_posfiles={}
    coherence_test_nodes={}
    coherence_test_files={}
    coherence_test_outfiles={}
    coherence_test_bfiles={}
    combine_nodes_by_time={}
    snr_files_by_time={}
    coherence_test_dir=os.path.join(opts.run_path,'coherencetest')
    checkDir(coherence_test_dir)
    # Loop over all the times
    for time in times:
	if timeslides:
		if time in timeslides.keys():
			ts=timeslides[time]
	else:
		ts=None
        ifos=segments_by_time[time].keys()
        if len(ifos)<2:
            print 'Unable to setup coherence test for time %.3f, only %s data available.\n'%(time,ifos[0])
            continue
        coherence_test_node=OddsPipeline.CoherenceTestNode(coherence_test_job)
        coherence_test_node.add_var_arg(nest_files[time]+'_B.txt')
        coherence_test_files[time]={}
        coherence_test_posfiles[time]={}
        coherence_test_bfiles[time]={}
        coherence_inspnest_nodes[time]={}
        combine_nodes_by_time[time]={}
        snr_files_by_time[time]={}
        for ifo in segments_by_time[time].keys():
            combine_node=OddsPipeline.CombineZNode(combine_each_time_job)
            combine_nodes_by_time[time][ifo]=combine_node
            coherence_test_posfiles[time][ifo]=os.path.join(combine_dir,'posfile_%.3f_%s.dat'%(time,ifo))
            coherence_test_bfiles[time][ifo]=os.path.join(combine_dir,'bayesfactor_%.3f_%s.txt'%(time,ifo))
            snr_files_by_time[time][ifo]=os.path.join(snr_dir,'snr_%s_%9.1f.dat'%(ifo,time))
            combine_node.add_file_opt('bayesfactor',coherence_test_bfiles[time][ifo],file_is_output_file=True),
            combine_node.add_file_opt('outpos',coherence_test_posfiles[time][ifo],file_is_output_file=True)
            if time_event is not None and not opts.disable_inject:
                thisevent=time_event[time]
            else:
                thisevent=None
            if nparallel>1:
                (merge_node,inspnest_nodes)=OddsPipeline.setup_parallel_nest(cp,inspnest_job,merge_job,time,segments_by_time[time],nest_path,ifos=[ifo],event=thisevent,nodeclass=nestclass,timeslides=ts)
                dag.add_node(merge_node)
                coherence_test_files[time][ifo]=merge_node.get_output_files()[0]+'_B.txt'
                combine_node.add_parent(merge_node)
                coherence_test_node.add_parent(merge_node)
                coherence_inspnest_nodes[time][ifo]=merge_node
            else:
                inspnest_node=OddsPipeline.setup_single_nest(cp,inspnest_job,time,segments_by_time[time],nest_path,ifos=[ifo],event=thisevent,nodeclass=nestclass,timeslides=ts)
                inspnest_nodes=[inspnest_node]
                coherence_test_files[time][ifo]=inspnest_node.get_output_files()[0]+'_B.txt'
                combine_node.add_parent(inspnest_node)
                coherence_inspnest_nodes[time][ifo]=inspnest_node
            combine_node.add_var_arg(coherence_test_files[time][ifo])
            combine_node.add_file_opt('headers',coherence_test_files[time][ifo]+'_params.txt')
            coherence_test_node.add_var_arg(coherence_test_files[time][ifo])
            map(coherence_test_node.add_parent, inspnest_nodes)
            map(dag.add_node,inspnest_nodes)
            dag.add_node(combine_node)
        coherence_test_outfiles[time]=os.path.join(coherence_test_dir,'coherence_test_%.3f.txt'%(time))
        coherence_test_node.add_file_opt('outfile',coherence_test_outfiles[time],file_is_output_file=True)
        coherence_test_node.add_parent(inspnest_nodes_by_time[time])
        coherence_test_nodes[time]=coherence_test_node
        dag.add_node(coherence_test_node)


times_by_seg={}
seg_block={}
inspnest_by_seg={}
nestfiles_by_seg={}
co_nests_by_seg={}
co_nestfiles_by_seg={}
# Also want to combine all jobs within each segment
for time in times:
    ifos=segments_by_time[time].keys()
    smallest_segment=segments.segment(max([segments_by_time[time][ifo][0][0] for ifo in ifos]), min([segments_by_time[time][ifo][0][1] for ifo in ifos]))
#   smallest_segment=reduce(lambda a,b:a and b, [segments_by_time[time][ifo][0] for ifo in segments_by_time[time].keys()])
    if smallest_segment not in times_by_seg.keys():
        times_by_seg[smallest_segment]=[]
        inspnest_by_seg[smallest_segment]=[]
        nestfiles_by_seg[smallest_segment]=[]
        co_nests_by_seg[smallest_segment]={}
        co_nestfiles_by_seg[smallest_segment]={}
        if len(segments_by_time[time].keys())>1:
            for ifo in segments_by_time[time].keys():
                co_nests_by_seg[smallest_segment][ifo]=[]
                co_nestfiles_by_seg[smallest_segment][ifo]=[]
    times_by_seg[smallest_segment].append(time)
    inspnest_by_seg[smallest_segment].append(inspnest_nodes_by_time[time])
    nestfiles_by_seg[smallest_segment].append(nest_files[time])
    # Append the coherence test files
    if opts.coherence_test and len(segments_by_time[time].keys())>1:
        for ifo in segments_by_time[time].keys():
            co_nests_by_seg[smallest_segment][ifo].append(coherence_inspnest_nodes[time][ifo])
            co_nestfiles_by_seg[smallest_segment][ifo].append(coherence_test_files[time][ifo])

print 'Broke into %i chunks for recombination'%(len(times_by_seg))

combine_sub_file=os.path.join(opts.run_path,'combine.sub')
combine_job=OddsPipeline.CombineZJob(cp,combine_sub_file,opts.jobs_log_path)
combine_nodes={}
combine_dir=os.path.join(opts.run_path,'combine')
checkDir(combine_dir)

if opts.disable_pages is False:
    # Set up web page producing nodes
    results_page_sub_file=os.path.join(opts.run_path,'results_page.sub')
    results_page_job=OddsPipeline.ResultsPageJob(cp,results_page_sub_file,opts.jobs_log_path)
    if(len(times)==1):
        results_path=os.path.join(cp.get('results','basedir'),str(time))
    else:
        results_path=os.path.join(cp.get('results','basedir'),'%.3f-%.3f'%(min(times),max(times)))
    if(opts.notify is not None):
        results_page_job.set_notification('Complete')
        results_page_job.add_condor_cmd('notify_user',opts.notify)

# Analyse each segment

for seg in times_by_seg.keys():
    combine_node=OddsPipeline.CombineZNode(combine_job)
    posfile=os.path.join(combine_dir,'posterior_samples_%.3f-%.3f.dat'%(seg[0],seg[1]))
    bffile=os.path.join(combine_dir,'bayesfactor_%.3f-%.3f.txt'%(seg[0],seg[1]))
    combine_node.add_file_opt('outpos',posfile)
    combine_node.add_file_opt('bayesfactor',bffile)
    outsamp=os.path.join(combine_dir,'nested_samples_%.3f-%.3f.dat'%(seg[0],seg[1]))
    combine_node.add_file_opt('outsamp',outsamp)
    map(combine_node.add_var_arg,nestfiles_by_seg[seg])
    map(lambda a: combine_node.add_file_opt('headers',a+'_params.txt'),nestfiles_by_seg[seg])
    map(combine_node.add_parent,inspnest_by_seg[seg])
    dag.add_node(combine_node)
    combine_nodes[seg]=combine_node
    
    if opts.coherence_test and len(co_nests_by_seg[seg].keys())>1:
        co_combine_node={}
        co_pos_file={}
        co_bf_file={}
        co_outsamp={}
        for ifo in co_nests_by_seg[seg].keys():
            co_combine_node[ifo]=OddsPipeline.CombineZNode(combine_job)
            co_pos_file[ifo]=posfile.replace('.dat','_'+ifo+'.dat')
            co_bf_file[ifo]=bffile.replace('.txt','_'+ifo+'.txt')
            co_outsamp[ifo]=outsamp.replace('.dat','_'+ifo+'.dat')
            co_combine_node[ifo].add_file_opt('outpos',co_pos_file[ifo],file_is_output_file=True)
            co_combine_node[ifo].add_file_opt('bayesfactor',co_bf_file[ifo],file_is_output_file=True)
            co_combine_node[ifo].add_file_opt('outsamp',co_outsamp[ifo],file_is_output_file=True)
            map(co_combine_node[ifo].add_parent,co_nests_by_seg[seg][ifo])
            map(lambda a:co_combine_node[ifo].add_file_opt('headers',a+'_params.txt'),co_nestfiles_by_seg[seg][ifo])
            map(co_combine_node[ifo].add_var_arg,co_nestfiles_by_seg[seg][ifo])
            dag.add_node(co_combine_node[ifo])
        coherence_test_node=OddsPipeline.CoherenceTestNode(coherence_test_job)
        cotest_outfile_seg=os.path.join(coherence_test_dir,'coherence_test_%.3f-%.3f.txt'%(seg[0],seg[1]))
        coherence_test_node.add_var_opt('outfile',cotest_outfile_seg)
        coherence_test_node.add_var_arg(outsamp)
        map(coherence_test_node.add_var_arg, co_outsamp.values())
        map(coherence_test_node.add_parent, co_combine_node.values())
        coherence_test_node.add_parent(combine_node)
        dag.add_node(coherence_test_node)
        # Generate results page for each IFO too
        if opts.disable_pages is False:
            for ifo in co_nests_by_seg[seg].keys():
                results_node=OddsPipeline.ResultsPageNode(results_page_job)
                results_node.add_file_opt('data',co_pos_file[ifo])
                results_node.add_parent(co_combine_node[ifo])
                seg_results_path=os.path.join(results_path,'%.3f-%.3f'%(seg[0],seg[1]),ifo)
                results_node.add_var_opt('outpath',seg_results_path)
                results_node.add_var_opt('bsn',co_bf_file[ifo])
                dag.add_node(results_node)


    if opts.disable_pages is False:
        # Produce web page output for each segment
        results_node=OddsPipeline.ResultsPageNode(results_page_job)
        results_node.add_file_opt('data',posfile)
        seg_results_path=os.path.join(results_path,'%.3f-%.3f'%(seg[0],seg[1]))
        checkDir(seg_results_path)
        results_node.add_var_opt('outpath',seg_results_path)
        results_node.add_var_opt('bsn',bffile)
        if opts.coherence_test and len(co_nests_by_seg[seg].keys())>1:
            results_node.add_parent(coherence_test_node)
            results_node.add_var_arg('--bci '+cotest_outfile_seg)
        results_node.add_parent(combine_node)
        dag.add_node(results_node)

if opts.enable_pages_alltimes is True:
    # Also produce plots for each time chunk
    for time in times:
        results_page_node=OddsPipeline.ResultsPageNode(results_page_job)
        results_page_node.add_file_opt('data',posfiles_by_time[time])
        results_dir=os.path.join(results_path,'timebins',str(time))
        checkDir(results_dir)
        results_page_node.add_var_opt('outpath',results_dir)
        if time_event is not None:   
            if time_event[time] is not None:
                results_page_node.set_event_number(time_event[time])
                results_page_node.add_var_arg('--inj '+os.path.abspath(opts.inj))
        results_page_node.add_parent(combine_each_time_nodes[time])
        results_page_node.add_var_opt('bsn',bayesfiles_by_time[time])
        results_page_node.add_var_arg('--snr '+snrs_by_time[time])
	dag.add_node(results_page_node)

        if opts.coherence_test and (time in coherence_inspnest_nodes.keys()):
            results_page_node.add_var_arg('--bci '+coherence_test_outfiles[time])
            results_page_node.add_parent(coherence_test_nodes[time])
            if opts.coherence_test and (time in coherence_inspnest_nodes.keys()):
                for ifo in coherence_inspnest_nodes[time].keys():
                    results_node=OddsPipeline.ResultsPageNode(results_page_job)
                    results_node.add_file_opt('data',coherence_test_posfiles[time][ifo])
                    if time_event is not None:   
                        if time_event[time] is not None:
                            results_node.set_event_number(time_event[time])
                            results_node.add_var_arg('--inj '+os.path.abspath(opts.inj))
                    results_node.add_parent(combine_nodes_by_time[time][ifo])
                    ifo_results_path=os.path.join(results_dir,ifo)
                    results_node.add_var_opt('outpath',ifo_results_path)
                    results_node.add_var_opt('bsn',coherence_test_bfiles[time][ifo])
                    results_node.add_var_arg('--snr '+snr_files_by_time[time][ifo])
                
                    dag.add_node(results_node)


# Write the DAG file
dag.write_sub_files()
dag.write_dag()
dag.write_script()

# End of program
print 'Successfully created DAG file.'
print 'Now run condor_submit_dag %s\n'%(dag.get_dag_file())

if opts.condor_submit:
    import subprocess
    from subprocess import Popen
    
    x = subprocess.Popen(['condor_submit_dag',dag.get_dag_file()])
    x.wait()
    if x.returncode==0:
        print 'Submitted DAG file'
    else:
        print 'Unable to submit DAG file'

