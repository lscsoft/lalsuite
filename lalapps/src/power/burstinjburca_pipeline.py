#!/usr/bin/env python2.2
"""
burst_pipeline.py - standalone burst pipeline driver script
This has been created by modifying Duncan Brown's original script for running
the standalone inspiral code

$Id$

This script produced the necessary condor submit and dag files to run
the standalone burst code on LIGO data
"""

__author__ = 'Saikat Ray Majumder <saikat@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import os
import sys
import string
import getopt
import re
import tempfile
import ConfigParser

def isplay(t):
  if ((t - 729273613) % 6370) < 600:
    return 1
  else:
    return 0

class BurstChunk:
  def __init__(self,start,end):
    self.start = start
    self.end = end
    self.length = end - start

class ScienceSegment:
  def __init__(self,id,start,end,duration,segpad):
    self.id = id
    self.start = start
    self.end = end
    self.duration = duration
    self.used = 0
    self.chunks = []
    self.segpad = segpad
    self.startpad = self.start - self.segpad
    self.endpad = self.end + self.segpad

  def createchunks(self,length,overlap,play_only):
    dur = self.duration
    start = self.start
    seg_incr = length - overlap
    while dur >= 0:
      if dur >= length:
        middle = start + length/2
        end = start + length
      else:
          middle = start + dur/2
          end = start + dur
      if (not play_only) or ( play_only 
        and ( isplay(start) or isplay(middle) or isplay(end) ) ):
        self.chunks.append(BurstChunk(start,end))
      start += seg_incr
      dur -= seg_incr
    self.used = self.duration - dur


class BurstPipeline:
  """
  Contains a dictionary of science segments and chunks that
  define the data to analyze and how it should be analyzed
  """
  def __init__(self,config_file):
    """
    configfile = the name of the file containing configuration parameters
    """
    # parse the configuration file into a dictionary
    self.config = {}
    cp = ConfigParser.ConfigParser()
    cp.read(config_file)
    for sec in cp.sections():
      name = string.lower(sec)
      self.config[name] = {}
      for opt in cp.options(sec):
        self.config[name][opt] = string.strip(cp.get(sec,opt))
    self.basename = config_file.split('.')[0]
    tempfile.template = self.basename
    self.logfile = tempfile.mktemp()

  def parsesegs(self):
    self.segments = []
    segpad = int(self.config['burst']['pad-data'])
    # lines that start with an octothorpe are comments
    comment_line = re.compile(r'\A#')
    for line in open(self.config['input']['segments']):
      if not comment_line.match(line):
        segpars = tuple(map(int,line.split()))
        # change this line if the format of the science segment file changes
        self.segments.append( 
          ScienceSegment(segpars[0],segpars[1],segpars[2],segpars[3],segpad) )
  
  def createchunks(self,play_only):
    # compute the chunk and overlap length in seconds
    chunk_length = int(self.config['datacond']['chunk-length'])
    overlap = int(self.config['datacond']['overlap'])
    for seg in self.segments:
      seg.createchunks(chunk_length,overlap,play_only)

  def status(self,play_only):
    log_fh = open( self.basename + '.pipeline.log', 'w' )
    print >> log_fh, """$Id$"""
    total_science_used = 0
    total_science = 0
    num_chunks = 0
    for sciseg in self.segments:
      total_science += sciseg.duration
      total_science_used += sciseg.used
      num_chunks += len(sciseg.chunks)
    print >> log_fh, """\
read %d science segments
got %d seconds of science data of which %d seconds is usable
created %d burst chunks for analysis\
""" % ( len(self.segments), total_science, total_science_used, num_chunks )
    if play_only:
      print >> log_fh, "filtering only segments that overlap with playground"
    print >> log_fh, "condor log file is", self.logfile
    log_fh.close()

  def frcachesub(self):
    sub_fh = open( self.basename + '.frcache.sub', 'w' )
    print >> sub_fh, """\
universe = scheduler
executable = %s
arguments = --lal-cache \\
  --instrument $(site) --type %s \\
  --start $(frstart) --end $(frend)
environment = LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH);PYTHONPATH=$ENV(PYTHONPATH)
log = %s
error = datafind/frcache-$(site)-$(frstart)-$(frend).$(cluster).$(process).err
output = cache/frcache-$(site)-$(frstart)-$(frend).out
notification = never
queue
""" % (self.config['condor']['datafind'],
       self.config['input']['datatype'],
       self.logfile)
    sub_fh.close()

  def binjsub(self):
    sub_fh = open( self.basename + '.binj.sub', 'w' )
    print >> sub_fh, """\
universe = standard
executable = %s
arguments = --gps-start-time $(frstart) --gps-end-time $(frend) --user-tag $(usertag) \\
""" % (self.config['condor']['binj']),
    for sec in ['binjarg']:
      for arg in self.config[sec].keys():
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """    
log = %s
error = binj/binj-$(frstart)-$(frend).$(cluster).$(process).err
output = binj/binj-$(frstart)-$(frend).$(cluster).$(process).out
notification = never
queue
"""%  (self.logfile)
    

  def burstsub(self):

    boolargs = re.compile(r'(printSpectrum|verbose)')
    sub_fh = open( self.basename + '.burst.sub', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments =  --nseg $(nseg) \\
  --npts $(npts) \\
  --numpts $(numpts) \\
  --channel $(channel) \\
  --threshold $(alphath) \\
  --start_time $(start) \\ 
  --framecache cache/frcache-$(site)-$(frstart)-$(frend).out \\
  --calcache /ldas_outgoing/calibration/cache_files/$(ifo)-CAL-V03-729273600-734367600.cache \\
  --injfile HL-INJECTIONS_1_$(usertag)-$(frstart)-$(frduration).xml \\
  --comment $(usertag) \\
  --cluster \\
 """ % (self.config['condor']['universe'],self.config['condor']['burst']),
    for sec in ['arg']:
      for arg in self.config[sec].keys():
        if boolargs.match(arg):
          if re.match(self.config[sec][arg],'true'): 
            print >> sub_fh, "--" + arg,
        else:
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """
log = %s
error = burst/burst-$(ifo)-$(start)-$(end).$(cluster).$(process).err
output = burst/burst-$(ifo)-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.logfile)

    boolargs = re.compile(r'(printSpectrum|verbose)')
    sub_fh = open( self.basename + '.burstH2.sub', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments =  --nseg $(nseg) \\
  --npts $(npts) \\
  --numpts $(numpts) \\
  --channel $(channel) \\
  --threshold $(alphath) \\
  --start_time $(start) \\ 
  --framecache cache/frcache-$(site)-$(frstart)-$(frend).out \\
  --calcache $(calfile) \\
  --injfile HL-INJECTIONS_1_$(usertag)-$(frstart)-$(frduration).xml \\
  --comment $(usertag) \\
  --cluster \\
 """ % (self.config['condor']['universe'],self.config['condor']['burst']),
    for sec in ['arg']:
      for arg in self.config[sec].keys():
        if boolargs.match(arg):
          if re.match(self.config[sec][arg],'true'): 
            print >> sub_fh, "--" + arg,
        else:
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """
log = %s
error = burst/burst-$(ifo)-$(start)-$(end).$(cluster).$(process).err
output = burst/burst-$(ifo)-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.logfile)

    

  def coincsub(self):

    sub_fh = open( self.basename + '.coinc1.sub', 'w' )
    print >> sub_fh, """\
universe = standard
executable = %s
arguments = --ifo-a H1-$(usertag)-POWER-$(start)-$(duration).xml \\
  --ifo-b H2-$(usertag)-POWER-$(start)-$(duration).xml \\
  --outfile coincidence-$(usertag)-POWER-$(start)-$(duration).xml --dt 250
log = %s
error = coin/coinc-$(start)-$(end).$(cluster).$(process).err
output = coin/coinc-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.config['condor']['coinc'],self.logfile)

    sub_fh = open( self.basename + '.coinc2.sub', 'w' )
    print >> sub_fh, """\
universe = standard
executable = %s
arguments = --ifo-a L1-$(usertag)-POWER-$(start)-$(duration).xml \\
  --ifo-b coincidence-$(usertag)-POWER-$(start)-$(duration).xml \\
  --outfile H1H2L1-$(usertag)-POWER-$(start)-$(duration).xml --dt 250
log = %s
error = coin/H1H2L1-$(start)-$(end).$(cluster).$(process).err
output = coin/H1H2L1-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.config['condor']['coinc'],self.logfile)

    



  def builddag(self,cache,burst,simu,coin):
    chanH1 = self.config['input']['channel-name']
    chanH2 = self.config['input']['channelh-name']
    chanL1 = self.config['input']['channell-name']
    siteH1 = chanH1[0]
    siteH2 = chanH2[0]
    siteL1 = chanL1[0]
    ifoH1  = chanH1[0:2]
    ifoH2  = chanH2[0:2]
    ifoL1  = chanL1[0:2]
    t      = int(self.config['datacond']['t'])
    srate  = int(self.config['arg']['srate'])
    olap   = int(self.config['arg']['olap'])
    usertag = self.config['input']['user-tag']
    cala   = self.config['calib']['cala']
    calb   = self.config['calib']['calb']
    alphathL = float(self.config['input']['alphal'])
    alphathH = float(self.config['input']['alphah'])

    print alphathL
    
    dag_fh = open( self.basename + ".dag", "w" )
    
    # jobs to generate the frame cache files
    for seg in self.segments:
      jobname = 'frcache_%s_%d_%d' % (siteH1,seg.startpad,seg.endpad)
      print >> dag_fh, 'JOB %s %s.frcache.sub' % (jobname,self.basename),
      if cache: print >> dag_fh, 'DONE',
      print >> dag_fh, '\nVARS %s site="%s" frstart="%s" frend="%s"' % (
      jobname, siteH1, seg.startpad, seg.endpad )
    for i in range(1,len(self.segments)):
      print >> dag_fh, 'PARENT frcache_%s_%s_%s CHILD frcache_%s_%s_%s' % (
        siteH1,self.segments[i-1].startpad,self.segments[i-1].endpad,
        siteH1,self.segments[i].startpad,self.segments[i].endpad)
    for seg in self.segments:
      jobname = 'frcache_%s_%d_%d' % (siteL1,seg.startpad,seg.endpad)
      print >> dag_fh, 'JOB %s %s.frcache.sub' % (jobname,self.basename),
      if cache: print >> dag_fh, 'DONE',
      print >> dag_fh, '\nVARS %s site="%s" frstart="%s" frend="%s"' % (
      jobname, siteL1, seg.startpad, seg.endpad )
    for i in range(1,len(self.segments)):
      print >> dag_fh, 'PARENT frcache_%s_%s_%s CHILD frcache_%s_%s_%s' % (
        siteL1,self.segments[i-1].startpad,self.segments[i-1].endpad,
        siteL1,self.segments[i].startpad,self.segments[i].endpad)


   # jobs to generate the injection files
    for seg in self.segments:
      jobname = 'binj_%d_%d' % (seg.startpad,seg.endpad)
      print >> dag_fh, 'JOB %s %s.binj.sub' % (jobname,self.basename),
      if simu: print >> dag_fh, 'DONE',
      print >> dag_fh, '\nVARS %s frstart="%s" frend="%s" usertag="%s" ' % (
      jobname, seg.startpad, seg.endpad, usertag )
    for i in range(1,len(self.segments)):
      print >> dag_fh, 'PARENT binj_%s_%s CHILD binj_%s_%s' % (
        self.segments[i-1].startpad,self.segments[i-1].endpad,
        self.segments[i].startpad,self.segments[i].endpad)


    
    # jobs to run the burst code
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'frcache_%s_%s_%s' % (siteH1,seg.startpad,seg.endpad)
        parentB = 'binj_%s_%s' % (seg.startpad,seg.endpad)
        jobname = 'burst_%s_%s_%s' % (ifoH1,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.burst.sub' % (jobname,self.basename),
        if burst: print >> dag_fh, 'DONE',
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d" frduration="%s" usertag="%s" alphath="%g"\
""" % ( jobname,siteH1,ifoH1,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanH1,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)),seg.endpad-seg.startpad, usertag, alphathH)
        print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'frcache_%s_%s_%s' % (siteL1,seg.startpad,seg.endpad)
        parentB = 'binj_%s_%s' % (seg.startpad,seg.endpad)
        jobname = 'burst_%s_%s_%s' % (ifoL1,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.burst.sub' % (jobname,self.basename),
        if burst: print >> dag_fh, 'DONE',        
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d" frduration="%s" usertag="%s" alphath="%g"\
""" % ( jobname,siteL1,ifoL1,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanL1,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)),seg.endpad-seg.startpad, usertag, alphathL)
        print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'frcache_%s_%s_%s' % (siteH2,seg.startpad,seg.endpad)
        parentB = 'binj_%s_%s' % (seg.startpad,seg.endpad)
        jobname = 'burst_%s_%s_%s' % (ifoH2,chunk.start,chunk.end)
        if chunk.end < 731849040:
          print >> dag_fh, 'JOB %s %s.burstH2.sub' % (jobname,self.basename),
          if burst: print >> dag_fh, 'DONE',          
          print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d" frduration="%s" usertag="%s" calfile = "%s" alphath="%g"\
""" % ( jobname,siteH2,ifoH2,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanH2,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)),seg.endpad-seg.startpad, usertag, cala, alphathH)
          print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)
        else:
          print >> dag_fh, 'JOB %s %s.burstH2.sub' % (jobname,self.basename),
          if burst: print >> dag_fh, 'DONE',          
          print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d" frduration="%s" usertag="%s" calfile = "%s" alphath="%g"\
""" % ( jobname,siteH2,ifoH2,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanH2,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)),seg.endpad-seg.startpad, usertag, calb, alphathH)
          print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)



     # jobs to run the coinc code
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'burst_%s_%s_%s' % (ifoH1,chunk.start,chunk.end)
        parentB = 'burst_%s_%s_%s' % (ifoH2,chunk.start,chunk.end)
        jobname = 'coinc_1_%s_%s' % (chunk.start,chunk.end)
        if coin:
          print >> dag_fh, 'JOB %s %s.coinc1.sub' % (jobname,self.basename),
          print >> dag_fh, """
VARS %s ifo="%s" ifo="%s" start="%d" end="%d" chunklen="%d" duration="%d" usertag="%s"\
""" % ( jobname,ifoH2,ifoH1,chunk.start,chunk.end,chunk.length,(chunk.end-chunk.start+1),usertag)
          print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'burst_%s_%s_%s' % (ifoL1,chunk.start,chunk.end)
        parentB = 'coinc_1_%s_%s' % (chunk.start,chunk.end)
        jobname = 'coinc_2_%s_%s' % (chunk.start,chunk.end)
        if coin:
          print >> dag_fh, 'JOB %s %s.coinc2.sub' % (jobname,self.basename),
          print >> dag_fh, """
VARS %s ifo="%s" start="%d" end="%d" chunklen="%d" duration="%d" usertag="%s"\
""" % ( jobname,ifoL1,chunk.start,chunk.end,chunk.length,(chunk.end-chunk.start+1),usertag)
          print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)


    dag_fh.close()
     
def usage():
  msg = """\
Usage: bursttriplifo_pipeline.py [OPTIONS] 

   -f, --config-file FILE   use configuration file FILE
   -p, --play               only create chunks that overlap with playground
   -v, --version            print version information and exit
   -h, --help               print help information
   -c, --cache              flag the frame cache query as done
   -s, --simu               flag the simulated injections as done
   -b, --burst              flag the burst jobs as done
   -l, --log-path           directory to write condor log file

This program generates a DAG to run the burst code. The configuration file 
should specify the parameters needed to run the jobs and must be specified 
with the --config-file (or -f) option. 

A directory which condor uses to write a log file must be specified with the
--log-file (or -l) option. This must be a non-NFS mounted directory. The name
of the log file is automatically created and will be unique for each
invocation of bursttriplifo_pipeline.py.

A file containing science segments to be analyzed should be specified in the 
[input] section of the configuration file with a line such as

segments = S2H1H2L1Science_3.txt

This should contain four whitespace separated coulumns:

  segment_id    gps_start_time  gps_end_time    duration

that define the science segments to be used. Lines starting with # are ignored.

The chunks start and stop times are computed from the science segment times
and used to build the DAG.

If the --play (or -p) flag is specified, then only chunks that overlap the S2
playground defined by:

  ((t - 729273613) % 6370) < 600

are included in the DAG.

If the --cache (or -c) option is specifed, the generation of the frame cache 
files is marked as done in the DAG.

If the --simu (or -s) option is specified, the generation of the injection
files are marked as done in the DAG.

If the --burst (or -b) option is specified, the excess power jobs are marked
as done in the DAG.

\
"""
  print msg

shortop = "f:vhcsbkpl:"
longop = [
  "config-file=",
  "version",
  "help",
  "cache",
  "simu",
  "burst",
  "coin",
  "play",
  "log-path="
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

config_file = None
no_cache = None
no_simu = None
no_burst = None
no_coin = None
play_only = None
log_path = None

for o, a in opts:
  if o in ("-v", "--version"):
    print """$Id$"""
    sys.exit(0)
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  if o in ("-f", "--config-file"):
    config_file = a
  if o in ("-c", "--cache"):
    no_cache = 1
  if o in ("-b", "--burst"):
    no_burst = 1
  if o in ("-k", "--coin"):
    no_coin = 1
  if o in ("-s", "--simu"):
    no_simu = 1  
  if o in ("-p", "--play"):
    play_only = 1
  if o in ("-l", "--log-path"):
    log_path = a

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use -f FILE or --config-file FILE to specify location."
  sys.exit(1)

if not log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use -l PATH or --log-path PATH to specify a location."
  sys.exit(1)
tempfile.tempdir = log_path

try: os.mkdir('cache')
except: pass
try: os.mkdir('datafind')
except: pass
try: os.mkdir('burst')
except: pass
try: os.mkdir('coin')
except: pass
try: os.mkdir('binj')
except: pass

pipeline = BurstPipeline(config_file)
pipeline.parsesegs()
pipeline.createchunks(play_only)
pipeline.status(play_only)
pipeline.frcachesub()
pipeline.burstsub()
pipeline.coincsub()
pipeline.binjsub()
pipeline.builddag(no_cache,no_burst,no_simu,no_coin)
