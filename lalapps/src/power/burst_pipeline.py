#!/usr/bin/env python2.2
"""
burst_pipeline.py - standalone burst pipeline driver script
This has been created by modifying Duncan Brown's original script for running
the standalone inspiral code

$Id$

This script produced the necessary condor submit and dag files to run
the standalone burst code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
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
  define the data to analyze and hos it should be analyzed
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
environment = LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH)
log = %s
error = datafind/frcache-$(site)-$(frstart)-$(frend).$(cluster).$(process).err
output = cache/frcache-$(site)-$(frstart)-$(frend).out
notification = never
queue
""" % (self.config['condor']['datafind'],
       self.config['input']['datatype'],
       self.logfile)
    sub_fh.close()



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
  --start_time $(start) \\ 
  --framecache cache/frcache-$(site)-$(frstart)-$(frend).out \\
  --comment realdata-$(start) \\
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

    sub_fh = open( self.basename + '.coinc.sub', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments = --ifo-a L-realdata-$(start)-POWER-$(start)-$(duration).xml \\
  --ifo-b H-realdata-$(start)-POWER-$(start)-$(duration).xml \\
  --outfile coincidence-$(start)-$(end).xml --dt 10
log = %s
error = coin/coinc-$(start)-$(end).$(cluster).$(process).err
output = coin/coinc-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.config['condor']['universe'],self.config['condor']['coinc'],self.logfile),



  def builddag(self,cache):
    chanH1 = self.config['input']['channel-name']
    chanL1 = self.config['input']['channell-name']
    siteH1 = chanH1[0]
    siteL1 = chanL1[0]
    ifoH1  = chanH1[0:2]
    ifoL1  = chanL1[0:2]
    t      = int(self.config['datacond']['t'])
    srate  = int(self.config['arg']['srate'])
    olap   = int(self.config['arg']['olap'])
    
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


    # jobs to run the burst code
    for seg in self.segments:
      for chunk in seg.chunks:
        parent = 'frcache_%s_%s_%s' % (siteH1,seg.startpad,seg.endpad)
        jobname = 'burst_%s_%s_%s' % (ifoH1,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.burst.sub' % (jobname,self.basename),
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d"\
""" % ( jobname,siteH1,ifoH1,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanH1,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)))
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)
    for seg in self.segments:
      for chunk in seg.chunks:
        parent = 'frcache_%s_%s_%s' % (siteL1,seg.startpad,seg.endpad)
        jobname = 'burst_%s_%s_%s' % (ifoL1,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.burst.sub' % (jobname,self.basename),
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" nseg = "%d" npts = "%d" numpts = "%d"\
""" % ( jobname,siteL1,ifoL1,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chanL1,(2*chunk.length-1),(t*srate),(((t*srate)-olap)*(2*chunk.length-1)+(3*olap)))
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)



     # jobs to run the coinc code
    for seg in self.segments:
      for chunk in seg.chunks:
        parentA = 'burst_%s_%s_%s' % (ifoH1,chunk.start,chunk.end)
        parentB = 'burst_%s_%s_%s' % (ifoL1,chunk.start,chunk.end)
        jobname = 'coinc_%s_%s' % (chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.coinc.sub' % (jobname,self.basename),
        print >> dag_fh, """
VARS %s ifo="%s" ifo="%s" start="%d" end="%d" chunklen="%d" duration="%d"\
""" % ( jobname,ifoH1,ifoL1,chunk.start,chunk.end,chunk.length,(chunk.end-chunk.start+1))
        print >> dag_fh, 'PARENT %s %s CHILD %s' % (parentA, parentB, jobname)    


    dag_fh.close()
     
def usage():
  msg = """\
Usage: burst_pipeline.py [OPTIONS]

   -f, --config-file FILE   use configuration file FILE
   -p, --play               only create chunks that overlap with playground
   -v, --version            print version information and exit
   -h, --help               print help information
   -c, --cache              flag the frame cache query as done
   -b, --bank               flag the bank generation and cache query as done
   -l, --log-path           directory to write condor log file

This program generates a DAG to run the inspiral code. The configuration file 
should specify the parameters needed to run the jobs and must be specified 
with the --config-file (or -f) option. See the LALapps documentation for more
information on the syntax of the configuation file.

A directory which condor uses to write a log file must be specified with the
--log-file (or -l) option. This must be a non-NFS mounted directory. The name
of the log file is automatically created and will be unique for each
invocation of inspiral_pipeline.py.

A file containing science segments to be analyzed should be specified in the 
[input] section of the configuration file with a line such as

segments = S2TripleCoincidetScienceSegments.txt

This should contain four whitespace separated coulumns:

  segment_id    gps_start_time  gps_end_time    duration

that define the science segments to be used. Lines starting with # are ignored.

The length of the number of inspiral segments, their overlap and length is 
determined from the config file and the length of an inspiral chunk is
computed (typically this is 1024 seconds for S2).

The chunks start and stop times are computed from the science segment times
and used to build the DAG.

If the --play (or -p) flag is specified, then only chunks that overlap the S2
playground defined by:

  ((t - 729273613) % 6370) < 600

are included in the DAG.

If the --cache (or -c) option is specifed, the generation of the frame cache 
files is marked as done in the DAG. The cache files are expected to exist by
the bank generation and inspiral codes.

If the --bank (or -b) option is specified, both the frame cache query and
generation of the template bank are marked as done. The cache files and
template banks are expected to exist by the inspiral code.
\
"""
  print msg

shortop = "f:vhcbipl:"
longop = [
  "config-file=",
  "version",
  "help",
  "cache",
  "burst",
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

pipeline = BurstPipeline(config_file)
pipeline.parsesegs()
pipeline.createchunks(play_only)
pipeline.status(play_only)
pipeline.frcachesub()
pipeline.burstsub()
pipeline.coincsub()
pipeline.builddag(no_cache)
