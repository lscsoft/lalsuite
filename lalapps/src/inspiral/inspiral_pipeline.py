#!/usr/bin/env python2.2
"""
inspiral_pipeline.py - standalone inspiral pipeline driver script

$Id$

This script produced the necessary condor submit and dag files to run
the standalone inspiral code on LIGO data
"""

__author__ = 'Duncan Brown <duncan@gravity.phys.uwm.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import os
import sys
import string
import getopt
import re
import ConfigParser

class InspiralChunk:
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
    while dur >= length:
      middle = start + length / 2
      end = start + length
      if (not play_only) or ( play_only 
        and ( isplay(start) or isplay(middle) or isplay(end) ) ):
        self.chunks.append(InspiralChunk(start,end))
      start += seg_incr
      dur -= seg_incr
    self.used = self.duration - dur

  def isplay(self,t):
    if ((t - 729273613) % 6370) < 600:
      return 1
    else:
      return 0

class InspiralPipeline:
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

  def parsesegs(self):
    self.segments = []
    segpad = int(self.config['inspiral']['pad-data'])
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
    numpoints = int(self.config['datacond']['segment-length'])
    numseg = int(self.config['datacond']['number-of-segments'])
    overlap = int(self.config['inspiral']['segment-overlap'])
    srate = int(self.config['datacond']['sample-rate'])
    chunk_length = (numpoints * numseg - ( numseg - 1 ) * overlap ) / srate
    chunk_overlap = overlap / srate
    for seg in self.segments:
      seg.createchunks(chunk_length,chunk_overlap,play_only)

  def status(self):
    log_fh = open( self.basename + '.pipeline.log', 'w' )
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
created %d inspiral chunks for analysis
""" % ( len(self.segments), total_science, total_science_used, num_chunks )
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
log = %s.log
error = datafind/frcache-$(site)-$(frstart)-$(frend).$(cluster).$(process).err
output = cache/frcache-$(site)-$(frstart)-$(frend).out
notification = never
queue
""" % (self.config['condor']['datafind'],
       self.config['input']['datatype'],
       self.basename)
    sub_fh.close()

  def banksub(self):
    boolargs = re.compile(r'(disable-high-pass|write-strain-spectrum|verbose)')
    sub_fh = open( self.basename + '.tmpltbank.sub', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments = --gps-start-time $(start) --gps-end-time $(end) \\
  --channel-name $(channel) --calibration-cache $(calcache) \\
  --frame-cache cache/frcache-$(site)-$(frstart)-$(frend).out \\
 """ % (self.config['condor']['universe'],self.config['condor']['tmpltbank']),
    for sec in ['datacond','bank']:
      for arg in self.config[sec].keys():
        if boolargs.match(arg):
          if re.match(self.config[sec][arg],'true'): 
            print >> sub_fh, "--" + arg,
        else:
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """
log = %s.log
error = bank/tmpltbank-$(ifo)-$(start)-$(end).$(cluster).$(process).err
output = bank/tmpltbank-$(ifo)-$(start)-$(end).$(cluster).$(process).out
notification = never
queue""" % self.basename
    sub_fh.close()

  def inspiralsub(self):
    boolargs = re.compile(r'(disable-high-pass|enable-event-cluster|verbose|enable-output|disable-output)')
    sub_fh = open( self.basename + '.inspiral.sub', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments = --gps-start-time $(start) --gps-end-time $(end) \\
  --channel-name $(channel) --calibration-cache $(calcache) \\
  --frame-cache cache/frcache-$(site)-$(frstart)-$(frend).out \\
  --bank-file $(site)-TMPLTBANK-$(start)-$(chunklen).xml \\
 """ % (self.config['condor']['universe'],self.config['condor']['inspiral']),
    for sec in ['datacond','inspiral']:
      for arg in self.config[sec].keys():
        if boolargs.match(arg):
          if re.match(self.config[sec][arg],'true'): 
            print >> sub_fh, "--" + arg,
        else:
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """
log = %s.log
error = inspiral/inspiral-$(ifo)-$(start)-$(end).$(cluster).$(process).err
output = inspiral/inspiral-$(ifo)-$(start)-$(end).$(cluster).$(process).out
notification = never
queue
""" % (self.basename)

  def builddag(self,cache,bank,inspiral):
    chan = self.config['input']['channel-name']
    site = chan[0]
    ifo  = chan[0:2]
    dag_fh = open( self.basename + ".dag", "w" )
    
    # jobs to generate the frame cache files
    for seg in self.segments:
      jobname = 'frcache_%s_%d_%d' % (site,seg.startpad,seg.endpad)
      print >> dag_fh, 'JOB %s %s.frcache.sub' % (jobname,self.basename),
      if cache: print >> dag_fh, 'done',
      print >> dag_fh, '\nVARS %s site="%s" frstart="%s" frend="%s"' % (
      jobname, site, seg.startpad, seg.endpad )
    for i in range(1,len(self.segments)):
      print >> dag_fh, 'PARENT frcache_%s_%s_%s CHILD frcache_%s_%s_%s' % (
        site,self.segments[i-1].startpad,self.segments[i-1].endpad,
        site,self.segments[i].startpad,self.segments[i].endpad)
    
    # jobs to generate the template banks
    for seg in self.segments:
      parent = 'frcache_%s_%s_%s' % (site,seg.startpad,seg.endpad)
      for chunk in seg.chunks:
        jobname = 'tmpltbank_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.tmpltbank.sub' % (jobname,self.basename),
        if bank: print >> dag_fh, 'done',
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" channel="%s" calcache="%s"\
""" % ( jobname,site,ifo,seg.startpad,seg.endpad,chunk.start,chunk.end,chan,
self.config['input'][string.lower(ifo) + '-cal'])
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)

    # jobs to run the inspiral code
    for seg in self.segments:
      for chunk in seg.chunks:
        parent = 'tmpltbank_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        jobname = 'inspiral_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.inspiral.sub' % (jobname,self.basename),
        if inspiral: print >> dag_fh, 'done',
        print >> dag_fh, """
VARS %s site="%s" ifo="%s" frstart="%s" frend="%s" start="%d" end="%d" chunklen="%d" channel="%s" calcache="%s"\
""" % ( jobname,site,ifo,seg.startpad,seg.endpad,chunk.start,chunk.end,
chunk.length,chan,self.config['input'][string.lower(ifo) + '-cal'])
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)
    dag_fh.close()

def usage():
  msg = """\
Usage: inspiral_pipeline.py [OPTIONS]

   -f, --config-file FILE       use configuration file FILE
   -v, --version                print version information and exit
   -h, --help                   print help information
   -c, --cache                  flag the LALdataFind query as done
   -b, --bank                   flag the template bank generation as done
   -i, --inspiral               flag the inspiral code as done
   -p, --play                   only create chunks that overlap with playground

This program generates a DAG to run the inspiral code. The configuration
file should specify the parameters needed to run the jobs and must be
specified with the --config-file (or -f) option.

Unless the --cache (or -c) option is specified the DAG will generate frame
cache files by querying LDR, otherwise these will be expected to exist.

The DAG will template bank generation job marked as done unless the
--bank (or -b) option is given. Similarly the DAG will have the inspiral jobs
marked as done unless the --inspiral (or -i) option is given. If the bank
generation is not specified and the inspiral is specifed, the DAG will 
expect the banks to have been previously generated with the correct names.
\
"""
  print msg

shortop = "f:vhcbip"
longop = [
  "config-file=",
  "version",
  "help",
  "cache",
  "bank",
  "inspiral",
  "play"
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

config_file = None
do_cache = None
do_bank = None
do_inspiral = None
play_only = None

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
    do_cache = 1
  if o in ("-b", "--bank"):
    do_bank = 1
  if o in ("-i", "--inspiral"):
    do_inspiral = 1
  if o in ("-p", "--play"):
    play_only = 1

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use -f FILE or --config-file FILE to specify location."
  sys.exit(1)

try: os.mkdir('cache')
except: pass
try: os.mkdir('datafind')
except: pass
try: os.mkdir('bank')
except: pass
try: os.mkdir('inspiral')
except: pass

pipeline = InspiralPipeline(config_file)
pipeline.parsesegs()
pipeline.createchunks(play_only)
pipeline.status()
pipeline.frcachesub()
pipeline.banksub()
pipeline.inspiralsub()
pipeline.builddag(do_cache,do_bank,do_inspiral)
