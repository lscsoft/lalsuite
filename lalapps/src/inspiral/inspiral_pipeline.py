#!/usr/bin/python2.2
"""
pipeline.py - standalone inspiral pipeline driver script

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

class ScienceSegment:
  def __init__(self,id,start,end,duration):
    self.id = id
    self.start = start
    self.end = end
    self.duration = duration
    self.chunks = []
  def createchunks(self,length,overlap):
    dur = self.duration
    start = self.start
    incr = dur - overlap
    while dur >= length:
      self.chunks.append(InspiralChunk(start,start+length))
      dur =- incr
      start += incr

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
    # lines that start with an octothorpe are comments
    comment_line = re.compile(r'\A#')
    for line in open(self.config['input']['segments']):
      if not comment_line.match(line):
        segpars = tuple(map(int,line.split()))
        # change this line if the format of the science segment file changes
        self.segments.append( 
          ScienceSegment(segpars[0],segpars[1],segpars[2],segpars[3]) )
  
  def createchunks(self):
    # compute the chunk and overlap length in seconds
    numpoints = int(self.config['datacond']['segment-length'])
    numseg = int(self.config['datacond']['number-of-segments'])
    overlap = int(self.config['inspiral']['segment-overlap'])
    srate = int(self.config['datacond']['sample-rate'])
    chunk_length = (numpoints * numseg - ( numseg - 1 ) * overlap ) / srate
    chunk_overlap = overlap / srate
    for seg in self.segments:
      seg.createchunks(chunk_length,chunk_overlap)

  def frcachesub(self):
    sub_fh = open( self.basename + '.frcache.condor', 'w' )
    print >> sub_fh, """\
universe = scheduler
executable = %s
arguments = --lal-cache \\
  --instrument $(site) --type %s \\
  --start $(frstart) --end $(frend)
log = %s.log
error = frcache-$(site)-$(frstart)-$(frend).err
output = frcache-$(site)-$(frstart)-$(frend).out
notification = never
queue
""" % (self.config['condor']['datafind'],
       self.config['input']['datatype'],
       self.basename)
    sub_fh.close()

  def banksub(self):
    boolargs = re.compile(r'(disable-high-pass|write-strain-spectrum)')
    sub_fh = open( self.basename + '.tmpltbank.condor', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments = --gps-start-time $(start) --gps-end-time $(end) \\
  --channel-name $(channel) --calibration-cache $(calcache) \\
  --frame-cache frcache-$(site)-$(frstart)-$(frend).out \\
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
error = tmpltbank-$(ifo)-$(start)-$(end).err
output = tmpltbank-$(ifo)-$(start)-$(end).out
notification = never
queue""" % self.basename
    sub_fh.close()

  def inspiralsub(self):
    boolargs = re.compile(r'(disable-high-pass|enable-event-cluster)')
    sub_fh = open( self.basename + '.inspiral.condor', 'w' )
    print >> sub_fh, """\
universe = %s
executable = %s
arguments = --gps-start-time $(start) --gps-end-time $(end) \\
  --channel-name $(channel) --calibration-cache $(calcache) \\
  --frame-cache frcache-$(site)-$(frstart)-$(frend).out \\
 """ % (self.config['condor']['datafind'],self.config['condor']['tmpltbank']),
    for sec in ['datacond','inspiral']:
      for arg in self.config[sec].keys():
        if boolargs.match(arg):
          if re.match(self.config[sec][arg],'true'): 
            print >> sub_fh, "--" + arg,
        else:
          print >> sub_fh, "--" + arg, self.config[sec][arg], 
    print >> sub_fh, """
log = %s.log
error = inspiral-$(ifo)-$(start)-$(end).err
output = inspiral-$(ifo)-$(start)-$(end).out
notification = never
queue
""" % (self.basename)

  def builddag(self,cache,bank,inspiral):
    chan = self.config['input']['channel-name']
    site = chan[0]
    ifo  = chan[0:2]
    dag_fh = open( self.basename + ".dag", "w" )
    print >> dag_fh, "DOT %s.dot UPDATE OVERWRITE" % self.basename
    
    # jobs to generate the frame cache files
    for seg in self.segments:
      jobname = 'frcache_%s_%d_%d' % (site,seg.start,seg.end)
      print >> dag_fh, 'JOB %s %s.frcache.condor' % (jobname,self.basename),
      if not cache: print >> dag_fh, 'done',
      print >> dag_fh, '\nVARS %s site="%s"' % (jobname,site)
      print >> dag_fh, 'VARS %s frstart="%s"' % (jobname,seg.start)
      print >> dag_fh, 'VARS %s frend="%s"' % (jobname,seg.end)
    for i in range(1,len(self.segments)):
      print >> dag_fh, 'PARENT frcache_%s_%s_%s CHILD frcache_%s_%s_%s' % (
        site,self.segments[i-1].start,self.segments[i-1].end,
        site,self.segments[i].start,self.segments[i].end)
    
    # jobs to generate the template banks
    for seg in self.segments:
      parent = 'frcache_%s_%s_%s' % (site,seg.start,seg.end)
      for chunk in seg.chunks:
        jobname = 'tmpltbank_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.tmpltbank.condor' % (jobname,self.basename),
        if not cache: print >> dag_fh, 'done',
        print >> dag_fh, '\nVARS %s site="%s"' % (jobname,site)
        print >> dag_fh, 'VARS %s ifo="%s"' % (jobname,ifo)
        print >> dag_fh, 'VARS %s frstart="%s"' % (jobname,seg.start)
        print >> dag_fh, 'VARS %s frend="%s"' % (jobname,seg.end)
        print >> dag_fh, 'VARS %s start="%d"' % (jobname,chunk.start)
        print >> dag_fh, 'VARS %s end="%d"' % (jobname,chunk.end)
        print >> dag_fh, 'VARS %s channel="%s"' % (jobname,chan)
        print >> dag_fh, 'VARS %s calcache="%s"' % (jobname,
          self.config['input'][string.lower(ifo) + '-cal'])
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)

    # jobs to run the inspiral code
    for seg in self.segments:
      for chunk in seg.chunks:
        parent = 'tmpltbank_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        jobname = 'inspiral_%s_%s_%s' % (ifo,chunk.start,chunk.end)
        print >> dag_fh, 'JOB %s %s.inspiral.condor' % (jobname,self.basename),
        if not cache: print >> dag_fh, 'done',
        print >> dag_fh, '\nVARS %s site="%s"' % (jobname,site)
        print >> dag_fh, 'VARS %s ifo="%s"' % (jobname,ifo)
        print >> dag_fh, 'VARS %s frstart="%s"' % (jobname,seg.start)
        print >> dag_fh, 'VARS %s frend="%s"' % (jobname,seg.end)
        print >> dag_fh, 'VARS %s start="%d"' % (jobname,chunk.start)
        print >> dag_fh, 'VARS %s end="%d"' % (jobname,chunk.end)
        print >> dag_fh, 'VARS %s channel="%s"' % (jobname,chan)
        print >> dag_fh, 'VARS %s calcache="%s"' % (jobname,
          self.config['input'][string.lower(ifo) + '-cal'])
        print >> dag_fh, 'PARENT %s CHILD %s' % (parent, jobname)
    dag_fh.close()

def usage():
  msg = """\
Usage: pipeline.py [OPTIONS]

   -f, --config-file FILE       use configuration file FILE
   -v, --version                print version information and exit
   -h, --help                   print help information
   -c, --cache                  query LDR to generate cache files
   -b, --bank                   run the bank generation code in the DAG
   -i, --inspiral               run the inspiral code in the DAG

This program generates a DAG to run the inspiral code. The configuration
file should specify the parameters needed to run the jobs and must be
specified with the --config-file (or -f) option.

If the --cache (or -c) option is specified the DAG will generate frame
cache files by querying LDR, otherwise these will be expected to exist.

The DAG will have each template bank generation job marked as done unless the
--bank (or -b) option is given. Similarly the DAG will have the inspiral jobs
marked as done unless the --inspiral (or -i) option is given. If the bank
generation is not specified and the inspiral is specifed, the DAG will 
expect the banks to have been previously generated with the correct names.
\
"""
  print msg

shortop = "f:vhcbi"
longop = [
  "config-file=",
  "version",
  "help",
  "cache",
  "bank",
  "inspiral"
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

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use -f FILE or --config-file FILE to specify location."
  sys.exit(1)

pipeline = InspiralPipeline(config_file)
pipeline.parsesegs()
pipeline.createchunks()
pipeline.frcachesub()
pipeline.banksub()
pipeline.inspiralsub()
pipeline.builddag(do_cache,do_bank,do_inspiral)
