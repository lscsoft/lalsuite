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
    # compute the chunk and overlap length in seconds
    numpoints = int(self.config['datacond']['segment-length'])
    numseg = int(self.config['datacond']['number-of-segments'])
    overlap = int(self.config['datacond']['segment-overlap'])
    srate = int(self.config['datacond']['sample-rate'])
    self.chunk = (numpoints * numseg - ( numseg - 1 ) * overlap ) / srate
    self.overlap = int(self.config['datacond']['segment-overlap']) / srate

  def parsesegs(self):
    self.segments = []
    comment_line = re.compile(r'\A#')
    for line in open(self.config['input']['segments']):
      if not comment_line.match(line):
        self.segments.append(
          { 'segment' : tuple(map(int,line.split())), 'chunks' : [] } )

  def buildchunks(self):
    for seg in self.segments:
      (id, start, end, length) = seg['segment']
      while length >= self.chunk:
        seg['chunks'].append(tuple([start,start+self.chunk]))
        start += self.chunk - self.overlap
        length -= self.chunk - self.overlap

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
    args = ()
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
    print >> sub_fh, """\
log = %s.log
error = tmpltbank-$(start)-$(end).err
output = tmpltbank-$(start)-$(end).out
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
error = tmpltbank-$(start)-$(end).err
output = tmpltbank-$(start)-$(end).out
notification = never
queue
""" % (self.basename)

  def builddag(self,cache,bank,inspiral):
    chan = self.config['input']['channel-name']
    site = chan[0]
    ifo  = chan[0:2]
    dag_fh = open( self.basename + ".dag", "w" )
    print >> dag_fh, "dot %s.dot update overwrite" % self.basename
    
    # jobs to generate the frame cache files
    for seg in self.segments:
      (id, start, end, length) = seg['segment']
      jobname = 'frcache_%s_%s_%s' % (site,start,end)
      print >> dag_fh, 'job %s %s.frcache.condor' % (jobname,self.basename),
      if not cache: print >> dag_fh, 'done',
      print >> dag_fh, '\nvars %s site="%s"' % (jobname,site)
      print >> dag_fh, 'vars %s frstart="%s"' % (jobname,start)
      print >> dag_fh, 'vars %s frend="%s"' % (jobname,end)
    for i in range(1,len(self.segments)):
      (id, start_p, end_p, length) = self.segments[i-1]['segment']
      (id, start_c, end_c, length) = self.segments[i]['segment']
      print >> dag_fh, 'parent frcache_%s_%s_%s child frcache_%s_%s_%s' % (
        site,start_p,end_p,site,start_c,end_c)
    
    # jobs to generate the template banks
    for seg in self.segments:
      (id, frstart, frend, length) = seg['segment']
      frparent = 'frcache_%s_%s_%s' % (site,frstart,frend)
      for start,end in seg['chunks']:
        jobname = 'tmpltbank_%s_%s_%s' % (ifo,start,end)
        print >> dag_fh, 'job %s %s.tmpltbank.condor' % (jobname,self.basename),
        if not cache: print >> dag_fh, 'done',
        print >> dag_fh, '\nvars %s site="%s"' % (jobname,site)
        print >> dag_fh, 'vars %s frstart="%s"' % (jobname,start)
        print >> dag_fh, 'vars %s frend="%s"' % (jobname,end)
        print >> dag_fh, 'vars %s start="%d"' % (jobname,start)
        print >> dag_fh, 'vars %s end="%d"' % (jobname,end)
        print >> dag_fh, 'vars %s channel="%s"' % (jobname,chan)
        print >> dag_fh, 'vars %s calcache="%s"' % (jobname,
          self.config['input'][string.lower(ifo) + '-cal'])
        print >> dag_fh, 'parent %s child %s' % (frparent, jobname)
    dag_fh.close()

  def status(self):
    print "writing to files with basename", self.basename
    print "generating", self.chunk, "second chunks",
    print "with", self.overlap, "second overlap"
    try:
      print self.segments
    except AttributeError:
      print "no science segments defined"

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
pipeline.buildchunks()
pipeline.frcachesub()
pipeline.banksub()
pipeline.inspiralsub()
pipeline.builddag(do_cache,do_bank,do_inspiral)
pipeline.status()
