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

def LoadConfig(file):
  """
  returns a dictionary with key's of the form
  <section>.<option> and the values 
  """
  config = {}
  cp = ConfigParser.ConfigParser()
  cp.read(file)
  for sec in cp.sections():
    name = string.lower(sec)
    config[name] = {}
    for opt in cp.options(sec):
      config[name][opt] = string.strip(cp.get(sec,opt))
  return config

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

# load the configuration file into a dictionary
config = LoadConfig(config_file)

# load the science segment file into a dictionary
sci_seg = {}
cmt_ln = re.compile(r'\A#')
for line in open(config['input']['segments']):
  if not cmt_ln.match(line):
    sci_seg.setdefault(tuple(map(int,line.split())), [])

# compute the overlap and chunk length in seconds
numpoints = int(config['datacond']['segment-length'])
numseg = int(config['datacond']['number-of-segments'])
numovrlap = int(config['datacond']['segment-overlap'])
srate = int(config['datacond']['sample-rate'])
ovrlap_len = numovrlap / srate
chunk_len = (numpoints * numseg - ( numseg -1 ) * numovrlap) / srate

# generate the list of chunk start/end tuples
for key in sci_seg.keys():
  (id, start, end, length) = key
  while length >= chunk_len:
    sci_seg[key].append(tuple([start, start + chunk_len]))
    start += (chunk_len - ovrlap_len)
    length -= (chunk_len - ovrlap_len)

ldrcache_sub_fh = open('ldrcache.condor','w')
print >> ldrcache_sub_fh, """\
universe = scheduler
executable = %s
arguments = --lal-cache \\
            --instrument $(ifo) --type %s \\
            --start $(start) --end $(end)
log = inspiral_dag.log
error = ldrcache-$(ifo)-$(start)-$(end).err
output = ldrcache-$(ifo)-$(start)-$(end).out
notification = never
queue
""" % (config['condor']['datafind'],config['input']['datatype'])
ldrcache_sub_fh.close()

bank_sub_fh = open('tmpltbank.condor','w')
print >> bank_sub_fh, """\
universe = %s
executable = %s
arguments = """,
bank_sub_fh.close()

dag_fh = open('inspiral_pipeline.dag','w')

print >> dag_fh, "# inspiral pipeline dag"
print >> dag_fh, "dot inspiral_pipeline.dot update overwrite"

child = None
last  = (None,None)
for key in sci_seg.keys():
  interval = (1,2)
  job = 'ldrcache-%d-%d' % interval
  print >> dag_fh, 'job %s ldrcache.condor' % job
  print >> dag_fh, 'vars %s start="%d"' % (job, interval[0])
  print >> dag_fh, 'vars %s stop="%d"' % (job, interval[1])
  if not child:
    child = 1
  else:
    print >> dag_fh, 'parent %s child %s' % (lastjob,job)
  lastjob = job

dag_fh.close()
