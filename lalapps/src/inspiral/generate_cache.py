#!/usr/bin/env @PYTHONPROG@
"""
qscan.in - simple dag generator for q scans

$Id$

This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
#import inspiral

class getCache(UserDict):
  """
  An instance of a lal cache
  """
  def __init__(self, options):
    UserDict.__init__(self)
    self.dir = os.listdir(options.cache_path)
    self.options = options
  def getCacheType(self, type):
    self[type] = []
    p = re.compile(type)
    m = re.compile("-")
    x = re.compile(".xml")
    for fname in  self.dir:
      if p.search(fname):
        ifo = m.split(fname)[0]
        start = m.split(fname)[-2]
        dur = x.split(m.split(fname)[-1])
        entry = lal.CacheEntry(ifo+" "+self.options.science_run+" " \
                               +start+" "+dur[0]+" "+fname)
        self[type].append(entry)

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-p", "--cache-path",action="store",type="string",\
    metavar=" PATH",help="directory to find XML files")

parser.add_option("-r", "--science-run", action="store",type="string",\
    metavar=" RUN", help="name of science run")


command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#################################
# if --version flagged
if opts.version:
  print "$Id$"
  sys.exit(0)

#################################
if not opts.cache_path:
  print >> sys.stderr, "No cache path specified."
  print >> sys.stderr, "Use --cache-path PATH to specify a location."
  sys.exit(1)

if not opts.science_run:
  print >> sys.stderr, "No science run specified."
  print >> sys.stderr, "Use --science-run RUN to specify a run."
  sys.exit(1)

cache = getCache(opts)

cache.getCacheType("TMPLTBANK")
cache.getCacheType("INSPIRAL_")
cache.getCacheType("INSPIRAL-")
cache.getCacheType("TRIGBANK")
cache.getCacheType("THINCA-")
cache.getCacheType("THINCA_")

# Write bank cache file
bank = open('bank.cache','w')
for fname in cache["TMPLTBANK"]:
  bank.write(str(fname)+"\n")
bank.close()

# Write trigbank cache file
trigbank = open('trigbank.cache','w')
for fname in cache["TRIGBANK"]:
  trigbank.write(str(fname)+"\n")
trigbank.close()

# Write first inspiral cache file
first_inspiral = open('first_inspiral.cache','w')
for fname in cache["INSPIRAL-"]:
  first_inspiral.write(str(fname)+"\n")
first_inspiral.close()

# Write second inspiral cache file
second_inspiral = open('second_inspiral.cache','w')
for fname in cache["INSPIRAL-"]:
  second_inspiral.write(str(fname)+"\n")
second_inspiral.close()

# Write first thinca cache file
first_thinca = open('first_thinca.cache','w')
for fname in cache["THINCA-"]:
  first_thinca.write(str(fname)+"\n")
first_thinca.close()

# Write second thinca cache file
second_thinca = open('second_thinca.cache','w')
for fname in cache["THINCA_"]:
  second_thinca.write(str(fname)+"\n")
second_thinca.close()


