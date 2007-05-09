#!/usr/bin/env @PYTHONPROG@
"""
Something

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
from pylab import *
from glue import segments 
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal.fu_utils import *

########## CLASS TO WRITE LAL CACHE FROM HIPE OUTPUT #########################
class getCache(UserDict):
  """
  An instance of a lal cache
  """
  def __init__(self, options):
    UserDict.__init__(self)
    self.dir = os.listdir(options.cache_path)
    self.options = options
    self.types = ['TMPLTBANK', 'INSPIRAL_', 'INSPIRAL-', \
                 'TRIGBANK', 'THINCA-', 'THINCA_']
    self.oNames = ['bank.cache', 'trigbank.cache', 'first_inspiral.cache' \
         'second_inspiral.cache', 'first_thinca.cache', 'second_thinca.cache'] 
    self.nameMaps = map(None, self.oNames, self.types)
    self.ifoTypes = ['H1','H2','L1','H1H2','H1L1','H2L1','H1H2L1']
    self.ifoDict = {'H1':[],'H2':[],'L1':[],'H1H2':[], \
                    'H1L1':[],'H2L1':[],'H1H2L1':[]}
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
        cache_path = os.path.abspath(self.options.cache_path)
        try:
          entry = lal.CacheEntry(ifo+" "+self.options.science_run+" " \
              +start+" "+dur[0]+" "+"file://localhost"
              +cache_path+"/"+fname)
          self[type].append(entry)
        except:
          pass

  def getCacheAll(self):
    for type in self.types:
      self.getCacheType(type)

  def writeCacheType(self,oName,type):
    cName = open(oName,'w')
    for fname in self[type]:
      cName.write(str(fname)+"\n")
      cName.close()

  def writeCacheAll(self):
    for oName, type in self.nameMaps:
      self.writeCacheType(str(oName),type)
  
  def filesMatchingGPS(self, time, type):
    cacheSubSet = self.ifoDict
    for cache in self[type]:
      for ifo in self.ifoTypes:
       try:
         if (time >= cache.to_segmentlistdict()[ifo][0][0]) and \
            (time < cache.to_segmentlistdict()[ifo][0][1]):
           cacheSubSet[ifo].append(cache)
       except:
         pass
    return cacheSubSet
      

##############################################################################
# redefine the SimInspiral columns of interest
##############################################################################
lsctables.SimInspiralTable.loadcolumns = [
    "waveform",
    "geocent_end_time",
    "geocent_end_time_ns",
    "h_end_time",
    "h_end_time_ns",
    "l_end_time",
    "l_end_time_ns",
    "source",
    "mass1",
    "mass2",
    "mchirp",
    "eta",
    "distance",
    "spin1x",
    "spin1y",
    "spin1z",
    "spin2x",
    "spin2y",
    "spin2z",
    "eff_dist_h",
    "eff_dist_l",
    "eff_dist_g",
    "eff_dist_t",
    "eff_dist_v"]

##############################################################################
# redefine the SnglInspiral columns of interest
##############################################################################
lsctables.SnglInspiralTable.loadcolumns = [
    "ifo",
    "end_time",
    "end_time_ns",
    "eff_distance",
    "mass1",
    "mass2",
    "mchirp",
    "eta",
    "snr",
    "chisq",
    "chisq_dof",
    "sigmasq",
    "event_id"]

      
##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

###### OPTION PARSING AND SANITY CHECKS #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-p", "--cache-path",action="store",type="string",\
    metavar=" PATH",help="directory to find  all hipe XML files")

parser.add_option("-r", "--science-run", action="store",type="string",\
    metavar=" RUN", help="name of science run")

parser.add_option("-g","--xml-glob",action="store",type="string",\
    default=None, metavar=" XML_GLOB", \
    help="GLOB of coire xml files to read (Loudest events)" )

parser.add_option("-K","--statistic",action="store",type="string",\
    default="effective_snrsq",metavar=" STAT",\
    help="coincident statistic (default = effective_snr)")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

if opts.version:
  print "$Id$"
  sys.exit(0)

if not opts.cache_path:
  print >> sys.stderr, "No cache path specified."
  print >> sys.stderr, "Use --cache-path PATH to specify a location."
  sys.exit(1)

if not opts.science_run:
  print >> sys.stderr, "No science run specified."
  print >> sys.stderr, "Use --science-run RUN to specify a run."
  sys.exit(1)

#if not opts.xml_glob:
#  print >> sys.stderr, "Must specify a GLOB of xmls to read"
#  sys.exit(1)

#if not opts.statistic:
#  print >> sys.stderr, "Must specify a statistic to use"
#  sys.exit(1)


############# TURN THE HIPE OUTPUT INTO LAL CACHE FILES #######################

cache = getCache(opts)
cache.getCacheAll()
cache.writeCacheAll()

############# READ IN THE COIRE FILES #########################################

if opts.xml_glob:
  pass
#  found, coincs = readFiles(opts.xml_glob,getstatistic(opts))
#  followuptrigs = getfollowuptrigs(opts,coincs,missed)

sys.exit(0)

