"""
followup utilities

This
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'

import sys
import os, shutil
import urllib
try:
  import sqlite3
except ImportError:
  # pre 2.5.x
  from pysqlite2 import dbapi2 as sqlite3

from subprocess import *
import copy
import re
from StringIO import StringIO
import exceptions
import glob
import fileinput
import linecache
import string
import random
import numpy
import cPickle
import gzip
from scipy import interpolate
from commands import getstatusoutput
import math
import fnmatch
from optparse import *
from types import *
import matplotlib
matplotlib.use('Agg')
import operator
from UserDict import UserDict

from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue.ligolw import dbtables
from pylal import CoincInspiralUtils
from pylal import frutils
from glue import iterutils
from glue import pipeline
from glue.lal import *
from glue import lal
from glue import markup
from lalapps import inspiralutils
from glue.segmentdb import segmentdb_utils
from glue.segmentdb import query_engine
from pylal.xlal import date as xlaldate
#Part of bandaid
from xml import sax
from pylal import db_thinca_rings
from pylal import git_version
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

#####################################################
## Misc Constants used by multiple classes/methods ##
#####################################################
interferometers=["H1","H2","L1","V1"]
defaultsegmentserver="https://segdb.ligo.caltech.edu"
runEpochs={"L1":{
  "s6d_current":(961545687,999999999),
  "s6c":(949449543,961545687),
  "s6b":(937800015,947260815),
  "s6a":(931035296,935798487),
  "s5":(815155213,875232014)
  },
           "H1":{
  "s6d_current":(961545687,999999999),
  "s6c":(949449543,961545687),
  "s6b":(937800015,947260815),
  "s6a":(931035296,935798487),
  "s5":(815155213,875232014)
  },
           "H2":{
  "s6d_current":(961545687,999999999),
  "s6c":(949449543,961545687),
  "s6b":(937800015,947260815),
  "s6a":(931035296,935798487),
  "s5":(815155213,875232014)
  },
           "V1":{
  "vsr1":(863557214,875250014),
  "vsr2":(931035615,947023215),
  "vsr3":(965520016,999999999)
  }
           }
# method to return the science run start and stop
def getRunTimes(runName=None,ifo=None):
  """
  This method returns the start and stop of a LIGO
  run given a gpstime.  (start,stop)
  """
  if runName == None or ifo == None:
    return ()
  (a,b)=runEpochs[ifo.upper().strip()][runName.lower().strip()]
  return (a,b)

def getRunTimesAsLIGOSegment(runName=None,ifo=None):
  """
  This method returns the start and stop of a LIGO
  run given a gpstime.  (start,stop)
  """
  try:
    (a,b)=getRunTimes(runName,ifo)
  except:
    return ()
  return segments.segment(a,b)
  
def getRunEpoch(runTime=None, ifo=None):
  """
  This method returns the science run epoch name for LIGO,
  given a gps time.
  S5,S6...
  """
  if runTime == None or ifo == None:
    return ()
  outputEpoch=""
  for myEpoch in runEpochs[ifo.upper().strip()].keys():
    (start,stop)=runEpochs[ifo.upper().strip()][myEpoch]
    if (start<=runTime) and (runTime<=stop):
      outputEpoch=myEpoch
  return outputEpoch

def getKnownIfoEpochs(ifoName=None):
  """
  Simple method if given an ifo string gives a list
  of know epochs for that ifo.
  """
  if ifoName==None:
    return []
  return runEpochs[ifoName].keys()

def getKnownIfos():
  """
  This method dumps the names as a list of know ifos.
  """
  return interferometers

#Double defined home_dir funct from stfu_pipe
def home_dir():
	return os.path.split(os.environ['HOME'])[0]

########## CLASS TO WRITE LAL CACHE FROM HIPE OUTPUT #########################
class getCache(UserDict):
  """
  An instance of a lal cache
  """
  def __init__(self, options,cp, currentPath):
    UserDict.__init__(self)
    self.dir = os.listdir(os.getcwd())
    self.options = options
    self.types = ['TMPLTBANK', 'TRIGBANK', 'INSPIRAL_FIRST', \
                 'INSPIRAL_SECOND', 'THINCA_FIRST', 'THINCA_SECOND']
    self.iniNames = ['tmpltbank-path', 'trigbank-path', 'first-inspiral-path', \
         'second-inspiral-path', 'first-coinc-path', 'second-coinc-path']
    self.iniNameMaps = map(None, self.iniNames, self.types)
    self.oNames = ['bank.cache', 'trigbank.cache', 'first_inspiral.cache', \
         'second_inspiral.cache', 'first_thinca.cache', 'second_thinca.cache']
    self.nameMaps = map(None, self.oNames, self.types)
    self.ifoTypes = ['H1','H2','L1','V1','H1H2','H1L1','H2L1','H1V1', \
                          'H2V1','L1V1','H1H2L1','H1H2V1','H1L1V1', \
                          'H2L1V1','H1H2L1V1']


    if options.generate_fu_cache:
      print >> sys.stderr, "\nOption --generate-fu-cache is specified, it overwrites the hipe cache files  which already exist"
      self.getCacheAll(cp)
      self.writeCacheAll('fu_hipe.cache')
      print >> sys.stderr, "\nHIPE CACHE FILE WRITTEN TO: fu_hipe.cache"
    try:
      os.chdir("cache")
      os.chdir(currentPath)
    except:
      if len(string.strip(cp.get('followup-hipe-cache','hipe-cache-path'))) > 0:
        try:
          os.symlink(string.strip(cp.get('followup-hipe-cache','hipe-cache-path')), 'cache')
        except: print >> sys.stderr, "WARNING: cache directory exists, cannot link"
      else: pass


    # Since we are continuing get useful stuff from the ini file.
    if options.generate_fu_cache or not cp.has_option('followup-triggers','hipe-output-cache'):
      cacheString = 'fu_hipe.cache'
    else:
      cacheString = string.strip(cp.get('followup-triggers','hipe-output-cache'))
    if cp.has_option('followup-triggers','triggers-tag'):
      self.triggerTag = string.strip(cp.get('followup-triggers','triggers-tag'))
    else:
      self.triggerTag = ''

    self.cache = Cache.fromfile(open(cacheString))



  def ifoDict(self):
    return {'H1':[],'H2':[],'L1':[],'V1':[],'H1H2':[], 'H1L1':[],'H2L1':[], \
                   'H1V1':[],'H2V1':[],'L1V1':[],'H1H2L1':[],'H1H2V1':[], \
                   'H1L1V1':[],'H2L1V1':[],'H1H2L1V1':[]}

  def getCacheType(self, iniName, type, cp=None):
    self[type] = []
    p = re.compile(type)
    f = re.compile("FOLLOWUP")
    m = re.compile("-")
    x = re.compile(".xml")
    g = re.compile(".gz")
    try:
      dir = os.listdir(string.strip(cp.get('followup-hipe-cache',iniName)))
      cache_path = os.path.abspath(string.strip(cp.get('followup-hipe-cache',iniName)))
    except:
      dir = self.dir
      cache_path = os.path.abspath(self.options.cache_path)
    for fname in dir:
      if f.search(fname): continue
      if p.search(fname):
        fnamebis =  g.split(fname)[0]
        ifo = m.split(fnamebis)[0]
        tag = m.split(fnamebis)[1]
        start = m.split(fnamebis)[-2]
        dur = x.split(m.split(fnamebis)[-1])
        try:
          # scirun = string.strip(cp.get('followup-hipe-cache','science-run'))
          # tmpentry = ifo+" "+scirun+" "+start+" "+dur[0]+" "+"file://localhost" +cache_path+"/"+fname
          tmpentry = ifo+" "+tag+" "+start+" "+dur[0]+" "+"file://localhost" +cache_path+"/"+fname
          entry = lal.CacheEntry(tmpentry)
          self[type].append(entry)
        except: pass

  def getCacheAll(self,cp=None):
    for iniName, type in self.iniNameMaps:
      self.getCacheType(iniName,type,cp)


  def writeCacheType(self,oName,type,cName):
    #cName = open(oName,'w')
    for fname in self[type]:
      cName.write(str(fname)+"\n")
    #cName.close()

  def writeCacheAll(self,outputfile):
    cName = open(outputfile,'w')
    for oName, type in self.nameMaps:
      self.writeCacheType(str(oName),type,cName)
    cName.close()


  def getListFromCache(self,cache):
    listInstance = []
    for ifo in cache:
      listInstance = listInstance + cache[ifo].pfnlist()
    return listInstance

  def filesMatchingGPSinCache(self, time=None, cacheType=None, ifo_tag=None, ifo_in_coinc=None, ):
    cacheSubSet = self.ifoDict()

    if self.cache:
      for ifo in self.ifoTypes:
        try:
          if time:
            if time[ifo]:
              time_ifo = time[ifo]
            else:
              #if ifo found in the analysed times, but not in coinc...
              if ifo_tag and len(ifo)==2 and re.search(ifo,ifo_tag):
                time_ifo = time[ifo_in_coinc[0]]
              else:
                continue
            seg1 = segments.segment(time_ifo,time_ifo+1)
            seg2 = segments.segment(time_ifo-1,time_ifo)
          else:
            seg1 = None
            seg2 = None          
	  listInstance = self.cache.sieve(ifo,cacheType,None,True)
          listInstance = listInstance.sieve(None,None,seg1)
          listInstance = listInstance.sieve(None,None,seg2)
          cacheSubSet[ifo] = listInstance
        except:
          continue #print "couldn't get sub cache for " + ifo 
    else: print "no valid cache file"
    return(cacheSubSet)


  def getProcessParamsFromCache(self, opts, cp, subCache, tag, time, getInsp = False, ifo_in_coinc=None):

    process = self.ifoDict()
    inspTable = self.ifoDict()
    for ifo in subCache:
      for f in subCache[ifo]:
        # As soon as process[ifo] is filled we can stop iterating in this loop.
        if process[ifo]:
          break

        path = f.path

        # if the option "--convert-evenid" is called, a copy of the xml file
        # is made under LOCAL_XML_COPY, and the event_id is converted 
        # from int_8s to ilwd:char
        if opts.convert_eventid:
          path = self.doFileCopyAndEventIdConvert(cp,[path],True)[0]
        elif opts.create_localcopy:
          path = self.doFileCopyAndEventIdConvert(cp,[path])[0]
        doc = utils.load_filename(path)
        proc = table.get_table(doc, lsctables.ProcessParamsTable.tableName)
        if getInsp:
          insp = table.get_table(doc, lsctables.SnglInspiralTable.tableName)
        for row in proc:          
          if str(row.param).find("--ifo-tag") >= 0:
             ifoTag = row.value
          # The default is to assume that the file matches if there is no tag
          else: ifoTag = tag
          # this is a hack to handle the "-userTag" string in some xml files...
          if str(row.param).find("-userTag") >= 0:
             row.param = "--user-tag"
          # end of the hack...

        # We check that the ifo names appearing in the tag (input argument) are 
        # found in the ifoTag of the process param table. But we don't require 
        # an exact match. We allow a "tag" being double times, and an "ifoTag"
        #  being triple times. This is to handle the case of candidates which 
        # were identified when one ifo was vetoed while still in science mode.
        # Notice that this method is a hack. However implementing
        # this loose condition on the ifo-tag was required in order to be able
        # to run the followup robustly over the candidates of the 1st calendar
        # of S5 (after the IFO TIME SEPARATION).
        # THIS SHOULD BE IMPROVED IN FUTURE ANALYSES #

        ifoTagTest = True
        for j in range(0,len(tag)-1,2):
          ifostring = tag[j:j+2]
          if ifostring not in ifoTag:
            ifoTagTest = False
            break
        if ifoTagTest:

          if time[ifo]:
            time_ifo = time[ifo]
          else:
            # if time[ifo] is not defined (happens for a double coinc in triple
            #  times) then use the end_time in the first ifo of the coinc
            time_ifo = time[ifo_in_coinc[0]]
          search = table.get_table(doc, lsctables.SearchSummaryTable.tableName)
          for row in search:
            out_start_time = float(row.out_start_time)
            out_start_time_ns = float(row.out_start_time_ns)/1000000000
            out_end_time = float(row.out_end_time)
            out_end_time_ns = float(row.out_end_time_ns)/1000000000
            if ( (time_ifo >= (out_start_time+out_start_time_ns)) and (time_ifo <= (out_end_time+out_end_time_ns)) ):
              process[ifo] = proc
              if getInsp: inspTable[ifo] = insp
              break          

    if getInsp:
      return process, inspTable
    else:
      return process

  def processFollowupCache(self, cp, opts, trig, type="INSPIRAL_SECOND_"):

    if type == "INSPIRAL_SECOND_":
      # get the list of interferometers from the ifoTag
      listIfoTime = []
      for j in range(0,len(trig.ifoTag)-1,2):
        ifo = trig.ifoTag[j:j+2]
        listIfoTime.append(ifo)

      # get an ifo tag with wild cards (example: *H1*H2*)
      wildIfoTag = "*" + "*".join(listIfoTime) + "*"
      # now we define a wildIfoTag type, for example '*H1*H2*'
      # This type will be used to sieve (in filesMatchingGPSinCache())
      wildType = type + wildIfoTag

      try:
        inspiral_process_params = self.getProcessParamsFromCache( \
                       opts,cp, self.filesMatchingGPSinCache(\
                       trig.gpsTime, wildType, trig.ifoTag, trig.ifolist_in_coinc), \
                       trig.ifoTag, trig.gpsTime, False, trig.ifolist_in_coinc)
      except:
        print "couldn't get inspiral process params for " + str(trig.eventID)
        print "gpsTime:",trig.gpsTime
        print "wildType:",wildType
        print "ifoTag:",trig.ifoTag
        print "ifolist_in_coinc:",trig.ifolist_in_coinc
        inspiral_process_params = []
      return inspiral_process_params
    else:
      try:
        inspiral_process_params, sngl = self.getProcessParamsFromCache( \
                       opts,cp, self.filesMatchingGPSinCache(\
                       trig.gpsTime, type), \
                       trig.ifoTag, trig.gpsTime,True)
      except:
        print "couldn't get "+type+" process params for " + str(trig.eventID)
        inspiral_process_params = []
        sngl = []
      return inspiral_process_params, sngl

  def readTriggerFiles(self,cp,opts):
    #Patch CVT 
    numtrigs=0
    if cp.has_option('followup-triggers','sqlite-triggers'):
      sqliteTriggerFile=cp.get('followup-triggers','sqlite-triggers')
      if cp.has_option('followup-triggers','num-trigs'):
        numtrigs=cp.get('followup-triggers','num-trigs')
      else:
        numtrigs=100
      excludeTags=""
      found, coincs, search = _bandaid_(sqliteTriggerFile,numtrigs,getstatistic('far',"",""),excludeTags)
    else:
      if not self.cache or self.triggerTag == "":
        #If using multiple XML files (CVT)
        xml_glob = string.strip(cp.get('followup-triggers','xml-glob'))
        triggerList = glob.glob(xml_glob) 
      else:
        #Search for XML matching string self.triggerTag (CVT)
        triggerCache = self.filesMatchingGPSinCache(None,self.triggerTag)
        triggerList = self.getListFromCache(triggerCache)
        #Variable triggerList is a list of XML files to search(CVT)

      # if the option "--convert-eventid" is called, a copy of the xml file 
      # is made under LOCAL_XML_COPY, and the event_id is converted
      # from int_8s to ilwd:char
      if opts.convert_eventid:
        triggerList = self.doFileCopyAndEventIdConvert(cp,triggerList,True)
      elif opts.create_localcopy:
        triggerList = self.doFileCopyAndEventIdConvert(cp,triggerList)

      numtrigs = string.strip(cp.get('followup-triggers','num-trigs'))
      statistic =  string.strip(cp.get('followup-triggers','statistic'))
      bla =  string.strip(cp.get('followup-triggers','bitten-l-a'))
      blb =  string.strip(cp.get('followup-triggers','bitten-l-b'))
      #triggerList 
      if cp.has_option('followup-triggers','exclude-tags'):
        excludedTags = string.strip(cp.get('followup-triggers','exclude-tags'))
      else:
        excludedTags = None      
      found, coincs, search = readFiles(triggerList,getstatistic(statistic,bla,blb),excludedTags)
    print "COINC VARIABLE (fu_utils):",coincs,type(coincs)
    if coincs==None:
      print "COINCS object is None!\n"
    return numtrigs, found, coincs, search


  def doFileCopyAndEventIdConvert(self,cp,inputxmlfilelist,convert=False):
    newfilelist = []
    for inputxmlfile in inputxmlfilelist:
      if not os.path.isfile("LOCAL_XML_COPY/" + inputxmlfile.split('/')[-1]):
        shutil.copy(inputxmlfile,'LOCAL_XML_COPY')
        if convert:
          convert_process = call(cp.get('condor','pylal_conv_eventid') + " LOCAL_XML_COPY/" + inputxmlfile.split('/')[-1], shell=True)
          if convert_process != 0:
            print >> sys.stderr, "ligolw_conv_inspid could not be run on file " + inputxmlfile.split('/')[-1]
            sys.exit(1)
      else:
        print "The file " + inputxmlfile.split('/')[-1] + " already exist in LOCAL_XML_COPY. It will not be overwritten"
      newfilelist.append("LOCAL_XML_COPY/" + inputxmlfile.split('/')[-1])
    return newfilelist

##############################################################################
# functions to parse qscan cache files (used in analyseQscan)
##############################################################################

def getPathFromCache(fileName,type,ifo=None,time=None):
  qscanList = []
  cacheList = listFromFile(fileName)
  if len(cacheList) == 0:
    return qscanList
  for line in cacheList:
    test_line = True
    if not re.search(type,line.split(' ')[1]):
      test_line = False
    if ifo:
      if not ifo == line.split(' ')[0]:
        test_line = False
      else: pass
    if time:
      if not time == line.split(' ')[2]:
        test_line = False
      else: pass
    if test_line:
      path_output = line.split(' ')[-1]
      time_output = line.split(' ')[2]
      type_output = line.split(' ')[1]
      ifo_output = line.split(' ')[0]
      qscanList.append([path_output,time_output,type_output,ifo_output])
    else: continue
  return qscanList

##############################################################################
# class for qscan intermediate table
##############################################################################

class QscanIntermediateTable(table.Table):
  tableName = "qscan:intermediate:table"
  validcolumns = {
    "type": "lstring",
    "param": "lstring",
    "ifo": "lstring",
    "qscan_time": "lstring",
    "qscan_dir": "lstring",
  }

##############################################################################
# Function to get qscan foreground.
##############################################################################

def getForegroundTimes(cp,opts,ifo):
  times = []
  fileName = ""

  if cp.has_option('followup-triggers',ifo+'times') and string.strip(cp.get('followup-triggers',ifo+'times')):
    fileName = string.strip(cp.get('followup-triggers',ifo+'times'))
    if not os.path.exists(fileName):
      print >> sys.stderr, "File " + fileName + " does not exist"
      sys.exit(1)
    timeList = listFromFile(fileName)
    if not timeList:
      print >> sys.stderr, "the list of times in file " + fileName + " could not be found"
      sys.exit(1)
    times = stringToFloatList(timeList)

  else:
    print >> sys.stderr, "Field "+ifo+"times in [triggers] section is not specified. The ifo " + ifo + " will be ignored"

  return times,fileName

##############################################################################
# Function to get qscan background. 
##############################################################################

def getQscanBackgroundTimes(cp, opts, ifo, segFile):
    times = []
    fileName = ''
    segmentListLength = 0
  
    if cp.has_option('followup-background-qscan-times',ifo+'range'):
      rangeString = string.strip(cp.get('followup-background-qscan-times',ifo+'range'))
    else:
      rangeString = ''
    if cp.has_option('followup-background-qscan-times',ifo+'segment-list'):
      segmentListFile = string.strip(cp.get('followup-background-qscan-times',ifo+'segment-list'))
    else:
      segmentListFile = ''
    if cp.has_option('followup-background-qscan-times',ifo+'time-list'):
      timeListFile = string.strip(cp.get('followup-background-qscan-times',ifo+'time-list'))
    else:
      timeListFile = ''

    if len(rangeString) == 0 and len(segmentListFile) == 0 and len(timeListFile) == 0:
      print "No qscan background specified for " + ifo
    else:

      # Generate the list of science segments (excluding cat 1 vetoes) if a time range is provided in the ini file
      if rangeString:
        epochStart = rangeString.split(',')[0]
        epochEnd = rangeString.split(',')[1]
        #opts.gps_start_time = int(epochStart)
        #opts.gps_end_time = int(epochEnd)
        #opts.use_available_data = False
        #opts.run_data_quality = False
        # overwrite the ini file if the field "analyze" in section [segments] exist...
        #cp.set("segments", "analyze", "Science")

        #inspiralutils.findSegmentsToAnalyze(cp,opts,ifo,dq_url_pattern,segFile)
        #segmentListFile = segFile[ifo]

      # Use the segment list if provided, and generate a list of random times
      if rangeString or segmentListFile:
        segmentList = pipeline.ScienceData()
        segmentMin = cp.getint('followup-background-qscan-times','segment-min-len')
        segmentPading = cp.getint('followup-background-qscan-times','segment-pading')
        if segmentListFile:
          segmentList.read(segmentListFile,segmentMin)
        elif rangeString:
          segmentString="DMT-SCIENCE"
          if ifo.lower() == "v1":
            segmentString="ITF_SCIENCEMODE"
          segmentList = getSciSegs(ifo,int(epochStart),int(epochEnd),True,None,segmentString,segmentMin,segmentPading)
        segmentListLength = segmentList.__len__()

        if segmentListLength == 0:
          print "No Science Segments found for " + ifo

        else:          
          segmentListStart = segmentList.__getitem__(0).start()
          segmentListEnd = segmentList.__getitem__(segmentListLength - 1).end()

          seed = cp.getint('followup-background-qscan-times','random-seed')
          statistics = cp.getint('followup-background-qscan-times','background-statistics')

          random.seed(seed)
          counter = 0
          times = []

          while counter < statistics:
            gps = float(segmentListStart) + float(segmentListEnd - segmentListStart)*random.random()
            testList = copy.deepcopy(segmentList)
            secondList = [pipeline.ScienceSegment(tuple([0,int(gps),int(gps)+1,1]))]        
            if testList.intersection(secondList):
              times.append(gps)
              counter = counter + 1
            else:
              continue
          # Save the list of times in a file (for possibly re-using it later)
          timeList = floatToStringList(times)
          fileName = "timeList-" + ifo + "-" + repr(seed) + "-" + repr(statistics) + "-" + repr(segmentListStart) + "-" + repr(segmentListEnd) + ".txt"
          saveRandomTimes(timeList,fileName)

      # Use the time-list file if provided
      if segmentListLength == 0 and not len(timeListFile) == 0:
        timeList = listFromFile(timeListFile)
        fileName = timeListFile
        if not timeList:
          print >> sys.stderr, "the list of times in file " + timeListFile + " could not be found"
          sys.exit(1)
        times = stringToFloatList(timeList)

    return times, fileName

##############################################################################
# function to read/write a list of strings in a file
##############################################################################
def listFromFile(fileName):
  listInstance = []
  try:
    file = open(fileName,"r")
  except:
    print >> sys.stderr, "could not open file " + fileName
    return listInstance
  list_in_file = file.readlines()
  file.close()
  if not len(list_in_file):
    print >> sys.stderr, "No lines found in file " + fileName
    print >> sys.stderr, "Is the first line blank ?"
    return listInstance
  for line in list_in_file:
    if not line[0] == '#':
      listInstance.append(string.strip(line))
    else: pass
  return listInstance

def saveRandomTimes(timeList,fileName):
  file = open(fileName,"w")
  for time in timeList:
    file.write(time + '\n')
  file.close()

def stringToFloatList(listin):
  listout = []
  for line in listin:
    listout.append(float(line))
  return listout

def floatToStringList(listin):
  listout = []
  for line in listin:
    listout.append(repr(line))
  return listout

##############################################################################
# function to read in a list of files and extract the simInspiral tables
# and sngl_inspiral tables
##############################################################################

##############################################################################
#
# This is a section of code the "BANDAID"
#
##############################################################################

def _bandaid_(sqlFile,triggerCap=100,statistic=None,excludeTags=None):
  """
  Sqlite bandaid put in place until pipeline re-write.
  Cristina Valeria Torres: Tue-Sep-08-2009:200909081422
  hacked from the get_xml method.
  """
  sims=None
  search=None
  coincs=None
  #Create emtpy table object
  #Build by Hand
  snglInspiralTable=lsctables.SnglInspiralTable(sax.xmlreader.AttributesImpl({u"Name":lsctables.SnglInspiralTable.tableName}))
  colnamefmt = u":".join(snglInspiralTable.tableName.split(":")[:-1]) + u":%s"
  for key, value in snglInspiralTable.validcolumns.items():
    snglInspiralTable.appendChild(table.Column(sax.xmlreader.AttributesImpl({u"Name": colnamefmt % key, u"Type": value})))
  snglInspiralTable._end_of_columns()
  dataStream=table.TableStream(sax.xmlreader.AttributesImpl({u'Name':snglInspiralTable.tableName})).config(snglInspiralTable)
  #Pull in sqlite data
  try:
    sqlDB=sqlite3.connect(sqlFile)
    sqlDBsock=sqlDB.cursor()
    #Query the SQL Coincs table to get the gpstimes of the top N events
    liveTimeProgram="thinca"
    rawSeglists = db_thinca_rings.get_thinca_zero_lag_segments(sqlDB, program_name = liveTimeProgram)
    playground_segs = segmentsUtils.S2playground(rawSeglists.extent_all())
    sqlDB.create_function("is_playground", 2, lambda seconds, nanoseconds: LIGOTimeGPS(seconds, nanoseconds) in playground_segs)
    oString00=""" SELECT coinc_inspiral.end_time + coinc_inspiral.end_time_ns * 1.0e-9, coinc_event.coinc_event_id   FROM coinc_inspiral JOIN coinc_event ON (coinc_event.coinc_event_id  == coinc_inspiral.coinc_event_id) WHERE NOT  is_playground(coinc_inspiral.end_time, coinc_inspiral.end_time_ns)   AND NOT EXISTS(SELECT * FROM time_slide WHERE   time_slide.time_slide_id == coinc_event.time_slide_id AND   time_slide.offset != 0) ORDER BY combined_far LIMIT ? """%(triggerCap)
    gpsTimesToFetch=sqlDBsock.execute(oString00).fetchall()
    #Get sngl ID from coinc_map
    coincList=list()
    for item in gpsTimesToFetch:
      coincID=item[1]
      coincList.append(coincID)
    idString="SELECT event_id FROM coinc_event_map WHERE "
    for cID in coincList:
      idString=idString+""" (coinc_event_id == "%s") OR """%(cID)
    idString=idString.rstrip("OR ")
    snglEventIDs=sqlDBsock.execute(idString).fetchall()
    # Create piece of query string 
    oString01=oString02=""
    for eventInfo in snglEventIDs:
      eventID=eventInfo[0]
      oString01=oString01+""" (event_id == "%s") OR """%(eventID)
    oString01=oString01.rstrip("OR ")
    for col in snglInspiralTable.columnnames:
      oString02="%s,%s"%(str(oString02),str(col))
    oString02=oString02.lstrip(",")
    oString03="SELECT %s FROM sngl_inspiral WHERE %s"%(oString02,oString01,)
    sqlDBsock.execute(oString03)
    #Qusery coinc table using cut and paste of chads script lalapps_cbc_plotsummary
    #select %s from sngl_inspiral table where END_TIME-1 <= x+ns <= END_TIME+1
    results=sqlDBsock.fetchall()
    #Create tmp sngl table to speed things up
    oString07="""CREATE TEMPORARY TABLE sngl_inspiral_temp AS %s"""%(oString03)
    sqlDBsock.execute(oString07)
    #Grab mapping of coinc event to sngl event ids
    #Remap sngl event ids
    newResults=list()
    for elem in results:
      newID=sqlDBsock.execute("""SELECT coinc_event_id FROM coinc_event_map WHERE event_id == '%s'"""%(elem[14])).fetchall()
      newRow=list()
      for item in elem:
        newRow.append(item)
      newRow[14]=newID[0][0].replace("coinc_event:coinc_event_id:","sngl_inspiral:event_id:")
      newResults.append(newRow)
    oldResults=results
    results=newResults
    delim=dataStream.getAttribute("Delimiter")
    dataString=""
    for elem in results:
      for elem2 in elem:
        dataString="%s%s%s"%(dataString,delim,str(elem2))
    dataString=str(dataString+",").lstrip(delim)
    snglInspiralTable.appendChild(dataStream)  
    dataStream.appendData(dataString)  
    fp=open('completeSNGLINSPIRALTABLE.xml','w')
    snglInspiralTable.write(file=fp)
    fp.close()
    # SETUP SEARCH SUMMARY TABLE
    #Create a search summary table to return here
    SearchSummaryTable=lsctables.SearchSummaryTable(sax.xmlreader.AttributesImpl({u"Name":lsctables.SearchSummaryTable.tableName}))
    colnamefmt = u":".join(SearchSummaryTable.tableName.split(":")[:-1]) + u":%s"
    for key, value in SearchSummaryTable.validcolumns.items():
      SearchSummaryTable.appendChild(table.Column(sax.xmlreader.AttributesImpl({u"Name": colnamefmt % key, u"Type": value})))
      SearchSummaryTable._end_of_columns()
      SearchSummaryStream=table.TableStream(sax.xmlreader.AttributesImpl({u'Name':SearchSummaryTable.tableName})).config(SearchSummaryTable)
    oString05=""
    for col in SearchSummaryTable.columnnames:
      oString05="%s,search_summary.%s"%(str(oString05),str(col))
    oString05=oString05.lstrip(",")
    oString04=""" SELECT %s FROM search_summary LEFT OUTER JOIN sngl_inspiral_temp ON (search_summary.process_id == sngl_inspiral_temp.process_id) WHERE """%(oString05)
    for eventInfo in snglEventIDs:
      eventID=eventInfo[0]
      oString04=oString04+""" ( sngl_inspiral_temp.event_id == "%s" ) OR """%(eventID)
    oString04=oString04.rstrip("OR ")
    #Take SQL result and put into the searchsummarystream
    searchSummaryRecords=sqlDBsock.execute(oString04).fetchall()
    dataString=""
    for elem in searchSummaryRecords:
      for elem2 in elem:
                dataString="%s%s%s"%(dataString,delim,str(elem2))
    dataString=str(dataString+",").lstrip(delim)
    SearchSummaryTable.appendChild(SearchSummaryStream)
    SearchSummaryStream.appendData(dataString)
    fp2=open('completeSearchSummaryTable.xml','w')
    SearchSummaryTable.write(file=fp2)
    fp2.close()
    search=SearchSummaryTable
  except:
    sys.stderr.write("Error converting SQLITE file %s \n"%(sqlFile))
    sys.stderr.write("Generating empty sngl_inspiral table xmldoc object.\n")
    sys.stderr.flush()
    dataString=''
    coincs=snglInspiralTable
    raise
  #Mimicking method readfiles from here on out
  if statistic == None:
    coincs=snglInspiralTable
  else:
    coincInspiralTable = \
                       CoincInspiralUtils.coincInspiralTable(snglInspiralTable,statistic)
    coincs=coincInspiralTable

  return sims,coincs,search

##############################################################################
#
# End the BANDAID
#
##############################################################################

def readFiles(fileGlob,statistic=None,excludedTags=None):
  """
  read in the Sngl and SimInspiralTables from a list of files
  if Sngls are found, construct coincs, add injections (if any)
  also return Sims (if any)
  @param fileGlob: glob of input files
  @param statistic: statistic to use in creating coincs
  """
  #if fileGlob is empty return empty structures
  if not fileGlob:
    print "Warning: No glob specified, returning empty structures..."
    return None, CoincInspiralUtils.coincInspiralTable(), None

  fList = []
  for thisFile in fileGlob:
    if excludedTags:
      for thisTag in excludedTags.split(","):
        if fnmatch.fnmatch(thisFile.split('/')[-1],thisTag.strip()):
          print "WARNING: the following file will be excluded:"
          print thisFile
          continue
    fList.append(thisFile)

  if len(fList) == 0:
    print "Warning: After removing forbidden tags, no remaining files in glob. Returning empty structures..."
    return None, CoincInspiralUtils.coincInspiralTable(), None

  sims = None
  coincs = None
  search = None
  for thisFile in fList:
    doc = utils.load_filename(thisFile)
    # extract the sim inspiral table
    try:
      simInspiralTable = \
          table.get_table(doc, lsctables.SimInspiralTable.tableName)
      if sims: sims.extend(simInspiralTable)
      else: sims = simInspiralTable
    except: simInspiralTable = None

    # extract the sngl inspiral table, construct coincs
    try: snglInspiralTable = \
      table.get_table(doc, lsctables.SnglInspiralTable.tableName)
    except: 
      snglInspiralTable = None
      searchSumTable = None
    if snglInspiralTable:
      coincInspiralTable = \
        CoincInspiralUtils.coincInspiralTable(snglInspiralTable,statistic)
      if simInspiralTable:
        coincInspiralTable.add_sim_inspirals(simInspiralTable)
      # extract the search_summary table only if a sngl inspiral table is found
      searchSumTable = table.get_table(doc,lsctables.SearchSummaryTable.tableName)
      if coincs: 
        coincs.extend(coincInspiralTable)
        search.extend(searchSumTable)         
      else: 
        coincs = coincInspiralTable
        search = searchSumTable    
  print "processed " + str(len(fList)) + " files..." 
  #for ck in coincs:
  #  print getattr(ck, 'event_id'), getattr(ck, 'numifos')
  return sims,coincs,search


#############################################################################
# function to set up directories
#############################################################################

def createdir(newDir):
  try:
    os.mkdir(newDir)
  except:
    pass


def setupdirs():

  try:
    os.chdir("clustered")
    os.chdir('..')
  except: os.mkdir("clustered")

  try:
    os.chdir("found")
    os.chdir('..')
  except: os.mkdir("found")

  try:
    os.chdir("missed")
    os.chdir('..')
  except: os.mkdir("missed")

  try:
    os.chdir("followuptrigs")
    os.chdir('..')
  except: os.mkdir("followuptrigs")

  try:
    os.chdir("followupfound")
    os.chdir('..')
  except: os.mkdir("followupfound")

  try:
    os.chdir("followupmissed")
    os.chdir('..')
  except: os.mkdir("followupmissed")

  try:
    os.chdir("logs")
    os.chdir('..')
  except: os.mkdir("logs")

  try:
    os.chdir("datafind_cache")
    os.chdir('..')
  except: os.mkdir("datafind_cache")

#############################################################################
# function to return the number of slides in a file (as a string)
#############################################################################
def getslidenum(fName):
  command = "grep -m 1 -e '--num-slides' " + fName + " | sed -e 's@\"@\'@g"
  fStr = os.popen(command).readlines()
  if fStr:
    fStrlist = fStr[0].rsplit(",")
    return fStrlist[-2]
  else:
   return "0"

#############################################################################
# function to glob files
#############################################################################
def globxmlandinj(opts):

  #Glob for the files both injections and trigger files
  fList = glob.glob(opts.xml_glob)
  if not fList:
    print >>sys.stderr, "The glob for " + opts.xml_glob + " returned no files"
    sys.exit(1)

  if opts.injection_glob:
    iList = glob.glob(opts.injection_glob)
    if not iList:
      print >> sys.stderr, "The glob for " + opts.injection_glob + " returned no files"
      sys.exit(1)
  else:
    iList = None

  return fList,iList

##############################################################################
# function to coire the input files
##############################################################################
def coire(opts,fList,iList=None):
  """
  do coires
  """
  sims = None
  coincs = None

  # we need to get slide num and data type from the files
  # this is just hard coded now!!!
  for xmls in fList:
    if opts.verbose:
      print "running lalapps_coire on " + xmls
    numslides = getslidenum(xmls)
    command = "lalapps_coire " + \
              "--glob '" + xmls + "'" + \
              " --data-type " + str(opts.data_type) + \
              " --num-slides " + numslides + " --coinc-stat " + \
              opts.statistic + " --cluster-time " + \
              str(1000*opts.cluster_time)

    if iList:
      command += " --injection-file " + iList[0]
      command += " --injection-window " + str(opts.injection_coinc_time)
      command += " --missed-injections missed/MISSED" + xmls
      command += " --output found/FOUND" + xmls
    else:
      command += " --output clustered/CLUSTER"
      command += str(opts.cluster_time) + "s_" + xmls
    if opts.bitten_l_a and opts.bitten_l_b:
      command += " --h1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h2-bittenl-a " + str(opts.bitten_l_a)
      command += " --l1-bittenl-a " + str(opts.bitten_l_a)
      command += " --h1-bittenl-b " + str(opts.bitten_l_b)
      command += " --h2-bittenl-b " + str(opts.bitten_l_b)
      command += " --l1-bittenl-b " + str(opts.bitten_l_b)

    os.system(command)

##############################################################################
# function to extract the statistic information
##############################################################################
def getstatistic(stat, bla, blb):

  if stat == "effective_snrsq":
    newstat = "effective_snr"
  else:
    newstat = stat

  statistic=CoincInspiralUtils.coincStatistic( newstat, bla, blb )
  return statistic

##############################################################################
# function to query segment server looking for science segments
##############################################################################
def getSciSegs(ifo=None,
               gpsStart=None,
               gpsStop=None,
               cut=bool(False),
               serverURL=None,
               segName="DMT-SCIENCE",
               seglenmin=None,
               segpading=0
):
  """
  This method is designed to query the server specified by SERVERURL
  if not specified the method will use the environment variable
  S6_SEGMENT_SERVER to determine who to query.
  The method will return the segments that are between and overlaping
  with the variable gpsStart and gpsStop.  If the flag cut is
  specified to be True then the returned science segment list will be
  cut so that the 
  times are between gpsStart and gpsStop inclusive.  In addition to
  these required arguments you must also specify in a text string the
  IFO of interest.  Valid entries are L1 H1 V1 , but only one IFO at a
  time can be specified.  You can call this method by specifying
  specific keyswords ifo,gpsStart,gpsStop,cut,serverURL.  For example
  to call using no segment cuts and the default URL try:
  x=getSciSegs(gpsStart=987654321,gpsStop=876543210)
  A query failure will give an error but no records found for the
  options specified will return an empty glue.pipeline.ScienceData()
  segment list.  We restrict the max version query to operate only
  over the segment interval requested not all know segments!
  Returns a data structure of type glue.pipeine.ScienceData()
  """
  if sum([x==None for x in (ifo,gpsStart,gpsStop)])>0:
    sys.stderr.write("Invalid arguments given to getSciSegs.\n")
    return None
  ifo=ifo.strip()
  query01 ="""SELECT \
  segment.start_time, \
  segment.end_time \
  FROM segment, segment_definer \
  WHERE \
  segment.segment_def_id  = segment_definer.segment_def_id AND \
  segment.segment_def_cdb = segment_definer.creator_db AND \
  segment_definer.name = '%s' AND \
  segment_definer.ifos = '%s' AND \
  NOT (segment.start_time > %s OR  %s > segment.end_time) AND \
  segment_definer.version = (SELECT MAX(sd.version) FROM \
  segment_definer AS sd, segment_summary AS ss WHERE \
  sd.name = segment_definer.name AND \
  sd.segment_def_id = ss.segment_def_id AND \
  NOT (ss.start_time > %s OR %s > ss.end_time)) \
  ORDER BY segment.start_time \
  """
  #Determine who to query if not specified.
  if serverURL == None:
    serverURL=os.getenv('S6_SEGMENT_SERVER')
    if serverURL == None:
      serverURL=defaultsegmentserver
  try:
    connection=None
    connection=segmentdb_utils.setup_database(serverURL)
  except Exception, errMsg:
    sys.stderr.write("Error connection to %s\n"\
                     %(serverURL))
    sys.stderr.write("Error Message :\t %s \n"%(errMsg))
    return None
  try:
    sqlQuery=query01%(segName,ifo,gpsStop,gpsStart,gpsStop,gpsStart)
    engine=query_engine.LdbdQueryEngine(connection)
    queryResult=engine.query(sqlQuery)
  except Exception, errMsg:
    sys.stderr.write("SciSeg query failed %s\n"%(serverURL))
    sys.stdout.write("Error fetching sci segs %s : %s\n"%(gpsStart,gpsStop))
    sys.stderr.write("Error message seen: %s\n"%(str(errMsg)))
    sys.stderr.write("Query Tried: \n %s \n"%(sqlQuery))
    return
  engine.close()
  queryResult.sort()
  #Take segment information and turn into
  #ScienceData() object
  segListTemp = pipeline.ScienceData()
  #Append the raw data into the ScienceData class()
  segIndex=0
  for rawStart,rawStop in queryResult:
    if cut:
      if int(rawStart)<int(gpsStart):
        rawStart=gpsStart
      if int(rawStop)>int(gpsStop):
        rawStop=gpsStop
    segListTemp.append_from_tuple((segIndex,rawStart,rawStop,rawStop-rawStart))
    segIndex=+1
    segListTemp.coalesce()
  if not seglenmin: return segListTemp
  else:
    if segpading and 2*segpading >= seglenmin:
      sys.stderr.write("segpading must be smaller than seglenmin/2\n")
      sys.exit(1)
    segList = pipeline.ScienceData()
    for indice in range(0,segListTemp.__len__()):
      segTemp = segListTemp.__getitem__(indice)
      if segTemp.dur() >= seglenmin:
        segList.append_from_tuple((segTemp.id(),segTemp.start()+segpading,segTemp.end()-segpading,segTemp.dur()-2*segpading))
    return segList
#End getSciSegs()
#

#############################################################################
# Follow up list class definition
#############################################################################
class followUpList:
  """
  Class to hold gps and ifo pairs to send to subsequent functions
  It also holds an instance of the coinc class which contains
  All of the relevant xml information
  """
  def __init__(self,Coincs = None, Missed = None ):
    self.gpsTime = {"H1" : None, "H2" : None, "L1" : None,
                  "G1" : None, "V1" : None, "T1" : None}
    self.ifos = '' # added to construct a new field in "followups", see function getfollowuptrigs
    self.ifolist_in_coinc = None
    self.ifoTag = None
    self.coincs = Coincs
    self.missed = Missed
    self.eventID = None
    self.statValue = None
    self.far = None
    self.page = None
    self.summarydir = None
    self.summarypage = None
    self.types = ['H1H2', 'H1H2L1', 'H1H2L1V1', 'H1H2V1', 'H1L1', 'H1L1V1', 'H1V1', 'H2L1', 'H2L1V1', 'H2V1', 'L1V1']
    self.rank = -1
    #self.magic_number = 250.0
    self.magic_number = None
  def get_coinc_type(self):
    ifostring = ""
    for t in ['H1','H2','L1','V1']:
      if self.gpsTime[t]: ifostring+=t
    return ifostring

  def add_coincs(self,Coincs):
    setattr(self,"coincs",Coincs)
    # convert the ilwdchar "Coincs.event_id" into an integer. This removes
    # the extra table_name and column_name and keeps only the ID string
    self.eventID = int(Coincs.event_id)
    self.statValue = Coincs.stat
    self.ifos, self.ifolist_in_coinc = Coincs.get_ifos()
    if self.is_trigs():
      self.summarydir = "followuptrigs"
    if self.is_found():
      self.summarydir = "followupfound"
  
  def add_missed(self,Missed):
    setattr(self,"missed",Missed)
    self.summarydir = "followupmissed"
  
  def is_trigs(self):
    if isinstance(self.coincs,CoincInspiralUtils.coincInspiralTable.row):
      return 1
  
  def is_found(self):
    sim = None
    try: 
      sim = isinstance(self.coincs.sim,lsctables.SimInspiral)
    except: return 0
    if sim:
      return 1
  def is_missed(self):
    if isinstance(self.missed,lsctables.SimInspiralTable):
      return 1
  def add_page(self,page):
    self.page = page

  def write_trigger_info(self, fobj):
    fobj.write("Rank:"+str(self.rank)+",ID:"+str(self.eventID)+",Stat:"+str(self.statValue)+",FAR:"+str(self.far)+",Type:"+str(self.get_coinc_type())+",IfoTime:"+str(self.ifoTag)+",")
    for ifo in ['H1','H2','L1','V1']:
      if self.gpsTime[ifo]:
        if self.magic_number:
          eff_snr = getattr(self.coincs,ifo).get_effective_snr(self.magic_number)
        else:
          eff_snr = getattr(self.coincs,ifo).get_effective_snr()
        fobj.write(ifo+"time:"+repr(self.gpsTime[ifo])+","+ifo+"mchirp:"+str(getattr(self.coincs,ifo).mchirp)+","+ifo+"eta:"+str(getattr(self.coincs,ifo).eta)+","+ifo+"mass1:"+str(getattr(self.coincs,ifo).mass1)+","+ifo+"mass2:"+str(getattr(self.coincs,ifo).mass2)+","+ifo+"snr:"+str(getattr(self.coincs,ifo).snr)+","+ifo+"chisq:"+str(getattr(self.coincs,ifo).chisq)+","+ifo+"chisq_dof:"+str(getattr(self.coincs,ifo).chisq_dof)+","+ifo+"duration:"+str(getattr(self.coincs,ifo).template_duration)+","+ifo+"eff_snr:"+str(eff_snr)+",")
    fobj.write("\n")

    

#############################################################################
# Function to generate a trigbank xml file
#############################################################################
def generateXMLfile(cp,ckey,ifo,outputPath=None,type='plot',use_max_template=None,use_col_from_installation=None):

  if outputPath:
    try:
      os.mkdir(outputPath) 
    except: pass

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
 
  # if we want our trigbank files to use the loudest template
  maxSNR = 0
  if use_max_template:
    maxSNR = 0
    maxIFO = ""
    for t in ckey:
      snr = t.snr
      if snr > maxSNR:
        maxSNR = snr
        maxIFO = t.ifo
    trig = getattr(ckey,maxIFO)
  else:
    trig = getattr(ckey,ifo)

  trigger = copy.deepcopy(trig)

  if type == "coh":
    if cp.has_option("followup-coh-trigbank",ifo+"_channel"):
      trigger.channel = cp.get("followup-coh-trigbank",ifo+"_channel")
    else:
      print >> sys.stderr, "the section [followup-coh-trigbank] in the .ini file should contain a field \""+ifo+"_channel\""
      sys.exit(1)
  else:
    properChannelTrig = getattr(ckey,ifo)
    trigger.channel = properChannelTrig.channel
  trigger.ifo = ifo
  #print ifo, getattr(ckey,ifo).ifo
  #print trig.channel, trig.ifo, maxSNR, getattr(ckey, 'event_id')
  # BEFORE WE MAKE A NEW TABLE FIGURE OUT WHAT COLUMNS ARE VALID !!!
  valid_columns = trigger.__slots__
  columns = []
  notcolumns = []
  for col in valid_columns:
    try: 
      getattr(trigger,col)
      columns.append(col)
    except:
      notcolumns.append(col)
  # IF "use_col_from_installation" IS TRUE, ADD NEW COLUMNS TO THE SNGL_INSPIRAL TABLE
  if use_col_from_installation:
    for col in notcolumns:
      print "\n adding column " + col
      columns.append(col)
      setattr(trigger,col,0)

  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table) 

  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable,columns)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
  sngl_inspiral_table.append(trigger)

  fileName = ifo + '-TRIGBANK_FOLLOWUP_' + type + str(int(ckey.event_id)) + ".xml.gz"
  if outputPath:
    fileName = outputPath + '/' + fileName
  utils.write_filename(xmldoc, fileName, verbose = False, gz = True)   

def generateBankVetoBank(fuTrig, ifo,sngl,subBankSize,outputPath=None):

  sngl_copy = copy.deepcopy(sngl)
  trig =  getattr(fuTrig.coincs,ifo)
  mass = trig.mass1 + trig.mass2
  if outputPath:
    try:
      os.mkdir(outputPath)  
    except: pass
  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table)
  valid_columns = lsctables.SnglInspiralTable.validcolumns
  columns = []
  notcolumns = []
  for col in valid_columns:
    try:
      getattr(sngl_copy[0],col)
      columns.append(col)
    except:
      notcolumns.append(col) 
  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
  sngl_copy.sort(lambda a, b: cmp(a.mtotal, b.mtotal))
  index = sngl_copy.getColumnByName('mtotal').index(mass)
  fromEnd = len(sngl_copy)-index-int(subBankSize/2)
  if fromEnd < 0:
    sngl_sub = sngl_copy[index+fromEnd-int(subBankSize/2):-1]
  else: 
    sngl_sub = sngl_copy[index-int(subBankSize/2):index+int(subBankSize/2)]
 
  for row in sngl_sub:
    for col in notcolumns:
      setattr(row,col,0)
    sngl_inspiral_table.append(row)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
 
  fileName = ifo + '-BANKVETO_FOLLOWUP_' + str(fuTrig.eventID) + ".xml.gz"
  if outputPath:
    fileName = outputPath + '/' + fileName
  utils.write_filename(xmldoc, fileName, verbose = True, gz = True)
  return fileName

#############################################################################
# Function to return the follow up list of coinc triggers
#############################################################################
def getfollowuptrigs(cp,numtrigs,trigtype=None,page=None,coincs=None,missed=None,search=None,trigbank_test=None,ifar=True,add_columns=False):

  followups = []
  xmfilenum = 0
  trigtime = 0
  cnt = 0
  # get segments to followup
  seglistname = string.strip(cp.get('followup-triggers','segment-list'))
  if cp.has_option('followup-triggers','eff-snr-denom-fac'):
    magic_number = float(string.strip(cp.get('followup-triggers','eff-snr-denom-fac')))
  else:
    magic_number = None
  if seglistname: 
    seglist = segmentsUtils.fromsegwizard(open(seglistname,'r'))
    seglist.coalesce()
  else: seglist = []
  if seglist: print "WARNING: restricting triggers to specified segment list"

  # print out the trigger info
  #trigInfo = open("trigger_info.txt","w")

 
  if coincs:
      if not ifar:
          sim = None
          try:
              sim = isinstance(coincs[0].sim,lsctables.SimInspiral)
          except: pass
          if sim:
	      print "following up injections..."
              coincs.sort() # This does an descending sort even for found inj !!! CHANGE THIS? It USED TO BE ASCENDING
          else:
              coincs.sort()
          numTrigs = 0
          # the loop over coincident triggers
          for ckey in coincs:
              cnt += 1
              fuList = followUpList()
              fuList.add_coincs(ckey)
              if page:
                  fuList.add_page(page)
              ifo_list = ['H1','H2','L1','G1','V1','T1']
              for ifo in ifo_list:
                  try:
                      getattr(ckey,ifo)
		      #print getattr(ckey, 'numifos'), getattr(ckey, 'event_id')
                      fuList.gpsTime[ifo] = (float(getattr(ckey,ifo).end_time_ns)/1000000000)+float(getattr(ckey,ifo).end_time)
                      trigtime = fuList.gpsTime[ifo]
                  except: fuList.gpsTime[ifo] = None
        
	      # if the trigger type is wrong don't count it: continue
              if trigtype: 
	        stop = 0
	        for t in trigtype.split(','):
		  if t.strip() == fuList.get_coinc_type(): stop = 1
		if not stop: continue  

	      # if trigger is not in requested time move on 
              if seglist:
                if (trigtime in seglist): print "trigger " + str(trigtime) + " in seglist ranked " + str(cnt)
                else: continue

              numTrigs += 1
              if numTrigs > int(numtrigs):
                break
              # write the trig bank files	      
              for ifo in ifo_list:
                if fuList.gpsTime[ifo] and trigbank_test:
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"plot",None,add_columns)
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"notrig",None,add_columns)
                      # Also make a trigbank xml with the max template
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"coh",True,add_columns)
                  xmfilenum += 3

              # now, find the ifoTag associated with the triggers, 
              # using the search summary tables...
              if fuList.ifolist_in_coinc:
                  firstIfo = fuList.ifolist_in_coinc[0]
                  triggerTime = fuList.gpsTime[firstIfo]
                  if search:  
                      for chunk in search:
                          out_start_time = float(chunk.out_start_time)
                          out_start_time_ns = float(chunk.out_start_time_ns)/1000000000
                          out_end_time = float(chunk.out_end_time)
                          out_end_time_ns = float(chunk.out_end_time_ns)/1000000000
                          if ( (triggerTime >= (out_start_time+out_start_time_ns)) and (triggerTime <= (out_end_time+out_end_time_ns)) ):
                              fuList.ifoTag = chunk.ifos
                              break
              fuList.rank = cnt
              fuList.magic_number = magic_number
              followups.append(fuList)
              #write info to file
              #fuList.write_trigger_info(trigInfo,cnt)
              if search:
                # generate a cohbank xml file for this coinc trigger
                generateCohbankXMLfile(ckey,fuList.gpsTime[firstIfo],fuList.ifoTag,fuList.ifolist_in_coinc,search,'trigTemplateBank',"coh",add_columns)
                # generate a trigbank with max templates for ifo in science but not found in coincidence
                if trigbank_test:
                  for j in range (0,len(fuList.ifoTag)-1,2):
                    itf = fuList.ifoTag[j:j+2]
                    if not itf in fuList.ifolist_in_coinc:
                      generateXMLfile(cp,ckey,itf,'trigTemplateBank',"coh",True,add_columns)

      else:

          ifarList=list()
          # the loop over coincident triggers with ifar
          for ckey in coincs:
              cnt += 1
              myIFAR = getattr(ckey,ckey.get_ifos()[1][0]).alpha
              myKey = ckey
              ifarList.append([myIFAR,myKey])
          ifarList.sort()
          ckeyList=list()
          for key in ifarList:
              ckeyList.append(key[1])
          numTrigs = 0
          cnt = 0
          for ckey in ckeyList:
              cnt += 1

              fuList = followUpList()
              fuList.far = getattr(ckey,ckey.get_ifos()[1][0]).alpha
              fuList.add_coincs(ckey)
              if page:
                  fuList.add_page(page)
              ifo_list = ['H1','H2','L1','G1','V1','T1']
              for ifo in ifo_list:
                  try:
                      getattr(ckey,ifo)
                      fuList.gpsTime[ifo] = (float(getattr(ckey,ifo).end_time_ns)/1000000000)+float(getattr(ckey,ifo).end_time)
                      trigtime = fuList.gpsTime[ifo]
                  except: fuList.gpsTime[ifo] = None

              # if the trigger type is wrong don't count it
              if trigtype:
                stop = 0
                for t in trigtype.split(','):
                  if t.strip() == fuList.get_coinc_type(): stop = 1
                if not stop: continue
              # if trigger is not in requested time move on 
              if seglist:
                if (trigtime in seglist): print "trigger " + str(trigtime) + " in seglist ranked " + str(cnt)
                else: continue

              numTrigs += 1
              if numTrigs > int(numtrigs):
                break

              # make the trig bank files
              for ifo in ifo_list:
                if fuList.gpsTime[ifo] and trigbank_test:
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"plot",None,add_columns)
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"notrig",None,add_columns)
                  # Also make a trigbank xml with the max template
                  generateXMLfile(cp,ckey,ifo,'trigTemplateBank',"coh",True,add_columns)
                  xmfilenum += 3

              # now, find the ifoTag associated with the triggers, 
              # using the search summary tables...
              if fuList.ifolist_in_coinc:
                  firstIfo = fuList.ifolist_in_coinc[0]
                  triggerTime = fuList.gpsTime[firstIfo]
                  if search:  
                      for chunk in search:
                          out_start_time = float(chunk.out_start_time)
                          out_start_time_ns = float(chunk.out_start_time_ns)/1000000000
                          out_end_time = float(chunk.out_end_time)
                          out_end_time_ns = float(chunk.out_end_time_ns)/1000000000
                          if ( (triggerTime >= (out_start_time+out_start_time_ns)) and (triggerTime <= (out_end_time+out_end_time_ns)) ):
                              fuList.ifoTag = chunk.ifos
                              break 
              fuList.rank = cnt
              fuList.magic_number = magic_number
              followups.append(fuList)
              
              #write info to file
              #fuList.write_trigger_info(trigInfo,cnt)
              # generate a cohbank xml file for this coinc trigger
              if search:
                generateCohbankXMLfile(ckey,fuList.gpsTime[firstIfo],fuList.ifoTag,fuList.ifolist_in_coinc,search,'trigTemplateBank',"coh",add_columns)
                if trigbank_test:
                  for j in range (0,len(fuList.ifoTag)-1,2):
                    itf = fuList.ifoTag[j:j+2]
                    if not itf in fuList.ifolist_in_coinc:
                      generateXMLfile(cp,ckey,itf,'trigTemplateBank',"coh",True,add_columns)


  # the missed stuff doesnt work yet!!!
  print "produced " + str(xmfilenum) + " trig bank files..."

  # compute the trigger types found
  typesdict = {"H1H2" : 0, "H1L1" : 0, "H1V1" : 0, "H2L1" : 0, "H2V1" : 0, "L1V1" : 0, "H1H2L1" : 0, "H1H2V1" : 0, "H2L1V1" : 0, "H1L1V1" : 0, "H1H2L1V1" : 0}
  for f in followups:
    typesdict[f.get_coinc_type()] += 1
  for f in typesdict:
    print "found " + str(typesdict[f]) + " " + str(f) +  " triggers"
  #trigInfo.close()
  if missed:
    followups
  return followups


#############################################################################
# Method to write a table of inspiral parameters
#############################################################################
def writeParamTable(trigger,opts):
  page = markup.page(mode="strict_html")
  page._escape = False
  #doctype="""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">"""
  #doctype+="""\n<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">"""

  title = "Inspiral Parameters for trigger " + str(trigger.eventID)
  #page.init(title=title, doctype=doctype) 
  page.init(title=title)

  paramList = [['SNR','snr'],['CHISQ','chisq'],['Chirp Mass','mchirp'],['Eta','eta'],['Mass 1','mass1'],['Mass 2','mass2'],['Eff Dist (Mpc)', 'eff_distance']]

  page.add('<table bgcolor=cyan border=1px>')
  page.tr();
  page.td('Ifo');
  page.td('End Time')
  for param in paramList:
    page.td("<b>"+param[0]+"</b>");
  page.td('<b>Combined Stat</b>');
  if not opts.disable_ifarsorting:
    page.td('<b>FAR</b>');
  page.tr.close()
  for ifo in trigger.ifolist_in_coinc:
    page.tr();
    page.td(ifo);
    page.td(repr(trigger.gpsTime[ifo]));
    for param in paramList:
      page.td("%0.3f"%(eval("getattr(trigger.coincs,ifo)."+param[1]))); 
    page.td("%0.5f"%(trigger.statValue));
    if not opts.disable_ifarsorting:
      page.td("%0.4f"%(trigger.far))
    page.tr.close()
  page.table.close()

  if not os.access('PARAM_TABLES',os.F_OK):
    os.mkdir('PARAM_TABLES')
  else: pass
  html_file = file("PARAM_TABLES/table_" + str(trigger.eventID) + ".txt","w")
  html_file.write(page(False))
  html_file.close()

#############################################################################
# Class to hold summary HTML information for all of the functions
#############################################################################
class summaryHTMLTable:

  def __init__(self,trig):
    if trig.is_trigs() and not trig.is_found():
      self.summarypath =  "followuptrigs/"
    if trig.is_trigs() and trig.is_found():
      self.summarypath = "followupfound/"
    self.eventID = trig.eventID
    self.statValue = trig.statValue
    self.H1time = trig.gpsTime["H1"]
    self.H2time = trig.gpsTime["H2"]
    self.L1time = trig.gpsTime["L1"]
    self.G1time = trig.gpsTime["G1"]
    self.V1time = trig.gpsTime["V1"]
    self.T1time = trig.gpsTime["T1"]
    self.containers = []

class HTMLcontainer:

  def __init__(self,trig,name,alt_web=None):
    # The missed injections dont work yet!!!
    self.name = name.rsplit(".")[-1]
    self.detailpath = ""
    if trig.is_trigs() and not trig.is_found():   
      os.chdir("followuptrigs")
      try: 
        os.chdir(self.name)
        os.chdir("../../")
      except: 
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = trig.page + "/followuptrigs/" + self.name + "/"
      self.localdetailpath = "followuptrigs/" + self.name + "/"

    if trig.is_trigs() and trig.is_found():
      os.chdir("followupfound")
      try:
        os.chdir(self.name)
        os.chdir("../../")
      except:
        os.mkdir(self.name)
        os.chdir("../")
      self.detailpath = trig.page + "/followupfound/" + self.name + "/"
      self.localdetailpath = "followupfound/" + self.name + "/"

    self.image = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".png"
    self.localimage = self.localdetailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + self.name + ".png" 
    self.text = "click here"
    if alt_web: name = alt_web
    else: name = self.name
    self.link = self.detailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + name + ".html"
    self.locallink = self.localdetailpath + str(trig.statValue) + "_" + str(trig.eventID) + "_" + name + ".html"



##############################################################################
# Function to write the HTML tables to pages
##############################################################################
def writeIULHeader(file):
  file.write('<%method title>Follow Up Report</%method><%method headline>Follow Up Report</%method><%method cvsid>' + git_version.id + '</%method>\n')

def beginSummaryTable(file, table):
  file.write("<h3>Trigger [" + str(table.eventID) + 
             "] with combined statistic = " + str(table.statValue) + "</h3>")
  file.write('\n<br><table width=800 border=1>')
  
def endSummaryTable(file, table):
  file.write('\n</table>')

def writeModule(file, container):
  file.write('\n<tr><td width=400><font color="red" size=5>'+
             container.name + '</font>')
  file.write('\n<br><a href= "' + container.link + '">' + container.text 
             + "</a></td>")
  file.write('\n<td><a href= "' + container.image + '">\n<img src="' + 
              container.image + '" width=400 alt="No Image"></a></td></tr>')

def writeHTMLTables(summaryHTMLlist):
  for table in summaryHTMLlist:
    tableFile = open(table.summarypath + str(table.statValue) + "_" + 
                      str(table.eventID) + "_summary.html" ,'w')
    beginSummaryTable(tableFile,table)
    for container in table.containers:
      writeModule(tableFile, container)
    endSummaryTable(tableFile,table)
    writeIULHeader(tableFile)
    tableFile.close()

class HTMLTable:
  def __init__(self):
    self.columns = []
    self.headers = []

  def add_column(self,rows,header):
    self.columns.append(rows)
    self.headers.append(header)

  def write(self,file):
    file.write('\n<br><table><tr>')
    cnt = 0
    for i in self.columns:
      file.write('<td><b>'+ self.headers[cnt]+'</b><table>')
      cnt +=1
      for j in i:
        file.write('<tr><td>'+str(j)+'</td></tr>\n')
      file.write('</table></td>')
    file.write('\n</tr></table><br>\n')

##############################################################################
# Function to publish the web tree
##############################################################################
def publishOnHydra(page):
  indexFile = open("index.html","w")
  patt = re.compile('.*summary.html')
  
  # First do the found injections
  files = os.chdir('followupfound')
  files = os.listdir('.')
  table = HTMLTable()
  fList = []
  stat = []
  ID = []
  if files:
    for f in files:
      temp = patt.match(f)
      if temp:
        fList.append((float(temp.group().rsplit('_')[0]), \
                  temp.group().rsplit('_')[1]))
    indexFile.write('<h3>Found Injection follow ups</h3>\n')
    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=False)
    for i in sortedF:
      stat.append(str(i[0]))
      ID.append('<a href="' +page+ '/followupfound/'+str(i[0])+'_'+str(i[1])+
                '_summary.html">'+i[1]+'</a>')
    table.add_column(stat,'Stat Value')
    table.add_column(ID,'ID')
    table.write(indexFile)  
  os.chdir('..')

# then do trigs
  files = os.chdir('followuptrigs')
  files = os.listdir('.')
  table = HTMLTable()
  fList = []
  stat = []
  ID = []
  if files:
    for f in files:
      temp = patt.match(f)
      if temp:
        fList.append((float(temp.group().rsplit('_')[0]), \
                  temp.group().rsplit('_')[1]))
    indexFile.write('<h3>Found Trigger follow ups</h3>\n')
    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=True)
    for i in sortedF:
      stat.append(str(i[0]))
      ID.append('<a href="' +page+ '/followuptrigs/'+str(i[0])+'_'+str(i[1])+
                '_summary.html">'+i[1]+'</a>')
    table.add_column(stat,'Stat Value')
    table.add_column(ID,'ID')
    table.write(indexFile)
  os.chdir('..')



#  # Then do the triggers 
#  files = os.chdir('followuptrigs')
#  files = os.listdir('.')
#  table = HTMLTable()
#  fList = []
#  stat = []
#  ID = []
#  if files:
#    for f in files:
#      temp = patt.match(f)
#      if temp:
#        fList.append((float(temp.group().rsplit('_')[0]), \
#                  temp.group().rsplit('_')[1]))
#    indexFile.write('<h3>Trigger follow ups</h3>\n')
#    sortedF = sorted(fList,key=operator.itemgetter(0),reverse=True)
#    for i in sortedF:
#      stat.append(str(i[0]))
#      ID.append('<a href="' +page+ '/followuptrigs/'+str(i[0])+'_'+str(i[1])+
#                '_summary.html">'+i[1]+'</a>')    
#    table.add_column(stat,'Stat Value')
#    table.add_column(ID,'ID')
#    table.write(indexFile)
#  
#  os.chdir('..')

  writeIULHeader(indexFile)
  indexFile.close()
  #This needs to be done so that it doesn't keep writing over this 
  #directory... ;)
  os.system('scp -r index.html followuptrigs followupfound followupmissed ' +
            'hydra.phys.uwm.edu:/home/htdocs/uwmlsc/root/'+
            page + '.')

############################################################
# Class for checking follow up GPS times to Nelson's veto lists.
############################################################
class nVeto:
  """
  Author: Cristina Valeria Torres
  Provides a method for checking S5y1 and S5y2 follow up candidates
  against locally cached static veto lists produced by Nelson Christensen
  """
  def __init__(self):
    """
    Initialize this class.
    """
    self.database=list([["L1","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/L1/VetoLists",list()],
                        ["H1","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/H1/VetoLists",list()],
                        ["H2","/archive/home/romain/Projects/LowMassCBC/20051104-20061114/triggers/upperlimits_v99_july08/nelsonVetoes/H2/VetoLists",list()]])
      #Form [[filename,list(gpsIntervals)],[filename2,list(gpsIntervals)]...
    self.vetoExtension=".veto"
    self.tolWin=float(0.01)
    self.filesLoaded=bool(False)
  #End self.__init__()

  def __loadFiles__(self):
    """
    This method is invoked automatically if required.
    It creates the dynamically structures to hold the veto lists.
    """
    #Check that the paths specified exist

    for ifoName,ifoPath,vetoData in self.database:
        if not os.path.exists(ifoPath):
          raise IOError("%s Path Invalid: %s"%(ifoName,ifoPath))
    for index in range(self.database.__len__()):
      self.__processPath__(index)
    self.filesLoaded=True
  #End self.__loadFiles__()

  def __processPath__(self,index=""):
    """
    This method is not to be called explicity.  It is called by
    self.__loadFiles__(). The role of this method is to fill up 
    the instance with the data required to check the veto intervals.
    """
    path=self.database[index][1]
    listOfFiles=os.listdir(path)
    listOfVetoFiles=list()
    for file in listOfFiles:
      if file.endswith(self.vetoExtension):
        listOfVetoFiles.append(file)
    listOfFiles=[os.path.normpath(path+"/"+file) for file in listOfVetoFiles]
    for filename in listOfFiles:
      fileData=list()
      fp=open(str(filename),'r')
      txtData=fp.readlines()
      fp.close()
      for row in txtData:
        a,b,c=map(float,row.split())
        fileData.append([a,b,c])
      self.database[index][2].append([filename,fileData])
  #End self.__processPath__()

  def __checkIntervals__(self,vList,gpsTime):
    """
    This method is not to be called explicity. It is called by
    self.findInterval(IFO,gpsTime).  This method returns a 
    text string giving the veto list names and intervals that
    the gpsTime intersects with.
    """
    vetoMatch=list()
    for vetoNameLong,vetoList in vList:
      for startT,stopT,KWSig in vetoList:
        if ((startT-self.tolWin)<=gpsTime<=(stopT+self.tolWin)):
          vetoMatch.append([os.path.basename(vetoNameLong).split(self.vetoExtension,1)[0],startT,stopT])
    tmpList=list()
    for a,b,c in vetoMatch:
      tmpList.append("%s %s %s\n"%(str(a),str(b),str(c)))
    outputString=str().join(tmpList)
    if vetoMatch.__len__()==0:
      outputString="%s %s %s\n"%("NONE","0","0")
    return outputString
  #End self.__checkIntervals__()

  def setIfoPath(self,ifoName="",ifoPath=""):
    """
    If you want to use something other than the hardwired path
    defaults to the veto files invoke this method with the IFO
    for which you need to adjust the path to those veto files.
    ie To adjust L1 paths do
    self.setIfopath("L1","/path/to/new/veto/files")
    """
    myIndex=-1
    for i in range(self.database.__len__()):
      if (self.database[i][0].lower()==ifoName.lower()):
        myIndex=i
    if myIndex==-1:
      raise IOError("%s Invalid IFO Name Path Unchanged."%(ifoName))
    self.database[myIndex][1]=ifoPath
  #End self.setIfoPath()

  def findInterval(self,ifoName,gpsTime):
    """
    Given an IFO name and GPS time it searches the veto lists to
    determine which if any of the veto intervals interset this
    gps time.
    """
    if not self.filesLoaded:
      self.__loadFiles__()
    myIndex=-1
    for i in range(self.database.__len__()):
      if (self.database[i][0].lower()==ifoName.lower()):
        myIndex=i
    if myIndex==-1:
      raise IOError("%s Invalid IFO Name"%(ifoName))
    vList=self.database[myIndex][2]
    outputString=self.__checkIntervals__(vList,gpsTime)
    return outputString
  #End self.findInterval()
#
#End nVeto() Class Definition
############################################################

#############################################################################
#Class to integrated SNR Ratio Testing into CBC follow up
#Amber Stuver and Cristina V. Torres
#Thu-Apr-23-2009:200904231541 
#############################################################################
# RATIO TEST CLASS
class ratioTest:
  """
  Apply the SNR ratio test to a signal coincident at two
  interferometers.  This test uses the time of flight and SNR
  strengths measured at two seperate sites to determine if the ratio
  of observed signal strengths is physically possible.  If the ratio
  is phyiscally possible but unlikely the result of this test is a
  probability that the seen SNR ratio for a given time of flight is
  physically possible.
  Contacts: Cristina Valeria Torres and Amber Stuver
  """
  def __init__(self):
    """
    This checks for information required to conduct this test.  The
    initialization of this class requires the user to input a pointer
    to the pickle file that holds all the necessary data to perform
    this test.  If no pickle file is specified then the class object
    looks for the pickle file at a predetermined http server and
    attempts to download that file into the home directory of the user
    invoking this class.
    """
    self.pickleLoaded=bool(False)
    self.pickleURL="https://ldas-jobs.ligo.caltech.edu/~ctorres/DQstuff/ratioTest.pickle"
    self.picklePath=os.path.normpath(os.getenv("HOME")+"/")
    self.pickleName="ratioTest.pickle"
    self.localPickle=os.path.normpath(self.picklePath+self.pickleName)
    self.ifoLambda={
      'LHO':{'LHO':None,'LLO':float(0.726479)},
      'LLO':{'LHO':float(0.726479),'LLO':None},
      'VIRGO':{'LHO':float(1.1066),'LLO':float(1.0753)},
      'GEO':{'LHO':float(1.0602),'LLO':float(1.1291)},
      'TAMA':{'LHO':float(1.1089),'LLO':float(1.1221)}
      }
    self.urlPattern=str("https://ldas-jobs.ligo.caltech.edu/~ctorres/DQstuff/ratioMinMax_%s_%s_hires.jpg")
    self.ifoURL={
      'LHO':{'LHO':None,'LLO':self.urlPattern%("LHO","LLO")},
      'LLO':{'LHO':self.urlPattern%("LHO","LLO"),'LLO':None},
      'VIRGO':{'LHO':self.urlPattern%("LHO","VIRGO"),'LLO':self.urlPattern%("LLO","VIRGO")},
      'GEO':{'LHO':self.urlPattern%("LHO","GEO"),'LLO':self.urlPattern%("LLO","GEO")},
      'TAMA':{'LHO':self.urlPattern%("LHO","TAMA"),'LLO':self.urlPattern%("LLO","TAMA")}
      }

    self.pickleData=dict()
  #End __init__()

  def __loadPickle__(self,path2Pickle=None):
    """
    The code required to download the pickle from a specified location
    and install the pickle into the default location for use.  This
    method is only called if the class can not find the pickle requred
    to run.  If we pass it a alternative complete path it tries to load it from
    there.  
    """
    if path2Pickle==None:
      path2Pickle=self.localPickle
    else:
      self.localPickle=path2Pickle
    self.__openPickle__()
    if not self.pickleLoaded:
      print "In future we will try to fetch a web posted pickle!\n"
      print "URL of pickle is %s\n"%self.pickleURL
      print "Download and install this pickle by hand.\n"
  #End __loadPickle__():

  def __openPickle__(self):
    """
    This method opens the pre-existing pickle file and assigns the
    information in this pickle file into variables already declared in
    the __init__() method of this class object.  
    """
    #Load pickle of this structure
    try:
      inputFP=gzip.open(self.localPickle,'rb')
      self.pickleData=cPickle.load(inputFP)
      inputFP.close()
      self.pickleLoaded=bool(True)
    except:
      self.pickleLoaded=bool(False)
  #End __openPickle__(self):

  def __createPickle__(self):
    """
    This method is used to create the needed pickle file locally from
    input text files.  This is a hard wired method that should likely
    never need to be invoked by anyone other than the authors of this
    class.
    PICKLE FORMAT for inBound determination
    X={'ifo1':{'ifo2':[list(t),list(min Ratio),list(max Ratio)]},
    ...
    }
    """
    detectorPairs={
      'LHO_GEO':'LHO_GEO_ratio.txt',
      'LHO_LLO':'LHO_LLO_ratio.txt',
      'LHO_TAMA':'LHO_TAMA_ratio.txt',
      'LHO_VIRGO':'LHO_VIRGO_ratio.txt',
      'LLO_GEO':'LLO_GEO_ratio.txt',
      'LLO_TAMA':'LLO_TAMA_ratio.txt',
      'LLO_VIRGO':'LLO_VIRGO_ratio.txt'
      }
    for key in detectorPairs.keys():
      #Load a single file.
      tVector=list()
      minRVector=list()
      maxRVector=list()
      (ifo1Name,ifo2Name)=str(key).strip().split("_")
      fp=open(detectorPairs[key],"r")
      rawData=fp.readlines()
      #print "IFO1: %s \t  IFO2: %s"%(ifo1Name,ifo2Name)
      for row in rawData:
        (a,b,c)=row.split()
        tVector.append(float(a))
        minRVector.append(float(b))
        maxRVector.append(float(c))
      if not self.pickleData.has_key(ifo1Name):
        self.pickleData[ifo1Name]=dict()
      self.pickleData[ifo1Name][ifo2Name]=[tVector,minRVector,maxRVector]
      if not self.pickleData.has_key(ifo2Name):
        self.pickleData[ifo2Name]=dict()
      self.pickleData[ifo2Name][ifo1Name]=[tVector,minRVector,maxRVector]
    #Save pickle of this structure
    outputFP=gzip.open(self.pickleName,'wb')
    cPickle.dump(self.pickleData,outputFP,2)
    outputFP.close()
  #End __createPickle__()

  def mapToObservatory(self,ifo=None):
    """
    Expects H1, V1 etc and maps it to LHO VIRGO etc.
    """
    obsMap={}
    obsMap["L1"]="LLO"
    obsMap["H1"]="LHO"
    obsMap["H2"]="LHO"
    obsMap["V1"]="VIRGO"
    obsMap["G1"]="GEO"
    obsMap["T1"]="TAMA"
    if ifo != None:
      try:
        return obsMap[ifo.upper()]
      except:
        return None
    return None
  #End mapToObservatory()

  def findURL(self,ifo1=None,ifo2=None):
    """
    Manipulates the keys as needed determine URL of relative figure from
    self.ifoURL for given detector pair. This returns a URL that
    can be used to point to a precomposed plot.
    """
    ifo1=ifo1.strip().upper()
    ifo2=ifo2.strip().upper()
    firstKeyElements=self.ifoURL.keys()
    if firstKeyElements.__contains__(ifo1):
      secondKeyElements=self.ifoLambda[ifo1].keys()
      if secondKeyElements.__contains__(ifo2):
        firstKey=ifo1
        secondKey=ifo2
      else:
        firstKey=ifo2
        secondKey=ifo1
    try:
      output=str(self.ifoURL[firstKey][secondKey])
      return output
    except:
      return None
    #End findURL()
    
  def fetchLambda(self,ifo1=None,ifo2=None):
    """
    Manipulates the keys as needed to load values from self.ifoLambda
    properly.
    """
    ifo1=ifo1.strip().upper()
    ifo2=ifo2.strip().upper()
    firstKeyElements=self.ifoLambda.keys()
    if firstKeyElements.__contains__(ifo1):
      secondKeyElements=self.ifoLambda[ifo1].keys()
      if secondKeyElements.__contains__(ifo2):
        firstKey=ifo1
        secondKey=ifo2
      else:
        firstKey=ifo2
        secondKey=ifo1
    try:
      output=self.ifoLambda[firstKey][secondKey]
      return output
    except:
      return None
    #End fetchLambda()
  
  def setPickleLocation(self,newSpot=None):
    """
    Resets the location of the pickle file to some place that is non
    standard.
    """
    if newSpot==None:
      return
    else:
      #Break path up and store parts and entire piece
      self.picklePath=os.path.dirname(newSpot)
      self.pickleName=os.path.basename(newSpot)
      self.localPickle=newSpot
    return
  #End setPickleLocation()

  def testRatio(self,ifo1="NULL",ifo2="NULL",timeOfFlight=float(1.0),\
                  myRatio=None):
    """
    A call to this method performs the check using information about
    the IFOs and the time of flight.  A float is returned by the
    method.  It contains the probablity that this SNR ratio give a
    particular time of flight is physically probably. If an error is
    encountered method returns possible values of -98,-97,-96,-95,-94.
    """
    errCode={
      -99:"ifo1 == ifo2",
      -98:"Pickle file not loaded",
      -97:"IFO lambda not found",
      -96:"Time Delay unphysical",
      -95:"Interpolation function call failure",
      -94:"TOF not found!",
      -93:"Unknown problem"
      }
    if (ifo1=="NULL") or (ifo2=="NULL") or (myRatio == None) or (myRatio==0):
      return(errCode[-99])
    if not self.pickleLoaded:
      self.__loadPickle__()
    if not self.pickleLoaded:
      return(errCode[-98])
    ifo1=ifo1.strip().upper()
    ifo2=ifo2.strip().upper()
    LV=None
    LV=self.fetchLambda(ifo1,ifo2)
    if LV == None and ifo1 != ifo2:
      return errCode[-97]
    #Check bound set ratio TO 1 return
    #Extract 3 data vectors
    (t,minR,maxR)=self.pickleData[ifo1][ifo2]
    if not(min(t)<=timeOfFlight<=max(t)):
      return float(0)
    rPrimeFunc=interpolate.interp1d(t,[minR,maxR],kind='linear')
    (newMinR,newMaxR)=rPrimeFunc(timeOfFlight)
    if (newMinR.__len__() > 1 or newMaxR.__len__() > 1):
      return errCode[-95]
    else:
      newMinR=newMinR[0]
      newMaxR=newMaxR[0]
    if (newMinR<=myRatio<=newMaxR):
      return float(1.0)
    else:
      myOutput=float(0.0)
      #P(ln(SNR)) = 1-exp(-|ln(SNR)|*0.726479)
      #This above expression from Amber's webpage a tad misleading
      #We mean 1-expcdf(abs(SNR),MU) which is 
      #exp(-SNR/mu)
      myOutput=numpy.exp(-numpy.fabs(numpy.log(myRatio)*LV))
      return myOutput
    return errCode[-94]
  #End testRatio()

  def checkPairs(self,ifoPairs=None):
    """
    Give a list of ifo pairs like [["LLO",SNR,End_Time],["LHO",SNR,End_Time]...]]
    Return the results of the Ratio test as 2D list in form
    [[IFO1,IFO2,URL,%Likely],[...],....]
    """
    pairingList=list()
    for A,a in enumerate(ifoPairs):
      for B,b in enumerate(ifoPairs):
        if (A!=B) and not pairingList.__contains__([B,b,A,a]):
          pairingList.append([A,a,B,b])
    outputList=list()
    for index1,data1,index2,data2 in pairingList:
      ifoA=self.mapToObservatory(data1[0])
      ifoB=self.mapToObservatory(data2[0])
      if ifoA != ifoB:
        gpsA=numpy.float64(data1[2])
        gpsB=numpy.float64(data2[2])
        snrA=float(data1[1])
        snrB=float(data2[1])
        try:
          snrRatio=snrA/snrB
        except:
          snrRatio=0
        gpsDiff=gpsA-gpsB
        result=self.testRatio(ifoA,ifoB,gpsDiff,snrRatio)
        myURL=self.findURL(ifoA,ifoB)
        outputList.append([ifoA,ifoB,gpsDiff,snrRatio,myURL,result])
    return outputList

  def generateHTMLTable(self,outputList=None):
    """
    Creates a small snippet of HTML for display on web browser.
    """
    resultString=(" <table border=1px>\
     <tr><th>IFO:IFO</th><th>ToF</th><th>Deff Ratio</th><th>Prob</th><th>Figure</th></tr>")
    for ifoA,ifoB,gpsDiff,snrRatio,pairURL,result in outputList:
          myURL=str('<a href="%s"><img height=150px src="%s"></a>'%(pairURL,pairURL))
          myString="<tr><td>%s:%s</td><td>%2.4f</td><td>%5.2f</td><td>%1.3f</td><td>%s</td></tr>"%\
              (ifoA,ifoB,gpsDiff,snrRatio,result,myURL)
          resultString="%s %s"%(resultString,myString)
    resultString=" %s </table>"%(resultString)
    return resultString
  
  def generateMOINMOINTable(self,outputList=None):
    """
    Creates a small snippet of MOINMOIN wiki syntax for display on web browser.
    """
    #[[ImageLink(image,target[,width=width[,height=height]][,alt=alttag])]]    
    resultString="|| IFO:IFO || ToF || Deff Ratio || Probability || Figure ||\n"
    for ifoA,ifoB,gpsDiff,snrRatio,pairURL,result in outputList:
      myURL=str('[[%s|{{%s|RatioTestPlot|width=300}}]]')%(pairURL,pairURL)
      myString="|| %s:%s || %2.4f || %5.2f || %1.3f || %s ||\n"%\
                (ifoA,ifoB,gpsDiff,snrRatio,result,myURL)
      resultString="%s%s"%(resultString,myString)
    resultString="%s \n"%(resultString)
    return resultString
# End Class ratioTest()
#############################################################################


#############################################################################
# Function to generate a coherentbank xml file
#############################################################################
def generateCohbankXMLfile(ckey,triggerTime,ifoTag,ifolist_in_coinc,search,outputPath=None,type='plot',use_col_from_installation=None):

  if outputPath:
    try:
      os.mkdir(outputPath)
    except: pass

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())

  # Cohbank files use the loudest single-ifo template for all ifos in a coinc
  maxSNR = 0;
  maxIFO = ""
  for t in ckey:
    snr = t.snr
    if snr > maxSNR:
      maxSNR = snr
      maxIFO = t.ifo
  trig = getattr(ckey,maxIFO)

  trigcopy = copy.deepcopy(trig)

  #This is a hack since data channel can differ among ifos
  #properChannelTrig = getattr(ckey,maxIFO)
  #trig.channel = properChannelTrig.channel
  # BEFORE WE MAKE A NEW TABLE FIGURE OUT WHAT COLUMNS ARE VALID !!!
  valid_columns = trigcopy.__slots__
  columns = []
  notcolumns = []
  for col in valid_columns:
    try:
      getattr(trigcopy,col)
      columns.append(col)
    except:
      notcolumns.append(col)
  # IF "use_col_from_installation" IS TRUE, ADD NEW COLUMNS TO THE SNGL_INSPIRAL TABLE
  if use_col_from_installation:
    for col in notcolumns:
      print "\n adding column " + col
      columns.append(col)
      setattr(trigcopy,col,0)

  process_table = lsctables.New(lsctables.ProcessTable)
  xmldoc.childNodes[-1].appendChild(process_table)

  process_params_table = lsctables.New(lsctables.ProcessParamsTable)
  xmldoc.childNodes[-1].appendChild(process_params_table)

  search_summary_table = lsctables.New(lsctables.SearchSummaryTable)
  xmldoc.childNodes[-1].appendChild(search_summary_table)
  for chunk in search:
    out_start_time = float(chunk.out_start_time)
    out_start_time_ns = float(chunk.out_start_time_ns)/1000000000
    out_end_time = float(chunk.out_end_time)
    out_end_time_ns = float(chunk.out_end_time_ns)/1000000000
    if ( (triggerTime >= (out_start_time+out_start_time_ns)) and (triggerTime <= (out_end_time+out_end_time_ns)) ):
      search_summary_table.append(chunk)
      break

  sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable,columns)
  xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
  # Each coherent bank file should have trig rows for all ifos in a coinc
  #for ifo in ifolist_in_coinc:
  # Each coherent bank file should have trig rows for all ifos in science mode

  for j in range(0,len(ifoTag)-1,2):
    itf = ifoTag[j:j+2]
    trigger = copy.deepcopy(trigcopy)
    trigger.ifo = itf
    sngl_inspiral_table.append(trigger)

  fileName = ifoTag + '-COHBANK_FOLLOWUP_' + str(int(ckey.event_id)) + '-' + str(int(triggerTime)) + "-2048.xml.gz"
  if outputPath:
    fileName = outputPath + '/' + fileName
  utils.write_filename(xmldoc, fileName, verbose = False, gz = True)

  #Also write input file for clustering code "cohire"
  chiaFileName = 'followUpChiaJob/' + ifoTag + '-CHIA_1_' + str(int(ckey.event_id)) + '-' + str(int(triggerTime)-1) + "-2.xml.gz"
  cohireInputFile = ifoTag + '-COHIRE_FOLLOWUP_' + str(int(ckey.event_id)) + '-' + str(int(triggerTime)-1) + "-2.txt"
  if outputPath:
    cohireInputFile = outputPath + '/' + cohireInputFile
  ff = open(cohireInputFile,'w')
  print >>ff, chiaFileName

  return maxIFO


class followupDQV:
  """
  This class is intended to provide a mechanism to access DQ segment
  information and veto segment information put into the segment
  database.  This class will replace the previously defined class of
  followupdqdb.  This class also provides access to the needed code
  to create a local background of DQ statistics for improved
  background studies.
  Contact: Cristina Valeria Torres
  """
  def __init__(self,LDBDServerURL=None,quiet=bool(False),pickle=None,blinded=False):
    """
    This class setups of for connecting to a LDBD server specified at
    command line to do segment queries as part of the follow up
    pipeline.  If the user does not specify the LDBD server to use the
    method will use the environment variable S6_SEGMENT_SERVER to
    determine who to query.  The LDBD URL should be in the following form
    ldbd://myserver.domain.name:808080. You can specify the path to a
    background DQ pickle if the path is valid the class opens it
    otherwise it queries the segment DB and builds the pickle.
    Warning! recreating the DQ background is VERY slow.
    """
    self.__connection__= None
    self.__engine__= None
    self.__installPath__=home_dir()+"/ctorres/followupbackgrounds/dq/"
    self.__blinded__=blinded
    self.__blindFlags__=[\
      "DMT-INJECTION_INSPIRAL",
      "DMT-INJECTION"\
      ]
    if pickle==None:
      self.__backgroundPickle__=None
    else:
      self.__backgroundPickle=pickle
      self.__backgroundPickle__=os.path.expanduser(self.__backgroundPickle__)
    self.__haveBackgroundDict__=bool(False)
    self.__havecategories__=bool(False)
    #A dict for a dict of GPStimes and list all flags seen
    self.__backgroundDict__=dict()
    self.__category__=dict()
    #Access a dict of dicts for a flag names with % stored
    self.__backgroundResults__=dict()
    self.__backgroundPoints__=int(1000)
    self.__columns__=["Ifo","Flag","Ver","Start","Offset",\
                      "Stop","Offset","Size","DQ Rank","Cat(s)"]
    #Dict should be a dict of lists
    self.__backgroundTimesDict__=dict()
    self.ifos=interferometers
    self.ifos.sort()
    self.triggerTime=int(-1)
    self.serverURL=defaultsegmentserver
    if LDBDServerURL==None:
      envServer=None
      envServer=os.getenv('S6_SEGMENT_SERVER')
      if envServer!=None:
        self.serverURL=envServer
      sys.stderr.write("Warning no LDBD Server URL specified \
defaulting to %s\n"%(self.serverURL))
    else:
      self.serverURL=LDBDServerURL
    self.resultList=list()
    #New more robust query
    self.dqvQuery="""SELECT \
    segment_definer.ifos, \
    segment_definer.name, \
    segment_definer.version, \
    segment_definer.comment, \
    segment.start_time, \
    segment.end_time \
    FROM segment,segment_definer \
    WHERE \
    segment_definer.segment_def_id = segment.segment_def_id AND \
    segment_definer.creator_db = segment.segment_def_cdb AND \
    NOT (segment.start_time > %s OR %s > segment.end_time) \
    ORDER BY segment_definer.ifos,segment_definer.name,segment.start_time \
    """
    tmpBlind=""
    for blindFlag in self.__blindFlags__:
      tmpBlind+=""" AND NOT (segment_definer.name = '%s') """%blindFlag
    self.dqvQueryBlinded=self.dqvQuery.split("ORDER")[0]+\
                          tmpBlind+\
                          " ORDER "+\
                          self.dqvQuery.split("ORDER")[1]
  #End __init__()


  def __merge__(self,inputList=None):
    """
    Takes an input list of tuples representing start,stop and merges
    them placing them in time order when returning the coalesced list of
    tuples.
    """
    outputList=list()
    if type(inputList) != type(list()):
      sys.stderr.write("Wrong variable type passed as argument in\
 followupDQV.__merge__()\n")
      return None
    if inputList.__len__() < 1:
      return  inputList
    inputList.sort()
    while inputList:
        segA=inputList.pop()
        overlap=True
        #Assume next segment overlaps segA
        while overlap:
            #Pop of next segment if available
            if inputList.__len__() > 0:
                segB=inputList.pop()
            else:
                #No overlap possible no segs left!
                segB=(-1,-1)
                overlap=False
            #Three cases of intersection
            #Overlap Left
            if (
                (segB[0]<= segA[0] <= segB[1])
                and
                (segA[1] >= segB[1])
                ):
                segA=(segB[0],segA[1])
            #Overlap Right
            elif (
                  (segB[0]<= segA[1] <= segB[1])
                  and
                  (segA[1] <= segB[0])
                 ):
                segA=(segA[0],segB[1])
            #Bridge over
            elif (
                (segB[0]<=segA[0])
                and
                (segB[1]>=segA[1])
                ):
                segA=(segB[0],segB[1])
            else:
                #Put segment back there was no overlap!
                if not((-1,-1)==segB):
                  inputList.append(segB)
                overlap=False
        outputList.append(segA)
        outputList.sort()
    return outputList
  #End __merge__() method

  def getInstallPath(self):
    """
    Returns a string pointing to the local install location of where
    the DQ background pickles should be found.
    """
    return self.__installPath__
    
  def figure_out_pickle(self,ifoEpochList=None):
    fileMask="followup_background_%s.pickle"
    installPath=self.__installPath__
    if ifoEpochList==None:
      return None
    elif type(ifoEpochList) == type(str()):
      if ifoEpochList=="automatic":
        ifoEpochList=[(x,getRunEpoch(self.triggerTime,x)) \
                      for x in self.ifos]
      else:
        sys.stderr.write("Given input of %s can't parse it!\n"%ifoEpochList)
        return None
    ifoEpochList.sort()
    tmpString=""
    for A,B in ifoEpochList:
      tmpString+="%s:%s_"%(A,B)
    tmpString.rstrip("_")
    return installPath+fileMask%(tmpString)
  #End figure_out_pickle

  def resetBackgroundPoints(self,pointCount=None):
    """
    Used to manipulate the number of points to create background with.
    """
    if pointCount==None:
      return
    else:
      self.__backgroundPoints__=int(pointCount)

  def getPicklePointer(self):
    """
    Returns the pickle to be used by function.
    """
    return self.__backgroundPickle__

  def resetPicklePointer(self,filename=None):
    """
    If you called the class definition with the wrong pickle path.
    You can reset it with this method.
    """
    if filename==None:
      sys.stdout.write("Path information to background pickle unchanged.\n")
    elif filename.__contains__("~"):
      self.__backgroundPickle__=os.path.expanduser(filename)
    else:
      self.__backgroundPickle__=filename
  #End resetPicklePointer
  
  def query(self,queryString=None):
    """
    Simple wrapper method to do query and return result
    """
    closeConnection=False
    if queryString == None:
      return []
    if (self.__engine__ == None) or \
           (self.__connection__ == None):
      self.__connectToSegmentDB__()
      closeConnection=True
    try:
      myOutput=self.__engine__.query(queryString)
    except Exception, errMsg:
      sys.stderr.write("Query failed %s \n"%(self.serverURL))
      sys.stdout.write("Error fetching query results.\n")
      sys.stderr.write("Error message seen: %s\n"%(str(errMsg)))
      sys.stderr.write("Query Tried: \n %s \n"%(queryString))
      return []
    if closeConnection:
      self.__disconnectFromSegmentDB__()      
    #Clean up extra space in query results
    myFinalOutput=list()
    for row in myOutput:
      newRow=list()
      for col in row:
        if type(col)==str:
          newRow.append(col.strip())
        else:
          newRow.append(col)
      myFinalOutput.append(newRow)
    return myFinalOutput
  #End query

  def fetchInformation(self,triggerTime=None,window=300):
    """
    Wrapper for fetchInformationDualWindow that mimics original
    behavior
    """
    return self.fetchInformationDualWindow(triggerTime,window,window,ifoList='DEFAULT')

  def __connectToSegmentDB__(self,serverURL=None):
    """
    Private method to execute connection to segment DB
    """
    if serverURL==None:
      serverURL=self.serverURL
    try:
      self.__connection__=segmentdb_utils.setup_database(serverURL)
    except Exception, errMsg:
      sys.stderr.write("Error connecting to %s\n"\
                         %(serverURL))
      sys.stderr.write("Error Message :\t %s\n"%(str(errMsg)))
      self.__connection__=None
      self.resultList=list()
    try:
      self.__engine__=query_engine.LdbdQueryEngine(self.__connection__)
    except Exception, errMsg:
      sys.stderr.write("Error building query engine using %s\n"\
                         %(serverURL))
      sys.stderr.write("Error Message :\t %s\n"%(str(errMsg)))
      self.__engine__=None
      self.resultList=list()
    return

  def __disconnectFromSegmentDB__(self):
    """
    Private method to close query engine.
    """
    self.__engine__.close()
    self.__engine__=None
    return
  
  def fetchInformationDualWindow(self,triggerTime=None,frontWindow=300,\
                                 backWindow=150,ifoList='DEFAULT'):
    """
    This method is responsible for queries to the data server.  The
    results of the query become an internal list that can be converted
    into an HTML table.  The arguments allow you to query with trigger
    time of interest and to change the window with each call if
    desired. The version argument will fetch segments with that
    version or higher.
    """
    if ifoList=="DEFAULT":
      ifoList=interferometers
    if (ifoList == None) or \
       (len(ifoList) < 1):
      sys.stderr.write("Ifolist passed is malformed! : %s\n"%ifoList)
      return
    #Set the internal class variable self.ifos
    self.ifos=ifoList
    if sum([x.upper() in interferometers for x in ifoList]) < 1:
      sys.stderr.write("Valid ifos not specified for DQ lookups. %s\n"%ifoList)
      return
    triggerTime=float(triggerTime)
    if triggerTime==int(-1):
      sys.stdout.write("Specify trigger time please.\n")
      return
    else:
      self.triggerTime = float(triggerTime)
    gpsEnd=int(triggerTime)+int(backWindow)
    gpsStart=int(triggerTime)-int(frontWindow)
    #If this is a blinded check use blinded query
    if self.__blinded__:
      self.resultList=self.query(self.dqvQueryBlinded%(gpsEnd,gpsStart))
    else:
      self.resultList=self.query(self.dqvQuery%(gpsEnd,gpsStart))
    if len(self.resultList) < 1:
      sys.stdout.write("Query Completed, Nothing Returned for time %s.\n"%(triggerTime))
    #Coalesce the segments for each DQ flag
    #Reparse the information
    newDQSeg=list()
    if self.resultList.__len__() > 0:
      #Obtain list of all flags, ignore IFOs not specified
      uniqSegmentName=list()
      for ifo,name,version,comment,start,end in self.resultList:
        if (not uniqSegmentName.__contains__((ifo,name,version,comment))) and \
               (ifo.strip().upper() in ifoList):
          uniqSegmentName.append((ifo,name,version,comment))
        #Add the SCIENCE segment no matter which IFOs are specified!
        if ((name.lower().__contains__('science')) and \
            not (ifo.strip().upper() in ifoList)):
          uniqSegmentName.append((ifo,name,version,comment))
      #Save textKey for all uniq segments combos
      for uifo,uname,uversion,ucomment in uniqSegmentName:
        segmentIntervals=list()
        #Extra segments based on uniq textKey
        for ifo,name,version,comment,start,end in self.resultList:
          if (uifo,uname,uversion,ucomment)==(ifo,name,version,comment):
            segmentIntervals.append((start,end))
        segmentIntervals.sort()
        #Coalesce those segments
        newSegmentIntervals=self.__merge__(segmentIntervals)
        #Write them to the object which we will return
        for newStart,newStop in newSegmentIntervals:
          newDQSeg.append([uifo,uname,uversion,ucomment,newStart,newStop])
        newDQSeg.sort()
        del segmentIntervals
    #Reset the result list to the IFO restricted set
    self.resultList=newDQSeg
    return newDQSeg
  #End method fetchInformation()

  def createDQbackground(self,ifoEpochList=list(),pickleLocale=None):
    """
    Two inputs a list of tuples (ifo,epochname) for each instrument.
    Also a place to save the potential pickle to for quick access
    later.    """
    if type(ifoEpochList) != type(list()):
      raise Exception, \
            "Invalid input argument ifoEpochList,%s type(%s)"\
            %(ifoEpochList,type(ifoEpochList))
    #Make sure epoch exists for reach ifo
    for ifo,epoch in ifoEpochList:
      if ifo not in runEpochs.keys():
        raise Exception, "Bad ifo specified, %s"%ifo
      if epoch not in runEpochs[ifo].keys():
        raise Exception, "Bad ifo epoch specified, %s:%s"%(ifo,epoch)
    #If pickle location given try to load that pickle first.
    backgroundPickle=False
    if pickleLocale!=None:
      #If pickle file exists read it if not make sure we can
      #generate it properly otherwise skip creating background
      if os.path.isfile(pickleLocale):
        try:
          self.__backgroundDict__=cPickle.load(file(pickleLocale,'r'))
          self.__haveBackgroundDict__=True
          backgroundPickle=True
        except:
          backgroundPickle=False
          self.__haveBackgroundDict__=False
          sys.stderr.write("Error importing the pickle file! %s\n"\
                           %(pickleLocale))
          return
        for (ifo,epoch) in ifoEpochList:
            if (ifo.upper().strip(),epoch.upper().strip()) \
                   not in self.__backgroundDict__["ifoepoch"]:
              raise Exception,\
                    "Invalid ifo and epoch information in \
generated background expected %s got %s"%(\
                      self.__backgroundDict__["ifoepoch"],
                      ifoEpochList)
        return
    else:
      #Assume an automatic pickle creation
      autoPath=self.figure_out_pickle("automatic")
      (dirname,filename)=os.path.split(autoPath)
      if not os.access(dirname,os.W_OK):
        sys.stderr.write("Warning: Insufficient disk \
permissions to create DQ background pickle file:%s.\n"%(autoPath))
        backgroundPickle=False
        self.__haveBackgroundDict=False
        return
    #Setup for large volumne of queries
    #Determine random background times for each IFO
    sciSegments=dict()
    for myIfo,myEpoch in ifoEpochList:
      #Find epoch start and stop
      if myIfo.strip().lower() == "v1":
        mySegName="ITF_SCIENCEMODE"
      else:
        mySegName="DMT-SCIENCE"
      (epochStart,epochStop)=getRunTimes(myEpoch,myIfo)
      tmpResults=getSciSegs(myIfo.strip(),\
                            epochStart,\
                            epochStop,\
                            cut=True,\
                            serverURL=self.serverURL,\
                            segName=mySegName,\
                            seglenmin=None,\
                            segpading=0)
      sciSegments[myIfo]=[(x.start(),x.end()) for x in tmpResults]
      if len(sciSegments[myIfo])<1:
        sys.stderr.write("Warning:No segments found for %s %s\n"%(myIfo,myEpoch))
    #Create list of random times in segments of sciSegments[ifo]
    for myIfo in sciSegments.keys():
      #Background Timepoint
      mySampleCount=0
      emergencyStop=int(1.50*self.__backgroundPoints__)
      while mySampleCount < self.__backgroundPoints__:
        if mySampleCount>emergencyStop:
          sys.stderr.write("Aborting creation of DQ background timelist!\n")
          os.abort()
        myPoint=numpy.random.uniform(epochStart,epochStop,1)
        nearestIndex=numpy.searchsorted([x for x,y in sciSegments[myIfo]],\
                                        int(myPoint),\
                                        side='left')
        if nearestIndex > 1:
          nearestIndex=nearestIndex-1
          #Index is next segment past time in question
          (start,end)=sciSegments[myIfo][nearestIndex]
        else:
          #Exception is first segment might hold point
          (start,end)=sciSegments[myIfo][nearestIndex]          
        if (start<=myPoint) and (myPoint<=end):
          if myIfo not in self.__backgroundTimesDict__.keys():
            self.__backgroundTimesDict__[myIfo]=list()
          self.__backgroundTimesDict__[myIfo].append("%10.0f"%myPoint)
          mySampleCount=mySampleCount+1
    #Take each list of ifo times and look up the DQ flags
    self.__connectToSegmentDB__()
    for myIfo,myTimes in self.__backgroundTimesDict__.iteritems():
      for myTime in myTimes:
        myQuery=self.dqvQuery%(myTime,myTime,myTime,myTime)
        #Insert befor "ORDER BY"
        myQueryA,myQueryB=myQuery.split("ORDER BY",1)
        myQuery=myQueryA+\
                 """ AND segment_definer.ifos = '%s' """%myIfo+\
                 """ ORDER BY """+myQueryB
        if myIfo not in self.__backgroundDict__.keys():
          self.__backgroundDict__[myIfo]=dict()
        self.__backgroundDict__[myIfo]["%s"%myTime]=self.query(myQuery)
    self.__disconnectFromSegmentDB__()
    #Integrate epoch listing into dict()
    self.__backgroundDict__["ifoepoch"]=[(myIfo.strip().upper(),
                                          myEpoch.strip().upper()) \
                                         for myIfo,myEpoch in ifoEpochList]
    self.__haveBackgroundDict__=True
    #Save the created DQ background to a pickle, skip saving on error!
    #That is assuming we didn't get our data from a pickle already!
    if not backgroundPickle:
      try:
        cPickle.dump(self.__backgroundDict__,file(pickleLocale,'w'))
      except:
        sys.stdout.write("Problem saving pickle of DQ information.")
        sys.stdout.write("Trying to place pickle in your home directory.")
        try:
          cPickle.dump(self.__backgroundDict__,
                       file(home_dir()+"/"+os.path.basename(pickleLocale),'w'))
        except:
          sys.stdout.write("Really ignoring pickle generation now!\n")
  #End createDQbackground

  def estimateDQbackground(self):
    """
    This method looks at the self.resultlist inside the instance.
    Using this and 1000 generated time stamp it tabulates a
    ranking of flag prevelance, binomial probability 'p'
    """
    if len(self.resultList) < 1:
      self.__backgroundResults__=list()
      self.__backgroundTimesDict__=dict()
      self.__backgroundDict__=dict()
      self.__haveBackgroundDict__=bool(False)
      return
    #Create DQ background, specify pickle locale to try and load first
    #Determine which IFOs from DQ listing
    uniqIfos=list()
    for ifo,b,c,d,e,f in self.resultList:
      if ifo not in uniqIfos:
        uniqIfos.append(ifo)
    ifoEpochList=[(x,getRunEpoch(self.triggerTime,x)) for x in self.ifos]
    self.createDQbackground(ifoEpochList,self.__backgroundPickle__)
    if not self.__haveBackgroundDict__:
      sys.stderr.write("Could not either open or save DQ \
background in %s no background data available!\n"%(self.__backgroundPickle__))
      self.__backgroundResults__=list()
      self.__backgroundTimesDict__=dict()
      self.__backgroundDict__=dict()
      self.__haveBackgroundDict__=bool(False)
      return
    #Calculate the binomial 'p' value for the flags in the table.
    if self.resultList < 1:
          sys.stderr.write("Aborting tabulation of binomial P\n")
          os.abort()
    seenFlags=dict()
    for ifo,name,version,comment,start,stop in self.resultList:
      if ifo.strip() not in seenFlags.keys():
        seenFlags[ifo]=list()
      seenFlags[ifo].append(name)
    for myIfo in seenFlags.keys():
      tmpFlags=list()
      #Handles flags seen but not present in pickle structure
      if myIfo.strip() not in self.__backgroundDict__.keys():
        if myIfo not in self.__backgroundResults__.keys():
          self.__backgroundResults__[myIfo]=dict()
        for outsideFlag in seenFlags[myIfo]:
          self.__backgroundResults__[myIfo][outsideFlag]=float(-0.0)
      else:
        #Handles all flags seen compared to background pickle
        for backgroundTime,backgroundFlags in \
                self.__backgroundDict__[myIfo.strip()].iteritems():
          tmpFlags.extend([name for ifo,name,ver,com,start,stop in backgroundFlags])
        if myIfo not in self.__backgroundResults__.keys():
          self.__backgroundResults__[myIfo]=dict()
        for myFlag in seenFlags[myIfo]:
          self.__backgroundResults__[myIfo][myFlag]=tmpFlags.count(myFlag)/float(self.__backgroundPoints__)
    self.__haveBackgroundDict__=True
    #Return background estimating
    
                                            
  def getDQbackgroundList(self):
    """
    Simple wrapper to return a list of DQ background information.
    """
    return self.__backgroundResults__
  
  def getResultList(self):
    """
    Simple calling function to create a list object of the results.
    """
    return self.resultList
  
  def generateResultList(self):
    """
    Simple calling function to create a list object of the results.
    """
    return self.getResultList()
  #End generateResultList
  
  def generateHTMLTable(self,tableType="BOTH"):
    """
    Return a HTML table already formatted using the module MARKUP to
    keep the HTML tags complient.  This method does nothing but return
    the result of the last call to self.fetchInformation() The flag
    names associated with LIGO will have links to the channel wiki in
    them also.
    Types that will invoke a not everything behaviour are
    DQ and VETO
    """
    ligo=["L1","H1","H2","V1"]
    channelWiki="https://ldas-jobs.ligo.caltech.edu/cgi-bin/chanwiki?%s"
    if self.triggerTime==int(-1):
      return ""
    myColor="grey"
    tableString="<table bgcolor=grey border=1px>"
    titleString="<tr>"
    tableEmptyString="<tr bgcolor=%s>"%myColor
    rowString="<tr bgcolor=%s> "
    for col in self.__columns__:
      titleString+="<th>%s</th>"%col
      rowString+="<td>%s</td>"
      tableEmptyString+="<td>None</td>"
    titleString+="</tr>\n"
    tableEmptyString+="</tr>\n"
    rowString+="</tr>\n"
    tableString+=titleString
    if len(self.resultList) == 0:
      tableString+=tableEmptyString
    for ifo,name,version,comment,start,stop in self.resultList:
      #If we have background information fetch it
      if self.__haveBackgroundDict__:
        myBackgroundRank=str("%3.1f"%(100.0*self.__backgroundResults__[ifo][name])).rjust(5)
      else:
        myBackgroundRank="None"        
      if self.__havecategories__:
        myCategory=self.__category__[ifo][name]
      else:
        myCategory="None"
      offset1=start-self.triggerTime
      offset2=stop-self.triggerTime
      size=int(stop-start)
      if (offset1>=0) and (offset2>=0):
        myColor="green"
      if (offset1<=0) and (offset2<=0):
        myColor="yellow"
      if (offset1<=0) and (offset2>=0):
        myColor="red"
      if name.lower().__contains__('science'):
        myColor="skyblue"
      if tableType.upper().strip() == "DQ":
        if not name.upper().startswith("UPV"):
          tableString+=rowString%(myColor,ifo,name,version,\
                                  start,offset1,stop,offset2,\
                                  size,myBackgroundRank,myCategory)
      elif tableType.upper().strip() == "VETO":
        if name.upper().startswith("UPV"):
          tableString+=rowString%(myColor,ifo,name,version,\
                                  start,offset1,stop,offset2,\
                                  size,myBackgroundRank,myCategory)
      elif tableType.upper().strip() not in ["VETO","DQ"]:
        tableString+=rowString%(myColor,ifo,name,version,\
                                start,offset1,stop,offset2,size,\
                                myBackgroundRank,myCategory)
    tableString+="</table>"
    return tableString
  #End method generateHTMLTable()

  def generateMOINMOINTable(self,tableType="BOTH"):
    """
    Return a MOINMOIN table.
    """
    ligo=["L1","H1","H2","V1"]
    channelWiki="https://ldas-jobs.ligo.caltech.edu/cgi-bin/chanwiki?%s"
    if self.triggerTime==int(-1):
      return ""
    myColor="grey"
    tableString=""
    titleString=""
    emptyRowString=""
    rowString=""
    for i,col in enumerate(self.__columns__):
      if i == 0:
        titleString+="""||<rowbgcolor="%s"> %s """%(myColor,col)
        rowString+="""||<rowbgcolor="%s"> %s """
        emptyRowString+="""||<rowbgcolor="%s"> None """%myColor
      else:
        titleString+="""|| %s """%col
        rowString+="""|| %s """
        emptyRowString+="""|| None """
    titleString+="""||\n"""
    rowString+="""||\n"""
    emptyRowString+="""||\n"""
    tableString+=titleString
    #Extract only DQ row or only VETO rows
    tmpResultList=list()
    for myRow in self.resultList:
      ifo,name,version,comment,start,stop=myRow
      #Select base on table type
      if ((tableType.upper() == "DQ") and \
          (not name.strip().upper().startswith("UPV"))):
        tmpResultList.append(myRow)
      elif ((tableType.upper() == "VETO") and \
            (name.strip().upper().startswith("UPV"))):
        tmpResultList.append(myRow)
      elif tableType.upper().strip() not in ["VETO","DQ"]:
        tmpResultList.append(myRow)
    if len(tmpResultList) == 0:
      tableString+=emptyRowString
    for ifo,name,version,comment,start,stop in tmpResultList:
      #If we have background information fetch it
      if self.__haveBackgroundDict__:
        myBackgroundRank=str("%3.1f"%(100.0*self.__backgroundResults__[ifo][name])).rjust(5)
      else:
        myBackgroundRank="None"        
      if self.__havecategories__:
        myCategory=self.__category__[ifo][name]
      else:
        myCategory="None"
      offset1=start-self.triggerTime
      offset2=stop-self.triggerTime
      size=int(stop-start)
      if (offset1>=0) and (offset2>=0):
        myColor="green"
      if (offset1<=0) and (offset2<=0):
        myColor="yellow"
      if (offset1<=0) and (offset2>=0):
        myColor="red"
      if name.lower().__contains__('science'):
        myColor="skyblue"
      tableString+=rowString%(myColor,str(ifo).strip(),name,version,\
                              start,offset1,stop,offset2,size,\
                              myBackgroundRank,myCategory)
    tableString+="\n"
    return tableString

  
#End class followupDQV

######################################################################
#New Class Definition for determining DQ segments active for a given
#GPS time 
#
class followupdqdb:
    """
    This class provides a method for doing DQ flag queries required
    by the follow up pipeline when creating the automatic checklists.
    Invoking this class requires a boolean option if you 
    need to override the default behavior for this class
    initialization.
    To override the default initialization call like:
    dq=followupDQdb(True)
    else just do not worry about it and everything should work.
    """
    def __init__(self,default=False):
        """
        The __init__ method which can be overrridden using 
        other methods defined in this class.
        """
        self.defaultVersion=99
        self.activeRecords="1"
        self.db=None
        self.dbSocket=None
        self.myHome=None
        self.myHome=os.getenv("HOME")
        self.urlPattern="http://ldas-cit.ligo.caltech.edu/segments/S5/%s/dq_segments.txt"
        self.pathPattern="%s/%s/dq_segments.txt"
        self.sqlFile="/followupDQ.sqlite"
        self.sqlPath=self.myHome+"/"
        self.ifoList=["L1","H1","H2","V1"]
        self.dbFiles=[]
        for ifo in self.ifoList:
            thisFile=self.pathPattern%(self.myHome,ifo)
            thisUrl=self.urlPattern%(ifo)
            thisIfo=ifo
            thisFileExists=os.path.isfile(thisFile)
            self.dbFiles.append([thisIfo,thisFile,thisUrl,thisFileExists])
        #Default Behavior Code Below
        if not(default) and os.path.isfile(self.sqlPath+self.sqlFile):
            self.__connectDB__()
            try:
                x=self.__createRawCursor__()
                for table in self.ifoList:
                    x.execute("select * from %s limit 100"%(table))
                del x
            except:
                sys.stderr.write("Sqlite database at %s seems\
 corrupted. It will be rebuilt.\n"%(self.sqlPath+self.sqlFile))
                self.db.close()
                self.db=None
                os.unlink(self.sqlPath+self.sqlFile)
                self.createDB()
                self.__connectDB__()
        elif not(default) and not(os.path.isfile(self.sqlPath+self.sqlFile)):
            self.createDB()
            self.__connectDB__()
        else:
            self.db=None
    #End __init__()

    def __connectDB__(self):
        """
        Simple method to wrap connections to the DB to start queries.
        """
        self.db=sqlite.connect(self.sqlPath+self.sqlFile)
    #End connectDB()

    def __disconnectDB__(self):
        """
        Wrapper method to close the database file.
        """
        self.db.commit()
        self.db.close()
    #End __disconnectDB__()

    def close(self):
      """
      This explicity allows the user to close the database when done
      making all the queries we want.  Closing the DB helps prevent
      corruption and clears existing locks to the DQ so others can add
      elements to the DB.
      """
      self.__disconnectDB__()
    #End close()

    def __createRawCursor__(self):
        """
        A method that is only used for testing.  It return a cursor 
        object to the database currently open.
        """
        return self.db.cursor()
    #End __createRawCursor__()

    def createDB(self):
        """
        Method to create the sqlite database and install it to 
        location sqlPath.
        """
        sys.stderr.write("Trying to create SQLite database.\n")
        self.__connectDB__()
        self.dbSocket=self.db.cursor()
        for table,file,url,exists in self.dbFiles:
            sys.stderr.write("Adding table for %s\n"%(table))
            if not exists:
                sys.stderr.write("Missing %s fetching it from %s \n"%(file,url))
                sys.stderr.write("Downloading... please wait.\n")
                myPath=os.path.dirname(file)
                if not os.path.exists(myPath):
                    sys.stderr.write("File path not found creating\
path %s\n"%(myPath))
                    os.makedirs(myPath)
                urllib.urlretrieve(url,filename=file)
                sys.stderr.write("Downloaded. Saved to %s \n"%(file))
                sys.stderr.write("Converting to SQL syntax and inserting.\n")
            commandString="create table %s (flag text,version\
 integer,start integer,stop integer, active integer)"%(table)
            self.dbSocket.execute(commandString)
            fp=open(file,"r")
            for row in fp.readlines():
                token=None
                if len(row) > 0 and row[0] != "#":
                    token=row.split()
                    if token.__len__() < 5:
                        print "Problems with file %s"%(file)
                        self.dbSocket.commit()
                        db.close()
                        os.unlink(self.sqlPath+self.sqlFile)
                        os.abort()
                    else:
                        commandString2="insert into %s \
(flag,version,start,stop,active) values ('%s',%i,%i,%i,%i)"\
%(table,token[0],int(token[1]),int(token[2]),int(token[3]),int(token[4]))
                        self.dbSocket.execute(commandString2)
            self.db.commit()
            fp.close()
        self.db.close()
        sys.stderr.write("Created SQLite database and closed it.\n")
    #End method createDB            


    def setURL(self,ifo=None,url=None):
        """
        This method allows use to specify a ifo L1 etc and the URL
        to use to go onto the web fetching the appropriate DQ segments
        file. If the method is called with invalid options it does
        nothing.
        """
        if not(ifo.upper() in self.ifoList) or url == None:
            return
        for i in range(0,self.dbFiles.__len__()):
            if self.dbFiles[i][0] == ifo:
                self.dbFiles[i][2]=url
        return
    #End setURL()

    def setFile(self,ifo=None,file=None):
        """
        This method allows you to specify an override location
        specifying where or where to install the ASCII dq_segments
        type file for a given IFO.
        """
        if not(ifo.upper() in self.ifoList) or file == None:
            return
        for i in range(0,self.dbFiles.__len__()):
            if self.dbFiles[i][0] == ifo:
                self.dbFiles[i][1]=file
                if os.path.isfile(self.dbFiles[i][1]):
                    self.dbFiles[i][3]=True
                else:
                    self.dbFiles[i][3]=False
        return
    #End setFile()

    def setDB(self,file=None):
        """
        This method set the location to install the sqlite file to
        or the location of a preinstalled sqlite filed derived from
        ASCII dq_segments files.
        """
        if file==None:
            return
        filename=os.path.normpath(file)
        self.sqlFile=os.path.basename(filename)
        self.sqlPath=os.path.dirname(filename)
        if os.path.isdir(self.sqlPath):
          self.sqlPath=self.sqlPath+"/"
        return

    def setVersion(self,version=None):
        """
        This method sets the version number or greater to fetch from
        the segment lists sqlite database.
        """
        if (version == None):
           return
        self.defaultVersion=version
        return

    def queryDB(self,gps=None,window=int(30)):
        """
        Takes two arguments the gps time as an integer and the window
        to search around this time with default to 30 seconds. It
        returns a dictionary variable assesible via
        X["ifo"] which will be a list of [flag,version,start,stop,active]
        
        """
        try:
            dbSocket=self.__createRawCursor__()
        except:
            self.__connectDB__()
            dbSocket=self.__createRawCursor__()
        results=dict()
        for table in self.ifoList:
            results[table]=list()
            iStart=int(gps)-int(window)
            iStop=int(gps)+int(window)
            commandString=\
                "select * from %s where\
     (((start<=%i) and (stop >= %i)) or\
    ((start>=%i) and (stop <= %i)) or\
    ((start<=%i) and (stop >= %i))) and (active == 1)\
     and version == %i \
     order by flag,version desc"%\
    (table,\
    iStart,iStart,\
    iStart,iStop,\
    iStop,iStop,self.defaultVersion)
            dbSocket.execute(commandString)            
            results[table]=dbSocket.fetchall()
        return results
    #End def queryDB()

#End followupDQdb class()
######################################################################

#A loose method to retrieve the iLog url given a integer for of
#GPStime
def getiLogURL(time=None,ifo=None):
  """
  This method returns a URL string to point you to ilog day page for
  specified IFO and GPStime. Valid IFO labels are V1, L1, H1 or H2.
  """
  #Example URL to give to query between two dates, still need to click
  #search button but this is a start(Dec 4, 2009 -> Dec 24,2009)
  #https://pub3.ego-gw.it/logbook/index.php?area=logbook&ref=detailedsearch&datefrom=04/12/2009&dateto=24/12/2009
  #Due to difficulty in navigating virgo ilog search (t0-1 day,t0+1 day)
  dateString="%s/%s/%s"
  urls={
    'default':"http://www.ligo.caltech.edu/~pshawhan/scilinks.html",
    'V1':"https://pub3.ego-gw.it/logbook/",
    'L1':"http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s\
&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to=",
    'H1':"http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s\
&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to=",
    'H2':"http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s\
&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to="
    }
  outputURL=urls['default']
  if ((ifo==None) or (time==None)):
    return urls['default']
  gpsTime=LIGOTimeGPS(time)
  Y,M,D,doy,h,m,s,ns,junk=xlaldate.XLALGPSToUTC(gpsTime)
  gpsStamp=dateString%(str(M).zfill(2),str(D).zfill(2),str(Y).zfill(4))
  if ('H1','H2','L1').__contains__(ifo.upper()):
    outputURL=urls[ifo.upper()]%gpsStamp
  if ('V1').__contains__(ifo.upper()):
    outputURL=urls[ifo.upper()]
  return outputURL
#End def getiLogURL

def getGlitchReportURL(time=None):
  """
  This method is esentially a wrapper method until we have a better
  approach to linking directly to a specific glitch report. The method
  expects an interger respresentation of GPS time.
  """
  stopS5=int(875232014)
  defaultURL="https://www.lsc-group.phys.uwm.edu/twiki/bin/view/DetChar/GlitchStudies"
  s5URL="http://lancelot.mit.edu/~dicredic/S5scimon.html"
  if time==None:
    return defaultURL
  if int(time) <= stopS5:
    return s5URL
  else:
    return defaultURL
#End getGlitchReportURL

def getDailyStatsURL(time=None):
  """
  This method points you to the right URL to look at the daily stats
  pages.
  """
  stopS5=int(875232014)
  defaultURL="http://blue.ligo-wa.caltech.edu/scirun/S6/DailyStatistics/"
  s5Link="http://blue.ligo-wa.caltech.edu/scirun/S5/DailyStatistics/"
  if time==None:
    return defaultURL
  if int(time) <= stopS5:
    return s5Link
  gpsTime=LIGOTimeGPS(time)
  Y,M,D,doy,h,m,s,ns,junk=xlaldate.XLALGPSToUTC(gpsTime)
  linkText="%s/%s/%s/"%(str(Y).zfill(4),str(M).zfill(2),str(D).zfill(2))
  outputLink=defaultURL+linkText
  return outputLink
#End getDailyStatsURL


class omega_config_parser(object):
	def __init__(self, file):
		self.data = []
		self.parse = False
		for l in open(file).readlines():
			#FIXME imposes that curly braces are on there own line, I don't know if that is required
			if l.strip() == "{":
				self.start_parse()
			if l.strip() == "}":
				self.end_parse()
			self.sanitize(l)
			self.parse_line(self.sanitize(l))

	def sanitize(self,l):
		return l.replace("'","").replace("[","").replace("]","")

	def start_parse(self):
		self.dict = {}
		self.parse = True

	def end_parse(self):
		self.data.append(copy.deepcopy(self.dict))
		self.parse = False

	def parse_line(self,l):
		if not self.parse: return
		line = l.split()
		# FIXME could fail if there are nulls
		if len(line) > 1: self.dict[line[0].rstrip(':')] = line[1:]

	def to_channel_dict(self):
		dict = {}
		for d in self.data:
			dict[d['channelName'][0]] = d
		return dict

	def _spec_name(self,dict):
		return '*%s_%.2f_spectrogram_whitened.png' % (dict['channelName'][0], float(dict['plotTimeRanges'][0]))

	def to_plot_tuple(self):
		out = []
		for d in self.data:
			out.append((d['channelName'][0], self._spec_name(d),''))
		return out

####################################################################
# Class to discover and get data for reconstructing the FOMs 
####################################################################
class getFOMdata:
  """
  This class will be a conglomerate of methods to rebuild custom
  version of the top N control room plots for H1 and L1.
  Eventually we try to integrate V1 also.
  """
  def __init__(self,verbose=False):
    """
    This object does not need any input arguments.
    """
    self.gpsTime=None
    self.ifoList=interferometers
    self.preWindow=int(10*3600)
    self.postWindow=int(2*3600)
    self.dataDict=dict()
    self.haveData=False
    self.__verbose__=verbose
    self.channelDict={"Inspiral Range":{
      "H1":{"frametype":"SenseMonitor_H1_M",
            "channels":["H1:DMT-SNSM_EFFECTIVE_RANGE_MPC.rms"]},
      "L1":{"frametype":"SenseMonitor_L1_M",
            "channels":["L1:DMT-SNSM_EFFECTIVE_RANGE_MPC.rms"]}
      },
                      "L1 STS 0.03-0.1Hz":{
      "L1":{"frametype":"STS_Blrms_M",
            "channels":["L1:DMT-BRMS_SEI_ETMX_STS2_X_0p03-0p1Hz.rms",
                        "L1:DMT-BRMS_SEI_ETMY_STS2_Y_0p03-0p1Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_X_0p03-0p1Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_Y_0p03-0p1Hz.rms"]}
      },
                      "L1 STS 0.1-0.35Hz":{
      "L1":{"frametype":"STS_Blrms_M",
            "channels":["L1:DMT-BRMS_SEI_ETMX_STS2_X_0p1-0p2Hz.rms",
                        "L1:DMT-BRMS_SEI_ETMY_STS2_Y_0p1-0p2Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_X_0p1-0p2Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_Y_0p1-0p2Hz.rms",
                        "L1:DMT-BRMS_SEI_ETMX_STS2_X_0p2-0p35Hz.rms",
                        "L1:DMT-BRMS_SEI_ETMY_STS2_Y_0p2-0p35Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_X_0p2-0p35Hz.rms",
                        "L1:DMT-BRMS_SEI_LVEA_STS2_Y_0p2-0p35Hz.rms"]}
      },
                      "L1 1-3Hz":{
      "L1_Seis":{"frametype":"Seis_Blrms_M",
                 "channels":["L0:DMT-BRMS_PEM_EY_SEISX_1_3Hz.rms",
                             "L0:DMT-BRMS_PEM_EY_SEISZ_1_3Hz.rms",
                             "L0:DMT-BRMS_PEM_LVEA_SEISZ_1_3Hz.rms"]},
      "L1_STS":{"frametype":"STS_Blrms_M",
                "channels":["L1:DMT-BRMS_SEI_ETMX_STS2_X_1-3Hz.rms",
                            "L1:DMT-BRMS_SEI_ETMY_STS2_Y_1-3Hz.rms",
                            "L1:DMT-BRMS_SEI_LVEA_STS2_X_1-3Hz.rms",
                            "L1:DMT-BRMS_SEI_LVEA_STS2_Y_1-3Hz.rms"]}
      },
                      "L1 Gurlap 3-10Hz":{
      "L1":{"frametype":"Seis_Blrms_M",
            "channels":["L0:DMT-BRMS_PEM_EX_SEISX_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_EY_SEISY_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_LVEA_SEISX_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_LVEA_SEISY_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_EX_SEISZ_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_EY_SEISZ_3_10Hz.rms",
                        "L0:DMT-BRMS_PEM_LVEA_SEISZ_3_10Hz.rms"]}
      },
                      "H1 0.03-0.3Hz":{
      "H1":{"frametype":"Seis_Blrms_M",
            "channels":["H0:DMT-BRMS_PEM_VAULT_SEISZ_0.03-0.1Hz.rms",
                        "H0:DMT-BRMS_PEM_EX_SEISZ_0.1_0.3Hz.rms",
                        "H0:DMT-BRMS_PEM_EY_SEISZ_0.1_0.3Hz.rms"]}
      },
                      "H1 3-10Hz":{
      "H1":{"frametype":"Seis_Blrms_M",
            "channels":["H0:DMT-BRMS_PEM_EX_SEISZ_3_10Hz.rms",
                        "H0:DMT-BRMS_PEM_LVEA_SEISZ_3_10Hz.rms",
                        "H0:DMT-BRMS_PEM_EY_SEISZ_3_10Hz.rms",
                        "H0:DMT-BRMS_PEM_VAULT_SEISZ_3_10Hz.rms"]}
      },
                      "H1 1-3Hz":{
      "H1":{"frametype":"Seis_Blrms_M",
            "channels":["H0:DMT-BRMS_PEM_EX_SEISZ_1_3Hz.rms",
                        "H0:DMT-BRMS_PEM_LVEA_SEISZ_1_3Hz.rms",
                        "H0:DMT-BRMS_PEM_EY_SEISZ_1_3Hz.rms",
                        "H0:DMT-BRMS_PEM_VAULT_SEISZ_1_3Hz.rms"]}
      },
                      "H1 10-30Hz":{
      "H1":{"frametype":"Seis_Blrms_M",
            "channels":["H0:DMT-BRMS_PEM_EX_SEISZ_10_30Hz.rms",
                        "H0:DMT-BRMS_PEM_LVEA_SEISZ_10_30Hz.rms",
                        "H0:DMT-BRMS_PEM_EY_SEISZ_10_30Hz.rms"]}
      }
                      }
    self.graphList=self.channelDict.keys()

  def getGraphKeys(self):
    """
    Simple method which returns a list of know graphs
    that can be generated.
    """
    return self.graphList
  
  def __getTimeSeries__(self,
                        frameType=None,
                        observatory=None,
                        channel=None,
                        start=None,
                        stop=None):
    """
    A method to simplify life to query ligo_data_find then open
    a frame time series object and hand that object back.  This method
    is meant to be internal to this class.
    """
    if (frameType==None) or \
           (channel==None) or \
           (start==None) or \
           (stop==None) or \
           (observatory == None):
      return None
    if len(observatory) > 1:
      observatory=observatory[0].upper()
    else:
      raise Exception, "Unexpected ifo option encountered, %s"%(observatory)
    #Just in case time stamps are backwards!
    if start>stop:
      tmp=start
      start=stop
      stop=tmp
    #Query for the data
    myEmptyCommand="""ligo_data_find --gaps --type=%s --observatory=%s \
--gps-start-time=%s --gps-end-time=%s --url-type=%s \
--lal-cache --no-proxy"""
    myCommand=myEmptyCommand%(frameType,\
                              observatory,\
                              start,\
                              stop,\
                              "file")
    if self.__verbose__:
      sys.stdout.write("Executing data query:%s\n"%myCommand)
    (errorCode,cmdOutput)=getstatusoutput(myCommand)
    if self.__verbose__:
      sys.stdout.write("Query Completed!\n")
    if errorCode != 0:
      if self.__verbose__:
        sys.stderr.write("Error seen! : Code %s: Message %s\n"%(errorCode,cmdOutput))
      return None
    if self.__verbose__:
      sys.stdout.write("Opening cache files.\n")
    cacheInRam=StringIO(cmdOutput)
    dataStream=frutils.FrameCache(lal.Cache.fromfile(cacheInRam))
    cacheInRam.close()
    if self.__verbose__:
      sys.stdout.write("Reading data stream.\n")
    dataVector=dataStream.fetch(channel,
                                start,
                                stop)
    if self.__verbose__:
      sys.stdout.write("Done reading in trend data.\n")
    return dataVector
                                  
  def setWindows(self,preWindow=None,postWindow=None):
    """
    Set a leading window in seconds and a following window in
    seconds for easy in graph centering.  Use integer seconds!
    """
    if preWindow!=None:
      self.preWindow=int(float(preWindow))
    if postWindow!=None:
      self.postWindow=int(float(postWindow))
      
  def setGPS(self,trigger=None):
    """
    Specify as a string (integer) a gps time of interest.
    """
    if trigger==None:
      raise Exception, "None variable specified as trigger time."
    self.gpsTime=str(int(trigger))

  def setGraphs(self,graphList=None):
      """
      Set the name of the graphs to create. To see the available
      name look at self.channelDict.keys()
      """
      if  graphList==None:
        raise Exception, "None variable specified at graphs to plot."
      for myKey in graphList:
        if myKey not in self.getGraphKeys():
          sys.stderr.write("Expected any of these %s, but got \
these %s\n"%(self.getGraphKeys(),graphList))
          raise Exception, "Invalid graph key specified."
      self.graphList=graphList

  def setIfos(self,myList=None):
    """
    Specify a list object with ifo names ie:
    ["L1",...,"H1"]
    """
    if myList==None:
      raise Exception, "None variable passed to method setIfos."
    for myIfo in myList:
      if not myIfo.strip().upper() in interferometers:
        raise Exception, "Invalid IFO specified, got %s\n"%myIfo
    self.ifoList=myList

  def getGraphList(self):
    """
    Method returns a list of keys to the know graphs to make.
    """
    return self.channelDict.keys()
    
  def getData(self):
    """
    Argue-less method that returns the data associated with
    the ifos and time in question.
    """
    if self.haveData:
      return self.dataDict
    else:
      self.queryGraphData()
      return self.dataDict

  def deleteGraphData(self):
    """
    Simple method to explicity kill the variable.
    """
    del self.dataDict
    self.dataDict=dict()
      
  def queryGraphData(self):
    """
    Method that does the appropriate data finding and
    reads the frame data into memory.
    """
    for myGraphName in self.channelDict.keys():
      self.dataDict[myGraphName]=dict()
      for graphIFO in self.channelDict[myGraphName].keys():
        graphFrameType=self.channelDict[myGraphName][graphIFO]["frametype"]
        for channel in self.channelDict[myGraphName][graphIFO]["channels"]:
            myStart=int(float(self.gpsTime)) - self.preWindow
            myStop =int(float(self.gpsTime)) + self.postWindow
            try:
              self.dataDict[myGraphName][channel]=self.__getTimeSeries__(graphFrameType,
                                                                       graphIFO,
                                                                       channel,
                                                                       myStart,
                                                                       myStop)
            except:
              sys.stderr.write("Error accessing frame data. Skipping \
trace for graph %s.\n"%(myGraphName))
              sys.stderr.write("Channel    :%s\n"%channel)
              sys.stderr.write("Frame Type :%s\n"%graphFrameType)
              sys.stderr.write("GPS Start  :%s\n"%myStart)
              sys.stderr.write("GPS Stop   :%s\n"%myStop)
