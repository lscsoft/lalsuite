#!/usr/bin/env python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
__author__ = 'Cristina Valeria Torres <cristina.torres@ligo.org>'

import sys
import numpy
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
from StringIO import StringIO

#####################################################
## Function related to handling and processing of  ##
## frame data manipulation as part of followup     ##
## procedures                                      ##
#####################################################

class connectToFrameData:
  """
  A class that uses classes defined in pylal.frutils.
  This handles searching for data via ligo_data_find
  and creating an object which can fetch specific
  channels from that data epoch given a frame type.
  """
  def __init__(self,\
               gpsStart=None,\
               gpsEnd=None,\
               frameType="R",\
               observatory=None,\
               urlType="local",\
               noGaps=True,\
               verbose=False):
      self.__emptyquery__="""ligo_data_find --gaps \
      --type=%s --observatory=%s --gps-start-time=%s \
      --gps-end-time=%s --url-type=%s  --lal-cache --no-proxy"""
      if gpsStart == None or gpsEnd == None:
          print "Error with gps arguments given are NoneType!"
      self.globalStart=int(numpy.floor(gpsStart))
      self.globalEnd=int(numpy.ceil(gpsEnd))
      self.observatory=observatory
      self.frameType=frameType
      self.urlType=urlType
      self.query=self.__emptyquery__%(self.frameType,\
                                      self.observatory,\
                                      self.globalStart,\
                                      self.globalEnd,\
                                      self.urlType)
      #
      # Create a connection to the data
      #
      (errorCode,cmdOutput)=getstatusoutput(self.query)
      self.dataStream=None
      if errorCode != 0:
        sys.stderr.write("%s\n"%cmdOutput)
        raise Exception, "Error access ligo_data_find"
      memFP=StringIO(cmdOutput)
      self.dataStream=frutils.FrameCache(lal.Cache.fromfile(memFP))
      memFP.close()

  def getDataStream(self,\
                       channel=None,\
                       gpsStart=None,\
                       gpsEnd=None,\
                       verbose=None):
      """
      Fetch data inside the bounds specified during
      the init call to this class. If gpsStart and gpsEnd
      are None then we fetch data around gpsZero with
      a window of pm gpsWindow.
      """
      if gpsStart == None:
        gpsStart=self.globalStart
      if gpsEnd == None:
        gpsEnd=self.globalEnd
      gpsStart=int(numpy.floor(gpsStart))
      gpsEnd=int(numpy.ceil(gpsEnd))
      if channel == None:
        raise Exception, "Channel name passed to method as None"
      return  self.dataStream.fetch(channel,gpsStart,gpsEnd)

  def saveAsTextDataStream(self,\
                           channel=None,\
                           gpsStart=None,\
                           gpsEnd=None,\
                           filename=None,
                           verbose=None):
    """
    Write the data from the stream to an ascii text file
    """
    if gpsStart==None:
      gpsStart=self.globalStart
    if gpsEnd==None:
      gpsEnd=self.globalEnd
    if filename==None:
      filename="default_data_file.txt"
    if channel==None:
      raise Exception, "No channel name given to method."
    tmpData=self.getDataStream(channel,gpsStart,gpsEnd)
    tmpDataDict=tmpData.metadata.todict()
    myFP=file(filename,'w')
    for myKey,myVal in tmpMetaData.iteritems():
                  myFP.write("#%s:%s\n"%(myKey,myVal))
    myFP.write("#t0:%s\n"%(gpsStart))
    for datum in tmpData.tolist():
      myFP.write("%s\n"%datum)
    myFP.close()

  def convertFrObjectToXYLists(self,dataStream):
    """
    Returns two lists as (tStamps,dataPoints) the method
    expects you to give is an object FrameCache from
    frutils.FrameCache(xxxx).fetch() or use this getDataStream()
    method to get a variable of this type.
    """
    metaDataDict=dataStream.metadata.todict()
    dT=float(metaDataDict["dt"])
    segmentInfo=metaDataDict["segments"]
    dataPoints=dataStream.tolist()
    myStart=segmentInfo[0][0]
    myStop=myStart+len(dataPoints)*dT
    tStamps=arange(myStart,myStop,dT)
    return (tStamps,dataPoints)



