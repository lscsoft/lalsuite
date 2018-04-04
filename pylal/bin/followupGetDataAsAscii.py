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
__author__ = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__   = 'followupGetDataAsAscii.py'
"""
This script lets you specify a channel and gps markers, from this a
ascii file is generated and written to disk.  This is a simple tool
to allow novice persons to analyse data without understanding the
complete LIGO data handling infrastructure.
"""
from StringIO import StringIO
from commands import getstatusoutput
from glue import lal
from numpy import mean,std,median,min,max,abs,array
from pylal import frutils
from sys import stderr,stdout,exit
import optparse

usage = """usage: %prog [options]"""
validIFOList=["L","H","V"]
parser = optparse.OptionParser(usage)
#Specify all needed options
parser.add_option("-c","--channel-name",action="store",type="string",\
                   metavar="CHANNELNAME",default=None,\
                   help="Specify the channel name you want to fetch \
data for.  This option can only take one channel at a time.")
parser.add_option("-s","--gps-start",action="store",type="string",\
                   metavar="gpsStart",default=None,\
                   help="Specify the GPS integer start time of the \
data interval of interest")
parser.add_option("-e","--gps-end",action="store",type="string",\
                   metavar="gpsEnd",default=None,\
                   help="Specify the GPS integer end time of the \
data interval of interest.")
parser.add_option("-d","--gps-delta",action="store",type="string",\
                   metavar="gpsDeltaT",default=None,\
                   help="Specify a relative offset from specified \
GPS start time.  This argument takes precidence over specifying \
the explict gps end of the interval of interest.")
parser.add_option("-i","--ifo",action="store",type="string",\
                   metavar="IFO",default=None,help="Specify the \
IFO of interest as one of the following from the set of valid IFO \
args %s"%validIFOList)
parser.add_option("-t","--type",action="store",type="string",\
                   metavar="FRAMETYPE",default=None,\
                   help="Specify the type of frame to make a query \
for.  See the ligo_data_find help for a complete list of \
frame types.")
parser.add_option("-u","--url-type",action="store",type="string",\
                   metavar="URLTYPE",default="file",\
                   help="The default URL type to query for is local \
file.  If you wish to change it to another availble type see \
the ligo_data_find help for the available types.")
parser.add_option("-o","--output-filename",action="store",\
                   type="string",metavar="FILENAME",default=None,\
                   help="Set this flag to override the automagic \
naming of files to something other than channel_gpsstart_deltaT.ascii")
parser.add_option("-v","--verbose",action="store_true",\
                   default=False,\
                   help="Set this flag to make the software \
output progress information to the shell.")
#
# Process command arguments
#
(opts,args) = parser.parse_args()
if (opts.channel_name == None) or \
   (opts.url_type == None) or \
   (opts.type == None) or \
   (opts.gps_start == None) or \
   ((opts.gps_end == None) and (opts.gps_delta == None)):
    stderr.write("Invalid arguments passed or insufficient \
arguments given!\n")
    exit(1)
#
# Seek out data
#
myEmptyCommand="""ligo_data_find --gaps --type=%s --observatory=%s \
--gps-start-time=%s --gps-end-time=%s --url-type=%s \
--lal-cache --no-proxy"""
myStartTime=int(opts.gps_start)
if opts.gps_delta:
    myEndTime=myStartTime+int(opts.gps_delta)
elif opts.gps_end:
    myEndTime=int(opts.gps_end)
else:
    stderr.write("Error specifying the GPS interval \
options --gps-end or --gps-delta.\n")
myFullCommand=myEmptyCommand%(opts.type,\
                              opts.ifo,\
                              myStartTime,\
                              myEndTime,\
                              opts.url_type)
if opts.verbose:
    stdout.write("Querying for data location.\n")
    stdout.write("Query Used:%s\n"%myFullCommand)
(errorCode,cmdOutput)=getstatusoutput(myFullCommand)
if errorCode != 0:
    stderr.write("Error querying for data location!\n")
    stderr.write("%s\n%s"%(errorCode,cmdOutput))
    exit(1)
elif opts.verbose:
    stdout.write("Query completed.\n")
#
# Load data into memory for reparsing
#
if opts.verbose:
    stdout.write("Reading in frame data for channel %s.\n"%opts.channel_name)
memoryFP=StringIO(cmdOutput)
dataStream=frutils.FrameCache(lal.Cache.fromfile(memoryFP))
memoryFP.close()
dataVector=dataStream.fetch(opts.channel_name,myStartTime,myEndTime)
if opts.verbose:
    stdout.write("Done reading in data.\n")
#
# Write out metadata to ascii header
#
if opts.verbose:
    stdout.write("Writing metadata and data to file.\n")
if opts.output_filename:
    myFilename=str(opts.output_filename)
else:
    myFilename="%s_%s_%s.ascii"%(opts.channel_name,\
                                 myStartTime,\
                                 myEndTime-myStartTime)
fp=open(myFilename,'w')
myMetaDataDict=dataVector.metadata.todict()
for myKey,myVal in myMetaDataDict.iteritems():
    fp.write("#%s:%s\n"%(myKey,myVal))
fp.write("#t0:%s\n"%(myStartTime))
#
# Write out data
#
for datum in dataVector.tolist():
    fp.write("%s\n"%datum)
#
# Close out file and
#
fp.close()
if opts.verbose:
    stdout.write("Data written to file and file closed.\n")
