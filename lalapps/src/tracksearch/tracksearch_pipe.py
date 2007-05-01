#!/usr/bin/env python
"""
  * Copyright (C) 2004, 2005 Cristina V. Torres
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with with program; see the file COPYING. If not, write to the
  *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  *  MA  02111-1307  USA
"""
"""
This is a python script that will setup a big segment run with given
arguements and a combination from the ini file specified
"""
import sys,os
import getopt

from optparse import OptionParser
import ConfigParser
import time
from glue import pipeline
import tracksearch

#Parse the command line options and create help features
parser = OptionParser()
parser.add_option("-f","--file",dest="iniFile",
                  help="The path to the tracksearch initialization file to use for this run",
                  default='tracksearch.ini',
                  metavar="FILE")
parser.add_option("-s","--segment_file",dest="segmentList",
                  help="The filename of the segments list created by segwizard",
                  default='segments.txt',
                  metavar="FILE")
parser.add_option("-i","--invoke_inject",dest="invokeInject",
                  help="This option will create a pipeline with injections invoked if there is a section in the config file for doing injections. You must still specify name and path to 2 column text data to inject.",
                  default=False,action="store_true",
                  metavar="injectFile")

(options,args) = parser.parse_args()

iniFile=os.path.abspath(str(options.iniFile))
segmentList=os.path.abspath(str(options.segmentList))
injectFileFlag=bool(options.invokeInject)
if not os.path.exists(iniFile):
    print 'Can not find iniFile: ',os.path.basename(iniFile)
    os.abort()
if not os.path.exists(segmentList):
    print 'Can not find segentlist: ',os.path.basename(segmentList)
    os.abort()
#Read in the configuration file options
cp = ConfigParser.ConfigParser()
cp.read(iniFile)

#Run iniFile partial option sanity check
testIniFile=tracksearch.tracksearchCheckIniFile(cp)
#Find Errors
testIniFile.checkOpts()
#Print Errors Found
if testIniFile.numberErrors() > 0 :
    testIniFile.printErrorList()
    os.abort()

#Check ini file for an injection block heading.
if (injectFileFlag == True):
    print "Checking injection options."
    if not(testIniFile.hasInjectSec()):
        print "Double check the paramter file's injection section!"
        os.abort()

#We assume the input segment list is > 100 so we will try to loop it
dataBlockSize=int(float(str.strip(cp.get('layerconfig','layerTopBlockSize'))))

#Convert the segment list to smaller blocks
reformatSegList=tracksearch.tracksearchConvertSegList(segmentList,dataBlockSize)
#Setup logfile mask
logFilePath=cp.get('filelayout','logpath')
logFileMask=logFilePath+'/logFile_'
if not(os.path.isdir(logFilePath)):
       print 'Log path does not exist!'
       print 'Expected to find:',logFilePath
       tracksearch.buildDir(logFilePath)
       print 'Created...'

reformatSegList.writeSegList()
segmentList=reformatSegList.getSegmentName()
allData=pipeline.ScienceData()
allData.read(segmentList,dataBlockSize)
allData.make_chunks(dataBlockSize)
# Loop over our new collection of 100s chunks to set up each DAG
indexA=0
for block in allData:
    indexA+=1
    logfile=logFileMask+str(indexA)+'.log'
    tracksearchBlock=tracksearch.tracksearch(cp,block,injectFileFlag,logfile)
    tracksearchBlock.createJobs()
    tracksearchBlock.writePipelineDAG()
        
print "We prepared "+str(allData.__len__())+" search DAGs."
if allData.__len__() == 0:
    print "Please check your segment list."
    os.abort()
sys.exit(0)


