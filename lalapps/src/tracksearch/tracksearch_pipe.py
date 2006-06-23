#!/usr/bin/env python
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
(options,args) = parser.parse_args()

iniFile=os.path.abspath(str(options.iniFile))
segmentList=os.path.abspath(str(options.segmentList))

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
    
#We assume the input segment list is > 100 so we will try to loop it
dataBlockSize=int(float(str.strip(cp.get('layerconfig','layerTopBlockSize'))))

#Convert the segment list to smaller blocks
reformatSegList=tracksearch.tracksearchConvertSegList(segmentList,dataBlockSize)
#Setup logfile mask
logFilePath=cp.get('fileLayout','workpath')
logFileMask=logFilePath+'/logFile_'
if not os.path.isdir(logFilePath):
            try: os.mkdir(logFilePath)
            except: pass
reformatSegList.writeSegList()
segmentList=reformatSegList.getSegmentName()
allData=pipeline.ScienceData()
allData.read(segmentList,dataBlockSize)
allData.make_chunks(dataBlockSize)
# Loop over our new collection of 100s chunks to set up each DAG
indexA=0
for block in allData:
    """This code is only for 1 100sec stretch"""
    indexA+=1
    logfile=logFileMask+str(indexA)+'.log'
    tracksearchBlock=tracksearch.tracksearch(cp,block,logfile)
    tracksearchBlock.createJobs()
    tracksearchBlock.writePipelineDAG()
        
print "We prepared "+str(allData.__len__())+" search DAGs."
sys.exit(0)


