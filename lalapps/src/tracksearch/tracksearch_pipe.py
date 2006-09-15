#!/usr/bin/env python2.3
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
                  default="",
                  metavar="injectFile")

(options,args) = parser.parse_args()

iniFile=os.path.abspath(str(options.iniFile))
segmentList=os.path.abspath(str(options.segmentList))
injectFile=str(options.invokeInject)
if not os.path.exists(iniFile):
    print 'Can not find iniFile: ',os.path.basename(iniFile)
    os.abort()
if not os.path.exists(segmentList):
    print 'Can not find segentlist: ',os.path.basename(segmentList)
    os.abort()
if injectFile != "":
    if not os.path.exists(injectFile):
        print 'Can not find the file to inject: ',os.path.basename(injectFile)
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
if ((injectFile != "") and (testIniFile.injectSecTest() == False)):
    print "Ini file doesn't have an injection section!"
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
    """This code is only for 1 100sec stretch"""
    indexA+=1
    logfile=logFileMask+str(indexA)+'.log'
    tracksearchBlock=tracksearch.tracksearch(cp,block,injectFile,logfile)
    tracksearchBlock.createJobs()
    tracksearchBlock.writePipelineDAG()
        
print "We prepared "+str(allData.__len__())+" search DAGs."
sys.exit(0)


