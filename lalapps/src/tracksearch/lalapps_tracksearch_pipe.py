#
# Copyright (C) 2004, 2005 Cristina V. Torres
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA
#
"""
This is a python script that will setup a big segment run with given
arguements and a combination from the ini file specified
"""
import sys
import os
import getopt
from optparse import OptionParser
import ConfigParser
import time
from glue import pipeline
try:
    from lalapps import tracksearch
except:
    import tracksearch
try:
    from lalapps.tracksearchutils import progressSpinner
except:
    from tracksearchutils import progressSpinner
try:
    from lalapps.tracksearchutils import autoCreateSegmentList
except:
    from tracksearchutils import autoCreateSegmentList
try:
    from lalapps.tracksearchutils import determineDataPadding
except:
    from tracksearchutils import determineDataPadding
    

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
parser.add_option("-r","--remove_pipeline",dest="removepipeline",
                  default=False,
                  action="store_true",
                  help="This option is used to remove an installed pipeline with the specified configuration INI and segment file.  It deciphers the location of where the pipes were installed and removes the contents and directories. Use this option CAUTIOUSLY!"
                  )
parser.add_option("-b","--rubber_block_length",dest="topBlockFloat",
                  default=False,
                  action="store_true",
                  help="If this flag is specified then don't adjust the entries in the segment list, instead prepare a simple 1 layer investigation with no interval restriction on the layerTopBlockSize field in the config file.  We have an implicit minimum block size to use which is equal to a single tracksearchTime job call given the configuration options in INI file.  If the INI file contains configuration options for more than just a 1 layer run then pipe builder will exited with an error. We build chunks equal to the length of layerTopBlockSize plus smaller ones to maximize the data we check.  The lower bound on chunk size is a 1 tracksearchTime job worth of data."
                  )
parser.add_option("-u","--unique_path",dest="uniquePath",
                  default='',
                  help="Mechanism to reuse a ini file where UNIQUE_PATH is appended to paths from work_path log_path and dag_path, these options are found in [filelayout] section of ini file.")
parser.add_option("-v","--verbose",dest="verbose",
                  help="Invoke a more verbose setting with progress indicator.",
                  action="store_true",
                  default=False
                  )
parser.add_option("-e","--specific_event",dest="specificGPS",
                  default='0',
                  help="Mechanism to setup a tracksearch pipeline for a specific GPS time stamp.  This method uses the INI file specified to setup an appropriate segment list for the run such that the GPS time specified is approximately centered in the middle of the data to be analyzed.  The amount of surrounding data (durations) etc are all defined in the INI file. Using this option will implicity set the flag associated with --unique_path unless it has already been set to something by the user.",
                  )
parser.add_option("-d","--override_burn",dest="overrideBurn",
                  default=False,
                  action="store_true",
                  help="Specify this flag to override the burn data option in ini file.   You may want to invoke this if the input segment list is the left over data from a previous pipe construction.  The previous pipe built a data used list with the edges already burned.  This does not affect the actual analysis pipeline.  It still uses the data padding (burn durations).  This keeps from discarding good data multiple time from same segment."
                  )

(options,args) = parser.parse_args()

iniFile=os.path.abspath(os.path.expanduser(str(options.iniFile)))
segmentList=os.path.abspath(os.path.expanduser(str(options.segmentList)))
injectFileFlag=bool(options.invokeInject)
removepipeline=bool(options.removepipeline)
verbose=bool(options.verbose)
overrideBurn=bool(options.overrideBurn)
topBlockFloat=bool(options.topBlockFloat)
uniquePath=str(options.uniquePath)
specificGPS=float(options.specificGPS)

if not os.path.exists(iniFile):
    print 'Can not find iniFile: ',os.path.basename(iniFile)
    os.abort()
if (not os.path.exists(segmentList) and specificGPS == 0):
    print 'Can not find segentlist: ',os.path.basename(segmentList)
    os.abort()
    
#Read in the configuration file options
cp = ConfigParser.ConfigParser()
cp.read(iniFile)

#If specificGPS flag is used to start pipe builder
#It is placed in same path as INI file can be found also.
if (specificGPS != 0):
    segmentList=autoCreateSegmentList(cp,iniFile,specificGPS)
    if uniquePath == '':
        uniquePath='/event_'+str(specificGPS)

#Adjust the paths using unique_path option
newWorkPath=os.path.normpath(os.path.expanduser(str(cp.get('filelayout','workpath')))+"/"+uniquePath)
newDagPath=os.path.normpath(os.path.expanduser(str(cp.get('filelayout','dagpath')))+"/"+uniquePath)
newLogPath=os.path.normpath(os.path.expanduser(str(cp.get('filelayout','logpath')))+"/"+uniquePath)
cp.set('filelayout','workpath',newWorkPath)
cp.set('filelayout','dagpath',newDagPath)
cp.set('filelayout','logpath',newLogPath)

#Use config information to remove files if requested
if removepipeline:
    print "Removing this particular pipeline setup!"
    #print "Waiting 3 seconds..."
    #time.sleep(4)
    pathList=[]
    for option in cp.options('filelayout'):
        workpath=cp.get('filelayout',option)
        if os.path.isdir(workpath):
            print "Removing:",workpath
            print "Feature not available yet."
            pathList.append(workpath)
    print "For now cut and paste the following."
    for entry in pathList:
        print "rm -rf ",entry
    sys.exit(0)
    
#Run iniFile partial option sanity check
testIniFile=tracksearch.tracksearchCheckIniFile(cp)
#Find Errors
testIniFile.checkOpts()
#Print Errors Found
if testIniFile.numberErrors() > 0 :
    testIniFile.printErrorList()
    os.abort()
#Print Warnings
if testIniFile.numberWarnings() > 0:
    testIniFile.printWarningList()

#Print memory use estimate!
memoryEstimate=testIniFile.getMemoryUseEstimate()
psdSmoothBandEstimate=testIniFile.getBandwithSmoothedEstimate()
print "***"
print "For this dag configuration the peak memory use estimate is: "+str(memoryEstimate)+" MBs."
print "As a rule of thumb this value should be about one third of the systems available memory or less."
print "***"
if verbose:
	print "For this dag setup if smoothing is used we estimate that a band of "+str(psdSmoothBandEstimate)+" Hz will be smoothed during the whitening process."

#Check ini file for an injection block heading.
if (injectFileFlag == True):
    print "Checking injection options."
    if not(testIniFile.hasInjectSec()):
        print "Double check the parameter file's injection section!"
        os.abort()

#We assume the input segment list has entries exceeding layerTopBlockSize
#so we will try to loop it.  If the pipe builder was invoked with a FLOATING
#top block size then we will issue an error IFF there is more than 1 layer configured
segmentListName=segmentList
dataBlockSize=int(float(str.strip(cp.get('layerconfig','layerTopBlockSize'))))
if not(topBlockFloat):
    #Convert the segment list to smaller blocks
    reformatSegList=tracksearch.tracksearchConvertSegList(segmentList,dataBlockSize,cp,topBlockFloat,overrideBurn)
    reformatSegList.writeSegList()
    segmentListName=reformatSegList.getSegmentName()
    allData=pipeline.ScienceData()
    allData.read(segmentListName,dataBlockSize)
    allData.make_chunks(dataBlockSize)
else:
    #Do optimized floating blocks
    #Check for layer2
    setSize=int(cp.get('layerconfig','layer1SetSize'))
    timeScale=float(cp.get('layerconfig','layer1TimeScale'))
    minSize=setSize*timeScale
    print "Building pipe with rubber block size option enabled."
    print "Minimum duration block: "+str(minSize)
    print "Maximum duration block: "+str(dataBlockSize)
    for opt in cp.options('layerconfig'):
        if str(opt).lower().__contains__(str('layer2TimeScale').lower()):
            print "Error found additional layerconfig options for multi-resolution search."
            print "Aborting pipeline construction."
            os.abort
    rubberSegList=tracksearch.tracksearchConvertSegList(segmentList,minSize,cp,topBlockFloat,overrideBurn)
    rubberSegList.writeSegList()
    rubberSegListName=rubberSegList.getSegmentName()
    allData=pipeline.ScienceData()
    allData.read(rubberSegListName,minSize)
    allData.make_optimised_chunks(0,dataBlockSize)
    
#Setup logfile mask
logFilePath=cp.get('filelayout','logpath')
logFileMask=logFilePath+'/logFile_'
if not(os.path.isdir(logFilePath)):
       print 'Log path does not exist!'
       print 'Expected to find:',logFilePath
       tracksearch.buildDir(logFilePath)
       print 'Created...'

indexA=0
# The pipeline needs a method change for handling large segment lists
# We should switch from 1 DAG per analysis segment to
# a single DAG with the proper relationships to analyze each segment
# each segment should be logically independant and thus one dag
# can handle a search on a limitless amount of data using a single DAG
# to manage the workflow.  The results of each segment should still be
# stored independently

#Write out the data chunks actually configured for search.
print "Writing the actual segment list to disk for reference."
tracksearch.writeChunkListToDisk(allData,str(segmentListName+".dataUsed"))

##### New initialization method
tracksearchLogFile="tracksearchDag.log"
tracksearchBlock=tracksearch.tracksearch(cp,injectFileFlag,tracksearchLogFile)
progress=progressSpinner(verbose,2)
progress.setTag('Building ')
if allData.__len__() == 0:
    print "Please check your segment list."
    print "Also check INI files option in section [layerconfig]"
    os.abort()
else:
    print "Building :"+str(allData.__len__())+" blocks into pipe."
for block in allData:
    progress.updateSpinner()
    tracksearchBlock.createJobs(block)
progress.closeSpinner
tracksearchBlock.writePipelineDAG()
        
print "\n We prepared "+str(allData.__len__())+" analysis segments."
print "The dag info should be found in :",cp.get('filelayout','dagpath')
sys.exit(0)


 
