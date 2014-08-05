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
from lalapps import tracksearch
try:
    from lalapps.tracksearch import buildDir
except:
    from tracksearch import buildDir
try:
    from lalapps.tracksearchutils import progressSpinner
except:
     from tracksearchutils import progressSpinner
try:
    from lalapps.tracksearchutils import generateFileList
except:
    from tracksearchutils import generateFileList


#Parse in command line inputs
# 
# Use an INI file
#  Ini file from the original pipeline configuration options
#  A new custom ini file with only cut related sections
#  All possible INI sections for cutting are
#  [graphicsthreshold],[candidatethreshold],[pylibraryfiles],[condor]
#
# We allow overriding these options via command line or defaulting via 
# hardware.
#

# Dag builder method
# Sub builder method

parser = OptionParser()
parser.add_option("-f","--file",dest="iniFile",
                  help="This is the ini file required to create a DAG for thresholding a batch of candidate files.  You can borrow the original search INI file or use a custom INI file specific to the way you want to cut your candidates files.",
                  default="tracksearchCut.ini"
                  )
parser.add_option("-s","--single_file_list",dest="singleList",
                  help="This is a plain text file with each candidate list one per line using an absolute path or path relative to the execution location for lalapps_tracksearch_cutpipe. To cut a single file, just specify the candidate file using this same flag, the code will see if this is a list of files or a single file to process.",
                  default="candidateFileList.txt"
                  )
parser.add_option("-b","--build_figures",dest="buildFigures",
                  help="Setting this flag will create a set of histograms which should be included INI file.  If this option is not present in the INI file then an error is returned unless we use the option --figure_override (Override not yet implemented)",
                  default=False,action="store_true"
                  )
parser.add_option("-p","--plot_triggers",dest="plotTriggers",
                  help="Setting this flag will cause the resulting candidate object or file to be used to create a trigger plot (line drawing) of the triggers which meet the threshold criteria.",
                  default=False,action="store_true"
                  )
parser.add_option("-o","--output_name",dest="outputName",
                  help="Setting this will create a dag via the form  MYPIPE.dag and a matching submit file via MYJOB.sub. These two files will be placed in the current working directory. These files should be launched from a local disk area, but this can be overridden with the appropriate condor flag.",
                  default="MYPIPE"
                  )
parser.add_option("-u","--output_path",dest="outputPath",
                  default="./",
                  help="Setting this will set the output the cut pipe to write the resulting candidate files and figures to the path specified by OUTPUTPATH, creating it if needed.  If you don't set this then it is assumed that you will want OUTPUTPATH to be instead in the current directory and named MYPIPE.RESULTDIR"
                  )
parser.add_option("-d","--dag_locks",dest="dagLockPath",
                  default="/tmp/dagLocks",
                  help="This is an optional location (local disk) for the lock files Condor needs to keep the log files on a non-local disk.  The default is /tmp/dagLocks"
                  )

(options,args)=parser.parse_args()
dagLocks=str(os.path.normpath(options.dagLockPath))
buildDir(dagLocks)
outputName=str(os.path.normpath(options.outputName))
outputPath=str(os.path.normpath(options.outputPath))
outputResultsPath=str(os.path.normpath(outputPath+"/"+outputName+".RESULTS"))
singleList=os.path.normpath(options.singleList)
buildFigures=options.buildFigures
plotTriggers=options.plotTriggers
iniFile=os.path.normpath(options.iniFile)


listOfFiles=generateFileList(singleList)

if not os.path.exists(iniFile):
    print 'Can not find iniFile: ',os.path.basename(iniFile)
    os.abort()
if not os.path.exists(singleList):
    print 'Can not find segentlist: ',os.path.basename(singleList)
    os.abort()

#Load up the iniFile
cp= ConfigParser.ConfigParser()
cp.read(iniFile)

#Setup DAG
outputPath=os.path.abspath(os.path.normpath(outputPath))
dagLog=os.path.normpath(dagLocks+"/"+outputName+".LOG")
myDag=pipeline.CondorDAG(os.path.normpath(dagLog))
myDag.set_dag_file(os.path.normpath(str(outputName)))

#Setup SUB
tsHandler=os.path.expanduser(cp.get('condor','clustertool'))
tsUniverse=str(cp.get('condor','clustertool_universe')).lower()
if not os.path.exists(str(tsHandler)):
    print "ERROR: Can't find tracksearch handler executable."
    os.abort()
    

myJob=pipeline.CondorDAGJob(tsUniverse,tsHandler)
myJob.set_sub_file(str(outputName)+".sub")
logDir=os.path.normpath(outputPath+"/logs/")
buildDir(logDir)
myJob.set_stdout_file(os.path.abspath(os.path.normpath(logDir+"/log_$(cluster)_$(process).out")))
myJob.set_stderr_file(os.path.abspath(os.path.normpath(logDir+"/log_$(cluster)_$(process).err")))


#All all the proper options and macro handles
if not cp.has_section('candidatethreshold'):
    print "NO [candidatethreshold] section!\n"
    os.abort()
else:
    if cp.has_option('candidatethreshold','expression-threshold'):
        newVal=cp.get('candidatethreshold','expression-threshold')
        if (newVal.__contains__('"') and not newVal.__contains__('\\')):
            newVal=str(newVal).replace('"','\\"')
        cp.set('candidatethreshold','expression-threshold',newVal)
    myJob.add_ini_opts(cp,'candidatethreshold')

if not (not buildFigures and not plotTriggers):
    if not cp.has_section('graphicsthreshold'):
        print "NO [graphicsthreshold] section!\n"
        os.abort()
    else:
        if not plotTriggers:
            cp.remove_option('graphicsthreshold','line-plot')
            cp.remove_option('graphicsthreshold','timescale')
        if not buildFigures:
            cp.remove_option('graphicsthreshold','histogram')
            cp.remove_option('graphicsthreshold','trigger_property')
        if not plotTriggers and not buildFigures:
            cp.remove_option('graphicsthreshold','print')
        myJob.add_ini_opts(cp,'graphicsthreshold')

#Add options for filename to process
myJob.add_opt('file','$(macroFileToProcess)')

#Use condor file transfer mechanism to move candidate files into
#position
myJob.add_condor_cmd('transfer_input_files','$(macroFileList)')
myJob.add_condor_cmd('when_to_transfer_output','on_exit')
myJob.add_condor_cmd('initialdir',outputResultsPath)
buildDir(outputResultsPath)

#Setup dag nodes
#Loop over files to process
if not cp.has_section('pylibraryfiles'):
    print "NO [pylibraryfiles] section!\n"
    os.abort()
else:
    libraryFile=os.path.expanduser(cp.get('pylibraryfiles','pyutilfile'))
    if not os.path.exists(str(libraryFile)):
                         print "ERROR: Library file not found."
                         os.abort()

for thisFile in listOfFiles:
    myFile=thisFile.strip("\n")
    thisNode=pipeline.CondorDAGNode(myJob)
    thisNode.add_macro('macroFileList',str(myFile)+","+libraryFile)
    if tsUniverse == 'local':
        thisNode.add_macro('macroFileToProcess',myFile)
    else:
        thisNode.add_macro('macroFileToProcess',os.path.basename(myFile))
    myDag.add_node(thisNode)

#Write out the files that constitute the dag
myDag.write_sub_files()
myDag.write_dag()
