#!/usr/bin/env python2.3

__author__ = 'Charlie Torres <charlie@phys.utb.edu>'
__date__ = '$Date$'
__version__ = ''

from optparse import OptionParser
import ConfigParser
import getopt
import math
import os
import string
import sys
import time
from candidateUtils import *

#Begin MAIN part of this python code to actually manipulate the curve
#lists


#################### Internal Methods #################

def buildCandidateGlob(fileList):
    #Build glob var from filelist.  Memory friendly
    #Returns copy of results
    canList=fileList
    if canList.__len__() == 1:
        print "Globbing is not possible only one file found!"
        print "Returning the result as that entry!"
        tmpObject=candidateList()
        tmpObject.loadfile(canList.pop(0))
        newCandidateObject=copy.deepcopy(tmpObject)
    else:
        newCandidateObject=candidateList()
        newCandidateObject.loadfile(canList.pop(0))
    for entry in canList:
            print " "
            print "Loading candidate list file: ",entry
            tmpCandidate=candidateList()
            tmpCandidate.loadfile(entry)
            print "Globing file:",entry
            newCandidateObject=newCandidateObject.globList(tmpCandidate,True)
            del tmpCandidate
    return copy.deepcopy(newCandidateObject)
#End method

def gnuplotScriptFile(filename):
    """
    Invokes gnuPlot in an unusual hackish was to create a PS scatter
    plot to view the pixels found as part of a curve(s). Add a .plt
    extension to the filename specified
    """
    txtTemplate='plot "%s"\n set size 1.0,0.6\n set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 14\n set output "%s.ps"\n replot\n set terminal x11\n set size 1,1\n'
    output_fp=open(filename+'.plt','w')
    output_fp.write(txtTemplate%(filename,filename))
    output_fp.close()
#End method gnuplotScriptFile()

###################End internal methods #################

parser = OptionParser()

parser.add_option("-f","--file",dest="filename",
                  default="",
                  help="This specifies the list of curve files to process contained in FILE, either rethresholding, globbing.  If the input is a directory then it is assumed that we will be working on all the seen files.  /PATH/*.ext is allowed.  Only one of the other following options should be specified at run-time.  Combining the following options will yield unpredictable results!",
                  metavar="FILE"
                  )
parser.add_option("-g","--glob",dest="glob",
                  default=False,
                  action="store_true",
                  help="This option tells the code to glob all the lists together into one larger curve candidate listing.  The output be a single or sequence of files with GLOB:Start:00000000,0:Stop:00000000,0:TF:000:000:.candidates. This option can not be used with clobber",
                  )
parser.add_option("-o","--outfile",dest="outfile",
                  default="",
                  help="This overrides the automagic output file nameing characteristics of this utility.  Please invoke the option to write to the output to a definite file location of your choice."
                  )
parser.add_option("-c","--clobber",dest="clobberFilename",
                   default="",
                   help="This is a file or directory/mask of all files that should be globbed together then checked against the globbed list from --file option.  Entries candidates in FILE with reside in list CLOBBER will be deleted from list FILE and written to disk along with an globbed list which is not clobbered.The filename will be like CLOBBER:Start:00000000,0:Stop:00000000,0:TF:000:000:AND:Start:00000000,0:Stop:00000000,0:TF:000:000:.candidates.  This option can not be used with glob, we assume specified --file are globs of candidates.   We glob the --clobber option if neccessary to perform the clobber."
                   )
parser.add_option("-t","--threshold",dest="thresholdString",
                  default="",
                  help="This is the the thresholding options that can be done to an already complete candidate list.  We can reprocess the list via the following syntax, 1.0,>,and,=,7  so this would be written like -t 1.0,>,and,=,7 which would logically translated to: The power for the curve should be greater than 1.0 and the length of the curve should equal 7 to pass our threshold.  Please read the comments in this python script for a better explaination."
                  )
parser.add_option("-T","--expression_threshold",dest="expThreshold",
                  default="",
                  help="This is a more flexible thresholding interface.  We are allowed to manipulate variables to build expressions. Variables allowed are:\n P ,power \n L ,pixel length \n T , time duration in seconds \n F , frequency bandwith of kurve \n  See the help for the candidateList class for better explaination of the valid syntax allowed. ENCLOSE THE EXPRESSION IN DOUBLE QUOTES!",
                  )
parser.add_option("-s","--write_summary",dest="dumpSummaryDisk",
                  default=False,
                  action="store_true",
                  help="This will write the summary information for the candidate file(s) opened. The filename scheme replaces candidate with summary. The data store is Length Pixels,Power,Duration Sec,Bandwidth Hz"
                  )
parser.add_option("-d","--display_summary",dest="dumpSummaryScreen",
                  default=False,
                  action="store_true",
                  help="This will display the summary information for the candidate file(s) opened."
                  )
parser.add_option("-p","--print",dest="print2file",
                  default=False,
                  action="store_true",
                  help="This option states that files specified with the --file option should be printed to a similiar list of coordinate pairs that may be plotted using some external X interface, such as gnuplot. To show you a visual representation of the candidates found.  THIS OPTION NOT YET IMPLEMENTED!"
                  )
(options,args)=parser.parse_args()
filename=str(options.filename)
glob=options.glob
clobberFilename=str(options.clobberFilename)
threshold=str(options.thresholdString)
expThreshold=str(options.expThreshold)
printFile=options.print2file
outfile=str(options.outfile)
dumpSummaryScreen=bool(options.dumpSummaryScreen)
dumpSummaryDisk=bool(options.dumpSummaryDisk)

if filename == "":
    print "Filename argument either not specified or invalid!"
    os.abort()

#Load the file/dir/single file specified by --file option
canList=[]
canList=generateFileList(filename)

#SECTION TO DO THE GLOBBING OF MANY CANDIDATE FILES
if (glob and (canList.__len__() >= 1)):
    newGlobFile=candidateList()
    if outfile != "":
        outName=outfile
    else:
        outName=newGlobFile.__filemaskGlob__()
    #If file preexists erase it first!
    if os.path.isfile(outName):
        print "Prexistinf file found:",outName
        print "Removing!"
        os.unlink(outName)
    for entry in canList:
        newGlobFile.globListFile(outfile,entry)

    #SECTION TO DO THE CLOBBERING OF A CANDIDATE FILE WITH ANOTHER
elif (clobberFilename != '') and (canList.__len__() == 1):
    #Create clobberer
    #Load file(s) to clobber with
    clobberList=generateFileList(clobberFilename)
    #If there is more than one file we need to glob them first
    newCandidateClobberObject=buildCandidateGlob(clobberList)
    #Create candidate list to clobber.
    clobberVictim=candidateList()
    clobberVictim.loadfile(canList[0])
    #Clobber the Victim
    newClobberedList=clobberVictim.clusterClobberWith(newCandidateClobberObject)
    #Write the results to the disk!
    if outfile != "":
        newClobberedList.writefile(outfile)
    else:
        newClobberedList.writefile(clobberVictim.__filemaskClob__(newCandidateClobberObject))

    #SECTION TO APPLY THRESHOLD DO NOT USE THIS ANYMORE
elif ((threshold != "") and (canList.__len__() >= 1)):
    #Apply specified thresholds to --file argument!
    args=str(threshold).split(',')
    for entry in canList:
        candidateObject=candidateList()
        candidateObject.loadfile(entry)
        candidateObject.applyNewThresholds(args[0],args[1],args[2],args[3],args[4])
    pathName=''
    if outfile != "":
        print "Sorry can not just use 1 filename for multiple file saves!"
        print "Taking path information and saving collection of files there."
        pathName=os.path.basename(outfile)
    saveFiles=pathName+'Threshold:'+str(threshold)+':'+candidateObject.filename[0]
    candidateObject.writefile(saveFiles)

    #SECTION 4 PRINTING
elif ((canList.__len__() >= 1) and (printFile)):
    #Iterate of files creating plotable graphic for each!
    for entry in canList:
        candidateObject=candidateList()
        candidateObject.loadfile(entry)
        if (outfile != "") and (canList.__len__() == 1):
            candidateObject.writePixelList(outfile,'tf')
        else:
            pathName=os.path.dirname(candidateObject.filename[0])
            saveFiles=pathName+'/ScatterPlot:'+os.path.basename(candidateObject.filename[0])
            print "Writing scatter plotable data :",saveFiles
            candidateObject.writePixelList(saveFiles,'tf')        
            gnuplotScriptFile(saveFiles)
        del entry
        del candidateObject

    #SECTION APPLY ABITRARY THRESHOLDS
elif ((expThreshold != "") and (canList.__len__() >=1)):
    #Carry out thresholding
    for entry in canList:
        candidateObject=candidateList()
        candidateObject.loadfile(entry)
        candidateResults=candidateObject.applyArbitraryThresholds(expThreshold)
        expThresholdName=str(expThreshold).replace('(','--').replace(')','--')
        pathName=''
        if (outfile != "") and (canList.__len__() == 1):
            candidateResults.writefile(outfile)
        else:
            pathName=os.path.dirname(entry)
            saveFiles=pathName+'/Threshold:'+str(expThresholdName)+':'+os.path.basename(entry)
            print "Writing file :",saveFiles
            candidateResults.writefile(saveFiles)        
        del entry
        del candidateResults
        del candidateObject
    #SECTION TO DUMP SUMMARY TO DISK OR SCREEN
elif ((canList.__len__() >=1) and dumpSummaryDisk):
    for entry in canList:
        candidateObject=candidateList()
        candidateObject.loadfile(entry)
        candidateObject.writeSummary()
        del candidateObject
        
elif ((canList.__len__() >=1) and dumpSummaryScreen):
    for entry in canList:
        candidateObject=candidateList()
        candidateObject.loadfile(entry)
        candidateObject.printSummary()
        del entry

    #IF there are no files to work on found.    
elif (canList.__len__() < 0):
      print "It appears there are no files to process!"

    #THIS SECTION SHOULD NEVER HAPPEN
else:
    print "Error with combination of arguments given!"
    print options
    print args
    os.abort()
