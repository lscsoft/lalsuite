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
__author__ = 'Cristina Torres <cristina@phys.utb.edu>'
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
disableGraphics=False
try:
    from pylab import *
except RuntimeError,ImportError:
    disableGraphics=True

# Need to add import for pylab and if the import fails disable the
# printing options available from the command line with error messages
# to the user.

#Begin MAIN part of this python code to actually manipulate the curve
#lists


#################### Internal Methods #################

def buildCandidateGlob(fileList,verboseSwitch=False):
    #Build glob var from filelist.  Memory friendly
    #Returns copy of results
    canList=fileList
    if canList.__len__() == 1:
        if verboseSwitch:
            print "Globbing is not possible only one file found!"
            print "Returning the result as that entry!"
        tmpObject=candidateList()
        #tmpObject.loadfile(canList.pop(0))
        tmpObject.__loadfileQuick__(canList.pop(0))
        newCandidateObject=copy.deepcopy(tmpObject)
    else:
        newCandidateObject=candidateList()
        #newCandidateObject.loadfile(canList.pop(0))
        newCandidateObject.__loadfileQuick__(canList.pop(0))
    for entry in canList:
        if verboseSwitch:
            print " "
            print "Loading candidate list file: ",entry
        tmpCandidate=candidateList()
        tmpCandidate.__loadfileQuick__(entry)
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
    txtTemplate='plot "%s" with lines\n set size 1.0,0.6\n set terminal postscript portrait enhanced mono dashed lw 1 "Helvetica" 14\n set output "%s.ps"\n replot\n set terminal x11\n set size 1,1\n set title "%s"\n'
    output_fp=open(filename+'.plt','w')
    output_fp.write(txtTemplate%(filename,filename,os.path.basename(filename)))
    output_fp.close()
#End method gnuplotScriptFile()

###################End internal methods #################
usage=''
if disableGraphics:
    usage = "X windows not detected graphing flags not available.\ncandidateHandler.py [args]"
else:
    usage = "candidateHandler.py [args]"
parser = OptionParser(usage)
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
parser.add_option("-T","--expression_threshold",dest="expThreshold",
                  default="",
                  help="This is a flexible thresholding interface.  We are allowed to manipulate variables to build expressions. Variables allowed are:\n P ,power \n L ,curve length \n D , time duration in seconds \n B , frequency bandwith of kurve \n  T, start time float\n S, stop time float\n F, start Freq \n, G, stop Freq \n V, bright pixel frequency\n H, bright pixel time(float)\n J, power of bright pixel \n M, mean pixel power \n C, variance pixel power \n See the help for the candidateList class for better explaination of the valid syntax allowed. You can impose multiple threshold expresssions for one file operation.  This saves time lost to IO of running candidateHandler.py multiple times.  Use the comma as a delimiter B>10,L<100 would apply to thresholds and save them in two appropriately marked files. ENCLOSE THE EXPRESSION IN DOUBLE QUOTES",
                  )
parser.add_option("-s","--write_summary",dest="dumpSummaryDisk",
                  default=False,
                  action="store_true",
                  help="This will write the summary information for the candidate file(s) opened. The filename scheme replaces candidate with summary. The data store is Length Pixels,Power,Duration Sec,Bandwidth Hz.  This option can be invoked along with the thresholding option to write out summary files of the results immediately after thresholding the data."
                  )
parser.add_option("-d","--display_summary",dest="dumpSummaryScreen",
                  default=False,
                  action="store_true",
                  help="This will display the summary information for the candidate file(s) opened."
                  )
parser.add_option("-p","--print",dest="imageFilemask",
                  default='',
                  help="This option states that files specified with the --file option should be printed to a similiar list of coordinate pairs that may be plotted using some external X interface, such as gnuplot. To show you a visual representation of the candidates found.  THIS OPTION NOT YET IMPLEMENTED! IT IS BEING REDONE DO NOT USE THIS OPTION! Should eventually take a filename as output arg!b"
                  )
parser.add_option("-i","--image",dest="imageCount",
                  default='0',
                  help="This switch allows you to tell the candidateHandler to graph the data, if the input library is clobbered together from many candidate lists then specify the maximum number of figures you are willing to create.  The default is currently 1")
parser.add_option("-x","--stats",dest="stats",
                  default='',
                  help="Determine the specified stats for all candidates in input file for variable in question.  See the -T options variable labels for explaination.  This result is written out to the screen as mean,min,max,std.  To measure stats of more than one trait use a comma seperated list P,L as input args. FUNCTIONALITY NOT IMPLEMENTED YET!"
                  )
parser.add_option("-v","--verbose",dest="verbose",
                  default=False,
                  action="store_true",
                  help="Sets up the verbose features. Prints diagnostic messages and progress meters where appropriate.")
(options,args)=parser.parse_args()
filename=str(options.filename)
glob=options.glob
clobberFilename=str(options.clobberFilename)
expThreshold=str(options.expThreshold)
graph2file=printFile=options.imageFilemask
graph2screen=int(options.imageCount)
outfile=str(options.outfile)
dumpSummaryScreen=bool(options.dumpSummaryScreen)
dumpSummaryDisk=bool(options.dumpSummaryDisk)
verboseSwitch=bool(options.verbose)
measureStats=str(options.stats)
if verboseSwitch:
    print "Setting verbose mode to candidateHandler call."
    
if filename == "":
    print "Filename argument either not specified or invalid!"
    os.abort()

#Load the file/dir/single file specified by --file option
canList=[]
canList=generateFileList(filename)

#SECTION TO DO THE GLOBBING OF MANY CANDIDATE FILES
if (glob and (canList.__len__() >= 1)):
    if outfile != "":
        outName=outfile
    else:
        outName=newGlobFile.__filemaskGlob__()
    #If file preexists erase it first!
    if os.path.isfile(outName):
        if verboseSwitch:
            print "Prexisting file found:",outName
            print "Removing!"
        os.unlink(outName)
    #Create new globbed candidate structure to write to disk!
    newGlobFile=buildCandidateGlob(canList,verboseSwitch)
    #SECTION TO DO THE CLOBBERING OF A CANDIDATE FILE WITH ANOTHER
elif (clobberFilename != '') and (canList.__len__() == 1):
    #Create clobberer
    #Load file(s) to clobber with
    clobberList=generateFileList(clobberFilename)
    #If there is more than one file we need to glob them first
    newCandidateClobberObject=buildCandidateGlob(clobberList,verboseSwitch)
    #Create candidate list to clobber.
    clobberVictim=candidateList(verboseSwitch)
    clobberVictim.__loadfileQuick__(canList[0])
    #Clobber the Victim
    newClobberedList=clobberVictim.clusterClobberWith(newCandidateClobberObject)
    #Write the results to the disk!
    if outfile != "":
        newClobberedList.writefile(outfile)
    else:
        newClobberedList.writefile(clobberVictim.__filemaskClob__(newCandidateClobberObject))

    #SECTION 4 PRINTING
elif ((canList.__len__() >= 1) and (printFile)):
    #Iterate of files creating plotable graphic for each!
    if verboseSwitch:
        print "Preparing for printing :",canList.__len__()," files."
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        candidateObject.__loadfileQuick__(entry)
        if (outfile != "") and (canList.__len__() == 1):
            candidateObject.writePixelList(outfile,'tf+time')
        else:
            pathName=os.path.dirname(candidateObject.filename[0])
            saveFiles=pathName+'/ScatterPlot:'+os.path.basename(candidateObject.filename[0])
            candidateObject.writePixelList(saveFiles,'tf+time')        
            gnuplotScriptFile(saveFiles)
        del entry
        del candidateObject
    #Create a scatterplot on screen of the curves in the candidateList
elif ((canList.__len__() >=1) and ((graph2screen>0) or (graph2file!=''))):
    if (graph2screen < canList.__len__()):
        print "Number of images to create exceeds your request!"
        print "Estimate images in this data structure are ",canList.__len__()
        #Needs continue anyway option here!
        os.abort()
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        candidateObject.__loadfileQuick__(entry)
        candidateObject.graphdata(graph2file)
        del candidateObject
    #SECTION APPLY ABITRARY THRESHOLDS
elif ((expThreshold != "") and (canList.__len__() >=1)):
    #Carry out thresholding on all entries in canList
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        candidateObject.__loadfileQuick__(entry)
        #Build Threshold list
        expThresholdLIST=str(expThreshold).split(',')
        for singleThreshold in expThresholdLIST:
            candidateResults=candidateObject.applyArbitraryThresholds(singleThreshold)
            singleThresholdName=str(singleThreshold).replace('(','--').replace(')','--')
            pathName=''
            if (outfile != "") and (canList.__len__() == 1):
                #Prepend Threshold label to outfile
                if expThresholdLIST.__len__()>1:
                    outfile=singleThreshold+outfile
                candidateResults.writefile(outfile)
                if dumpSummaryDisk:
                    candidateResults.writeSummary(outfile)
            else:
                pathName=os.path.dirname(entry)
                saveFiles=pathName+'/Threshold:'+str(singleThresholdName)+':'+os.path.basename(entry)
                candidateResults.writefile(saveFiles)
                if verboseSwitch:
                    print "Wrote file :",saveFiles
            if dumpSummaryDisk:
                candidateResults.writeSummary(saveFiles)
            del candidateResults
        del entry
    del candidateObject

    #SECTION TO DUMP SUMMARY TO DISK WHEN NO THRESHOLDING REQD
elif ((canList.__len__() >=1) and dumpSummaryDisk and (expThreshold =="") ):
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        if verboseSwitch:
            print "Processing:",entry
        candidateObject.__loadfileQuick__(entry)
        candidateObject.writeSummary()
        del candidateObject
    #SECTION TO DUMP SUMMARY TO SCREEN        
elif ((canList.__len__() >=1) and dumpSummaryScreen):
    for entry in canList:
        candidateObject=candidateList(verboseSwitch)
        candidateObject.__loadfileQuick__(entry)
        candidateObject.printSummary()
        del entry

    #IF there are no files to work on found.    
elif (canList.__len__() < 0):
      print "It appears there are no files to process!"

    #THIS SECTION SHOULD NEVER HAPPEN
else:
    print "Error with combination of arguments given!"
    print "Legitimate command line arguments don't get this error!"
    print "Options in parser data structure."
    print options
    print "Corresponding argument values."
    print args
    print "Candidates files found :",canList.__len__()
    os.abort()
 
