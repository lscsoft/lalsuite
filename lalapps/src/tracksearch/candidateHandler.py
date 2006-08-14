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

parser = OptionParser()

parser.add_option("-f","--file",dest="filename",
                  default="",
                  help="This specifies the list of curve files to process contained in FILE, either rethresholding, globbing.  If the input is a directory then it is assumed that we will be working on all the seen files.  /PATH/*.ext is allowed.",
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
parser.add_option("-p","--print",dest="print2file",
                  default=False,
                  action="store_true",
                  help="This option states that files specified with the --file option should be printed to a similiar list of files that may be plotted using some external X interface, such as gnuplot. To show you a visual representation of the candidates found.  THIS OPTION NOT YET IMPLEMENTED!"
                  )
(options,args)=parser.parse_args()
filename=str(options.filename)
glob=options.glob
clobberFilename=str(options.clobberFilename)
threshold=str(options.thresholdString)
printFile=options.print2file
outfile=str(options.outfile)

#Load the file/dir/single file specified by --file option
canList=[]
canList=generateFileList(filename)

#What are we doing with it?
if (glob and (canList.__len__() > 1)):
    canObjects=[]
    for entry in canList:
        tmpCandidate=candidateList()
        tmpCandidate.loadfile(entry)
        canObjects.append(tmpCandidate)
    #Glob the input files and exit writing to disk!
    newCandidateObject=copy.deepcopy(canObjects.pop(0))
    for entry in canObjects:
        newCandidateObject=newCandidateObject.globList(entry)
    if outfile != "":
        newCandidateObject.writefile(outfile)
    else:
        newCandidateObject.writefile(newCandidateObject.__filemaskGlob__())
elif (clobberFilename != '') and (canList.__len__() == 1):
    #Create clobberer
    #Load file(s) to clobber with
    clobberList=generateFileList(clobberFilename)
    #If there is more than one file we need to glob them first
    tmpObjClobber=candidateList()
    for entry in clobberList:
        object=candidateList()
        object.loadfile(entry)
        tmpObjClobber.append(object)
    #Do the glob
    newCandidateClobberObject=candidateList()
    for entry in tmpObjClobber:
        newCandidateClobberObject=newCandidateClobberObject.globList(entry)
    #Create candidate list to clobber.
    clobberVictim=candidateList()
    clobberVictim.loadfile(canList[0])
    #Clobber the Victim
    newClobberedList=clobberVictim.clusterCheckAgainst(newCandidateClobberObject)
    #Write the results to the disk!
    if outfile != "":
        newClobberedList.writefile(outfile)
    else:
        newClobberedList.writefile(clobberVictim.__filemaskClob__(newCandidateClobberObject))
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
elif ((canList.__len__() >= 1) and (printFile)):
    #Iterate of files creating plotable graphic for each!
    print "Printing not ready yet!!"
else:
    print "Error with combination of arguments given!"
    os.abort()

