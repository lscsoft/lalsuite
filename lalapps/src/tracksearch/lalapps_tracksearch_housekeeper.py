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
__author__ = 'Cristina Torres <cristina@phys.utb.edu>'
__date__ = '$Date$'
__version__ = ''

import sys
from optparse import OptionParser
import shutil
import ConfigParser
import getopt
import math
import os
import string
import time
import commands

# Begin MAIN part of this simple script
usage="tshousekeeper [args]"
parser = OptionParser()

parser.add_option("-m","--masks",dest="removeMasks",
                  default="",
                  help="This allows a specification of file names to remove X.dat Y.candidates etc from a comma seperate list dat,candidates .  This way we can specify what to clean out.",
                  metavar="MASKS"
                  )
parser.add_option("-a","--action",dest="action",
                  default="no",
                  help="This option allows us to selectively ignore cleaning, archive files associated with MASKS or delete those files all togher.  The valid options here are yes,no, and archive. If you choose archive please note you need a GNU-type tar installation and the tar.gz archive will be placed in PARENT_DIR",
                  metavar="OPTION"
                  )
parser.add_option(
                 "-i","--ignore-path",dest="ignore",
                 default="",
                 help="This allows a simple string to specify paths to not housekeep.  By nature this program will recurse the directory tree clearing out the files specified by MASK.  If you want to skip a type of directory with a command name such as RESULTS.  Specifying --ignore-path=RESULS will mean any file with path /home/user/dataPipe/OTHERTXT_RESULTS/file.MASK will be spared from the cleaning process.  All that is required is a simple check the RESULTS string must be somewhere in the path to the file.MASK",
                 metavar="SPARE_DIR"
                 )
parser.add_option(
                 "-p","--parent-dir",dest="parent",
                 default="",
                 help="This specifies the path (PARENT) where cleaning will start moving down the directory structure.  We will spare directories with RESULTS in their path names.  This ideally should be an absolute path def to avoid ambiguity.",
                 metavar="PARENT_DIR"
                 )
parser.add_option(\
                  "-g","--migrate-to",dest="migrateto",\
                  default=None,
                  help="Specify a complete path when setting this\
option to move the spared files to the migrate location.  If this\
location doesn't exist we will create it.",
                  metavar="MIGRATE_TO"
                  )
parser.add_option(\
                 "-e","--migrate-from",dest="migratefrom",
                 default=None,
                 help="Defaults to option --ignore-path if set\
 otherwise no migration will happen.  Set this to the directory\
 you want surviving files placed into.  Filename collisions are\
 possible we don't check to see what exists before the\
 migration.",
                 metavar="MIGRATE_FROM"
                 )    
(options,args)=parser.parse_args()
masks=str(options.removeMasks).split(',');
if (masks == ""):
    print "You didn't specify a file extension mask to remove!"
    os.abort()
else:
    print "Processing the following mask(s)."
    print masks

action=str(options.action).lower().lstrip().rstrip()
actionKeywords=['yes','no','archive']
if not(actionKeywords.__contains__(action)):
    print "Error with action specified! ",action," keyword not defined."
    os.abort()
if action == 'archive':
    print "Sorry the tar feature is not implemented yet!"
    print "Rerun with yes/no options."
    os.abort()

parentDir=os.path.abspath(str(options.parent))
if (parentDir == '/'):
    print "You requested cleaning starting and the file system root / ARE YOU SURE?"
    print "This is unsafe. Sorry I'm aborting."
    os.abort()

ignorePathStr=str(options.ignore)
if (ignorePathStr == ""):
    print "You are choosing to process all subdirectories of ",parentDir
    print "No files matching MASKS will be spared!"
else:
    print "Ignoring paths containing the string :",ignorePathStr

migrateTo=options.migrateto
migrateFrom=options.migratefrom
if ((migrateFrom == None) and (ignorePathStr != "")):
    migrateFrom=ignorePathStr

#Create the structure containing all the files down from the parent directory.
files2clean=[]
allFilesSeen=[]
for root,dir,files in os.walk(parentDir):
    for entry in files:
        thisFile=os.path.join(root,entry)
        allFilesSeen.append(thisFile)
        for maskKey in masks:
            if os.path.basename(thisFile).__contains__(maskKey):
                    files2clean.append(thisFile)

#Keep only files which are not in the ignorePathStr
files2remove=[]
if (ignorePathStr != ""):
    print "Checking for files to keep."
    for entry in files2clean:
        if not(entry.__contains__(ignorePathStr)):
            files2remove.append(entry)
else:
   files2remove=files2clean
    
# Report number of files 2 clean
print "We found ",files2remove.__len__()," files in ",parentDir," and its subdirectories to clean."
print "We spared ",files2clean.__len__() - files2remove.__len__()," files from cleaning."
fileSizeSum=sum([os.path.getsize(entry) for entry in files2remove])
fileSizeSpared=sum([os.path.getsize(entry) for entry in files2clean])-fileSizeSum
print "Matching files to remove consume ",fileSizeSum," bytes of disk space."
print "Clearing additional spared file(s) will restore ",fileSizeSpared," more bytes of disk space freed."

#Remove the files we are supposed to remove!
if (action == 'yes'):
    print "Removing files listed below: "
    for file in files2remove:
        print file
        os.remove(file)
    print "Done."
else:
    print "Leaving file system intact!"
 
#Migrate files not cleaned out of the directory structure
if (migrateTo != None) and (migrateFrom != None):
    #Set paths
    migrateTo=os.path.abspath(migrateTo)
    if ignorePathStr != migrateFrom:
        migrateFrom=os.path.abspath(migrateFrom)
    #Create migrateTo location if not there!
    if not os.path.isdir(migrateTo):
        os.makedirs(migrateTo)
    #Loop over the files to be migrated.
    if os.path.isdir(migrateTo):
        for file in allFilesSeen:
            #If file lives in migrateFrom path
            if file.__contains__(migrateFrom):
                print "Moving file ",file," to ",migrateTo
                shutil.move(file,migrateTo)

