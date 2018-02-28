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
"""
This script is responsible for creating a pickle background file for
calculating the followup background percentages that will appear in
the checklists.  Pre-built background files rely on epoch definitions
from fu_utils and significantly decrease the execution time of veto
and dq queries from the followup codes.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__    = 'followupGenerateDQBackground.py'

import optparse
import sys
import os
from pylal import git_version
from pylal.fu_utils import followupDQV,\
     getKnownIfos,\
     getKnownIfoEpochs,\
     getRunTimesAsLIGOSegment

sys.path.append('@PYTHONLIBDIR@')

usage = """usage: %prog [options]"""
parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
#Add all options to setup the query
parser.add_option("-i","--ifo-epoch-list",action="store",type="string",\
                  metavar="X1:EPOCH",default="",\
                  help="Specify the IFOs and respective epochs to \
build the DQ background for.  The proper format for this argument \
is as follows: L1:S6A,H1:S6A,V1:VSR2.  You can make incoherent \
combinations, so be aware of the arguments given."\
                  )
parser.add_option("-o","--output-location",action="store",type="string",\
                  metavar="abspath",default=None,\
                  help="Specify the path to the place the background \
DQ pickle into when complete.  This should probably be the same \
location as the Qscan background directory tree but in \DQ instead \
of \Qscan.  The actual filename is controlled via a fixed filemask \
this filemask composes a filename based of ifos and epoch \
combinations chosen with --ifo-epoch-list."\
                  )
parser.add_option("-s","--show-ifo-epoch",action="store_true",\
                  default=False,\
                  help="Specify this flag to have the script \
dump out all know IFO names and named search epochs.  These are \
defined in fu_utils module."\
                  )
parser.add_option("-a","--available-backgrounds",action="store_true",\
                  default=False,\
                  help="Specify this flag to query the default \
location of background pickles to see which are available on \
this system.")
parser.add_option("-v","--verbose",action="store_true",\
                  default=False,\
                  help="Turns on verbose status updates to stdout."\
                  )
######################################################################
(opts,args) = parser.parse_args()
#Create a DQ querying object
dqService=followupDQV()
if opts.show_ifo_epoch:
    for myIfo in getKnownIfos():
        validEpochs=getKnownIfoEpochs(myIfo)
        sys.stdout.write("%s has valid epochs %s\n"%(myIfo,validEpochs))
    sys.exit(0)
if opts.available_backgrounds:
    sys.stdout.write("These are all possible combinations of ifo and \
epoch. Be Warned! Some combinations may be nonsense.\n")
    #Determine all background combinations
    allIfos=getKnownIfos()
    #Ifo singles
    sys.stdout.write("Singles\n")
    for j in allIfos:
        for myEpoch in getKnownIfoEpochs(j):
            mySegment=getRunTimesAsLIGOSegment(myEpoch,j)
            myPickle=dqService.figure_out_pickle([(j,myEpoch)])
            #Check for existing pickle
            if os.path.isfile(myPickle) and \
               not os.access(myPickle,os.R_OK):
                sys.stderr.write(
                    "%s:%s found without \
 read permissions!\n"%(j,myEpoch))
            elif not os.path.isfile(myPickle):
                sys.stdout.write(
                    "%s:%s background file \
 is not found!\n"%(j,myEpoch))
            else:
                sys.stdout.write(
                    "%s:%s found at %s\n"\
                    %(j,myEpoch,myPickle))
    #Ifo doubles
    sys.stdout.write("Doubles\n")
    doubleList=list()
    for j in allIfos:
        for k in allIfos:
            pair=[j,k]
            pair.sort()
            (a,b)=pair
            if ((a,b) not in doubleList) and (a!=b):
                doubleList.append((a,b))
    for a,b in doubleList:
        epochsA=getKnownIfoEpochs(a)
        epochsB=getKnownIfoEpochs(b)
        for epochA in epochsA:
            for epochB in epochsB:
                segmentA=getRunTimesAsLIGOSegment(epochA,a)
                segmentB=getRunTimesAsLIGOSegment(epochB,b)
                if segmentA.intersects(segmentB):
                    epochCombo=[(a,epochA),(b,epochB)]
                    expectedPickle=dqService.figure_out_pickle(epochCombo)
                    #Check to see if pickle is there
                    if os.path.isfile(expectedPickle) and \
                       not os.access(expectedPickle,os.R_OK):
                        sys.stderr.write(
                            "%s:%s,%s:%s found without \
read permissions!\n"%(a,epochA,b,epochB))
                    elif not os.path.isfile(expectedPickle):
                        sys.stdout.write(
                            "%s:%s,%s:%s background file \
is not found!\n"%(a,epochA,b,epochB))
                    else:
                        sys.stdout.write(
                            "%s:%s,%s:%s found at %s\n"\
                            %(a,epochA,b,epochB,expectedPickle))
    #Ifo triples
    sys.stdout.write("Triples.\n")
    tripleList=list()
    for j in allIfos:
        for k in allIfos:
            for z in allIfos:
                pair=[j,k,z]
                pair.sort()
                (a,b,c)=pair
                if ((a,b,c) not in tripleList) and (a!=b) and (b!=c) and (a!=c):
                    tripleList.append((a,b,c))
    for a,b,c in tripleList:
        epochsA=getKnownIfoEpochs(a)
        epochsB=getKnownIfoEpochs(b)
        epochsC=getKnownIfoEpochs(c)
        for epochA in epochsA:
            for epochB in epochsB:
                for epochC in epochsC:
                    segmentA=getRunTimesAsLIGOSegment(epochA,a)
                    segmentB=getRunTimesAsLIGOSegment(epochB,b)
                    segmentC=getRunTimesAsLIGOSegment(epochC,c)
                    if segmentA.intersects(segmentB) and \
                       segmentB.intersects(segmentC) and \
                       segmentC.intersects(segmentA):
                        epochCombo=[(a,epochA),(b,epochB),(c,epochC)]
                        expectedPickle=dqService.figure_out_pickle(epochCombo)
                        #Check to see if pickle is there
                        if os.path.isfile(expectedPickle) and \
                           not os.access(expectedPickle,os.R_OK):
                            sys.stderr.write(
                                "%s:%s,%s:%s,%s:%s found without \
read permissions!\n"%(a,epochA,b,epochB,c,epochC))
                        elif not os.path.isfile(expectedPickle):
                            sys.stdout.write(
                                "%s:%s,%s:%s,%s:%s background file \
is not found!\n"%(a,epochA,b,epochB,c,epochC))
                        else:
                            sys.stdout.write(
                                "%s:%s,%s:%s,%s:%s found at %s\n"\
                                %(a,epochA,b,epochB,c,epochC,expectedPickle))
    sys.exit(0)
#
#Check the form of the ifo epoch option
#
epochInfo=opts.ifo_epoch_list.upper().strip()
ifoEpochList=list()
for myPair in epochInfo.split(","):
    myIfo,myEpoch=myPair.split(":")
    myIfo=myIfo.upper()
    myEpoch=myEpoch.lower()
    if opts.verbose:
        sys.stdout.write("Checking %s %s\n"%(myIfo,myEpoch))
    if myIfo not in getKnownIfos():
        sys.stderr.write(" %s not valid ifo.\n"%myIfo)
        sys.stderr.write("Valid ifos are %s\n"%getKnownIfos())
        sys.exit(-1)
    if myEpoch not in getKnownIfoEpochs(myIfo):
        sys.stderr.write("%s not valid epoch for %s ifo.\n"%(myEpoch,myIfo))
        sys.stderr.write("valid epochs for %s ifo are %s\n."%\
                         (myIfo,getKnownIfos()))
        sys.exit(-1)
    ifoEpochList.append((myIfo,myEpoch))
#Check for proper permissions etc to the output location etc
outputFullPath=dqService.figure_out_pickle(ifoEpochList)
(systemPath,filename)=os.path.split(outputFullPath)
if opts.output_location:
    systemPath=os.path.expanduser(os.path.normpath(opts.output_location))
if opts.verbose:
    sys.stdout.write("Checking Path for write permission \
%s\n"%systemPath)
    if not os.path.isdir(systemPath):
        sys.stderr.write("Output path does not exist.\n")
        try:
            os.makedirs(systemPath)
        except:
            sys.stderr.write("Error creating path!\n")
            sys.exit(-2)
    if not os.access(systemPath,os.W_OK):
        sys.stderr.write("User does not have write permission to \
 %s\n"%systemPath)
outputFullPath=os.path.normpath(systemPath+"/"+filename)
if opts.verbose:
    sys.stdout.write("Storing background to %s\n"%outputFullPath)
#Begin the string of queries
if opts.verbose:
    sys.stdout.write("Beginning query of DQ server.\n")
    ifoCount=len(ifoEpochList)
    queryCount=dqService.__backgroundPoints__
    sys.stdout.write("A total of %i queries will be made to \
server.\n"%(ifoCount*queryCount))
dqService.createDQbackground(ifoEpochList,outputFullPath)
if opts.verbose:
    sys.stdout.write("Done conducting queries.")
sys.exit(0)

    
                     
                     
        
