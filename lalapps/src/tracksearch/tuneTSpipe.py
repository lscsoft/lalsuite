#!/usr/bin/env python2.3
"""
This is a complete standalone python script which will acess the
search ini file and call the tracksearch_pipe.py program to setup all
the pipes needed to tune the search algorithm.  This script will also
require its own ini file to setup the analysis tuning.  The convention
we want to follow for tuning ini files to give them the extens .tun
"""

from candidateUtils import *
from glue import pipeline
from optparse import OptionParser
from tracksearch import buildDir
import ConfigParser
import commands
import copy
import getopt
import math
import os
import string
import sys
import time
import tracksearch

#
#
# The following are all class/methods defined to actually carry out components
# of the tuning process
#
#
class tuneObject:
    def __init__(self,tunFileCP):
        self.cpTun=tunFileCP
        self.masterIni=cp.get('all','masterini')
        self.cpIni=ConfigParser.ConfigParser()
        self.cpIni.read(self.masterIni)
        self.lambaList=[]
        self.lamHopts=cp.get('all','LH')
        self.lamLopts=cp.get('all','LL')
        if self.lamHopts.count(";")!=2:
            print "Error with LH ini file delimiters!"
            os.abort()
        if self.lamLopts.count(";")!=2:
            print "Error with LL ini file delimiters!"
            os.abort()
        self.LH=self.lamHopts.split(";")
        self.LL=self.lamLopts.split(";")
        self.myIni=cp.get('all','masterini')
        self.batchMask=cp.get('all','iniBatchLabel')
        self.home=cp.get('all','tuningHome')
        self.installPipes=self.home+'/pipes'
        self.log=cp.get('all','tuningLogs')
        self.seglist=cp.get('all','seglist')
        self.FApipeNames=[]
        self.DEpipeNames=[]
        self.pipeBuilder=cp.get('all','pipeProgram')
        #Create directory to place all pipes if not already present!
        buildDir(self.home)
    #End __init__ method 

    def performFAsetup(self):
        h=float(self.LH[2])
        l=float(self.LL[2])
        deltah=float(self.LH[1])
        deltal=float(self.LL[1])
        stoph=float(self.LH[0])
        stopl=float(self.LL[0])
        if h>stoph:
            print "Error in config, inconsistent LH options."
            os.abort()
        if l>stopl:
            print "Error in config, inconsistent LL options."
            os.abort()
        i=0
        while (h <= stoph):
            l=float(self.LL[2])
            while (l <= stopl):
                coord=[h,h*l]
                pH=str(str(coord[0]).__getslice__(0,int(str(coord[0]).index(".")+4))).zfill(8)
                pL=str(str(coord[1]).__getslice__(0,int(str(coord[1]).index(".")+4))).zfill(8)
                pipeIniName=self.home+'/'+'FAsetup_'+self.batchMask+':'+pH+':'+pL+':'+'.ini'
                self.FApipeNames.append([pipeIniName,pH,pL,coord[0],coord[1]])
                l=l+deltal
                i=i+1
            h=h+deltah
        print "Creating each pipe configuration script."
        for pipe in self.FApipeNames:
            self.FAconfigure(pipe[3],pipe[4],pipe[0])
        print "Building ",self.FApipeNames.__len__()," pipes. This may take awhile."
        for pipe in self.FApipeNames:
            self.createPipe(pipe[0])
        print "Ok. Finished"
    #end performFAsetup method

    def FAconfigure(self,LH,LL,iniName):
        """
        Takes the cp and manipulates it to make a config file object.
        This object is suited for the FA rate investigation and written
        to disk with specified filename
        """
        newCP=copy.deepcopy(self.cpIni)
        # New/Rewritten ini options
        uniqStr='_'+str(LH)+'_'+str(LL)+'_'
        newCP.set('filelayout','workpath',self.installPipes+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        # If injection section present strip it out!
        if newCP.has_section('tracksearchinjection'):
            newCP.remove_section('tracksearchinjection')
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End FAconfigure method

    def createPipe(self,pipeIniName):
        """
        This method actually invokes a shell to create the pipe in a unique
        directory under the main directory specified in [all] tuningHome
        which tailors each pipe ini file entry [filelayout] workpath
        """
        seglist=self.seglist
        path2bin=self.pipeBuilder
        commandLine=path2bin+' --file='+pipeIniName+' --segment_file='+seglist
        commandOut=commands.getstatusoutput(commandLine)
        if commandOut[0] != 0:
            print " "
            print "There was a problem installing a pipe for ini file :",pipeIniName
            print " "
            print commandOut[1]
            print " "
    # End method createPipe

    def performDEsetup(self):
        #1) Read in lambaH lambaL P L result file
        #2) Create corresponding ini file
        #3) Run pipe creation script to create injection pipes
        myIni=self.myIni
        results=self.readFAresults()
        for entry in results:
            h=entry[0]
            l=entry[1]
            p=entry[2]
            l=entry[3]
            pipeIniName='DEsetup_'+self.home+self.batchMask+str(h)+str(float(h*l))+'.ini'
            self.DEpipeNames.append(pipeIniName)
            self.DEeditIniFile(pipeIniName)
            self.createPipe(pipeIniName)
    #End performDEsetup

    def DEeditIniFile(self,LH,LL,P,L,iniName):
        """
        This method will edit the original ini file in such a way as to
        prepare the detection efficiency runs to be launched.
        """
        newCP=copy.deepcopy(self.cp)
        #New/Rewritten ini options
        newCP.set('filelayout','workpath','NEW_WORKPATH')
        newCP.set('logpath','logpath','NEW_LOGPATH')
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        newCP.set('tracksearchbase','length_threshold',L)
        newCP.set('tracksearchbase','power_threshold',P)
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End DEeditIniFile

    def performFAcalc(self):
        #create 2 list of longest p,l each pipe
        #calculate the specified p,l for tun file percentile
        #write a rank file to disk 4c form

        print 'hi'
    #End FAcalc method

    def performDEcalc(self):
        #count num of can files with 1 entry+
        #determine percentage succsessful
        #write 2 ranked 4c list lh ll #map percentage

        print 'there'
    #End performDEcalc
#
#
# This is where the main body of the python script will reside.
#
#
parser = OptionParser()

parser.add_option("-f","--file",dest="tunFile",
                  default="",
                  help="This specifies the complete path to the initialization file which will be used to carry out the tracksearch pipeline tuning process.",
                  metavar="FILE.tun"
                  )
parser.add_option("-a","--false-alarm-setup",dest="FAsetup",
                  default=False,
                  action="store_true",
                  help="Specify this option to create all pipeline jobs required to determine the optimal search parameters for all possible \lambda_{h,l} pairs")
parser.add_option("-d","--detect-efficiency-setup",dest="DEsetup",
                  default=False,
                  action="store_true",
                  help="Specify this option after determining our various false alarm rates.  It is in this step that we perform software injections to determine the efficiency of the parameter options.  By this we mean what is the recovery rate of the injected signals.")
parser.add_option("-c","--false-alarm-calculate",dest="FAcalc",
                  default=False,
                  action="store_true",
                  help="Specify this option after executing the false-alarm-setup(FAS) jobs successfully.  The script will then parse the results from the FAS jobs and determine the values of integrated curve power (IP) and curve length (CL) to achieve the false alarm rates requested in the tun file.  The solution will be written to the parent directory which contains all the pipeline tuning jobs.")
parser.add_option("-e","--detect-effeciency-calculate",dest="DEcalc",
                  default=False,
                  action="store_true",
                  help="Specify this option after executing the detect-efficiency-setup (DES) jobs successfully.  The script will then parse through the results from the DES jobs and rank the parameter space point investigate by detection efficiency createing two files a simple rank list and secondary file which can be plotted into efficiency contours with a program such as Matlab.")

(options,args)=parser.parse_args()

#Check the selected options make sure only one setup/calculate flag is called
tunFile=os.path.abspath(str(options.tunFile))
FAsetup=bool(options.FAsetup)
DEsetup=bool(options.DEsetup)
FAcalc=bool(options.FAcalc)
DEcalc=bool(options.DEcalc)
if options.tunFile == "":
    print "The tun file required to perform any step of the tuning process is not specified!"
    os.abort()
elif not os.path.exists(tunFile):
    print 'Can not find iniFile: ',os.path.basename(tunFile),' in path ',os.path.dirname(tunFile)
    os.abort()
countOpts=0

if FAsetup:
    countOpts=countOpts+1
if DEsetup:
    countOpts=countOpts+1
if FAcalc:
    countOpts=countOpts+1
if DEcalc:
    countOpts=countOpts+1
if countOpts != 1:
    print "There is a problem with the options specified at the command line.  Only one action at a time can be specified!"
    os.abort()
#
# Setup the configuration file (tun) variable which points to the needed vals
#
cp = ConfigParser.ConfigParser()
cp.read(tunFile)
#
# The main conditional branch to the individual methods
#
tuner=tuneObject(cp)
if FAsetup:
    tuner.performFAsetup()
elif DEsetup:
    tuner.performDEsetup()
elif FAcalc:
    tuner.performFAcalc()
elif DEcalc:
    tuner.performDEcalc()
else:
    print "This should never happen!"
    os.abort()
#
# End main program
#
