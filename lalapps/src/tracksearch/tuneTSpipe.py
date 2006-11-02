#!/usr/bin/env python
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
import pickle

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
        self.installPipes=self.home+'/FA_pipes'
        self.installIni=self.home+'/FA_ini'
        self.installPipes2=self.home+'/DE_pipes'
        self.installIni2=self.home+'/DE_ini'
        self.log=cp.get('all','tuningLogs')
        self.seglist=cp.get('all','seglist')
        self.mySigmaFA=float(cp.get('false-alarm-calculate','FAR'))
        self.FApipeNames=[]
        self.DEpipeNames=[]
        self.pipeBuilder=cp.get('all','pipeProgram')
        #Create directory to place all pipes if not already present!
        buildDir(self.home)
    #End __init__ method 

    def __createLHLL__(self):
        """
        This method makes a simple 5 element listing for later use.
        The data will
        look like [pipeName,prettyH,prettyL,H,L].  It will span the options in
        [all],LH and LL
        """
        i=0
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
        while (h <= stoph):
            l=float(self.LL[2])
            while (l <= stopl):
                coord=[h,h*l]
                pH=str(str(coord[0]).__getslice__(0,int(str(coord[0]).index(".")+4))).zfill(8)
                pL=str(str(coord[1]).__getslice__(0,int(str(coord[1]).index(".")+4))).zfill(8)
                pipeIniName=self.installIni+'/'+self.batchMask+':'+pH+':'+pL+':'+'.ini'
                self.FApipeNames.append([pipeIniName,pH,pL,coord[0],coord[1]])
                l=l+deltal
                i=i+1
            h=h+deltah
    #End __createLHLL__

    def __findFiles__(self,searchPath,keywordList):
        """
        Method that will involk walk to scan a path hierarchy and check that each
        abs path matches the entries in the list keywordList ['a','b','c'] in filename.
        All keyword must match the filename somewhere to be returned as a result.
        """
        filesinPath=[]
        for root,dir,files in os.walk(searchPath):
            for entry in files:
                thisFile=os.path.join(root,entry)
                hits=0
                for key in keywordList:
                    if thisFile.__contains__(key):
                        hits=hits+1
                if hits >= keywordList.__len__():
                    filesinPath.append(thisFile)
        return filesinPath
    #End method __findFiles__

    
    def __determineDEinPath__(self,candPath):
        """
        This uses the shell and grep to count the number of candidate files with
        a trigger in some particular path
        """
        cmdString='grep -h -i "Total Curves" '+str(candPath)+'/* | grep -c -v 0'
        cmdString2='grep -h -i "Total Curves" '+str(candPath)+'/* | wc -l '
        errs=commands.getstatusoutput(cmdString)
        errs2=commands.getstatusoutput(cmdString2)
        filecount=int(errs2[1])
        trigFound=int(errs[1])
        return [trigFound,filecount]
    # End __determineDEinPath__

    def __buildDAGfrom__(self,dagFilename,root2dags):
        """
        Scans all subdirectories to root2dag for dags to create a master dag file for tuning execution.
        This module uses the shell and find utilities to locate all dags for inclusion.
        """
#        cmdString='find '+root2dags+'/* | grep ".dag"'
#        errs=commands.getstatusoutput(cmdString)
 #       if errs[0] != 0:
  #          print "There was an error finding the individual dag submit files for each pipe!"
  #          print errs[1]
  #          print " "
  #      fileList=errs[1].split('\n')
        fileList=self.__findFiles__(root2dags,['.dag'])
        print "The master DAG will have ",fileList.__len__()," nodes/subdags."
        print "Preparing the master dag may take awhile."
        for dagEntry in fileList:
            cmdString2='condor_submit_dag -no_submit '+dagEntry
            errs=commands.getstatusoutput(cmdString2)
            if errs[0] != 0:
                print "There is a problem creating Master DAG node :",dagEntry
                print errs[1]
                print " "
                print " "
        i=0
        masterDAG_fp=file(dagFilename,'w')
        masterDAG_fp.write('DOT MasterDag.dot\n')
        for entry in fileList:
            masterDAG_fp.write('Job '+str(i)+' '+str(entry)+'.condor.sub \n')
            i=i+1
        masterDAG_fp.close()
        print "Master dag is ready."
    #End __buildDAGfrom__

    def performFAsetup(self):
        buildDir(self.installPipes)
        buildDir(self.installIni)
        self.__createLHLL__()
        print "Creating each pipe configuration script."
        countMe=0
        modValue=20
        for pipe in self.FApipeNames:
            self.FAconfigure(pipe[3],pipe[4],pipe[0])
            if countMe%modValue==0:
                print ".",
            countMe=countMe+1
        print " "
        print "Building ",self.FApipeNames.__len__()," pipes. This may take awhile."
        for pipe in self.FApipeNames:
            #False turns off possible interpretations.
            self.createPipe(pipe[0],False)
            if countMe%modValue==0:
                print ".",
            countMe=countMe+1
        print " "
        print "Ok. Finished"
        #Search the workspace to construct a huge parent DAG for submitting all the DAGs
        root2dags=self.installPipes
        dagFilename=self.home+'/FA_Tuning.dag'
        self.__buildDAGfrom__(dagFilename,root2dags)
    #end performFAsetup method

    def FAconfigure(self,LH,LL,iniName):
        """
        Takes the cp and manipulates it to make a config file object.
        This object is suited for the FA rate investigation and written
        to disk with specified filename
        """
        newCP=copy.deepcopy(self.cpIni)
        # New/Rewritten ini options
        uniqStr='FA_'+str(LH)+'_'+str(LL)+'_'
        newCP.set('filelayout','workpath',self.installPipes+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        # If injection section present strip it out!
        if newCP.has_section('tracksearchinjection'):
            newCP.remove_section('tracksearchinjection')
        # If candidatethreshold section present strip it out!
        if newCP.has_section('candidatethreshold'):
            newCP.remove_section('candidatethreshold')
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End FAconfigure method

    def createPipe(self,pipeIniName,injectPipeFlag):
        """
        This method actually invokes a shell to create the pipe in a unique
        directory under the main directory specified in [all] tuningHome
        which tailors each pipe ini file entry [filelayout] workpath
        True mean leave 4 injectionsection
        False means remove injectionsection
        """
        seglist=self.seglist
        path2bin=self.pipeBuilder
        if injectPipeFlag == True:
            commandLine=path2bin+' --file='+pipeIniName+' --segment_file='+seglist+' --invoke_inject'
        else:
            commandLine=path2bin+' --file='+pipeIniName+' --segment_file='+seglist
        commandOut=commands.getstatusoutput(commandLine)
        if commandOut[0] != 0:
            print " "
            print "There was a problem installing a pipe for ini file :",pipeIniName
            print " "
            print commandOut[1]
            print " "
    # End method createPipe

    def performFAcalc(self):
        """
        This method will parse through all the pipes candidate files to find
        the individual 
        P and L values which would alone yield the specified P,L value for
        the given parameter FAR=3 imples 3sigma or about 99% FA rate.  This
        option is in the tun file under section [false-alarm-calculate] FAR
        """
        #1) Load the RESULTS File for level 1 tseries data
        #2) Invoke the candidate class method to give P,L list for given
        #   sigma value
        #3) Record the LH,LL,P,L quadruple to a file FAR_Results.dat
        myFAR=self.mySigmaFA
        
        # Determine the pipe that should be installed based on tun file.
        # resultFiles=self.__findResultFiles__(self.installPipes)
        resultFiles=self.__findFiles__(self.installPipes,['_1.candidates','RESULTS'])
        print "Checking ",resultFiles.__len__()," total candidate files. This may take a while."
        auxoutData=[]
        outputData=[]
        countMe=0
        modValue=20
        for entry in resultFiles:
            if countMe%modValue==0:
                print ".",
            countMe=countMe+1
            #Revising to avoid opening data files
            candidate=candidateList()
            myStat == candidate.candidateStatFromFile(entry)
            if myStat == []:
                meanP=varP=stdP=meanL=stdL=0
                threshP=0
                threshL=3
            else:
                meanP=myStat[3]
                varP=myStat[4]
                stdP=math.sqrt(varP)
                meanL=myStat[5]
                varL=myStat[6]
                stdL=math.sqrt(varL)
                #Threshold FAR values
                threshP=meanP+(myFAR*stdP)
                threshL=int(round(float((meanL+(myFAR*stdL)))))
            myOpts=str(str(entry).replace(self.installPipes,'').split('/')[1]).split('_')
            myLH=myOpts[1]
            myLL=myOpts[2]
            auxoutData.append([float(myLH),float(myLL),meanL,stdL,meanP,stdP])
            outputData.append([float(myLH),float(myLL),float(threshP),float(threshL)])
        auxout_fp=open(self.home+'/FA_results.aux','w')
        auxout_fp.write("LH,LL,Mean L,Stddev L,Mean P,Stddev P\n")
        for entry in auxoutData:
            auxout_fp.write(str(entry)+'\n')
        auxout_fp.close()
        output_fp=open(self.home+'/FA_results.dat','w')
        for entry in outputData:
            output_fp.write(str(entry)+'\n')
        output_fp.close()
        pickleout_fp=open(self.home+'/FA_results.pickle','w')
        pickle.dump(outputData,pickleout_fp)
        pickleout_fp.close()
    #End performFAcalc method

    def readFAresults(self):
        """
        Parse the results file to setup the appropriate P,L thresholds to do
        a detection efficiency study.  It uses pickle module to open up the file.
        """
        rawText=[]
        importFile=self.home+'/FA_results.pickle'
        if not(os.path.isfile(importFile)):
            print "Did you remember to run the false alarm calculate feature? ./tuneTSpipe.py --file=file.tun -c"
            os.abort()
        import_fp=open(importFile,'r')
        FAresults=pickle.load(import_fp)
        import_fp.close()
        return FAresults
    #End readFAresults method
    
    def performDEsetup(self):
        #1) Read in lambaH lambaL P L result file
        #2) Create corresponding ini file
        #3) Run pipe creation script to create injection pipes
        buildDir(self.installPipes2)
        buildDir(self.installIni2)
        myIni=self.myIni
        Orig_results=self.readFAresults()
        #Re-encode results to get one ini file for H,L with P and then H,L with Len
        results=[]
        for entry in Orig_results:
            h=entry[0]
            l=entry[1]
            p=entry[2]
            len=entry[3]
            #Setup P entry
            results.append([h,l,p,3])
            #Setup L entry
            results.append([h,l,0,len])
        print "The detection efficiency pipes to be created are",results.__len__()," this may take a while."
        countMe=0
        modValue=20
        for entry in results:
            if countMe%modValue==0:
                print ".",
            countMe=countMe+1
            h=entry[0]
            l=entry[1]
            p=entry[2]
            len=entry[3]
            pipeIniName=self.installIni2+'/'+self.batchMask+':'+str(h)+':'+str(float(l))+':'+str(float(p))+':'+str(int(len))+':'+'.ini'
            self.DEpipeNames.append(pipeIniName)
            self.DEeditIniFile(h,l,p,len,pipeIniName)
            #True states use the ini file to set injection into the pipeline.
            self.createPipe(pipeIniName,True)
        #Search the workspace to construct a huge parent DAG for submitting all the DAGs
        root2dags=self.installPipes2
        dagFilename=self.home+'/DE_Tuning.dag'
        self.__buildDAGfrom__(dagFilename,root2dags)
    #End performDEsetup

    def DEeditIniFile(self,LH,LL,P,L,iniName):
        """
        This method will edit the original ini file in such a way as to
        prepare the detection efficiency runs to be launched.
        """
        newCP=copy.deepcopy(self.cpIni)
        #New/Rewritten ini options
        uniqStr='DE_'+str(LH)+'_'+str(LL)+'_'+str(P)+'_'+str(L)+'_'
        newCP.set('filelayout','workpath',self.installPipes2+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        newCP.set('tracksearchbase','length_threshold',L)
        newCP.set('tracksearchbase','power_threshold',P)
        # The master Ini must have an injection section to continue
        if not newCP.has_section('tracksearchinjection'):
            print "Error the master ini in our tun file has no injection section!"
            os.abort
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End DEeditIniFile

    def NEW_DEeditIniFile(self,LH,LL,P,L,iniName):
        """
        This method will edit the original ini file in such a way as to
        prepare the detection efficiency runs to be launched.
        """
        newCP=copy.deepcopy(self.cpIni)
        #New/Rewritten ini options
        uniqStr='DE_'+str(LH)+'_'+str(LL)+'_'+str(P)+'_'+str(L)+'_'
        newCP.set('filelayout','workpath',self.installPipes2+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        #We adjust the options associated with the python post processing routines.
        expString="(P>%f)and(L>%i)"%(P,L)
        newCP,set('candidatethreshold','expression_threshold',expString)
        # The master Ini must have an injection section to continue
        if not newCP.has_section('tracksearchinjection'):
            print "Error the master ini in our tun file has no injection section!"
            os.abort
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End DEeditIniFile
    
    def performDEcalc(self):
        #count num of can files with 1 entry+
        #determine percentage succsessful
        #write 2 ranked 4c list lh ll #map percentage
#        globFiles=self.__findLayer1GlobFiles__(self.installPipes2)
        globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates'])
        outputPickle=[]
        outputP=[]
        outputL=[]
        numInjections=int(self.cpIni.get('tracksearchinjection','inject_count'))
        print "Calculating efficiencies for ",globFiles.__len__()," files."
        for entry in globFiles:
            candObj=candidateList()
            candObj.loadfile(entry)
            fileName=entry
            curveCount=candObj.totalCount
            myOpts=str(str(entry).replace(self.installPipes2,'').split('/')[1]).split('_')
            h=myOpts[1]
            l=myOpts[2]
            p=myOpts[3]
            len=myOpts[4]
            eff=float(candObj.totalCount)/float(numInjections)
            outputPickle.append([h,l,p,len,eff,candObj.totalCount,numInjections,entry])
            outputP.append([float(h),float(l),float(p),float(eff),int(candObj.totalCount),int(numInjections)])
            outputL.append([float(h),float(l),int(float(len)),float(eff),int(candObj.totalCount),int(numInjections)])
        #Write out results files to disk
        format4CL = "%10.5f  %10.5f  %15i    %10.5f %10i %10i\n"
        format4CP = "%10.5f  %10.5f  %10.5f  %10.5f %10i %10i\n"
        format6C = "%10.5f %10.5f  %10.5f %10i  %10.5f %10i %10i\n"
        file_fp=open(self.home+'/DE_results_L.dat','w')
        file_fp.write("H L P Len Eff Count \n")
        for entry in outputL:
            file_fp.write(format4CL % (entry[0],entry[1],entry[2],entry[3],entry[4],entry[5]))
        file_fp.close()
        file_fp=open(self.home+'/DE_results_P.dat','w')
        file_fp.write("H L P Len Eff Count \n")
        for entry in outputP:
            file_fp.write(format4CP % (entry[0],entry[1],entry[2],entry[3],entry[4],entry[5]))
        file_fp.close()
        file_fp=open(self.home+'/DE_results_Mix.dat','w')
        file_fp.write("H L P Len Eff Count \n")
        for entry in outputPickle:
            file_fp.write(format6C % (float(entry[0]),
                                      float(entry[1]),
                                      float(entry[2]),
                                      int(float(entry[3])),
                                      float(entry[4]),
                                      int(float(entry[5])),
                                      int(float(entry[6]))))
        file_fp.close()
        pickle_fp=open(self.home+'/DE_results.pickle','w')
        pickle.dump(outputPickle,pickle_fp)
        pickle_fp.close()
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

if FAsetup: #ready
    countOpts=countOpts+1
if DEsetup: #NOT READY
    countOpts=countOpts+1
if FAcalc: #NOT READY
    countOpts=countOpts+1
if DEcalc: #NOT READY
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
