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
This is a complete standalone python script which will acess the
search ini file and call the tracksearch_pipe.py program to setup all
the pipes needed to tune the search algorithm.  This script will also
require its own ini file to setup the analysis tuning.  The convention
we want to follow for tuning ini files to give them the extens .tun
"""
import sys
from glue import pipeline
from optparse import OptionParser
import ConfigParser
import commands
import copy
import getopt
import math
import os
import string
import time
import pickle
try:
    from lalapps import tracksearch
except:
    import tracksearch
try:
    from lalapps.tracksearchutils import *
except:
    from tracksearchutils import *
try:
    from lalapps.tracksearch import buildDir
except:
    from tracksearch import buildDir


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
        if not os.path.isfile(self.masterIni):
            print "Error with masterIni in tun configuration file!"
            print self.masterIni
            os.abort()
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
        self.installPipes=os.path.normpath(self.home+'/FA_pipes')
        self.installIni=os.path.normpath(self.home+'/FA_ini')
        self.installPipes2=os.path.normpath(self.home+'/DE_pipes')
        self.installIni2=os.path.normpath(self.home+'/DE_ini')
        self.log=cp.get('all','tuningLogs')
        self.dagpath=cp.get('all','tuningDags')
        self.seglist=cp.get('all','seglist')
        self.mySigmaFA=float(cp.get('false-alarm-calculate','FAR'))
        self.FApipeNames=[]
        self.DEpipeNames=[]
        self.pipeBuilder=cp.get('all','pipeProgram')
        #Set the pickle file from Curves Found in FA_pipe run.
        self.curveFoundPickle=[]
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

    def __findDirs__(self,searchPath,dirName,shortCircuit):
        """
        Finds any path that terminates its name with dirName.  This is an
        exact search ie find != Find and no partial text searches as in
        __findFiles__. If a path can be built like /searchPath/dirName  if shortCircuit
        is true we stop the downward recursion.  Can not fine /searchPath/dirName/Layer/dirName
        directories with shortCircuit equal True!
        """
        countMe=0
        modValue=1
        dirsInPath=[]
        for root,dir,files in os.walk(searchPath):
            for subDir in dir:
                tmpPath=os.path.abspath(os.path.normpath(root+'/'+subDir))
                if tmpPath.endswith(dirName):
                    dirsInPath.append(tmpPath)
                    if countMe%modValue==0:
                        sys.stdout.writelines('.')
                        sys.stdout.flush()
                        countMe=countMe+1
            if shortCircuit == True:
                while dir.__contains__(dirName):
                    dir.remove(dirName)
        return dirsInPath
    #End __findDirs__(self,searchPath,dirName)

    def __scanDirectoryCandidateFiles__(self,entry):
        """
        Looks at all *.candidate files tallying up the num of files with an
        entry
        """
        pathString=os.path.normpath(entry+'/*.candidates')
        cmdString='cd '+entry+' ;grep -m 1 "Curve number" *.candidates | wc -l'
        [err,fileCountTXT]=commands.getstatusoutput(cmdString)
        if fileCountTXT.isdigit():
            triggerCount=int(fileCountTXT)
        else:
            triggerCount=0
        return triggerCount
    #End __scanDirectoryCandidateFiles__(entry)

    def __determineDEinPath__(self,candPath):
        """
        This uses the shell and grep to count the number of candidate files
        with a trigger in some particular path.  NOT CONSIDERED RELIABLE!
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
        fileList=self.__findFiles__(root2dags,['.dag'])
        print "The master DAG will have ",fileList.__len__()," nodes/subdags."
        print "You should find the master dag contained in:",root2dags
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

    def __writeFAfiles__(self,auxoutData,outputData):
        """
        Common code to write out FA_results file.  This facilitates ease
        of use for FA routines. Expects two input lists [] an empty list
        will not be written out.
        """
        if auxoutData !=[]:
            auxout_fp=open(self.home+'/FA_results.aux','w')
            auxout_fp.write("LH,LL,Mean L,Stddev L,Mean P,Stddev P Count\n")
            txt=""
            for entry in auxoutData:
                for element in entry:
                    txt=txt+str("%15.5f")%(float(element))
                entryTXT=txt
                txt=""
                #entryTXT=str(entry).replace('[','').replace(']','').replace("'",'').replace(',',' ')
                auxout_fp.write(str(entryTXT)+'\n')
            auxout_fp.close()
        if outputData != []:
            output_fp=open(self.home+'/FA_results.dat','w')
            for entry in outputData:
                output_fp.write(str(entry)+'\n')
            output_fp.close()
            pickleout_fp=open(self.home+'/FA_results.pickle','w')
            pickle.dump(outputData,pickleout_fp)
            pickleout_fp.close()
            output_fp=open(self.home+'/FA_results_MatlabPlot.dat','w')
            formatContour = "%10.5f %10.5f %10.5f %10i\n"
            contourData=outputData
            for entry in contourData:
                output_fp.write(formatContour % (float(entry[0]),
                                                 float(entry[1]),
                                                 float(entry[2]),
                                                 int(entry[3])))
            output_fp.close()
            formatContour = "%10.5f %10.5f %10.5f %10i\n"
            output_fp=open(self.home+'/FA_results_CurvesFound_MatlabPlot.dat','w')
            lineData=[]
            for entry in auxoutData:
                lineData.append([entry[0],entry[1],entry[6],entry[6]])
                output_fp.write(formatContour % (float(entry[0]),
                                                 float(entry[1]),
                                                 float(entry[6]),
                                                 int(entry[6])))
        pickleout_fp=open(self.home+'/FA_results_Found.pickle','w')
        pickle.dump(lineData,pickleout_fp)
        pickleout_fp.close()
        output_fp.close()
    #End __writeFAfiles__()
    
    def getBackground(self,LH,LL,curvePickleFile):
        """
        Uses specified pickle file to try and find the curves found at a given
        LH and LL.  This is used as the off-source background rate.
        """
        if self.curveFoundPickle == []:
            print "Pickle unopen trying to load:",curvePickleFile
            print "Returning 1% of seen noise background rate, due to expected Z P,L choices."
            if not(os.path.isfile(curvePickleFile)):
                print "Did you remember to run the false alarm calculate feature? ./tuneTSpipe.py --file=file.tun -c"
                os.abort()
            input_fp=file(curvePickleFile,'r')
            self.curveFoundPickle=pickle.load(input_fp)
            self.curveFoundPickle.sort()
            input_fp.close()
            print "Done loading pickle jar."
        start=0
        stop=self.curveFoundPickle.__len__()
        current=int(round(stop/2))
        found=False
        tries=0
        #Ini file has z=4.5 is .99999 or a background of 0.00001
        percentage=0.000001
        while not(found):
            #if tries > int(self.curveFoundPickle.__len__()/4):
            if tries > 10:
                "Error looking up parameters."
                os.abort()
            tries=tries+1
            entry=self.curveFoundPickle[current]
            [lh,ll]=[entry[0],entry[1]]
#            print "C,F,H,==,<,>",current,"[",LH,LL,"]","[",lh,ll,"]",[LH,LL]==[lh,ll],[LH,LL]<[lh,ll],[LH,ll]>[lh,ll]
            if [LH,LL]==[lh,ll]:
                found=True
                return int(percentage*entry[2])
            elif ((current==start) or (current==stop)):
                return int(0)
            elif [LH,LL]<[lh,ll]:
                tmp=current
                current=current-int(round((current-start)/2.0))
                stop=tmp
            elif [LH,LL]>[lh,ll]:
                tmp=current
                current=current+int(round((stop-current)/2.0))
                start=tmp
            if current<start or current>stop:
                return int(0)
        print "This error should never happen!"
        return int(0)
    #End getBackground()
    
    def performFAsetup(self):
        buildDir(self.installPipes)
        buildDir(self.installIni)
        self.__createLHLL__()
        print "Creating each pipe configuration script."
        countMe=0
        modValue=10
        for pipe in self.FApipeNames:
            self.FAconfigure(pipe[3],pipe[4],pipe[0])
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
        print " "
        print "Building ",self.FApipeNames.__len__()," pipes. This may take awhile."
        countMe=0
        for pipe in self.FApipeNames:
            #False turns off possible injections.
            self.createPipe(pipe[0],False)
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
        print " "
        print "Ok. Finished"
        #Search the workspace to construct a huge parent DAG for submitting all the DAGs
        root2dags=os.path.normpath(self.dagpath+'/FA/')
        dagFilename=os.path.normpath(self.dagpath+'/FA/FA_tuning.dag')
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
        newCP.set('filelayout','workpath',os.path.normpath(self.installPipes+'/'+uniqStr))
        newCP.set('filelayout','logpath',os.path.normpath(self.log+'/'+uniqStr))
        newCP.set('filelayout','dagpath',os.path.normpath(self.dagpath+'/FA/'+uniqStr))
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
        contourData=[]
        countMe=0
        modValue=10
        for entry in resultFiles:
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
            #Revising to avoid opening data files
            candidate=candidateList()
            #print "Checking ",entry
            myStat=candidate.candidateStatsOnDisk(entry)
            if myStat == []:
                curves=meanP=varP=stdP=meanL=stdL=0
                threshP=0
                threshL=3
            else:
                curves=myStat[2]
                meanP=myStat[3]
                varP=myStat[4]
                stdP=math.sqrt(varP)
                meanL=myStat[5]
                varL=myStat[6]
                stdL=math.sqrt(varL)
                maxP=myStat[7]
                maxL=myStat[8]
                #Threshold FAR values
                threshP=meanP+(myFAR*stdP)
                threshL=int(round(float((meanL+(myFAR*stdL)))))
                if threshP < maxP:
                    #print "Threshold seems inappropriate:",threshP,maxP
                    #print "Resetting IP threshold to 1 event in data set!"
                    threshP=maxP
                if threshL < maxL:
                    #print "Threshold seems inappropriate:",threshL,maxL
                    #print "Resetting CL threshold to 1 event in data set!"
                    threshL=maxL
            myOpts=str(str(entry).replace(self.installPipes,'').split('/')[1]).split('_')
            myLH=myOpts[1]
            myLL=myOpts[2]
            #Ignore lines where both Pval and Lval are zero
            if ((threshP !=0) and (threshL != 0)):
                auxoutData.append([float(myLH),float(myLL),float(meanL),float(stdL),float(meanP),float(stdP),int(curves)])
                outputData.append([float(myLH),float(myLL),float(threshP),float(threshL)])
                contourData.append([float(myLH),float(myLL),float(threshP),float(threshL)])
        print " "
        self.__writeFAfiles__(auxoutData,outputData)
    #End performFAcalc method

    def rewriteFARfile(self):
        """
        This method uses the output of performFAcalc to re-write the threshold files this is in case we want to change the FAR value and rerun the detection efficiency tests
        """
        #Load up the statistic aux file.
        auxout_fp=open(self.home+'/FA_results.aux','r')
        auxTxt=auxout_fp.readlines()
        auxout_fp.close()
        outputData=[]
        myFAR=self.mySigmaFA
        headerTxt=auxTxt.pop(0)
        #Recalculate stats        
        for entry in auxTxt:
            (myLH,myLL,meanL,stdL,meanP,stdP,curves)=entry.split()
            threshP=float(meanP)+(myFAR*float(stdP))
            threshL=int(round(float((float(meanL)+(myFAR*float(stdL))))))
            outputData.append([float(myLH),float(myLL),float(threshP),float(threshL)])
        #Rewrite the threshold pickle file and dat file
        #
        self.__writeFAfiles__([],outputData)
    #End rewriteFARfile
    
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
        modValue=10
        for entry in results:
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
            h=entry[0]
            l=entry[1]
            p=entry[2]
            len=entry[3]
            pipeIniName=self.installIni2+'/'+self.batchMask+':'+str(h)+':'+str(float(l))+':'+str(float(p))+':'+str(int(len))+':'+'.ini'
            self.DEpipeNames.append(pipeIniName)
            #Check to see if the threshold section exits...
            #If it does we implicity keep all triggers for FA and DE
            #Else the DE runs try to save disk space and allow tracksearch
            #to perform our simple L and P thresholds. Check MasterIni
            if self.cpIni.has_section('candidatethreshold'):
                self.keepTrigs_DEeditIniFile(h,l,p,len,pipeIniName)
            else:
                self.saveDisk_DEeditIniFile(h,l,p,len,pipeIniName)
            #True states use the ini file to set injection into the pipeline.
            self.createPipe(pipeIniName,True)
        #Search the workspace to construct a huge parent DAG for submitting all the DAGs
        root2dags=os.path.normpath(self.dagpath+'/DE/')
        dagFilename=os.path.normpath(self.dagpath+'/DE/DE_tuning.dag')
        self.__buildDAGfrom__(dagFilename,root2dags)
    #End performDEsetup

    def saveDisk_DEeditIniFile(self,LH,LL,P,L,iniName):
        """
        This method will edit the original ini file in such a way as to
        prepare the detection efficiency runs to be launched.
        Do not use this method!
        """
        newCP=copy.deepcopy(self.cpIni)
        #New/Rewritten ini options
        uniqStr='DE_'+str(LH)+'_'+str(LL)+'_'+str(P)+'_'+str(L)+'_'
        newCP.set('filelayout','workpath',self.installPipes2+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('filelayout','dagpath',os.path.normpath(self.dagpath+'/DE/'+uniqStr))
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

    def keepTrigs_DEeditIniFile(self,LH,LL,P,L,iniName):
        """
        This method will edit the original ini file in such a way as to
        prepare the detection efficiency runs to be launched.
        """
        newCP=copy.deepcopy(self.cpIni)
        #New/Rewritten ini options
        uniqStr='DE_'+str(LH)+'_'+str(LL)+'_'+str(P)+'_'+str(L)+'_'
        newCP.set('filelayout','workpath',self.installPipes2+'/'+uniqStr)
        newCP.set('filelayout','logpath',self.log+uniqStr)
        newCP.set('filelayout','dagpath',os.path.normpath(self.dagpath+'/DE/'+uniqStr))
        newCP.set('tracksearchbase','start_threshold',LH)
        newCP.set('tracksearchbase','member_threshold',LL)
        #We adjust the options associated with the python post processing routines.
        expString="(P>%f)and(L>%i)"%(P,L)
        newCP.set('candidatethreshold','expression_threshold',expString)
        # The master Ini must have an injection section to continue
        if not newCP.has_section('tracksearchinjection'):
            print "Error the master ini in our tun file has no injection section!"
            os.abort
        # Write out this modified ini file to disk
        fp=open(iniName,'w')
        newCP.write(fp)
        fp.close()
    #End DEeditIniFile

    def __performShellDEcalc__(self):
        """
        Relys on quickly check /1/ directions and counting number of canidate
        files. Then it counts the number of files with 1 or more entries.
        From this we determine the detection efficiency.  Though we won't know
        anything about the triggers themselves.  To work successfully we need
        to execute the search with masterIni [housekeeper] to keep *.candidate
        files in the /1/ directory!  We assume that the trials are configured
        at one injection per map! This is critical, see masterIni
        [tracksearchinjection] to verify this!
        """
        print "Remember we assumed that masterIni configured with 1 injection per map!"
        print "Options in [tracksearchinjection]"
        for entry in self.cpIni.options('tracksearchinjection'):
            print entry," ",self.cpIni.get('tracksearchinjection',entry)
        print ""
        if self.cpIni.has_section('layerconfig'):
            topBlockSize=int(self.cpIni.get('layerconfig','layerTopBlockSize'))
            timeScale=float(self.cpIni.get('layerconfig','layer1TimeScale'))
            deltaT=timeScale
            mapCounts=topBlockSize/timeScale
        else:
            print "Error missing required section LAYERCONFIG"
            os.abort
        if self.cpIni.has_section('candidatethreshold'):
            print "Ok.  There is a problem. Please review your tuning setup.  You have a candidateThreshold section in the masterIni file.  We can't do python thresholding with this tuning subroutine."
            print "Trying alternate tuning routine: self.__performMapDEcalc__"
            print "BE WARNED THIS ALTERNATIVE ROUTINE IS REALLY REALLY SLOW!!!"
            print "Please see tuneTSpipe.py comments for explaination."
            [outputPickle,outputP,outputL,outputContour]=self.__performMapDEcalc__()
            return [outputPickle,outputP,outputL,outputContour]
        else:
            globFiles=[]
            #Waste of CPU cycles to guess the number of files to process
            #globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates','DE_'])
        #Scan for directories /1/ labeled
        print "Checking for output files in tuning directory."
        shortCircuit=True #Stops tranversing once a path ending in arg2 is determined!
        foundDirs=self.__findDirs__(self.installPipes2,'1',shortCircuit)
        print "\n Calculating efficiencies for approximately ",foundDirs.__len__()," trials."
        countMe=0
        modValue=10
        outputPickle=[]
        outputL=[]
        outputP=[]
        outputContour=[]
        refTime=time.time()
        debugFile=file(str(os.path.basename(tunFile))+'_debug.info','w');
        for entry in foundDirs:
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
            myOpts=str(str(entry).replace(self.installPipes2,'').split('/')[1]).split('_')
            #Write each directory scanned to a file.  Debug purposes...
            deltaTime="%3.2f"%(time.time()-refTime)
            debugFile.write(str(deltaTime)+' '+str(entry)+'\n')
            debugFile.flush()
            if myOpts.__len__() > 4:
                triggerCount=self.__scanDirectoryCandidateFiles__(entry)
                h=myOpts[1]
                l=myOpts[2]
                p=myOpts[3]
                len=myOpts[4]
                eff=float(triggerCount/mapCounts)
                curveCount=int(triggerCount)
                outputPickle.append([h,l,p,len,eff,mapCounts,0,0])
                if float(p)>0:
                    outputP.append([float(h),float(l),float(p),float(eff),int(curveCount),int(mapCounts)])
                if float(len)>3:
                    outputL.append([float(h),float(l),int(float(len)),float(eff),int(curveCount),int(mapCounts)])
                outputContour.append([float(h),float(l),float(p),int(curveCount)])
            else:
                print "Omit::Problem checking results for :",entry
                print "Found opts::",myOpts
        #End the loop
        print " "
        debugFile.close()
        return [outputPickle,outputP,outputL,outputContour]
    #end __performShellDEcalc__()
        
    def __performMapDEcalc__(self):
        """
        Checks for triggers using a map definition from GWDAW9 paper.
        Takes the ini file information to count the number of maps that
        were triggered using the Glob file for this information.
        Take summary information from file and process it using the masterini file information
        Assume that there is only 1 valid injection per map/interval of interest.
        """
        if self.cpIni.has_section('layerconfig'):
            topBlockSize=int(self.cpIni.get('layerconfig','layerTopBlockSize'))
            timeScale=float(self.cpIni.get('layerconfig','layer1TimeScale'))
            deltaT=timeScale
            mapCounts=topBlockSize/timeScale
        else:
            print "Error missing required section LAYERCONFIG"
            os.abort
        if self.cpIni.has_section('candidatethreshold'):
            globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates','DE_','Threshold:'])
        else:
            globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates','DE_'])
        print "Calculating efficiencies for ",globFiles.__len__()," files."
        countMe=0
        modValue=10
        outputPickle=[]
        outputL=[]
        outputP=[]
        outputContour=[]
        numInjections=int(self.cpIni.get('tracksearchinjection','inject_count'))
        numInjections=0
        for entry in globFiles:
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            myOpts=str(str(entry).replace(self.installPipes2,'').split('/')[1]).split('_')
            globStartTime=float(str(os.path.basename(entry)).split(':')[2])
            h=myOpts[1]
            l=myOpts[2]
            p=myOpts[3]
            len=myOpts[4]
            #Load glob file
            candObj=candidateList()
            candObj.loadfile(entry)
            fileSummary=candObj.dumpCandidateKurveSummary()
            fileSummary.sort()
            del candObj
            #Create Summary
            injectCount=0
            currentTime=globStartTime
            #Check through Summary and record answer
            indexLow=0
            indexHigh=0
            while currentTime < globStartTime+topBlockSize:
                while fileSummary[indexLow]<currentTime:
                    indexLow=indexLow+1
                while fileSummary[indexHigh]<(currentTime+deltaT):
                    indexHigh=indexHigh+1
                #Take the difference this is num of triggers between currentTime and cT+dT
                trigFound=indexHigh-indexLow
                if trigFound>0:
                    injectCount=injectCount+1
                currentTime=currentTime+deltaT
            eff=injectCount/mapCounts
            curveCount=int(injectCount)
            outputPickle.append([h,l,p,len,eff,mapCounts,0,0])
            if float(p)>0:
                outputP.append([float(h),float(l),float(p),float(eff),int(curveCount),int(numInjections)])
            if float(len)>3:
                outputL.append([float(h),float(l),int(float(len)),float(eff),int(curveCount),int(numInjections)])
            outputContour.append([float(h),float(l),float(p),int(curveCount)])
        #End the loop
        print " "
        return [outputPickle,outputP,outputL,outputContour]
    #end __performMapDEcalc__()

    def __performNonMapDEcalc__(self):
        """
        This counts just triggers from the glob not considering our orignal trigger definition.
        """
        #count num of can files with 1 entry+
        #determine percentage succsessful
        #write 2 ranked 4c list lh ll #map percentage
        #globFiles=self.__findLayer1GlobFiles__(self.installPipes2)
        if self.cpIni.has_section('candidatethreshold'):
            globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates','DE_','Threshold:'])
        else:
            globFiles=self.__findFiles__(self.installPipes2,['Glob','_1.candidates','DE_'])
        outputPickle=[]
        outputP=[]
        outputL=[]
        outputContour=[]
        numInjections=int(self.cpIni.get('tracksearchinjection','inject_count'))
        numInjections=0
        print "Calculating efficiencies for ",globFiles.__len__()," files."
        countMe=0
        modValue=10
        for entry in globFiles:
            if countMe%modValue==0:
                sys.stdout.writelines('.')
                sys.stdout.flush()
            countMe=countMe+1
            candObj=candidateList()
            candStats=candObj.candidateStatsFromFile(entry)
            fileName=entry
            curveCount=candStats[2]
            myOpts=str(str(entry).replace(self.installPipes2,'').split('/')[1]).split('_')
            h=myOpts[1]
            l=myOpts[2]
            p=myOpts[3]
            len=myOpts[4]
            eff=curveCount
            outputPickle.append([h,l,p,len,eff,curveCount,numInjections,entry])
            if float(p)>0:
                outputP.append([float(h),float(l),float(p),float(eff),int(curveCount),int(numInjections)])
            if float(len)>3:
                outputL.append([float(h),float(l),int(float(len)),float(eff),int(curveCount),int(numInjections)])
            outputContour.append([float(h),float(l),float(p),int(curveCount)])
        print " "
        return [outputPickle,outputP,outputL,outputContour]
    #End __performNonMapDEcalc__()
        
    def performDEcalc(self):
        #[outputPickle,outputP,outputL,outputContour]=self.__performNonMapDEcalc__()
        [outputPickle,outputP,outputL,outputContour]=self.__performShellDEcalc__()
        #Write out results files to disk
        format4CL = "%10.5f  %10.5f  %15i    %10.5f %10i %10i\n"
        format4CP = "%10.5f  %10.5f  %10.5f  %10.5f %10i %10i\n"
        format6C = "%10.5f %10.5f  %10.5f %10i  %10.5f %10i %10i\n"
        formatContour = "%10.5f %10.5f %10.5f %10.5f\n"
        file_fp=open(self.home+'/DE_results_L.dat','w')
        file2_fp=open(self.home+'/DE_results_Matlab_Length.dat','w')
        file3_fp=open(self.home+'/DE_results_Matlab_Length_NoBack.dat','w')
        file_fp.write("H L P Len Eff Count \n")
        for entry in outputL:
            file_fp.write(format4CL % (entry[0],entry[1],entry[2],entry[3],entry[4],entry[5]))
            background=self.getBackground(entry[0],entry[1],self.home+'/FA_results_Found.pickle')
            #Omit zero result entries!
            if float(entry[4])>0:
                file2_fp.write(formatContour % (entry[0],entry[1],entry[2],entry[3]))
            file3_fp.write(formatContour % (entry[0],entry[1],entry[2],entry[3]-background))
        file_fp.close()
        file2_fp.close()
        file3_fp.close()
        file_fp=open(self.home+'/DE_results_P.dat','w')
        file1_fp=open(self.home+'/DE_results_Matlab_Power.dat','w')
        file3_fp=open(self.home+'/DE_results_Matlab_Power_NoBack.dat','w')
        file_fp.write("H L P Len Eff Count \n")
        for entry in outputP:
            file_fp.write(format4CP % (entry[0],entry[1],entry[2],entry[3],entry[4],entry[5]))
            background=self.getBackground(entry[0],entry[1],self.home+'/FA_results_Found.pickle')
            background=int(background)
            if float(entry[4])>0:
                file1_fp.write( formatContour %(entry[0],entry[1],entry[2],entry[3]))
            file3_fp.write( formatContour %(entry[0],entry[1],entry[2],entry[3]-background))
        file_fp.close()
        file1_fp.close()
        file3_fp.close()
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
        newOutputContour=list()
        lList=[]
        pList=[]
        tmpList=copy.deepcopy(outputContour)
        tmpList.sort()
        file_fp=open(self.home+'/DE_results_MatlabPlot.dat','w')
        for entry in tmpList:
            file_fp.write(formatContour % (float(entry[0]),
                                           float(entry[1]),
                                           float(entry[2]),
                                           int(float(entry[3]))))
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
usage="tstunepipe [args]"
parser = OptionParser(usage)

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
parser.add_option("-r","--false-alarm-recalculate",dest="FAagain",
                  default=False,
                  action="store_true",
                  help="Invoke this option to recalculate IP and CL thresholds if the false alarm rate has been altered in the tun file.")
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
FAagain=bool(options.FAagain)
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
if FAagain:
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
elif FAagain:
    tuner.rewriteFARfile()
else:
    print "This should never happen!"
    os.abort()
#
# End main program
#
 
