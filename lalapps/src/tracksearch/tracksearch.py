"""
Classes and methods for the tracksearch pipeline
"""

__author__ = 'Charlie Torres <charlie@phys.utb.edu>'
__date__ = '$Date$'
__version__ = ''

import os
import re
import time
import string
import math
import exceptions
import ConfigParser
from glue import pipeline


#Method to Build directory which is not present
def buildDir(dirpath):
    pathArray=[]
    stopVar=0
    #Check for trailing slash
    oldDir=str(os.path.abspath(dirpath)).__add__('/')
    while stopVar == 0:
        splitContainer=os.path.split(oldDir)
        newSec=splitContainer[0]
        pathArray.append(newSec)
        oldDir=newSec
        if ((splitContainer[1] == '') and (splitContainer[0] == '/')):
            stopVar=1
    #Now reverse this list and recurse creating the pathing as needed
    pathArray.reverse()
    for pathBlock in pathArray:
        if os.path.isfile(pathBlock):
            print 'We have a problem the path requested may clobber an existing file!'
            print 'File :',str(pathBlock)
            os.abort()
        if not os.path.isdir(pathBlock):
           print 'Making path:',pathBlock
           try: os.mkdir(pathBlock)
           except: pass
    
    #Now double check that directory exists else flag error
    if not os.path.isdir(dirpath):
        print '***'
        print 'SHIT, directory still not ready contain(s) ',pathArray.__len__(),' branches!'
        print dirpath
        print '+++'
        print pathArray
        print '***'
#End def

def isPowOf2(input):
    inNum=int(input)
    ans=inNum&(inNum-1)
    if ans== 0:
        return 1
    else:
        return 0
#End def

def powOf2Floor(input):
    ans=0
    inNum=int(input)
    if inNum < 0:
        print "ERROR Unexpected negative value found!"
    newPow2=1
    expPow=1
    if not isPowOf2(inNum):
        while newPow2 < inNum:
            newPow2=math.pow(2,expPow)
            expPow+=1
        expPow+=-1
        ans=math.pow(2,expPow)
    else:
        ans=inNum
    return ans

def determineLayerPath(cp,blockID,layerID):
    layerPath=str('%s/%s/%s/'%(cp.get('fileLayout','workpath'),blockID,layerID))
    return layerPath
#End def

def determineBlockPath(cp,blockID):
    blockPath=str('%s/%s/'%(cp.get('fileLayout','workpath'),blockID))
    return blockPath
#End def

class tracksearchCheckIniFile:
    """
    This class will check the ini file for appropriate arguements
    where possible.  The errors will be listed upon request
    """

    def __init__(self,cp):
        self.iniOpts=cp
        self.errList=[]
    #End init

    def checkOpts(self):
        #Check for existance of [condor] section files
        condorOptList=self.iniOpts.options('condor')
        for entry in condorOptList:
            optValue=self.iniOpts.get('condor',entry)
            if str(optValue).__contains__('/'):
                if not os.path.exists(str(optValue)):
                    self.errList.append('Can not find :'+str(entry)+':'+str(optValue))

        #Check [tracksearchbase] section
        lambdaH=0
        lambdaL=1
        lambdaH=self.iniOpts.get('tracksearchbase','start_threshold')
        lambdaL=self.iniOpts.get('tracksearchbase','member_threshold')
        if lambdaL > lambdaH :
            self.errList.append('Error invalid choice for lambda parameters')
        CL=float(self.iniOpts.get('tracksearchbase','length_threshold'))
        if CL < 3:
            self.errList.append('Error invalid choice for curve length threshold')
        IP=float(self.iniOpts.get('tracksearchbase','power_threshold'))
        if IP <= 0:
            self.errList.append('Error invalid choice for integrated power threshold')
        #Check [tracksearchtime] section
        WS=1
        NOFB=0
        WS=int(self.iniOpts.get('tracksearchtime','window_size'))
        NOFB=int(self.iniOpts.get('tracksearchtime','number_of_freq_bins'))
        if WS > NOFB:
            self.errList.append('Error window length inconsistent!')
        #Check [layerconfig] section
        LTBS=self.iniOpts.get('layerconfig','layerTopBlockSize')
        layerOpts=self.iniOpts.options('layerconfig')
        layerOpts.sort()
        layerTimes=[]
        layerTimesOrig=[]
        for entry in layerOpts:
            optValue=float(self.iniOpts.get('layerconfig',entry))
            if str(entry).__contains__('TimeScale'):
                layerTimes.append(optValue)
                layerTimesOrig.append(optValue)
        layerTimes.sort()
        if layerTimesOrig != layerTimes:
            self.errList.append('Error inconsistent layer time scales!')

    #end checkOpts def

    def numberErrors(self):
        return self.errList.__len__()
    #end numberErrors def

    def printErrorList(self):
        print self.numberErrors(),' Errors found!'
        for error in self.errList:
            print error
    #end printErrorList
#end Class

            
class tracksearchConvertSegList:
    """
    This class declaration is responsible for opening some segment
    file and reparsing it into block sizes of exactly n seconds each
    from data in the original segments file.  It is this file from where we
    construct an analysis DAG for each line item.
    """
    def __init__(self,segmentFileName,newDuration):
        tracksearchConvertSegList.origSegObject=pipeline.ScienceData()
        tracksearchConvertSegList.newSegObject=pipeline.ScienceData()
        tracksearchConvertSegList.newSegFilename=segmentFileName+'.revised'
        tracksearchConvertSegList.segFilename=segmentFileName
        tracksearchConvertSegList.duration=newDuration

    def writeSegList(self):
        #Read segents from file
        tracksearchConvertSegList.origSegObject.read(tracksearchConvertSegList.segFilename,tracksearchConvertSegList.duration)
        #Read the file header
        input_fp=open(tracksearchConvertSegList.segFilename,'r')
        newHeading=[]
        for headline in input_fp.readlines():
            if str(headline)[0] == '#':
                 newHeading.append(headline)
        newHeading.append('# This file is a reprocessed list at 100s intervals.\n')
        newHeading.append('# We drop sections were there is not 100s of data\n')
        newHeading.append('# This file drawn from '+tracksearchConvertSegList.segFilename+'\n')
        output_fp=open(tracksearchConvertSegList.newSegFilename,'w')
        for newHeadline in newHeading:
            output_fp.write(newHeadline)
        #Write redivided list to file
        index=0
        for bigChunks in tracksearchConvertSegList.origSegObject:
            bigChunks.make_chunks(tracksearchConvertSegList.duration)
            for newChunk in bigChunks:
                index+=1
                label="%d %d %d %d\n"%(index,newChunk.start(),newChunk.end(),newChunk.end()-newChunk.start())
                output_fp.write(label)
        output_fp.close

    def getSegmentName(self):
        #Return a string containing the name of the revised segment list
        return(tracksearchConvertSegList.newSegFilename)


class tracksearchTimeJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    The class running a TIME instance of the tracksearch code
    """
    def __init__(self,cp,block_id,layer_id,dagDir):
        self.dagDirectory=dagDir
        #ConfigParser object -> cp
        self.__executable = cp.get('condor','tracksearch')
        self.__universe = cp.get('condor','universe');
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp)
        blockID=block_id
        layerID=layer_id
        layerPath=determineLayerPath(cp,blockID,layerID)
        blockPath=determineBlockPath(cp,blockID)
        #THERE IS A PROBLEM WITH GPS START ARGUEMENT TO THE CODE
        #Parse all the options from BASE and Timeseries Config sections
        for sec in ['tracksearchbase']:
            self.add_ini_opts(cp,sec)
        for sec in ['tracksearchtime']:
            self.add_ini_opts(cp,sec)
        #Read expected job sampling rate
        sampleRate=float(cp.get('tracksearchtime','sample_rate'))
        #Read expected TF overlapping percentage
        overlapPercentage=float(cp.get('layerconfig','layerOverlapPercent'))
        #Set each trials total_time_point
        setSize=int(cp.get('layerconfig','layer1SetSize'))
        timeScale=float(cp.get('layerconfig','layer1TimeScale'))
        TTP=int(sampleRate*setSize*timeScale)
        self.add_opt('total_time_points',str(TTP))
        #Set segment size for each TF map
        STP=powOf2Floor(timeScale*sampleRate)
        self.add_opt('segment_time_points',str(STP))
        #Set overlap size for each TF map
        overlap=int(math.floor(STP*overlapPercentage))
        self.add_opt('overlap',str(overlap))

        #Check for needed directories if not present make them!
        buildDir(blockPath)
        buildDir(layerPath)
        buildDir(blockPath+'/logs')

        #Set the condor submit initial directory option this
        #should point to the block,layer location for output files
        self.add_condor_cmd('initialdir',layerPath)
        #Set log
        channelName=string.strip(cp.get('tracksearchtime','channel_name'))
        self.set_stdout_file(os.path.normpath(blockPath+'/logs/tracksearchTime-'+channelName+'-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.normpath(blockPath+'/logs/tracksearchTime-'+channelName+'-$(cluster)-$(process).err'))
        self.set_sub_file(self.dagDirectory+'tracksearchTime.sub')
    #End init
#End Class

class tracksearchMapJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    The class running a MAP instance of the tracksearch code
    """
    def __init__(self,cp,block_id,layer_id,dagDir):
        self.dagDirectory=dagDir
        #ConfigParser object -> cp
        self.__executable = cp.get('condor','tracksearch')
        self.__universe = cp.get('condor','universe');
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp)
        blockID=block_id
        layerID=layer_id
        ### WE NEED A CONDOR DAG VARIABLE TO THE INITIAL DIR ARGUMENT
        layerPath=determineLayerPath(cp,blockID,layerID)
        blockPath=determineBlockPath(cp,blockID)
        #Parse all the options from BASE and Timeseries Config sections
        for sec in ['tracksearchbase']:
            self.add_ini_opts(cp,sec)
        for sec in ['tracksearchmap']:
            self.add_ini_opts(cp,sec)
        #Check for needed directories if not present make them!
        buildDir(blockPath)
        buildDir(layerPath)
        buildDir(blockPath+'/logs/')
        #Set the condor submit initial directory option this
        #should point to the block,layer location for output files
        self.add_condor_cmd('initialdir','$(macroStartDir)')
        #Set log
        self.set_stdout_file(os.path.normpath(blockPath+'/logs/tracksearchMap-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.normpath(blockPath+'/logs/tracksearchMap-$(cluster)-$(process).err'))
        self.set_sub_file(self.dagDirectory+'tracksearchMap.sub')
    #End init
#End Class

class tracksearchAveragerJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    The class that defines how we will perform an averager job.  The averager
    jobs, are tasks that merge two TF maps into one.  This is performed on a
    set of maps.  The set is then a new image of reduced dimesions specified
    via input arguments
    """
    def __init__(self,cp,block_id,layer_id,dagDir):
        self.dagDirectory=dagDir
    #Setup job options to take a cache of caches and build new maps
            #ConfigParser object -> cp
        self.__executable = cp.get('condor','averager')
        self.__universe = cp.get('condor','universe');
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp)
        blockID=block_id
        layerID=layer_id
        layerPath=determineLayerPath(cp,blockID,layerID)
        blockPath=determineBlockPath(cp,blockID)
        #Parse all the options from BASE and Timeseries Config sections
        for sec in ['averagerbase']:
            self.add_ini_opts(cp,sec)
        for sec in ['tracksearchmap']:
            self.add_ini_opts(cp,sec)
        #Check for needed directories if not present make them!
        buildDir(blockPath)
        buildDir(layerPath)
        buildDir(blockPath+'/logs/')
        #Set the condor submit initial directory option this
        #should point to the block,layer location for output files
        self.add_condor_cmd('initialdir','$(macroStartDir)')
        #Set log
        self.set_stdout_file(os.path.normpath(blockPath+'/logs/tracksearchAverager-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.normpath(blockPath+'/logs/tracksearchAverager-$(cluster)-$(process).err'))
        self.set_sub_file(self.dagDirectory+'tracksearchAverager.sub')
    #End init
#End Class

class tracksearchMapCacheBuildJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    This is the class definition which will drive the compiled python file name
    parsing code to create the layers of map building caches required
    to complete the pipeline
    """
    def __init__(self,cp,block_id,layer_id,dagDir):
        self.dagDirectory=dagDir
    #Setup job options to take a cache of caches and build new maps
            #ConfigParser object -> cp
        self.__executable = cp.get('condor','cachebuilder')
        self.__universe = cp.get('condor','universe');
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp)
        blockID=block_id
        layerID=layer_id
        #WE NEED A DAG VARIABLE THE LAYERid FOR PROPER INITIAL DIR
        layerPath=determineLayerPath(cp,blockID,layerID)
        blockPath=determineBlockPath(cp,blockID)
        #From ini sections [averagerbase] and [layerconfig] we will determine
        # the proper arguments to give the mapCacheBuildJob
        layerconfigOpts=[]
        section='layerconfig'
        for opt in cp.options(section):
            arg=string.strip(cp.get(section,opt))
            if str(opt).__contains__(str('TimeScale').lower()):
                layerconfigOpts.append(arg)
        if layerconfigOpts.__len__() <= 1:
            print 'Error with section Layerconfig!'
        if layerconfigOpts.__len__() < layerID:
            print 'Error invoking mapBuilderObject: layerID problem!',layerconfigOpts,layerID
        layerconfigOpts.sort()
        #Error condition layerID >= len(layerList)
        #Determine overlap conditions from [tracksearchbase] as percentage
        percentOverlap=float(string.strip(cp.get('layerconfig','layerOverlapPercent')))
        #Get the overlap of new cache files.
        layerSizeLabel='layer%sSetSize'%layerID
        layerTimeLabel='layer%sTimeScale'%layerID
        self.jobSetSize=cp.get('layerconfig',layerSizeLabel)
        self.expectedTotalDuration=cp.get('layerconfig','layerTopBlockSize')
        self.overlapTime=float(string.strip(cp.get('layerconfig',layerTimeLabel)))*percentOverlap
        self.mapTime=float(string.strip(cp.get('layerconfig',layerTimeLabel)))
        #Work on all files listed
        #OK THIS IS INCORRECT SETTING OF MACROS
        self.add_opt('file','$macroFile')
        self.add_opt('start_time','$macroStartTime')
        self.add_opt('map_set_duration','$macroMapSetDuration')
        self.add_opt('new_map_duration','$macroNewMapDuration')
        self.add_opt('overlap_maps','$macroOverlapMaps')
        #Parse all the options from BASE and Timeseries Config sections
        for sec in ['cachebuilder']:
            self.add_ini_opts(cp,sec)
        #Check for needed directories if not present make them!
        buildDir(blockPath)
        buildDir(layerPath)
        buildDir(blockPath+'/logs/')
        #Adjust initial dir so output ends up in proper place!
        self.add_condor_cmd('initialdir','$(macroStartDir)')
        #Add a condor macro variable for specialized value passing
        #Set log
        self.set_stdout_file(os.path.normpath(blockPath+'/logs/tracksearchCacheBuild-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.normpath(blockPath+'/logs/tracksearchCacheBuild-$(cluster)-$(process).err'))
        self.set_sub_file(self.dagDirectory+'tracksearchCacheBuild.sub')
    #End init

    def getJobsetList(self):
        #Generate the theoretical list of Jobset names to process based
        #on the configuration options given to the job
        outputList=[]
        A=float(self.expectedTotalDuration)
        B=float(self.mapTime)
        C=float(self.overlapTime)
        F=float(self.jobSetSize)
        jobCacheNumEstimate=int(int((A-C)/(B-C))/int(F))
        if (jobCacheNumEstimate == 0):
            print 'Error check object initialization.'
            return outList
        for num in range(1,jobCacheNumEstimate+1):
            outputList.append('JobSet_'+str(num)+'.cacheSet')
        return outputList
    #End getJobsetList

    def getJobTSAList(self):
        #Generate the theoretical list of Jobset names to process based
        #on the configuration options given to the job
        outputList=[]
        A=float(self.expectedTotalDuration)
        B=float(self.mapTime)
        C=float(self.overlapTime)
        F=float(self.jobSetSize)
        jobCacheNumEstimate=int(int((A-C)/(B-C))/int(F))
        if (jobCacheNumEstimate == 0):
            print 'Error check object initialization.'
            return outList
        for num in range(1,jobCacheNumEstimate+1):
            outputList.append('JobSet_'+str(num)+'.cacheTSA')
        return outputList
    #End getJobTSAList
#End Class

class tracksearchDataFindJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    THIS CODE IS NOT USED PLEASE IGNORE THIS CLASS IS DEAD!
    This class wraps up the calls LSCdataFindJob and sets relationships
    """
    def __init__(self,cp,block_id,dagDir):
        self.dagDirectory=dagDir
        blockID=block_id
        layerID=1
        blockPath=cp.get('fileLayout','workpath')+'/'+blockID+'/'+layerID+'/'
        cachePath=blockPath+'/'+'cache'
        logPath=blockPath
        self.df_job=pipeline.LSCDataFindJob(cachePath,logPath,cp)
        pad = 0
        length=int(cp.get('layerconfig','layer1TimeScale'))
        overlap=0
        self.df_job.set_sub_file(self.dagDirectory+'datafind.sub')
        self.df_job.add_condor_cmd('initialdir',determineLayerPath())
    #End init
#End Class
                   
class tracksearchTimeNode(pipeline.CondorDAGNode,
                          pipeline.AnalysisNode):
    """
    The class acting as a generic template to do a Time search
    """
    def __init__(self,job):
        #Expects to run CondorDAGJob object for tracksearch
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)
    #End init
#End Class

class tracksearchMapNode(pipeline.CondorDAGNode,
                         pipeline.AnalysisNode):
    """
    The class acting as a generic template to do a MAP search
    """
    def __init__(self,job):
        #Expects to run CondorDAGJob object for tracksearch
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)
    #End init
#End Class

class tracksearchAveragerNode(pipeline.CondorDAGNode,
                         pipeline.AnalysisNode):
    """
    The class acting as a generic template to do a MAP merging from specified
    input map cache file to process
    """
    def __init__(self,job):
        #Expects to run CondorDAGJob object for tracksearch
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)
    #End init
#End Class      

class tracksearchMapCacheBuildNode(pipeline.CondorDAGNode,
                         pipeline.AnalysisNode):
    """
    The class acting as a generic template to do a MAP cache build.
    It managaes a compiled python script to process a dir list of files.
    """
    def __init__(self,job):
        #Expects to run CondorDAGJob object for tracksearch
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)
    #End init
#End Class

class tracksearch:
    """
    The class which wraps all 100 second blocks into a uniq DAG
    we expect an input of a Science Segment of exactly n seconds
    specified via the ini file
    """
    def __init__(self,cparser,scienceSegment,logfile):
        #Expects a fixed size PsuedoScience Segment
        #The dag for this run will be then run from this data
        self.sciSeg=scienceSegment
        self.runStartTime=self.sciSeg.start()
        #cp -> ConfigParser object
        self.cp=cparser
        #The variable to link a previous jobblock into next block
        self.blockJobLinkList=[]
        self.blockID=str(self.sciSeg.start())+':'+str(self.sciSeg.dur())
        self.dag = pipeline.CondorDAG(logfile)
        self.dagName='tracksearchDAG_'+str(self.sciSeg.start())+'_duration_'+str(self.sciSeg.dur())
        self.resultPath=self.cp.get('fileLayout','workpath')
        self.dagDirectory=self.getDagDirectory()
        self.dagFilename=self.dagDirectory+'/'+self.dagName
        buildDir(self.dagDirectory)
        self.dag.set_dag_file(self.dagFilename)
        #Variables that are common to all search layers
        self.percentOverlap=float(string.strip(self.cp.get('layerconfig','layerOverlapPercent')))
        self.layerTimeScaleInfo=[]
        tempLib=[]
        for opt in self.cp.options('layerconfig'):
            section='layerconfig'
            arg=string.strip(self.cp.get(section,opt))
            if str(opt).__contains__(str('TimeScale').lower()):
                entry=[]
                entry.append(str(opt))
                entry.append(str(arg))
                tempLib.append(entry)
        tempLib.sort()
        self.layerTimeScaleInfo=tempLib
        tempLib=[]
        self.layerSetSizeInfo=[]
        for opt in self.cp.options('layerconfig'):
            section='layerconfig'
            arg=string.strip(self.cp.get(section,opt))
            if str(opt).__contains__(str('SetSize').lower()):
                entry=[]
                entry.append(str(arg))
                tempLib.append(entry)
        tempLib.sort()
        self.layerSetSizeInfo=tempLib
        if self.layerTimeScaleInfo.__len__()!=self.layerSetSizeInfo.__len__():
            print 'ERROR with section [layerconfig]'
        self.layerCount=self.layerTimeScaleInfo.__len__()

    #End Init

    def getDagDirectory(self):
        #The blockID is just the GPSstart:Duration
        dagDirectory=str(self.resultPath+'/DAG_files/'+self.blockID+'/')
        return dagDirectory
    #end def
    
    def startingSearchLayer(self,layerID):
        #This layer is responsible for initial data find jobs and the
        #tracksearch time jobs in a single layer search this would be
        #the only layer to actually be used
        #RETURNS the last node of this layer to satisify the
        #parent child relationship to next execution layer
        #layerID is expected to be a 1 Duh, first layer!
        print 'Preparing layer ',layerID,' of block ',self.blockID
        #We need to make chunks of data from this sciSeg which are
        #layer1SetSize long in time
        mapTime=float(self.cp.get('layerconfig','layer1TimeScale'))
        overlapTime=mapTime*self.percentOverlap
        setSize=float(self.cp.get('layerconfig','layer1SetSize'))
        ifo=str(self.cp.get('datafind','observatory'))
        #Create the chunk list that we will makeing data find jobs for
        #Give allowance for 1 extra map worth of data just in case
        self.sciSeg.make_chunks(setSize*mapTime+math.ceil(overlapTime),math.ceil(overlapTime))
        #Setup time marker for entire run
        runStartTime=self.runStartTime
        #What is name of dagDir location of dags and sub files
        dagDir=self.dagDirectory
        #Create dataFindJob
        #This path is relative to the jobs initialDir arguments!
        dataFindInitialDir=determineLayerPath(self.cp,self.blockID,layerID)
        outputDir=dataFindInitialDir
        dataFindLogPath=os.path.normpath(determineBlockPath(self.cp,self.blockID)+'/logs/')
        df_job = pipeline.LSCDataFindJob(outputDir,dataFindLogPath,self.cp)
        df_job.set_sub_file(self.dagDirectory+'datafind.sub')
        df_job.add_condor_cmd('initialdir',str(dataFindInitialDir))
        prev_df = None
        for chunk in self.sciSeg:
            df_node=pipeline.LSCDataFindNode(df_job)
            df_node.set_start(chunk.start())
            df_node.set_end(chunk.end())
            df_node.set_observatory(ifo[0])
            if prev_df:
                df_node.add_parent(prev_df)
            self.dag.add_node(df_node)
            prev_df = df_node
        #Setup the jobs after queries to LSCdataFind as equal children
        block_id=self.blockID
        layer_id=layerID
        tracksearchTime_job=tracksearchTimeJob(self.cp,block_id,layer_id,dagDir)
        prevLayerJobList=[]
        #THE FOLLOWIN LOOP NEEDS TO BE REVISED
        #WE ARE LOOSING THE CACHE NAMES IN THE PREVIOUS LOOP
        for chunk in self.sciSeg:
            tracksearchTime_node=tracksearchTimeNode(tracksearchTime_job)
            #Fix format of FILENAME
            #Inserted a SETB Pool Hack for testing
            tracksearchTime_node.add_var_opt('cachefile','/home/charlie/pipeTest/testFrame.cache')
            #tracksearchTime_node.add_var_opt('frame_cache',df_node.get_output())

            #Set the node job time markers from chunk!
            tracksearchTime_node.add_var_opt('gpsstart_seconds',chunk.start())

            #We want these nodes to be a child of last data find job
            tracksearchTime_node.add_parent(df_node)
            self.dag.add_node(tracksearchTime_node)
            prevLayerJobList.append(tracksearchTime_node)
        return prevLayerJobList
    #Finished with getting data and layer 1 of analysis submit files
    #end def

    def finalSearchLayer(self,layerID,nodeLinks):
        #This layer will setup last nodes of dag.  These nodes
        #should perform a final map search.  No additional cache
        #building needs to be done here
        #RETURNS no node linkage
        closeme=0
    #end def

    def clusterSearchLayer(self):
        #clustering and other "post jobs" go here
        closeme=0
    #end def

    def intermediateSearchLayers(self,layerID,nodeLinks):
        #This layer performs the required map cache building from
        #previous ith layer.  It is these caches which are joined into
        #new maps and placed as part of the i+1 layer and analyzed
        #RETURNS node linkage to possibly perform another search layer
        layerNum=layerID
        layerNumPrevious=layerID-1
        layer_id = layerNum 
        block_id=self.blockID
        dagDir=self.dagDirectory
        if (layerNumPrevious < 1):
            print 'Error calling the intermediate search layer method!'
        print 'Preparing layer ',layerID,' of block ',self.blockID
        prevLayerJobList=nodeLinks
        # Setup each additonal individual layer
        # The cache build node list clear after first pass through loop
        if (layerNum < 2):
            prevLayerJobList=[]
        cacheBuild_job=tracksearchMapCacheBuildJob(self.cp,block_id,layerNum,dagDir)
        cacheBuild_node=tracksearchMapCacheBuildNode(cacheBuild_job)
        #Add directory to process code expects file but behavior will
        #adjust to a directory listing if that is the case
        #Specify directory contain map files to setup
        #Should be previous layer already processed!
        cacheBuildMapDir=determineLayerPath(self.cp,block_id,layerNumPrevious)
        cacheBuildWorkDir=determineLayerPath(self.cp,block_id,layerNum)
        cacheBuild_node.add_macro('macroStartDir',cacheBuildWorkDir)
        cacheBuild_node.add_macro('macroFile',cacheBuildMapDir)            
        #Set the time of this run to start preping caches for
        cacheBuild_node.add_macro('macroStartTime',self.runStartTime)
        #Lookup the proper layer options
        layerMapNewDur=float(self.cp.get('layerconfig','layer'+str(layerNum)+'TimeScale'))
        layerMapSetDur=float(layerMapNewDur*float(self.cp.get('layerconfig','layer'+str(layerNum)+'SetSize')))
        layerMapOverlap=float(self.percentOverlap*layerMapNewDur)
        cacheBuild_node.add_macro('macroMapSetDuration',layerMapSetDur)
        cacheBuild_node.add_macro('macroNewMapDuration',layerMapNewDur)
        cacheBuild_node.add_macro('macroOverlapMaps',layerMapOverlap)
        #Make this process the child of all frame analysis jobs
        for parentJob in prevLayerJobList:
            cacheBuild_node.add_parent(parentJob)
        self.dag.add_node(cacheBuild_node)
        
        # The merge map
        tracksearchAverager_job=tracksearchAveragerJob(self.cp,block_id,layerNum,dagDir)

        jobSetList=cacheBuild_job.getJobsetList()
        jobTSAList=cacheBuild_job.getJobTSAList()
        #Var to store copies of these objects to get right parent relation
        averagerJobListing=[]
        #Loop over all theoretical jobsets for map making
        averagerWorkDir=determineLayerPath(self.cp,block_id,layerNum)
        for cacheSet in jobTSAList:
            tracksearchAverager_node=tracksearchAveragerNode(tracksearchAverager_job)
            tracksearchAverager_node.add_var_opt('multi_cache',cacheSet)
            tracksearchAverager_node.add_macro('macroStartDir',averagerWorkDir)
            tracksearchAverager_node.add_parent(cacheBuild_node)
            self.dag.add_node(tracksearchAverager_node)
            averagerJobListing.append(tracksearchAverager_node)
        #Run this new layer map analysis
        nextLayerJobList=[]
        #The results of this analysis should be stored in layerNum layer
        tracksearchMap_job=tracksearchMapJob(self.cp,block_id,layerNum,dagDir)
        #The entries in the JobSet file will be the input cache sets
        #for this layer of the Analysis
        for cacheSet in jobSetList:
            tracksearchMap_node=tracksearchMapNode(tracksearchMap_job)
            tracksearchMap_node.add_var_opt('inject_map_cache',cacheSet)
            tracksearchMap_node.add_macro('macroStartDir',determineLayerPath(self.cp,block_id,layer_id+1))
            for parent in averagerJobListing:
                tracksearchMap_node.add_parent(parent)
            self.dag.add_node(tracksearchMap_node)
            nextLayerJobList.append(tracksearchMap_node)
        prevLayerJobList=nextLayerJobList
        return nextLayerJobList
    #end def
    
    def createJobs(self):
        #This method will call all the needed jobs/node objects to create
        #a coherent job dag
        #CHECK TO MAKE SURE LAYERCONFIG IS VALID HAS ENTRIES
        #We need a consistency check between above and [tracksearchbase]
        layerID=1

        nodeLinkage=self.startingSearchLayer(layerID)

        layerCount=self.layerCount
        #set of numbers [i,j)
        for layerNum in range(2,layerCount+1):
            nodeLinkage=self.intermediateSearchLayers(layerNum,nodeLinkage)

        self.finalSearchLayer(layerCount,nodeLinkage)


        """
        We still need the pre and post script setups for each job
        then the correct config section in the ini file
        also need to the final track clustering code
        There is an implied fix to mapBuild.py code so Jobset uses
        GPS markers and not floats in the filenames!!
        """
    #End createJobs
    def writePipelineDAG(self):
        #This method lays out the dag files to a simple submission to condor
        self.dag.write_sub_files()
        self.dag.write_dag()
    #End writeTracksearchTopBlockPipe
#End Class

