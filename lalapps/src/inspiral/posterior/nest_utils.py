# DAG Class definitions for Odds Pipeline
# (C) 2011 John Veitch, Kiersten Ruisard, Kan Wang

import glue
from glue import pipeline
import os


class InspNestJob(pipeline.CondorDAGJob):
    """
    Class defining the Condor Job for inspnest to be run as part of a pipeline
    Input arguments:
        cp          - A ConfigParser object containing the inspnest section
        submitFile     - Path to store the submit file
        logdir         - A directory to hold the stderr, stdout files of inspnest
    """
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','inspnest')
        pipeline.CondorDAGJob.__init__(self,"standard",exe)
        # Set the options which are always used
        self.set_sub_file(submitFile)
        self.add_ini_opts(cp,'inspnest')
        self.__cp=cp
        self.set_stdout_file(os.path.join(logdir,'inspnest-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'inspnest-$(cluster)-$(process).err'))
    get_cp = lambda self: self.__cp

class LALInferenceJob(pipeline.CondorDAGJob):
    """
    Class defining the Condor Job for inspnest to be run as part of a pipeline
    Input arguments:
        cp          - A ConfigParser object containing the inspnest section
        submitFile     - Path to store the submit file
        logdir         - A directory to hold the stderr, stdout files of inspnest
    """
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','lalinference')
        pipeline.CondorDAGJob.__init__(self,"standard",exe)
        # Set the options which are always used
        self.set_sub_file(submitFile)
        self.add_ini_opts(cp,'lalinference')
        self.__cp=cp
        self.set_stdout_file(os.path.join(logdir,'lalinference-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'lalinference-$(cluster)-$(process).err'))
    get_cp = lambda self: self.__cp

class LALInferenceNode(pipeline.CondorDAGNode):
    """
    Class defining a condor DAG node for lalinferencenest jobs
    """
    def __init__(self,li_job):
        pipeline.CondorDAGNode.__init__(self,li_job)
    def set_seed(self,seed):
        self.add_var_opt('randomseed',seed)
    def set_dataseed(self,seed):
        self.add_var_opt('dataseed',seed)
    def add_ifo_data(self,data_tuples,ifos=None,shift_time_dict=None):
        """
        Add list of IFOs and data to analyse.
        data_tuples is a dictionary of available (cache,channel) tuples
        ifos is an optional list of IFOs to include. If not specified analyse all available
        """
        if shift_time_dict is not None:
            raise Exception('ERROR: Time slides not supported by lalinference yet')
        cp = self.job().get_cp()
        allifos=data_tuples.keys()
        if ifos is None:
            self.__ifos=allifos
        else:
            self.__ifos=ifos

        ifostring='['
        cachestring='['
        channelstring='['
        first=True
        for ifo in self.__ifos:
            if first:
                delim=''
                first=False
            else: delim=','
            if data_tuples[ifo][1] is None:
                cache=self.job().get_cp().get('data',ifo.lower()+'-channel')
            else:
                cache=data_tuples[ifo][1].get_df_node().get_output_files()[0]
                self.add_parent(data_tuples[ifo][1].get_df_node())
            ifostring=ifostring+delim+ifo
            cachestring=cachestring+delim+cache
            channelstring=channelstring+delim+self.job().get_cp().get('data',ifo.lower()+'-channel')
        ifostring=ifostring+']'
        cachestring=cachestring+']'
        channelstring=channelstring+']'
        self.add_var_arg('--IFO '+ifostring)
        self.add_var_arg('--channel '+channelstring)
        self.add_var_arg('--cache '+cachestring)
        # Start at earliest common time
        # NOTE: We perform this arithmetic for all ifos to ensure that a common data set is
        # Used when we are running the coherence test.
        # Otherwise the noise evidence will differ.
        starttime=max([int(data_tuples[ifo][0][0]) for ifo in allifos])
        endtime=min([int(data_tuples[ifo][0][1]) for ifo in allifos])
        self.__GPSstart=starttime
        self.__GPSend=endtime
        length=endtime-starttime
    
        # Now we need to adjust the start time and length to make sure the maximum data length
        # is not exceeded.
        trig_time=self.get_trig_time()
        maxLength=float(cp.get('analysis','analysis-chunk-length'))
        if(length > maxLength):
            while(self.__GPSstart+maxLength<trig_time and self.__GPSstart+maxLength<self.__GPSend):
                    self.__GPSstart+=maxLength/2.0
        # Override calculated start time if requested by user in ini file
        if self.job().get_cp().has_option('lalinference','psdstart'):
		self.__GPSstart=self.job().get_cp().getfloat('lalinference','psdstart')
		print 'Over-riding start time to user-specified value %f'%(self.__GPSstart)
		if self.__GPSstart<starttime or self.__GPSstart>endtime:
			print 'ERROR: Over-ridden time lies outside of science segment!'
			raise Exception('Bad psdstart specified')
	else: 
		self.add_var_opt('psdstart',str(self.__GPSstart))
        if self.job().get_cp().has_option('lalinference','psdlength'):
		length=self.job().get_cp().getfloat('lalinference','psdlength')
		print 'Over-riding PSD length to user-specified value %f'%(length)
	else:
		length=self.__GPSend-self.__GPSstart
	        if(length>maxLength):
        	    length=maxLength
	        self.add_var_opt('PSDlength',str(int(length)))
        self.add_var_opt('seglen',self.job().get_cp().get('analysis','psd-chunk-length'))


    def get_ifos(self):
        return ''.join(map(str,self.__ifos))

    def set_output(self,outfile):
        """
        Set output file command line argument and add outfile
        to the list of nodes output files
        """
        self.add_file_opt('outfile',outfile,file_is_output_file=True)

    def set_trig_time(self,time):
        """
        Set the end time of the signal for the centre of the prior in time
        """
        self.__trigtime=float(time)
        self.add_var_opt('trigtime',str(time))

    def set_event_number(self,event):

        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_opt('event',str(event))

    get_trig_time = lambda self: self.__trigtime

class InspNestNode(pipeline.CondorDAGNode):
    """
    Class defining a Condor DAG node for inspnest jobs.
    Input arguments:
    inspnest_job - A InspNestJob object
    """
    def __init__(self,inspnest_job):
        pipeline.CondorDAGNode.__init__(self,inspnest_job)
    def finalize(self):
        pipeline.CondorDAGNode.finalize(self)
        cp=self.job().get_cp()
        maxLength=float(cp.get('analysis','analysis-chunk-length'))
        trigtime=self.__trigtime
        if(self.__GPSend-self.__GPSstart > maxLength):
            while(self.__GPSend-self.__GPSstart >maxLength):
                if self.__GPSend-trigtime > trigtime-self.__GPSstart:
                    self.__GPSend-=1
                else:
                    self.__GPSstart+=1
        length=self.__GPSend-self.__GPSstart
        if self.__trigtime<self.__GPSstart or self.__trigtime>self.__GPSend:
            print 'ERROR in generating Inspnest Node, trig time %f not inside data segment (%f,%f)'%(self.__trigtime,self.__GPSstart,self.__GPSend)
            raise Exception('Bad data segment or trig time')
        self.add_var_opt('GPSstart',str(self.__GPSstart))
        self.add_var_opt('length',str(int(length)))
        self.add_var_opt('Nsegs',str(int(length/float(self.job().get_cp().get('analysis','psd-chunk-length')))))
    def set_seed(self,seed):
        self.add_var_opt('seed',seed)
    def set_dataseed(self,seed):
        self.add_var_opt('dataseed',seed)
    def add_ifo_data(self,data_tuples,ifos=None,shift_time_dict=None):
        """
        Add list of IFOs and data to analyse.
        data_tuples is a dictionary of available (cache,channel) tuples
        ifos is an optional list of IFOs to include. If not specified analyse all available
        """
        cp = self.job().get_cp()
        allifos=data_tuples.keys()
        if ifos is None:
            self.__ifos=allifos
        else:
            self.__ifos=ifos

        for ifo in self.__ifos:
            self.add_var_arg('--IFO '+ifo)
            self.add_var_arg('--channel '+self.job().get_cp().get('data',ifo.lower()+'-channel'))
            if data_tuples[ifo][1] is None:
                self.add_var_arg('--cache '+self.job().get_cp().get('data',ifo.lower()+'-channel'))
            else:
                self.add_var_arg('--cache '+data_tuples[ifo][1].get_df_node().get_output_files()[0])
                self.add_parent(data_tuples[ifo][1].get_df_node())
        # Start at earliest common time
        # NOTE: We perform this arithmetic for all ifos to ensure that a common data set is
        # Used when we are running the coherence test.
        # Otherwise the noise evidence will differ.
	offsets={}
	for k in allifos:
        	if shift_time_dict is not None and k[0] in allifos:
			offsets[k]=shift_time_dict[k[0]]
		else:
			offsets[k]=0.0
	starttime=max([int(data_tuples[ifo][0][0])+offsets[ifo] for ifo in allifos])
        endtime=min([int(data_tuples[ifo][0][1])+offsets[ifo] for ifo in allifos])
        self.__GPSstart=starttime
        self.__GPSend=endtime
        length=endtime-starttime
    
        if shift_time_dict:
        	# Make sure that all trig times are in our time interval.
        	# shift_times is a list of all time shifts.
        	# shift_times = shift_time_dict.values()
        	# mid_shift_time is the middle of shift times.
        	#max_trig_time=max(shift_times)+trig_time
        	#min_trig_time=min(shift_times)+trig_time
        	#length=min(length,maxLength)
        	#self.__GPSstart=max(self.__GPSstart,(min_trig_time+max_trig_time-length)/2) 
        	#self.add_var_opt('GPSstart',str(self.__GPSstart))
        	for ifo in shift_time_dict:
        		self.add_var_arg('--'+ifo+'GPSshift '+str(shift_time_dict[ifo]))

    def get_ifos(self):
        return ''.join(map(str,self.__ifos))

    def set_output(self,outfile):
        """
        Set output file command line argument and add outfile
        to the list of nodes output files
        """
        self.add_file_opt('out',outfile,file_is_output_file=True)

    def set_trig_time(self,time):
        """
        Set the end time of the signal for the centre of the prior in time
        """
        self.__trigtime=float(time)
        self.add_var_opt('end_time',str(time))

    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_opt('event',str(event))

    get_trig_time = lambda self: self.__trigtime
    

class CoherenceTestJob(pipeline.CondorDAGJob):
    """
    Class defining the coherence test job to be run as part of a pipeline.
    """
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','coherencetest')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        self.add_opt('coherent-incoherent-noise','')
        self.add_condor_cmd('getenv','True')
        self.set_stdout_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'coherencetest-$(cluster)-$(process).err'))
        self.set_sub_file(submitFile)

class CoherenceTestNode(pipeline.CondorDAGNode):
    """
    Class defining the node for the coherence test
    """
    def __init__(self,coherencetest_job):
        pipeline.CondorDAGNode.__init__(self,coherencetest_job)

class CombineZJob(pipeline.CondorDAGJob):
    """
    Class defining a combineZ script job to be run as part of a pipeline
    This job combine runs with adjacent prior areas and produces the posterior samples
    Input Arguments:
        cp        - A ConfigParser object containing the combinez section
        submitFile    - Path to store the submit file
        logdir        - A directory to hold the stderr, stdout files of combineZ
    """
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','combinez')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        self.add_opt('Nlive',str(int(cp.get('analysis','nlive'))*int(cp.get('analysis','nparallel'))))
        self.set_stdout_file(os.path.join(logdir,'combineZ-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'combineZ-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        self.set_sub_file(submitFile)

class CombineZNode(pipeline.CondorDAGNode):
    """
    Class defining a Condor DAG Node for combineZ jobs
    Input Arguments:
        combine_job - A CombineZJob object
    """
    def __init__(self,combine_job):
        pipeline.CondorDAGNode.__init__(self,combine_job)

class MergeJob(pipeline.CondorDAGJob):
    """
    Class defining a job which merges several parallel nested sampling jobs into a single file
    Input arguments:
        cp        - A configparser object containing the setup of the analysis
        submitFile    - Path to store the submit file
        logdor        - A directory to hold the stderr, stdout files of the merge runs
    """
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','mergescript')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'merge-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'merge-$(cluster)-$(process).err'))
        self.add_opt('Nlive',cp.get('analysis','nlive'))
        self.add_condor_cmd('getenv','True')

class MergeNode(pipeline.CondorDAGNode):
    """
    Class defining the DAG node for a merge job
    Input arguments:
        merge_job = A MergeJob object
    """
    def __init__(self,merge_job):
        pipeline.CondorDAGNode.__init__(self,merge_job)

#

class ResultsPageJob(pipeline.CondorDAGJob):
    def __init__(self,cp,submitFile,logdir):
        exe=cp.get('condor','resultspage')
        pipeline.CondorDAGJob.__init__(self,"vanilla",exe)
        self.set_sub_file(submitFile)
        self.set_stdout_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).out'))
        self.set_stderr_file(os.path.join(logdir,'resultspage-$(cluster)-$(process).err'))
        self.add_condor_cmd('getenv','True')
        # self.add_opt('Nlive',cp.get('analysis','nlive'))
        
        if cp.has_option('results','skyres'):
            self.add_opt('skyres',cp.get('results','skyres'))


class ResultsPageNode(pipeline.CondorDAGNode):
    def __init__(self,results_page_job):
        pipeline.CondorDAGNode.__init__(self,results_page_job)
    def set_event_number(self,event):
        """
        Set the event number in the injection XML.
        """
        if event is not None:
            self.__event=int(event)
            self.add_var_arg('--eventnum '+str(event))
# Function definitions for setting up groups of nodes

def setup_single_nest(cp,nest_job,end_time,data,path,ifos=None,event=None,nodeclass=InspNestNode,timeslides=None):
    """
    Setup nodes for analysing a single time
    cp - configparser object
    NestJob - Job for the inspnest nodes
    time - gps time of the centre of the prior
    data - dictionary of (seg,science_segment) tuple for ifos
    ifos - optional list of IFOs to analyse. If not specified analyse all available.
    nodeclass - custom class to use, must implement the InspNestNode above
    """
    nest_node=nodeclass(nest_job)
    nest_node.set_trig_time(end_time)
    nest_node.set_event_number(event)
    if cp.has_option('analysis','seed'):
        ini_chain_seed=int(cp.get('analysis','seed'))
    else:
        ini_chain_seed=100
    if cp.has_option('analysis','dataseed') and event is not None:
        ini_noise_seed=int(cp.get('analysis','dataseed'))
    nest_node.set_seed(str(ini_chain_seed))
    nest_node.set_dataseed(str(event+ini_noise_seed))
    nest_node.add_ifo_data(data,ifos,shift_time_dict=timeslides)
    outfile_name=os.path.join(path,'outfile_%f_%s.dat'%(end_time,nest_node.get_ifos()))
    nest_node.set_output(outfile_name)
    return nest_node


def setup_parallel_nest(cp,nest_job,merge_job,end_time,data,path,ifos=None,event=None,nodeclass=InspNestNode,timeslides=None):
    """
    Setup nodes for analysing a single time using
    parallel runs
    cp - configparser objection
    nest_job - job for the inspnest nodes
    merge_job - job for the merge nodes
    end_time - time to centre prior on
    data - dictionary of (seg,science_segment) tuple for ifos
    ifos - optional list of IFOs to analyse. If not specified analyse all available.
    nodeclass - custom class to use, must implement the InspNestNode above
    timeslides - dictionary of {ifo: shift} to specify time slides
    """
    nparallel=int(cp.get('analysis','nparallel'))
    merge_node=MergeNode(merge_job)
    merge_node.add_var_opt('Nlive',cp.get('analysis','nlive'))
    nest_nodes=[]
    if cp.has_option('analysis','seed'):
        ini_chain_seed=int(cp.get('analysis','seed'))
    else:
        ini_chain_seed=100
    if cp.has_option('analysis','dataseed') and event is not None:
        noise_seed=str(int(event)+int(cp.get('analysis','dataseed')))
    else: 
        noise_seed=0
    for i in range(nparallel):
        nest_node=nodeclass(nest_job)
        nest_node.set_trig_time(end_time)
        nest_nodes.append(nest_node)
        nest_node.add_ifo_data(data,ifos,shift_time_dict=timeslides)
        nest_node.set_event_number(event)
        p_outfile_name=os.path.join(path,'outfile_%f_%i_%s.dat'%(end_time,i,nest_node.get_ifos()))
        nest_node.set_seed(str(i+ini_chain_seed))
        nest_node.set_dataseed(str(noise_seed))
        merge_node.add_parent(nest_node)
        merge_node.add_file_arg(p_outfile_name)
        nest_node.set_output(p_outfile_name)
    outfile_name=os.path.join(path,'outfile_%f_%s.dat'%(end_time,nest_node.get_ifos()))
    merge_node.add_file_opt('out',outfile_name,file_is_output_file=True)
    merge_node.add_file_opt('headers',p_outfile_name+'_params.txt',file_is_output_file=False)
    return (merge_node,nest_nodes)
