"""
Classes needed for the binary_pulsar analysis pipeline.

This script produces the necessary condor submit and dag files to run
the binary pulsar upper limit analysis.
"""


__author__ = 'Chris Messenger <cm@star.sr.bham.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

import string
import exceptions
import os
from glue import pipeline


class binarypulsarError(exceptions.Exception):
  def __init__(self, args=None):
    self.args = args

#############################################################################    
#############################################################################
class sensitivityJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_CalculateSensitivity job used by the binary_pulsar pipeline.
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini
  file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','sensitivity')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    Static options are read from the sensitivity, ephemeris, and source sections
    """
    for sec in ['sensitivity']:
      self.add_ini_opts(cp,sec)
    
    for sec in ['ephemeris']:
      self.add_ini_opts(cp,sec)

    for sec in ['source']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/sensitivity-$(cluster)-$(process).out')
    self.set_stderr_file('logs/sensitivity-$(cluster)-$(process).err')
    self.set_sub_file('sensitivity.sub')


class sensitivityNode(pipeline.CondorDAGNode,pipeline.AnalysisNode):
  """
  An sensitivityNode runs an instance of lalapps_CalculateSensitivity in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A CondorDAGJob that can run an instance of lalapps_CalculateSensitivity.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__datadir = None
    self.__detector = None
    self.__tstart = None
    self.__tend = None
    self.__outdir = None

  def set_datadir(self, datadir):
    """
    Set the data directory.
    """
    self.add_var_opt('datadir',datadir)
    self.__datadir = datadir

  def set_detector(self, detector):
    """
    Set the code to use as detector.
    detector = code (e.g. LLO, LHO or GEO).
    """
    self.add_var_opt('detector', detector)
    self.__detector = detector

  def set_tstart(self, tstart):
    """
    Set the GPS start time of this analysis.
    """
    self.add_var_opt('tstart',tstart)
    self.__tstart = tstart

  def set_tend(self, tend):
    """
    Set the GPS end time of the sensitivity analysis.
    """
    self.add_var_opt('tend', tend)
    self.__tend = tend

  def set_outdir(self, outdir):
    """
    Set the output directory of the sensitivity analysis.
    """
    self.add_var_opt('outdir', outdir)
    self.__outdir = outdir

#############################################################################    
#############################################################################  
class makemeshJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_GenerateBinaryMesh job used by the binary_pulsar pipeline. 
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','makemesh')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    Static options are read from the mesh, ephemeris, and source sections
    """
    for sec in ['mesh']:
      self.add_ini_opts(cp,sec)
    
    for sec in ['ephemeris']:
      self.add_ini_opts(cp,sec)

    for sec in ['source']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/makemesh-$(cluster)-$(process).out')
    self.set_stderr_file('logs/makemesh-$(cluster)-$(process).err')
    self.set_sub_file('makemesh.sub')


class makemeshNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An makemeshNode runs an instance of lalapps_GenerateBinaryMesh in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A makemeshNode that can run an instance of lalapps_GenerateBinaryMesh.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__datadir = None
    self.__mismatch = None
    self.__outdir = None

  def set_datadir(self,datadir):
    """
    Set the data directory.
    """
    self.add_var_opt('datadir',datadir)
    self.__datadir = datadir

  def set_outdir(self,outdir):
    """
    Set the output mesh directory.
    """
    self.add_var_opt('outdir',outdir)
    self.__outdir = outdir

#############################################################################    
#############################################################################
class searchJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_ComputeFStatisticBinary job used by the binary_pulsar pipeline. 
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','search')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    Static options are read from the search, ephemeris, and source sections
    """
    for sec in ['search']:
      self.add_ini_opts(cp,sec)
    
    for sec in ['ephemeris']:
      self.add_ini_opts(cp,sec)

    for sec in ['source']:
      self.add_ini_opts(cp,sec)

    self.set_stdout_file('logs/search-$(cluster)-$(process).out')
    self.set_stderr_file('logs/search-$(cluster)-$(process).err')
    self.set_sub_file('search.sub')


class searchNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An searchNode runs an instance of lalapps_ComputeFStatisticBinary in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A searchNode that can run an instance of lalapps_ComputeFStatisticBinary.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__f_min = None
    self.__f_band = None
    self.__datadir = None
    self.__detector = None
    self.__binarytemplatebank = None
    self.__workdir = None

  def set_datadir(self,datadir):
    """
    Set the data directory.
    """
    self.add_var_opt('datadir',datadir)
    self.__datadir = datadir

  def set_f_min(self,f_min):
    """
    Set the minimum search frequency.
    """
    self.add_var_opt('f_min',f_min)
    self.__f_min = f_min

  def set_f_band(self,f_band):
    """
    Set the search frequency band.
    """
    self.add_var_opt('f_band',f_band)
    self.__f_band = f_band

def set_detector(self,detector):
    """
    Set the detector code as described earlier.
    """
    self.add_var_opt('detector',detector)
    self.__detector = detector

def set_binarytemplatebank(self,binarytemplatebank):
    """
    Set the location of the binary template bank for the search.
    """
    self.add_var_opt('binarytemplatebank',binarytemplatebank)
    self.__binarytemplatebank = binarytemplatebank

def set_workdir(self,workdir):
    """
    Set the search output directory.
    """
    self.add_var_opt('workdir',workdir)
    self.__workdir = workdir

#############################################################################    
#############################################################################
class coincidenceJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_FindCoincidence job used by the binary_pulsar pipeline. 
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','coincidence')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    No Static options for the coincidence analysis
    """
       
    self.set_stdout_file('logs/coincidence-$(cluster)-$(process).out')
    self.set_stderr_file('logs/coincidence-$(cluster)-$(process).err')
    self.set_sub_file('coincidence.sub')


class coincidenceNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An coincidenceNode runs an instance of lalapps_FindCoincidence in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A coincidenceNode that can run an instance of lalapps_FindCoincidence.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__presults = None
    self.__sresults = None
    self.__f_min = None
    self.__f_max = None
    self.__searchfreqmeshfile = None
    self.__sdataparamsfile = None
    self.__outdir = None

  def set_presults(self,presults):
    """
    Set the primary search results directory.
    """
    self.add_var_opt('presults',presults)
    self.__presults = presults

  def set_sresults(self,sresults):
    """
    Set the secondary search results directory.
    """
    self.add_var_opt('sresults',sresults)
    self.__sresults = sresults

  def set_f_min(self,f_min):
    """
    Set the minimum primary results frequency.
    """
    self.add_var_opt('f_min',f_min)
    self.__f_min = f_min

  def set_f_max(self,f_max):
    """
    Set the maximum primary results frequency.
    """
    self.add_var_opt('f_max',f_max)
    self.__f_max = f_max

  def set_searchfreqmeshfile(self,searchfreqmeshfile):
    """
    Set the name of the file containing search frequency - mesh information.
    """
    self.add_var_opt('searchfreqmeshfile',searchfreqmeshfile)
    self.__searchfreqmeshfile = searchfreqmeshfile

  def set_sdataparamsfile(self,sdataparamsfile):
    """
    Set the name of the file containing the secondary data set parameters.
    """
    self.add_var_opt('sdataparamsfile',sdataparamsfile)
    self.__sdataparamsfile = sdataparamsfile

  def set_outdir(self,outdir):
    """
    Set the output coincidence directory.
    """
    self.add_var_opt('outdir',outdir)
    self.__outdir = outdir

#############################################################################    
#############################################################################
class injectionsJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_Injections job used by the binary_pulsar pipeline. 
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','injections')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    No Static options for the injections analysis
    """
       
    self.set_stdout_file('logs/injections-$(cluster)-$(process).out')
    self.set_stderr_file('logs/injections-$(cluster)-$(process).err')
    self.set_sub_file('injections.sub')


class injectionsNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An injectionsNode runs an instance of lalapps_Injections in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A injectionsNode that can run an instance of lalapps_Injections.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__configfile = None
    self.__f_min = None
    self.__f_max = None
    self.__pdatadir = None
    self.__sdatadir = None
    self.__gw_amplitude = None
    self.__ntrials = None

  def set_configfile(self,configfile):
    """
    Set the input config file for the injections.
    """
    self.add_var_opt('configfile',configfile)
    self.__configfile = configfile

  def set_f_min(self,f_min):
    """
    Set the minimum injection frequency
    """
    self.add_var_opt('f_min',f_min)
    self.__f_min = f_min

  def set_f_max(self,f_max):
    """
    Set the maximum injection frequency
    """
    self.add_var_opt('f_max',f_max)
    self.__f_max = f_max

  def set_pdatadir(self,pdatadir):
    """
    Set the directory containing the primary search data.
    """
    self.add_var_opt('pdatadir',pdatadir)
    self.__pdatadir = pdatadir

  def set_sdatadir(self,sdatadir):
    """
    Set the directory containing the secondary search data.
    """
    self.add_var_opt('sdatadir',sdatadir)
    self.__sdatadir = sdatadir

  def set_gw_amplitude(self,gw_amplitude):
    """
    Set the injection signal amplitude.
    """
    self.add_var_opt('gw_amplitude',gw_amplitude)
    self.__gw_amplitude = gw_amplitude

  def set_ntrials(self,ntrials):
    """
    Set the number of signal injections.
    """
    self.add_var_opt('ntrials',ntrials)
    self.__ntrials = ntrials

#############################################################################    
#############################################################################    
class upperlimitJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  A lalapps_ComputeUL job used by the binary_pulsar pipeline. 
  The stdout and stderr from the job are directed to the logs directory.
  The path to the executable and the universe is determined from the ini file.
  """
  def __init__(self,cp):
    """
    cp = ConfigParser object from which options are read.
    """
    self.__executable = cp.get('condor','upperlimit')
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp)

    """
    Static options are read from the upperlimit section
    """
    for sec in ['upperlimit']:
      self.add_ini_opts(cp,sec)
    
    self.set_stdout_file('logs/upperlimit-$(cluster)-$(process).out')
    self.set_stderr_file('logs/upperlimit-$(cluster)-$(process).err')
    self.set_sub_file('upperlimit.sub')


class upperlimitNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  An upperlimitNode runs an instance of lalapps_ComputeUL in a Condor DAG.
  """
  def __init__(self,job):
    """
    job = A upperlimitNode that can run an instance of lalapps_ComputeUL.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

    """
    initialise the job variables
    """
    self.__injectionsdir = None
    self.__coresultsdir = None
    self.__injfreqmeshfile = None
    self.__injh0file = None
    self.__outfile = None

  def set_injectionsdir(self,injectionsdir):
    """
    Set the location of the directory containing the injection results.
    """
    self.add_var_opt('injectionsdir',injectionsdir)
    self.__injectionsdir = injectionsdir

  def set_coresultsdir(self,coresultsdir):
    """
    Set the location of the directory containing the coincidence results.
    """
    self.add_var_opt('coresultsdir',coresultsdir)
    self.__coresultsdir = coresultsdir

  def set_injfreqmeshfile(self,injfreqmeshfile):
    """
    Set the location of the file containing the injection frequencies.
    """
    self.add_var_opt('injfreqmeshfile',injfreqmeshfile)
    self.__injfreqmeshfile = injfreqmeshfile

  def set_injh0file(self,injh0file):
    """
    Set the location of the file containing the gw injection amplitudes.
    """
    self.add_var_opt('injh0file',injh0file)
    self.__injh0file = injh0file

  def set_outfile(self,outfile):
    """
    Set the ouotput file. 
    """
    self.add_var_opt('outfile',outfile)
    self.__outfile = outfile
    
