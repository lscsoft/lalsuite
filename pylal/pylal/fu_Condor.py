"""
This module contains condor jobs / node classes for the followup dag

This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
import math
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict
sys.path.append('@PYTHONLIBDIR@')

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from glue import segments
from glue import segmentsUtils
from pylal.webUtils import *
from pylal.webCondor import *
from lalapps import inspiral
from pylal import fu_utils
from glue.ligolw import lsctables

###### WRAPPER FOR CONDOR DAG - TO MAKE THE FOLLOWUP DAG WEBIFIABLE ###########
###############################################################################

class followUpDAG(pipeline.CondorDAG, webTheDAG):

  def __init__(self, config_file, log_path):
    self.basename = re.sub(r'\.ini',r'', config_file) 
    tempfile.tempdir = log_path
    tempfile.template = self.basename + '.dag.log.'
    logfile = tempfile.mktemp()
    fh = open( logfile, "w" )
    fh.close()
    pipeline.CondorDAG.__init__(self,logfile)
    self.set_dag_file(self.basename)
    self.jobsDict = {}
    # The list remote_nodes will contain the list of nodes run remotely 
    # (such as V1 qscans)
    self.remote_nodes = []

########################################################
#### Methods common to several followup classes ########
########################################################

def checkHipeCachePath(cp):
  try:
    if len(string.strip(cp.get('followup-hipe-cache','hipe-cache-path'))) > 0:
      hipeCachePath = string.strip(cp.get('followup-hipe-cache','hipe-cache-path'))
    else:
      hipeCachePath = None
    return(hipeCachePath)
  except:
    print >> sys.stderr, "ERROR: failure in checkHipeCachePath()"
    return None

###############################################################################
### Following two methods are to check CP objects
### If CP object missing section a temporary default valued
### config parser (CP) object is generated
###############################################################################

def verifyCP(cp,defaults):
  """
  This method takes in a cp object and give a set of defaults check to
  make sure the section in question exists and at least the options in
  the default are specified with some value.  It return TRUE or FALSE
  depending if the cp object contains the sections and options
  specified by the input DEFAULTS.
  """
  return cp.has_section(defaults["section"]) and \
  all(cp.has_option(defaults["section"], opt) for opt in defaults["options"])
  # End verifyCP

def modifyCP(cp,defaults):
  """
  Appended the configuration information in defaults into the config
  parser (cp) object and return a copy of this newly update cp object.
  """
  if not(cp.has_section(defaults["section"])):
    cp.add_section(defaults["section"])
  for key, val in defaults["options"].iteritems():
    if not cp.has_option(defaults["section"], key):
      cp.set(defaults["section"], val)
  #End modifyCP


###############################################################################
#### A CLASS TO DO FOLLOWUP INSPIRAL JOBS ####################################
###############################################################################
class followUpInspJob(inspiral.InspiralJob,webTheJob):
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "inspiral_head":"lalapps_inspiral"
      }
    }
  def __init__(self,cp,type='plot'):
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    inspiral.InspiralJob.__init__(self,cp)
    if type == 'head':
      self.set_executable(string.strip(cp.get('condor','inspiral_head')))
    self.name = 'followUpInspJob' + type
    self.setupJobWeb(self.name)


class followUpInspNode(inspiral.InspiralNode,webTheNode):
  
  def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag, datafindCache, d_node, datafindCommand, type='plot', sngl_table = None):
    self.sample_rate = string.strip(cp.get('coh-inspiral','sample-rate'))
    if 1:#try:
      self.output_file_name = ""
      #inspiral.InspiralNode.__init__(self, inspJob) 
      # the use of this class would require some reorganisation in fu_Condor.py
      # and webCondor.py in order to set up the jobs following the same scheme
      # as the way it is done for the Inspiral pipeline...
      pipeline.CondorDAGNode.__init__(self,inspJob)
      injFile = self.checkInjections(cp)      
      hipeCache = checkHipeCachePath(cp)

      if type == "plot" or type == "notrig" or type == "coh" or type == "chia":
        # Here we define the trig-start-time and the trig-end-time;
        # The difference between these two times should be kept to 2s
        # Otherwise change the clustering window also
        hLengthAnalyzed = 1
	if type == "coh" or type == "chia": hLengthAnalyzed = 1.0
        self.set_trig_start( int(trig.gpsTime[ifo] - hLengthAnalyzed + 0.5) )
        self.set_trig_end( int(trig.gpsTime[ifo] + hLengthAnalyzed + 0.5) )

      if type == "plot" or type == "notrig" or type == "coh":
        self.add_var_opt("write-snrsq","")
        self.add_var_opt("write-chisq","")
        self.add_var_opt("write-spectrum","")
        self.add_var_opt("write-template","")
      if type == "chia" or type == "notrig" or type == "coh":
        self.add_var_opt("write-cdata","")

      if injFile: 
        self.set_injections( injFile )

      skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']
      if not hipeCache:
        skipParams.append('frame-cache')
        self.add_var_opt('frame-cache',datafindCache)        

      # initialize the extension of the output file. If the option 
      # write_compress is found in procParams the extension will be overwritten
      # later as .xml.gz
      extension = ".xml"
      for row in procParams:
        param = row.param.strip("-")
        value = row.value
	# override the options for coherent jobs (useful for skymaps so that
	# we can bump up the sample rate)
	if type == "coh" and cp.has_option("coh-inspiral",param):
	  value = cp.get("coh-inspiral",param)
        if type == "chia" and cp.has_option("coh-inspiral",param):
          value = cp.get("coh-inspiral",param)
        if param == 'bank-file':
          bankFile = value
        if type == "notrig" or type == "coh" or type == "chia":
        # if forceTrigger is true, we loose the thresholds to
        # make sure to get a trigger
          if param == 'snr-threshold': value = "0.1"
          # rsq veto must be disabled
          if param == 'do-rsq-veto': continue
          if param == 'enable-rsq-veto': continue
          # chisq veto is disabled by loosing its threshold 
          # we still want to generate the chisq time-series
          if param == 'chisq-threshold': value = "1.0e+06"
          # Using a window of 1s for clustering will allow us to always get
          # at least one trigger
          if param == 'cluster-method': value = 'window'
          if param == 'cluster-window': continue
          pass
        if param in skipParams: continue
        self.add_var_opt(param,value)
        # The attributes _AnalysisNode__end, _AnalysisNode__start,
        # _InspiralAnalysisNode__pad_data need to be defined before calling the
        # method "writeAll" of "class webTheDAG". This method calls
        # "write_sub_files()" in pipeline.py, which itself relies on the
        # "finalize()" method of "class InspiralAnalysisNode" in inspiral.py .
        # This is where all these attributes are being used. This hack is 
        # required because "inspiral.InspiralNode.__init__(self, inspJob)"
        # currently does not work within "class followUpInspNoDE"
        if param == 'gps-end-time':
          self.__end = value
          self._AnalysisNode__end = int(value)
        if param == 'gps-start-time':
          self.__start = value
          self._AnalysisNode__start = int(value)
        if param == 'pad-data': 
          self._InspiralAnalysisNode__pad_data = int(value)
        if param == 'ifo-tag':
          self.__ifotag = value
        if param == 'channel-name': self.inputIfo = value[0:2]
        if param == 'write-compress':
          extension = '.xml.gz'

      if type == "notrig" or type == "coh" or type == "chia":
        self.add_var_opt('cluster-window',str(hLengthAnalyzed/2.))
        self.add_var_opt('disable-rsq-veto',' ')

      # add the arguments that have been specified in the section 
      # [inspiral-extra] of the ini file (intended for 12-18 month analysis)
      if cp.has_section("followup-inspiral-extra"):
        for (name,value) in cp.items("followup-inspiral-extra"):
          self.add_var_opt(name,value)

      if type == "plot" or type == "coh":
        bankFile = 'trigTemplateBank/' + self.inputIfo + '-TRIGBANK_FOLLOWUP_' + type + str(trig.eventID) + '.xml.gz'
      if type == "chia":
        bankFile = 'trigTemplateBank/' + self.inputIfo + '-TRIGBANK_FOLLOWUP_coh' + str(trig.eventID) + '.xml.gz'
      if type == "notrig":
        bankFile = 'trigTemplateBank/' + ifo + '-TRIGBANK_FOLLOWUP_' + type + str(trig.eventID) + '.xml.gz'
      self.set_bank(bankFile)

      if not ifo == self.inputIfo and not type == "coh" and not type == "chia":
        second_user_tag = "_" + ifo + "tmplt"
      else:
        second_user_tag = ""
      self.set_user_tag("FOLLOWUP_" + str(trig.eventID) + second_user_tag)
      self.__usertag = "FOLLOWUP_" + str(trig.eventID) + second_user_tag


      # THIS IS A HACK FOR NOW, THERE IS PROBABLY A BETTER WAY TO DO THIS
      if (type == 'head'): 
        subBankSize = string.strip(cp.get('followup-inspiral-head','bank-veto-subbank-size'))
        if opts.inspiral_head:
          bankFileName = fu_utils.generateBankVetoBank(trig, ifo, str(trig.gpsTime[ifo]), sngl_table[ifo],int(subBankSize),'BankVetoBank')
        else: bankFileName = 'none'      
        self.add_var_opt("bank-veto-subbank-size", string.strip(cp.get('followup-inspiral-head','bank-veto-subbank-size')))
        self.add_var_opt("order", string.strip(cp.get('followup-inspiral-head','order')))
        self.set_bank(bankFileName)


      # the output_file_name is required by the child job (plotSNRCHISQNode)
      if type == "plot" or type == "notrig" or type == "coh" or type == "chia":
        self.output_file_name = inspJob.outputPath + self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)) + extension

      self.set_id(self.inputIfo + "-INSPIRAL_" + self.__ifotag + "_" + self.__usertag + "-" + self.__start + "-" + str(int(self.__end)-int(self.__start)))

      self.outputCache = self.inputIfo + ' ' + 'INSPIRAL' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name  + '\n' + self.inputIfo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.__start) + ' ' + str(int(self.__end)-int(self.__start)) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

      self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
      self.add_var_opt("output-path",inspJob.outputPath)

      if not opts.disable_dag_categories:
        self.set_category(inspJob.name.lower())

      try:
        if d_node.validNode and eval('opts.' + datafindCommand):
          self.add_parent(d_node)
      except: 
        print >> sys.stderr, "Didn't find a datafind job, I'll assume I don't need it"

      if type == "plot" or type == "notrig":
        if opts.inspiral:
          dag.addNode(self,'inspiral')
          self.validate()
        else: self.invalidate()

      if type == 'head':
        if opts.inspiral_head:
          dag.addNode(self,'inspiral-head')
          self.validate()
        else: self.invalidate()

      if type == 'coh':
        if opts.coh_inspiral:
          dag.addNode(self,'coh-inspiral')
          self.validate()
        else: self.invalidate()

      if type == "chia":
        if opts.plot_chia:
          dag.addNode(self,'chia-inspiral')
          self.validate()
        else: self.invalidate()

    else: #except:
      try:
        print "couldn't add inspiral job for " + self.inputIfo + "@ "+ str(trig.gpsTime[ifo])
        # if self.inputIfo does not exist (happens when inspiral cache and xml files not available), then use ifo in the string.
      except:
        print "couldn't add inspiral job for " + ifo + "@ "+ str(trig.gpsTime[ifo])

  def checkInjections(self,cp):
    try:
      if len(string.strip(cp.get('followup-triggers','injection-file'))) > 0:
        injectionFile = string.strip(cp.get('followup-triggers','injection-file'))
      else:
        injectionFile = None
      return(injectionFile)
    except:
      print >> sys.stderr, "ERROR: failure in followUpInspNode.checkInjections()"
      return None

########## PLOT SNR  CHISQ TIME SERIES ########################################
###############################################################################

##############################################################################
# jobs class for plot snr chisq 

class plotSNRCHISQJob(pipeline.CondorDAGJob,webTheJob):
  """
  A followup plotting job for snr and chisq time series
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "plotsnrchisq":"plotsnrchisq_pipe"
      }
    }
  def __init__(self, options, cp, tag_base='PLOT_FOLLOWUP'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'plotSNRCHISQJob'
    self.__executable = string.strip(cp.get('condor','plotsnrchisq'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)

##############################################################################
# node class for plot snr chisq

class plotSNRCHISQNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a plotSNRCHISQ followup job
  """
  def __init__(self,job,ifo,fileName,trig,page,dag,inspiralNode,opts,ifoString=None):
    """
    job = A CondorDAGJob that can run an instance of plotSNRCHISQ followup.
    """
    if ifoString:
      time = trig.gpsTime[ifoString]
    else:
      time = trig.gpsTime[ifo]
    self.friendlyName = 'Plot SNR/CHISQ/PSD'
    try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.output_file_name = ""
      self.add_var_opt("frame-file",fileName.replace(".xml",".gwf").strip(".gz"))
      self.add_var_opt("inspiral-xml-file",fileName)

      duration = 2.0 # width of the time series to be displayed
      self.add_var_opt("plot-width",duration)

      self.add_var_opt("gps",time)
      self.add_var_opt("gps-start-time",time-duration*.5)
      self.add_var_opt("gps-end-time",time+duration*.5)

      self.add_var_opt("ifo-times",ifo)
      self.add_var_opt("ifo-tag","FOLLOWUP_" + ifo)

      if ifoString:
        self.add_var_opt("user-tag",ifoString+'tmplt_'+str(trig.eventID))
        self.id = job.name + '-' + ifo + '-' + ifoString + 'tmplt' + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      else:
        self.add_var_opt("user-tag",str(trig.eventID))
        self.id = job.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.setupNodeWeb(job,True, dag.webPage.lastSection.lastSub,page,None,None)

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      if inspiralNode.validNode: self.add_parent(inspiralNode)
      if opts.plots:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()
    except: 
      self.invalidate()
      print "couldn't add plot job for " + str(ifo) + "@ "+ str(time)

##############################################################################
# job class for producing the skymap

class lalapps_skyMapJob(pipeline.CondorDAGJob,webTheJob):
  """
  Generates sky map data
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "lalapps_skymap":"lalapps_skymap"
      }
    }
  def __init__(self, options, cp, tag_base='SKY_MAP'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'lalapps_skyMapJob'
    self.__executable = string.strip(cp.get('condor','lalapps_skymap'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)
    self.ra_res = string.strip(cp.get('skymap','ra-res'))
    self.dec_res = string.strip(cp.get('skymap','dec-res'))
    self.sample_rate = string.strip(cp.get('coh-inspiral','sample-rate'))

##############################################################################
# job class for producing the skymap

class pylal_skyPlotJob(pipeline.CondorDAGJob,webTheJob):
  """
  Plots the sky map output of lalapps_skymap
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "pylal_skyPlotJob":"pylal_plot_inspiral_skymap"
      }
    }
  def __init__(self, options, cp, tag_base='SKY_PLOT'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'pylal_skyPlotJob'
    self.__executable = string.strip(cp.get('condor','pylal_skyPlotJob'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)
    self.ra_res = string.strip(cp.get('skymap','ra-res'))
    self.dec_res = string.strip(cp.get('skymap','dec-res'))
    self.sample_rate = string.strip(cp.get('coh-inspiral','sample-rate'))


##############################################################################
# job class for producing the skymap

class lalapps_skyMapNode(pipeline.CondorDAGNode,webTheNode):
  """
  A C code for computing the sky map
  An example command line is:

lalapps_skymap --h1-frame-file H1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.gwf --l1-frame-file L1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.gwf --v1-frame-file V1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088205-2048.gwf --event-id 866088314000001908 --ra-res 512 --dec-res 256 --h1-xml-file H1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.xml.gz --l1-xml-file L1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088022-2048.xml.gz --v1-xml-file V1-INSPIRAL_SECOND_H1H2L1V1_FOLLOWUP_866088314000001908-866088205-2048.xml.gz --output-file chad.txt
  """
  def __init__(self,job,trig,opts):
    self.ifo_list = ["H1","L1","V1"]
    #self.already_added_ifo_list = []

    self.ra_res = job.ra_res
    self.dec_res = job.dec_res
    self.sample_rate = job.sample_rate
    pipeline.CondorDAGNode.__init__(self,job)
    self.friendlyName = 'Produce sky map of event'    
    self.id = job.name + '-skymap-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(job)
    # required by pylal_skyPlotNode
    # this program now gzips its files (otherwise they are really huge)
    self.output_file_name = job.outputPath + self.id+".txt.gz"
    self.add_var_opt("output-file",self.output_file_name)
    self.add_var_opt("ra-res",self.ra_res)
    self.add_var_opt("dec-res",self.dec_res)
    self.add_var_opt("event-id",trig.eventID)
    self.add_var_opt("h1-frame-file","none")
    self.add_var_opt("h1-xml-file","none")
    self.add_var_opt("h2-frame-file","none")
    self.add_var_opt("h2-xml-file","none")
    self.add_var_opt("l1-frame-file","none")
    self.add_var_opt("l1-xml-file","none")
    self.add_var_opt("v1-frame-file","none")
    self.add_var_opt("v1-xml-file","none")
    self.add_var_opt("sample-rate",self.sample_rate)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

  def append_insp_node(self,inspNode,ifo):
    if ifo in self.ifo_list:
      fileName = str(inspNode.output_file_name)
      self.add_var_opt(ifo.lower()+"-frame-file",str(fileName.replace(".xml",".gwf").strip(".gz")))
      self.add_var_opt(ifo.lower()+"-xml-file",str(fileName))
      if inspNode.validNode: self.add_parent(inspNode)
      
    else: pass #print >> sys.stderr, "WARNING: Already added " + ifo


  def add_node_to_dag(self,dag,opts,trig):
    if opts.sky_map:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: 
      self.invalidate()
      #print "couldn't add sky map job for " + str(trig.eventID)



##############################################################################
# job class for producing the skymap

class pylal_skyPlotNode(pipeline.CondorDAGNode,webTheNode):
  """
  A python code for plotting the sky map
  An example command line is

  /pylal_plot_inspiral_skymap --event-id 866088314000001908 --ra-res 512 --dec-res 256 --output-path . --page-rel-path . --output-web-file test.html --page . --injection-right-ascension 0 --injection-declination 0 --map-data-file chad.txt 
  """
  def __init__(self,job,trig,skyMapNode,dag,page,opts):
    # Always initialize the CondorDAGNode
    pipeline.CondorDAGNode.__init__(self,job)
    
    self.friendlyName = 'Produce a plot of the sky map of an event'
    
    self.id = job.name + '-skymap-plot' + str(trig.statValue) + '_' + str(trig.eventID)
    # This node outputs pretty pictures, so we need to tell the setupNodeWeb()
    # method where to put these things.  We'll put it in the last section 
    # not the last subsection since this is an "event" plot not a single ifo
    # trigger plot
    self.setupNodeWeb(job,True, dag.webPage.lastSection,page,None,None)
    # this is the output of the skyMapNode.  It contains the data to make a
    # sky map plot.  (RA,DEC,Probability)
    # an example sky map plotting command line is:
    #

    self.add_var_opt("map-data-file",skyMapNode.output_file_name)
    self.add_var_opt("user-tag",str(trig.eventID))
    self.add_var_opt("ifo-tag",trig.ifos)
    self.add_var_opt("ifo-times",trig.ifos)
    self.add_var_opt("ra-res",str(skyMapNode.ra_res))
    self.add_var_opt("dec-res",str(skyMapNode.dec_res))
    self.add_var_opt("stat-value", str(trig.statValue))
    # if this is a software injection pass along the information to the
    # plotting code so that it can make a mark where the injection should have
    # been :)
    if trig.is_found():
      inj_ra = trig.coincs.sim.longitude
      inj_dec = trig.coincs.sim.latitude
      self.add_var_opt("injection-right-ascension",str(inj_ra))
      self.add_var_opt("injection-declination",str(inj_dec))

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    try:
      if skyMapNode.validNode: self.add_parent(skyMapNode)
    except: pass
    if opts.sky_map_plot:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()
 

############### DATAFIND CLASSES ##############################################
###############################################################################

class followupDataFindJob(pipeline.LSCDataFindJob,webTheJob):
  defaults={
    "section":"condor",
    "options":
      {
      "universe":"vanilla",
      "datafind":"ligo_data_find"
      }
    }

  def __init__(self, config_file, source):

    if source == 'futrig':
      self.name = 'qscanDataFindJob'
    if source == 'inspiral':
      self.name = 'inspiralDataFindJob'

    # unfortunately the logs directory has to be created before we call LSCDataFindJob
    try:
      os.mkdir(self.name)
      os.mkdir(self.name + '/logs')
    except: pass
    pipeline.LSCDataFindJob.__init__(self, self.name, self.name + '/logs', config_file)
    if source == 'futrig':
      self.setup_cacheconv(config_file)
    self.setupJobWeb(self.name) # this will overwrite the stdout and stderr set up by LSCDataFindJob

  def setup_cacheconv(self,cp):
    # create a shell script to call convertlalcache.pl if the value of $RETURN is 0
    convert_script = open(self.name + '/cacheconv.sh','w')
    convert_script.write("""#!/bin/bash
    if [ ${1} -ne 0 ] ; then
      exit 1
    else
      %s ${2} ${3}
    fi
    """ % string.strip(cp.get('condor','convertcache')))
    convert_script.close()
    os.chmod(self.name + '/cacheconv.sh',0755)


class followupDataFindNode(pipeline.LSCDataFindNode,webTheNode):
 
  def __init__(self, job, source, type, cp, time, ifo, opts, dag, datafindCommand, procParams=None):
    try:
      self.outputFileName = ""
      pipeline.LSCDataFindNode.__init__(self,job)
      self.id = str(ifo) + '-' + repr(time) + '-' + str(type)
      self.setupNodeWeb(job,False,None,None,None,dag.cache)
      if source == 'futrig':
        self.outputFileName = self.setup_fu_trig(job, cp, time, ifo, type)
        nodeName = "qscan data find"
      if source == 'inspiral':
        self.outputFileName = self.setup_inspiral(cp,ifo,type,procParams)
        nodeName = "inspiral data find"

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      # if the selected "ifo" needs to be done remotely (this the case for 
      # Virgo qscan datafind) do not add the node to the dag
      if eval('opts.' + datafindCommand) and \
        not( cp.has_option("followup-"+type,"remote-ifo") and \
        cp.get("followup-"+type,"remote-ifo")==ifo ):
          dag.addNode(self,nodeName)
          self.validNode = True
      else: self.validNode = False
    except:
      self.validNode = False
      print >> sys.stderr, "could not set up the datafind jobs for " + type

  def setup_inspiral(self,cp,ifo,type,procParams):
    for row in procParams:
      param = row.param.strip("-")
      value = row.value
      if param == 'gps-start-time': startTime = value
      if param == 'gps-end-time': endTime = value
      if param == 'pad-data': paddataTime = value
    self.set_observatory(ifo[0])
    self.set_start(int(startTime) - int(paddataTime))
    self.set_end(int(endTime) + int(paddataTime))
    self.set_type(cp.get("followup-"+type,ifo + '_type'))
    lalCache = self.get_output()
    return(lalCache)

  def setup_fu_trig(self, job, cp, time, ifo, type):
    # 1s is substracted to the expected startTime to make sure the window
    # will be large enough. This is to be sure to handle the rouding to the
    # next sample done by qscan.
    self.q_time = cp.getint("followup-"+type,'search-time-range')/2
    self.set_observatory(ifo[0])
    self.set_start(int( time - self.q_time - 1))
    self.set_end(int( time + self.q_time + 1))
    if cp.has_option("followup-"+type, ifo + '_type'): 
      self.set_type( cp.get("followup-"+type, ifo + '_type' ))
    else:
      if not( cp.has_option("followup-"+type,"remote-ifo") and \
      cp.get("followup-"+type,"remote-ifo")==ifo ):
        self.set_type( cp.get("followup-"+type, 'type' ))
      else: self.set_type("dummy")
    lalCache = self.get_output()
    qCache = lalCache.rstrip("cache") + "qcache"
    self.set_post_script(job.name + "/cacheconv.sh $RETURN %s %s" %(lalCache,qCache) )
    return(qCache)

##############################################################################
# qscan class for qscan jobs

class qscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A qscan job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "qscan":"wpipeline"
      }
    }

  def __init__(self, opts, cp, tag_base='QSCAN'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__executable = string.strip(cp.get('condor','qscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(tag_base,None,cp)
    self.setup_checkForDir()

  def setup_checkForDir(self):
    # create a shell script to check for the existence of the qscan output directory and rename it if needed
    checkdir_script = open(self.name + '/checkForDir.sh','w')
    checkdir_script.write("""#!/bin/bash
    if [ -d $1/$2 ]
    then
      matchingList=$(echo $(find $1 -name $2.bk*))
      COUNTER=1
      for file in $matchingList
      do
        let COUNTER=COUNTER+1
      done
      mv $1/$2 $1/$2.bk.$COUNTER
    fi
    """)
    checkdir_script.close()
    os.chmod(self.name + '/checkForDir.sh',0755)

##############################################################################
# qscan class for qscan Node

class qscanNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,time,cp,qcache,ifo,name, opts, d_node, dag, datafindCommand, qscanCommand, trig=None,qFlag=None):
    """
    job = A CondorDAGJob that can run an instance of qscan.
    """
    self.friendlyName = name
    self.id = ifo + '-' + name + '-' + repr(time)

    pipeline.CondorDAGNode.__init__(self,job)
    if name.split('-')[0]=='background':
      self.add_var_arg('scan')
    else:
      self.add_var_arg('scan -r')
    qscanConfig = string.strip(cp.get("followup-"+name, ifo + 'config-file'))
    self.add_var_arg("-c "+qscanConfig)
    self.add_var_arg("-f "+qcache)

    if cp.has_option("followup-"+name, ifo + 'output') and string.strip(cp.get("followup-"+name, ifo + 'output')):
      output = string.strip(cp.get("followup-"+name, ifo + 'output'))
    else:
      #output = dag.publish_path + '/' + job.name + '/' + name + '/' + ifo
      output = job.name + '/' + name + '/' + ifo
    if not os.access(output,os.F_OK):
      os.makedirs(output)
    else:
      if not os.access(output,os.W_OK):
        print >> sys.stderr, 'path '+output+' is not writable'
        sys.exit(1)

    self.add_var_arg("-o "+output+"/"+repr(time))
    self.add_var_arg(repr(time))

    self.set_pre_script(job.name + "/checkForDir.sh %s %s" \
    %(output, repr(time)))

    #get the absolute output path whatever the path might be in the ini file
    #absoutput = os.path.abspath(output)

    #self.outputName = absoutput + '/' + repr(time) # redirect output name
    self.outputName = output + '/' + repr(time)

    #prepare the string for the output cache
    self.outputCache = ifo + ' ' + name + ' ' + repr(time) + ' ' + self.outputName + '\n'

    #extract web output from the ini file if the job is QSCAN
    if job.name == 'QSCAN':
      if cp.has_option("followup-"+name,ifo+'web') and string.strip(cp.get("followup-"+name,ifo+'web')):
        pageOverride = string.strip(cp.get("followup-"+name,ifo+'web'))+'/'+repr(time)
      else:
        #pageOverride = dag.page + '/' + job.name + '/' + name + '/' + ifo + '/' + repr(time)
        pageOverride = job.name + '/' + name + '/' + ifo + '/' + repr(time)
      self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,dag.page,pageOverride,dag.cache)

    else:
      self.setupNodeWeb(job,False,None,None,None,dag.cache) 

    # This command will force Condor to see the qscan jobs successful even
    # they fail. This is useful when the followups are rerun on candidates 
    # already analysed, since when a qscan directory exists, the qscan job
    # will fail. By setting the post_script to true the qscan job will
    # still be reported as successful, so that an analyseQscan job can be run
    # immediately after. 
    # self.set_post_script("/bin/true")

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    # only add a parent if it exists
    try:
      if d_node.validNode and eval('opts.' + datafindCommand):
        self.add_parent(d_node)
    except: pass

    # if the selected "ifo" needs to be done remotely (this the case for 
    # Virgo qscans) do not add the node to the dag
    if eval('opts.' + qscanCommand):
      if not(cp.has_option("followup-"+name,"remote-ifo") and \
      cp.get("followup-"+name,"remote-ifo")==ifo):
        dag.addNode(self,self.friendlyName)
        self.validNode = True
      else:
        dag.remote_nodes.append(self)
    else: self.validNode = False
 #   except: 
 #     self.validNode = False
 #     print >> sys.stderr, "could not set up the qscan job for " + self.id


##############################################################################
# class for remote qscan jobs

class remoteQscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A remote qscan job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "submit_remote_scan":"submit_remote_scan.py"
      }
    }

  def __init__(self, opts, cp, tag_base='REMOTESCAN'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    
    if not os.path.exists(tag_base):
       os.mkdir(tag_base)
    self.setup_executable(tag_base)
    self.__executable = tag_base + '/remote_scan_wrapper.sh'
    self.__universe = "scheduler"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(tag_base,None)

  def setup_executable(self,tag_base):
    starter_script = open(tag_base + '/remote_scan_wrapper.sh','w')
    starter_script.write("""#!/bin/bash
    dotbashrc=$1
    executable=$2
    gpstime=$3
    configfile=$4
    qscantype=$5
    remoteoutput=$6
    remotereceiver=$7
    outputpath=$8
    shift 8
    source $dotbashrc
    $executable --gps-time $gpstime --config-file $configfile --qscan-type $qscantype --remote-output $remoteoutput --remote-receiver $remotereceiver --output-path $outputpath
    """)
    starter_script.close()
    os.chmod(tag_base + '/remote_scan_wrapper.sh',0755)


##############################################################################
# class for remote qscan Node

class remoteQscanFgNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a remote qscan job
  """
  def __init__(self,job,time,cp,ifo,name,opts,dag,qscanCommand):
    """
    job = A CondorDAGJob that can run an instance of remote qscan.
    """
    self.friendlyName = name
    self.id = ifo + '-' + name + '-' + repr(time)
    pipeline.CondorDAGNode.__init__(self,job)
    self.add_macro("macroid", self.id)
    self.jobName = job.name

    self.add_var_arg(string.strip(cp.get("followup-remote-scan","virgo-env-path")))
    self.add_var_arg(string.strip(cp.get("condor","submit_remote_scan")))
    self.add_var_arg(repr(time))
    self.add_var_arg(string.strip(cp.get("followup-"+name,ifo+"config-file")))
    self.add_var_arg("_".join(name.split("-")[1:len(name.split("-"))]))
    self.add_var_arg(string.strip(cp.get("followup-remote-scan","remote-output")))
    self.add_var_arg(string.strip(cp.get("followup-remote-scan","remote-server")))

    outputdir = 'QSCAN' + '/' + name + '/' + ifo + '/' + repr(time)
    self.add_var_arg(outputdir)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    if eval('opts.' + qscanCommand):
        dag.addNode(self,"Remote " + self.friendlyName)
        self.validate()
    else: self.invalidate()


##############################################################################
# distributeQscanJob class: the job

class distributeQscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A job to distribute the results of the qscans that have been run remotely (for LV search)
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "distribute_q":"distrib_fu_qscan_results.py"
      }
    }

  def __init__(self,cp):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'distributeQscanJob'
    self.__executable = string.strip(cp.get('condor','distribute_q'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__)

##############################################################################
# distributeQscanNode class: the node

class distributeQscanNode(pipeline.CondorDAGNode, webTheNode):
  """
  A node to distribute the results of the qscans that have been run remotely (for LV search)
  """
  def __init__(self,job,foregroundCache,backgroundCache,ifo,inputFile,opts,dag):

    self.friendlyName = "distributeQscanResults"

    pipeline.CondorDAGNode.__init__(self,job)
    self.add_var_opt('qscan-input-file',inputFile)
    self.add_var_opt('qscan-cache-background',backgroundCache)
    self.add_var_opt('qscan-cache-foreground',foregroundCache)
    self.add_var_opt('remote-ifo',ifo)

    typeList=""
    for type in ["qscan","seismic-qscan"]:
      typeList += type + ","
    self.add_var_opt('qscan-type-list',typeList.strip(','))

    if opts.distrib_remote_q:
      dag.addNode(self,self.friendlyName)
      self.validNode = True
    else: self.validNode = False

##############################################################################
# analyse qscan class: the job

class analyseQscanJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup analyseQscan job to interprete the qscans
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "analyseQscan":"analyseQscan.py"
      }
    }
  def __init__(self,options,cp,tag_base='ANALYSE_QSCAN'):
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'analyseQscanJob'
    self.__executable = string.strip(cp.get('condor','analyseQscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)

##############################################################################
# analyse qscan class: the node

class analyseQscanNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a followup analyseQscan job
  """
  def __init__(self,job,time,ifo,name,foregroundCache,backgroundCache,cp,opts,dag,command):
    """
    job = A CondorDAGJob that can run an instance of analyseQscan followup.
    """
    self.friendlyName = 'analyse ' + name
    self.id = ifo + '-' + name + '-' + repr(time)

    nameList = name.split('-')[1:len(name.split('-'))]
    shortName = ''
    for word in nameList:
      shortName = shortName + word + '-'

    try:
      pipeline.CondorDAGNode.__init__(self,job)
      if cp.has_option('followup-analyse-qscan','generate-qscan-xml'):
        self.add_var_opt('generate-qscan-xml','')
      self.add_var_opt('z-threshold',cp.getfloat('followup-analyse-qscan','z-threshold'))
      if cp.has_option('followup-analyse-qscan','plot-z-distribution'):
        self.add_var_opt('plot-z-distribution','')
        self.add_var_opt('z-min',cp.getfloat('followup-analyse-qscan','z-min'))
        self.add_var_opt('z-max',cp.getfloat('followup-analyse-qscan','z-max'))
        self.add_var_opt('z-bins',cp.getfloat('followup-analyse-qscan','z-bins'))
      if cp.has_option('followup-analyse-qscan','plot-dt-distribution'):
        self.add_var_opt('plot-dt-distribution','')
        self.add_var_opt('dt-min',cp.getfloat('followup-analyse-qscan',shortName + 'dt-min'))
        self.add_var_opt('dt-max',cp.getfloat('followup-analyse-qscan',shortName + 'dt-max'))
        self.add_var_opt('dt-bins',cp.getfloat('followup-analyse-qscan','dt-bins'))
      if cp.has_option('followup-analyse-qscan','plot-z-scattered'):
        self.add_var_opt('plot-z-scattered','')
      if cp.has_option('followup-analyse-qscan','plot-z-scattered') or cp.has_option('followup-analyse-qscan','plot-dt-distribution'):
        if not ifo=='V1':
          refChannel = cp.get('followup-analyse-qscan',shortName + 'ref-channel').split(',')[0].strip()
        else:
          refChannel = cp.get('followup-analyse-qscan',shortName + 'ref-channel').split(',')[1].strip()
        self.add_var_opt('ref-channel',refChannel)
      self.add_var_opt('ifo-times',ifo)
      self.add_var_opt('type',name)
      self.add_var_opt('gps-string',repr(time))
      self.add_var_opt('ifo-tag',ifo)
      self.add_var_opt('user-tag',repr(time).replace('.','_') + "_" + shortName.replace('-','_').strip("_"))

      self.add_var_opt('qscan-cache-foreground',foregroundCache)
      self.add_var_opt('qscan-cache-background',backgroundCache)

      self.setupNodeWeb(job,True,None,dag.page,None,None)
      # get the table of the qscan job associated to this trigger
      #if not(cp.has_option("followup-"+name,"remote-ifo") and cp.get("followup-"+name,"remote-ifo")==ifo):
        #for node in dag.get_nodes():
          #if isinstance(node,qscanNode):
            #if node.id == self.id:
              # link the analyseQscan output page to the qscan table
              #node.webTable.row[0].cell[0].linebreak()
              #node.webTable.row[0].cell[0].link(self.webLink,"qscan background vs qscan foreground")
              #break
      # if remote-ifo is analysed, find the associated qscan jobs in dag.remote_nodes
      #else:
        #for node in dag.remote_nodes:
          #if isinstance(node,qscanNode):
            #if node.id == self.id:
              #node.webTable.row[0].cell[0].linebreak()
              #node.webTable.row[0].cell[0].link(self.webLink,"qscan background vs qscan foreground")
              #break

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())
      
      # add the parents to this node
      for node in dag.get_nodes():
        # if node distributeQscanNode is valid and if remote ifo is analysed,
        # add distributeQscanNode as parent
        if isinstance(node,distributeQscanNode):
          if cp.has_option("followup-"+name,"remote-ifo") and cp.get("followup-"+name,"remote-ifo")==ifo:
            if node.validNode:
              self.add_parent(node)
        # if node remoteQscanFgNode is valid and if remote ifo is analysed,
        # add remoteQscanFgNode as parent
        if isinstance(node,remoteQscanFgNode):
          if cp.has_option("followup-"+name,"remote-ifo") and cp.get("followup-"+name,"remote-ifo")==ifo:
            if node.friendlyName == name and node.validNode:
              self.add_parent(node)
        # add all qscan nodes of the same type as parents
        if isinstance(node,qscanNode): 
          if node.validNode:
            if (node.friendlyName == name or \
            node.friendlyName.replace('background','foreground') == name) \
            and node.id.split('-')[0] == ifo:
              self.add_parent(node)

      if eval('opts.' + command):
        dag.addNode(self,self.friendlyName)
        self.validNode = True
      else: self.validNode = False

    except:
      self.validNode = False
      print "couldn't add " + name + " analyseQscan job for " + ifo + "@ "+ repr(time)

##############################################################################
# class for h1h2 qevent jobs

class h1h2QeventJob(pipeline.CondorDAGJob, webTheJob):
  """
  A h1h2 qevent job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "qscan":"wpipeline"
      }
    }
  def __init__(self, opts, cp):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.name = 'h1h2QeventJob'
    self.__executable = string.strip(cp.get('condor','qscan'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)    
    self.setupJobWeb(self.name)
    #self.setup_cachecat()
    self.setup_cacheAndCheckDir()

#  def setup_cachecat(self):
#    # create a shell script to cat all the required cache files
#    cat_script = open(self.name + '/cachecat.sh','w')
#    cat_script.write("""#!/bin/bash
#    cat ${1} ${2} > ${3}
#    """)
#    cat_script.close()
#    os.chmod(self.name + '/cachecat.sh',0755)

  def setup_cacheAndCheckDir(self):
    # create a shell script to cat all the required cache files
    # create a shell script to check for the existence of the qscan output directory and rename it if needed
    checkdir_script = open(self.name + '/checkForDir.sh','w')
    checkdir_script.write("""#!/bin/bash
    cat ${1} ${2} > ${3}
    if [ -d $4/$5 ]
    then
      matchingList=$(echo $(find $4 -name $5.bk*))
      COUNTER=1
      for file in $matchingList
      do
        let COUNTER=COUNTER+1
      done
      mv $4/$5 $4/$5.bk.$COUNTER
    fi
    """)
    checkdir_script.close()
    os.chmod(self.name + '/checkForDir.sh',0755)

#############################################################################
# class for h1h2 qevent Node

class h1h2QeventNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a qscan job
  """
  def __init__(self,job,dNode,times,ifoList,name,cp,opts,dag,qeventCommand):
    """
    job = A CondorDAGJob that can run an instance of H1H2 qevent.
    """

    ifoString = ''
    for ifo in ifoList:
      ifoString = ifoString + ifo

    self.friendlyName = name
    self.id = ifoString + '-' + name + '-' + str(times[ifoList[0]])

    pipeline.CondorDAGNode.__init__(self,job)

    cache_type_temp = dNode[ifoList[0]].outputFileName.split('-')[1]
    cache_type = cache_type_temp[3:len(cache_type_temp)]
    cache_start = []
    cache_end = []
    for ifo in ifoList:
      cache_temp = dNode[ifo].outputFileName.split('.')[0]
      cache_start.append(cache_temp.split('-')[2])
      cache_end.append(cache_temp.split('-')[-1])
    cache_start_time = max(cache_start)

    qeventcache = job.name + '/' + ifoString + '_' + cache_type + '-' + \
    str(max(cache_start)) + '-' + str(min(cache_end)) + '.qcache'


    if cp.has_option("followup-"+name, ifoString + '-output') and string.strip(cp.get("followup-"+name, ifoString + '-output')):
      output = string.strip(cp.get("followup-"+name, ifoString + '-output'))
    else:
      #output = dag.publish_path + '/' + job.name + '/' + name + '/' + ifoString
      output = job.name + '/' + name + '/' + ifoString
    if not os.access(output,os.F_OK):
      os.makedirs(output)
    else:
      if not os.access(output,os.W_OK):
        print >> sys.stderr, 'path '+output+' is not writable'
        sys.exit(1)

    self.add_var_arg('event')
    qeventConfig = string.strip(cp.get("followup-"+name, ifoString + '-config-file'))
    self.add_var_arg('-p '+qeventConfig)
    self.add_file_arg('-f '+qeventcache)
    self.add_var_arg('-o '+output+'/'+repr(times[ifoList[0]]))
    self.add_var_arg(repr(times[ifoList[0]]))
    eventDuration = string.strip(cp.get("followup-"+name, 'duration'))
    self.add_var_arg(eventDuration)

    #self.set_pre_script(job.name + "/cachecat.sh %s %s %s" \
    #%(dNode[ifoList[0]].outputFileName, dNode[ifoList[1]].outputFileName, \
    #qeventcache))
    self.set_pre_script(job.name + "/checkForDir.sh %s %s %s %s %s" \
    %(dNode[ifoList[0]].outputFileName, dNode[ifoList[1]].outputFileName, \
    qeventcache, output, repr(times[ifoList[0]])))

    #get the absolute output path whatever the path might be in the ini file
    absoutput = os.path.abspath(output)
    self.outputName = absoutput + '/' + repr(times[ifoList[0]]) # redirect output name

    #prepare the string for the output cache
    self.outputCache = ifoString + ' ' + name + ' ' + repr(times[ifoList[0]]) + ' ' + self.outputName + '\n'

    if cp.has_option("followup-"+name,ifoString+'-web') and string.strip(cp.get("followup-"+name,ifoString+'-web')):
      pageOverride = string.strip(cp.get("followup-"+name,ifoString+'-web'))+'/'+repr(times[ifoList[0]])
    else:
      #pageOverride = dag.page + '/' + job.name + '/' + name + '/' + ifoString + '/' + repr(times[ifoList[0]])
      pageOverride = job.name + '/' + name + '/' + ifoString + '/' + repr(times[ifoList[0]])
    self.setupNodeWeb(job,False,dag.webPage.lastSection.lastSub,dag.page,pageOverride,dag.cache)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    for ifo in ifoList:
      if dNode[ifo].validNode: self.add_parent(dNode[ifo])
      else: pass

    if eval('opts.' + qeventCommand):
      dag.addNode(self,self.friendlyName)
      self.validNode = True
    else: self.validNode = False


###############################################################################
# FrCheck Jobs and Nodes

class FrCheckJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup job for checking frames
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "frame_check":"frame_check"
      }
    }
  def __init__(self, options, cp, tag_base='FRCHECK'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'FrCheckJob'
    self.__executable = string.strip(cp.get('condor','frame_check'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)


class FrCheckNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a FrCheck followup job
  """
  def __init__(self, FrCheckJob, procParams, ifo, trig, cp, opts, dag, datafindCache, d_node, datafindCommand):

    try:
      hipeCache = checkHipeCachePath(cp)

      if not hipeCache:
        cacheFile = datafindCache
      else:
        for row in procParams:
          param = row.param.strip("-")
          value = row.value
          if param == 'frame-cache': cacheFile = value

      self.friendlyName = 'Frame Check'
 
      pipeline.CondorDAGNode.__init__(self,FrCheckJob)
      self.add_var_opt("frame-cache", cacheFile)
      self.add_var_opt("frame-check-executable", string.strip(cp.get('followup-frameCheck','executable')))
      self.add_var_opt("ifo-times",ifo)
      self.add_var_opt("ifo-tag","FOLLOWUP_"+ifo)
      self.add_var_opt("user-tag",trig.eventID)
      self.id = FrCheckJob.name + '-' + ifo + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.setupNodeWeb(FrCheckJob,True, dag.webPage.lastSection.lastSub,dag.page,None,None)

      if not opts.disable_dag_categories:
        self.set_category(FrCheckJob.name.lower())

      try:
        if d_node.validNode and eval('opts.' + datafindCommand):
          self.add_parent(d_node)
      except: pass

      if opts.frame_check:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add frame check job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])

class IFOstatus_checkJob(pipeline.CondorDAGJob, webTheJob):
  """
  A followup job for downloading summary plots
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "IFOstatus_check":"IFOstatus_check"
      }
    }
  def __init__(self, options, cp, tag_base='IFOSTATUS'):
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'IFOstatus_checkJob'
    self.__executable = string.strip(cp.get('condor','IFOstatus_check'))
    self.__universe = "local"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)
    
class IFOstatus_checkNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a FrCheck followup job
  """
  def __init__(self, IFOstatus_checkJob, ifo, trig, cp,opts,dag):

    self.friendlyName = 'IFO status summary plots'
    pipeline.CondorDAGNode.__init__(self,IFOstatus_checkJob)
    self.add_var_opt("ifo-times", ifo)
    self.add_var_opt("gps-time", trig.gpsTime[ifo])
    self.add_var_opt("ifo-tag", "FOLLOWUP_"+ifo)
    self.add_var_opt("user-tag", str(trig.eventID))
    self.id = IFOstatus_checkJob.name + '-' + str(ifo) + '-' + str(trig.statValue) + '_' + str(trig.eventID)
    self.setupNodeWeb(IFOstatus_checkJob,True, dag.webPage.lastSection.lastSub,dag.page,None,None)

    if not opts.disable_dag_categories:
      self.set_category(IFOstatus_checkJob.name.lower())

    if opts.ifo_status_check:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()

##############################################################################

class followupoddsJob(pipeline.CondorDAGJob, webTheJob):
  """
  A model selection job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "followupodds":"lalapps_inspnest"
      }
    }
  def __init__(self,options,cp,tag_base='FOLLOWUPODDS'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)    
    self.__prog__='followupoddsjob'
    self.__executable=string.strip(cp.get('condor','followupodds'))
    self.__universe="standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__prog__,tag_base)

class followupoddsNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of the model selection followup job
  """
  def __init__(self,followupoddsJob,procParamsTable,trig,randomseed,cp,opts,dag):
    try:
      IFOs = trig.ifolist_in_coinc
      time_prior = string.strip(cp.get('followup-odds','time_prior'))
      #Nlive = string.strip(cp.get('followup-odds','live-points'))
      Nlive = string.strip(cp.get('followup-odds','min-live'))
      Nmcmc = string.strip(cp.get('followup-odds','Nmcmc'))
      srate = string.strip(cp.get('followup-odds','sample_rate'))
      Approximant = string.strip(cp.get('followup-odds','approximant'))
      self.friendlyName = 'Odds followup job'
      pipeline.CondorDAGNode.__init__(self,followupoddsJob)
      cacheFiles=[]
      GPSstarts=[]
      GPSends=[]
      for ifo in IFOs:
        for row in procParamsTable[ifo]:
          param=row.param.strip("-")
          value=row.value
          if param == 'frame-cache': cacheFile=value
          if param == 'gps-start-time':
            GPSstarts.append(float(value))
          if param == 'gps-end-time':
            GPSends.append(float(value))
        self.add_var_arg("--IFO "+str(ifo))
        self.add_var_arg("--cache " +str(cacheFile))
      #Check the start and end times are OK
      
      GPSstart=str(max(GPSstarts)+64)
      GPSend=str(min(GPSends)-64)

      #Check the start and end times are OK
      GPSstart=str(max(GPSstarts)+64)
      GPSend=str(min(GPSends)-64)

      outputname = followupoddsJob.name + '/'+followupoddsJob.name+'-' \
                   +trig.ifos+'-'+str(trig.statValue)+'_'+str(trig.eventID)+'_'+randomseed[0]+'.dat'
      self.add_var_opt("Nlive",Nlive)
      self.add_var_opt("GPSstart",GPSstart)
      self.add_var_opt("length",str(float(GPSend)-float(GPSstart)))
      self.add_var_opt("approximant",Approximant)
      self.add_var_opt("out",outputname)
      self.add_var_opt("Nsegs",str((int(float(GPSend))-int(float(GPSstart)))/8))
      self.add_var_opt("dt",time_prior)
      self.add_var_opt("end_time",trig.gpsTime[ifo])
      self.add_var_opt("Mmin",2.8)
      self.add_var_opt("Mmax",30)
      self.add_var_opt("srate",srate)
      self.add_var_opt("seed",randomseed[0])
      #self.add_var_opt("randomseed","[" + randomseed[0] + "," + randomseed[1] + "]")
      self.id = followupoddsJob.name + '-' + trig.ifos + '-' + str(trig.statValue) + '_' + str(trig.eventID) + '_' + randomseed[0]
      self.outputCache = trig.ifos + ' ' + followupoddsJob.name + ' ' +\
                         self.id.split('-')[-1]+' '+outputname+'\n'
      self.add_var_opt("channel",string.strip(cp.get("followup-coh-trigbank",trig.ifos[0:2]+"_channel")))

      #print "Using IFOs " + str(IFOs)

      self.setupNodeWeb(followupoddsJob,False,None,None,None,dag.cache)

      #print "Arguments: " + str(self.get_cmd_line())

      if opts.odds and float(GPSend)-float(GPSstart)>=24:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "Couldn't add followupOdds job for " + str(trig.gpsTime[ifo])
  

###########################################################################

class followupOddsPostJob(pipeline.CondorDAGJob,webTheJob):
  """
  The post-processing of odds jobs
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "oddsPostScript":"OddsPostProc.py"
      }
    }
  def __init__(self,options,cp,tag_base='FOLLOWUPODDSPOST'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__='followupOddsPostJob'
    self.__executable=string.strip(cp.get('condor','oddsPostScript'))
    self.__universe="vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__prog__,tag_base)


##############################################################################

class followupOddsPostNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs the post-processing script
  """
  def __init__(self,oddsPostJob,procParams,trig,oddsjoblist,cp,opts,dag):
    try:
      self.friendlyName = 'Odds plotting job'
      pipeline.CondorDAGNode.__init__(self,oddsPostJob)

      # get the list of odds .txt files to be used as input
      for oddsjobId in oddsjoblist:
        oddsfile = oddsjobId.split('-')[0]+'/' + oddsjobId + '.dat' # here we assume that the directory name is mcmcId.split('-')[0]
        self.add_var_arg("--data " + oddsfile)

      # Get the number of live points in each run
      Nlive = string.strip(cp.get('followup-odds','min-live'))
      self.add_var_opt("Nlive",Nlive)

      if cp.has_option('followup-odds','web') and string.strip(cp.get('followup-odds','web')):
        outputpath = string.strip(cp.get('followup-odds','web'))
      else:
        outputpath = oddsPostJob.name + "/" + str(trig.eventID)
      if not os.access(outputpath,os.F_OK):
        os.mkdir(outputpath)

      self.add_var_opt("outpath",outputpath)

      self.id = oddsPostJob.name + '-' + trig.ifos + '-' + str(trig.statValue) + '_' + str(trig.eventID)

      #output_page = self.id
      self.outputCache = self.id.replace('-',' ') + " " + os.path.abspath(outputpath) + "/" + self.id + "\n"

      self.setupNodeWeb(oddsPostJob,False,dag.webPage.lastSection.lastSub,None,None,dag.cache)
 
      # only add a parent if it exists
      for node in dag.get_nodes():
        if isinstance(node,followupoddsNode):
          if not node.id.find(trig.ifos + '-' + str(trig.statValue) + '_' + str(trig.eventID)) == -1:
            try:
              if node.validNode: self.add_parent(node)
            except: pass

      if opts.odds:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else:
        self.invalidate()

    except:
      self.invalidate()
      print "couldn't add odds post job for " + str(trig.ifos) + "@ "+ str(trig.gpsTime[trig.ifolist_in_coinc[-1]])


##############################################################################

class followupmcmcJob(pipeline.CondorDAGJob, webTheJob):
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "followupmcmc":"lalapps_followupMcmc"
      }
    }
  """
  An mcmc job
  """
  def __init__(self, options, cp, tag_base='FOLLOWUPMCMC'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'followupmcmcJob'
    self.__executable = string.strip(cp.get('condor','followupmcmc'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__prog__,tag_base)

###############################################################################

class followupmcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of an mcmc followup job
  """
  #def __init__(self, followupmcmcJob, procParams, ifo, trig, randomseed, cp,opts,dag):
  def __init__(self, followupmcmcJob, procParams, trig, randomseed, cp, opts, dag, ifo=None):

    try:
      time_margin = string.strip(cp.get('followup-mcmc','prior-coal-time-marg'))
      iterations = string.strip(cp.get('followup-mcmc','iterations'))
      tbefore = string.strip(cp.get('followup-mcmc','tbefore'))
      tafter = string.strip(cp.get('followup-mcmc','tafter'))
      massmin = string.strip(cp.get('followup-mcmc','massmin'))
      massmax = string.strip(cp.get('followup-mcmc','massmax'))
      dist90 = string.strip(cp.get('followup-mcmc','dist90'))
      dist10 = string.strip(cp.get('followup-mcmc','dist10'))

      self.friendlyName = 'MCMC followup'
      pipeline.CondorDAGNode.__init__(self,followupmcmcJob)

      if ifo:
        IFOs = [ifo]
        self.ifonames = ifo
      else:
        IFOs = trig.ifolist_in_coinc
        self.ifonames = trig.ifos

      cacheFiles = ""
      channelNames = ""
      chunk_end_list = {}
      chunk_start_list = {}
      for itf in IFOs:
        for row in procParams[itf]:
          param = row.param.strip("-")
          value = row.value
          if param == 'frame-cache': cacheFile = value
          if param == 'channel-name': channel = value
          if param == 'gps-end-time': chunk_end = value
          if param == 'gps-start-time': chunk_start = value
        cacheFiles += cacheFile + ","
        channelNames += channel + ","
        chunk_end_list[itf] = chunk_end
        chunk_start_list[itf] = chunk_start

      if len(IFOs) > 1:
        maxSNR = 0
        maxIFO = ""
        for trigger in trig.coincs:
          snr = trigger.snr
          if snr > maxSNR:
            maxSNR = snr
            maxIFO = trigger.ifo
      else:
        maxIFO = IFOs[0]
      trig_tempo = getattr(trig.coincs,maxIFO)
      triggerRef = copy.deepcopy(trig_tempo)
      self.ifoRef = maxIFO
        
      self.add_var_opt("template",string.strip(cp.get('followup-mcmc','template')))
      self.add_var_opt("iterations",iterations)
      self.add_var_opt("randomseed",randomseed)
      self.add_var_opt("tcenter","%0.3f"%trig.gpsTime[maxIFO])
      self.add_var_opt("tbefore",tbefore)
      self.add_var_opt("tafter",tafter)

      tmin = trig.gpsTime[maxIFO] - float(time_margin)
      tmax = trig.gpsTime[maxIFO] + float(time_margin)
      self.add_var_opt("priorparameters","[" + massmin + "," + massmax + "," + str(tmin) + "," + str(tmax) + "," + dist90 + "," + dist10 + "]")

      param_mchirp = triggerRef.mchirp
      param_eta = triggerRef.eta
      param_distance = triggerRef.eff_distance
      self.add_var_opt("guess","[" + str(param_mchirp) + "," + str(param_eta) + "," + str(trig.gpsTime[maxIFO]) + "," + str(param_distance) + "]")
      #self.add_var_opt("fixed","[altitude=0,azimuth=0,inclination=0,polarisation=0]")

      #self.add_var_opt("readdata","")

      ########################################################################
      # GET THE FRAME FILE INFO - THIS NEEDS TO BE CHANGED !!!
      # THE MCMC CODE SHOULD TAKE THE SAME PSD AS LALAPPS_INSPIRAL.
      ########################################################################

      #self.add_var_opt("cachefiledata",cacheFile)
      #self.add_var_opt("cachefilenoise",cacheFile)
      self.add_var_opt("cachefile","["+cacheFiles.strip(",")+"]")
      self.add_var_opt("filechannel","["+channelNames.strip(",")+"]")

      psdEstimateStart = ""
      psdEstimateEnd = ""
      for itf in IFOs:
        datainchunk_before = int(trig.gpsTime[maxIFO]) - 75 - 64 - int(chunk_start_list[itf])
        datainchunk_after = int(chunk_end_list[itf]) - 64 - int(trig.gpsTime[maxIFO]) - 32
        if datainchunk_after > datainchunk_before:
          psdEstimateStart += str(int(trig.gpsTime[maxIFO]) + 32) + ","
          psdEstimateEnd += str(int(chunk_end_list[itf]) - 64) + ","
        else:
          psdEstimateStart += str(int(chunk_start_list[itf]) + 64) + ","
          psdEstimateEnd += str(int(trig.gpsTime[maxIFO]) - 75) + ","

      self.add_var_opt("psdestimatestart","["+psdEstimateStart.strip(",")+"]")
      self.add_var_opt("psdestimateend","["+psdEstimateEnd.strip(",")+"]")

      self.add_var_opt("importanceresample",10000)

      self.id = followupmcmcJob.name + '-' + self.ifonames + '-' + str(trig.statValue) + '_' + str(trig.eventID) + '_' + randomseed
      outputName = followupmcmcJob.name+'/'+self.id
      self.outputCache = self.ifonames + ' ' + followupmcmcJob.name + ' ' + self.id.split('-')[-1] + ' ' + outputName + '.csv\n'

      self.setupNodeWeb(followupmcmcJob,False,None,None,None,dag.cache)
      self.add_var_opt("outfilename",outputName)

      if not opts.disable_dag_categories:
        self.set_category(followupmcmcJob.name.lower())

      if opts.mcmc:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add followupmcmc job for " + self.ifonames + "@ "+ str(trig.gpsTime[maxIFO])

##############################################################################

class followupspinmcmcJob(pipeline.CondorDAGJob, webTheJob):
  """
  A spinning MCMC job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "followupspinspiral":"SPINspiral"
      }
    }
  def __init__(self, options, cp, tag_base='FOLLOWUPSPINMCMC'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'followupspinmcmcJob'
    self.__executable = string.strip(cp.get('condor','followupspinspiral'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__prog__,tag_base)


###############################################################################

class followupspinmcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of spinning mcmc followup job
  """
  def __init__(self, followupspinmcmcJob, procParams, trig, cp, opts, dag, chain_number):

    try:

      self.friendlyName = 'SPIN MCMC followup'
      pipeline.CondorDAGNode.__init__(self,followupspinmcmcJob)

      IFOs = trig.ifolist_in_coinc

      maxSNR = 0
      maxIFO = ""
      for trigger in trig.coincs:
        snr = trigger.snr
        if snr > maxSNR:
          maxSNR = snr
          maxIFO = trigger.ifo
      trig_tempo = getattr(trig.coincs,maxIFO)
      triggerRef = copy.deepcopy(trig_tempo)
      self.ifoRef = maxIFO

      #self.add_var_opt("gps-time", str(trig.gpsTime[maxIFO]))
      self.add_var_arg(string.strip(cp.get("followup-spin-mcmc","input_file")))

      self.id = followupspinmcmcJob.name + '-' + trig.ifos + '-' + str(trig.statValue) + '_' + str(trig.eventID) + '_' + str(chain_number)

      self.setupNodeWeb(followupspinmcmcJob,False,None,None,None,dag.cache)

      if not opts.disable_dag_categories:
        self.set_category(followupspinmcmcJob.name.lower())

      if opts.spin_mcmc:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()

    except:
      self.invalidate()
      print "couldn't add followupspinmcmc job for " + trig.ifos + "@ "+ str(trig.gpsTime[maxIFO])

###############################################################################

class plotmcmcJob(pipeline.CondorDAGJob, webTheJob):
  """
  A plot mcmc job
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "plotmcmc":"plotmcmc.py"
      }
    }
  def __init__(self, options, cp, tag_base='PLOTMCMC'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'plotmcmcJob'
    self.__executable = string.strip(cp.get('condor','plotmcmc'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.setupJobWeb(self.__prog__,tag_base)

###############################################################################

class plotmcmcNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of  plotmcmc job
  """
  def __init__(self, plotmcmcjob, trig, mcmcIdList, cp, opts, dag, ifo, ifonames):

    try:
      self.friendlyName = 'plot MCMC'
      pipeline.CondorDAGNode.__init__(self,plotmcmcjob)

      if cp.has_option('followup-plotmcmc','burnin'):
        burnin = string.strip(cp.get('followup-plotmcmc','burnin'))
        if burnin.strip():
          self.add_var_opt("burnin",burnin)

      plot_routine = string.strip(cp.get('followup-plotmcmc','plot_routine'))
      executable = string.strip(cp.get('followup-plotmcmc','executable'))
      sim = None
      try:
        sim = isinstance(trig.coincs.sim,lsctables.SimInspiral)
      except:
        pass
      if sim:
        time = eval("trig.coincs.sim." + ifo[0:1].lower() + "_end_time")
        time_ns = eval("trig.coincs.sim." + ifo[0:1].lower() + "_end_time_ns")
        gps = float(time) + float(time_ns)/1000000000.
        mchirp = trig.coincs.sim.mchirp
        eta = trig.coincs.sim.eta
        distance = trig.coincs.sim.distance
        phi = trig.coincs.sim.phi0
      else:
        gps = trig.gpsTime[ifo]
        mchirp = getattr(trig.coincs,ifo).mchirp
        eta = getattr(trig.coincs,ifo).eta
        distance = getattr(trig.coincs,ifo).eff_distance
        phi = "0.0"

      self.add_var_opt("plot-routine",plot_routine)
      self.add_var_opt("executable",executable)
      self.add_var_opt("reference-time",gps)
      self.add_var_opt("reference-mchirp",mchirp)
      self.add_var_opt("reference-eta",eta)
      self.add_var_opt("reference-distance",distance)
      self.add_var_opt("reference-phi",phi)
      # get the list of MCMC .txt files to be used as input
      mcmcfilelist = ""
      for mcmcId in mcmcIdList:
        mcmcfilelist += mcmcId.split('-')[0]+'/' + mcmcId + '.csv,' # here we assume that the directory name is mcmcId.split('-')[0]
      self.add_var_opt("mcmc-file",mcmcfilelist.strip(','))

      self.id = plotmcmcjob.name + '-' + ifonames + '-' + str(trig.statValue) + '_' + str(trig.eventID)
      self.add_var_opt("identity",self.id)

      if cp.has_option('followup-plotmcmc', 'output') and string.strip(cp.get('followup-plotmcmc', 'output')):
        outputpath = string.strip(cp.get('followup-plotmcmc', 'output'))
      else:
        #outputpath = dag.publish_path + '/' + plotmcmcjob.name
        outputpath = plotmcmcjob.name
      if not os.access(outputpath,os.F_OK):
        os.makedirs(outputpath)
      else:
        if not os.access(outputpath,os.W_OK):
          print >> sys.stderr, 'path '+outputpath+' is not writable'
          sys.exit(1)

      if cp.has_option('followup-plotmcmc','web') and string.strip(cp.get('followup-plotmcmc','web')):
        webpath = string.strip(cp.get('followup-plotmcmc','web'))
      else:
        #webpath = dag.page + '/' + plotmcmcjob.name
        webpath = plotmcmcjob.name

      output_page = webpath + '/' + self.id
      self.outputCache = self.id.replace('-',' ') + " " + os.path.abspath(outputpath) + "/" + self.id + "\n"
      self.setupNodeWeb(plotmcmcjob,False,dag.webPage.lastSection.lastSub,None,output_page,dag.cache)

      self.add_var_opt("output-path",outputpath)

      if not opts.disable_dag_categories:
        self.set_category(plotmcmcjob.name.lower())

      # only add a parent if it exists
      for node in dag.get_nodes():
        if isinstance(node,followupmcmcNode):
          if not node.id.find(ifonames + '-' + str(trig.statValue) + '_' + str(trig.eventID)) == -1:
            try:
              if node.validNode: self.add_parent(node)
            except: pass

      if opts.plot_mcmc:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else:
        self.invalidate()

    except:
      self.invalidate()
      print "couldn't add plot mcmc job for " + ifonames + "@ "+ str(trig.gpsTime[ifo])


#### A CLASS TO DO FOLLOWUP COHERENT INSPIRAL JOBS ############################
###############################################################################
class followUpChiaJob(inspiral.ChiaJob,webTheJob):
  """
  Generates coherent inspiral data
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "chia":"lalapps_coherent_inspiral"
      }
    }
  def __init__(self, options, cp, tag_base='CHIA'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'followUpChiaJob'
    self.__executable = string.strip(cp.get('condor','chia'))
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self._InspiralAnalysisNode__pad_data = 0
    self.setupJobWeb(self.__prog__,tag_base)

##############################################################################

class followUpChiaNode(inspiral.ChiaNode,webTheNode):
  """
  A C code for computing the coherent inspiral statistic.
  An example command line is:
lalapps_coherent_inspiral --segment-length 1048576 --dynamic-range-exponent 6.900000e+01 --low-frequency-cutoff 4.000000e+01 --bank-file H1H2-COHBANK_COHERENT_H1H2_PLAYGROUND-823269333-600.xml --sample-rate 4096 --cohsnr-threshold 5.500000e+00 --ifo-tag H1H2 --frame-type LSC-STRAIN --H1-framefile H1-INSPIRAL_COHERENT_H1H2_PLAYGROUND-823269286-2048.gwf --H2-framefile H2-INSPIRAL_COHERENT_H1H2_PLAYGROUND-823268952-2048.gwf --gps-end-time 823269933 --gps-start-time 823269333 --write-cohsnr --write-cohnullstat --write-cohphasediff --write-events --verbose
  """

  #def __init__(self, chiaJob, procParams, trig, cp,opts,dag, trig_node, notrig_node ):

  def __init__(self,job,trig,opts,dag,cp):

    if 1:#try:
      # the use of this class would require some reorganisation in fu_Condor.py
      # and webCondor.py in order to set up the jobs following the same scheme
      # as the way it is done for the Inspiral pipeline...
      pipeline.CondorDAGNode.__init__(self,job)
      self.friendlyName = 'Produce coherent inspiral plots of event'
      self.id = job.name + '-CHIA-' + str(trig.statValue) + '_' + str(trig.eventID)
      #self.setupNodeWeb(job)
      self.output_file_name = ""
      self.add_var_opt("segment-length",string.strip(cp.get('chia','segment-length')))
      self.add_var_opt("dynamic-range-exponent",string.strip(cp.get('chia','dynamic-range-exponent')))
      self.add_var_opt("low-frequency-cutoff",string.strip(cp.get('chia','low-frequency-cutoff')))
      self.add_var_opt("sample-rate",string.strip(cp.get('chia','sample-rate')))
      self.add_var_opt("cohsnr-threshold",string.strip(cp.get('chia','cohsnr-threshold')))
      self.add_var_opt("ra-step",string.strip(cp.get('chia','ra-step')))
      self.add_var_opt("dec-step",string.strip(cp.get('chia','dec-step')))
      self.add_var_opt("numCohTrigs",string.strip(cp.get('chia','numCohTrigs')))
      self.add_var_opt("user-tag",str(trig.eventID))
      self.add_var_opt("ifo-tag",trig.ifoTag)
      # required by followUpChiaPlotNode
      bankFile = 'trigTemplateBank/' + trig.ifoTag + '-COHBANK_FOLLOWUP_' + str(trig.eventID) + '-' + str(int(trig.gpsTime[trig.ifolist_in_coinc[0]])) + '-2048.xml.gz'

      self.add_var_opt("write-events","")
      self.add_var_opt("write-compress","")
      self.add_var_opt("write-cohsnr","")
      self.add_var_opt("write-cohnullstat","")
      self.add_var_opt("write-h1h2nullstat","")
      self.add_var_opt("write-cohh1h2snr","")
      #self.add_var_opt("verbose","")
      self.set_bank(bankFile)
      self._InspiralAnalysisNode__pad_data = 0

      #CHECK: needed here? self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
      self.setupNodeWeb(job,False,None,None,None,dag.cache)
      self.add_var_opt("output-path",job.outputPath)

      # Here we define the trig-start-time and the trig-end-time;
      # The difference between these two times should be kept to 2s
      # Otherwise change the clustering window also
      hLengthAnalyzed = 1
      self.add_var_opt("gps-start-time",int(trig.gpsTime[trig.ifolist_in_coinc[0]]) - int(hLengthAnalyzed) )
      self.add_var_opt("gps-end-time",int(trig.gpsTime[trig.ifolist_in_coinc[0]]) + int(hLengthAnalyzed) )
      skipParams = ['minimal-match', 'bank-file', 'injection-file', 'trig-start-time', 'trig-end-time']

      if opts.plot_chia:
        dag.addNode(self,'chia')
        self.validate()
      else: self.invalidate()

    else: #except:
      print >> sys.stderr, "Didn't find a coherent inspiral job, I'll assume I don't need it"
      # initialize the extension of the output file. If the option
      # write_compress is found in procParams the extension will be overwritten

  def append_insp_node(self,inspNode,ifo):
    fileName = str(inspNode.output_file_name)
    #self.add_var_opt("H1-xml-file",str(fileName))
    #self.add_var_opt(ifo+"-framefile",str(fileName.replace(".xml",".gwf").strip(".gz")))
    self.add_var_arg("--"+ifo+"-framefile "+str(fileName.replace(".xml",".gwf").strip(".gz")))
    #self.add_var_opt(ifo.lower()+"-xml-file",str(fileName))
    if inspNode.validNode: self.add_parent(inspNode)


  def add_node_to_dag(self,dag,opts,trig):
    if opts.plot_chia:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else:
      self.invalidate()
      print "couldn't add coherent-inspiral job for " + str(trig.eventID)

########## CLUSTER coherent triggers ##########################################
###############################################################################
    
##############################################################################
# jobs class for clustering coherent triggers

class followUpCohireJob(pipeline.CondorDAGJob,webTheJob):
  """
  A clustering job for coherent inspiral triggers
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "cohire":"lalapps_cohire"
      }
    }
  def __init__(self, options, cp, tag_base='COHIRE'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'followUpCohireJob'
    self.__executable = string.strip(cp.get('condor','cohire'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)


##############################################################################
# node class for clustering coherent triggers 

class followUpCohireNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a cohire (coherent trigger clustering) job
  """
  def __init__(self,job,chiaXmlFilePath,trig,chiaNode,dag,page,opts):
    """
    job = A CondorDAGJob that can run an instance of COHIRE clustering
    """
    self.friendlyName = 'Produce xml file of clustered coherent triggers'
    if 1: #try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.output_file_name = ""
      inputFileName = 'trigTemplateBank/' + trig.ifoTag + '-COHIRE_FOLLOWUP_' + str(trig.eventID) + '-' + str(int(trig.gpsTime[trig.ifolist_in_coinc[0]]-1)) + '-2.txt'
      outputXmlFile = chiaXmlFilePath + trig.ifoTag + '-COHIRE_FOLLOWUP_' + str(trig.eventID) + '-' + str(int(trig.gpsTime[trig.ifolist_in_coinc[0]]-1)) + '-2.xml.gz'
      summaryFileName = chiaXmlFilePath + trig.ifoTag + '-COHIRE_SUMMARY_FOLLOWUP_' + str(trig.eventID) + '-' + str(int(trig.gpsTime[trig.ifolist_in_coinc[0]]-1)) + '-2.txt'
      self.add_var_opt("input",inputFileName)
      self.add_var_opt("data-type","all_data")
      self.add_var_opt("output",outputXmlFile)
      self.add_var_opt("summary-file",summaryFileName)
      self.add_var_opt("cluster-algorithm","snr")
      self.add_var_opt("sort-triggers","")
      self.add_var_opt("snr-threshold",5.0)
      self.add_var_opt("cluster-time",4000)
      self.id = job.name + '-' + str(trig.statValue) + '-' + str(trig.eventID)
      self.setupNodeWeb(job,False,None,None,None,dag.cache)
      skipParams = ['enable-output','output-path']

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      #try: 
      if chiaNode.validNode: self.add_parent(chiaNode)
      #except: pass
      if opts.plot_chia:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()
    #except: 
    #  self.invalidate()
    #  print "couldn't add cohire clustering job for event " + str(trig.eventID)

###############################################################################

########## PLOT Coherent search and null statistics TIME SERIES ###############
###############################################################################

##############################################################################
# jobs class for plotting coherent inspiral search and null stat timeseries

class plotChiaJob(pipeline.CondorDAGJob,webTheJob):
  """
  A followup plotting job for coherent inspiral search and null stat timeseries
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "plotchiatimeseries":"plotchiatimeseries"
      }
    }

  def __init__(self, options, cp, tag_base='PLOT_CHIA'):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'plotChiaJob'
    self.__executable = string.strip(cp.get('condor','plotchiatimeseries'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__,tag_base)

##############################################################################
# node class for plotting coherent inspiral search and null stat timeseries

class plotChiaNode(pipeline.CondorDAGNode,webTheNode):
  """
  Runs an instance of a plotChia followup job
  """
  def __init__(self,job,chiaXmlFilePath,trig,cohireNode,dag,page,opts,cp):
    """
    job = A CondorDAGJob that can run an instance of plotChiaJob followup.
    """
    self.friendlyName = 'Plot CHIA time-series'
    try:
      pipeline.CondorDAGNode.__init__(self,job)
      self.output_file_name = ""
      chiaXmlFileName = chiaXmlFilePath + trig.ifoTag + '-COHIRE_FOLLOWUP_' + str(trig.eventID) + '-' + str(int(trig.gpsTime[trig.ifolist_in_coinc[0]]-1)) + '-2.xml.gz'
      self.add_var_opt("chiaXmlFile",chiaXmlFileName)
      self.add_var_opt("gps-start-time",int(trig.gpsTime[trig.ifolist_in_coinc[0]]-1))
      self.add_var_opt("gps-end-time",int(trig.gpsTime[trig.ifolist_in_coinc[0]]+1))
      self.add_var_opt("sample-rate",string.strip(cp.get('chia','sample-rate')))
      self.add_var_opt("user-tag",str(trig.eventID))
      self.add_var_opt("ifo-tag",trig.ifoTag)
      self.add_var_opt("ifo-times",trig.ifoTag)
      self.id = job.name + '-' + str(trig.statValue) + '-' + str(trig.eventID)
      self.setupNodeWeb(job,True,None,None,None,dag.cache)

      if not opts.disable_dag_categories:
        self.set_category(job.name.lower())

      #try: 
      if cohireNode.validNode: self.add_parent(cohireNode)
      #except: pass
      if opts.plot_chia:
        dag.addNode(self,self.friendlyName)
        self.validate()
      else: self.invalidate()
    except: 
      self.invalidate()
      print "couldn't add chia plotting job for event " + str(trig.eventID)

  def append_insp_node(self,inspNode,ifo):
    fileName = str(inspNode.output_file_name)
    self.add_var_arg("--"+ifo+"-framefile "+str(fileName.replace(".xml",".gwf").strip(".gz")))
    if inspNode.validNode: self.add_parent(inspNode)

###############################################################################

############# Generate the checklists for candidate followups #################
###############################################################################

##############################################################################
# job class for making candidate checklists

class makeCheckListJob(pipeline.CondorDAGJob,webTheJob):
  """
  A job to prepare the checklist of a candidate
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "makechecklist":"makeCheckList.py"
      }
    }
  def __init__(self, options, cp):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'CHECKLIST'
    self.__executable = string.strip(cp.get('condor','makechecklist'))
    self.__universe = "local"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__)

##############################################################################
# node class for making candidate checklists

class makeCheckListNode(pipeline.CondorDAGNode,webTheNode):
  """
  A node to prepare the checklist of a candidate
  """
  def __init__(self,job,trig,cp,opts,dag):
    """
    """
    self.friendlyName = 'Make checklist'
    pipeline.CondorDAGNode.__init__(self,job)

    self.id = job.name + "-" + str(trig.eventID)

    self.add_var_opt("trigger-id",str(trig.eventID))
    gpsList = ""
    ifolist = ""
    for ifo in trig.ifolist_in_coinc:
      ifolist += ifo
      gpsList += repr(trig.gpsTime[ifo]) + ","
    self.add_var_opt("trigger-gps",gpsList.strip(","))
    self.add_var_opt("ifolist-in-coinc",ifolist)
    self.add_var_opt("user-tag",str(trig.eventID)+"_"+ifolist+"_"+str(int(trig.gpsTime[trig.ifolist_in_coinc[0]])))
    self.add_var_opt("ifo-times",trig.ifoTag)
    self.add_var_opt("ifo-tag",trig.ifoTag)
    if cp.has_option("followup-dq","input-sql"):
      self.add_var_opt("data-quality-database",cp.get("followup-dq","input-sql"))
    elif cp.has_option("followup-dq","server-url"):
      self.add_var_opt("segment-url",cp.get("followup-dq","server-url"))
    if cp.has_option("followup-ratiotest","input-pickle"):
      self.add_var_opt("SNR-ratio-test",cp.get("followup-ratiotest","input-pickle"))
    if cp.has_section("followup-analyse-qscan"):
      if cp.has_option("followup-analyse-qscan","hoft-qscan-ref-channel"):
        self.add_var_opt("hoft-channel-ref",cp.get("followup-analyse-qscan","hoft-qscan-ref-channel"))

    if cp.has_option("followup-foreground-qscan","remote-ifo"):
      remote_ifo = cp.get("followup-foreground-qscan","remote-ifo")
      self.add_var_opt("remote-qscan-web",remote_ifo+","+string.strip(cp.get("followup-foreground-qscan",remote_ifo+"web")))
    if cp.has_option("followup-foreground-seismic-qscan","remote-ifo"):
      remote_ifo = cp.get("followup-foreground-seismic-qscan","remote-ifo")
      self.add_var_opt("remote-seismic-qscan-web",remote_ifo+","+string.strip(cp.get("followup-foreground-seismic-qscan",remote_ifo+"web")))

    self.setupNodeWeb(job,True,None,None,None,dag.cache)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    for node in dag.get_nodes():
      if isinstance(node,IFOstatus_checkNode) or isinstance(node,FrCheckNode) or isinstance(node,plotSNRCHISQNode) or isinstance(node,pylal_skyPlotNode) or isinstance(node,plotChiaNode) or isinstance(node,plotmcmcNode) or isinstance(node,followupTriggerNode):
        if str(trig.eventID) in node.id and node.validNode:
          self.add_parent(node)
    
    if opts.make_checklist:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()


##############################################################################
# job class for plotting triggers in the chunk versus time

class followupTriggerJob(pipeline.CondorDAGJob,webTheJob):
  """
  A job to plot the triggers in the chunk
  """
  defaults={
    "section":"condor",
    "options":{
      "universe":"vanilla",
      "fu_triggers":"fup_triggers.py"
      }
    }
  def __init__(self, options, cp):
    """
    """
    if not(verifyCP(cp,self.defaults)):
      modifyCP(cp,self.defaults)
    self.__prog__ = 'followUpTriggers'
    self.__executable = string.strip(cp.get('condor','fu_triggers'))
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('getenv','True')
    self.setupJobWeb(self.__prog__)

##############################################################################
# node class for plotting triggers in the chunk versus time

class followupTriggerNode(pipeline.CondorDAGNode,webTheNode):
  """
  A node to plot triggers in the chunk
  """
  def __init__(self,job,trig,cp,opts,dag):
    """
    """
    self.friendlyName = 'plot triggers'
    pipeline.CondorDAGNode.__init__(self,job)

    self.id = job.name + "-" + str(trig.eventID)

    if opts.convert_eventid:
      self.add_var_opt("old-document",True)

    if opts.generate_fu_cache or not cp.has_option('followup-triggers','hipe-output-cache'):
      cacheString = 'fu_hipe.cache'
    else:
      cacheString = string.strip(cp.get('followup-triggers','hipe-output-cache'))

    followupTag = string.strip(cp.get("followup-triggersInChunk","tag"))

    if cp.has_option("followup-triggersInChunk","exttrig"):
      self.add_var_opt("followup-exttrig",True)

    if cp.has_option("followup-triggersInChunk","sned"):
      self.add_var_opt("followup-sned",string.strip(cp.get("followup-triggersInChunk","sned")))

    self.add_var_opt("gps-time",float(trig.gpsTime[trig.ifolist_in_coinc[0]]))
    self.add_var_opt("ifo-times",trig.ifoTag)
    self.add_var_opt("followup-tag",followupTag)
    self.add_var_opt("windows",string.strip(cp.get('followup-triggersInChunk','windows')))
    self.add_var_opt("event-id",str(trig.eventID))
    self.add_var_opt("cache-file",cacheString)
    self.setupNodeWeb(job,True,None,None,None,dag.cache)

    if not opts.disable_dag_categories:
      self.set_category(job.name.lower())

    if opts.followup_triggers:
      dag.addNode(self,self.friendlyName)
      self.validate()
    else: self.invalidate()

