"""
Something


This program creates cache files for the output of inspiral hipe
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'

##############################################################################
# import standard modules and append the lalapps prefix to the python path
import sys, os, copy, math
import socket, time
import re, string
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict

##############################################################################
# import the modules we need to build the pipeline
from glue import pipeline
from glue import lal
from pylab import *
from glue import segments 
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from pylal import CoincInspiralUtils
from pylal.fu_utils import *
from pylal import Fr      
##############################################################################
# redefine the columns of interest
##############################################################################
lsctables.SimInspiralTable.loadcolumns = ["waveform","geocent_end_time",
    "geocent_end_time_ns","h_end_time","h_end_time_ns","l_end_time",
    "l_end_time_ns","source","mass1","mass2","mchirp","eta","distance",
    "spin1x","spin1y","spin1z","spin2x","spin2y","spin2z","eff_dist_h",
    "eff_dist_l","eff_dist_g","eff_dist_t","eff_dist_v"]

lsctables.SnglInspiralTable.loadcolumns = ["ifo","end_time","end_time_ns",
    "eff_distance","mass1","mass2","mchirp","eta","snr","chisq","chisq_dof",
    "sigmasq","event_id"]

################ TRIG BANK FROM SIRE FILE CONDOR DAG JOB ######################

class trigBankFollowUpJob(pipeline.CondorDAGJob):
  """
  A followup trig bank job
  """
  def __init__(self, options, tag_base='TRIGBANK_FOLLOWUP'):
    """
    """
    self.__executable = "lalapps_trigbank"
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
#    self.tag_base = tag_base
#    self.options = options
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/futrigbank-$(macrogpsstarttime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/futrigbank-$(macrogpsstarttime)-$(cluster)-$(process).err')
    self.set_sub_file('futrigbank.sub')


class trigBankFollowUpNode(pipeline.CondorDAGNode):
  """
  Runs an instance of a trig bank followup job
  """
  def __init__(self,job,time,proc,opts):
    """
    job = A CondorDAGJob that can run an instance of trigbank followup.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.output_file_name = ""
    self.output_user_tag = "FOLLOWUP_" + str(time)

    for row in proc:
      self.add_var_opt(row.param[2:],row.value)
      if row.param[2:] == 'gps-start-time':
        self.add_macro("macrogpsstarttime",row.value)
        self.start_time = row.value
      if row.param[2:] == 'input-ifo':
        self.input_ifo = row.value
      if row.param[2:] == 'ifo-tag':
        self.ifo_tag = row.value
      if row.param[2:] == 'gps-end-time':
        self.end_time = row.value

    self.add_var_opt("user-tag",self.output_user_tag)
    self.add_file_arg(opts.xml_glob)
    self.output_file_name = self.input_ifo + "-TRIGBANK_" + self.ifo_tag + \
      "_" + self.output_user_tag + "-" + self.start_time + "-" + \
      str(eval(self.end_time)-eval(self.start_time)) + ".xml"

class inspiralFollowUpJob(pipeline.CondorDAGJob):
  """
  A followup trig bank job
  """
  def __init__(self, options, tag_base='INSPIRAL_FOLLOWUP'):
    """
    """
    self.__executable = "lalapps_inspiral"
    self.__universe = "standard"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.tag_base = tag_base
    self.options = options
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")

    self.set_stdout_file('logs/fuinspiral-$(macrogpsstarttime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/fuinspiral-$(macrogpsstarttime)-$(cluster)-$(process).err')
    self.set_sub_file('fuinspiral.sub')


class inspiralFollowUpNode(pipeline.CondorDAGNode):
  """
  Runs an instance of a trig bank followup job
  """
  def __init__(self,job,time,proc,opts,bankFile):
    """
    job = A CondorDAGJob that can run an instance of trigbank followup.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.output_file_name = ""
    self.output_user_tag = "FOLLOWUP_" + str(time)
    self.add_rsq_veto = 1
    for row in proc:
      if row.param[2:] == 'bank-file':
        continue
      self.add_var_opt(row.param[2:],row.value)
      if row.param[2:] == 'gps-start-time':
        self.add_macro("macrogpsstarttime",row.value)
        self.start_time = row.value
      if row.param[2:] == 'channel-name':
        self.input_ifo = row.value[0:2]
      if row.param[2:] == 'ifo-tag':
        self.ifo_tag = row.value
      if row.param[2:] == 'gps-end-time':
        self.end_time = row.value
      if row.param[2:] == 'enable-rsq-veto':
        self.add_rsq_veto = 0
      if row.param[2:] == 'disable-rsq-veto':
        self.add_rsq_veto = 0
    self.add_var_opt("user-tag",self.output_user_tag)
    self.add_var_opt("write-snrsq","")
    self.add_var_opt("write-chisq","")
    self.add_var_opt("write-spectrum","")
    self.add_var_opt("bank-file",bankFile)
    if self.add_rsq_veto:
      self.add_var_opt("enable-rsq-veto","")
    self.output_file_name = self.input_ifo + "-INSPIRAL_" + self.ifo_tag + \
      "_" + self.output_user_tag + "-" + self.start_time + "-" +\
      str(eval(self.end_time)-eval(self.start_time)) + ".xml"
    
class plotSNRCHISQJob(pipeline.CondorDAGJob):
  """
  A followup trig bank job
  """
  def __init__(self, options, tag_base='PLOT_FOLLOWUP'):
    """
    """
    self.__executable = "plotsnrchisq_pipe.py"
    self.__universe = "vanilla"
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
    self.set_stdout_file('logs/fuplotsnrchisq-$(macrogpsstarttime)-$(cluster)-$(process).out')
    self.set_stderr_file('logs/fuplotsnrchisq-$(macrogpsstarttime)-$(cluster)-$(process).err')
    self.set_sub_file('fuplotsnrchisq.sub')


class plotSNRCHISQNode(pipeline.CondorDAGNode):
  """
  Runs an instance of a trig bank followup job
  """
  def __init__(self,job,time,fileName):
    """
    job = A CondorDAGJob that can run an instance of trigbank followup.
    """
    pipeline.CondorDAGNode.__init__(self,job)
    self.output_file_name = ""
    self.add_var_opt("frame-file",fileName.replace(".xml",".gwf"))
    self.add_var_opt("gps",time)
    self.add_var_opt("inspiral-xml-file",fileName)
    self.add_var_opt("output-path",".")

##############################################################################
#
#  MAIN PROGRAM
#
##############################################################################

###### OPTION PARSING AND SANITY CHECKS #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser( usage )

parser.add_option("-v", "--version",action="store_true",default=False,\
    help="print version information and exit")

parser.add_option("-c", "--cache-path",action="store",type="string",\
    metavar=" PATH",help="directory to find  all hipe XML files")

parser.add_option("-r", "--science-run", action="store",type="string",\
    metavar=" RUN", help="name of science run")

parser.add_option("-g","--xml-glob",action="store",type="string",\
    default=None, metavar=" XML_GLOB", \
    help="GLOB of coire xml files to read (Loudest events)" )

parser.add_option("-K","--statistic",action="store",type="string",\
    default="effective_snrsq",metavar=" STAT",\
    help="coincident statistic (default = effective_snr)")

parser.add_option("-a","--bitten-l-a",action="store",type="float",\
    default=3,metavar=" BLA",\
    help="bitten l a parameter")

parser.add_option("-b","--bitten-l-b",action="store",type="float",\
    default=3,metavar=" BLB",\
    help="bitten l b parameter")

parser.add_option("-p","--page",action="store",type="string",\
    default="investigations/s5/people/followups/",metavar=" PAGE",\
    help="web page path (default 'investigations/s5/people/followups/'")

parser.add_option("-n","--num-trigs",action="store",type="float",\
    default=1,metavar=" NUM_TRIGS",\
    help="number of triggers to follow up (default = 1)")

parser.add_option("-l", "--log-path",action="store",type="string",\
    metavar=" PATH",help="directory to write condor log file")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

if opts.version:
  print "$Id$"
  sys.exit(0)

if not opts.cache_path:
  print >> sys.stderr, "No cache path specified."
  print >> sys.stderr, "Use --cache-path PATH to specify a location."
  sys.exit(1)

if not opts.science_run:
  print >> sys.stderr, "No science run specified."
  print >> sys.stderr, "Use --science-run RUN to specify a run."
  sys.exit(1)

#if not opts.xml_glob:
#  print >> sys.stderr, "Must specify a GLOB of xmls to read"
#  sys.exit(1)

#if not opts.statistic:
#  print >> sys.stderr, "Must specify a statistic to use"
#  sys.exit(1)

print opts.log_path
############# TURN THE HIPE OUTPUT INTO LAL CACHE FILES #######################

cache = getCache(opts)
cache.getCacheAll()
cache.writeCacheAll()
print >> sys.stderr, "\nHIPE CACHE FILES WRITTEN TO:"
for n, t in cache.nameMaps:
  print >> sys.stderr, " * " + n + " [" + str(len(cache[t])) + "]"

##############################################################################
# create a log file that the Condor jobs will write to
basename = 'followup'
tempfile.tempdir = opts.log_path
tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

##############################################################################
# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file( basename )
subsuffix = '.sub'

############# READ IN THE COIRE FILES #########################################
if not opts.xml_glob:
  sys.exit(0)

found, coincs = readFiles(opts.xml_glob,getstatistic(opts))
missed = None
followuptrigs = getfollowuptrigs(opts,coincs,missed)

print "\n.......Found " + str(len(coincs)) + " trigs to follow up" 

trigJob = trigBankFollowUpJob(opts)
inspJob = inspiralFollowUpJob(opts)
plotJob = plotSNRCHISQJob(opts)

trigJobCnt = 0
inspJobCnt = 0
plotJobCnt = 0

print "\n.......Setting up pipeline jobs"
for trig in followuptrigs:
  try:
    trig_process_params = cache.getProcessParamsFromCache( \
                     cache.filesMatchingGPS(trig.gpsTime,'TRIGBANK'), \
                     trig.gpsTime)  
  except: "couldn't get trigbank process params for " + trig.gpsTime
  try: 
    inspiral_process_params = cache.getProcessParamsFromCache( \
                     cache.filesMatchingGPS(trig.gpsTime,'INSPIRAL_'), \
                     trig.gpsTime)
  except: "couldn't get trigbank process params for " + trig.gpsTime
  for ifo in trig_process_params:
    try: trig.gpsTime[ifo]
    except: continue
    if trig.gpsTime[ifo]:
      try:
        trigNode = trigBankFollowUpNode(trigJob,trig.gpsTime[ifo],trig_process_params[ifo],opts)
        dag.add_node(trigNode)
        trigJobCnt+=1
      except:
        trigNode = None
        print "couldn't add trigbank job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])
      if trigNode:
        try:
          inspiralNode = inspiralFollowUpNode(inspJob,trig.gpsTime[ifo],inspiral_process_params[ifo],opts,trigNode.output_file_name)
          inspiralNode.add_parent(trigNode)
          dag.add_node(inspiralNode)
          inspJobCnt+=1
        except: 
          inspiralNode = None
          print "couldn't add inspiral job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])
      if inspiralNode:
        try:
          plotNode = plotSNRCHISQNode(plotJob,trig.gpsTime[ifo],inspiralNode.output_file_name)
          plotNode.add_parent(inspiralNode)
          dag.add_node(plotNode)
          plotJobCnt+=1
        except:
          print "couldn't add plot job for " + str(ifo) + "@ "+ str(trig.gpsTime[ifo])


print "\nFound " + str(trigJobCnt)+ " trig bank jobs and "+ str(inspJobCnt) +" inspiral Jobs " \
      + str(plotJobCnt) + " plot jobs"

print "\n.......Writing DAG"
dag.write_sub_files()
dag.write_dag()

sys.exit(0)


#  for trig in followuptrigs:
#    print getattr(trig.coincs,'H1').snr

############# SET UP TRIG BANK JOBS   #########################################

############# SET UP INSPIRAL JOBS    #########################################

############# SET UP PLOTTING JOBS    #########################################



sys.exit(0)

#print coincs.getsngls('H1')

#for coinc in coincs:
#  print getattr(coinc, 'H1').snr

#    print cache.filesMatchingGPS(trig.gpsTime,'TRIGBANK')

