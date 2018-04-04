#!/usr/bin/env python

__prog__ = "wscan_background.py"
__title__ = "Generate bakground of omega scans"

##############################################################################

import os, sys, subprocess, string, socket, shutil, tarfile
from optparse import *
import ConfigParser
import time

from glue import pipeline
from glue import gpstime
from glue.ligolw import dbtables
from pylal import fu_utils
from pylal import git_version
from pylal import stfu_pipe

##############################################################################
# Useful methods

class create_default_config(object):
  def __init__(self,home_base):
    cp = ConfigParser.ConfigParser()

    cp.add_section("condor")
    cp.set("condor","datafind",self.which("ligo_data_find"))

    cp.add_section("fu-condor")
    cp.set("fu-condor","datafind",self.which("ligo_data_find"))
    cp.set("fu-condor","convertcache",self.which("convertlalcache.pl"))
    #cp.set("fu-condor","qscan",home_base+"/romain/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")
    cp.set("fu-condor","setuplogfile",self.which("wscan_bg_setup_log.py"))

    cp.add_section("datafind")

    cp.add_section("fu-q-rds-datafind")
    for ifo in ["H1","H2","L1"]:
      cp.set("fu-q-rds-datafind",ifo+"-search-time-range","1024")
    cp.set("fu-q-rds-datafind","V1-search-time-range","2048")

    cp.add_section("fu-q-hoft-datafind")
    for ifo in ["H1","H2","L1","V1"]:
      cp.set("fu-q-hoft-datafind",ifo+"-search-time-range","128")

    cp.add_section("followup-background-qscan-times")
    cp.set("followup-background-qscan-times","H1range","")
    cp.set("followup-background-qscan-times","L1range","")
    cp.set("followup-background-qscan-times","V1range","")
    cp.set("followup-background-qscan-times","segment-min-len","2048")
    cp.set("followup-background-qscan-times","segment-pading","64")
    cp.set("followup-background-qscan-times","random-seed","1")
    cp.set("followup-background-qscan-times","background-statistics","20")

    # fu-bg-ht-qscan SECTION
    cp.add_section("fu-bg-ht-qscan")
    for config in ["H1config","H2config","L1config","V1config"]:
      cp.set("fu-bg-ht-qscan",config,self.__qscan_config("s5_background_" + self.__config_name(config[:2],'hoft') + ".txt"))

    # fu-bg-rds-qscan SECTION
    cp.add_section("fu-bg-rds-qscan")
    for config in ["H1config","H2config","L1config","V1config"]:
      cp.set("fu-bg-rds-qscan",config,self.__qscan_config("s5_background_" + self.__config_name(config[:2],'rds') + ".txt"))

    # fu-bg-seismic-qscan SECTION
    cp.add_section("fu-bg-seismic-qscan")
    for config in ["H1config","H2config","L1config","V1config"]:
      cp.set("fu-bg-seismic-qscan",config,self.__qscan_config("s5_background_" + self.__config_name(config[:2],'seismic') + ".txt"))

    # OUTPUT SECTION
    cp.add_section("fu-output")
    cp.set("fu-output","log-path","/usr1/" + os.getenv("USER"))

    # REMOTE JOBS SECTION
    cp.add_section("fu-remote-jobs")
    remoteIfos,remoteJobs = self.get_remote_jobs()
    cp.set('fu-remote-jobs','remote-ifos',remoteIfos)
    cp.set('fu-remote-jobs','remote-jobs',remoteJobs)

    # CONDOR MAX JOBS SECTION
    cp.add_section("condor-max-jobs")
    cp.set("condor-max-jobs","ligo_data_find_Q_HT_","3")
    cp.set("condor-max-jobs","ligo_data_find_Q_RDS_","3")
    cp.set("condor-max-jobs","ligo_data_find_Q_HT_","3")
    cp.set("condor-max-jobs","ligo_data_find_Q_RDS_","3")
    cp.set("condor-max-jobs","wpipeline_BG_HT_","150")
    cp.set("condor-max-jobs","wpipeline_BG_RDS_","150")
    cp.set("condor-max-jobs","wpipeline_BG_SEIS_RDS_","150")

    self.cp = cp
    self.set_qscan_executable()

  def set_qscan_executable(self):
    host = self.__get_hostname()
    if 'phy.syr.edu' in host:
      self.cp.set("fu-condor","qscan",self.__home_dirs()+"/rgouaty/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")
    else:
      self.cp.set("fu-condor","qscan",self.__home_dirs()+"/romain/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")

  def __home_dirs(self):
    return os.path.split(os.environ['HOME'])[0]

  def __qscan_config(self,config):
    #FIXME why isn't there an environment variable for things in lalapps share?
    path = self.which('lalapps_inspiral')
    if path: path = os.path.split(path)[0]
    else:
      print >>sys.stderr, "COULD NOT FIND QSCAN CONFIG FILE %s IN %s, ABORTING" % (config, path)
      raise ValueError
      sys.exit(1)
    out = path.replace('bin','share/lalapps') + '/' + config
    if not os.path.isfile(out):
      print >>sys.stderr, "COULD NOT FIND QSCAN CONFIG FILE %s IN %s, ABORTING" % (config, out)
      raise ValueError
      sys.exit(1)
    return out

  def __config_name(self,ifo,type):
    fileMap={
            "L1":{"hoft":"L1_hoft_cbc","rds":"L0L1-RDS_R_L1-cbc","seismic":"L0L1-RDS_R_L1-seismic-cbc"},
            "H1":{"hoft":"H1_hoft_cbc","rds":"H0H1-RDS_R_L1-cbc","seismic":"H0H1-RDS_R_L1-seismic-cbc"},
            "H2":{"hoft":"H2_hoft_cbc","rds":"H0H2-RDS_R_L1-cbc","seismic":"H0H2-RDS_R_L1-seismic-cbc"},
            "V1":{"hoft":"V1_hoft_cbc","rds":"V1-raw-cbc","seismic":"V1-raw-seismic-cbc"}
            }
    return fileMap[ifo][type]

  def get_remote_jobs(self):
    host = self.__get_hostname()
    #FIXME add more hosts as you need them
    if 'ligo.caltech.edu' or 'ligo-la.caltech.edu' or 'ligo-wa.caltech.edu' or 'phys.uwm.edu' or 'aei.uni-hannover.de' or 'phy.syr.edu' in host:
      remote_ifos = "V1"
      remote_jobs = "ligo_data_find_Q_RDS_,wpipeline_BG_RDS_,wpipeline_BG_SEIS_RDS_"
      return remote_ifos, remote_jobs
    return '', ''

  def __get_hostname(self):
    host = socket.getfqdn()
    return host

  def which(self,prog):
    which = subprocess.Popen(['which',prog], stdout=subprocess.PIPE)
    out = which.stdout.read().strip()
    if not out: print >>sys.stderr, "WARNING: could not find %s in your path, unless you have an ini file to overide the path to %s the DAG will fail" % (prog,prog)
    return out

  def get_cp(self):
    return self.cp

def overwrite_config(cp,config):
  for section in config.sections():
    if not cp.has_section(section): cp.add_section(section)
    for option in config.options(section):
      cp.set(section,option,config.get(section,option))
  return cp

class setupLogFileJob(pipeline.CondorDAGJob,stfu_pipe.FUJob):

	def __init__(self,opts,cp):
		self.__executable = string.strip(cp.get('fu-condor','setuplogfile'))
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,cp=cp,dir='',tag_base='')

class setupLogFileNode(pipeline.CondorDAGNode,stfu_pipe.FUNode):

	def __init__(self,dag,job,cp,time_range,tag='start'):
		pipeline.CondorDAGNode.__init__(self,job)
		self.add_var_arg("--" + tag + "-run")
		self.add_var_opt("log-name",time_range+".log")
		outputString = "omega/" + stfu_pipe.science_run(int(time_range.split("_")[0])).upper() + "/background"
		if cp.has_option('fu-output','output-dir') and cp.get('fu-output','output-dir'):
			output = cp.get('fu-output','output-dir') + '/' + outputString
		else:
			output = outputString
		self.add_var_opt("output-path",output)

		if tag == 'terminate':
			for node in dag.get_nodes():
				if isinstance(node,stfu_pipe.fuQscanNode):
					self.add_parent(node)
		dag.add_node(self)
		self.validate()

##############################################################################
# METHOD TO PREPARE FILES FOR REMOTE VIRGO SCANS AT LYON CCIN2P3
##############################################################################

class prepareLyonRemoteScans(object):
  def __init__(self,cp,ifo,timeref):
    depIfoIniConfig = ifo+'config'
    self.depIfoDir = ifo+'_qscans_config'
    depQscanList = ['bg-rds-qscan', 'bg-seismic-qscan']

    # Build directories
    if not os.path.isdir(self.depIfoDir):
      os.makedirs(self.depIfoDir)
    if not os.path.isdir(self.depIfoDir+'/CONFIG'):
      os.makedirs(self.depIfoDir+'/CONFIG')
    if not os.path.isdir(self.depIfoDir+'/RESULTS'):
      os.makedirs(self.depIfoDir+'/RESULTS')
    for depQscan in depQscanList:
      if not os.path.isdir(self.depIfoDir+'/RESULTS/results_'+depQscan):
        os.makedirs(self.depIfoDir+'/RESULTS/results_'+depQscan)
    if not os.path.isdir(self.depIfoDir+'/SCRIPTS'):
      os.makedirs(self.depIfoDir+'/SCRIPTS')
    if not os.path.isdir(self.depIfoDir+'/TIMES'):
      os.makedirs(self.depIfoDir+'/TIMES')

    # Copy the qscan configuration files
    for depQscan in depQscanList:
      if cp.has_option("fu-"+depQscan, depIfoIniConfig):
        qscanConfig = self.fix_config_for_science_run( cp.get("fu-"+depQscan, depIfoIniConfig).strip(), timeref )
        if qscanConfig!='':
          print 'copy '+qscanConfig+' -----> '+self.depIfoDir+'/CONFIG/'+depQscan+'_config.txt'
          shutil.copy(qscanConfig,self.depIfoDir+'/CONFIG/'+depQscan+'_config.txt')

    # Copy the scripts used in the remote computing center
    #   first, get the path to the scripts
    scriptPath = subprocess.Popen(["which", "analyseQscan.py"], stdout=subprocess.PIPE).communicate()[0]

    scriptPath = scriptPath.strip('analyseQscan.py\n')
    shutil.copy(scriptPath+'/qsub_wscanlite.sh',self.depIfoDir+'/SCRIPTS/')
    shutil.copy(scriptPath+'/wscanlite_in2p3.sh',self.depIfoDir+'/SCRIPTS/')
    shutil.copy(scriptPath+'/prepare_sendback.py',self.depIfoDir)
    shutil.copy(scriptPath+'/virgo_qscan_in2p3.py',self.depIfoDir)
    os.chmod(self.depIfoDir+'/virgo_qscan_in2p3.py',0755)

  def fix_config_for_science_run(self, config, time):
    run = stfu_pipe.science_run(time)
    config_path = os.path.split(config)
    out = "/".join([config_path[0], config_path[1].replace('s5',run).replace('s6',run)])
    return out

##############################################################################
#MAIN PROGRAM
##############################################################################

home_dir = os.getenv("HOME")
home_base = "/".join(home_dir.split("/")[0:-1])

######################## OPTION PARSING  #####################################
usage = """usage: %prog [options]
"""

parser = OptionParser(usage, version=git_version.verbose_msg)

parser.add_option("-f","--config-file",action="store",type="string",\
    default="",help="configuration file is optional")

parser.add_option("-i","--ifos",action="store",type="string",\
    default="H1L1V1",help="list of requested ifos, expected format is of " \
    "the kind \"H1H2L1\" ")

parser.add_option("", "--disable-dag-categories",action="store_true",\
    default=False,help="disable the internal dag category maxjobs")

parser.add_option("","--no-ht-qscan", action="store_true",\
    default=False,help="disable hoft qscan nodes")

parser.add_option("","--no-rds-qscan", action="store_true",\
    default=False,help="disable rds qscan nodes")

parser.add_option("","--no-seismic-qscan", action="store_true",\
    default=False,help="disable seismic qscan nodes")

parser.add_option("","--no-htQscan-datafind", action="store_true",\
    default=False,help="disable hoft qscan datafind nodes")

parser.add_option("","--no-rdsQscan-datafind", action="store_true",\
    default=False,help="disable rds qscan datafind nodes")

parser.add_option("","--prepare-scan-ccin2p3", action="store_true",\
    default=False,help="prepare the scripts to submit Virgo omega scans at " \
    "Lyon CC")

command_line = sys.argv[1:]
(opts,args) = parser.parse_args()

#############################################################################

default_cp = create_default_config(home_base)
cp = default_cp.get_cp()
if opts.config_file: 
 config = ConfigParser.ConfigParser()
 config.read(opts.config_file)
 cp = overwrite_config(cp,config)

ifos_list = []
for j in range(0,len(opts.ifos)-1,2):
  ifo = opts.ifos[j:j+2]
  ifos_list.append(ifo)

#Get the start-end times of yesterday...
ifo_range = ",".join(stfu_pipe.get_day_boundaries(int(gpstime.GpsSecondsFromPyUTC(time.time())) - 86400))
# print "Start time : " + ifo_range.split(",")[0] + "   End Time : " + ifo_range.split(",")[-1]

range_string = ""
#Check the time ranges for each ifo in the ini file and , if they are left empty fill them with yesterday's start-end times.
for ifo_index,ifo in enumerate(ifos_list):
    if cp.has_option("followup-background-qscan-times",ifo+"range"):
      if not cp.get("followup-background-qscan-times",ifo+"range"):
        cp.set("followup-background-qscan-times",ifo+"range",ifo_range)
      range_string = string.strip(cp.get("followup-background-qscan-times",ifo+"range")).replace(',','_')

#Get current UTC time to be used in the ini file name
time_now = "_".join([str(i) for i in time.gmtime()[0:6]])

#Initialize dag
dag = stfu_pipe.followUpDAG(time_now + "-" + range_string + ".ini",cp,opts)

# CONDOR JOB CLASSES
htdataJob	= stfu_pipe.fuDataFindJob(cp,tag_base='Q_HT',dir='')
rdsdataJob	= stfu_pipe.fuDataFindJob(cp,tag_base='Q_RDS',dir='')
htQscanBgJob	= stfu_pipe.qscanJob(opts,cp,tag_base='BG_HT',dir='')
rdsQscanBgJob	= stfu_pipe.qscanJob(opts,cp,tag_base='BG_RDS',dir='')
seisQscanBgJob	= stfu_pipe.qscanJob(opts,cp,tag_base='BG_SEIS_RDS',dir='')
setupLogJob	= setupLogFileJob(opts,cp)

start_node = setupLogFileNode(dag,setupLogJob,cp,range_string,'start')

for ifo in ifos_list:

    # FIX ME: input argument segFile is not needed any more
    segFile = {}
    times, timeListFile = fu_utils.getQscanBackgroundTimes(cp,opts,ifo,segFile)

    #Prepare files for remote scans at Lyon CC...
    if opts.prepare_scan_ccin2p3 and ifo in cp.get("fu-remote-jobs","remote-ifos").strip().split(","):
        if times:
          CCRemoteScans = prepareLyonRemoteScans(cp,ifo,times[0])
        else:
          CCRemoteScans = None

    for qtime in times:
      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (REGULAR DATA SET)
      dNode = stfu_pipe.fuDataFindNode(dag,rdsdataJob,cp,opts,ifo,sngl=None,qscan=True,trigger_time=qtime,data_type='rds',p_nodes=[start_node])

      # SETUP DATAFIND JOBS FOR BACKGROUND QSCANS (HOFT)
      dHoftNode = stfu_pipe.fuDataFindNode(dag,htdataJob,cp,opts,ifo,sngl=None,qscan=True,trigger_time=qtime,p_nodes=[start_node])

      # SETUP BACKGROUND QSCAN JOBS
      qHtBgNode = stfu_pipe.fuQscanNode(dag,htQscanBgJob,cp,opts,qtime,ifo,dHoftNode.output_cache.path,p_nodes=[dHoftNode],type="ht",variety="bg")

      qRdsBgNode = stfu_pipe.fuQscanNode(dag,rdsQscanBgJob,cp,opts,qtime,ifo,dNode.output_cache.path,p_nodes=[dNode],type="rds",variety="bg")

      qSeisBgNode = stfu_pipe.fuQscanNode(dag,seisQscanBgJob,cp,opts,qtime,ifo,dNode.output_cache.path,p_nodes=[dNode],type="seismic",variety="bg")

      if opts.prepare_scan_ccin2p3 and ifo in cp.get("fu-remote-jobs","remote-ifos").strip().split(","):
        if hasattr(qRdsBgNode.output_cache,'__iter__'):
          dag.output_cache.extend(qRdsBgNode.output_cache)
        else:
          dag.output_cache.append(qRdsBgNode.output_cache)
        if hasattr(qSeisBgNode.output_cache,'__iter__'):
          dag.output_cache.extend(qSeisBgNode.output_cache)
        else:
          dag.output_cache.append(qSeisBgNode.output_cache)

    # WRITE TIMES FOR REMOTE (DEPORTED)CALCULATIONS
    if opts.prepare_scan_ccin2p3:
      if ifo in cp.get("fu-remote-jobs","remote-ifos").strip().split(",") and timeListFile and CCRemoteScans:
        shutil.copy(timeListFile,CCRemoteScans.depIfoDir+'/TIMES/background_qscan_times.txt')

end_node = setupLogFileNode(dag,setupLogJob,cp,range_string,'terminate')

#### ALL FINNISH ####

# Prepare for moving the deported qscan directory (tar-gzip)
if opts.prepare_scan_ccin2p3 and CCRemoteScans:

  tar = tarfile.open(CCRemoteScans.depIfoDir+'.tar.gz','w:gz')
  tar.add(CCRemoteScans.depIfoDir)
  tar.close()

cp.write(open(time_now + "-" + range_string + ".ini","w"))
dag.write_all()

