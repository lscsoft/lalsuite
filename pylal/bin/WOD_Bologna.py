#!/usr/bin/python

from optparse import OptionParser

try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3
import sys, os, copy, math, re, subprocess, string, tempfile, socket
import time as time_method
import ConfigParser

from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue import pipeline
from glue.ligolw import ilwd
from glue import lal
from pylal import db_thinca_rings
from pylal import stfu_pipe
from pylal import fu_utils
from lalapps import inspiral

#from pylal.fu_Condor import *
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
dbtables.lsctables.LIGOTimeGPS = LIGOTimeGPS

###############################################################################
###### UTILITY FUNCTIONS ######################################################
###############################################################################

def write_user_info(cp):

	outputdir = cp.get('fu-output','output-dir').strip()
	print >>sys.stdout, "\n********* Information for user **********\n"
	print >>sys.stdout, "Before submitting your dag, update your grid-proxy by running:\n$ grid-proxy-init\n"
	print >>sys.stdout, "Then submit the dag to condor by running:\n$ condor_submit_dag WOD_Bologna.dag\n"
	print >>sys.stdout, "To monitor the progress of your dag please run:\n$ tail -f WOD_Bologna.dag.dagman.out \n"
	print >>sys.stdout, "The results will be written under the directory:\n %s" % (outputdir)

def exportGPSEventToDisk(tevent, dir, cnt, dag, filename=None):
        """
        """
        #If directy coincHeadings not there make the directory if
        #filename is None
        headingDir= dir + "/coinc_info"
        ifos = tevent.instruments
        instruments =tevent.instruments
        time = tevent.time

        idKeyStr="%s_%s" % (str(time), instruments)
        if filename==None:
                filename="coincEvent.info"
                stfu_pipe.mkdir(headingDir)
                filename=os.path.normpath(headingDir+'/'+idKeyStr+'_'+ filename)
        fp=open(filename,'w')
        fp.write("#DIR\t\t\tRANK\tFAR\t\tSNR\tIFOS\tINSTRUMENTS\tTIME\t\tMASS\n")
        fp.write("%-16s\t%d\t%0.2e\t%.2f\t%s\t%s\t\t%.3f\t%.2f\n" % (dir, cnt, 0, 0, ifos, instruments, float(time), 0) )
        fp.write("#DIR\t\t\tIFO\tTIME\t\tSNR\tCHISQ\tMASS1\tMASS2\n")
        rowString="%-16s\t%s\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\n"
        content=list()
        for ifo in tevent.ifos_list:
                content.append(rowString%(dir,
                                          ifo,
                                          float(time),
                                          float(0),
                                          float(0),
                                          float(0),
                                          float(0)))
        cache = lal.CacheEntry(instruments, "COINC_INFO_"+dir.upper(), segments.segment(float(time), float(time)), "file://localhost/"+os.path.abspath(filename))
        dag.output_cache.append(cache)
        fp.writelines(content)
        fp.close()
        return os.path.split(filename)[1]

class time_only_event(object):
	def __init__(self, ifostr, time):
		self.ifos_list = dbtables.lsctables.instrument_set_from_ifos(ifostr)
		self.ifos = ifostr
		self.instruments = ifostr
		self.time = float(time)

class time_only_events(object):
	def __init__(self, timestr):
		self.events = []
		if not timestr: return
		values = timestr.split(',')
		for value in values:
			ifos = value.split(':')[0]
			time = value.split(':')[1]
			self.events.append(time_only_event(ifos,time))

class extractTimesFromFile(object):
	def __init__(self, inputfile):
		self.events = []
		trigger_list = fu_utils.listFromFile(inputfile)
		if not trigger_list: return
		for value in trigger_list:
			ifos = value.split(' ')[0]
			time = value.split(' ')[1]
			self.events.append(time_only_event(ifos,time))

class create_default_config_wod(object):
	def __init__(self, configfile=None):
		cp = ConfigParser.ConfigParser()
		self.cp = cp
		self.time_now = "_".join([str(i) for i in time_method.gmtime()[0:6]])
		self.ini_file=self.time_now + ".ini"
		home_base = stfu_pipe.home_dirs()

		# CONDOR SECTION NEEDED BY THINGS IN INSPIRAL.PY
                cp.add_section("condor")
		cp.set("condor","datafind",self.which("ligo_data_find"))
		cp.set("condor","universe","standard")
		# SECTIONS TO SHUT UP WARNINGS
                #cp.add_section("data")

		# DATAFIND SECTION
		cp.add_section("datafind")

		# FU-CONDOR SECTION
		cp.add_section("fu-condor")
		cp.set("fu-condor","datafind",self.which("ligo_data_find"))
		cp.set("fu-condor","convertcache",self.which("convertlalcache.pl"))
		self.set_qscan_executable()

		# fu-q-hoft-datafind SECTION
		cp.add_section("fu-q-hoft-datafind")
		for ifo in ["V1"]:
		#for ifo in ["H1","H2","L1","V1"]:
			cp.set("fu-q-hoft-datafind",ifo+"-search-time-range","128")
		# fu-q-rds-datafind SECTION
		cp.add_section("fu-q-rds-datafind")
		#for ifo in ["H1","H2","L1"]:
		#	cp.set("fu-q-rds-datafind",ifo+"-search-time-range","1024")
		cp.set("fu-q-rds-datafind","V1-search-time-range","2048")

		# fu-fg-ht-qscan SECTION
		cp.add_section("fu-fg-ht-qscan")
		for config in ["V1config"]:
		#for config in ["H1config","H2config","L1config","V1config"]:
			cp.set("fu-fg-ht-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'hoft') + ".txt","QSCAN CONFIG"))

		# fu-fg-rds-qscan SECTION
		cp.add_section("fu-fg-rds-qscan")
		#for config in ["H1config","H2config","L1config"]:
		#	cp.set("fu-fg-rds-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'rds') + ".txt","QSCAN CONFIG"))
		cp.set("fu-fg-rds-qscan","V1config","/storage/gpfs_virgo3/virgo/omega/configurations/s6_foreground_V1-raw-cbc.txt")

		# fu-fg-seismic-qscan SECTION
		cp.add_section("fu-fg-seismic-qscan")
		#for config in ["H1config","H2config","L1config"]:
		#	cp.set("fu-fg-seismic-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'seismic') + ".txt","QSCAN CONFIG"))
		cp.set("fu-fg-seismic-qscan","V1config","/storage/gpfs_virgo3/virgo/omega/configurations/s6_foreground_V1-raw-seismic-cbc.txt")

		# REMOTE JOBS SECTION
		cp.add_section("fu-remote-jobs")
		remoteIfos,remoteJobs = self.get_remote_jobs()
		cp.set('fu-remote-jobs','remote-ifos',remoteIfos)
		cp.set('fu-remote-jobs','remote-jobs',remoteJobs)

		# FU-OUTPUT SECTION
		cp.add_section("fu-output")
		cp.set("fu-output","log-path",self.log_path())
		cp.set("fu-output","output-dir",self.web_dir())
		cp.set("fu-output","web-url", self.web_url())

		# CONDOR MAX JOBS SECTION
		cp.add_section("condor-max-jobs")

		# if we have an ini file override the options
		if configfile:
			user_cp = ConfigParser.ConfigParser()
			user_cp.read(configfile)
		else:
			# otherwise see if a file with the standard ini file exists in the directory, the user probably intends to use it
			try:
				user_cp = ConfigParser.ConfigParser()
				user_cp.read('WOD_Bologna.ini')
			except: pass
		# override the default options
		if user_cp: self.overwrite_config(user_cp,cp)

	def overwrite_config(self,config,cp):
		for section in config.sections():
			if not cp.has_section(section): cp.add_section(section)
			for option in config.options(section):
				cp.set(section,option,config.get(section,option))

	def log_path(self):
		host = stfu_pipe.get_hostname()
		#FIXME add more hosts as you need them
		if 'phy.syr.edu' in host: return '/usr1/' + os.environ['USER']
		if 'caltech.edu' in host: return '/usr1/' + os.environ['USER']
		if 'phys.uwm.edu' in host: return '/people/' + os.environ['USER']
		if 'aei.uni-hannover.de' in host: return '/local/user/' + os.environ['USER']

	def which(self,prog):
		which = subprocess.Popen(['which',prog], stdout=subprocess.PIPE)
		out = which.stdout.read().strip()
		if not out: print >>sys.stderr, "WARNING: could not find %s in your path, unless you have an ini file to overide the path to %s the DAG will fail" % (prog,prog)
		return out

	def web_dir(self):
		host = stfu_pipe.get_hostname()
		#FIXME add more hosts as you need them
		if 'caltech.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'phy.syr.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return os.path.abspath(os.environ['HOME']) + '/WWW/LSC/followups/' + self.time_now
		print sys.stderr, "WARNING: could not find web directory, returning empty string"
		return ''

	def web_url(self):
		host = stfu_pipe.get_hostname()
		#FIXME add more hosts as you need them
		if 'ligo.caltech.edu' in host: return "https://ldas-jobs.ligo.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-la.caltech.edu' in host: return "https://ldas-jobs.ligo-la.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-wa.caltech.edu' in host: return "https://ldas-jobs.ligo-wa.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return "https://ldas-jobs.phys.uwm.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'phy.syr.edu' in host: return "https://sugar-jobs.phy.syr.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return "https://atlas3.atlas.aei.uni-hannover.de/~" + os.environ['USER'] + '/LSC/followups/' + self.time_now
		print sys.stderr, "WARNING: could not find web server, returning empty string"
                return ''

	def __config_name(self,ifo,type):
		fileMap={
			"L1":{"hoft":"L1_hoft_cbc","rds":"L0L1-RDS_R_L1-cbc","seismic":"L0L1-RDS_R_L1-seismic-cbc"},
			"H1":{"hoft":"H1_hoft_cbc","rds":"H0H1-RDS_R_L1-cbc","seismic":"H0H1-RDS_R_L1-seismic-cbc"},
			"H2":{"hoft":"H2_hoft_cbc","rds":"H0H2-RDS_R_L1-cbc","seismic":"H0H2-RDS_R_L1-seismic-cbc"},
			"V1":{"hoft":"V1_hoft_cbc","rds":"V1-raw-cbc","seismic":"V1-raw-seismic-cbc"}
			}
		return fileMap[ifo][type]

	def __find_config(self,config,description):
		#FIXME why isn't there an environment variable for things in lalapps share?
		path = self.which('lalapps_inspiral')
		if path: path = os.path.split(path)[0]
		else:
			print >>sys.stderr, "COULD NOT FIND " + description + " FILE %s IN %s, ABORTING" % (config, path)
			raise ValueError
			sys.exit(1)
		out = path.replace('bin','share/lalapps') + '/' + config
		if not os.path.isfile(out):
			print >>sys.stderr, "COULD NOT FIND " + description + " FILE %s IN %s, ABORTING" % (config, out)
			raise ValueError
			sys.exit(1)
		return out

	def get_remote_jobs(self):
		host = stfu_pipe.get_hostname()
		#FIXME add more hosts as you need them
		if 'ligo.caltech.edu' or 'ligo-la.caltech.edu' or 'ligo-wa.caltech.edu' or 'phys.uwm.edu' or 'aei.uni-hannover.de' or 'phy.syr.edu' in host:
			remote_ifos = "V1"
			remote_jobs = "ligo_data_find_Q_RDS_full_data,wpipeline_FG_RDS_full_data,wpipeline_FG_SEIS_RDS_full_data,ligo_data_find_Q_RDS_playground,wpipeline_FG_RDS_playground,wpipeline_FG_SEIS_RDS_playground,ligo_data_find_Q_RDS_gps_only,wpipeline_FG_RDS_gps_only,wpipeline_FG_SEIS_RDS_gps_only,ligo_data_find_Q_RDS_time_slides,wpipeline_FG_RDS_time_slides,wpipeline_FG_SEIS_RDS_time_slides"
			return remote_ifos, remote_jobs
		return '', ''

	def set_qscan_executable(self):
		host = stfu_pipe.get_hostname()
		if 'phy.syr.edu' in host:
			self.cp.set("fu-condor","qscan",stfu_pipe.home_dirs()+"/rgouaty/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")
		else:
			self.cp.set("fu-condor","qscan",stfu_pipe.home_dirs()+"/romain/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")

	def get_cp(self):
		return self.cp

	def write(self):
		self.get_cp().write(open(self.ini_file,"w"))

def parse_command_line():
	parser = OptionParser(
		version = "%prog",
		description = "Pipeline to setup Remote Wscans On Demand"
	)
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")	
	parser.add_option("-f", "--config-file", default="WOD_Bologna.ini", help="the config file, default looks for stfu_pipe.ini in path, if none is found it makes one from your environment (only provide a config file if you know you must override something)")
	parser.add_option("-g", "--gps-times", default='', help="Specify gps times to follow up. Format --gps-times=ifos:time,ifos:time (e.g. --gps-times=H1L1:888888888.999,H2L1:787787787.987,H1H2L1:999999999.999). No segment validation is done. If there is no data for these times it will crash.")
	parser.add_option("-i", "--input-file", default='', help="Specify gps times to follow up inside a text file. Format --gps-times=myfile.txt. No segment validation is done. If there is no data for these times it will crash.")

	parser.add_option("","--disable-dag-categories",action="store_true",\
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

	parser.add_option("","--do-remoteScans", action="store_true",\
	default=True,help="enable the remote scans through condor flocking." \
	" This option should be deprecated as soon as the condor flocking is" \
	" set up on every LIGO cluster.")

	options, filenames = parser.parse_args()

	if not filenames: filenames = []

	return options, (filenames or [])

###############################################################################
##### MAIN ####################################################################
###############################################################################

# Parse options and config files
options, filenames = parse_command_line()
default_cp = create_default_config_wod(options.config_file)
cp = default_cp.get_cp()

# Initialize dag
dag = stfu_pipe.followUpDAG(options.config_file,cp,options)

# if the "do_remoteScans" option is valid, we need to prepare a job to check grid proxy path at the starting of the dag.
setup_proxy_job         = stfu_pipe.setupProxyJob(options,cp)
setup_proxy_node        = stfu_pipe.setupProxyNode(dag,setup_proxy_job,cp,options)

search='gps_only'
q_ht_data_find_job      = stfu_pipe.fuDataFindJob(cp, tag_base='qdatafind', dir=search)
q_rds_data_find_job     = stfu_pipe.fuDataFindJob(cp, tag_base='Q_RDS', dir=search)
ht_qscan_job            = stfu_pipe.qscanJob(options,cp, dir=search, tag_base='FG_HT')

remote_datafind_job     = stfu_pipe.remoteDatafindJob(options,cp, tag_base='Q_RDS', dir=search)
remote_rds_qscan_job    = stfu_pipe.remoteQscanJob(options,cp, dir=search, tag_base='FG_RDS')
remote_seis_qscan_job   = stfu_pipe.remoteQscanJob(options,cp, dir=search, tag_base='FG_SEIS_RDS')
distrib_remote_rds_qscan_job = stfu_pipe.distribRemoteQscanJob(options,cp, dir=search, tag_base='FG_RDS')
distrib_remote_seis_qscan_job = stfu_pipe.distribRemoteQscanJob(options,cp, dir=search, tag_base='FG_SEIS_RDS')

if options.gps_times:
	gpsevents = time_only_events(options.gps_times)
elif options.input_file:
	gpsevents = extractTimesFromFile(options.input_file)
else:
	print >> sys.stderr, "an argument is missing in the command:\n You need to use one of the options --gps-times or --input-file"
	sys.exit(1)

for cnt, event in enumerate(gpsevents.events):
	for ifo in event.ifos_list:
		if options.verbose:
			print >>sys.stdout, "following up %s @ %s" % (ifo, event.time)
		# h(t) QSCAN datafind Nodes
		ht_qscan_data_find_node = stfu_pipe.fuDataFindNode(dag, q_ht_data_find_job, cp, options, ifo, trigger_time=event.time, qscan=True)
		# RDS QSCAN datafind Nodes
		rds_qscan_data_find_node = stfu_pipe.fuDataFindNode(dag,q_rds_data_find_job, cp, options, ifo, trigger_time=event.time, qscan=True, data_type="rds")
		# h(t) QSCAN Nodes
		ht_qscan_node = stfu_pipe.fuQscanNode(dag, ht_qscan_job, cp, options, event.time, ifo, ht_qscan_data_find_node.outputFileName, p_nodes=[ht_qscan_data_find_node],type="ht")
		if cp.has_option('fu-remote-jobs','remote-ifos') and ifo in cp.get('fu-remote-jobs','remote-ifos'):
			remote_rds_qscan_node = stfu_pipe.fuRemoteQscanNode(dag, remote_rds_qscan_job, cp, options, event.time, ifo, p_nodes=[rds_qscan_data_find_node,setup_proxy_node], type="rds")
			remote_seis_qscan_node = stfu_pipe.fuRemoteQscanNode(dag, remote_seis_qscan_job, cp, options, event.time, ifo, p_nodes=[rds_qscan_data_find_node,setup_proxy_node], type="seismic")
			distrib_remote_rds_qscan_node = stfu_pipe.distribRemoteQscanNode(dag, distrib_remote_rds_qscan_job, cp, options, ifo, event.time, p_nodes=[remote_rds_qscan_node,rds_qscan_data_find_node], type="rds")
			distrib_remote_seis_qscan_node = stfu_pipe.distribRemoteQscanNode(dag, distrib_remote_seis_qscan_job, cp, options, ifo, event.time, p_nodes=[remote_seis_qscan_node,rds_qscan_data_find_node], type="seismic")

	# write the event info to disk and cache for stfu_page
	exportGPSEventToDisk(event, search, cnt, dag)

#### ALL FINNISH ####
default_cp.write()
dag.write_all()
write_user_info(cp)
