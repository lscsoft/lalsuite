"""
This module contains condor jobs / node classes for the SQlite Triggered Follow Up dag
"""

__author__ = 'Chad Hanna <channa@phys.lsu.edu>, Cristina Torres <cristina.torres@ligo.org>, Romain Gouaty <gouaty@lapp.in2p3.fr>'

try:
	import sqlite3
except ImportError:
	# pre 2.5.x
	from pysqlite2 import dbapi2 as sqlite3

import sys, os, copy, math, math, subprocess, socket, re, string
import time as time_method
from optparse import *
import tempfile
import ConfigParser
import urlparse
from UserDict import UserDict

sys.path.append('@PYTHONLIBDIR@')

from glue import segments
from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
#from glue.ligolw import dbtables
from glue.ligolw import utils
from glue import pipeline
from glue import lal
#from pylal import db_thinca_rings
from lalapps import inspiral
from pylal import date
from pylal.xlal import date as xlaldate
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
#dbtables.lsctables.LIGOTimeGPS = LIGOTimeGPS
lsctables.LIGOTimeGPS = LIGOTimeGPS

###############################################################################
##### UTILITY FUNCTIONS #######################################################
###############################################################################


def mkdir(output):
	# MAKE SURE WE CAN WRITE TO THE OUTPUT DIRECTORY
	if not os.access(output,os.F_OK): os.makedirs(output)
	else:
		if not os.access(output,os.W_OK):
			print >> sys.stderr, 'path '+output+' is not writable'
			sys.exit(1)

def science_run(time):
	if time >= 815155213 and time <= 875232014: return 's5'
	if time >= 931035296 and time <= 999999999: return 's6'
	print >>sys.stderr, "COULD NOT DETERMINE SCIENCE RUN from %d" % (int(time),)
	sys.exit(1)

def get_day_boundaries(time):

  # determine the start Time : 00:00:00 UTC from the day before
  # and the end time, 00:00:00 UTC the current day

  gps = LIGOTimeGPS(time)
  start_gps = int(date.utc_midnight(gps))
  end_gps = start_gps + 86400
  return str(start_gps),str(end_gps)


def figure_out_type(time, ifo=None, data_type='hoft'):
	"""
Run boundaries (from CBC analyses):
VSR1: 863557214 - 875232014
S5:   815155213 - 875232014
VSR2/S6: 931035296 - ...
Frame types for S5/VSR1:
() For RDS_R_L1 data set:
type    channel_name
RDS_R_L1     H1:LSC-DARM_ERR
RDS_R_L1     H2:LSC-DARM_ERR
RDS_R_L1     L1:LSC-DARM_ERR
() For hoft data:
type    channel_name
H1_RDS_C03_L2  H1:LSC-STRAIN
H2_RDS_C03_L2  H2:LSC-STRAIN
L1_RDS_C03_L2  L1:LSC-STRAIN
HrecV2_16384Hz      V1:h_16384Hz
Frame types for S6/VSR2:
() For RDS_R_L1 data set:
type    channel_name
H1_RDS_R_L1   H1:LSC-DARM_ERR
L1_RDS_R_L1   L1:LSC-DARM_ERR
() For hoft data:
H1_LDAS_C00_L2  H1:LDAS-STRAIN
L1_LDAS_C00_L2  L1:LDAS-STRAIN
HrecOnline      V1:h_16384Hz
	"""
	L1HoftTypes=(
		("L1_RDS_C03_L2","L1:LSC-STRAIN",815155213,875232014),
		("L1_LDAS_C00_L2","L1:LDAS-STRAIN",931035296,941997600),
		("L1_LDAS_C02_L2","L1:LDAS-STRAIN",941997600,999999999)
		)
	H1HoftTypes=(
		("H1_RDS_C03_L2","H1:LSC-STRAIN",815155213,875232014),
		("H1_LDAS_C00_L2","H1:LDAS-STRAIN",931035296,941997600),
		("H1_LDAS_C02_L2","H1:LDAS-STRAIN",941997600,999999999)
		)
	H2HoftTypes=(
		("H2_RDS_C03_L2","H2:LSC-STRAIN",815155213,875232014),
		("H1_LDAS_C00_L2","H1:LDAS-STRAIN",931035296,999999999)
		)
	V1HoftTypes=(
		("HrecV2_16384Hz","V1:h_16384Hz",863557214,875232014),
		("HrecOnline","V1:h_16384Hz",931035296,999999999)
		)
	L1RdsTypes=(
		("RDS_R_L1","L1:LSC-DARM_ERR",815155213,875232014),
		("L1_RDS_R_L1","L1:LSC-DARM_ERR",931035296,999999999)
		)
	H1RdsTypes=(
		("RDS_R_L1","H1:LSC-DARM_ERR",815155213,875232014),
		("H1_RDS_R_L1","H1:LSC-DARM_ERR",931035296,999999999)
		)
	H2RdsTypes=(
		("RDS_R_L1","H2:LSC-DARM_ERR",815155213,875232014),
		("H1_RDS_R_L1","H1:LSC-DARM_ERR",931035296,999999999)
		)
	V1RdsTypes=(
		("raw","V1:Pr_B1_ACp",863557214,875232014),
		("raw","V1:Pr_B1_ACp",931035296,999999999)
		)
	channelMap={
		"L1":{"hoft":L1HoftTypes,"rds":L1RdsTypes},
		"H1":{"hoft":H1HoftTypes,"rds":H1RdsTypes},
		"H2":{"hoft":H2HoftTypes,"rds":H2RdsTypes},
		"V1":{"hoft":V1HoftTypes,"rds":V1RdsTypes}
		}
	#Use the IFO type to select the channel type
	foundType=""
	foundChannel=""
	if ifo == None:
		print time," ifo argument to figure_out_type should not be null!"
		os.abort()
		
	for type,channel,start,stop in channelMap[ifo][data_type]:
		if ((start<=time) and (time<=stop)):
			foundType=type
			foundChannel=channel
			break
	if foundType == "":
		print time,ifo + " not found in method stfu_pipe.figure_out_type"
		os.abort()
	return str(foundType), str(foundChannel)

def figure_out_cache(time,ifo):

	cacheList=(
		(home_dirs()+"/romain/followupbackgrounds/omega/S5/background/background_815155213_875232014.cache",815155213,875232014,"H1H2L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/S6a/background/background_931035296_935798415.cache",931035296,935798415,"H1L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/S6b/background/background_937800015_944587815.cache",935798415,944587815,"H1L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/S6b/background/background_944587815_947260815.cache",944587815,947260815,"H1L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/S6c/background/background_948672015_961545615.cache",948672015,961545687,"H1L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/S6d/background/background_961545607_968803223.cache",961545687,999999999,"H1L1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/VSR2aRerun/background/background_931035296_935798415.cache",931035296,935798415,"V1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/VSR2bRerun/background/background_937800015_947260815.cache",935798415,947116815,"V1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/VSR3preCommissioning/background/background_966124815_968025615.cache",964310415,968284815,"V1"),
		(home_dirs()+"/romain/followupbackgrounds/omega/VSR3postCommissioning/background/background_968544015_971568015.cache",968284815,999999999,"V1")
		)

	foundCache = ""
	for cacheFile,start,stop,ifos in cacheList:
		if ((start<=time) and (time<stop) and ifo in ifos):
			foundCache = cacheFile
			break

	if 'phy.syr.edu' in get_hostname():
		foundCache = foundCache.replace("romain","rgouaty")

	if foundCache == "":
		print ifo, time, " not found in method stfu_pipe.figure_out_cache"	
	else:
		if not os.path.isfile(foundCache):
			print "file " + foundCache + " not found"
			foundCache = ""

	return foundCache

def home_dirs():
	return os.path.split(os.path.abspath(os.environ['HOME']))[0]

def get_hostname():
	host = socket.getfqdn()
	return host

###############################################################################
##### USEFULL FUNCTIONS CALLED BY PYTHON JOBS
###############################################################################

def getParamsFromCache(fileName,type,ifo=None,time=None):
	qscanList = []
	cacheList = lal.Cache.fromfile(open(fileName))
	if not cacheList:
		return qscanList
	cacheSelected = cacheList.sieve(description=type,ifos=ifo)
	if time:
		if math.floor(float(time)) != math.ceil(float(time)):
			cacheSelected = cacheSelected.sieve(segment=segments.segment(math.floor(float(time)), math.ceil(float(time))))
		else:
			cacheSelected = cacheSelected.sieve(segment=segments.segment(math.floor(float(time))-0.5, math.floor(float(time))+0.5))

	for cacheEntry in cacheSelected:
		path_output = cacheEntry.path
		time_output = str(cacheEntry.segment[0])
		type_output = cacheEntry.description
		ifo_output = cacheEntry.observatory
		qscanList.append([path_output,time_output,type_output,ifo_output])

	return qscanList

###############################################################################
##### CONDOR JOB CLASSES ######################################################
###############################################################################

# DO SOME STANDARD STUFF FOR FU JOBS
class FUJob(object):
	"""
	"""

	def __init__(self):
		pass

	def __conditionalLoadDefaults__(self,defaults,cp):
		if not(cp.has_section(defaults["section"])):
			cp.add_section(defaults["section"])
		for key, val in defaults["options"].iteritems():
			if not cp.has_option(defaults["section"], key):
				cp.set(defaults["section"], key, val)

	def setupJob(self, name='job', dir= '', tag_base=None, cp=None):
		# Give this job a name.  Then make directories for the log files and such
		# This name is important since these directories will be included in
		# the web tree.
                if dir and not os.path.isdir(dir):
		  os.mkdir(dir)
		self.name = name + "_" + tag_base + "_" + os.path.split(dir.rstrip('/'))[1]
                if dir:
		  self.relPath = dir + '/' + self.name
                else:
                  self.relPath = self.name
		self.outputPath = os.getcwd() + '/' + self.relPath + '/'
		self.tag_base = tag_base
		if not os.path.isdir(self.relPath):
			os.mkdir(self.relPath)
		if not os.path.isdir(self.relPath+'/logs'):
			os.mkdir(self.relPath+'/logs')
		if not os.path.isdir(self.relPath+'/Images'):
			os.mkdir(self.relPath+'/Images')
		if not os.path.isdir(self.relPath+'/DataProducts'):
			os.mkdir(self.relPath+'/DataProducts')
		# Set up the usual stuff and name the log files appropriately
		try: self.add_condor_cmd('environment',"KMP_LIBRARY=serial;MKL_SERIAL=yes")
		except: pass
		self.set_sub_file(self.name+'.sub')
		self.set_stdout_file(self.outputPath+'/logs/'+self.name+'-$(macroid).out')
		self.set_stderr_file(self.outputPath+'/logs/'+self.name+'-$(macroid).err')
		if cp:
			if cp.has_section("condor-memory-requirement") and cp.has_option("condor-memory-requirement",name):
				requirement = cp.getint("condor-memory-requirement",name)
				self.add_condor_cmd("Requirements", "(Memory > " + str(requirement) + ")")

# QSCAN JOB CLASS
class qscanJob(pipeline.CondorDAGJob, FUJob):
	"""
	A qscan job
	"""
	def __init__(self, opts, cp, dir='', tag_base=''):
		"""
		"""
		self.__executable = string.strip(cp.get('fu-condor','qscan'))
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)
		self.setup_checkForDir()
		self.setup_rm_lock()

	def is_dax(self):
		return False

	def setup_checkForDir(self):
		# create a shell script to check for the existence of the qscan output directory and rename it if needed
		checkdir_script = open('checkForDir.sh','w')
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
		os.chmod('checkForDir.sh',0755)

	def setup_rm_lock(self):
		rm_lock_script = open('rmLock.sh','w')
		rm_lock_script.write("#!/bin/bash\nif [ -e $1 ]\nthen\n\trm $1\nfi")
		rm_lock_script.close()
		os.chmod('rmLock.sh',0755)

# CLASS TO SETUP THE PROXY FILE NEEDED FOR REMOTE JOBS (CONDOR FLOCKING)

class setupProxyJob(pipeline.CondorDAGJob, FUJob):
	"""
	"""
	def __init__(self, opts, cp, dir='', tag_base=''):
		self.setup_proxy_script()
		self.__executable = "getProxy.sh"
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "local"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

	def setup_proxy_script(self):
		proxy_script = open('getProxy.sh','w')
		proxy_script.write("""#!/bin/bash
if [ ! -e \"proxy.pem\" ]
then
	file=`grid-proxy-info -path`
	cp ${file} proxy.pem
fi
		""")
		proxy_script.close()
		os.chmod('getProxy.sh',0755)

# DISTRIBUTE REMOTE QSCANS CLASS

class distribRemoteQscanJob(pipeline.CondorDAGJob, FUJob):
	"""
	This class sets up a script to be run as child of the remote scans in order to distribute its results to the appropriate paths. It takes the qscan tarball as input, uncompress it and copy the results to the path specified in cache file.
	Moreover this job also deletes the temporary remote datafind cache files in order to clean up the followup directory.
	"""
	def __init__(self, opts, cp, dir='', tag_base=''):
		"""
		"""
		self.setup_distrib_script(dir,tag_base)
		self.__executable = 'distribRemoteScan_'+dir+'_'+tag_base+'.sh'
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

	def setup_distrib_script(self,dir,tag_base):
		distrib_script = open('distribRemoteScan_'+dir+'_'+tag_base+'.sh','w')
		distrib_script.write("""#!/bin/bash
currentPath=`pwd` ;
mv $1 $2/. ;
cd $2 ;
tar -xzvf $1 ;
cd $currentPath ;
for figPath in `find $2/$3 -name "*.png" -print` ; do
	echo $figPath ;
	thumbPath=`echo $figPath | sed s/.png/.thumb.png/g` ;
	figDPI=120; 
	fullSize='600x';
	thumbSize='300x';
	convert -resize $thumbSize -strip -depth 8 -colors 256 $figPath $thumbPath  ;
done
rm $2/$1 ;
if [ -e $4 ]
then
	rm $4 ;
fi
		""")
		distrib_script.close()
		os.chmod('distribRemoteScan_'+dir+'_'+tag_base+'.sh',0755)

# REMOTE DATAFIND JOB

class remoteDatafindJob(pipeline.CondorDAGJob, FUJob):
	"""
	"""
	def __init__(self, opts, cp, tag_base='', dir=''):

		self.setup_df_submission_script(dir,tag_base)
		self.__executable = 'remoteDatafind_'+dir+'_'+tag_base+'.sh'
		self.dir = dir
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.add_condor_cmd('input','proxy.pem')
		self.add_condor_cmd('should_transfer_files','yes')
		self.add_condor_cmd('when_to_transfer_output','ON_EXIT_OR_EVICT')
		self.add_condor_cmd('+RunOnEGEEGrid','True')
		self.add_condor_cmd("Requirements","(Arch == \"INTEL\" || Arch == \"X86_64\" ) && ( Pilot_SiteName == \"Bologna\")")
		self.add_condor_cmd('transfer_output_files','$(macrooutputfile)')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

	def setup_df_submission_script(self,dir,tag_base):
		submit_script = open('remoteDatafind_'+dir+'_'+tag_base+'.sh','w')
		submit_script.write("""#!/bin/bash
. /opt/exp_software/virgo/etc/virgo-env.sh
. /opt/glite/etc/profile.d/grid-env.sh
export X509_USER_PROXY=`pwd`/proxy.pem
/opt/exp_software/virgo/lscsoft/etc/LSCdataFind --observatory $1 --gps-start-time $2 --gps-end-time $3 --url-type file --lal-cache --type $4 --output $5
outputCache=$5
outputQcache=${outputCache/.cache/.qcache}
/storage/gpfs_virgo3/virgo/omega/omega_r3270_glnxa64_binary/bin/convertlalcache $5 %s-%s-$outputQcache
		"""%(dir,tag_base))
		submit_script.close()
		os.chmod('remoteDatafind_'+dir+'_'+tag_base+'.sh',0755)

# REMOTE QSCAN CLASS

class remoteQscanJob(pipeline.CondorDAGJob, FUJob):
	"""
	A remote qscan job
	"""
	def __init__(self, opts, cp, dir='', tag_base=''):
		"""
		"""
		self.setup_submission_script(dir,tag_base)
		self.__executable = "remoteScan_"+dir+"_"+tag_base+".sh"
		self.dir = dir
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.add_condor_cmd('input','proxy.pem')
		self.add_condor_cmd('should_transfer_files','yes')
		self.add_condor_cmd('when_to_transfer_output','ON_EXIT_OR_EVICT')
		self.add_condor_cmd('+RunOnEGEEGrid','True')
		self.add_condor_cmd("Requirements","(Arch == \"INTEL\" || Arch == \"X86_64\" ) && ( Pilot_SiteName == \"Bologna\")")
		self.add_condor_cmd('transfer_output_files','$(macrooutputfile)')
		self.add_condor_cmd('transfer_input_files','$(macroinputfile)')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

	def setup_submission_script(self,dir,tag_base):
		submit_script = open('remoteScan_'+dir+'_'+tag_base+'.sh','w')
		submit_script.write("""#!/bin/bash
. /opt/exp_software/virgo/etc/virgo-env.sh
. /opt/glite/etc/profile.d/grid-env.sh
. /storage/gpfs_virgo3/virgo/omega/omega_env.sh
export X509_USER_PROXY=`pwd`/proxy.pem
/storage/gpfs_virgo3/virgo/omega/omega_r3270_glnxa64_binary/bin/wpipeline scan -r -c $1 -f $2 -o $3 $4

tar -czf %s-%s-$4.tgz $3
		"""%(dir,tag_base))
		submit_script.close()
		os.chmod('remoteScan_'+dir+'_'+tag_base+'.sh',0755)

# A CLASS TO ANALYSE QSCAN RESULTS
class analyseQscanJob(pipeline.CondorDAGJob,FUJob):

	def __init__(self,opts,cp,dir='',tag_base=''):

		self.__executable = string.strip(cp.get('fu-condor','analyseQscan'))
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.name_for_background = self.name + "_" + tag_base
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)


# A CLASS TO DO FOLLOWUP INSPIRAL JOBS 
class followUpInspJob(inspiral.InspiralJob,FUJob):

	def __init__(self,cp,dir='', tag_base=''):

		inspiral.InspiralJob.__init__(self,cp)

		if tag_base == 'head':
			self.set_executable(string.strip(cp.get('fu-condor','inspiral_head')))
		#self.name = 'followUpInspJob' + type

		self.name = os.path.split(self.get_executable().rstrip('/'))[1]
		self.setupJob(name=self.name, dir=dir, cp=cp, tag_base=tag_base)

# JOB CLASS FOR PRODUCING A SKYMAP
class skyMapJob(pipeline.CondorDAGJob,FUJob):
	"""
	Generates sky map data
	"""
	def __init__(self, options, cp, dir='', tag_base=''):
		"""
		"""
		#self.__prog__ = 'lalapps_skyMapJob'
		self.__executable = string.strip(cp.get('fu-condor','lalapps_skymap'))
		self.__universe = "standard"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')

		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,dir=dir, tag_base=tag_base)

		self.ra_res = cp.get("fu-skymap", 'ra-res')
		self.dec_res = cp.get("fu-skymap", 'dec-res')
		self.sample_rate = cp.get("fu-skymap", 'sample-rate')
		
# JOB CLASS FOR PRODUCING SKYMAP PLOT
class skyMapPlotJob(pipeline.CondorDAGJob,FUJob):
	"""
	Plots the sky map output of lalapps_skymap
	"""
	def __init__(self, options, cp, dir='',tag_base=''):
		"""
		"""
		#self.__prog__ = 'pylal_skyPlotJob'
		self.__executable = string.strip(cp.get('fu-condor','pylal_skyPlotJob'))
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,dir=dir,tag_base=tag_base)

		self.ra_res = cp.get("fu-skymap", 'ra-res')
		self.dec_res = cp.get("fu-skymap", 'dec-res')
		self.sample_rate = cp.get("fu-skymap", 'sample-rate')

# JOB CLASS FOR PLOTTING SNR AND CHISQ
class plotSNRChisqJob(pipeline.CondorDAGJob,FUJob):
	"""
	A followup plotting job for snr and chisq time series
	"""
	def __init__(self, options, cp, dir='', tag_base=''):
		"""
		"""
		#self.__prog__ = 'plotSNRCHISQJob'
		self.__executable = string.strip(cp.get('fu-condor','plotsnrchisq'))
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')

		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base,dir=dir)

class fuDataFindJob(pipeline.LSCDataFindJob,FUJob):
	def __init__(self, config_file, dir='', tag_base=''):
    
		#self.name = name
    
		# unfortunately the logs directory has to be created before we call LSCDataFindJob
		#try:
		#	os.mkdir(self.name)
		#	os.mkdir(self.name + '/logs')
		#except: pass

		self.__executable = string.strip(config_file.get('condor','datafind'))

		# MUST DO THIS FIRST SO THE LOG DIRECTORIES ARE MADE
		self.name = os.path.split(self.__executable.rstrip('/'))[1]

		self.setupJob(name=self.name, tag_base=tag_base, dir=dir)

		pipeline.LSCDataFindJob.__init__(self, self.relPath, self.relPath + '/logs', config_file)
		self.setup_cacheconv(config_file)
    
	def setup_cacheconv(self,cp):
		# create a shell script to call convertlalcache.pl if the value of $RETURN is 0
		convert_script = open('cacheconv.sh','w')
		#FIXME changed convert cache script to not fail on previous error?
		convert_script.write("""#!/bin/bash
%s ${1} ${2}
if [ ${3} = \'y\' ]; then
	cp ${2} .
fi
		""" % string.strip(cp.get('fu-condor','convertcache')))
		convert_script.close()
		os.chmod('cacheconv.sh',0755)

#This class is responsible for running the default job for making our
#wiki content
class makeCheckListWikiJob(pipeline.CondorDAGJob,FUJob):
	"""
	This actually launches a default wiki creation job
	"""
	def __init__(self,opts,cp,dir='',tag_base=''):
	    """
	    """
	    self.__executable = string.strip(cp.get("fu-condor",
						    "makeCheckListWiki").strip())
	    self.name = os.path.split(self.__executable.rstrip('/'))[1]
	    self.__universe = string.strip(cp.get("makeCheckListWiki",
						  "universe").strip())
	    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
	    self.add_condor_cmd('getenv','True')
	    self.setupJob(name=self.name,dir=dir,cp=cp,tag_base='_all')
#End makeCheckListWikiJob class

#This class is responsible for the followup page
class makeFollowupPageJob(pipeline.CondorDAGJob,FUJob):
	"""
	This actually launches a followup page job
	"""
	def __init__(self,opts,cp,dir='',tag_base=''):
	    """
	    """
	    self.__executable = string.strip(cp.get("fu-condor",
						"lalapps_followup_page").strip())
	    self.name = os.path.split(self.__executable.rstrip('/'))[1]
	    self.__universe = "vanilla" 
	    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
	    self.add_condor_cmd('getenv','True')
	    self.setupJob(name=self.name,dir=dir,cp=cp,tag_base='_all')
#End makeFollowupPageJob class

	    
#The class responsible for running the data quality flag finding job
class findFlagsJob(pipeline.CondorDAGJob, FUJob):
	"""
	A job which queries the ldbd database holding segment
	information about the DQ flags.
	"""
	defaults={"section":"fu-condor",
		  "options":{"universe":"local",
			     "dqflags":"followupQueryDQ.py"}
		  }

	def __init__(self, opts, cp, dir='', tag_base=""):
		"""
		"""
		self.__conditionalLoadDefaults__(findFlagsJob.defaults,cp)
		#self.__prog__ = 'findFlagsJob'
		self.__executable = string.strip(cp.get('fu-condor','dqflags'))
		self.__universe = string.strip(cp.get('fu-condor','universe'))
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)

#The class responsible for checking for know veto flags
class findVetosJob(pipeline.CondorDAGJob,FUJob):
	"""
	A job instance that queries the segment database for known
	types of active veto intervals.
	"""
	defaults={"section":"fu-condor",
		  "options":{"universe":"local",
			     "vetoflags":"followupQueryVeto.py"}
		  }
	def __init__(self, opts,cp, dir='', tag_base=""):
		"""
		"""
		self.__conditionalLoadDefaults__(findVetosJob.defaults,cp)
		#self.__prog__ = 'findVetosJob'
		self.__executable = string.strip(cp.get('fu-condor','vetoflags'))
		self.__universe = string.strip(cp.get('fu-condor','universe'))
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)

#The class responsible for Job Object running the customFOM builder python
#script!
class customFOMPlotJob(pipeline.CondorDAGJob,FUJob):
	"""
	This is a job class which allows us to wrap up the script for
	creating customized figure of merit(FOM) plots.  The script,
	followupCustomFOM.py, acutally contains a call to
	ligo_data_find, via a subprocess.  This removes our need
	to have a datafind parent job.
	"""
	defaults={"section":"fu-condor",
		  "options":{"universe":"vanilla",
			     "customfom":"followupCustomFOM.py"}
		  }
	def __init__(self, opts, cp, dir='', tag_base=""):
		"""
		"""
		self.__conditionalLoadDefaults__(customFOMPlotJob.defaults,cp)
		#self.__prog__ = 'customFOMPlotJob'
		self.__executable = string.strip(cp.get('fu-condor','customfom'))
		self.__universe = string.strip(cp.get('fu-condor','universe'))
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)


#The class responsible for running one type of parameter consistency check
class effDRatioJob(pipeline.CondorDAGJob,FUJob):
	"""
	A job that performs parameter consitency check for a trigger
	being followed up.
	"""
	defaults={"section":"fu-condor",
		  "options":{"universe":"local",
			     "effDRatio":"followupRatioTest.py"}
		  }
	def __init__(self, opts, cp, dir='', tag_base=""):
		"""
		"""
		self.__conditionalLoadDefaults__(effDRatioJob.defaults,cp)
		#self.__prog__ = 'effDRatioTest'
		self.__executable = string.strip(cp.get('fu-condor','effDRatio'))
		self.__universe = string.strip(cp.get('fu-condor','universe'))
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)

# Follow up chia job
class followUpChiaJob(inspiral.ChiaJob,FUJob):
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

	def __init__(self, options, cp, dir='', tag_base=''):
		"""
		"""
		self.__conditionalLoadDefaults__(followUpChiaJob.defaults,cp)
		#self.__prog__ = 'followUpChiaJob'
		self.__executable = string.strip(cp.get('condor','chia'))
		self.__universe = "standard"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')

		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)


##############################################################################
# jobs class for plotting coherent inspiral search and null stat timeseries

class plotChiaJob(pipeline.CondorDAGJob,FUJob):
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

	def __init__(self, options, cp, dir, tag_base=''):
		"""
		"""
		#if not(verifyCP(cp,self.defaults)): modifyCP(cp,self.defaults)
		self.__executable = string.strip(cp.get('fu-condor','plotchiatimeseries'))
		self.__universe = "vanilla"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,tag_base=tag_base, dir=dir)

##############################################################################
# jobs class for setting a mcmc run

class mcmcJob(pipeline.CondorDAGJob, FUJob):
	"""
	A job to set up a mcmc run
	"""
	def __init__(self,opts,cp,dir='',tag_base=''):
		"""
		"""
		self.__executable = string.strip(cp.get('fu-condor','mcmc'))
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "standard"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

##############################################################################
# jobs class for setting a spinmcmc run

class spinmcmcJob(pipeline.CondorDAGJob, FUJob):
	"""
	A job to set up a spinmcmc run
	"""
	def __init__(self,opts,cp,dir='',tag_base=''):
		self.__executable = string.strip(cp.get('fu-condor','spinmcmc'))
		self.name = os.path.split(self.__executable.rstrip('/'))[1]
		self.__universe = "standard"
		pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)


##############################################################################
# jobs class for setting a the plotting of mcmc results

class plotmcmcJob(pipeline.CondorDAGJob, FUJob):
        """
        A job to set up a plotmcmc run
        """
        def __init__(self,opts,cp,dir='',tag_base=''):
                """
                """
                self.__executable = string.strip(cp.get('fu-condor','plotmcmc'))
                self.name = os.path.split(self.__executable.rstrip('/'))[1]
                self.__universe = "vanilla"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
		self.add_condor_cmd('getenv','True')
		self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)

##############################################################################
# jobs class for setting a the plotting of spinmcmc results

class plotspinmcmcJob(pipeline.CondorDAGJob, FUJob):
        """
        A job to set up a plotspinmcmc run
        """
        def __init__(self,opts,cp,dir='',tag_base=''):
                """
                """
                self.__executable = string.strip(cp.get('fu-condor','plotspinmcmc'))
                self.name = os.path.split(self.__executable.rstrip('/'))[1]
                self.__universe = "vanilla"
                pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
                self.add_condor_cmd('getenv','True')
                self.setupJob(name=self.name,dir=dir,cp=cp,tag_base=tag_base)


#############################################################################
###### CONDOR NODE CLASSES ##################################################
#############################################################################

class FUNode:
	"""
	"""

	def __init__(self):
		self.invalidate()

	def validate(self):
		self.validNode = True

	def invalidate(self):
		self.validNode = False

	def __conditionalLoadDefaults__(self,defaults,cp):
		if not(cp.has_section(defaults["section"])):
			cp.add_section(defaults["section"])
		for key, val in defaults["options"].iteritems():
			if not cp.has_option(defaults["section"], key):
				cp.set(defaults["section"], key, val)

	def setupPlotNode(self, job):
		self.add_var_opt("output-path",job.outputPath)
		self.add_var_opt("enable-output","")

# QSCAN NODE
class fuQscanNode(pipeline.CondorDAGNode,FUNode):
	"""
QScan node.  This node writes its output to the web directory specified in
the inifile + the ifo and gps time.  For example:
 
	/archive/home/channa/public_html/followup/htQscan/H1/999999999.999

The omega scan command line is 

	wpipeline scan -r -c H1_hoft.txt -f H-H1_RDS_C03_L2-870946612-870946742.qcache -o QSCAN/foreground-hoft-qscan/H1/870946677.52929688 870946677.52929688

	"""
	def __init__(self, dag, job, cp, opts, time, ifo, frame_cache, p_nodes=[], type="ht", variety="fg"):
		"""
		"""
		pipeline.CondorDAGNode.__init__(self,job)

		self.scan_type = variety.upper() + "_" + type.replace("seismic","seis").upper()
		self.scan_ifo = ifo
		preString="omega/"+ifo+"/%s/"+science_run(time).upper()+""
		if variety == "bg":
			self.add_var_arg('scan')
			preString = preString%("background")
			oldPreString="omega/" + science_run(time).upper() + "/background"
		else:
			self.add_var_arg('scan -r')
			preString = preString%("foreground")
			oldPreString="omega/" + science_run(time).upper() + "/foreground"
		config = self.fix_config_for_science_run( cp.get('fu-'+variety+'-'+type+'-qscan', ifo+'config').strip(), time )
		if cp.get('fu-'+variety+'-'+type+'-qscan', ifo+'config').strip() != config:
			cp.set('fu-'+variety+'-'+type+'-qscan',ifo+'config',config)
		self.add_var_arg("-c " + config )

		if type == "ht":
			dataString = figure_out_type(time, ifo, 'hoft')[0]
		else:
			dataString = figure_out_type(time, ifo, 'rds')[0]
		if type == "seismic":
			dataString = dataString + "_SEIS"
		if dataString[:2]!=ifo:
			dataString = ifo + "_" + dataString
		timeString = "-".join(get_day_boundaries(int(time)))
		if cp.has_option('fu-output','output-dir') and cp.get('fu-output','output-dir'):
		  output = cp.get('fu-output','output-dir') + '/' + preString + '/' + dataString + '/' + timeString
		else:
		  output = os.getcwd() + '/' + preString + '/' + dataString + '/' + timeString

		# CREATE AND MAKE SURE WE CAN WRITE TO THE OUTPUT DIRECTORY
		mkdir(output)

		output_path = output+"/"+str(time)
		self.add_var_arg("-o " + output_path)
		
		self.output_cache = lal.CacheEntry(ifo, job.name.upper(), segments.segment(float(time), float(time)), "file://localhost/"+output_path)
		# ADD FRAME CACHE FILE
		self.add_var_arg("-f "+frame_cache)
		
		self.add_var_arg(repr(time))

		self.set_pre_script( "checkForDir.sh %s %s" %(output, str(time)) )
		#FIXME is deleting the lock file the right thing to do?
		self.set_post_script( "rmLock.sh %s/%s/lock.txt" %(output, str(time)) )

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if not(cp.has_option('fu-remote-jobs','remote-jobs') and job.name in cp.get('fu-remote-jobs','remote-jobs') and cp.has_option('fu-remote-jobs','remote-ifos') and ifo in cp.get('fu-remote-jobs','remote-ifos')):
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			if not (type=="ht" and opts.no_ht_qscan) and not (type=="rds" and opts.no_rds_qscan) and not (type=="seismic" and opts.no_seismic_qscan):
				dag.add_node(self)
				self.validate()
			else:
				self.invalidate()
		else:
			self.invalidate()

	def fix_config_for_science_run(self, config, time):
		run = science_run(time)
		config_path = os.path.split(config)
		out = "/".join([config_path[0], config_path[1].replace('s5',run).replace('s6',run)])
		return out

# SETUP PROXY FOR REMOTE SCANS
class setupProxyNode(pipeline.CondorDAGNode,FUNode):

	def __init__(self, dag, job, cp, opts):

		pipeline.CondorDAGNode.__init__(self,job)
		dag.add_node(self)
		self.validate()

# DISTRIBUTE REMOTE QSCAN RESULTS
class distribRemoteQscanNode(pipeline.CondorDAGNode,FUNode):

	def __init__(self, dag, job, cp, opts, ifo, time, p_nodes=[], type=""):

		pipeline.CondorDAGNode.__init__(self,job)
		self.scan_type = type.replace("seismic","seis").upper()
		self.scan_ifo = ifo
		# WARNING: First element in p_nodes list is assumed to be the omega scan node
		self.add_var_arg(p_nodes[0].name_output_file)
		self.add_var_arg(p_nodes[0].output_path)
		self.add_var_arg(str(time))
		# WARNING: Second element in p_nodes list is assumed to be the datafind node
		self.add_var_arg(p_nodes[1].localFileName)

		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		dag.add_node(self)
		self.validate()

# REMOTE DATAFIND NODE
class remoteDatafindNode(pipeline.CondorDAGNode,FUNode):

	def __init__(self, dag, job, cp, opts, ifo, time, data_type="rds", p_nodes=[]):
		pipeline.CondorDAGNode.__init__(self,job)

		type, channel = figure_out_type(time,ifo,data_type)
		q_time = float(cp.get("fu-q-"+data_type+"-datafind",ifo+"-search-time-range"))/2.
		start_time = int( time - q_time - 1.)
		end_time = int( time + q_time + 1.)
		outputFileName = ifo[0]+'-'+type+'-'+str(start_time)+'-'+str(end_time)+'.qcache'

		# THIS IS THE DIRECTORY WHERE THE DATA WILL ULTIMATELY BE COPIED ONCE THE DATAPRODUCT ARE SENT BACK LOCALLY
		#self.output_cache = lal.CacheEntry(ifo, job.name.upper(), segments.segment(start_time, end_time), "file://localhost/"+job.outputPath+'/'+os.path.abspath(outputFileName))
	
		self.add_var_arg(ifo[0])
		self.add_var_arg(str(start_time))
		self.add_var_arg(str(end_time))
		self.add_var_arg(type)
		self.add_var_arg(outputFileName.replace('qcache','cache'))

		self.name_output_file = job.dir + "-" + job.tag_base + "-" + outputFileName
		self.add_macro("macrooutputfile", self.name_output_file)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())
		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		dag.add_node(self)
		self.validate()


# REMOTE QSCAN NODE
class fuRemoteQscanNode(pipeline.CondorDAGNode,FUNode):
	"""
	"""
	def __init__(self, dag, job, cp, opts, time, ifo, p_nodes=[], type="ht", variety="fg"):

		pipeline.CondorDAGNode.__init__(self,job)

		self.scan_type = variety.upper() + "_" + type.replace("seismic","seis").upper()
		self.scan_ifo = ifo
		preString="omega/"+ifo+"/%s/"+science_run(time).upper()+""
		if variety == "bg":
			preString = preString%("background")
		else:
			preString = preString%("foreground")
		config = cp.get('fu-'+variety+'-'+type+'-qscan', ifo+'config').strip()
		self.add_var_arg( config )

		if type == "ht":
			dataString = figure_out_type(time, ifo, 'hoft')[0]
		else:
			dataString = figure_out_type(time, ifo, 'rds')[0]
		if type == "seismic":
			dataString = dataString + "_SEIS"
		if dataString[:2]!=ifo:
			dataString = ifo + "_" + dataString
		timeString = "-".join(get_day_boundaries(int(time)))
		if cp.has_option('fu-output','output-dir') and cp.get('fu-output','output-dir'):
			output = cp.get('fu-output','output-dir') + '/' + preString + '/' + dataString + '/' + timeString
		else:
			output = os.getcwd() + '/' + preString + '/' + dataString + '/' + timeString

		# CREATE AND MAKE SURE WE CAN WRITE TO THE OUTPUT DIRECTORY
		mkdir(output)

		# THIS IS THE DIRECTORY WHERE THE DATA WILL ULTIMATELY BE COPIED ONCE THE DATAPRODUCT ARE SENT BACK LOCALLY
		self.output_path = output

		self.output_cache = lal.CacheEntry(ifo, job.name.replace("remoteScan_"+job.dir+"_"+job.tag_base+".sh","wpipeline").upper(), segments.segment(float(time), float(time)), "file://localhost/"+self.output_path+"/"+str(time))

		# ADD FRAME CACHE FILE
		#self.add_var_arg("/storage/gpfs_virgo3/virgo/omega/cbc/S6/foreground/RAW/V-raw-930000000-947260815.qcache")

		# The first parent node must be the cache file!
		input_cache_file = p_nodes[0].localFileName
		self.add_var_arg(input_cache_file)
		self.add_macro("macroinputfile", input_cache_file)

		# NOW WE NEED TO SET UP THE REMOTE OUTPUTPATH

		# String used in the naming of the omega scan directory
		self.add_var_arg(str(time))

		# Time at which the omega scan is performed
		self.add_var_arg(repr(time))

		self.name_output_file = job.dir + "-" + job.tag_base + "-" + repr(time) + ".tgz"
		self.add_macro("macrooutputfile", self.name_output_file)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		dag.add_node(self)
		self.validate()


# ANALYSEQSCAN NODE
class analyseQscanNode(pipeline.CondorDAGNode,FUNode):

	def __init__(self, dag, job, cp, opts, time, ifo):

		pipeline.CondorDAGNode.__init__(self,job)

		name = job.name

		if "SEIS" in name:
			data_type = "rds"
			shortName = "seis_rds"
		elif "HT" in name:
			data_type = "hoft"
			shortName = "ht"
		else:
			data_type = "rds"
			shortName = "rds"

		refChannel = figure_out_type(time, ifo, data_type)[1].split(":")[-1]
		self.add_var_opt('ref-channel',refChannel)
		if cp.has_option('fu-analyse-qscan','generate-qscan-xml'):
			self.add_var_opt('generate-qscan-xml','')
		self.add_var_opt('z-threshold',cp.getfloat('fu-analyse-qscan','z-threshold'))
		if cp.has_option('fu-analyse-qscan','plot-z-distribution'):
			self.add_var_opt('plot-z-distribution','')
			self.add_var_opt('z-min',cp.getfloat('fu-analyse-qscan','z-min'))
			self.add_var_opt('z-max',cp.getfloat('fu-analyse-qscan','z-max'))
			self.add_var_opt('z-bins',cp.getfloat('fu-analyse-qscan','z-bins'))
		if cp.has_option('fu-analyse-qscan','plot-dt-distribution'):
			self.add_var_opt('plot-dt-distribution','')
			self.add_var_opt('dt-min',cp.getfloat('fu-analyse-qscan',shortName.replace('_','-') + '-dt-min'))
			self.add_var_opt('dt-max',cp.getfloat('fu-analyse-qscan',shortName.replace('_','-') + '-dt-max'))
			self.add_var_opt('dt-bins',cp.getfloat('fu-analyse-qscan','dt-bins'))
		if cp.has_option('fu-analyse-qscan','plot-z-scattered'):
			self.add_var_opt('plot-z-scattered','')
		if cp.has_option('fu-analyse-qscan','plot-z-scattered') or cp.has_option('fu-analyse-qscan','plot-dt-distribution'):
			self.add_var_opt('ref-channel',refChannel)
		self.add_var_opt('ifo-times',ifo)
		self.add_var_opt('type',name.upper().replace("ANALYSEQSCAN.PY","WPIPELINE"))
		self.add_var_opt('short-type',job.name_for_background.upper().replace("ANALYSEQSCAN.PY","WPIPELINE")+'_')
		self.add_var_opt('gps-string',str(time))
		self.add_var_opt('ifo-tag',ifo)
		self.add_var_opt('user-tag',str(time).replace('.','_') + "_" + shortName)

		self.add_var_opt('qscan-cache-foreground',dag.basename+'.cache')
		
		if cp.has_option('fu-analyse-qscan',ifo+'-background-cache'):
			backgroundCache = cp.get('fu-analyse-qscan',ifo+'-background-cache').strip()
		else:
			backgroundCache = figure_out_cache(time,ifo)
			cp.set('fu-analyse-qscan',ifo+'-background-cache',backgroundCache)
		self.add_var_opt('qscan-cache-background',backgroundCache)

		self.output_file_name = "%s-analyseQscan_%s_%s-unspecified-gpstime.cache" % ( ifo, ifo, str(time).replace('.','_') + "_" + shortName)
		self.output_cache = lal.CacheEntry(ifo,job.name.upper(),segments.segment(float(time),float(time)),"file://localhost/"+job.outputPath+'/'+self.output_file_name)

		self.setupPlotNode(job)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		# add the parents to this node
		for node in dag.get_nodes():
			# if node distributeQscanNode is valid and if remote
			# ifo is analysed, add distributeQscanNode as parent
			#if isinstance(node,distributeQscanNode):
			#	if node.validNode:
			#		self.add_parent(node)
			# add all qscan nodes of the same type as parents
			if isinstance(node,fuQscanNode) or isinstance(node,distribRemoteQscanNode):
				if node.validNode:
					if (node.scan_type in name and node.scan_ifo == ifo):
						self.add_parent(node)
		if not (shortName=="ht" and opts.no_ht_analyseQscan) and not (shortName == "rds" and opts.no_rds_analyseQscan) and not (shortName == "seis_rds" and opts.no_seismic_analyseQscan) and backgroundCache:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()


# DATAFIND NODE
class fuDataFindNode(pipeline.LSCDataFindNode,FUNode):
    
	def __init__(self, dag, job, cp, opts, ifo, sngl=None, qscan=False, trigger_time=None, data_type="hoft", p_nodes=[]):

		self.outputFileName = ""
		pipeline.LSCDataFindNode.__init__(self,job)
		if qscan:
			if sngl: time = sngl.time
			else: time = trigger_time
			self.outputFileName = self.setup_qscan(job, cp, time, ifo, data_type)
		else:
			if not sngl:
				print >> sys.stderr, "argument \"sngl\" should be provided to class fuDataFindNode"
				sys.exit(1)
			self.outputFileName = self.setup_inspiral(job, cp, sngl, ifo)

		self.output_cache = lal.CacheEntry(ifo, job.name.upper(), segments.segment(self.get_start(), self.get_end()), "file://localhost/"+os.path.abspath(self.outputFileName))

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if not(cp.has_option('fu-remote-jobs','remote-jobs') and job.name in cp.get('fu-remote-jobs','remote-jobs') and cp.has_option('fu-remote-jobs','remote-ifos') and ifo in cp.get('fu-remote-jobs','remote-ifos')) or opts.do_remoteScans:
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			if not (data_type=="hoft" and not qscan and opts.no_insp_datafind) and not (data_type=="hoft" and qscan and opts.no_htQscan_datafind) and not (data_type=="rds" and opts.no_rdsQscan_datafind):
				dag.add_node(self)
				self.validate()
			else:
				self.invalidate()
		else:
			self.invalidate

	def setup_qscan(self, job, cp, time, ifo, data_type):
		# 1s is substracted to the expected startTime to make sure the window
		# will be large enough. This is to be sure to handle the rouding to the
		# next sample done by qscan.
		type, channel = figure_out_type(time,ifo,data_type)
		self.set_type(type)
		self.q_time = float(cp.get("fu-q-"+data_type+"-datafind",ifo+"-search-time-range"))/2.
		self.set_observatory(ifo[0])
		self.set_start(int( time - self.q_time - 1.))
		self.set_end(int( time + self.q_time + 1.))
		lalCache = self.get_output()
		qCache = lalCache.rstrip("lcf") + "qcache"

		if cp.has_option('fu-remote-jobs','remote-jobs') and job.name in cp.get('fu-remote-jobs','remote-jobs') and cp.has_option('fu-remote-jobs','remote-ifos') and ifo in cp.get('fu-remote-jobs','remote-ifos'):
			self.add_var_arg('--server ldr-bologna.phys.uwm.edu')
			postScriptTest = "y"
			self.localFileName = os.path.basename(qCache)
		else:
			self.add_var_arg('')
			postScriptTest = "n"

		self.set_post_script(os.getcwd()+"/cacheconv.sh %s %s %s" %(lalCache,qCache,postScriptTest) )
		return(qCache)

	def setup_inspiral(self, job, cp, sngl, ifo):
		# 1s is substracted to the expected startTime to make sure the window
		# will be large enough. This is to be sure to handle the rouding to the
		# next sample done by qscan.
		type, channel = figure_out_type(sngl.get_gps_start_time(),ifo)
		self.set_type(type)
		self.set_observatory(ifo[0])
		#FIXME use proper pad, not hardcode to 64
		self.set_start(sngl.get_gps_start_time()-64)
		self.set_end(sngl.get_gps_end_time()+64)
		self.add_var_arg('')
		lalCache = self.get_output()
		return(lalCache)


# INSPIRAL NODE
class followUpInspNode(inspiral.InspiralNode,FUNode):

  #def __init__(self, inspJob, procParams, ifo, trig, cp,opts,dag, datafindCache, d_node, datafindCommand, type='plot', sngl_table = None):
	def __init__(self, dag, job, cp, opts, sngl, frame_cache, chia, tag, p_nodes=[]):

		tlen = 1.0
		self.output_file_name = ""
		pipeline.CondorDAGNode.__init__(self,job)
		pipeline.AnalysisNode.__init__(self)

		#FIXME HANDLE INJECTION FILES AND datafind cache
		# injFile = self.checkInjections(cp)
		# self.set_injections( injFile )

		self.set_trig_start( int(sngl.time - tlen + 0.5) )
		self.set_trig_end( int(sngl.time + tlen + 0.5) )
                if not chia:
		  self.add_var_opt("write-snrsq","")
		  self.add_var_opt("write-chisq","")
		  self.add_var_opt("write-spectrum","")
		  self.add_var_opt("write-template","")
		self.add_var_opt("write-cdata","")

		skipParams = ['minimal-match', 'bank-file', 'user-tag', 'injection-file', 'trig-start-time', 'trig-end-time']

		extension = ".xml"
		for row in sngl.process_params:
			param = row.param.strip("-")
			value = row.value
			# override some options
			if param == 'frame-cache': value = frame_cache
			if param == 'snr-threshold': value = "0.1"
			if param == 'do-rsq-veto': continue
			if param == 'enable-rsq-veto': continue
			if param == 'chisq-threshold': value = "1.0e+06"
			if param == 'cluster-method': value = 'window'
			if param == 'cluster-window': continue
			if param == 'userTag': continue
			if param == 'user-tag': continue
			if param in skipParams: continue
			if param == 'channel-name':
				self.inputIfo = value[0:2]
				#HACK FOR GRB FOLLOWUPS: Channel names defined
				#in old GRB runs are obsolete. It is better to
				#figure out the channel name from the GPS time.
				if opts.do_grb:
					type,channel = figure_out_type(sngl.time,self.inputIfo)
					value = channel
			if param == 'injection-file': value = sngl.inj_file_name
			if param == 'gps-end-time':
				self.set_end(int(value))
				continue
			if param == 'gps-start-time':
				self.set_start(int(value))
				continue
			if param == 'ifo-tag':
				self.set_ifo_tag(value)
				continue
			self.add_var_opt(param,value)
			if param == 'pad-data':
				self.set_pad_data(int(value))
			if param == 'write-compress':
				extension = '.xml.gz'

		self.add_var_opt('cluster-window',str( tlen / 2.))
		self.add_var_opt('disable-rsq-veto',' ')
		bankFile = self.write_trig_bank(sngl, 'trig_bank/' + sngl.ifo + '-TRIGBANK_FOLLOWUP_' + repr(sngl.time) + '.xml.gz')
		self.set_bank(bankFile)

                if chia:
		  self.set_user_tag( tag.upper() + "_CHIA_FOLLOWUP_" + repr(sngl.time) )
                else:     
                  self.set_user_tag( tag.upper() + "_FOLLOWUP_" + repr(sngl.time) )

                self.output_file_name = job.outputPath + sngl.ifo + "-INSPIRAL_" + self.get_ifo_tag() + "_" + self.get_user_tag() + "-" + str(self.get_start()) + "-" + str(int(self.get_end())-int(self.get_start())) + extension
		self.outputCache = sngl.ifo + ' ' + 'INSPIRAL' + ' ' + str(self.get_start()) + ' ' + str(int(self.get_end())-int(self.get_start())) + ' ' + self.output_file_name  + '\n' + sngl.ifo + ' ' + 'INSPIRAL-FRAME' + ' ' + str(self.get_start()) + ' ' + str(int(self.get_end())-int(self.get_start())) + ' ' + self.output_file_name.replace(extension,".gwf") + '\n'

		self.add_var_opt("output-path",job.outputPath)
		self.output_cache = []
		self.output_cache.append(lal.CacheEntry(sngl.ifo, job.name.upper(), segments.segment(float(self.get_start()), float(self.get_end())), "file://localhost/"+self.output_file_name))
		self.output_cache.append(lal.CacheEntry(sngl.ifo, job.name.upper(), segments.segment(float(self.get_start()), float(self.get_end())), "file://localhost/"+self.output_file_name.replace(extension,'.gwf')))
		
		self.output_frame_file = self.output_file_name.replace(extension,'.gwf')

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

# 		# Wed-Aug-25-2010:201008251418 Added Pre & Post
# 		# scripts depends on value of output-path
# 		patchScript=create_default_config().which("followup_InspiralDataMover.sh")
# 		self.set_pre_script("%s %s"%(patchScript,job.outputPath))
# 		self.set_post_script("%s %s"%(patchScript,job.outputPath))
# 		# End temporary additions Wed-Aug-25-2010:201008251421 
# 		#add parents and put node in dag
		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_inspiral:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

	def write_trig_bank(self,sngl, name):
		try:
			os.mkdir('trig_bank')
		except: pass
		xmldoc = ligolw.Document()
		xmldoc.appendChild(ligolw.LIGO_LW())

		process_params_table = lsctables.New(lsctables.ProcessParamsTable)
		xmldoc.childNodes[-1].appendChild(process_params_table)

		sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
		xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
		sngl_inspiral_table.append(sngl.row)

		utils.write_filename(xmldoc, name, verbose=False, gz = True)
		return name

# Create checklist wiki files etc node
class makeCheckListWikiNode(pipeline.CondorDAGNode,FUNode):
	"""
	This class is responsible for running a final job which will
	create the default top 10 triggers for each trigger type.
	This will place these files into the publication directory so
	user can push wiki content onto the CBC wiki
	"""
	def __init__(self,dag,job,cp,opts):
		pipeline.CondorDAGNode.__init__(self,job)
		#Specify pipe location
		self.add_var_opt('followup-directory',cp.get("makeCheckListWiki",
							     "location").strip())
		#Specify pipe ini file
		self.add_var_opt('ini-file',cp.get("makeCheckListWiki",
						   "ini-file").strip())
		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())
		#Add this as child of all known jobs
		for parentNode in dag.get_nodes():
			if not isinstance(parentNode,makeFollowupPageNode):
				self.add_parent(parentNode)
		if not opts.no_makeCheckList:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

# Create followup page node
class makeFollowupPageNode(pipeline.CondorDAGNode,FUNode):
	"""
	This runs the followup page
	"""
	def __init__(self,dag,job,cp,opts):
		pipeline.CondorDAGNode.__init__(self,job)
		#FIXME Specify cache location (sortof hard coded)
		self.add_var_arg('followup_pipe.cache')
		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())
		#Add this as child of all known jobs
		for parentNode in dag.get_nodes():
			self.add_parent(parentNode)
		dag.add_node(self)
		self.validate()

# FIND FLAGS NODE 
class findFlagsNode(pipeline.CondorDAGNode,FUNode):
	"""
	This class is resposible for setting up a node to perform a
	query for the DQ flag for the trigger which under review.
	EXAMPLE
	followupQueryDQ.py --window=60,15 --trigger-time=929052945 --output-format=moinmoin --segment-url="https://segdb.ligo.caltech.edu:30015" --output-file=dqResults.wiki
	"""
	defaults={"section":"findFlags",
		  "options":{"window":"60,15",
			     "segment-url":"https://segdb.ligo.caltech.edu",
			     "output-format":"moinmoin",
			     "output-file":"dqResults.wiki",
			     "estimate-background":"",
			     "background-location":"automatic"}
		  }
	def __init__(self, dag, job, cp, opts, coincEvent=None):
		"""
		"""
		self.__conditionalLoadDefaults__(findFlagsNode.defaults,cp)
		pipeline.CondorDAGNode.__init__(self,job)
		self.add_var_opt("trigger-time",coincEvent.time)
		#Output filename
		oFilename="%s-findFlags_%s_%s.wiki"%(coincEvent.instruments,
						     coincEvent.ifos,
						     coincEvent.time)		
		self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
		self.add_var_opt("segment-url",cp.get('findFlags','segment-url'))
		self.add_var_opt("output-format",cp.get('findFlags','output-format'))
		self.add_var_opt("window",cp.get('findFlags','window'))
		if cp.has_option('findFlags','estimate-background'):
			self.add_var_opt("estimate-background",cp.get('findFlags','estimate-background'))
		if cp.has_option('findFlags','background-location'):
			self.add_var_opt("background-location",cp.get('findFlags','background-location'))
		if cp.has_option('findFlags','blind'):
			self.add_var_opt("blind",cp.get('findFlags','blind'))
		self.output_cache = lal.CacheEntry(coincEvent.ifos, job.name.upper(), segments.segment(float(coincEvent.time), float(coincEvent.time)), "file://localhost/"+job.outputPath+'/DataProducts/'+oFilename)

		#IFO arg string
		myArgString=""
		if hasattr(coincEvent, "sngl_inspiral"):
			for sngl in coincEvent.sngl_inspiral.itervalues():
				myArgString=myArgString+"%s,"%sngl.ifo
		elif hasattr(coincEvent, "ifos_list"):
			for ifo in coincEvent.ifos_list:
				myArgString=myArgString+"%s,"%ifo
		myArgString=myArgString.rstrip(",")
		self.add_var_opt("ifo-list",myArgString)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if not opts.no_findFlags:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

# FIND VETOS NODE 
class findVetosNode(pipeline.CondorDAGNode,FUNode):
	"""
	This class is responsible for creating a node in the dag which
	queries the segment database for veto segments active around
	the trigger time of the candidate.
	Command line example:
	followupQueryVeto.py --window=60,15 --trigger-time=929052945 --output-format=moinmoin --segment-url="https://segdb.ligo.caltech.edu:30015" --output-file=vetoResults.wiki
	"""
	defaults={"section":"findVetoes",
		  "options":{"window":"60,15",
			     "segment-url":"https://segdb.ligo.caltech.edu",
			     "output-format":"moinmoin",
			     "output-file":"vetoResults.wiki",
			     "estimate-background":"",
			     "background-location":"automatic"}
		  }
	def __init__(self, dag, job, cp, opts, coincEvent=None):
		"""
		"""
		self.__conditionalLoadDefaults__(findVetosNode.defaults,cp)
		pipeline.CondorDAGNode.__init__(self,job)
		self.add_var_opt("trigger-time",coincEvent.time)
		#Output filename
		oFilename="%s-findVetos_%s_%s.wiki"%(coincEvent.instruments,
						     coincEvent.ifos,
						     coincEvent.time)
		self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
		self.add_var_opt("segment-url",cp.get('findVetoes','segment-url'))
		self.add_var_opt("output-format",cp.get('findVetoes','output-format'))
		self.add_var_opt("window",cp.get('findVetoes','window'))
		if cp.has_option('findVetoes','estimate-background'):
			self.add_var_opt("estimate-background",cp.get('findVetoes','estimate-background'))
		if cp.has_option('findVetoes','background-location'):
			self.add_var_opt("background-location",cp.get('findVetoes','background-location'))
		if cp.has_option('findVetoes','blind'):
			self.add_var_opt("blind",cp.get('findVetoes','blind'))
		self.output_cache = lal.CacheEntry(coincEvent.ifos, job.name.upper(), segments.segment(float(coincEvent.time), float(coincEvent.time)), "file://localhost/"+job.outputPath+'/DataProducts/'+oFilename)

		#IFO arg string
		myArgString=""
		if hasattr(coincEvent, "sngl_inspiral"):
			for sngl in coincEvent.sngl_inspiral.itervalues():
				myArgString=myArgString+"%s,"%sngl.ifo
		elif hasattr(coincEvent, "ifos_list"):
			for ifo in coincEvent.ifos_list:
				myArgString=myArgString+"%s,"%ifo
		myArgString=myArgString.rstrip(",")
		self.add_var_opt("ifo-list",myArgString)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())
		if not opts.no_findVetoes:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

#The class responsible for Node Object running the customFOM builder python
#script!
class customFOMPlotNode(pipeline.CondorDAGNode,FUNode):
	"""
	This is a node that corresponds with the job class to whip up
	custom FOMs.   In general each node will have one condor
	changed variable, which is t0 (gps) of trigger.
	"""
	defaults={"section":"customfoms",
		  "options":{"plot-windows":"14400,7200",
			     "ifo-list":"L1,H1,V1"}
			  }
	def __init__(self, dag, job, cp, opts, coincEvent):
		"""
		Takes in a coincEvent object and prepares figure request.
		"""
		self.__conditionalLoadDefaults__(customFOMPlotNode.defaults,cp)
		pipeline.CondorDAGNode.__init__(self,job)
		if cp.has_option('customfoms','plot-windows'):
			self.add_var_opt('plot-windows',cp.get('customfoms','plot-windows'))
		if cp.has_option('customfoms','ifo-list'):
			self.add_var_opt('ifo-list',cp.get('customfoms','ifo-list'))
		self.add_var_opt("gps-time",coincEvent.time)
		self.add_var_opt("verbose","")
		self.add_var_opt("output-path",job.outputPath+'/DataProducts/')
		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())
		#FIX ME: if the cluster is not CIT do not enable these jobs
		if not opts.no_findVetoes and "ligo.caltech.edu" in get_hostname():
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()
		
# EFFECTIVE DISTANCE RATIO NODE 
class effDRatioNode(pipeline.CondorDAGNode,FUNode):
	"""
	This Node class performs a parameter consistency check using the
	sites claiming to detect the trigger and the observed
	effective distance at each site. A command line example is
	below:
	followupRatioTest.py -R /archive/home/ctorres/public_html/DQstuff/ratioTest.pickle -iL1 -jH1 -kV1 -I10 -J10 -K5 -A 1 -B 1 -C 1.0001 -pmoinmoin -o mmTable.wiki
	"""
	defaults={"section":"effDRatio",
		  "options":{"output-file":"effDRatio.wiki",
			     "output-format":"moinmoin",
			     "snr-ratio-test":"/archive/home/ctorres/public_html/DQstuff/ratioTest.pickle"}
		  }
	def __init__(self, dag, job, cp, opts, coincEvent=None):
		"""
		"""
		self.__conditionalLoadDefaults__(effDRatioNode.defaults,cp)
		pipeline.CondorDAGNode.__init__(self,job)
		oFilename="%s-effDRatio_%s_%s.wiki"%(coincEvent.instruments,
						     coincEvent.ifos,
						     coincEvent.time)
		self.add_var_opt("output-file",job.outputPath+'/DataProducts/'+oFilename)
		self.add_var_opt("output-format",cp.get('effDRatio','output-format'))
		self.add_var_opt("snr-ratio-test",cp.get('effDRatio','snr-ratio-test'))
		#Grab Sngl propteries from Coinc object
		index=1
		for ifo,snglEvent in coincEvent.sngl_inspiral.items():
			if ifo in coincEvent.ifos:
				myIFO=snglEvent.ifo
				mySNR=snglEvent.snr
				myTIME=snglEvent.time
				self.add_var_opt("ifo%i"%(index),myIFO)
				self.add_var_opt("snr%i"%(index),mySNR)
				self.add_var_opt("time%i"%(index),myTIME)
				index=index+1
		for rIndex in range(index,3+1):
			self.add_var_opt("ifo%i"%(rIndex),None)
			self.add_var_opt("snr%i"%(rIndex),None)
			self.add_var_opt("time%i"%(rIndex),None)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if not opts.no_effectiveRatio:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

##############################################################################
# job class for producing the skymap

class lalapps_skyMapNode(pipeline.CondorDAGNode,FUNode):
	"""
	A C code for computing the sky map
	"""
	def __init__(self,dag,job,cp, opts, coinc, sngl_node_dict, p_nodes=[]):
		self.ifo_list = ["H1","L1","V1"]
		#self.already_added_ifo_list = []

		self.ra_res = job.ra_res
		self.dec_res = job.dec_res
		self.sample_rate = job.sample_rate
		pipeline.CondorDAGNode.__init__(self,job)

		# this program now gzips its files (otherwise they are really huge)
		self.output_file_name = job.outputPath + str(coinc.time) + ".txt.gz"
		self.add_var_opt("output-file",self.output_file_name)
		self.add_var_opt("ra-res",self.ra_res)
		self.add_var_opt("dec-res",self.dec_res)

		# Initialize input files
		for ifo in ['h1','h2','l1','v1']:
			self.add_var_opt(ifo+"-frame-file","none")
			self.add_var_opt(ifo+"-xml-file","none")
			self.add_var_opt(ifo+"-channel-name","none")

		# Overide the sample rate
		self.add_var_opt("sample-rate",coinc.get_sample_rate())

		# Now add the data we actually have
		for ifo, sngl in sngl_node_dict.items():
			self.add_var_opt(ifo.lower()+"-frame-file",sngl.output_file_name.replace(".xml",".gwf").strip(".gz"))
			self.add_var_opt(ifo.lower()+"-xml-file",sngl.output_file_name)
		for ifo, sngl in coinc.sngl_inspiral_coh.items():
			self.add_var_opt( "%s-channel-name" % (ifo.lower(),), "%s:CBC-CData_%d" % (ifo.upper(), int(sngl.row.event_id)) )

		self.output_cache = lal.CacheEntry("".join(coinc.instruments.split(",")), job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+os.path.abspath(self.output_file_name))

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		# Add parents and put this node in the dag
		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_skymap:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

# job class for producing the skymap
class pylal_skyPlotNode(pipeline.CondorDAGNode,FUNode):
	"""
A python code for plotting the sky map
	"""
	def __init__(self, dag, job,cp, opts, coinc, skyMapNode, p_nodes = []):

		pipeline.CondorDAGNode.__init__(self,job)
		#self.setupNode(job,True, dag.webPage.lastSection,page,None,None)
		self.add_var_opt("map-data-file",skyMapNode.output_file_name)
		self.add_var_opt("user-tag",str(coinc.time))
		self.add_var_opt("ifo-tag",coinc.ifos)
		self.add_var_opt("ifo-times",coinc.instruments)
		self.add_var_opt("ra-res",str(skyMapNode.ra_res))
		self.add_var_opt("dec-res",str(skyMapNode.dec_res))
		self.add_var_opt("stat-value", str(coinc.combined_far))
		# setup default arguments for plot jobs
		self.setupPlotNode(job)
		# if this is a software injection pass along the information to the
		# plotting code so that it can make a mark where the injection should have
		# been :)
		if coinc.sim:
			inj_ra = coinc.sim.longitude
			inj_dec = coinc.sim.latitude
			self.add_var_opt("injection-right-ascension",str(inj_ra))
			self.add_var_opt("injection-declination",str(inj_dec))

		self.output_file_name = skyMapNode.output_file_name.replace('.txt','.png')

		self.output_file_name = "%s-plot_inspiral_skymap_%s_%s-unspecified-gpstime.cache" % ( coinc.instruments, coinc.ifos, str(coinc.time))

		self.output_cache = lal.CacheEntry("".join(coinc.instruments.split(",")), job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+job.outputPath + '/' + self.output_file_name)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_skymap:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

class followUpChiaNode(inspiral.ChiaNode,FUNode):
	"""
A C code for computing the coherent inspiral statistic.
An example command line is:
lalapps_coherent_inspiral --segment-length 1048576 --dynamic-range-exponent 6.900000e+01 --low-frequency-cutoff 4.000000e+01 --bank-file H1H2-COHBANK_COHERENT_H1H2_PLAYGROUND-823269333-600.xml --sample-rate 4096 --cohsnr-threshold 5.500000e+00 --ifo-tag H1H2 --frame-type LSC-STRAIN --H1-framefile H1-INSPIRAL_COHERENT_H1H2_PLAYGROUND-823269286-2048.gwf --H2-framefile H2-INSPIRAL_COHERENT_H1H2_PLAYGROUND-823268952-2048.gwf --gps-end-time 823269933 --gps-start-time 823269333 --write-cohsnr --write-cohnullstat --write-cohphasediff --write-events --verbose
	"""
	#def __init__(self, chiaJob, procParams, trig, cp,opts,dag, trig_node, notrig_node ):

	#def __init__(self,job,trig,opts,dag,cp):
	def __init__(self, dag, job, cp, opts, coinc, inspiral_node_dict, chia_node =None, p_nodes = []):

		# the use of this class would require some reorganisation in fu_Condor.py
		# and webCondor.py in order to set up the jobs following the same scheme
		# as the way it is done for the Inspiral pipeline...
		pipeline.CondorDAGNode.__init__(self,job)
		pipeline.AnalysisNode.__init__(self)
		self.output_file_name = ""
		sngl = coinc.sngl_inspiral_coh.values()[0]

                user_tag = "COHERENT-"+str(coinc.time)

		# These come from inspiral process param tables
		self.add_var_opt( "segment-length", sngl.get_proc_param('segment-length') )
		self.add_var_opt( "dynamic-range-exponent",sngl.get_proc_param('dynamic-range-exponent') )
		self.add_var_opt( "low-frequency-cutoff", sngl.get_proc_param('low-frequency-cutoff') )
		self.add_var_opt("sample-rate", sngl.get_proc_param('sample-rate') )
		# come from config file		
		self.add_var_opt("cohsnr-threshold",cp.get('chia','cohsnr-threshold'))
		self.add_var_opt("ra-step",cp.get('chia','ra-step'))
		self.add_var_opt("dec-step",cp.get('chia','dec-step'))
		self.add_var_opt("cdata-length",1.0)
		self.add_var_opt("user-tag",user_tag)
		self.add_var_opt("ifo-tag",coinc.instruments.replace(',',''))
		self.add_var_opt("write-events","")
		self.add_var_opt("write-compress","")
                self.add_var_opt("maximize-over-chirp","")
		self.add_var_opt("followup","")
		# required by followUpChiaPlotNode
		if chia_node:
			self.add_var_opt("exttrig","")
			self.add_var_opt("chia-file",chia_node.output_file_name)
			self.add_var_opt("write-cohsnr","")
			self.add_var_opt("write-cohnullstat","")
			self.add_var_opt("write-h1h2nullstat","")
			self.add_var_opt("write-cohh1h2snr","")
			

                hLengthAnalyzed = 1

		#CHECK: needed here? self.setupNodeWeb(inspJob,False,None,None,None,dag.cache)
		#self.setupNodeWeb(job,False,None,None,None,dag.cache)
		self.add_var_opt("output-path",job.outputPath)

		# Here we define the trig-start-time and the trig-end-time;
		# The difference between these two times should be kept to 2s
		# Otherwise change the clustering window also
		self.start = int(coinc.time) - int(hLengthAnalyzed)
		self.end = int(coinc.time) + int(hLengthAnalyzed)

		self.add_var_opt("gps-start-time",self.start)
		self.add_var_opt("gps-end-time",self.end)


                if chia_node:
		        self.output_file_name = "%s/%s-CHIA_%s-%d-%d.xml.gz" % (job.outputPath, coinc.instruments.replace(',',''), user_tag, self.start, self.end-self.start )
                else:
                        self.output_file_name = "%s/%s-CHIA_%s-%d-%d-ALLSKY.xml.gz" % (job.outputPath, coinc.instruments.replace(',',''), user_tag, self.start, self.end-self.start )
		self.output_frame_file = "%s/%s-CHIA_%s-%d-%d.gwf" % (job.outputPath, coinc.instruments.replace(',',''), user_tag, self.start, self.end-self.start )
		self.netnull_output_frame_file = "%s/%s-CHIA_NULL_STAT_%s-%d-%d.gwf" % (job.outputPath, coinc.instruments.replace(',',''), user_tag, self.start, self.end-self.start )

 		self.h1h2null_output_frame_file = "%s/H1H2-CHIA_NULL_STAT_%s-%d-%d.gwf" % (job.outputPath, user_tag, self.start, self.end-self.start )
 		self.h1h2coh_output_frame_file = "%s/H1H2-CHIA_COHSNR_%s-%d-%d.gwf" % (job.outputPath, user_tag, self.start, self.end-self.start )


		self.output_cache = []

		self.output_cache.append(lal.CacheEntry("".join(coinc.instruments.split(",")), job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+os.path.abspath(self.output_file_name)))

		self.output_cache.append(lal.CacheEntry("".join(coinc.instruments.split(",")), job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+os.path.abspath(self.output_frame_file)))

		self.output_cache.append(lal.CacheEntry("".join(coinc.instruments.split(",")), job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+os.path.abspath(self.netnull_output_frame_file)))


                bankname = 'trig_bank/%s-COHBANK_FOLLOWUP_%s-%d-%d.xml.gz' % (coinc.instruments.replace(',',''), str(coinc.time), int(coinc.time) - int(hLengthAnalyzed), 2 * int(hLengthAnalyzed))
		bankFile = self.write_trigbank(coinc, bankname)
		self.set_bank(bankFile)

		arg_str = ''
		for ifo,sngl in inspiral_node_dict.items():
			arg_str += " --" + ifo.lower()+"-framefile " + sngl.output_frame_file

		self.add_var_arg(arg_str)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_chia:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()
			
	def write_trigbank(self, coinc, name):
		try:
			os.mkdir('trig_bank')
		except: pass
		xmldoc = ligolw.Document()
		xmldoc.appendChild(ligolw.LIGO_LW())

		# ADD A PROCESS TABLE
		process_params_table = lsctables.New(lsctables.ProcessParamsTable)
		xmldoc.childNodes[-1].appendChild(process_params_table)

		# ADD A SEARCH SUMMARY TABLE
		search_summary_table = lsctables.New(lsctables.SearchSummaryTable)
		xmldoc.childNodes[-1].appendChild(search_summary_table)
		row = search_summary_table.RowType()
		#FIXME THIS IS PROBABLY NOT HOW TO DO IT
		row.process_id = "process:process_id:0"
		row.shared_object = None
		row.lalwrapper_cvs_tag = None
		row.lal_cvs_tag = None
		row.comment = "Awesome"
		row.ifos = coinc.instruments
		#FIXME adjust for what omega actually analyzed 
		row.set_in(segments.segment(LIGOTimeGPS(self.start,0), LIGOTimeGPS(self.end,0)))
		row.set_out(segments.segment(LIGOTimeGPS(self.start,0), LIGOTimeGPS(self.end,0)))
		row.nevents = None
		row.nnodes = None
		search_summary_table.append(row)

		sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
		xmldoc.childNodes[-1].appendChild(sngl_inspiral_table)
		for ifo, sngl in coinc.sngl_inspiral_coh.items():
			sngl_inspiral_table.append(sngl.row)
		
		utils.write_filename(xmldoc, name, verbose=False, gz = True)
		return name

##############################################################################
# node class for plot snr chisq

class plotSNRCHISQNode(pipeline.CondorDAGNode,FUNode):
	"""
Runs an instance of a plotSNRCHISQ followup job
	"""
	def __init__(self, dag, job, cp, opts, sngl, coinc, sngl_node, p_nodes=[]):
	#def __init__(self,job,ifo,fileName,trig,page,dag,inspiralNode,opts,ifoString=None):
		"""
job = A CondorDAGJob that can run an instance of plotSNRCHISQ followup.
		"""
		pipeline.CondorDAGNode.__init__(self,job)
		self.output_file_name = ""
		self.add_var_opt("frame-file",sngl_node.output_frame_file)
		self.add_var_opt("inspiral-xml-file",sngl_node.output_file_name)

		duration = 2.0 # width of the time series to be displayed
		self.add_var_opt("plot-width",duration)

		self.add_var_opt("gps",sngl.time)
		self.add_var_opt("gps-start-time",sngl.time-duration*.5)
		self.add_var_opt("gps-end-time",sngl.time+duration*.5)

		self.add_var_opt("ifo-times", coinc.instruments)
		self.add_var_opt("ifo-tag", sngl.ifo)

		self.add_var_opt("user-tag","FOLLOWUP_PLOTSNRCHISQ_" + str(sngl.time))

		self.output_file_name = "%s-plotsnrchisq_pipe_%s_%s-%d-%d.cache" % ( coinc.instruments, sngl.ifo, "FOLLOWUP_PLOTSNRCHISQ_" + str(sngl.time), int(sngl.time-duration*.5), math.ceil(sngl.time+duration*.5) - int(sngl.time-duration*.5) )

		self.output_cache = lal.CacheEntry(sngl.ifo, job.name.upper(), segments.segment(float(sngl.time), float(sngl.time)), "file://localhost/"+job.outputPath + '/' + self.output_file_name)

		self.setupPlotNode(job)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

                for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_plotsnrchisq:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()


##############################################################################
# node class for plotting coherent inspiral search and null stat timeseries

class plotChiaNode(pipeline.CondorDAGNode, FUNode):
	"""
Runs an instance of a plotChia followup job
	"""

	def __init__(self, dag, job, cp, opts, coinc, chia_node, insp_node_dict, p_nodes=[]):
	#def __init__(self,job,chiaXmlFilePath,trig,cohireNode,dag,page,opts,cp):
		"""
job = A CondorDAGJob that can run an instance of plotChiaJob followup.
		"""

		pipeline.CondorDAGNode.__init__(self,job)
		self.output_file_name = ""
		user_tag = "PLOT_CHIA_" + str(coinc.time)
		self.add_var_opt("chiaXmlFile",chia_node.output_file_name)
		self.add_var_opt("chiaFrameFile",chia_node.output_frame_file)
		self.add_var_opt("cohH1H2SNRFrameFile",chia_node.h1h2coh_output_frame_file)
		self.add_var_opt("H1H2NullStatFrameFile",chia_node.h1h2null_output_frame_file)
		self.add_var_opt("cohNullStatFrameFile",chia_node.netnull_output_frame_file)
		self.add_var_opt("gps-start-time",int(coinc.time-1))
		self.add_var_opt("gps-end-time",int(coinc.time+1))
		self.add_var_opt("sample-rate",str(coinc.get_sample_rate()))
		self.add_var_opt("user-tag",user_tag)
		ifos = "".join(coinc.ifos.split(","))
		instruments = "".join(coinc.instruments.split(","))
		self.add_var_opt("ifo-tag",ifos)
		self.add_var_opt("ifo-times",instruments)
		self.setupPlotNode(job)

		self.output_file_name = "%s-plotchiatimeseries_%s_%s-%d-%d.cache" % ( instruments, ifos, "PLOT_CHIA_" + str(coinc.time), int(coinc.time-1), math.ceil(int(coinc.time+1)) - int(coinc.time-1) )

		self.output_cache = lal.CacheEntry(instruments, job.name.upper(), segments.segment(float(coinc.time), float(coinc.time)), "file://localhost/"+job.outputPath + '/' + self.output_file_name)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		for node in p_nodes:
			if node.validNode:
				self.add_parent(node)
		if not opts.no_chia:
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

		for ifo, insp in insp_node_dict.items():
			self.add_var_arg("--"+ifo.lower()+"-framefile "+ insp.output_frame_file)


##############################################################################
# node class for running the mcmc code

class mcmcNode(pipeline.CondorDAGNode, FUNode):
	"""
	Runs a MCMC job
	"""
	def __init__(self,dag,job,cp,opts,coinc,frame_cache_list,randomseed,p_nodes,ifo_string=None):
		pipeline.CondorDAGNode.__init__(self,job)

		time_margin = string.strip(cp.get('fu-mcmc','prior-coal-time-marg'))
		iterations = string.strip(cp.get('fu-mcmc','iterations'))
		tbefore = string.strip(cp.get('fu-mcmc','tbefore'))
		tafter = string.strip(cp.get('fu-mcmc','tafter'))
		#FIXME: priors on masses and distances should depend on the parameters of the trigger
		massmin = string.strip(cp.get('fu-mcmc','massmin'))
		massmax = string.strip(cp.get('fu-mcmc','massmax'))
		dist90 = string.strip(cp.get('fu-mcmc','dist90'))
		dist10 = string.strip(cp.get('fu-mcmc','dist10'))

		if ifo_string:
			IFOs = frozenset([ifo_string])
			self.ifonames = ifo_string
			sngl_insp_string = "sngl_inspiral"
		else:
			IFOs = coinc.ifos_set
			self.ifonames = coinc.instruments
			sngl_insp_string = "sngl_inspiral_coh"

		channelNames = ""
		chunk_end_list={}
		chunk_start_list={}
		for itf in IFOs:
			sngl = eval("coinc." + sngl_insp_string + "[\'" + itf + "\']")
			for row in sngl.process_params:
          			param = row.param.strip("-")
          			value = row.value
          			if param == 'channel-name':
					channel = value
				if param == 'gps-end-time':
					chunk_end = value
				if param == 'gps-start-time':
					chunk_start = value
			channelNames += channel + ","
			chunk_end_list[itf] = int(chunk_end)
			chunk_start_list[itf] = int(chunk_start)

		if len(IFOs) > 1:
			self.ifoRef = coinc.max_trigger_ifo()
		else:
			self.ifoRef = ifo_string

		self.add_var_opt("template",string.strip(cp.get('fu-mcmc','template')))
		self.add_var_opt("iterations",iterations)
		self.add_var_opt("randomseed",randomseed)
		self.add_var_opt("tcenter","%0.3f"%coinc.sngl_inspiral[self.ifoRef].time)
		self.add_var_opt("tbefore",tbefore)
		self.add_var_opt("tafter",tafter)

		tmin = coinc.sngl_inspiral[self.ifoRef].time - float(time_margin)
		tmax = coinc.sngl_inspiral[self.ifoRef].time + float(time_margin)
		self.add_var_opt("priorparameters","[" + massmin + "," + massmax + "," + str(tmin) + "," + str(tmax) + "," + dist90 + "," + dist10 + "]")

		param_mchirp = coinc.sngl_inspiral[self.ifoRef].row.mchirp
		param_eta = coinc.sngl_inspiral[self.ifoRef].row.eta
		param_distance = coinc.sngl_inspiral[self.ifoRef].row.eff_distance
		self.add_var_opt("guess","[" + str(param_mchirp) + "," + str(param_eta) + "," + str(coinc.sngl_inspiral[self.ifoRef].time) + "," + str(param_distance) + "]")

		cacheFiles = ""
		for frameCache in frame_cache_list:
			cacheFiles += frameCache + ","
		self.add_var_opt("cachefile","["+cacheFiles.strip(",")+"]")
		self.add_var_opt("filechannel","["+channelNames.strip(",")+"]")

		psdEstimateStart = ""
		psdEstimateEnd = ""
		for itf in IFOs:
			datainchunk_before = int(coinc.sngl_inspiral[self.ifoRef].time) - 75 - 64 - chunk_start_list[itf]
			datainchunk_after = chunk_end_list[itf] - 64 - int(coinc.sngl_inspiral[self.ifoRef].time) - 32
			if datainchunk_after > datainchunk_before:
				psdEstimateStart += str(int(coinc.sngl_inspiral[self.ifoRef].time) + 32) + ","
				psdEstimateEnd += str(chunk_end_list[itf] - 64) + ","
			else:
				psdEstimateStart += str(chunk_start_list[itf] + 64) + ","
				psdEstimateEnd += str(int(coinc.sngl_inspiral[self.ifoRef].time) - 75) + ","

		self.add_var_opt("psdestimatestart","["+psdEstimateStart.strip(",")+"]")
		self.add_var_opt("psdestimateend","["+psdEstimateEnd.strip(",")+"]")

		self.add_var_opt("importanceresample",10000)

		self.id = job.name.upper() + '-' + self.ifonames.replace(",","") + '-' + str(int(coinc.coinc_event_id)) + '_' + randomseed
		self.outputName = job.outputPath + '/' + self.id
		self.add_var_opt("outfilename",self.outputName)

		self.start_time = min(chunk_start_list.values())
		self.end_time = max(chunk_end_list.values())
		self.output_cache = []
		self.output_cache.append(lal.CacheEntry(self.ifonames.replace(",",""), job.name.upper(), segments.segment(self.start_time,self.end_time), "file://localhost/"+self.outputName+".csv"))

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if opts.enable_bayesian:
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()

##############################################################################
# node class for running the spinmcmc code

class spinmcmcNode(pipeline.CondorDAGNode, FUNode):
	"""
	Runs a SPIN MCMC job
	"""
	def __init__(self,dag,job,cp,opts,coinc,frame_cache_list,p_nodes):
		pipeline.CondorDAGNode.__init__(self,job)

		iterations = string.strip(cp.get('fu-spinmcmc','iterations'))
		tbefore = string.strip(cp.get('fu-spinmcmc','tbefore'))
		tafter = string.strip(cp.get('fu-spinmcmc','tafter'))

		IFOs = coinc.ifos_set
		self.ifonames = coinc.instruments
		sngl_insp_string = "sngl_inspiral_coh"

		channelNames = ""
		ifoString = ""
		chunk_end_list={}
		chunk_start_list={}
		for itf in IFOs:
			sngl = eval("coinc." + sngl_insp_string + "[\'" + itf + "\']")
			for row in sngl.process_params:
				param = row.param.strip("-")
				value = row.value
				if param == 'channel-name':
					channel = value
				if param == 'gps-end-time':
					chunk_end = value
				if param == 'gps-start-time':
					chunk_start = value
			channelNames += channel + ","
			ifoString += itf + ","
			chunk_end_list[itf] = int(chunk_end)
			chunk_start_list[itf] = int(chunk_start)

		ifoString = ifoString.replace("H1","1")
		ifoString = ifoString.replace("L1","2")
		ifoString = ifoString.replace("V1","3")
		self.add_var_opt("network","["+ifoString.strip(",")+"]")

		self.ifoRef = coinc.max_trigger_ifo()

		self.add_var_opt("nIter",iterations)
                self.add_var_opt("tc","%0.3f"%coinc.sngl_inspiral[self.ifoRef].time)
		self.add_var_opt("beforetc",tbefore)
		self.add_var_opt("aftertc",tafter)

		param_mchirp = coinc.sngl_inspiral[self.ifoRef].row.mchirp
		param_eta = coinc.sngl_inspiral[self.ifoRef].row.eta
		param_distance = coinc.sngl_inspiral[self.ifoRef].row.eff_distance

		self.add_var_opt("mChirp",param_mchirp)
		self.add_var_opt("eta",param_eta)
		self.add_var_opt("dist",param_distance)

		cacheFiles = ""
		for frameCache in frame_cache_list:
			cacheFiles += frameCache + ","
		self.add_var_opt("cache","["+cacheFiles.strip(",")+"]")
		self.add_var_opt("channel","["+channelNames.strip(",")+"]")

#FIX ME: FOR NOW WE ARE LETTING THE CODE CHOSING AUTOMATICALLY THE DATA SEGMENT ON WHICH THE PSD IS COMPUTED
#		psdEstimateStart = ""
#		psdEstimateEnd = ""
#		for itf in IFOs:
#			datainchunk_before = int(coinc.sngl_inspiral[self.ifoRef].time) - 75 - 64 - chunk_start_list[itf]
#			datainchunk_after = chunk_end_list[itf] - 64 - int(coinc.sngl_inspiral[self.ifoRef].time) - 32
#			if datainchunk_after > datainchunk_before:
#				psdEstimateStart += str(int(coinc.sngl_inspiral[self.ifoRef].time) + 32) + ","
#				psdEstimateEnd += str(chunk_end_list[itf] - 64) + ","
#			else:
#				psdEstimateStart += str(chunk_start_list[itf] + 64) + ","
#				psdEstimateEnd += str(int(coinc.sngl_inspiral[self.ifoRef].time) - 75) + ","
#		self.add_var_opt("psdestimatestart","["+psdEstimateStart.strip(",")+"]")
#		self.add_var_opt("psdestimateend","["+psdEstimateEnd.strip(",")+"]")


		self.id = job.name.upper() + '-' + self.ifonames.replace(",","") + '-' + str(int(coinc.coinc_event_id))
		#FIXME: WHAT IS THE ACTUAL OUTPUT FILE?
 		self.outputName = job.outputPath + '/' + self.id
		self.add_var_opt("outputPath",job.outputPath)

                self.start_time = min(chunk_start_list.values())
                self.end_time = max(chunk_end_list.values())
                self.output_cache = []
                self.output_cache.append(lal.CacheEntry(self.ifonames.replace(",",""), job.name.upper(), segments.segment(self.start_time,self.end_time), "file://localhost/"+self.outputName))
		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if opts.enable_bayesian:
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()


##############################################################################
# node class for running the plotting of the mcmc results

class plotmcmcNode(pipeline.CondorDAGNode, FUNode):
        """
        Runs a plotmcmc job
        """
	def __init__(self,job,coinc,cp,opts,dag,ifo,ifonames,p_nodes):
		pipeline.CondorDAGNode.__init__(self,job)

		if job.tag_base=="sngl":
			sngl_insp_string = "sngl_inspiral"
		else:
			sngl_insp_string = "sngl_inspiral_coh"

		sngl = eval("coinc." + sngl_insp_string + "[\'" + ifo + "\']")

		if cp.has_option('fu-plotmcmc','burnin'):
			burnin = string.strip(cp.get('fu-plotmcmc','burnin'))
			if burnin.strip():
				self.add_var_opt("burnin",burnin)

		plot_routine = string.strip(cp.get('fu-plotmcmc','plot_routine'))
		executable = string.strip(cp.get('fu-plotmcmc','executable'))

		#FIXME: add a if statement to treat differently the injections. Reference values for injections should be the injected params.
		gps = sngl.time
		mchirp = sngl.row.mchirp
		eta = sngl.row.eta
		distance = sngl.row.eff_distance
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
		for node in p_nodes:
			mcmcfilelist += node.outputName + '.csv,'
		self.add_var_opt("mcmc-file",mcmcfilelist.strip(','))

		self.id = job.name.upper() + '-' + ifonames.replace(",","") + '-' + str(int(coinc.coinc_event_id))
		self.add_var_opt("identity",self.id)

		self.add_var_opt("output-path",job.outputPath)
		self.output_cache = lal.CacheEntry(ifonames.replace(",",""), job.name.upper(), segments.segment(p_nodes[0].start_time,p_nodes[0].end_time), "file://localhost/"+job.outputPath+"/"+self.id)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if opts.enable_bayesian:
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()


##############################################################################
# node class for running the plotting of the spin mcmc results

class plotspinmcmcNode(pipeline.CondorDAGNode, FUNode):

	def __init__(self,job,coinc,cp,opts,dag,ifo,ifonames,p_nodes):
		pipeline.CondorDAGNode.__init__(self,job)

		sngl_insp_string = "sngl_inspiral_coh"

		sngl = eval("coinc." + sngl_insp_string + "[\'" + ifo + "\']")

		plot_routine = string.strip(cp.get('fu-plotmcmc','plot_routine'))
		executable = string.strip(cp.get('fu-plotmcmc','executable'))

		#FIXME: add a if statement to treat differently the injections. Reference values for injections should be the injected params.
		gps = sngl.time
		mchirp = sngl.row.mchirp
		eta = sngl.row.eta
		distance = sngl.row.eff_distance
		#FIXME: HOW TO SETUP CORRECTLY THE FOLLOWING PARAMETERS?
		phi = "0.0"
		a_spin1 = "0.5"
		cs_th_sp1 = "0.1"
		phi_spin1 = "0.4"
		a_spin2 = "0.5"
		cs_th_sp2 = "0.5"
		phi_spin2 = "0.3"

		self.add_var_opt("plot-routine",plot_routine)
		self.add_var_opt("executable",executable)
		self.add_var_opt("reference-time",gps)
		self.add_var_opt("reference-mchirp",mchirp)
		self.add_var_opt("reference-eta",eta)
		self.add_var_opt("reference-distance",distance)
		self.add_var_opt("reference-phi",phi)
		self.add_var_opt("reference-a_spin1",a_spin1)
		self.add_var_opt("reference-a_spin2",a_spin2)
		self.add_var_opt("reference-phi_spin1",phi_spin1)
		self.add_var_opt("reference-phi_spin2",phi_spin2)
		self.add_var_opt("reference-cs_th_sp1",cs_th_sp1)
		self.add_var_opt("reference-cs_th_sp2",cs_th_sp2)

		# get the list of MCMC .txt files to be used as input
		mcmcfilelist = ""
		for node in p_nodes:
			mcmcfilelist += node.outputName
		self.add_var_opt("mcmc-file",mcmcfilelist.strip(','))

		self.id = job.name.upper() + '-' + ifonames.replace(",","") + '-' + str(int(coinc.coinc_event_id))
		self.add_var_opt("identity",self.id)

		self.add_var_opt("output-path",job.outputPath)
		self.output_cache = lal.CacheEntry(ifonames.replace(",",""), job.name.upper(), segments.segment(p_nodes[0].start_time,p_nodes[0].end_time), "file://localhost/"+job.outputPath+"/"+self.id)

		if not opts.disable_dag_categories:
			self.set_category(job.name.lower())

		if opts.enable_bayesian:
			for node in p_nodes:
				if node.validNode:
					self.add_parent(node)
			dag.add_node(self)
			self.validate()
		else:
			self.invalidate()


##############################################################################
###### CONDOR DAG THINGY #####################################################
##############################################################################

class followUpDAG(pipeline.CondorDAG):
	def __init__(self, config_file, cp, opts):
		log_path = cp.get('fu-output','log-path').strip()
		self.basename = re.sub(r'\.ini',r'', os.path.split(config_file)[1])
		tempfile.tempdir = log_path
		tempfile.template = self.basename + '.dag.log.'
		logfile = tempfile.mktemp()
		fh = open( logfile, "w" )
		fh.close()
		pipeline.CondorDAG.__init__(self,logfile)
		self.set_dag_file(self.basename)
		self.jobsDict = {}
		self.node_id = 0
		self.output_cache = []
		if not opts.disable_dag_categories:
			for cp_opt in cp.options('condor-max-jobs'):
					self.add_maxjobs_category(cp_opt,cp.getint('condor-max-jobs',cp_opt))

	def add_node(self,node):
		self.node_id += 1
		node.add_macro("macroid", self.node_id)
		pipeline.CondorDAG.add_node(self, node)
		try: self.output_cache.extend(node.output_cache)
		except:
			try: self.output_cache.append(node.output_cache)
			except: pass

	def write_all(self):
		self.write_sub_files()
		self.write_dag()
		self.write_script()
		self.write_output_cache()

	def write_output_cache(self):
		f = open(self.basename+".cache",'w')
		for c in self.output_cache:
			f.write(str(c)+'\n')
		f.close()

###############################################################################
###### CONFIG PARSER WRAPPING #################################################
###############################################################################
class create_default_config(object):
	def __init__(self, configfile=None):
		cp = ConfigParser.ConfigParser()
		self.cp = cp
		self.time_now = "_".join([str(i) for i in time_method.gmtime()[0:6]])
		self.ini_file=self.time_now + ".ini"
		home_base = home_dirs()
		
		# CONDOR SECTION NEEDED BY THINGS IN INSPIRAL.PY
		cp.add_section("condor")
		cp.set("condor","datafind",self.which("ligo_data_find"))
		cp.set("condor","inspiral",self.which("lalapps_inspiral"))
		cp.set("condor","chia", self.which("lalapps_coherent_inspiral"))
		cp.set("condor","universe","standard")
		# SECTIONS TO SHUT UP WARNINGS
		cp.add_section("inspiral")
		cp.add_section("data")	 
	
		# DATAFIND SECTION
		cp.add_section("datafind")

		# FU-CONDOR SECTION
		cp.add_section("fu-condor")
		cp.set("fu-condor","plotsnrchisq",self.which("plotsnrchisq_pipe"))
		cp.set("fu-condor","lalapps_skymap",self.which("lalapps_skymap"))
		cp.set("fu-condor","pylal_skyPlotJob",self.which("pylal_plot_inspiral_skymap"))
		cp.set("fu-condor","datafind",self.which("ligo_data_find"))
		cp.set("fu-condor","convertcache",self.which("convertlalcache.pl"))
		cp.set("fu-condor","chia", self.which("lalapps_coherent_inspiral"))
		cp.set("fu-condor","plotchiatimeseries", self.which("plotchiatimeseries"))
                cp.set("fu-condor","effDRatio", self.which("followupRatioTest.py"))
                cp.set("fu-condor","vetoflags", self.which("followupQueryVeto.py"))
                cp.set("fu-condor","customfom", self.which("followupCustomFOM.py"))
                cp.set("fu-condor","dqflags", self.which("followupQueryDQ.py"))
		cp.set("fu-condor","mcmc", self.which("lalapps_followupMcmc"))
		cp.set("fu-condor","spinmcmc", self.which("lalapps_spinspiral"))
		cp.set("fu-condor","plotmcmc", self.which("plotmcmc.py"))
		cp.set("fu-condor","plotspinmcmc", self.which("plotspinmcmc.py"))
		#FIXME SET THIS TO SOMETHING THAT WORKS
		#cp.set("fu-condor","qscan",home_base+"/romain/opt/omega/omega_r2062_glnxa64_binary/bin/wpipeline")
		self.set_qscan_executable()
		cp.set("fu-condor","analyseQscan", self.which("analyseQscan.py"))
		cp.set("fu-condor","makeCheckListWiki",self.which("makeCheckListWiki.py"))
		cp.set("fu-condor","lalapps_followup_page",self.which("lalapps_followup_page"))
		# makechecklistwiki SECTION
		cp.add_section("makeCheckListWiki")
		cp.set("makeCheckListWiki","universe","local")
		cp.set("makeCheckListWiki","location",os.getcwd())
		#Store full abs path in ini file!
		cp.set("makeCheckListWiki","ini-file",os.path.abspath(self.ini_file))
		
		# fu-q-hoft-datafind SECTION
		cp.add_section("fu-q-hoft-datafind")
		for ifo in ["H1","H2","L1","V1"]:
			cp.set("fu-q-hoft-datafind",ifo+"-search-time-range","128")

		# fu-q-rds-datafind SECTION
		cp.add_section("fu-q-rds-datafind")
		for ifo in ["H1","H2","L1"]:
			cp.set("fu-q-rds-datafind",ifo+"-search-time-range","1024")
		cp.set("fu-q-rds-datafind","V1-search-time-range","2048")
		
		# fu-fg-ht-qscan SECTION
		cp.add_section("fu-fg-ht-qscan")
		for config in ["H1config","H2config","L1config","V1config"]:
			cp.set("fu-fg-ht-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'hoft') + ".txt","QSCAN CONFIG"))

		# fu-fg-rds-qscan SECTION
		cp.add_section("fu-fg-rds-qscan")
		for config in ["H1config","H2config","L1config"]:
			cp.set("fu-fg-rds-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'rds') + ".txt","QSCAN CONFIG"))
		cp.set("fu-fg-rds-qscan","V1config","/storage/gpfs_virgo3/virgo/omega/configurations/s6_foreground_V1-raw-cbc.txt")

		# fu-fg-seismic-qscan SECTION
		cp.add_section("fu-fg-seismic-qscan")
		for config in ["H1config","H2config","L1config"]:
			cp.set("fu-fg-seismic-qscan",config,self.__find_config("s5_foreground_" + self.__config_name(config[:2],'seismic') + ".txt","QSCAN CONFIG"))
		cp.set("fu-fg-seismic-qscan","V1config","/storage/gpfs_virgo3/virgo/omega/configurations/s6_foreground_V1-raw-seismic-cbc.txt")

		# fu-analyse-qscan SECTION
		cp.add_section("fu-analyse-qscan")
		cp.set("fu-analyse-qscan","generate-qscan-xml","")
		cp.set("fu-analyse-qscan","z-threshold","0.0")
		cp.set("fu-analyse-qscan","z-min","0.0")
		cp.set("fu-analyse-qscan","z-max","30.0")
		cp.set("fu-analyse-qscan","z-bins","60")
		cp.set("fu-analyse-qscan","rds-dt-min","-0.6")
		cp.set("fu-analyse-qscan","rds-dt-max","0.6")
		cp.set("fu-analyse-qscan","ht-dt-min","-0.6")
		cp.set("fu-analyse-qscan","ht-dt-max","0.6")
		cp.set("fu-analyse-qscan","seis-rds-dt-min","-4.2")
		cp.set("fu-analyse-qscan","seis-rds-dt-max","4.2")
		cp.set("fu-analyse-qscan","dt-bins","120")
		cp.set("fu-analyse-qscan","plot-dt-distribution","")
		cp.set("fu-analyse-qscan","plot-z-scattered","")
		cp.set("fu-analyse-qscan","plot-z-distribution","")

		# FU-SKYMAP SECTION
		cp.add_section("fu-skymap")
		cp.set("fu-skymap","ra-res","1024")
		cp.set("fu-skymap","dec-res","512")
		cp.set("fu-skymap","sample-rate","4096")

		# FU-OUTPUT SECTION
		cp.add_section("fu-output")
		cp.set("fu-output","log-path",self.log_path())
		cp.set("fu-output","output-dir",self.web_dir())
		cp.set("fu-output","web-url", self.web_url())

		# CHIA SECTION
		cp.add_section("chia")
		cp.set('chia','cohsnr-threshold', "1")
		cp.set('chia','ra-step', "1")
		cp.set('chia','dec-step', "1")
		cp.set('chia','numCohTrigs', "2000")
		cp.set('chia', 'sample-rate', "4096")

		# EFFECTIVE DIST RATIO TEST SECTION
		cp.add_section("effDRatio")
		cp.set('effDRatio','snr-ratio-test',self.__find_config("ratioTest.pickle","RATIO TEST PICKLE"))

		# FU-MCMC SECTION
		cp.add_section("fu-mcmc")
		cp.set("fu-mcmc","chain_nb","6")
		cp.set("fu-mcmc","prior-coal-time-marg","0.050")
		cp.set("fu-mcmc","iterations","1000000")
		cp.set("fu-mcmc","tbefore","30")
		cp.set("fu-mcmc","tafter","1")
		cp.set("fu-mcmc","template","20SP")
		# FIXME: mass and distance priors should be choosen as a function of the parameters of the CBC trigger...
		cp.set("fu-mcmc","massmin","1.0")
		cp.set("fu-mcmc","massmax","15.0")
		cp.set("fu-mcmc","dist90","40.0")
		cp.set("fu-mcmc","dist10","80.0")

		# FU-PLOTMCMC SECTION
		cp.add_section("fu-plotmcmc")
		cp.set("fu-plotmcmc","plot_routine",self.__find_routine("mcmcsummary.R","R SCRIPT FOR MCMC PLOTS"))
		cp.set("fu-plotmcmc","executable","/usr/bin/R")

		# FU-SPINMCMC SECTION
                cp.add_section("fu-spinmcmc")
		cp.set("fu-spinmcmc","iterations","1000000")
		cp.set("fu-spinmcmc","tbefore","30")
		cp.set("fu-spinmcmc","tafter","1")

		# REMOTE JOBS SECTION
		cp.add_section("fu-remote-jobs")
		remoteIfos,remoteJobs = self.get_remote_jobs()
		cp.set('fu-remote-jobs','remote-ifos',remoteIfos)
		cp.set('fu-remote-jobs','remote-jobs',remoteJobs)

		# CONDOR MAX JOBS SECTION
		cp.add_section("condor-max-jobs")
		cp.set("condor-max-jobs","remoteScan_full_data_FG_RDS.sh_FG_RDS_full_data","20")
		cp.set("condor-max-jobs","remoteScan_full_data_FG_SEIS_RDS.sh_FG_SEIS_RDS_full_data","20")
		cp.set("condor-max-jobs","remoteScan_playground_FG_RDS.sh_FG_RDS_playground","20")
		cp.set("condor-max-jobs","remoteScan_playground_FG_SEIS_RDS.sh_FG_SEIS_RDS_playground","20")
		cp.set("condor-max-jobs","remoteScan_time_slides_FG_RDS.sh_FG_RDS_time_slides","20")
                cp.set("condor-max-jobs","remoteScan_time_slides_FG_SEIS_RDS.sh_FG_SEIS_RDS_time_slides","20")
		cp.set("condor-max-jobs","remoteDatafind_full_data_Q_RDS.sh_Q_RDS_full_data","10")
		cp.set("condor-max-jobs","remoteDatafind_full_data_Q_RDS.sh_Q_RDS_playground","10")
		cp.set("condor-max-jobs","remoteDatafind_full_data_Q_RDS.sh_Q_RDS_time_slides","10")
		cp.set("condor-max-jobs","ligo_data_find_HT_full_data","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_HT_full_data","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_RDS_full_data","3")
		cp.set("condor-max-jobs","ligo_data_find_HT_playground","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_HT_playground","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_RDS_playground","3")
		cp.set("condor-max-jobs","ligo_data_find_HT_time_slides","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_HT_time_slides","3")
		cp.set("condor-max-jobs","ligo_data_find_Q_RDS_time_slides","3")
		cp.set("condor-max-jobs","lalapps_followupmcmc_sngl_full_data","20")
		cp.set("condor-max-jobs","lalapps_followupmcmc_sngl_playground","20")
		cp.set("condor-max-jobs","lalapps_followupmcmc_sngl_time_slides","20")
		cp.set("condor-max-jobs","lalapps_followupmcmc_coh_full_data","20")
		cp.set("condor-max-jobs","lalapps_followupmcmc_coh_playground","20")
		cp.set("condor-max-jobs","lalapps_followupmcmc_coh_time_slides","20")
		cp.set("condor-max-jobs","lalapps_spinspiral_coh_full_data","20")
		cp.set("condor-max-jobs","lalapps_spinspiral_coh_playground","20")
		cp.set("condor-max-jobs","lalapps_spinspiral_coh_time_slides","20")

		# Following comments relate to default options
		# Generate by FUNode.__conditionalLoadDefaults__ method
		#findFlagsNode
		#findVetosNode
		
		# if we have an ini file override the options
		if configfile:
			user_cp = ConfigParser.ConfigParser()
			user_cp.read(configfile)
		else:
			# otherwise see if a file with the standard ini file exists in the directory, the user probably intends to use it
			try: 
				user_cp = ConfigParser.ConfigParser()
				user_cp.read('followup_pipe.ini')
			except: pass
		# override the default options
		if user_cp: self.overwrite_config(user_cp,cp)

	def write(self):
		self.get_cp().write(open(self.ini_file,"w"))

	def get_cp(self):
		return self.cp

	def set_qscan_executable(self):
		host = get_hostname()
		if 'phy.syr.edu' in host:
			self.cp.set("fu-condor","qscan",home_dirs()+"/rgouaty/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")
		else:
			self.cp.set("fu-condor","qscan",home_dirs()+"/romain/opt/omega/omega_r3270_glnxa64_binary/bin/wpipeline")		

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

	def __find_routine(self,script,description):
		path = self.which('lalapps_inspiral')
                if path: path = os.path.split(path)[0]
		else:
			print >>sys.stderr, "COULD NOT FIND " + description + " FILE %s IN %s, ABORTING" % (script, path)
			raise ValueError
                        sys.exit(1)
		out = path.replace('bin','share/lalapps') + '/' + script
		if not os.path.isfile(out):
			print >>sys.stderr, "COULD NOT FIND " + description + " FILE %s IN %s, ABORTING" % (script, out)
			raise ValueError
			sys.exit(1)
		return out

	def web_dir(self):
		host = get_hostname()
		#FIXME add more hosts as you need them
		if 'caltech.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'phy.syr.edu' in host: return os.path.abspath(os.environ['HOME']) + '/public_html/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return os.path.abspath(os.environ['HOME']) + '/WWW/LSC/followups/' + self.time_now
		print sys.stderr, "WARNING: could not find web directory, returning empty string"
		return ''

	def web_url(self):
		host = get_hostname()
		#FIXME add more hosts as you need them
		if 'ligo.caltech.edu' in host: return "https://ldas-jobs.ligo.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-la.caltech.edu' in host: return "https://ldas-jobs.ligo-la.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'ligo-wa.caltech.edu' in host: return "https://ldas-jobs.ligo-wa.caltech.edu/~" +os.environ['USER'] + '/followups/' + self.time_now
		if 'phys.uwm.edu' in host: return "https://ldas-jobs.phys.uwm.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'phy.syr.edu' in host: return "https://sugar-jobs.phy.syr.edu/~" + os.environ['USER'] + '/followups/' + self.time_now
		if 'aei.uni-hannover.de' in host: return "https://atlas3.atlas.aei.uni-hannover.de/~" + os.environ['USER'] + '/LSC/followups/' + self.time_now
		print sys.stderr, "WARNING: could not find web server, returning empty string"
		return ''

	def get_remote_jobs(self):
		host = get_hostname()
                #FIXME add more hosts as you need them
		if 'ligo.caltech.edu' or 'ligo-la.caltech.edu' or 'ligo-wa.caltech.edu' or 'phys.uwm.edu' or 'aei.uni-hannover.de' or 'phy.syr.edu' in host:
			remote_ifos = "V1"
			remote_jobs = "ligo_data_find_Q_RDS_full_data,wpipeline_FG_RDS_full_data,wpipeline_FG_SEIS_RDS_full_data,ligo_data_find_Q_RDS_playground,wpipeline_FG_RDS_playground,wpipeline_FG_SEIS_RDS_playground,ligo_data_find_Q_RDS_gps_only,wpipeline_FG_RDS_gps_only,wpipeline_FG_SEIS_RDS_gps_only,ligo_data_find_Q_RDS_time_slides,wpipeline_FG_RDS_time_slides,wpipeline_FG_SEIS_RDS_time_slides"
			return remote_ifos, remote_jobs
		return '', ''

	def log_path(self):
		host = get_hostname()
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

	def overwrite_config(self,config,cp):
		for section in config.sections():
			if not cp.has_section(section): cp.add_section(section)
			for option in config.options(section):
				cp.set(section,option,config.get(section,option))



#A get links to ifo FOMS[1,2,3]
def getFOMLinks(gpsTime=int(0),ifo=("default")):
	"""
	Simple method returns a list of links to FOMs ordered by FOM #
	The list is 2D ie:
	[['ifo,shift',LINKtoImage,LinktoThumb],['ifo,shift',LinktoImage,LinkToThumb]...]
	images marked [Eve,Owl,Day] via [p3,p2,p1] in filenames
	this methd only for S6 and later
	IFO naming start dates:
	There were three naming conventions mixed, then p1,p2,p3 and lastly Day,Eve,Owl
	LHO: 20090724 :: 932428815
	LLO: 20090708 :: 931046415
	It appears that the filenames are labeled by local times not
	utc??? We need to confirm this for this method CVT Fri-Jan-29-2010:201001291523 
	"""
	urls={
		"DEFAULT":"http://www.ligo.caltech.edu/~pshawhan/scilinks.html",
		"V1":"http://wwwcascina.virgo.infn.it/DetectorOperations/index.htm",
		"L1":"https://llocds.ligo-la.caltech.edu/scirun/S6/robofom/%s/%s%s_FOM%i%s.gif",
		"H1":"http://lhocds.ligo-wa.caltech.edu/scirun/S6/robofom/%s/%s%s_FOM%i%s.gif",
		"H2":"http://lhocds.ligo-wa.caltech.edu/scirun/S6/robofom/%s/%s%s_FOM%i%s.gif"
		}
	ifoTag=ifo.upper()
	shiftDuration=8;
	#Give the IFO and shift start hour as integer
	shiftStandardTime={'L1':{'day':14,'eve':22,'owl':6},
			   'H1':{'day':16,'eve':0,'owl':8},
			   'H2':{'day':16,'eve':0,'owl':8},
			   'V1':{'day':6,'eve':14,'owl':22}}
	shiftOrder=['day','eve','owl']
	shiftLabel={'day':'p1','eve':'p3','owl':'p2'}
	outputURLs=list()
	if ((ifo==None) or (gpsTime==None)):
		sys.stdout.write("getFOMLinks called incorrectly \
using default opts instead!\n")
		return [urls['DEFAULT']]
	outputURLs=[]
	#Just return immediately if V1 encoutered (HACK)
	if ifo.__contains__("V1"):
		return([['V1',urls[ifoTag],'']])
	#
	if shiftStandardTime.keys().__contains__(ifoTag):
		#Determine shift times n-1,n,n+1
		tOffset=3600*shiftDuration
		for thisTime in \
		[gpsTime-tOffset,gpsTime,gpsTime+tOffset]:
			Y,M,D,h,m,s,junk0,junk1,junk2=xlaldate.XLALGPSToUTC(LIGOTimeGPS(int(thisTime)))
			#Get shift label
			shiftString=''
			humanShiftLabel=''
			for shift,start in shiftStandardTime[ifoTag].iteritems():
				hours=[x%24 for x in range(start,start+shiftDuration)]
				if hours.__contains__(int(h)):
					shiftString=shiftLabel[shift]
					humanShiftLabel=shift
			#Need to back up one day in some cases
			#when the shift of interest started -1 day
			#before trigger time
			if (0 in hours) and (hours[0]<h):
				D=D-1;
			if D<1:
				D=1
				M=M-1
			if M<1:
				M=1
				Y=Y-1
			#Create txt string
			tString="%s%s%s"%(str(Y).zfill(4),str(M).zfill(2),str(D).zfill(2))
			if ('V1').__contains__(ifoTag):
				outputURLs.append(['V1',urls[ifoTag],''])
			else:
				sL=shiftString
				for fL in [1,2,3]:
					outputURLs.append(["%s,%s"%(ifoTag,humanShiftLabel),
							   urls[ifoTag]%(tString,tString,sL,fL,""),
							   urls[ifoTag]%(tString,tString,sL,fL,"Thumb")
						   ])
		return outputURLs
	else:
		return [urls['DEFAULT']]

#A simple method to convert GPS time to human readable for for
#checklist
def gpsTimeToReadableDate(gpsTime=float(0)):
	"""
	Pass in int form of gps time.
	"""
	lGTime=LIGOTimeGPS(int(gpsTime))
	Y,M,D,h,m,s,junk0,junk1,junk2=xlaldate.XLALGPSToUTC(lGTime)
	timeStamp=str("%s-%s-%s  %s:%s:%s UTC"%(str(Y).zfill(4),
						str(M).zfill(2),
						str(D).zfill(2),
						str(h).zfill(2),
						str(m).zfill(2),
						str(s).zfill(2)))
	return timeStamp

#A loose method to retrieve the iLog url given a integer for of
#GPStimeA
def getiLogURL(time=None,ifo=None):
	"""
	This method returns a URL string to point you to ilog day page for
	specified IFO and GPStime. Valid IFO labels are V1, L1, H1 or H2.
	"""
	time=int(float(time))
	dateString="%s/%s/%s"
	urls={
		'default':"http://www.ligo.caltech.edu/~pshawhan/scilinks.html",
		'V1':"https://pub3.ego-gw.it/logbook/index.php?area=logbook&ref=search&datefrom=%s&dateto=%s",
		'L1':"http://ilog.ligo-la.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to=",
		'H1':"http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to=",
		'H2':"http://ilog.ligo-wa.caltech.edu/ilog/pub/ilog.cgi?task=view&date_to_view=%s&group=detector&keywords_to_highlight=&text_to_highlight=&anchor_to_scroll_to="
		}
	outputURL=urls['default']
	if ((ifo==None) or (time==None)):
		return urls['default']
	gpsTime=LIGOTimeGPS(time)
	Y,M,D,h,m,s,junk0,junk1,junk2=xlaldate.XLALGPSToUTC(gpsTime)
	gpsStamp=dateString%(str(M).zfill(2),str(D).zfill(2),str(Y).zfill(4))
	if ('H1','H2','L1').__contains__(ifo.upper()):
		outputURL=urls[ifo.upper()]%gpsStamp
	if ('V1').__contains__(ifo.upper()):
		gpsTimePO=LIGOTimeGPS(time+(24*3600))		
		Y2,M2,D2,h2,m2,s2,junk0,junk1,junk2=xlaldate.XLALGPSToUTC(gpsTimePO)
		gpsStampPlusOne=dateString%(str(M2).zfill(2),str(D2).zfill(2),str(Y2).zfill(4))
		outputURL=urls[ifo.upper()]%(gpsStamp,gpsStampPlusOne)
	return outputURL

#End def getiLogURL

#Maps image paths to URLS for makeCheckListWiki.py
class filenameToURLMapper(object):
	  """
	  """
	  def __init__(self,publicationDirectory=None,publicationURL=None,verbose=False):
		    protocolTag="@PROTO@/"
		    self.verbose=verbose
		    self.validProtocols=["http://","https://"]
		    givenProtocol=""
		    if publicationDirectory == None or\
			   publicationURL == None:
			    sys.stderr.write("Error: Initializing filenameToURLMappe instance \
			    with None types.\n")
		    self.pDIR=publicationDirectory
		    self.pURL=publicationURL
		    for protocolCheck in self.validProtocols:
			if publicationDirectory.lower().startswith(protocolCheck):
				self.pDIR=publicationDirectory
				self.pURL=publicationURL
				raise Warning,"object initialized with publication directory and publication URL reversed\n"
		    for protocolCheck in self.validProtocols:
			    if self.pURL.lower().startswith(protocolCheck):
				    self.pURL="%s"%(self.pURL.replace(protocolCheck,protocolTag))
				    givenProtocol=protocolCheck
		    pd=self.pDIR.lstrip(os.path.sep).split(os.path.sep)
		    pu=self.pURL.split(os.path.sep)
		    self.pURL=publicationURL
		    pd.reverse()
		    pu.reverse()
		    cStringList=list()
		    cURLList=list()
                    #Seek matching path elements
		    mIndex=[pd[i]==pu[i] for i in range(min(len(pd),len(pu)))].index(False)
		    cURLList=pu[mIndex:]
		    cStringList=pd[mIndex:]
		    cStringList.reverse()
		    cURLList.reverse()
		    cURL=cString=""
		    for elem in cURLList:
			    cURL=cURL+"%s%s"%(os.path.sep,elem)
		    cURL=cURL+os.path.sep
		    if not self.pURL.startswith(os.path.sep):
			    cURL=cURL.lstrip(os.path.sep)
		    self.commonURL=os.path.normpath(cURL).replace(protocolTag,givenProtocol)
		    for elem in cStringList:
			    cString=cString+"%s%s"%(os.path.sep,elem)
		    cString=cString+os.path.sep
		    if not self.pDIR.startswith(os.path.sep):
			    cString=cString.lstrip(os.path.sep)
		    self.commonString=os.path.normpath(cString)

	  def publication_directory(self):
		  return self.pDIR

	  def publication_URL(self):
		  return self.pURL

	  def convert(self,filename=None):
		    #Strip of common path and create full blown URL
		    myURL=filename.replace(self.commonString,self.commonURL)
		    #Add a check to see if given filename is actually URL already!
		    if filename.strip() == "":
			    sys.stderr.write("Improper conversion for :%s\n"%filename)
			    raise Error, "object:filenameToURLMapper given empty string to convert!\n"
			    
		    if myURL == filename:
			    sys.stderr.write("Improper conversion for :%s\n"%filename)
			    sys.stderr.write("web-url        : %s\n"%self.pURL)
			    sys.stderr.write("publication dir: %s\n"%self.pDIR)
			    sys.stderr.write("Common String  : %s\n"%self.commonString)
			    sys.stderr.write("Common URL     : %s\n"%self.commonURL)
			    raise Warning, "object:filenameToURLMapper improperly initialized or given bad args\n"
		    if self.verbose:
			    sys.stdout.write("My URL         : %s\n"%myURL)
			    sys.stdout.write("My file        : %s\n"%filename)
			    sys.stdout.write("web-url        : %s\n"%self.pURL)
			    sys.stdout.write("publication dir: %s\n"%self.pDIR)
			    sys.stdout.write("Common String  : %s\n"%self.commonString)
			    sys.stdout.write("Common URL     : %s\n"%self.commonURL)
		    return myURL
