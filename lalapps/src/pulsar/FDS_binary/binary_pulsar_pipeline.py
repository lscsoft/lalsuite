#!/usr/local/bin/python2.3
"""

binary_pulsar_pipeline.py - binary pulsar pipeline driver script

This script produced the necessary condor submit and dag files to run
the binary pulsar pipeline on LIGO/GEO data

based on stochastic_pipe.in by Adam Mercer which itself was
based on inspiral_pipe.in by Duncan Brown
"""

# import required modules
import sys
import os
import getopt
import re
import math
import shutil
import string
import tempfile
import ConfigParser

# append the lalapps python path
sys.path.append('@PYTHONLIBDIR@')

# import the lalapps pipeline modules
from glue import pipeline
import binary_pulsar

# program usage
def usage():
  msg = """\
Usage: lalapps_binary_pulsar_pipeline [options]

  -h, --help               display this message
  -v, --version            print version information and exit

  -d, --selectdata         run LSCdataFind to create frame cache files
  -s, --search             run the F-statistic search
  -m, --makemesh           run the code to generate the orbital template banks
  -c, --coincidence        run coincidence analysis on F-statistic results
  -i, --injections         run monte-carlo injections
  -u, --upperlimits        run upperlimits analysis
  -w, --wipeclean          overwrite all previous results
  
  -C, --co-pipe            run coincidence pipeline
  -F, --usefulldata        use the full dataset as search dataset
  -P, --priority PRIO      run jobs with condor priority PRIO

  -f, --config-file FILE   use configuration file FILE
  -l, --log-path PATH      directory to write condor log file
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hvdsmcCFiuwP:f:l:"
longop = [
"help",
"version",
"selectdata",
"search",
"makemesh",
"coincidence",
"injections",
"upperlimits",
"wipeclean",
"co-pipe",
"usefulldata",
"priority",
"config-file=",
"log-path="
]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

# default options
config_file = None
do_selectdata = 0
do_search = 0
do_makemesh = 0
do_coincidence = 0
do_injections = 0
do_upperlimits = 0
do_wipeclean = 0
do_co_pipe = 0
do_usefulldata = 0
condor_prio = None
config_file = None
log_path = None

# process options
for o, a in opts:
  if o in ("-v", "--version"):
    print "$Id$"
    sys.exit(0)
  elif o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-d", "--selectdata"):
    do_selectdata = 1
  elif o in ("-s", "--search"):
    do_search = 1
  elif o in ("-m", "--makemesh"):
    do_makemesh = 1
  elif o in ("-c", "--coincidence"):
    do_coincidence = 1
  elif o in ("-i", "--injections"):
    do_injections = 1
  elif o in ("-u", "--upperlimits"):
    do_upperlimits = 1
  elif o in ("-w", "--wipeclean"):
    do_wipeclean = 1
  elif o in ("-C", "--co-pipe"):
    do_co_pipe = 1
  elif o in ("-F", "--usefulldata"):
    do_usefulldata = 1
  elif o in ("-P", "--priority"):
    condor_prio = a
  elif o in ("-f", "--config-file"):
    config_file = a
  elif o in ("-l", "--log-path"):
    log_path = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

if not config_file:
  print >> sys.stderr, "No configuration file specified."
  print >> sys.stderr, "Use --config-file FILE to specify location."
  sys.exit(1)

if not log_path:
  print >> sys.stderr, "No log file path specified."
  print >> sys.stderr, "Use --log-path PATH to specify a location."
  sys.exit(1)

if ((do_usefulldata)&(do_selectdata)):
  print >> sys.stderr, "Both --usefulldata and --selectdata specified."
  print >> sys.stderr, "Use only one or the other."
  sys.exit(1)

if not do_co_pipe:
  print >> sys.stderr, "Sorry, you have not selected the coincidence pipeline."
  print >> sys.stderr, "and the non-coincidence pipeline is not finished yet."
  print >> sys.stderr, "Use --co-pipe to do a coincidence analysis."
  sys.exit(1)

# try and make a directory to store the cache files and job logs
try: os.mkdir('cache')
except: pass
try: os.mkdir('logs')
except: pass

# create the config parser object and read in the ini file
cp = ConfigParser.ConfigParser()
cp.read(config_file)

# create a log file that the Condor jobs will write to
basename = re.sub(r'\.ini',r'',config_file)
tempfile.tempdir = log_path
tempfile.template = basename + '.dag.log.'
logfile = tempfile.mktemp()
fh = open( logfile, "w" )
fh.close()

# create the DAG writing the log to the specified directory
dag = pipeline.CondorDAG(logfile)
dag.set_dag_file(basename + '.dag')

# create the Condor jobs that will be used in the DAG
sensitivity_job = binary_pulsar.sensitivityJob(cp)
makemesh_job = binary_pulsar.makemeshJob(cp)
search_job = binary_pulsar.searchJob(cp)
coincidence_job = binary_pulsar.coincidenceJob(cp)
injections_job = binary_pulsar.injectionsJob(cp)
upperlimit_job = binary_pulsar.upperlimitJob(cp)

# set file names
subsuffix = '.sub'
sensitivity_job.set_sub_file(basename + '.sensitivity' + subsuffix)
makemesh_job.set_sub_file(basename + '.makemesh' + subsuffix)
search_job.set_sub_file(basename + '.search' + subsuffix)
coincidence_job.set_sub_file(basename + '.coincidence' + subsuffix)
injections_job.set_sub_file(basename + '.injections' + subsuffix)
upperlimit_job.set_sub_file(basename + '.upperlimit' + subsuffix)

# set the condor job priority
if condor_prio:
  sensitivity_job.add_condor_cmd('priority',condor_prio)
  makemesh_job.add_condor_cmd('priority',condor_prio)
  search_job.add_condor_cmd('priority',condor_prio)
  coincidence_job.add_condor_cmd('priority',condor_prio)
  injections_job.add_condor_cmd('priority',condor_prio)
  upperlimit_job.add_condor_cmd('priority',condor_prio)

# read in the config file parameters
sensitivity_code = cp.get('condor','sensitivity')
makemesh_code = cp.get('condor','makemesh')
search_code = cp.get('condor','search')
coincidence_code = cp.get('condor','coincidence')
injections_code = cp.get('condor','injections')
upperlimit_code = cp.get('condor','upperlimit')
getdataparams_code = cp.get('condor','getdataparams')
selectdata_code = cp.get('condor','selectdata')
freqmeshfile_code = cp.get('condor','freqmeshfile')
ephdir = cp.get('ephemeris','ephdir')
yr = cp.get('ephemeris','yr')
fulldata_primary = cp.get('data','fulldata_primary')
fulldata_secondary = cp.get('data','fulldata_secondary')
tsft = (int)(cp.get('data','tsft'))
resultsdir = cp.get('results','resultsdir')
det_primary = cp.get('detectors','det_primary')
det_secondary = cp.get('detectors','det_secondary')
tspan = cp.get('sensitivity','tspan')
tstep = cp.get('sensitivity','tstep')
sensitivity_nodes = (int)(cp.get('nodes','sensitivity_nodes'))
search_nodes = (int)(cp.get('nodes','search_nodes'))
coincidence_nodes = (int)(cp.get('nodes','coincidence_nodes'))
injections_nodes = (int)(cp.get('nodes','injections_nodes'))
sourcefile = cp.get('source','sourcefile')
source = cp.get('source','source')
dterms = cp.get('search','dterms')
windowsize = cp.get('search','windowsize')
dopplermax = cp.get('search','dopplermax')
overres = cp.get('search','overres')
binary = cp.get('search','binary')
thresh_primary = cp.get('thresholds','thresh_primary')
thresh_secondary = cp.get('thresholds','thresh_secondary')
mismatch = cp.get('makemesh','mismatch')
meshband = cp.get('makemesh','band')
co_nbins = cp.get('coincidence','nbins')
h0min = (float)(cp.get('injections','h0min'))
h0max = (float)(cp.get('injections','h0max'))
h0step = (float)(cp.get('injections','h0step'))
injection_band = cp.get('injections','injectionband')
ntrials = cp.get('injections','ntrials')
confidence = cp.get('upperlimit','confidence')
minres = cp.get('general','minres')

######################################################################
# define all directory names

generalresultsdir = resultsdir + '/generalresults'
sensitivitydir = resultsdir + '/sensitivityresults'
p_sensitivitydir = sensitivitydir + '/' + det_primary
s_sensitivitydir = sensitivitydir + '/' + det_secondary
optdatadir = resultsdir + '/optimaldata'
p_optdatadir = optdatadir + '/' + det_primary
s_optdatadir = optdatadir + '/' + det_secondary
templatebankdir = resultsdir + '/binarytemplatebanks'
p_templatebankdir = templatebankdir + '/' + det_primary
s_templatebankdir = templatebankdir + '/' + det_secondary
searchresultsdir = resultsdir + '/searchresults'
p_searchresultsdir = searchresultsdir + '/' + det_primary
s_searchresultsdir = searchresultsdir + '/' + det_secondary
coincidencedir = resultsdir + '/coincidenceresults'
injectionsdir = resultsdir + '/injections'

######################################################################
# create the results directory structure
# be careful if previous results already exist

# remove all previous results if requested
if do_wipeclean:
  if os.path.exists(resultsdir):
    print "WARNING : wiping results directory [%s] clean !!!" %(resultsdir)
    wipeclean_flag = input("Are you sure you want to delete EVERYTHING ? (y/n) ")
    if (wipeclean_flag=="y"):
      try: os.system('rm -rf ' + resultsdir)
      except:
        print "ERROR : failed to remove results directory [%s] !!!" %(resultsdir)
        sys.exit(1)
    elif (wipeclean_flag=="n"):
      print "EXITING : OK, will not delete directory [%s] !!!" %(resultsdir)
      sys.exit(1)
    else:
      print "ERROR : need to answer y or n !!!"
      sys.exit(1)
  else:
    print "WARNING : you have selected --wipeclean but directory [%s] does not exist !!!" %(resultsdir)
# else we make the main results directory if it doesn't exist
else:
  if os.path.exists(resultsdir):
    print "WARNING : results directory [%s] already exists !!!" %(resultsdir)
  else:
    try: os.mkdir(resultsdir)
    except:
      print "ERROR : failed to create results directory [%s] !!!" %(resultsdir)
      sys.exit(1)
      
# if we are not removing previous results then be careful !
# the policy is that if we have requested to do a part of the analysis
# then the directory associated with that part should not exist prior
# to creating this dag

# going to do this in an directory by directory basis
# ---------------------------------------------- sensitivity -----------------------------------------------

if ((os.path.exists(sensitivitydir))&(do_selectdata)):

  if os.path.exists(p_sensitivitydir):
    sensitivity_flag = (str)(raw_input("Are you sure you want to delete the primary sensitivity directory ? (y/n) "))
    if (sensitivity_flag=="y"):
      try: shutil.rmtree(p_sensitivitydir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(p_sensitivitydir)
        sys.exit(1)
      try: os.mkdir(p_sensitivitydir)
      except:
        print "ERROR : failed to create primary sensitivity directory [%s] !!!" %(p_sensitivitydir)
        sys.exit(1) 
    else:
      print "EXITING : OK, not deleting directory [%s]" %(p_sensitivitydir)
      sys.exit(1)
      
  if ((os.path.exists(s_sensitivitydir))&(do_co_pipe)):
    sensitivity_flag = (str)(raw_input("Are you sure you want to delete the secondary sensitivity directory ? (y/n) "))
    if (sensitivity_flag=="y"):
      try: shutil.rmtree(s_sensitivitydir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(s_sensitivitydir)
        sys.exit(1)
      try: os.mkdir(s_sensitivitydir)
      except:
        print "ERROR : failed to create secondary sensitivity directory [%s] !!!" %(s_sensitivitydir)
        sys.exit(1)  
    else:
      print "EXITING : OK, not deleting directory [%s]" %(s_sensitivitydir)
      sys.exit(1)

elif ((not os.path.exists(sensitivitydir))&(do_selectdata)):

  try: os.mkdir(sensitivitydir)
  except:
    print "ERROR : failed to create sensitivity directory [%s] !!!" %(sensitivitydir)
    sys.exit(1)
  try: os.mkdir(p_sensitivitydir)
  except:
    print "ERROR : failed to create primary sensitivity directory [%s] !!!" %(p_sensitivitydir)
    sys.exit(1)
  if (do_co_pipe):  
    try: os.mkdir(s_sensitivitydir)
    except:
      print "ERROR : failed to create secondary sensitivity directory [%s] !!!" %(s_sensitivitydir)
      sys.exit(1)   
    

# ---------------------------------------------- optimal -----------------------------------------------

if ((os.path.exists(optdatadir))&(do_makemesh)):

  if os.path.exists(p_optdatadir):
    optimal_flag = (str)(raw_input("Are you sure you want to delete the primary optimal data directory ? (y/n) "))
    if (optimal_flag=="y"):
      try: shutil.rmtree(p_optdatadir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(p_optdatadir)
        sys.exit(1)
      try: os.mkdir(p_optdatadir)
      except:
        print "ERROR : failed to create primary optimal data directory [%s] !!!" %(p_optdatadir)
        sys.exit(1) 
    else:
      print "EXITING : OK, not deleting directory [%s]" %(p_optdatadir)
      sys.exit(1)
      
  if ((os.path.exists(s_optdatadir))&(do_co_pipe)):
    optimal_flag = (str)(raw_input("Are you sure you want to delete the secondary optimal data directory ? (y/n) "))
    if (optimal_flag=="y"):
      try: shutil.rmtree(s_optdatadir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(s_optdatadir)
        sys.exit(1)
      try: os.mkdir(s_optdatadir)
      except:
        print "ERROR : failed to create secondary optimal data directory [%s] !!!" %(s_optdatadir)
        sys.exit(1)  
    else:
      print "EXITING : OK, not deleting directory [%s]" %(s_optdatadir)
      sys.exit(1)

elif ((not os.path.exists(optdatadir))&(do_makemesh)):

  try: os.mkdir(optdatadir)
  except:
    print "ERROR : failed to create optimal data directory [%s] !!!" %(optdatadir)
    sys.exit(1)
  try: os.mkdir(p_optdatadir)
  except:
    print "ERROR : failed to create primary optimal data directory [%s] !!!" %(p_optdatadir)
    sys.exit(1)
  if (do_co_pipe):  
    try: os.mkdir(s_optdatadir)
    except:
      print "ERROR : failed to create secondary optimal data directory [%s] !!!" %(s_optdatadir)
      sys.exit(1)   

# ---------------------------------------------- makemesh -----------------------------------------------

if ((os.path.exists(templatebankdir))&(do_makemesh)):

  if os.path.exists(p_templatebankdir):
    makemesh_flag = (str)(raw_input("Are you sure you want to delete the primary template bank directory ? (y/n) "))
    if (makemesh_flag=="y"):
      try: shutil.rmtree(p_templatebankdir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(p_templatebankdir)
        sys.exit(1)
      try: os.mkdir(p_templatebankdir)
      except:
        print "ERROR : failed to create primary template bank directory [%s] !!!" %(p_templatebankdir)
        sys.exit(1) 
    else:
      print "EXITING : OK, not deleting directory [%s]" %(p_templatebankdir)
      sys.exit(1)
      
  if ((os.path.exists(s_templatebankdir))&(do_co_pipe)):
    makemesh_flag = (str)(raw_input("Are you sure you want to delete the secondary template bank directory ? (y/n) "))
    if (makemesh_flag=="y"):
      try: shutil.rmtree(s_templatebankdir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(s_templatebankdir)
        sys.exit(1)
      try: os.mkdir(s_templatebankdir)
      except:
        print "ERROR : failed to create secondary template bank directory [%s] !!!" %(s_templatebankdir)
        sys.exit(1)  
    else:
      print "EXITING : OK, not deleting directory [%s]" %(s_templatebankdir)
      sys.exit(1)

elif ((not os.path.exists(templatebankdir))&(do_makemesh)):

  try: os.mkdir(templatebankdir)
  except:
    print "ERROR : failed to create template bank directory [%s] !!!" %(templatebankdir)
    sys.exit(1)
  try: os.mkdir(p_templatebankdir)
  except:
    print "ERROR : failed to create primary template bank directory [%s] !!!" %(p_templatebankdir)
    sys.exit(1)
  if (do_co_pipe):  
    try: os.mkdir(s_templatebankdir)
    except:
      print "ERROR : failed to create secondary template bank directory [%s] !!!" %(s_templatebankdir)
      sys.exit(1)   

# ---------------------------------------------- search -----------------------------------------------

if ((os.path.exists(searchresultsdir))&(do_search)):

  if os.path.exists(p_searchresultsdir):
    search_flag = (str)(raw_input("Are you sure you want to delete the primary search results directory ? (y/n) "))
    if (search_flag=="y"):
      try: shutil.rmtree(p_searchresultsdir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(p_searchresultsdir)
        sys.exit(1)
      try: os.mkdir(p_searchresultsdir)
      except:
        print "ERROR : failed to create primary search results directory [%s] !!!" %(p_searchresultsdir)
        sys.exit(1) 
    else:
      print "EXITING : OK, not deleting directory [%s]" %(p_searchresultsdir)
      sys.exit(1)
      
  if ((os.path.exists(s_searchresultsdir))&(do_co_pipe)):
    search_flag = (str)(raw_input("Are you sure you want to delete the secondary search results directory ? (y/n) "))
    if (search_flag=="y"):
      try: shutil.rmtree(s_searchresultsdir)
      except:
        print "ERROR : failed to delete directory [%s] !!!" %(s_searchresultsdir)
        sys.exit(1)
      try: os.mkdir(s_searchresultsdir)
      except:
        print "ERROR : failed to create secondary search results directory [%s] !!!" %(s_searchresultsdir)
        sys.exit(1)  
    else:
      print "EXITING : OK, not deleting directory [%s]" %(s_searchresultsdir)
      sys.exit(1)

elif ((not os.path.exists(searchresultsdir))&(do_search)):

  try: os.mkdir(searchresultsdir)
  except:
    print "ERROR : failed to create search results directory [%s] !!!" %(searchresultsdir)
    sys.exit(1)
  try: os.mkdir(p_searchresultsdir)
  except:
    print "ERROR : failed to create primary search results directory [%s] !!!" %(p_searchresultsdir)
    sys.exit(1)
  if (do_co_pipe):  
    try: os.mkdir(s_searchresultsdir)
    except:
      print "ERROR : failed to create secondary search results directory [%s] !!!" %(s_searchresultsdir)
      sys.exit(1)   

# ---------------------------------------------- coincidence -----------------------------------------------

if ((os.path.exists(coincidencedir))&(do_coincidence)):

  coincidence_flag = (str)(raw_input("Are you sure you want to delete the coincidence results directory ? (y/n) "))
  if (coincidence_flag=="y"):
    try: shutil.rmtree(coincidencedir)
    except:
      print "ERROR : failed to delete directory [%s] !!!" %(coincidencedir)
      sys.exit(1)
    try: os.mkdir(coincidencedir)
    except:
      print "ERROR : failed to create coincidence results directory [%s] !!!" %(coincidencedir)
      sys.exit(1) 
  else:
    print "EXITING : OK, not deleting directory [%s]" %(coincidencedir)
    sys.exit(1)
      
elif ((not os.path.exists(coincidencedir))&(do_coincidence)):

  try: os.mkdir(coincidencedir)
  except:
    print "ERROR : failed to create coincidence results directory [%s] !!!" %(coincidencedir)
    sys.exit(1)
 
# ---------------------------------------------- injections -----------------------------------------------

if ((os.path.exists(injectionsdir))&(do_injections)):

  injections_flag = (str)(raw_input("Are you sure you want to delete the injections results directory ? (y/n) "))
  if (injections_flag=="y"):
    try: shutil.rmtree(injectionsdir)
    except:
      print "ERROR : failed to delete directory [%s] !!!" %(injectionsdir)
      sys.exit(1)
    try: os.mkdir(injectionsdir)
    except:
      print "ERROR : failed to create injections results directory [%s] !!!" %(injectionsdir)
      sys.exit(1) 
  else:
    print "EXITING : OK, not deleting directory [%s]" %(injectionsdir)
    sys.exit(1)
      
elif ((not os.path.exists(injectionsdir))&(do_injections)):

  try: os.mkdir(injectionsdir)
  except:
    print "ERROR : failed to create injections results directory [%s] !!!" %(injectionsdir)
    sys.exit(1)
 
# ---------------------------------------------- general -----------------------------------------------

if os.path.exists(generalresultsdir):

  general_flag = (str)(raw_input("Do you want to delete the general results directory ? (y/n) "))
  if (general_flag=="y"):
    try: shutil.rmtree(generalresultsdir)
    except:
      print "ERROR : failed to delete directory [%s] !!!" %(generalresultsdir)
      sys.exit(1)
    try: os.mkdir(generalresultsdir)
    except:
      print "ERROR : failed to create general results directory [%s] !!!" %(generalresultsdir)
      sys.exit(1) 
  else:
    print "STATUS : OK, not deleting directory [%s]" %(generalresultsdir)
      
elif not os.path.exists(generalresultsdir):

  try: os.mkdir(generalresultsdir)
  except:
    print "ERROR : failed to create general results directory [%s] !!!" %(generalresultsdir)
    sys.exit(1)

 
######################################################################
# create the jobs here so we can sensibly do dependencies
# therefore we have to split up the times and frequencies here 

# get the primary data set parameters

p_fulldatasetparams = generalresultsdir + '/fulldatasetparameters_primary.data'
p_fullstampsfile = generalresultsdir + '/fulldatasettimestamps_primary.data'
p_optdatasetparams = generalresultsdir + '/optdatasetparameters_primary.data'
p_optstampsfile = generalresultsdir + '/optdatasettimestamps_primary.data'
try: os.system("%s --datadir %s --tsft %d --stampsfile %s --outfile %s" %(getdataparams_code, fulldata_primary, tsft, p_fullstampsfile, p_fulldatasetparams))
except:
  print "ERROR : failed to get primary dataset params, exiting."
  sys.exit(1)

if do_co_pipe:

  s_fulldatasetparams = generalresultsdir + '/fulldatasetparameters_secondary.data'
  s_fullstampsfile = generalresultsdir + '/fulldatasettimestamps_secondary.data'
  s_optdatasetparams = generalresultsdir + '/optdatasetparameters_secondary.data'
  s_optstampsfile = generalresultsdir + '/optdatasettimestamps_secondary.data'
  try: os.system("%s --datadir %s --tsft %d --stampsfile %s --outfile %s" %(getdataparams_code, fulldata_secondary, tsft, s_fullstampsfile, s_fulldatasetparams))
  except:
    print "ERROR : failed to get secondary dataset params, exiting."
    sys.exit(1)

######################################################################
# now divide the full dataset into smaller chunks for splitting over the nodes

if do_selectdata:

  # create the config parser object and read in the start and end times of the full dataset
  fdp = ConfigParser.ConfigParser()
  fdp.read(p_fulldatasetparams)
  p_full_tstart = (int)(fdp.get('dataparams','start'))
  p_full_tend = (int)(fdp.get('dataparams','end'))
  p_full_tobs = (int)(fdp.get('dataparams','tobs'))
  p_full_tspan = (int)(fdp.get('dataparams','tspan'))
  p_sens_nodes = sensitivity_nodes
  
  # if doing coincidence pipeline then read in params and split up nodes
  if do_co_pipe:
    sdp = ConfigParser.ConfigParser()
    sdp.read(s_fulldatasetparams)
    s_full_tstart = (int)(sdp.get('dataparams','start'))
    s_full_tend = (int)(sdp.get('dataparams','end'))
    s_full_tobs = (int)(sdp.get('dataparams','tobs'))
    s_full_tspan = (int)(sdp.get('dataparams','tspan'))
    
    # now do division of nodes based on estimated time to complete each dataset
    p_sens_nodes = math.floor(float(sensitivity_nodes)*((float(p_full_tobs))/float(p_full_tobs+s_full_tobs))+0.5)
    s_sens_nodes = math.floor(float(sensitivity_nodes)*((float(s_full_tobs))/float(s_full_tobs+p_full_tobs))+0.5)
    #print p_sens_nodes
    #print s_sens_nodes
  
  # generate the sensitivity primary start times
  p_sens_start = []
  p_sens_end = []
  p_timegap = (int)(math.ceil(float(p_full_tspan)/float(p_sens_nodes)))
  i=0
  while i < p_sens_nodes:
    p_sens_start.append((int)(p_full_tstart+(i*p_timegap)))
    p_sens_end.append((int)(p_full_tstart+((i+1)*p_timegap)))
    i = i+1

  
  if do_co_pipe:
    # generate the sensitivity secondary start times
    s_sens_start = []
    s_sens_end = []
    s_timegap = (int)(math.ceil(float(s_full_tspan)/float(s_sens_nodes)))
    i=0
    while i < s_sens_nodes:
      s_sens_start.append((int)(s_full_tstart+(i*s_timegap)))
      s_sens_end.append((int)(s_full_tstart+((i+1)*s_timegap)))
      i = i+1
               
######################################################################
# create the sensitivity jobs

if do_selectdata:

  sensitivity_primary = []
  i=0
  while i < p_sens_nodes:
  
    # fill in the primary sensitivity params
    sensitivity_primary.append(binary_pulsar.sensitivityNode(sensitivity_job))
    sensitivity_primary[i].set_datadir(fulldata_primary)
    sensitivity_primary[i].set_detector(det_primary)
    sensitivity_primary[i].set_tstart(p_sens_start[i])
    sensitivity_primary[i].set_tend(p_sens_end[i])
    sensitivity_primary[i].set_outdir(p_sensitivitydir)

    # add the primary sensitivity nodes
    dag.add_node(sensitivity_primary[i])
    
    i=i+1
    
  if do_co_pipe:

    sensitivity_secondary = []
    i=0
    while i < s_sens_nodes:
    
      # fill in the secondary sensitivity params
      sensitivity_secondary.append(binary_pulsar.sensitivityNode(sensitivity_job))
      sensitivity_secondary[i].set_datadir(fulldata_secondary)
      sensitivity_secondary[i].set_detector(det_secondary)
      sensitivity_secondary[i].set_tstart(s_sens_start[i])
      sensitivity_secondary[i].set_tend(s_sens_end[i])
      sensitivity_secondary[i].set_outdir(s_sensitivitydir)

      # add the primary sensitivity nodes
      dag.add_node(sensitivity_secondary[i])
      
      i=i+1


######################################################################
# create the makemesh jobs

if do_makemesh:

  # define some output file names
  p_optdataparams = generalresultsdir + '/optimal_data_params_primary.data'
  s_optdataparams = generalresultsdir + '/optimal_data_params_secondary.data'

  # fill in the primary makemesh params
  makemesh_primary = binary_pulsar.makemeshNode(makemesh_job)
  makemesh_primary.set_datadir(p_optdatadir)
  makemesh_primary.set_meshdir(p_templatebankdir)
  makemesh_primary.set_detector(det_primary)
  
  # add the neccessary prescript arguments
  makemesh_primary.set_pre_script(selectdata_code)
  makemesh_primary.add_pre_script_arg("--sensitivitydir %s" %(p_sensitivitydir))
  makemesh_primary.add_pre_script_arg("--fulldata %s" %(fulldata_primary))
  makemesh_primary.add_pre_script_arg("--optdatadir %s" %(p_optdatadir))
  makemesh_primary.add_pre_script_arg("--optdataparamsfile %s" %(p_optdataparams))
  makemesh_primary.add_pre_script_arg("--sourcefile %s" %(sourcefile))
  makemesh_primary.add_pre_script_arg("--source %s" %(source))

  # add the neccessary postscript arguments
  makemesh_primary.set_post_script(getdataparams_code)
  makemesh_primary.add_post_script_arg("--datadir %s" %(p_optdatadir))
  makemesh_primary.add_post_script_arg("--tsft %s" %(tsft))
  makemesh_primary.add_post_script_arg("--stampsfile %s" %(p_optstampsfile))
  makemesh_primary.add_post_script_arg("--outfile %s" %(p_optdatasetparams))
  
  # if we are not selecting the data then just copy the full data set
  if do_usefulldata:
    makemesh_primary.add_pre_script_arg("--justcopy")
  
  # add the primary makemesh dependencies
  if do_selectdata:
    i=0
    while i < p_sens_nodes:
      makemesh_primary.add_parent(sensitivity_primary[i])
      i=i+1
  
  # add the job to the dag
  dag.add_node(makemesh_primary)

  if do_co_pipe:

    # fill in the primary makemesh params
    makemesh_secondary = binary_pulsar.makemeshNode(makemesh_job)
    makemesh_secondary.set_datadir(s_optdatadir)
    makemesh_secondary.set_meshdir(s_templatebankdir)
    makemesh_secondary.set_detector(det_secondary)
    
    # add the neccessary prescript arguments
    makemesh_secondary.set_pre_script(selectdata_code)
    makemesh_secondary.add_pre_script_arg("--sensitivitydir %s" %(s_sensitivitydir))
    makemesh_secondary.add_pre_script_arg("--fulldata %s" %(fulldata_secondary))
    makemesh_secondary.add_pre_script_arg("--optdatadir %s" %(s_optdatadir))
    makemesh_secondary.add_pre_script_arg("--optdataparamsfile %s" %(s_optdataparams))
    makemesh_secondary.add_pre_script_arg("--sourcefile %s" %(sourcefile))
    makemesh_secondary.add_pre_script_arg("--source %s" %(source))

    # add the neccessary postscript arguments
    makemesh_secondary.set_post_script(getdataparams_code)
    makemesh_secondary.add_post_script_arg("--datadir %s" %(s_optdatadir))
    makemesh_secondary.add_post_script_arg("--tsft %s" %(tsft))
    makemesh_secondary.add_post_script_arg("--stampsfile %s" %(s_optstampsfile))
    makemesh_secondary.add_post_script_arg("--outfile %s" %(s_optdatasetparams))

      # if we are not selecting the data then just copy the full data set
    if do_usefulldata:
      makemesh_secondary.add_pre_script_arg("--justcopy")
    
    # add the primary makemesh dependencies
    if do_selectdata:
      i=0
      while i < s_sens_nodes:
        makemesh_secondary.add_parent(sensitivity_secondary[i])
        i=i+1
        
    # add the job to the dag
    dag.add_node(makemesh_secondary)
        
######################################################################
# setup parameters for the search

if ((do_search)|(do_coincidence)):

  # define freq-mesh search file name
  freqmesh_searchfile = generalresultsdir + '/search_freqmeshfile.data'
 
  # run script here to produce a search freq/mesh file
  try: os.system("%s --sourcefile %s --source %s --meshband %s --pdet %s --sdet %s --meshdir %s --nodes %s --minres %s --outfile %s"  %(freqmeshfile_code, sourcefile, source, meshband, det_primary, det_secondary, templatebankdir, search_nodes, minres, freqmesh_searchfile))
  except:
    print "ERROR : failed to generate a freq-mesh parameter file for the search !!!"
    sys.exit(1)

  # split up the frequency bands by reading in from the searchfreqmesh file
  f_min_search = []
  f_max_search = []
  f_band_search = []
  templatebank_primary = []
  templatebank_secondary = []

  # open the freqmesh file
  try: fp = open(freqmesh_searchfile,"r")
  except:
    print "ERROR : could not open search freq-mesh file [%s] !!!" %(freqmesh_searchfile)
    sys.exit(1)

  # loop over the contents of the file
  i=0
  for line in fp:
    temp = line.rstrip().split()
    f_min_search.append(float(temp[0]))
    f_max_search.append(float(temp[1]))
    f_band_search.append(float(temp[2]))
    templatebank_primary.append(temp[3])
    templatebank_secondary.append(temp[4])
    i=i+1

  # define the new number of search nodes
  co_search_nodes = i

  # close the file
  fp.close()

######################################################################
# create the search jobs

if do_search:

  search_primary = []
  i=0
  while i < co_search_nodes:

    # define output label
    p_searchlabel = "_%s_%.3f-%.3f.data" %(det_primary,f_min_search[i],f_max_search[i])
    
    search_primary.append(binary_pulsar.searchNode(search_job))
    search_primary[i].set_f_min(f_min_search[i])
    search_primary[i].set_f_band(f_band_search[i])
    search_primary[i].set_binarytemplatebank(templatebank_primary[i])
    search_primary[i].set_datadir(p_optdatadir)
    search_primary[i].set_detector(det_primary)
    search_primary[i].set_workdir(p_searchresultsdir)
    search_primary[i].set_ephdir(ephdir)
    search_primary[i].set_yr(yr)
    search_primary[i].set_Fthresh(thresh_primary)
    search_primary[i].set_label(p_searchlabel)

    # add primary search dependencies
    if do_makemesh:
      search_primary[i].add_parent(makemesh_primary)

    # add primary search dag
    dag.add_node(search_primary[i])
    i=i+1

  if do_co_pipe:

    search_secondary = []
    i=0
    while i < co_search_nodes:

      # define output label
      s_searchlabel = "_%s_%.3f-%.3f.data" %(det_secondary,f_min_search[i],f_max_search[i])
      
      search_secondary.append(binary_pulsar.searchNode(search_job))
      search_secondary[i].set_f_min(f_min_search[i])
      search_secondary[i].set_f_band(f_band_search[i])
      search_secondary[i].set_binarytemplatebank(templatebank_secondary[i])
      search_secondary[i].set_datadir(s_optdatadir)
      search_secondary[i].set_detector(det_secondary)
      search_secondary[i].set_workdir(s_searchresultsdir)
      search_secondary[i].set_ephdir(ephdir)
      search_secondary[i].set_yr(yr)
      search_secondary[i].set_Fthresh(thresh_secondary)
      search_secondary[i].set_label(s_searchlabel)

      # add secondary search dependencies
      if do_makemesh:
        search_secondary[i].add_parent(makemesh_secondary)

      # we are also adding the primary search jobs as dependencies because
      # depending on the datasets selected these can have dramatically different
      # running times.  So we will run them in sequence not in coincidence
      #j=0
      #while j < co_search_nodes:
        #search_secondary[i].add_parent(search_primary[j])
      
      # add secondary search dag
      dag.add_node(search_secondary[i])
      i=i+1

######################################################################
# setup the coincidence jobs frequencies

if ((do_co_pipe)&(do_coincidence)):

  # define freq-mesh search file name
  freqmesh_cofile = generalresultsdir + '/coincidence_freqmeshfile.data'
   
  # run script here to produce a search freq/mesh file (Note that this is only to split up frequencies)
  try: os.system("%s --sourcefile %s --source %s --meshband %s --pdet %s --sdet %s --meshdir %s --nodes %s --minres %s --outfile %s"  %(freqmeshfile_code, sourcefile, source, meshband, det_primary, det_secondary, templatebankdir, coincidence_nodes, minres, freqmesh_cofile))
  except:
    print "ERROR : failed to generate a freq-mesh parameter file for the coincidence !!!"
    sys.exit(1)

  # split up the frequency bands by reading in from the searchfreqmesh file
  f_min_co = []
  f_max_co = []
  f_band_co = []
  templatebank_primary_co = []
  templatebank_secondary_co = []

  # open the freqmesh file
  try: fp = open(freqmesh_cofile,"r")
  except:
    print "ERROR : could not open coincidence freq-mesh file [%s] !!!" %(freqmesh_cofile)
    sys.exit(1)

  # loop over the contents of the file
  i=0
  for line in fp:
    temp = line.rstrip().split()
    f_min_co.append(float(temp[0]))
    f_max_co.append(float(temp[1]))
    f_band_co.append(float(temp[2]))
    templatebank_primary_co.append(temp[3])
    templatebank_secondary_co.append(temp[4])
    i=i+1

  # define new coincidence nodes
  coincidence_nodes = i

  # close the file
  fp.close()

######################################################################
# create the coincidence jobs

if ((do_co_pipe)&(do_coincidence)):

  coincidence = []
  i=0
  while i < coincidence_nodes:
    coincidence.append(binary_pulsar.coincidenceNode(coincidence_job))
    coincidence[i].set_presultsdir(p_searchresultsdir)
    coincidence[i].set_sresultsdir(s_searchresultsdir)
    coincidence[i].set_fmin(f_min_co[i])
    coincidence[i].set_fmax(f_max_co[i])
    coincidence[i].set_freqmeshfile(freqmesh_searchfile)
    coincidence[i].set_coresultsdir(coincidencedir)
    coincidence[i].set_nbins(co_nbins)

    if do_selectdata:
      coincidence[i].set_sdataparamsfile(s_optdatasetparams)
    else:
      coincidence[i].set_sdataparamsfile(s_fulldatasetparams)

    # add secondary search dependencies (need to match the frequencies !!)
    if do_search:
      j=0
      while j < co_search_nodes:
        coincidence[i].add_parent(search_primary[j])
        coincidence[i].add_parent(search_secondary[j])
        j=j+1
    
    # add secondary search dag
    dag.add_node(coincidence[i])
    i=i+1

######################################################################
# create the injections jobs

if ((do_injections)|(do_upperlimits)):

  # need to loop over the frequencies, h0 and trials
  # do this so that we split things up over the requested
  # number of nodes

  # so minimum number of nodes is number of frequencies * number of h0 values

  # generate h0 values and save to file
  h0 = []
  #print "h0max is %e h0min is %e h0step is %e" %(h0max,h0min,h0step)
  nh0 = math.ceil(1.0 + (float(h0max-h0min)/float(h0step)))
  #print "nh0 is %f" %(nh0) 
  i=0
  
  # open the file
  h0_injfile = generalresultsdir + '/h0_injectionvalues.data'
  try: fp = open(h0_injfile, 'w')
  except:
    print "ERROR : strange, can't open file for outputting h0 values !!!"
    sys.exit(1)

  while i < nh0:
    h0.append("%.3e" %(h0min+(i*h0step)))
    fp.write("%.3e\n" %(float(h0[i])))
    i=i+1

  fp.close  
   
  # define output freqmeshfile name
  freqmesh_injfile = generalresultsdir + '/injections_freqmeshfile.data'

  # generate frequency values
  try: os.system("%s --sourcefile %s --source %s --splitband %s --meshband %s --pdet %s --sdet %s --meshdir %s --minres %s --outfile %s"  %(freqmeshfile_code, sourcefile, source, injection_band, meshband, det_primary, det_secondary, templatebankdir, minres, freqmesh_injfile))
  except:
    print "ERROR : failed to generate freq/mesh file for injections" 
    sys.exit(1)

  # read from freqmeshfile the frequencies and meshes to use
  f_min_injections = []
  f_max_injections = []
  templatebank_primary = []
  templatebank_secondary = []
  # open the file
  try: fp = open(freqmesh_injfile, 'r')
  except:
    print "ERROR : strange, can't open freqmesh file %s for injections !!!" %(freqmesh_injfile)
    sys.exit(1)
      
  # loop over each line in this file  
  nfreq=0
  for line in fp:
    temp = line.rstrip().split()
    f_min_injections.append(temp[0])
    f_max_injections.append(temp[1])
    templatebank_primary.append(temp[3])
    templatebank_secondary.append(temp[4])
    nfreq=nfreq+1 

  fp.close()

  # generate the ntrials value
  ntemp = nh0*nfreq
  n_ntrials = math.floor(injections_nodes/ntemp)
  if (n_ntrials==0):
    n_ntrials = 1
  trials = (int)(math.ceil(float(ntrials)/float(n_ntrials)))
  injection_nodes_new = nh0*nfreq*n_ntrials
  #print "nh0 is %d nfreq is %d n_ntrails is %d" %(nh0,nfreq,n_ntrials)
  #print "ntrials is %d trials is %d" %((int)(ntrials),(int)(trials))
  #print "injection nodes is %d" %(injection_nodes_new)

if do_injections:

  injections = []

  # start loop over frequency
  i=0
  count=0
  while i < nfreq:

    j=0
    while j < nh0:

      k=0
      while k < n_ntrials:

        # add the injection job params
        injections.append(binary_pulsar.injectionsNode(injections_job))
        injections[count].set_configfile(config_file)
        injections[count].set_f_min(f_min_injections[i])
        injections[count].set_f_max(f_max_injections[i])
        injections[count].set_pdatadir(p_optdatadir)
        injections[count].set_sdatadir(s_optdatadir)
        injections[count].set_ptemplatefile(templatebank_primary[i])
        injections[count].set_stemplatefile(templatebank_secondary[i])
        injections[count].set_gw_amplitude(h0[j])
        injections[count].set_ntrials(trials)
        injections[count].set_id(k)
        injections[count].set_outdir(injectionsdir)

        # add injections dependencies
        if do_search:
          injections[count].add_parent(makemesh_primary)
          injections[count].add_parent(makemesh_secondary)
                  
        # add injections job to dag
        dag.add_node(injections[count])

        # increment the count and the loop over ntrials
        count=count+1
        k=k+1

      # increment the loop over h0  
      j=j+1

    # increment the loop over freq bands  
    i=i+1  

  # record how many injection jobs were made
  injection_nodes = (int)(count)

######################################################################
# create the upperlimit job

if do_upperlimits:

  # define some output files
  maxoutfile = generalresultsdir + '/loudest_events.data'
  upperlimit_outfile = generalresultsdir + '/UL_results'

  # set the upperlimit parameters
  upperlimit = binary_pulsar.upperlimitNode(upperlimit_job)
  upperlimit.set_injectionsdir(injectionsdir)
  upperlimit.set_coresultsdir(coincidencedir)
  upperlimit.set_injfreqmeshfile(freqmesh_injfile)
  upperlimit.set_injh0file(h0_injfile)
  upperlimit.set_maxoutfile(maxoutfile)
  upperlimit.set_outfile(upperlimit_outfile)

  # add the upperlimit dependancies
  if do_injections:
    i=0
    while i < injection_nodes:
      upperlimit.add_parent(injections[i]);
      i=i+1
  if do_coincidence:
    i=0
    while i < coincidence_nodes:
      upperlimit.add_parent(coincidence[i]);
      i=i+1

  # add the upperlimits job to the dag
  dag.add_node(upperlimit)

######################################################################
# output the dag

# write out the DAG
dag.write_sub_files()
dag.write_dag()


sys.exit(0)

