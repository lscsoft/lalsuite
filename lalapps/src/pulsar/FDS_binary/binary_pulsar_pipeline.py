#!/usr/bin/env @PYTHONPROG@
"""
$Id$

binary_pulsar_pipeline.py - binary pulsar pipeline driver script

This script produced the necessary condor submit and dag files to run
the binary pulsar pipeline on LIGO/GEO data

based on stochastic_pipe.in by Adam Mercer which itself was
based on inspiral_pipe.in by Duncan Brown
"""

__author__ = 'Adam Mercer <ram@star.sr.bham.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import ConfigParser

# append the lalapps python path
sys.path.append('@PYTHONLIBDIR@')

# import the lalapps pipeline modules
from glue import pipeline
import stochastic

# program usage
def usage():
  msg = """\
Usage: lalapps_stochastic_pipe [options]

  -h, --help               display this message
  -v, --version            print version information and exit

  -d, --selectdata         run LSCdataFind to create frame cache files
  -s, --search             run the F-statistic search
  -c, --coincidence        run coincidence pipeline
  -i, --injections         run monte-carlo injections
  -u, --upperlimits        run upperlimits analysis

  -P, --priority PRIO      run jobs with condor priority PRIO

  -f, --config-file FILE   use configuration file FILE
  -l, --log-path PATH      directory to write condor log file
  """
  print >> sys.stderr, msg

# parse the command line options to figure out what we should do
shortop = "hvdstP:f:l:"
longop = [
"help",
"version",
"selectdata",
"search",
"coincidence",
"injections",
"upperlimits",
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
do_selectdata = None
do_search = None
do_coincidence = None
do_injections = None
do_upperlimits = None
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
  elif o in ("-c", "--coincidence"):
    do_coincidence = 1
  elif o in ("-i", "--injections"):
    do_injections = 1
  elif o in ("-u", "--upperlimits"):
    do_upperlimits = 1 
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
selectdata_job = binary_pulsar.selectdataJob(cp)
makemesh_job = binary_pulsar.makemeshJob(cp)
search_job = binary_pulsar.searchJob(cp)
coincidence = binary_pulsar.coincidencejob(cp)
injections = binary_pulsar.injectionsJob(cp)
upperlimit = binary_pulsar.upperlimitJob(cp)

# set file names
subsuffix = '.sub'
selectdata_job.set_sub_file(basename + '.selectdata'+ subsuffix)
makemesh_job.set_sub_file(basename + '.makemesh' + subsuffix)
search_job.set_sub_file(basename + '.search' + subsuffix)
coincidence_job.set_sub_file(basename + '.coincidence' + subsuffix)
injections_job.set_sub_file(basename + '.injections' + subsuffix)
upperlimit_job.set_sub_file(basename + '.upperlimit' + subsuffix)

# set the condor job priority
if condor_prio:
  selectdata_job.add_condor_cmd('priority',condor_prio)
  makemesh_job.add_condor_cmd('priority',condor_prio)
  search_job.add_condor_cmd('priority',condor_prio)
  coincidence_job.add_condor_cmd('priority',condor_prio)
  injections_job.add_condor_cmd('priority',condor_prio)
  upperlimit_job.add_condor_cmd('priority',condor_prio)

# read in the config file parameters
selectdata_code = cp.get('condor','selectdata')
makemesh_code = cp.get('condor','makemesh')
search_code = cp.get('condor','search')
coincidence_code = cp.get('condor','coincidence')
injections_code = cp.get('condor','injections')
upperlimit_code = cp.get('condor','uppelimit')
ephdir = cp.get('ephemeris','ephdir')
yr = cp.get('ephmeris','yr')
fulldata_primary = cp.get('data','fulldata_primary')
fulldata_secondary = cp.get('data','fulldata_secondary')
resultsdir = cp.get('results','resultsdir')
det_primary = cp.get('detectors','det_primary')
det_secondary = cp.get('detectors','det_secondary')
tspan = cp.get('selectdata','tspan')
tstep = cp.get('selectdata','tstep')
selectdata_nodes = (int)(cp.get('nodes','selectdata_nodes'))
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
meshband = cp.get('makemesh','meshband')
h0min = (float)(cp.get('injections','h0min'))
h0max = (float)(cp.get('injections','h0max'))
h0step = (floa)(cp.get('injections','h0step'))
injectionband = cp.get('injections','injectionband')
ntrials = cp.get('injections','ntrials')
confidence = cp.get('upperlimit','confidence')

######################################################################
# create the results directory structure
# be careful if previous results already exist



######################################################################
# create the jobs here so we can sensibly do dependencies
# therefore we have to split up the times and frequencies here 

# get the data set parameters
try: os.system(getdataparams %s %s %

if do_selectdata:

  # split up full data set 


######################################################################
# get the data set params

# get the primary dataset parameters
try: os.system(getdataparams %s %s %)
except:
  print ('failed to run %s on primary full data set', getdataparams) 
  sys.exit(1)

if do_coincidence:
  # get the secondary dataset parameters
  try: os.system(getdataparams %s %s %)
  except:
    print ('failed to run %s on secondary full data set', getdataparams)
    sys.exit(1)
               
######################################################################
# create the selectdata jobs

if do_selectdata:

  # read in the dataset config file
  pdcp = ConfigParser.ConfigParser()
  pdcp.read(fulldata_params)

  # read in primary dataset parameters
  start = (int)(pdcp.get('datasetparams','start'))
  end = (int)(pdcp.get('datasetparams','end'))

  # create the primary selectdata jobs
  start_primary = []
  selectdata_primary = []
  i=0
  time_gap = (int)((float)(end-start)/(float)(selectdata_nodes))
  while i < selectdata_nodes:
    start_primary.append(start+(i*time_gap))
    end_primary.append(start+((i+1)*time_gap))

    # fill in the primary selectdata params
    selectdata_primary.append(binary_pulsar.selectdataNode(selectdata_job))
    selectdata_primary.set_datadir(fulldata_primary)
    selectdata_primary.set_detector(det_primary)
    selectdata_primary.set_tstart(start_primary[i])
    selectdata_primary.set_tend(end_primary[i])
    selectdata_primary.set_outdir(selectdata_secondary_outdir)

    # set primary selectdata dependencies

    # add the primary selectdata nodes
    
    i=i+1
    
  if co_coincidence:

    # read in the dataset config file
    sdcp = ConfigParser.ConfigParser()
    sdcp.read(fulldata_params)

    # read in primary dataset parameters
    start = (int)(sdcp.get('datasetparams','start'))
    end = (int)(sdcp.get('datasetparams','end'))

    # split up secondary dataset 
    start_secondary = []
    selectdata_secondary = []
    i=0
    time_gap = (int)((float)(end-start)/(float)(selectdata_nodes))
    while i < selectdata_nodes:
      start_secondary.append(start+(i*time_gap))
      end_secondary.append(start+((i+1)*time_gap))

      # fill in the secondary selectdata params
      selectdata_secondary.append(binary_pulsar.selectdataNode(selectdata_job))
      selectdata_secondary.set_datadir(fulldata_secondary)
      selectdata_secondary.set_detector(det_secondary)
      selectdata_secondary.set_tstart(start_secondary[i])
      selectdata_secondary.set_tend(end_secondary[i])
      selectdata_secondary.set_outdir(selectdata_secondary_outdir)

      # set primary selectdata dependencies

      # add the primary selectdata nodes
      
      i=i+1


######################################################################
# create the makemesh jobs

if do_search:

  # fill in the primary makemesh params
  makemesh_primary = binary_pulsar.makemeshNode(makemesh_job)
  makemesh_primary.set_datadir = optadata_primary
  makemesh_primary.set_mismatch = mismatch
  makemesh_primary.set_outdir = makemesh_primary_outdir
  
  # add the neccessary prescript arguments 
  makemesh_primary.add_pre_script(getoptdata)
  makemesh_primary.add_pre_script_arg(fulldata_primary)
  makemesh_primary.add_pre_script_arg(optdata_primary)

  # if we are selecting the data use the optimal data params
  if do_selectdata:
    makemesh_primary.add_pre_script_arg(optdata_params_primary)
  # else we use the full data params
  else
    makemesh_primary.add_pre_script_arg(fulldata_params_primary)

  # add the primary makemesh dependencies

  # add the job to the dag

  if do_coincidence:

    # fill in the primary makemesh params
    makemesh_secondary = binary_pulsar.makemeshNode(makemesh_job)
    makemesh_secondary.set_datadir = optadata_secondary
    makemesh_secondary.set_mismatch = mismatch
    makemesh_secondary.set_outdir = makemesh_secondary_outdir
  
    # add the neccessary prescript arguments 
    makemesh_secondary.add_pre_script(getoptdata)
    makemesh_secondary.add_pre_script_arg(fulldata_secondary)
    makemesh_secondary.add_pre_script_arg(optdata_secondary)

    # if we are selecting the data use the optimal data params
    if do_selectdata:
      makemesh_secondary.add_pre_script_arg(optdata_params_secondary)
    # else we use the full data params
    else
      makemesh_secondary.add_pre_script_arg(fulldata_params_secondary)

    # add the secondary makemesh dependencies

    # add the job to the dag

######################################################################
# create the search jobs

if do_search:

  # run script here to produce a search freq/mesh file
  try: os.system(genfreqmeshfile %s %s %s)
  except:
    print 'failed to -----'
    sys.exit(1)

  # split up the frequency bands
  f_min_search = []
  f_band_search = []
  templatebank_primary = []
  templatebank_secondary = []

  search_primary = []
  i=0
  while i < search_nodes:
    search_primary.append(binary_pulsar.searchNode(search_job))
    search_primary.set_f_min(f_min_search[i])
    search_primary.set_f_band(f_band_search[i])
    search_primary.set_binarytemplatebank(templatebank_primary[i])
    search_primary.set_datadir(optdata_primary)
    search_primary.set_detector(det_primary)
    search_primary.set_workdir(search_primary_outdir)

    # add primary search dependencies
    search_primary.add_parent(makemesh_primary)

    # add primary search dag
    dag.add_node(search_primary[i])

  if do_coincidence:

    search_secondary = []
    i=0
    while i < search_nodes:
      search_secondary.append(binary_pulsar.searchNode(search_job))
      search_secondary.set_f_min(f_min_search[i])
      search_secondary.set_f_band(f_band_search[i])
      search_secondary.set_binarytemplatebank(templatebank_secondary[i])
      search_secondary.set_datadir(optdata_secondary)
      search_secondary.set_detector(det_secondary)
      search_secondary.set_workdir(search_secondary_outdir)

      # add secondary search dependencies
      search_secondary.add_parent(makemesh_secondary)

      # add secondary search dag
      dag.add_node(search_secondary[i])

######################################################################
# create the coincidence jobs

if do_coincidence:

  # run script here to produce a coincidence freq/mesh file
  try: os.system(genfreqmeshfile %s %s %s)
  except:
    print 'failed to -----'
    sys.exit(1)

  # split up the frequency bands
  f_min_coincidence = []
  f_band_coincidence = []
  templatebank_primary = []
  templatebank_secondary = []

  coincidence = []
  i=0
  while i < coincidence_nodes:
    coincidence.append(binary_pulsar.coincidenceNode(coincidence_job))
    coincidence.set_presults(search_primary_outdir)
    coincidence.set_sresults(search_secondary_outdir)
    coincidence.set_f_min(f_min_coincidence[i])
    coincidence.set_f_max(f_max_coincidence[i])
    coincidence.set_searchfreqmeshfile(searchfreqmeshfile)
    coincidence.set_sdataparamsfile(sdataparamsfile)
    coincidence.set_outdir(coincidence_outdir)

    # add secondary search dependencies (need to match the frequencies !!)
    
    # add secondary search dag
    dag.add_node(coincidence[i])

######################################################################
# create the injections jobs

if do_injections:

  # need to loop over the frequencies, h0 and trials
  # do this so that we split things up over the requested
  # number of nodes

  # so minimum number of nodes is number of frequencies * number of h0 values

  # generate h0 values
  h0 = []
  nh0 = 1+(math.floor)((h0_max-h0min)/h0step)
  while i < nh0:
    h0.append(h0min+(i*h0step)) # need to have correct formatting here !!!

  # generate frequency values
  try: os.system(genfreqmeshfile %s %s %s)
  except:
    print 'failed to -----'
    sys.exit(1)

  # split up the frequency bands
  f_min_injections = []
  f_band_injections = []
  templatebank_primary = []
  templatebank_secondary = []

  # generate the ntrials value
  ntemp = nh0*nfreq
  n_ntrials = 1 + (math.floor)(injections_nodes/ntemp)
  trials = ntrials/n_ntrials

  injections = []

  # start loop over frequency
  i=0
  while i < nfreq:

    j=0
    while j < nh0:

      k=0
      while k < n_ntrials:

        # add the injection job params
        injections.append(binary_pulsar.injectionsNode(injections_job))
        injections.set_configfile(config_file)
        injections.set_f_min(f_min_injections[i])
        injections.set_f_max(f_max_injections[i])
        injections.set_pdatadir(optdata_primary)
        injections.set_sdatadir(optdata_secondary)
        injections.set_gw_amplitude(h0[i])
        injections.set_ntrials(trials)
        injections.set_outdir(injections_outdir)

        # add injections dependencies

        # add injections job to dag
        dag.add_node(injections[i])

        k=k+1

      j=j+1

    i=i+1  

######################################################################
# create the upperlimit job

if do_upperlimits:

  # set the upperlimit parameters
  upperlimit = binary_pulsar.upperlimitNode(upperlimit_job)
  upperlimit.set_injectionsdir(injections_outdir)
  upperlimit.set_coresultsdir(coincidence_outdir)
  upperlimit.set_injfreqmeshfile(injfreqmeshfile)
  upperlimit.set_injh0file(injh0file)
  upperlimit.set_outfile(upperlimit_outfile)

  # add the upperlimit dependancies

  # add the upperlimits job to the dag
  dag.add_node(upperlimit)

######################################################################
# output the dag











prev_df = None



# write out the DAG
dag.write_sub_files()
dag.write_dag()


sys.exit(0)

