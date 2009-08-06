#!/usr/bin/python
"""
FollowUpPipeline.py - Main script which calls sub-scripts to perform a coherent follow up to continuous gravitational wave sources
"""

#Imports necessary modules
import os,sys,getopt
import numpy
import math

##############################################################
#Initialization of various default parameters for the follow up

v_max = 30290.0 #m/s
c = 299792458.0 #m/s

SFT_length = 1800 #seconds
time_per_SFT_per_template = 1e-6 #seconds

max_number_of_jobs = 500
number_of_events_to_follow_up = 1


valid_IFOs = ['H1','H2','L1','V1']

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FollowUpPipeline.py [options]

  -h, --help               display this message
  -o, --outputDirectory    Location where the output files will be placed, e.g. ./
  -C, --configFile         File which has default parameters for the follow-ups, i.e. ./followup_config.ini
  -d, --dataFile           File containing candidates to be examined, i.e. ./stage2.csv
  -s, --searchName         Label to be applied to the search files, i.e. Powerflux_S5_outliers
  -L, --logFileDirectory   Location where the condor .sub log file will be placed
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "ho:C:d:s:L:"
longop = [
  "help",
  "outputDirectory=",
  "configFile",
  "dataFile=",
  "searchName=",
  "logFileDirectory=",
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:],shortop,longop)
except getopt.GetoptError:
  print getopt.GetoptError
  usage()
  sys.exit(1)

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-o","--outputDirectory"):
    output_directory = str(a)
  elif o in ("-C","--configFile"):
    config_file = str(a)
  elif o in ("-d","--dataFile"):
    user_data_file = str(a)
  elif o in ("-s","--searchName"):
    search_name = str(a)
  elif o in ("-L","--LogFileDirectory"):
    log_file_directory = str(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

#################################################
#Makes a folder for this search
if os.path.exists(output_directory):
  print "Warning: " + output_directory + " exists but still continuing with DAG generation"
else:
  try:
    os.mkdir(output_directory)
  except:
    print "Error encountered trying to make directory: " + output_directory 
    sys.exit(1)

false_alarm_directory = output_directory.rstrip('/') + 'false_alarms'
if os.path.exists(false_alarm_directory):
  print "Warning: " + false_alarm_directory + " exists but still continuing"
else:
  try:
    os.mkdir(false_alarm_directory)
  except:
    print "Error encountered trying to make directory: " + false_alarm_directory
    sys.exit(1)
    
event_directory = output_directory.rstrip('/') + 'events'
if os.path.exists(event_directory):
  print "Warning: " + event_directory + " exists but still continuing"
else:
  try:
    os.mkdir(event_directory)
  except:
    print "Error encountered trying to make directory: " + event_directory
    sys.exit(1)
###########################################
#Create a log file with the input parameters
run_log_file = ''.join([output_directory.rstrip('/'),'/FollowUpPipeline.log'])
rlf = open(run_log_file,'w')
rlf.write(' '.join(sys.argv))
rlf.write('\n')
rlf.close()

####################################
#Parses the config file, if any
#The config file is overriden by command line parameters

config_file = os.path.abspath(config_file)
sys.path.append(os.path.dirname(config_file))
config = __import__(os.path.basename(config_file).strip('.py'));



# Variables Structure/Dictionary
Vars = {}

# Band in Right Ascension (radians)
try:
  Vars['alpha_band'] = config.alpha_band
except:
  print "alpha_band cannot be read"
  sys.exit(1)
  
# Band in Declination (radians)
try:
  Vars['delta_band'] = config.delta_band
except:
  print "delta_band cannot be read"
  sys.exit(1)
    
# Band in Frequency (Hz)
try:
  Vars['f0_band'] = config.f0_band
except:
  print "f0_band cannot be read"
  sys.exit(1)
      
# Band in 1st Spindown (f dot) (Hz/s)
try:
  Vars['f1_band'] = config.f1_band
except:
  print "f1_band cannot be read"
  sys.exit(1)
      
# Band in 2nd Spindown (f dot dot) (Hz/s^2)
try:
  Vars['f2_band'] = config.f2_band
except:
  print "f2_band cannot be read"
  sys.exit(1)

# Band in 3rd Spindown (f dot dot dot) (Hz/s^3)
try:
  Vars['f3_band'] = config.f3_band
except:
  print "f3_band cannot be read"
  sys.exit(1)

# First GPS time to be considered for data selection
try:
  Vars['start_time'] = config.start_time
except:
  print "start_time cannot be read"
  sys.exit(1)
    
# Last GPS time to be considered for data selection
try:
  Vars['end_time'] = config.end_time
except:
  print "end_time cannot be read"
  sys.exit(1)
      
# Time from which the frequency and spindown values come
try:
  Vars['param_time'] = config.param_time
except:
  print "param_time cannot be read"
  sys.exit(1)

# Maximum 1-D mismatch between template and possible signal
try:
  Vars['mismatch'] = config.mismatch
except:
  print "mismatch cannot be read"
  sys.exit(1)
        
# Earth and Sun Ephemeris directory
try:
  Vars['ephem_dir'] = config.ephem_dir
except:
  print "ephem_dir cannot be read"
  sys.exit(1)
      
# Earth and Sun Ephemeris year to use
try:
  Vars['ephem_year'] = config.ephem_year
except:
  print "ephem_year cannot be read"
  sys.exit(1)

# Number of candidates to keep per job 
try:
  Vars['num_candidates_to_keep'] = config.num_candidates_to_keep
except:
  print "num_candidates_to_keep cannot be read"
  sys.exit(1)

# Grid type used by CommputeFStatistic
try:
  Vars['grid_type'] = config.grid_type
except:
  print "grid_type cannot be read"
  sys.exit(1)

# Metric type used by ComputeFStatistic
try:
  Vars['metric_type'] = config.metric_type
except:
  print "metric_type cannot be read"
  sys.exit(1)    

# Max number of jobs to run per event per step
try:
  Vars['max_number_of_jobs'] = config.max_number_of_jobs
except:
  print "max_number_of_jobs cannot be read"
  sys.exit(1)

# Code to be used for search
try:
  Vars['code'] = config.code
except:
  print "code (the configuration variable) cannot be read"
  sys.exit(1)

# Directory containing the python scripts
try:
  Vars['python_dir'] = config.python_dir
except:
  print "python_dir cannot be read"
  sys.exit(1)

# Initial coherence time
try:
  Vars['start_coherence_time'] = config.start_coherence_time
except:
  print "start_coherence_time cannot be read"
  sys.exit(1)
  

#############################
#End user input section

  
###############################################################################
#This function parses the data file and generates the .dag file which is submitted on the cluster

def parse_Powerflux(data_file):
  main_dag_file_name = ''.join([output_directory.rstrip('/'),'/',search_name,'.dag'])
  main_sub_file_name = ''.join([output_directory.rstrip('/'),'/',search_name,'.sub'])

  dag_file = open(main_dag_file_name,'w')

  data = open(data_file,'r')
  
  dag_file.write('# Dag for base_name = ' + search_name + '\n')
  job = 0
  for line in data:
    if line[0] != '"':
      alpha = line.split('\t')[22]
      delta = line.split('\t')[23]
      f0 = line.split('\t')[20]
      f1 = line.split('\t')[21]
      f2 = 0
      f3 = 0
      outlier_number = line.split('\t')[0]
      base_name = search_name + '_' + outlier_number

      arg_list = ' '.join(['argList="',
                           '-o',base_name,
                           '-C',config_file,
                           '-i','0',
                           '-a',alpha,
                           '-d',delta,
                           '-z',str(Vars['alpha_band']),
                           '-c', str(Vars['delta_band']),
                           '-f', f0,
                           '-b',str(Vars['f0_band']),
                           '-T', str(Vars['param_time']),
                           '-s', f1,
                           '-m', str(Vars['f1_band']),
                           '--F2', str(f2),
                           '--F2Band', str(Vars['f2_band']),
                           '--F3', str(f3),
                           '--F3Band', str(Vars['f3_band']),
                           '-t',str(Vars['start_coherence_time']),
                           '-L', str(log_file_directory),
                           '"'])
      
      
      dag_file.write('JOB A' + str(job) + ' ' + main_sub_file_name +'\n')
      dag_file.write('VARS A' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + arg_list + '\n')
      dag_file.write('\n')
      job += 1

  dag_file.close()

  ########################################################
  #This step creates the .sub file(s) called by the .dag file

  sub_file = open(main_sub_file_name,'w')

  #This sub file is used to perform the actual search
  sub_file.write('universe= vanilla \n')
  sub_file.write('executable = ' + Vars['python_dir'].rstrip('/') + '/' + 'FstatFollowUp.py' +'\n')
  sub_file.write('output = node_' + search_name + '_A.out.$(JobID)\n')
  sub_file.write('error = node_' + search_name + '_A.err.$(JobID)\n')
  sub_file.write('log = ' + log_file_directory.rstrip('/') + '/' + search_name + '.log\n' )
  sub_file.write('\n')
  sub_file.write('arguments = $(arglist)\n')
  sub_file.write('queue\n')

  sub_file.close()


###################################
#Calls the correct function to parse the data file.
#Only one type of file and function at the moment

parse_Powerflux(user_data_file)


######################################
#This step submits the job to Condor


