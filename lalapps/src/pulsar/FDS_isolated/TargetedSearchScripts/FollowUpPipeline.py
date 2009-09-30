#!/usr/bin/python
"""
FollowUpPipeline.py - Main script which calls sub-scripts to perform a coherent follow up to continuous gravitational wave sources
"""

#Imports necessary modules
import os,sys,getopt
import numpy
import math
import parseConfigFile

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
  -C, --configFile         File which has default parameters for the follow-ups, i.e. ./followup_config.ini
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "hC:"
longop = [
  "help",
  "configFile",
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
  elif o in ("-C","--configFile"):
    config_file = str(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)


########
# Parse the config file

Vars = {}
Vars = parseConfigFile.parse_config_file(os.path.abspath(config_file),Vars)

Vars = parseConfigFile.generate_directory_names('',Vars)

Vars = parseConfigFile.generate_file_names(Vars['search_name'],Vars)
    

#################################################
#Makes a folder for this search
if os.path.exists(Vars['work_dir']):
  print "Warning: " + Vars['work_dir'] + " exists but still continuing with DAG generation"
else:
  try:
    os.mkdir(Vars['work_dir'])
  except:
    print "Error encountered trying to make directory: " + Vars['work_dir'] 
    sys.exit(1)

if os.path.exists(Vars['false_alarm_results_directory']):
  print "Warning: " + Vars['false_alarm_results_directory'] + " exists but still continuing"
else:
  try:
    os.mkdir(Vars['false_alarm_results_directory'])
  except:
    print "Error encountered trying to make directory: " + Vars['false_alarm_results_directory']
    sys.exit(1)
    
if os.path.exists(Vars['events_directory']):
  print "Warning: " + Vars['events_directory'] + " exists but still continuing"
else:
  try:
    os.mkdir(Vars['events_directory'])
  except:
    print "Error encountered trying to make directory: " + Vars['events_directory']
    sys.exit(1)

###########################################
#Create a log file with the input parameters

rlf = open(Vars['self_log_file_name'],'w')
rlf.write(' '.join(sys.argv))
rlf.write('\n')
rlf.close()
  

#############################
#End user input section

  
###############################################################################
#This function parses the data file and generates the .dag file which is submitted on the cluster

def parse_Powerflux(Vars):
  
  dag_file = open(Vars['self_dag_file_name'],'w')

  data = open(Vars['parameter_source_file'],'r')
  
  dag_file.write('# Dag for base_name = ' + Vars['search_name'] + '\n')
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
      base_name = Vars['search_name'] + '_' + outlier_number

      arg_list = ' '.join(['argList="',
                           '-o',base_name,
                           '-C',Vars['config_file'],
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
                           '-t',str(Vars['coherence_time']),
                           '"'])
      
      
      dag_file.write('JOB A' + str(job) + ' ' + Vars['self_sub_file_name'] +'\n')
      dag_file.write('VARS A' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + arg_list + '\n')
      dag_file.write('\n')

      dag_file_name = Vars['work_dir'].rstrip('/') + '/' + base_name + '/' + base_name + '_0.dag'
      resubmit_arg_list = ''.join(['argList=" ',
                                 '-d ',os.path.dirname(dag_file_name),
                                 ' -l ',Vars['submit_file_log_dir'],
                                 ' -f ',dag_file_name,
                                 '"'])

      dag_file.write('JOB B'+ str(job) + ' ' + Vars['self_resubmit_file_name'] +'\n')
      dag_file.write('VARS B' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + resubmit_arg_list + '\n')
      dag_file.write('\n')
      
      job += 1


  for z in range(0,job):
    dag_file.write('PARENT A' + str(z) + ' ' + 'CHILD B'+str(z) + '\n')

  dag_file.close()

  ########################################################
  #This step creates the .sub file(s) called by the .dag file

  sub_file = open(Vars['self_sub_file_name'],'w')

  #This sub file is used to perform the actual search
  sub_file.write('universe= vanilla \n')
  sub_file.write('executable = ' + Vars['python_dir'].rstrip('/') + '/' + 'FstatFollowUp.py' +'\n')
  sub_file.write('output = node_' + Vars['search_name'] + '_A.out.$(JobID)\n')
  sub_file.write('error = node_' + Vars['search_name'] + '_A.err.$(JobID)\n')
  sub_file.write('log = ' + Vars['self_sub_log_file_name'])
  sub_file.write('\n')
  sub_file.write('arguments = $(arglist)\n')
  sub_file.write('queue\n')

  sub_file.close()

  #This subfile handles the scripts run in the scheduler universe so more dags can be run
  resub_file = open(Vars['self_resubmit_file_name'],'w')
  
  resub_file.write('universe= scheduler \n')
  resub_file.write('getenv= true \n')
  resub_file.write('executable = ' + Vars['python_dir'].rstrip('/') + '/' + 'FstatFollowUpSubmit.py\n')
  resub_file.write('output = node_submit_' + Vars['search_name'] + '_B.out.$(JobID)\n')
  resub_file.write('error = node_submit_' + Vars['search_name'] + '_B.err.$(JobID)\n')
  resub_file.write('log = ' + Vars['self_sub_log_file_name'])
  resub_file.write('\n')
  resub_file.write('arguments = $(arglist)\n')
  resub_file.write('queue\n')

  resub_file.close()
  

###################################
#Calls the correct function to parse the data file.
#Only one type of file and function at the moment

parse_Powerflux(Vars)

sys.exit(0)


