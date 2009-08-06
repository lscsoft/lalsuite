#!/usr/bin/python
"""
ExamineFollowUpResults.py - Examines the results from a FStatFollowUp.py run, finds the largest events, makes a few plots, and calls FStatFollowUp.py if another iteration is still possible. 
"""

#Imports necessary modules
import os,sys,getopt
from numpy import *
import math

##############################################################
#Initialization of various default parameters for the follow up

v_max = 30290.0 #m/s
c = 299792458.0 #m/s

#Multiplies the estimated resolution to get parameter space to search
param_space_mult = 10.0

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUpExamine.py [options]

  -h, --help               display this message
  -C, --configFile         File which has default parameters for the follow-ups
  -o, --OutputDirectory    Base name for the directory
  -j, --TotalJobs          Total number of jobs in search
  -d, --InputDirectory     Directory where jobs are stored
  -b, --BaseName           Base file name of the jobs to be examined (no result/loudest, just base)
  -i, --Iteration          Indicates which run this is
  -t, --Tcoh               Coherent time of the jobs just run
  -T, --Time               Time of calculated ephemeris (start time of search)
  -L, --UserLogDirectory   Location to place condor log files
  -A, --AngularResolution  Spacing in Right Ascension or Declination (radians) of templates from previous iteration
  -N, --TemplateNumber     Total number of templates in the previous iteration
"""
  print >> sys.stderr, msg


#################################################
#This section reads in the command line arguments

shortop = "hC:o:j:d:b:i:t:T:L:A:N:"
longop = [
  "help",
  "ConfigFile=",
  "OutputDirectory=",
  "TotalJobs=",
  "InputDirectory=",
  "BaseName=",
  "Iteration=",
  "Tcoh=",
  "Time=",
  "UserLogDirectory=",
  "AngularResolution=",
  "TemplateNumber=",
  ]

#try:
opts, args = getopt.getopt(sys.argv[1:],shortop,longop)
#except getopt.GetoptError:
#  print "Getopt error"
#  usage()
#  sys.exit(1)

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-C","--ConfigFile"):
    config_file = str(a)   
  elif o in ("-o","--OutputDirectory"):
    output_dir = str(a)
  elif o in ("-j","--TotalJobs"):
    jobs = int(a)
  elif o in ("-d","--InputDirectory"):
    input_dir = str(a)
  elif o in ("-b","--BaseName"):
    base_name = str(a)
  elif o in ("-i","--Iteration"):
    run = int(a)
  elif o in ("-t","--Tcoh"):
    tcoh = float(a)
  elif o in ("-T","Time"):
    param_time = int(float(a))
  elif o in ("-L","UserLogDirectory"):
    user_log = str(a)
  elif o in ("-A","AngularResolution"):
    angular_resolution = float(a)
  elif o in ("-N","TemplateNumber"):
    total_templates = float(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

logFile = open(''.join([base_name,'_',str(run),'_examine.log']),'w')
logFile.write(' '.join(sys.argv))
logFile.write('\n')
logFile.close()


###################################
#Parses the config file, if any
#The config file is overriden by command line parameters

# Import the Configuration File into config
config_file = os.path.abspath(config_file)
sys.path.append(os.path.dirname(config_file))
config = __import__(os.path.basename(config_file.strip('.py')));

# First GPS time to be considered for data selection
try:
  config.start_time
except:
  print "start_time cannot be read"
  sys.exit(1)
    
# Last GPS time to be considered for data selection
try:
  config.end_time
except:
  print "end_time cannot be read"
  sys.exit(1)
      
# Maximum 1-D mismatch between template and possible signal
try:
  config.mismatch
except:
  print "mismatch cannot be read"
  sys.exit(1)
        
# Earth and Sun Ephemeris directory
try:
  config.ephem_dir
except:
  print "ephem_dir cannot be read"
  sys.exit(1)
      
# Earth and Sun Ephemeris year to use
try:
  config.ephem_year
except:
  print "ephem_year cannot be read"
  sys.exit(1)

# Number of candidates to keep per job 
try:
  config.num_candidates_to_keep
except:
  print "num_candidates_to_keep cannot be read"
  sys.exit(1)

# Grid type used by CommputeFStatistic
try:
  config.grid_type
except:
  print "grid_type cannot be read"
  sys.exit(1)

# Metric type used by ComputeFStatistic
try:
  config.metric_type
except:
  print "metric_type cannot be read"
  sys.exit(1)    

# Max number of jobs to run per event per step
try:
  config.max_number_of_jobs
except:
  print "max_number_of_jobs cannot be read"
  sys.exit(1)

# Location of data files, in comma seperated format
try:
  config.data_location_file = config.data_location_file.split(',')
except:
  print "data_location_file cannot be read"
  sys.exit(1)

# Location of the python files
try:
  config.python_dir
except:
  print "python_dir cannot be read"
  sys.exit(1)

# Which % of the largest 2F values to average to get the location of the signal
try:
  config.percent_largest_to_average
except:
  print "percent_largest_to_average cannot be read"
  sys.exit(1)

###################################
#Checks to see if all the results files exist
#If some results files are missing, throw an error and exit

###################################
#Find the loudest

input_dir = input_dir.rstrip('/')
input_dir = ''.join([input_dir,'/'])
output_dir = output_dir.rstrip('/') + '/'
new_base_name = ''.join([base_name,'_',str(run),'_loudest_'])
new_results_base_name = ''.join([base_name,'_',str(run),'_result_'])


loudest2F = 0
loudestFreq = 0
loudestF1dot = 0
loudestF2dot = 0
loudestF3dot = 0
loudestRA = 0
loudestDEC = 0

array2F = array([])
arrayFreq = array([])
arrayRA = array([])
arrayDEC = array([])
arrayF1dot = array([])
arrayF2dot = array([])
arrayF3dot = array([])

loudestResultFile = ''.join(['./' + 'loudestResult_',str(run),'.txt'])
lrf = open(loudestResultFile,'w')

for job in range(0,jobs):
    loudestFile = open(''.join([input_dir,new_base_name,str(job)]),'r')
    for line in loudestFile:
        if line.split('=')[0].strip() == 'twoF':
            temp2F = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Alpha':
            tempRA = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Delta':
            tempDEC = line.split('=')[1].strip().strip(';')
            
        elif line.split('=')[0].strip() == 'Freq':
            tempFreq = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f1dot':
            tempF1dot = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f2dot':
            tempF2dot = line.split('=')[1].strip().strip(';')

        elif line.split('=')[0].strip() == 'f3dot':
            tempF3dot = line.split('=')[1].strip().strip(';')

    loudestFile.close()

    lrf.write(' '.join([str(tempFreq),str(tempF1dot),str(tempF2dot),str(tempF3dot),str(tempRA),str(tempDEC),str(temp2F),"\n"]))

    if float(temp2F) > float(loudest2F):
        loudest2F = temp2F
        loudestRA = tempRA
        loudestDEC = tempDEC
        loudestFreq = tempFreq
        loudestF1dot = tempF1dot
        loudestF2dot = tempF2dot
        loudestF3dot = tempF3dot

#Finds all 2F values greater than or equal to a percentage (90%) of the loudest 2F value, and
#takes a weighted average to find the most likely source position in frequency, spindown and sky position
#parameter space

average_cutoff = float(loudest2F) * float(config.percent_largest_to_average)



for job in range(0,jobs):
  fullResultsFile = open(''.join([input_dir,new_results_base_name,str(job)]),'r')
  for line in fullResultsFile:
    if not (line[0] == '%'):
      #Tests to see if 2F greater than required cutoff, then records data if true
      if float(line.split(' ')[6].strip()) >= average_cutoff:
        array2F= append(array2F,float(line.split(' ')[6].strip()))
        arrayFreq = append(arrayFreq, float(line.split(' ')[0].strip()))
        arrayRA = append(arrayRA,float(line.split(' ')[1].strip()))
        arrayDEC = append(arrayDEC,float(line.split(' ')[2].strip()))
        arrayF1dot = append(arrayF1dot,float(line.split(' ')[3].strip()))
        arrayF2dot = append(arrayF2dot,float(line.split(' ')[4].strip()))
        arrayF3dot = append(arrayF3dot,float(line.split(' ')[5].strip()))

#Calculate the weighted mean of all the parameters (weighted by 2F value)
total2Fsum = sum(array2F)
meanFreq = sum(array2F*arrayFreq)/total2Fsum
meanRA = sum(array2F*arrayRA)/total2Fsum
meanDEC = sum(array2F*arrayDEC)/total2Fsum
meanF1dot = sum(array2F*arrayF1dot)/total2Fsum
meanF2dot = sum(array2F*arrayF2dot)/total2Fsum
meanF3dot = sum(array2F*arrayF3dot)/total2Fsum

#Use the calculated average parameters to generate the parameter space for the next longer coherent step

Tcoh = float(tcoh)
F_band = param_space_mult*1.0/Tcoh
F_dot_band = param_space_mult*1.0/Tcoh**2

RA_band = param_space_mult*angular_resolution
DEC_band = param_space_mult*angular_resolution

new_Tcoh = Tcoh * 4.0


F_c = float(loudest2F)/2.0
Prob_below = (1 - ((1 + F_c)/2.0)*math.exp(-F_c))
FA = 1 - Prob_below**total_templates



FA_file = open(''.join(['../false_alarms/FA_',base_name,'.txt']),'a')
if run == 0:
  FA_file.write('2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA Name Run\n')
FA_file.write(' '.join([loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(total_templates),str(FA),str(base_name),str(run)]))
FA_file.close()

sleep(1)

#Delete the previous run to save space
os.system(''.join(['rm -rf ',base_name,'_',str(run),'_run']))

#Delete condor log files to save space
os.system(''.join(['rm node_',base_name,'*']))


if FA < 0.01:
  interesting_event_file = open(''.join(['../events/',base_name,'_event_',str(run),'.txt']),'w')
  interesting_event_file.write('2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA\n')
  interesting_event_file.write(' '.join([loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(total_templates),str(FA)]))
  interesting_event_file.close()

  
  if float(Tcoh) >= (config.end_time - config.start_time):
    print "Last search used all available data - no further iterations"
    sys.exit(0)

  code = config.python_dir.rstrip('/') + '/FstatFollowUp.py '


  #Command string for the next iteration to run
  command_string = ' '.join([code,
                             '-o',base_name,
                             '-C',config_file,
                             '-i',str(run+1),
                             '-a',str(meanRA),
                             '-d',str(meanDEC),
                             '-z',str(RA_band),
                             '-c', str(DEC_band),
                             '-f', str(meanFreq),
                             '-b',str(F_band),
                             '-T', str(param_time),
                             '-s', str(meanF1dot),
                             '-m', str(F_dot_band),
                             '--F2', str(0),
                             '--F2Band', str(0),
                             '--F3', str(0),
                             '--F3Band', str(0),
                             '-t',str(new_Tcoh),
                             '-L', str(user_log)])
  
  os.system(command_string)
