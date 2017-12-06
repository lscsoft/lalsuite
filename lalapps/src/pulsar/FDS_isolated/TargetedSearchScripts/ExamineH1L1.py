#!/usr/bin/python

"""
ExamineFollowUpResults.py - Examines the results from a FStatFollowUp.py run, finds the largest events, makes a few plots, and calls FStatFollowUp.py if another iteration is still possible. 
"""

#Imports necessary modules
import os,sys,getopt
from numpy import *
import math
import time
import parseConfigFile

###################
#Initialization
Vars = {}

valid_IFOs = ['H1','L1']


################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUpExamine.py [options]

  -h, --help               display this message
  -C, --configFile         File which has default parameters for the follow-ups
  -o, --OutputName         Base name for the directory
  -j, --TotalJobs          Total number of jobs in search
  -i, --Iteration          Indicates which run this is
  -t, --CoherentTime       Coherent time of the jobs just run
  -T, --Time               Time of calculated ephemeris (start time of search)
  -A, --AngularResolution  Spacing in Right Ascension or Declination (radians) of templates from previous iteration
  -N, --TemplateNumber     Total number of templates in the previous iteration
"""
  print >> sys.stderr, msg


#################################################
#This section reads in the command line arguments

shortop = "hC:o:j:d:b:i:t:T:A:N:"
longop = [
  "help",
  "ConfigFile=",
  "OutputDirectory=",
  "TotalJobs=",
  "Iteration=",
  "CoherentTime=",
  "Time=",
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
    config_file_name = str(a)   
  elif o in ("-o","--OutputDirectory"):
    Vars['base_name'] = str(a)
  elif o in ("-j","--TotalJobs"):
    Vars['jobs'] = int(a)
  elif o in ("-i","--Iteration"):
    Vars['iteration'] = int(a)
  elif o in ("-t","--CoherentTime"):
    Vars['coherence_time'] = float(a)
  elif o in ("-T","Time"):
    Vars['param_time'] = int(float(a))
  elif o in ("-A","AngularResolution"):
    Vars['angular_resolution'] = float(a)
  elif o in ("-N","TemplateNumber"):
    Vars['total_templates'] = float(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

###################################
#Parses the config file, if any
#The config file is overriden by command line parameters

Vars = parseConfigFile.parse_config_file(os.path.abspath(config_file_name),Vars)
Vars = parseConfigFile.generate_directory_names(Vars['base_name'],Vars)

Vars = parseConfigFile.generate_file_names(Vars['run_name'],Vars)


###############################3
#Log the command line input

logFile = open(''.join([Vars['current_directory'].rstrip('/'),'/',Vars['run_name'],'_examine.log']),'w')
logFile.write(' '.join(sys.argv))
logFile.write('\n')
logFile.close()

###################################
#Checks to see if all the results files exist
#If some results files are missing, throw an error and exit

###################################
#Find the loudest

for this_IFO in valid_IFOs:

  if this_IFO == 'H1':
    job_offset = 0
  elif this_IFO == 'L1':
    job_offset = Vars['jobs']

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

  recordLoudestResultFile = ''.join([Vars['current_directory'],'/',this_IFO,'_loudestResult_',str(Vars['iteration']),'.txt'])
  record_loudest_results_file = open(recordLoudestResultFile,'w')

  for job in range(0+job_offset,Vars['jobs']+job_offset):
    loudestFile = open(''.join([Vars['output_run_directory'].rstrip('/'),'/',this_IFO,'_',Vars['loudest_base_prefix'],str(job)]),'r')
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

    record_loudest_results_file.write(' '.join([str(tempFreq),str(tempF1dot),str(tempF2dot),str(tempF3dot),str(tempRA),str(tempDEC),str(temp2F),"\n"]))

    if float(temp2F) > float(loudest2F):
        loudest2F = temp2F
        loudestRA = tempRA
        loudestDEC = tempDEC
        loudestFreq = tempFreq
        loudestF1dot = tempF1dot
        loudestF2dot = tempF2dot
        loudestF3dot = tempF3dot

  #Finds all 2F values greater than or equal to a percentage (i.e. 90%) of the loudest 2F value, and
  #takes a weighted average to find the most likely source position in frequency, spindown and sky position
  #parameter space

  average_cutoff = float(loudest2F) * float(Vars['percent_largest_to_average'])



  for job in range(0,Vars['jobs']):
    fullResultsFile = open(''.join([Vars['output_run_directory'],Vars['results_base_prefix'],str(job)]),'r')
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

  Vars['Coherence_time'] = float(Vars['coherence_time'])
  F_band = Vars['parameter_space_multiplier']*1.0/Vars['Coherence_time']
  F_dot_band = Vars['parameter_space_multiplier']*1.0/Vars['Coherence_time']**2

  RA_band = Vars['parameter_space_multiplier']*Vars['angular_resolution']
  DEC_band = Vars['parameter_space_multiplier']*Vars['angular_resolution']

  Vars['new_coherence_time'] = Vars['Coherence_time'] * 4.0


  F_c = float(loudest2F)/2.0
  Prob_below = (1 - ((1 + F_c)/2.0)*math.exp(-F_c))
  FA = 1 - Prob_below**Vars['total_templates']


  single_false_alarm_file = str(Vars['false_alarm_file']) + "." + str(this_IFO)

  if os.path.exists(single_false_alarm_file):
    false_alarm_file_exists = True
  else:
    false_alarm_file_exists = False
  
  FA_file = open(single_false_alarm_file,'a')
  if not false_alarm_file_exists:
    FA_file.write('Param_time Coherence_Time 2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA Name Run\n')
  FA_file.write(' '.join([str(Vars['param_time']),str(Vars['coherence_time']),loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(Vars['total_templates']),str(FA),str(Vars['base_name']),str(Vars['iteration']),'\n']))
  FA_file.close()

  time.sleep(1)

  #Cat result files for easy printing
  os.system('cat ' + Vars['output_run_directory'].rstrip('/') + '/' + str(this_IFO) + '_' + Vars['results_base_prefix'] + '* > ' + Vars['events_directory'].rstrip('/') + '/' + str(this_IFO) + '_plot_data_' + Vars['run_name'])

  #Delete the previous run to save space
  if not Vars['messy']:
    os.system(''.join(['rm -rf ',Vars['base_name'],'_',str(Vars['iteration']),'_run']))

  #Delete condor log files to save space
  if not Vars['messy']:
    os.system(''.join(['rm node_',Vars['base_name'],'*']))


sys.exit(0)
