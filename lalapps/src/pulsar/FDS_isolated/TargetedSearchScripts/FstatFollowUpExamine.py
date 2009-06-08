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
  -o, --OutputDirectory    Base name for the directory
  -j, --TotalJobs          Total number of jobs in search
  -d, --InputDirectory     Directory where jobs are stored
  -b, --BaseName           Base file name of the jobs to be examined (no result/loudest, just base)
  -i, --Iteration          Indicates which run this is
  -t, --Tcoh               Coherent time of the jobs just run
  -D, --DataFileLocation   Location of lists of sft files on nodes
  -I, --IFOs               IFOs used in search
  -T, --Time               Time of calculated ephemeris (start time of search)
  -e, --EphemerisDirectory Directory containing Sun and Earth ephemeris files
  -y, --EphemerisYear      Year of ephemeris files needed, i.e. "03-06"
  -L, --UserLogDirectory   Location to place condor log files
  -C, --SearchCode         Name and location of search code, i.e. ./lalapps_ComputeFStatisic_v2
  -A, --AngularResolution  Spacing in Right Ascension or Declination (radians) of templates from previous iteration
  -N, --TemplateNumber     Total number of templates in the previous iteration
"""
  print >> sys.stderr, msg


#################################################
#This section reads in the command line arguments

shortop = "ho:j:d:b:i:t:D:I:T:e:y:L:C:A:N:"
longop = [
  "help",
  "OutputDirectory=",
  "TotalJobs=",
  "InputDirectory=",
  "BaseName=",
  "Iteration=",
  "Tcoh=",
  "DataFileLocation=",
  "IFOs=",
  "Time=",
  "EphemerisDirectory=",
  "EphemerisYear=",
  "UserLogDirectory=",
  "SearchCode=",
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
  elif o in ("-D","DataFileLocation"):
    list_files = str(a)
  elif o in ("-I","IFOs"):
    ifos = str(a)
  elif o in ("-T","Time"):
    param_time = int(float(a))
  elif o in ("-e","EphemerisDirectory"):
    ephem_dir = str(a)
  elif o in ("-y","EphemerisYear"):
    ephem_year = str(a)
  elif o in ("-L","UserLogDirectory"):
    user_log = str(a)
  elif o in ("-C","SearchCode"):
    search_code = str(a)
  elif o in ("-A","AngularResolution"):
    angular_resolution = float(a)
  elif o in ("-N","TemplateNumber"):
    total_templates = float(a)
  else:
    print >> sys.stderr, "Unkonwon option:", o
    usage()
    sys.exit(1)

logFile = open(''.join([base_name,'_',str(run),'_examine.log']),'w')
logFile.write(' '.join(sys.argv))
logFile.write('\n')
logFile.close()




###################################
#Parses the config file, if any
#The config file is overriden by command line parameters


###################################
#Checks to see if all the results files exist
#If some results files are missing, throw an error and exit

###################################
#Find the loudest

input_dir = input_dir.rstrip('/')
input_dir = ''.join([input_dir,'/'])
output_dir = output_dir.rstrip('/')
output_dir = ''.join([output_dir,'/'])
new_base_name = ''.join([base_name,'_',str(run),'_loudest_'])


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

loudestResultFile = ''.join([output_dir,'loudestResult_',str(run),'.txt'])
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

Tcoh = float(tcoh)
F_band = param_space_mult*1.0/Tcoh
F_dot_band = param_space_mult*1.0/Tcoh**2

RA_band = param_space_mult*angular_resolution
DEC_band = param_space_mult*angular_resolution


F_c = float(loudest2F)/2.0
Prob_below = (1 - ((1 + F_c)/2.0)*math.exp(-F_c))
FA = 1 - Prob_below**total_templates

if FA < 0.01:
  event_file = open(''.join([base_name,'_event_',str(run),'.txt']),'w')
  event_file.write('2F RA DEC Freq F1dot F2dot F3dot Total-Templates FA\n')
  event_file.write(' '.join([loudest2F,loudestRA,loudestDEC,loudestFreq,loudestF1dot,loudestF2dot,loudestF3dot,str(total_templates),str(FA)]))
  event_file.close()
  


code = './FstatFollowUp.py '
start_time = 0
end_time = 999999999
mismatch = 0.15
cutoff2F = 30

os.system(''.join(['mkdir node-',str(run)]))
os.system(''.join(['mv node_',base_name,'* ./node-',str(run)]))

os.system('source /opt/lscsoft/lalapps/etc/lalapps-user-env.sh')

command_string = ''.join([code,' -o ',base_name,' -i ',str(run+1),
                          ' -a ',str(loudestRA),' -d ',str(loudestDEC),
                          ' -z ', str(RA_band), ' -c ', str(DEC_band),
                          ' -f ', str(loudestFreq), ' -b ', str(F_band),
                          ' -D ', str(list_files), ' -I ', str(ifos),
                          ' -S ', str(start_time), ' -E ', str(end_time),
                          ' -T ', str(param_time), ' -s ', str(loudestF1dot),
                          ' -m ', str(F_dot_band), ' -M ', str(mismatch),
                          ' -e ', str(ephem_dir), ' -y ', str(ephem_year),
                          ' -F ', str(cutoff2F), ' -C ', str(search_code),
                          ' -L ', str(user_log)])

os.system(command_string)
