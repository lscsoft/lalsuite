#!/usr/bin/python
"""
FstatFollowUp.py - Uses the coherent ComputeFStatistic_v2 code to follow up on a given parameter space.
"""

#Imports necessary modules
import os,sys,getopt
import numpy
import math
import parseConfigFile

print sys.version

##############################################################
#Initialization of various default parameters for the follow up

Vars = {}

valid_IFOs = ['H1','H2','L1','V1']

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUp.py [options]

  -h, --help               display this message
  -o, --OutputLabel        Base name for the directories and files to be generated
  -i, --Iteration          Iteration number, used to differentiate steps, generally starts at 0
  -C, --ConfigFile         Config file name
  -a, --Alpha              Sky position alpha (equatorial coordinates) in radians
  -d, --Delta              Sky position delta (equatorial coordinates) in radians 
  -z, --AlphaBand          Band in alpha (equatorial coordinates) in radians 
  -c, --DeltaBand          Band in delta (equatorial coordinates) in radians          
  -f, --F0                 Center frequency to be searched (in Hz)
  -b, --F0Band             Total Frequency band to be covered (from -band/2 to +band/2 centered on Freq)
  -s, --F1                 1st spin down central frequency/second to be searched (in Hz/s)
  -m, --F1Band             1st spin down band (from -band/2 to +band/2 centered on F1)
      --F2                 2nd spin down central frequency/second^2 to be searched (in Hz/s^2) 
      --F2Band             2nd spin down band (from -band/2 to +band/2 centered on F2)
      --F3                 3rd spin down central frequency/second^3 to be searched (in Hz/s^3)
      --F3Band             3rd spin down band (from -band/2 to +band/2 centered on F3)
  -T, --ParameterTime      Time which the Freq, dFreq, parameters come from (in GPS seconds)
  -t, --CoherenceTime      Desired span of data actually integrated to look for signals (in seconds)
  -F, --TwoF               Loudest 2F found
  
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "ho:i:C:I:a:d:z:c:f:b:s:m:T:t:F:"
longop = [
  "help",
  "OutputLabel=",
  "Iteration=",
  "ConfigFile",
  "Alpha=",
  "Delta=",
  "AlphaBand=",
  "DeltaBand=",
  "F0=",
  "F0Band=",
  "F1=",
  "F1Band=",
  "F2=",
  "F2Band=",
  "F3=",
  "F3Band=",
  "ParameterTime=",
  "CoherenceTime=",
  "TwoF=",
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:],shortop,longop)
except getopt.GetoptError:
  print getopt.GetoptError
  usage()
  sys.exit(-1)

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-o","--OutputLabel"):
    Vars['base_name'] = str(a)
  elif o in ("-i","--Iteration"):
    Vars['iteration'] = int(a)
  elif o in ("-C","--ConfigFile"):
    config_file_name = str(a)
  elif o in ("-a","--Alpha"):
    Vars['alpha'] = float(a)
  elif o in ("-d","--Delta"):
    Vars['delta'] = float(a)
  elif o in ("-z","--AlphaBand"):
    Vars['alpha_band'] = float(a)
  elif o in ("-c","--DeltaBand"):
    Vars['delta_band'] = float(a)
  elif o in ("-f","--F0"):
    Vars['f0'] = float(a)
  elif o in ("-b","--F0Band"):
    Vars['f0_band'] = float(a)
  elif o in ("-s","--F1"):
    Vars['f1'] = float(a)
  elif o in ("-m","--F1Band"):
    Vars['f1_band'] = float(a)
  elif o in ("--F2"):
    Vars['f2'] = float(a)
  elif o in ("--F2Band"):
    Vars['f2_band'] = float(a)
  elif o in ("--F3"):
    Vars['f3'] = float(a)
  elif o in ("--F3Band"):
    Vars['f3_band'] = float(a)
  elif o in ("-T","--ParameterTime"):
    Vars['param_time'] = int(a)
  elif o in ("-t","--CoherenceTime"):
    Vars['coherence_time'] = float(a)
  elif o in ("-F","--TwoF"):
    loudest2F = float(a)
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(-1)

print loudest2F

####################################
#Parses the config file, if any
#The config file is overriden by command line parameters

Vars = parseConfigFile.parse_config_file(os.path.abspath(config_file_name),Vars)
Vars = parseConfigFile.generate_directory_names(Vars['base_name'],Vars)

Vars = parseConfigFile.generate_file_names(Vars['run_name'],Vars)


  
#The following are derived from "Parameter space metric for combined diurnal and
# orbital motion" by Ian Jones, Ben Owen, and David Whitbeck
#which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls5directed&action=view&page=5
#and also from "Monte-Carlo tests of ComputeFStatistic for 4D parameter-spaces" 
#by Reinhard Prix which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls3knownpulsardemod&action=view&page=2

Vars['freq_resolution'] = 2*(3*Vars['mismatch']/((math.pi**2)*(Vars['coherence_time']**2)))**0.5
Vars['dfreq_resolution'] = 2*(4*5*9*Vars['mismatch']/((math.pi**2)*(Vars['coherence_time']**4)))**0.5

#Determine spacing in right ascension (alpha), declination (delta), 
#using the template grid of the actual code

###########

sky_template_count_file = ''.join([Vars['current_directory'],'/SkyTemplateCount_',Vars['run_name'],'.txt'])
find_total_templates = ' '.join([Vars['code'],
                           '--Alpha', str(Vars['alpha'] - Vars['alpha_band']/2.0),
                           '--Delta', str(Vars['delta'] - Vars['delta_band']/2.0),
                           '--AlphaBand', str(Vars['alpha_band']),
                           '--DeltaBand', str(Vars['delta_band']),
                           '--Freq', str(Vars['f0'] - Vars['f0_band']/2.0),
                           '--FreqBand', str(Vars['f0_band']),
                           '--dFreq', str(Vars['freq_resolution']),   
                           '--f1dot', str(Vars['f1']- Vars['f1_band']/2.0),
                           '--f1dotBand', str(Vars['f1_band']),
                           '--df1dot', str(Vars['dfreq_resolution']),
                           '--DataFiles', ''.join([Vars['final_sft_run_directory'].rstrip('/'),'/*']),
                           '--ephemDir', Vars['ephem_dir'],
                           '--ephemYear', Vars['ephem_year'],
                           '--refTime',str(Vars['param_time']),
                           '--gridType',str(Vars['grid_type']),
                           '--metricType',str(Vars['metric_type']),
                           '--metricMismatch',str(Vars['mismatch']),
                           '--countTemplates',
                           '>',sky_template_count_file])


os.system(find_total_templates)
templates_file = open(sky_template_count_file,'r')
for line in templates_file:
  if line.split(':')[0] == '%% Number of templates':
    sky_grid_templates = float(line.split(':')[1])

#Vars['angular_resolution'] = (Vars['alpha_band']*Vars['delta_band']/sky_grid_templates)**0.5

#Vars['total_templates'] = sky_grid_templates * Vars['f0_band'] * Vars['f1_band'] / (Vars['freq_resolution'] * Vars['dfreq_resolution'])

Vars['total_templates'] = sky_grid_templates

print "Number of templates is: " + str(sky_grid_templates)


F_c = float(loudest2F)/2.0
Prob_below = (1 - ((1 + F_c)/2.0)*math.exp(-F_c))
FA = 1 - Prob_below**Vars['total_templates']

print "The False Alarm probability for 2F of " + str(loudest2F) + " is " + str(FA)


sys.exit(0)
