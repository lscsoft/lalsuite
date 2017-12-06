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
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "ho:i:C:I:a:d:z:c:f:b:s:m:T:t:"
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
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(-1)


####################################
#Parses the config file, if any
#The config file is overriden by command line parameters

Vars = parseConfigFile.parse_config_file(os.path.abspath(config_file_name),Vars)
Vars = parseConfigFile.generate_directory_names(Vars['base_name'],Vars)

Vars = parseConfigFile.generate_file_names(Vars['run_name'],Vars)


#############################################
#Calculate some derived parameters and names
#Also performs some basic sanity checks

Vars['total_time'] = Vars['end_time'] - Vars['start_time']

if Vars['total_time'] < 0:
  print "No valid times: check your start and end times"
  sys.exit(-1)

#Need to double check this calculation
#Performs a coarse estimate

if Vars['coherence_time'] > Vars['total_time']:
  Vars['coherence_time'] = Vars['total_time']


##########################################################
#In this step, we loop through each IFO requested,
#examine the data available, use noise waiting to figure out
#which data to use, find the data on the cluster,
#strip out everything but a band around the source, and concatenate 
#the data into a single file per IFO (to ease I/O on the cluster) 

#Opens the data location file

Vars['time_array'] = numpy.array([])
Vars['psd_array'] = numpy.array([])
for this_IFO in valid_IFOs:
  found_file = False
  if this_IFO in Vars['IFOs']:
    for file_name in Vars['data_location_file']:
      if this_IFO in file_name.split('/')[-1]:
        try:
          data_file = open(file_name,'r')
        except:
          print "Unable to open " + file_name
          sys.exit(-1)
        found_file = True
    if found_file == False:
        print "I don't know which file to use for data location for " + this_IFO
        print "Try adding " + this_IFO + " somewhere to the file name"
        sys.exit(-1)

    #If a file containing the noise in band per sft does not exist,
    #create it
    psd_file_name = ''.join([Vars['current_directory'].rstrip('/'),'/Full_psd_',os.path.basename(Vars['current_directory']),'_',str(this_IFO),'.txt'])

    print psd_file_name
    if (not os.path.exists(psd_file_name)):
      print "Could not find " + psd_file_name
      
    data_file.close()
        
    #For each sft, note its start time in an array
    #Also note its psd for the frequency band
 
    try:
      psd_file = open(psd_file_name,'r')
    except:
      print "Unable to open " + psd_file_name
      sys.exit(-1)
    
    for line in psd_file:
      sft_start_time = int(line.split(' ')[0])
      if ((sft_start_time > Vars['start_time']) and
          (sft_start_time + Vars['SFT_length'] < Vars['end_time'])):
        Vars['time_array'] = numpy.append(Vars['time_array'],sft_start_time)
        psd_value = float(line.split(' ')[2])
        Vars['psd_array'] = numpy.append(Vars['psd_array'],psd_value)

    psd_file.close()
    #Determine the best section of data to use based on the number of SFTs
    #within a fixed duration, determined by the number of expected templates
      
############################################################################################
#This section determines the best set of SFTs to use given the noise of each SFT in the band

Vars['best_start_time'] = Vars['start_time']
best_noise_weight = 1
Vars['invert_psd_array'] = 1./Vars['psd_array']
first_time_pass = True
for x in Vars['time_array']:
  if numpy.sum((Vars['time_array'] > x) & (Vars['time_array'] < x + Vars['coherence_time'])) > 0:
    sum_invert_psd = numpy.sum(Vars['invert_psd_array'] * ((Vars['time_array'] > x) & (Vars['time_array'] < x + Vars['coherence_time'])))
    noise_weight = (1/sum_invert_psd)**0.5
    if (noise_weight < best_noise_weight) or (first_time_pass):
      best_noise_weight = noise_weight
      Vars['best_start_time'] = x
      first_time_pass = False

print best_noise_weight

sys.exit(0)
