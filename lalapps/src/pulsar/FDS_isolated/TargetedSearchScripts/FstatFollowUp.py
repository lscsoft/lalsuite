#!/usr/bin/python
"""
FstatFollowUp.py - Uses the coherent ComputeFStatistic_v2 code to follow up on a given parameter space.
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

grid_type = 2
metric_type = 1

valid_IFOs = ['H1','H2','L1','V1']

alpha = None
delta = None
alpha_band = None
delta_band = None
freq = None
f_band = None
final_iteration = False

################################
#Prints the usage of this script
def usage():
  msg = """\
Usage: FstatFollowUp.py [options]

  -h, --help               display this message
  -o, --OutputLabel        Base name for the directories and files to be generated
  -i, --Iteration          Iteration number, used to differentiate steps, generally starts at 0
  -C, --ConfigFile         Config file name
  -D, --DataFileLocation   Location(s) of data file indicating location of SFTs on the cluster, seperated by commas
                           i.e. "./H1ListOfDataOnNodes,./H2ListOfDataOnNodes,./L1ListOfDataOnNodes"
  -I, --IFOs               Which IFOs: 'H1','H2','L1','H1H2','H1L1','H2L1','H1H2L1', etc
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
  -S, --StartTime          Earliest time that a follow can use data from (in GPS seconds)
  -E, --EndTime            Latest time that a follow up can use data from (in GPS seconds)
  -T, --ParameterTime      Time which the Freq, dFreq, parameters come from (in GPS seconds)
  -t, --CoherenceTime      Desired span of data actually integrated to look for signals (in seconds)
  -e, --EphemerisDirectory Directory containing Sun and Earth ephemeris files
  -y, --EphemerisYear      Year of ephemeris files needed, i.e. "03-06"
  -N, --NumCandidatesToKeep  Set the total number of top list candidates to keep
  -M, --MetricMismatch     Maximal allowd mismatch for metric tiling (on a per dimension basis)
  -X, --ExecutableCode     Location of the exectuable to perform search, i.e. ./lalapps_ComputeFStatistic
  -L, --LogFileDirectory   Location where the condor .sub log file will be placed
"""
  print >> sys.stderr, msg

################################################
#This section reads in the command line arguments
shortop = "ho:i:C:D:I:a:d:z:c:f:b:s:m:S:E:T:t:e:y:N:M:X:L:"
longop = [
  "help",
  "OutputLabel=",
  "Iteration=",
  "ConfigFile",
  "DataFileLocation=",
  "IFOs=",
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
  "StartTime=",
  "EndTime=",
  "ParameterTime=",
  "CoherenceTime=",
  "EphemerisDirectory=",
  "EphemerisYear=",
  "NumCandidatesToKeep=",
  "MetricMismatch=",
  "ExectuableCode=",
  "LogFileDirectory=",
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
  elif o in ("-o","--OutputLabel"):
    output_label = str(a)
  elif o in ("-i","--Iteration"):
    iteration = int(a)
  elif o in ("-C","--ConfigFile"):
    config_file = str(a)
  elif o in ("-a","--Alpha"):
    alpha = float(a)
  elif o in ("-d","--Delta"):
    delta = float(a)
  elif o in ("-z","--AlphaBand"):
    alpha_band = float(a)
  elif o in ("-c","--DeltaBand"):
    delta_band = float(a)
  elif o in ("-f","--F0"):
    f0 = float(a)
  elif o in ("-b","--F0Band"):
    f0_band = float(a)
  elif o in ("-s","--F1"):
    f1 = float(a)
  elif o in ("-m","--F1Band"):
    f1_band = float(a)
  elif o in ("--F2"):
    f2 = float(a)
  elif o in ("--F2Band"):
    f2_band = float(a)
  elif o in ("--F3"):
    f3 = float(a)
  elif o in ("--F3Band"):
    f3_band = float(a)
  elif o in ("-D","--DataFileLocation"):
    data_location_file = str(a).split(',')
  elif o in ("-I","--IFOs"):
    IFOs = str(a)
  elif o in ("-S","--StartTime"):
    start_time = int(a)
  elif o in ("-E","--EndTime"):
    end_time = int(a)
  elif o in ("-T","--ParameterTime"):
    time_of_params = int(a)
  elif o in ("-t","--CoherenceTime"):
    coherence_time = int(a)
  elif o in ("-e","--EphemerisDirectory"):
    ephemeris_directory = str(a)
  elif o in ("-y","--EphemerisYear"):
    ephemeris_year = str(a)
  elif o in ("-N","--NumCandidatesToKeep"):
    num_candidates_to_keep = str(a)
  elif o in ("-M","--MetricMismatch"):
    metric_mismatch = float(a)
  elif o in ("-C","--ExecutableCode"):
    code = str(a)
  elif o in ("-L","--LogFileDirectory"):
    log_file_directory = str(a)
  else:
a    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

###########################################
#Create a log file with the input parameters
run_log_file = ''.join([output_label,'_commands_',str(iteration),'.log'])
rlf = open(run_log_file,'w')
rlf.write(' '.join(sys.argv))
rlf.write('\n')
rlf.close()

####################################
#Parses the config file, if any
#The config file is overriden by command line parameters




#############################
#End user input section

#################################################
#Makes a folder for this search iteration and moves to it
base_label = output_label + '_' + str(iteration)

run_directory = "./" + base_label + "_run/"
os.mkdir(run_directory)


#############################################
#Calculate some derived parameters and names
#Also performs some basic sanity checks

total_time = end_time - start_time

if total_time < 0:
  print "No valid times: check your start and end times"
  sys.exit(1)

#Need to double check this calculation
#Performs a coarse estimate


if coherence_time > total_time:
  coherence_time = total_time
  final_iteration = True


###########################################################
#Create a directory to store band passed, concatenated SFTs
final_sft_directory = run_directory + 'final_sft_' + base_label + '/'
os.mkdir(final_sft_directory)

##########################################################
#In this step, we loop through each IFO requested,
#examine the data available, use noise waiting to figure out
#which data to use, find the data on the cluster,
#strip out everything but a band around the source, and concatenate 
#the data into a single file per IFO (to ease I/O on the cluster) 

#Opens the data location file

time_array = numpy.array([])
psd_array = numpy.array([])
for this_IFO in valid_IFOs:
  found_file = False
  if this_IFO in IFOs:
    for file_name in data_location_file:
      if this_IFO in file_name.split('/')[-1]:
        try:
          data_file = open(file_name,'r')
        except:
          print "Unable to open " + file_name
          sys.exit(1)
        found_file = True
    if found_file == False:
        print "I don't know which file to use for data location for " + this_IFO
        print "Try adding " + this_IFO + " somewhere to the file name"
        sys.exit(1)

    #If a file containing the noise in band per sft does not exist,
    #create it
    psd_file_name = ''.join(['./Full_psd_',base_label,'_',str(this_IFO),'.txt'])
    
    if (not os.path.exists(psd_file_name)):
      node_directory = ''
      directory_count = 0
      cat_psd_files = ''
    
      for line in data_file:
        new_node_directory = '/'.join(line.split('/')[0:-1])
        if (new_node_directory != node_directory):
          node_directory = new_node_directory
          freq_avg_command = ''.join(['lalapps_FreqAverager_v2 -i "',
                                      node_directory,'/*',str(this_IFO),'*" -f ',
                                      str(freq - freq_band/2.0 - 0.1),
                                      ' -b ',str(freq_band + 0.2),
                                      ' -o ./psd_',str(this_IFO),'_',
                                      str(directory_count),'.txt'])
          directory_count = directory_count + 1
          print freq_avg_command
          os.system(freq_avg_command)

      cat_psd_files = ''.join(['cat ./psd_',str(this_IFO),'_* > ',psd_file_name])
      os.system(cat_psd_files)
      rm_temp_psd_files = 'rm psd_*_*.txt'
      os.system(rm_temp_psd_files)
      
    data_file.close()
        
    #For each sft, note its start time in an array
    #Also note its psd for the frequency band
 

    try:
      psd_file = open(psd_file_name,'r')
    except:
      print "Unable to open " + psd_file_name
      sys.exit(1)
    
    for line in psd_file:
      sft_start_time = int(line.split(' ')[0])
      if ((sft_start_time > start_time) and
          (sft_start_time + SFT_length < end_time)):
        time_array = numpy.append(time_array,sft_start_time)
        psd_value = float(line.split(' ')[2])
        psd_array = numpy.append(psd_array,psd_value)

    psd_file.close()
    #Determine the best section of data to use based on the number of SFTs
    #within a fixed duration, determined by the number of expected templates
      
############################################################################################
#This section determines the best set of SFTs to use given the noise of each SFT in the band

best_start_time = start_time
best_noise_weight = 1
invert_psd_array = 1./psd_array
first_time_pass = True
for x in time_array:
  if numpy.sum((time_array > x) & (time_array < x + coherence_time)) > 0:
    sum_invert_psd = numpy.sum(invert_psd_array * ((time_array > x) & (time_array < x + coherence_time)))
    noise_weight = (1/sum_invert_psd)**0.5
    if (noise_weight < best_noise_weight) or (first_time_pass):
      best_noise_weight = noise_weight
      best_start_time = x
      first_time_pass = False

#Determine what the frequency we should be band passing around is
parameter_time_difference = best_start time = time_or_params
freq_current = f0 + (f1 * parameter_time_difference) + (f2 * parameter_time_difference**2) + (f3 * parameter_time_difference**3) 

print "Best start time is = " + str(best_start_time)

for this_IFO in ['H1','H2','L1']:
  found_file = False
  if this_IFO in IFOs:
    for file_name in data_location_file:
      if this_IFO in file_name.split('/')[-1]:
        try:
          data_file = open(file_name,'r')
        except IOerror:
          print "Unable to open " + file_name
          sys.exit(1)
        found_file = True
    if found_file == False:
        print "I don't know which file to use for data location for " + this_IFO
        print "Try adding " + this_IFO + " somewhere to the file name"
        sys.exit(1)
        

    #Create a temporary directory containing the band passed data
    band_passed_directory = run_directory + 'band_passed_directory_' + base_label + '_' + this_IFO +'/'
    os.mkdir(band_passed_directory)
    
    #Use ConvertToSFTv2 to band pass the data and place in the just created directory    
    lowest_freq = freq_current - freq_band/2.0 - 1.0
    highest_freq = freq_current + freq_band/2.0 + 1.0
    print "Freq bandpass from: " + str(lowest_freq) + " to " + str(highest_freq)
    for line in data_file:
      sft_start_time = int(line.split('-')[-2])
      if (sft_start_time > best_start_time) and (sft_start_time < (best_start_time + Tcoh)):
        band_pass_command = ('lalapps_ConvertToSFTv2 ' + '-i ' + line.strip() + 
                             ' -I ' + this_IFO + ' -o ' + band_passed_directory + 
                             ' -f ' + str(lowest_freq) + ' -F ' + str(highest_freq))
        os.system(band_pass_command)
    data_file.close()
    
    
    #Merge stripped SFTs into a single file - remove any prior attempts
    cat_file_name = final_sft_directory + this_IFO + '-' + str(best_start_time) + '-' + str(Tcoh)
    print cat_file_name
    try:
      os.remove('cat_file_name')
    except OSError:
      pass

              
    for time_section in range(int((best_start_time / 1000000)),int((best_start_time+Tcoh)/ 1000000)+1):
      cat_file_name = final_sft_directory + this_IFO + '-' + base_name + '.sft'
      cat_command = ('cat ' + band_passed_directory + '*' + str(time_section) + 
                     ('?'*6) + '-' + str(SFT_length) + '.sft >> ' + cat_file_name)
      os.system(cat_command)

    #Removes the stripped SFTs since we now have the fully merged SFT
    remove_temporary = 'rm -rf ' + band_passed_directory
    os.system(remove_temporary)
  
###############################################################################
#This step generates the .dag file which is submitted on the cluster
#It first calculates the various step sizes and how the parameter space is 
#divided into the jobs

dag_file_name = base_label + '.dag'
sub_file_name = base_label + '.sub'
sub_file_name_examine = base_label + '_examine.sub'

output_results_dir = run_directory + 'output_results_' + base_label +'/'
os.mkdir(output_results_dir)

#The following are derived from "Parameter space metric for combined diurnal and
# orbital motion" by Ian Jones, Ben Owen, and David Whitbeck
#which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls5directed&action=view&page=5
#and also from "Monte-Carlo tests of ComputeFStatistic for 4D parameter-spaces" 
#by Reinhard Prix which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls3knownpulsardemod&action=view&page=2

freq_resolution = 2*(3*mismatch/((math.pi**2)*(Tcoh**2)))**0.5
dfreq_resolution = 2*(4*5*9*mismatch/((math.pi**2)*(Tcoh**4)))**0.5

#Determine spacing in right ascension (alpha), declination (delta), 
#using the template grid of the actual code

sky_template_count_file = ''.join(['SkyTemplateCount_',base_label,'.txt'])
find_angular_resolution = ' '.join([code,
                           '--Alpha', str(alpha - alpha_band/2.0),
                           '--Delta', str(delta - delta_band/2.0),
                           '--AlphaBand', str(alpha_band),
                           '--DeltaBand', str(delta_band),
                           '--Freq', str(freq),
                           '--FreqBand', str(0.0),
                           '--f1dot', str(f_dot),
                           '--f1dotBand', str(0.0),
                           '--DataFiles', ''.join([final_sft_directory,'*']),
                           '--ephemDir', ephemeris_directory,
                           '--ephemYear', ephemeris_year,
                           '--refTime',str(time_of_params),
                           '--gridType',str(grid_type),
                           '--metricType',str(metric_type),
                           '--metricMismatch',str(mismatch),
                           '--countTemplates',
                           '>',sky_template_count_file])
                           
os.system(find_angular_resolution)
ar_file = open(sky_template_count_file,'r')
for line in ar_file:
  if line.split(':')[0] == '%% Number of templates':
    sky_grid_templates = float(line.split(':')[1])

angular_resolution = (alpha_band*delta_band/sky_grid_templates)**0.5

total_templates = sky_grid_templates * freq_band * f_dot_band / (freq_resolution * dfreq_resolution)

alpha_steps = int(math.ceil(alpha_band/angular_resolution))
delta_steps = int(math.ceil(delta_band/angular_resolution))

if (alpha_steps*delta_steps < math.floor(max_number_of_jobs)/2):
  alpha_step_size = angular_resolution
  delta_step_size = angular_resolution
  dfreq_steps = int(math.floor(max_number_of_jobs/(alpha_steps*delta_steps)))
  dfreq_step_size = f_dot_band/dfreq_steps
elif (alpha_steps*delta_steps > max_number_of_jobs):
  alpha_steps = int(math.floor(max_number_of_jobs**0.5))
  delta_steps = alpha_steps
  alpha_step_size = alpha_band/alpha_steps
  delta_step_size = delta_band/delta_steps
  dfreq_steps = 1
  dfreq_step_size = f_dot_band
else:
  alpha_step_size = angular_resolution
  delta_step_size = angular_resolution
  dfreq_steps = 1
  dfreq_step_size = f_dot_band

dag_file_handle = open(dag_file_name,'w')

actual_number_of_jobs = alpha_steps*delta_steps*dfreq_steps

dag_file_handle.write('# Dag for base_name = ' + base_label + '\n')
job = 0
for alpha_count in range(0,alpha_steps):
  for delta_count in range(0,delta_steps):
    for dfreq_count in range(0,dfreq_steps):
      arg_list = ' '.join([' argList=" ',
                           '--Alpha', str(alpha - alpha_band/2.0 + alpha_step_size*alpha_count),
                           '--Delta', str(delta - delta_band/2.0 + delta_step_size*delta_count),
                           '--AlphaBand', str(alpha_step_size),
                           '--DeltaBand', str(delta_step_size),
                           '--Freq', str(freq - freq_band/2.0),
                           '--FreqBand', str(freq_band),
                           '--dFreq', str(freq_resolution),
                           '--f1dot', str(f_dot - f_dot_band/2.0 + dfreq_step_size*dfreq_count),
                           '--f1dotBand', str(dfreq_step_size),
                           '--df1dot', str(dfreq_resolution),
                           '--DataFiles', ''.join([final_sft_directory,'*']),
                           '--ephemDir', ephemeris_directory,
                           '--ephemYear', ephemeris_year,
                           '--NumCandidatesToKeep',str(num_candidates_to_keep),
                           '--refTime',str(time_of_params),
                           '--gridType',str(grid_type),
                           '--metricType',str(metric_type),
                           '--metricMismatch',str(mismatch),
                           '--outputFstat', ''.join([output_results_dir,
                                                     base_label,'_result_',str(job)]),
                           '--outputLoudest', ''.join([output_results_dir,
                                                       base_label,'_loudest_',str(job)]),
                           '"'])
                           
      dag_file_handle.write('JOB A' + str(job) + ' ' + sub_file_name +'\n')
      dag_file_handle.write('VARS A' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + arg_list + '\n')
      dag_file_handle.write('\n')
      job += 1

arg_list_examine = ' '.join([' argList=" ',
                             '--temp"'])

##############################################
#Call followup script after dag searches finish
dag_file_handle.write('JOB B0 ' + sub_file_name_examine + '\n')

arg_list_examine = ' '.join([' argList=" ',
                             '-o ./',
                             '-j', str(job),
                             '-b', str(output_label),
                             '-d', str(output_results_dir),
                             '-i', str(iteration),
                             '-t', str(Tcoh),
                             '-D', ','.join(data_location_file),
                             '-I', str(IFOs),
                             '-T', str(best_start_time),
                             '-e', str(ephemeris_directory),
                             '-y', str(ephemeris_year),
                             '-L', ''.join([log_file_directory.rstrip('/') + '/']),
                             '-C', str(code),
                             '-A', str(angular_resolution),
                             '-N', str(total_templates),
                             '"'])
                                           
dag_file_handle.write('VARS B0 ' + 'JobID="0" ' + arg_list_examine + '\n')
dag_file_handle.write('\n') 

dag_file_handle.write('PARENT ')
for z in range(0,job):
  dag_file_handle.write('A' + str(z) + ' ')
dag_file_handle.write('CHILD B0')

dag_file_handle.close()

########################################################
#This step creates the .sub file(s) called by the .dag file

#This sub file is used to perform the actual search
sub_file_handle = open(sub_file_name,'w')
sub_file_handle.write('universe = standard \n')
sub_file_handle.write('executable = ' + code +'\n')
sub_file_handle.write('output = node_' + base_label + '_A.out.$(JobID)\n')
sub_file_handle.write('error = node_' + base_label + '_A.err.$(JobID)\n')
sub_file_handle.write('log = ' + log_file_directory.rstrip('/') + '/' + base_label + '.log\n' )
sub_file_handle.write('\n')
sub_file_handle.write('arguments = $(arglist)\n')
sub_file_handle.write('queue\n')

sub_file_handle.close()

#This sub file handles examination of the output
sub_file_examine_handle = open(sub_file_name_examine,'w')
sub_file_examine_handle.write('universe = vanilla \n')
sub_file_examine_handle.write('executable = ./FstatFollowUpExamine.py \n')
sub_file_examine_handle.write('output = node_examine_' + base_label + '.out.$(JobID)\n')
sub_file_examine_handle.write('error = node_examine_' + base_label + '.out.$(JobID)\n')
sub_file_examine_handle.write('log = ' + log_file_directory.rstrip('/') + '/' + base_label + '.log\n')
sub_file_examine_handle.write('\n')
sub_file_examine_handle.write('arguments = $(arglist)\n')
sub_file_examine_handle.write('queue\n')

sub_file_examine_handle.close()


######################################
#This step submits the job to Condor
