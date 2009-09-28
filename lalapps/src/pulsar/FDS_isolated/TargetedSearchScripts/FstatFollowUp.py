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



#################################################
#Makes a folder for this search iteration and moves to it
if Vars['iteration'] == 0:
  os.system('rm -rf ' + Vars['current_directory'].rstrip('/') + '/' + Vars['search_name'] + '*')

if not (os.path.exists(Vars['current_directory'])):
  os.mkdir(Vars['current_directory'])
else:
  os.system('rm -rf ' + Vars['self_dag_file_name'] + '*')

try:
  os.mkdir(Vars['run_directory'])
except:
  print "Warning: Run directory already exists, but continuing with DAG generation anyways"

###########################################
#Create a log file with the input parameters
run_log_file = open(Vars['self_log_file_name'],'w')
run_log_file.write(' '.join(sys.argv))
run_log_file.write('\n')
run_log_file.close()


#############################
#End user input section


#def create_band_passed_sfts():

#def create_soft_linked_sfts():


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

###########################################################
#Create a directory to store band passed, concatenated SFTs
if not (os.path.exists(Vars['final_sft_run_directory'])):
  os.mkdir(Vars['final_sft_run_directory'])
else:
  os.system('rm -rf ' + Vars['final_sft_run_directory'])
  os.mkdir(Vars['final_sft_run_directory'])

###############################################
#Create a directory to store the final results of the search
if not (os.path.exists(Vars['output_run_directory'])):
  os.mkdir(Vars['output_run_directory'])
else:
  os.system('rm -rf ' + Vars['output_run_directory'])
  os.mkdir(Vars['output_run_directory'])

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
      node_directory = ''
      directory_count = 0
      cat_psd_files = ''
    
      for line in data_file:
        new_node_directory = '/'.join(line.split('/')[0:-1])
        if (new_node_directory != node_directory):
          node_directory = new_node_directory
          freq_avg_command = ''.join([Vars['python_dir'].rstrip('/'),'/lalapps_FreqAverager_v2 -i "',
                                      node_directory,'/*',str(this_IFO),'*" -f ',
                                      str(Vars['f0'] - Vars['f0_band']/2.0 - 0.1),
                                      ' -b ',str(Vars['f0_band'] + 0.2),
                                      ' -o ', Vars['current_directory'],'/psd_',str(this_IFO),'_',
                                      str(directory_count),'.txt'])
          directory_count = directory_count + 1
          os.system(freq_avg_command)

      cat_psd_files = ''.join(['cat ',Vars['current_directory'],'/psd_',str(this_IFO),'_* > ',psd_file_name])
      os.system(cat_psd_files)
      rm_temp_psd_files = ''.join(['rm ',Vars['current_directory'],'/psd_*_*.txt'])
      os.system(rm_temp_psd_files)
      
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

#Determine what the frequency we should be band passing around is
parameter_time_difference = Vars['best_start_time'] - Vars['param_time']
freq_current = Vars['f0'] + (Vars['f1'] * parameter_time_difference) + (Vars['f2'] * parameter_time_difference**2) + (Vars['f3'] * parameter_time_difference**3) 

print "Best start time is = " + str(Vars['best_start_time'])

for this_IFO in valid_IFOs:
  found_file = False
  if this_IFO in Vars['IFOs']:
    for file_name in Vars['data_location_file']:
      if this_IFO in file_name.split('/')[-1]:
        try:
          data_file = open(file_name,'r')
        except IOerror:
          print "Unable to open " + file_name
          sys.exit(-1)
        found_file = True
    if found_file == False:
        print "I don't know which file to use for data location for " + this_IFO
        print "Try adding " + this_IFO + " somewhere to the file name"
        sys.exit(-1)
        

    #Create a temporary directory containing the band passed data
    band_passed_directory = Vars['run_directory'] + 'band_passed_directory_' + Vars['run_name'] + '_' + this_IFO +'/'
    if not (os.path.exists(band_passed_directory)):
      os.mkdir(band_passed_directory)
    else:
      os.system('rm -rf ' + band_passed_directory)
      os.mkdir(band_passed_directory)
    
    #Use ConvertToSFTv2 to band pass the data and place in the just created directory    
    lowest_freq = freq_current - Vars['f0_band']/2.0 - Vars['frequency_band_wings']
    highest_freq = freq_current + Vars['f0_band']/2.0 + Vars['frequency_band_wings']
    print "Freq bandpass from: " + str(lowest_freq) + " to " + str(highest_freq)
    for line in data_file:
      sft_start_time = int(line.split('-')[-2])
      if (sft_start_time > Vars['best_start_time']) and (sft_start_time < (Vars['best_start_time'] + Vars['coherence_time'])):
        band_pass_command = (Vars['python_dir'].rstrip('/') + '/lalapps_ConvertToSFTv2 ' + '-i ' + line.strip() + 
                             ' -I ' + this_IFO + ' -o ' + band_passed_directory + 
                             ' -f ' + str(lowest_freq) + ' -F ' + str(highest_freq))
        os.system(band_pass_command)
    data_file.close()
    
    
    #Merge stripped SFTs into a single file - remove any prior attempts
    cat_file_name = Vars['final_sft_run_directory'] + this_IFO + '-' + str(Vars['best_start_time']) + '-' + str(Vars['coherence_time'])
    print cat_file_name
    try:
      os.remove('cat_file_name')
    except OSError:
      pass

              
    for time_section in range(int((Vars['best_start_time'] / 1000000)),int((Vars['best_start_time']+Vars['coherence_time'])/ 1000000)+1):
      cat_file_name = Vars['final_sft_run_directory'] + this_IFO + '-' + Vars['run_name'] + '.sft'
      cat_command = ('cat ' + band_passed_directory + '*' + str(time_section) + 
                     ('?'*6) + '-' + str(Vars['SFT_length']) + '.sft >> ' + cat_file_name)
      os.system(cat_command)

    #Removes the stripped SFTs since we now have the fully merged SFT
    remove_temporary = 'rm -rf ' + band_passed_directory
    os.system(remove_temporary)
  
###############################################################################
#This step generates the .dag file which is submitted on the cluster
#It first calculates the various step sizes and how the parameter space is 
#divided into the jobs

if not (os.path.exists(Vars['output_run_directory'])):
  os.mkdir(Vars['output_run_directory'])

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

#################

Vars['freq_step_size'] = Vars['f0_band']/Vars['max_number_of_jobs']


##################
dag_file_handle = open(Vars['self_dag_file_name'],'w')


dag_file_handle.write('# Dag for base_name = ' + Vars['run_name'] + '\n')
job = 0
for freq_count in range(0,Vars['max_number_of_jobs']):
  arg_list = ' '.join([' argList=" ',
                       '--Alpha', str(Vars['alpha'] - Vars['alpha_band']/2.0),
                       '--Delta', str(Vars['delta'] - Vars['delta_band']/2.0),
                       '--AlphaBand', str(Vars['alpha_band']),
                       '--DeltaBand', str(Vars['delta_band']),
                       '--Freq', str(Vars['f0'] - Vars['f0_band']/2.0 + Vars['freq_step_size']*freq_count),
                       '--FreqBand', str(Vars['freq_step_size']),
                       '--dFreq', str(Vars['freq_resolution']),
                       '--f1dot', str(Vars['f1'] - Vars['f1_band']/2.0),
                       '--f1dotBand', str(Vars['f1_band']),
                       '--df1dot', str(Vars['dfreq_resolution']),
                       '--DataFiles', ''.join([Vars['final_sft_run_directory'].rstrip('/'),'/*']),
                       '--ephemDir', Vars['ephem_dir'],
                       '--ephemYear', Vars['ephem_year'],
                       '--NumCandidatesToKeep',str(Vars['num_candidates_to_keep']),
                       '--refTime',str(Vars['param_time']),
                       '--gridType',str(Vars['grid_type']),
                       '--metricType',str(Vars['metric_type']),
                       '--metricMismatch',str(Vars['mismatch']),
                       '--TwoFthreshold',str(Vars['min_Fstat_to_keep']),
                       '--outputFstat', ''.join([Vars['output_run_directory'],
                                                 '/',Vars['results_base_prefix'],str(job)]),
                       '--outputLoudest', ''.join([Vars['output_run_directory'],'/',
                                                   Vars['loudest_base_prefix'],str(job)]),
                       '"'])
                           
  dag_file_handle.write('JOB A' + str(job) + ' ' + os.path.basename(Vars['self_sub_file_name']) +'\n')
  dag_file_handle.write('VARS A' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + arg_list + '\n')
  dag_file_handle.write('\n')
  job += 1

job = 0
for this_IFO in valid_IFOs:
  if this_IFO in Vars['IFOs']:
    for freq_count in range(0,Vars['max_number_of_jobs']):
      arg_list = ' '.join([' argList=" ',
                           '--Alpha', str(Vars['alpha'] - Vars['alpha_band']/2.0),
                           '--Delta', str(Vars['delta'] - Vars['delta_band']/2.0),
                           '--AlphaBand', str(Vars['alpha_band']),
                           '--DeltaBand', str(Vars['delta_band']),
                           '--Freq', str(Vars['f0'] - Vars['f0_band']/2.0 + Vars['freq_step_size']*freq_count),
                           '--FreqBand', str(Vars['freq_step_size']),
                           '--dFreq', str(Vars['freq_resolution']),
                           '--f1dot', str(Vars['f1'] - Vars['f1_band']/2.0),
                           '--f1dotBand', str(Vars['f1_band']),
                           '--df1dot', str(Vars['dfreq_resolution']),
                           '--DataFiles', ''.join([Vars['final_sft_run_directory'].rstrip('/'),'/',this_IFO,'*']),
                           '--ephemDir', Vars['ephem_dir'],
                           '--ephemYear', Vars['ephem_year'],
                           '--NumCandidatesToKeep',str(Vars['num_candidates_to_keep']),
                           '--refTime',str(Vars['param_time']),
                           '--gridType',str(Vars['grid_type']),
                           '--metricType',str(Vars['metric_type']),
                           '--metricMismatch',str(Vars['mismatch']),
                           '--TwoFthreshold',str(Vars['min_Fstat_to_keep']),
                           '--outputFstat', ''.join([Vars['output_run_directory'],
                                                     '/',this_IFO,'_',Vars['results_base_prefix'],str(job)]),
                           '--outputLoudest', ''.join([Vars['output_run_directory'],'/', this_IFO,'_',
                                                       Vars['loudest_base_prefix'],str(job)]),
                           '"'])
      
      dag_file_handle.write('JOB D' + str(job) + ' ' + os.path.basename(Vars['self_sub_file_name']) +'\n')
      dag_file_handle.write('VARS D' + str(job) + ' ' + 'JobID="' + str(job) + '" ' + arg_list + '\n')
      dag_file_handle.write('\n')
      job += 1

##############################################
#Call followup script after dag searches finish
dag_file_handle.write('JOB B0 ' + os.path.basename(Vars['self_sub_examine_file_name']) + '\n')

arg_list_examine = ' '.join([' argList=" ',
                             '-o', str(Vars['base_name']),
                             '-C', str(Vars['config_file']),
                             '-j', str(job),
                             '-i', str(Vars['iteration']),
                             '-t', str(Vars['coherence_time']),
                             '-T', str(Vars['best_start_time']),
                             '-A', str(Vars['angular_resolution']),
                             '-N', str(Vars['total_templates']),
                             '"'])
                                           
dag_file_handle.write('VARS B0 ' + 'JobID="0" ' + arg_list_examine + '\n')
dag_file_handle.write('\n')



Vars['dag_file_for_next_run'] = Vars['current_directory'] + '/' + Vars['base_name'] + '_' + str(Vars['iteration']+1) + '.dag'

arg_list_resubmit = ' '.join([' argList =" ',
                              '-l',Vars['submit_file_log_dir'],
                              '-f',Vars['dag_file_for_next_run'],
                              '-d',Vars['current_directory'],
                              '"'])

dag_file_handle.write('JOB C0 ' + os.path.basename(Vars['self_resubmit_file_name']) + '\n')
dag_file_handle.write('VARS C0 ' + 'JobID="0" ' + arg_list_resubmit + '\n')
dag_file_handle.write('\n')

dag_file_handle.write('PARENT ')
for z in range(0,Vars['max_number_of_jobs']):
  dag_file_handle.write('A' + str(z) + ' ')
dag_file_handle.write('CHILD B0')
dag_file_handle.write('\n')

dag_file_handle.write('PARENT B0 CHILD C0')

dag_file_handle.close()

########################################################
#This step creates the .sub file(s) called by the .dag file

#This sub file is used to perform the actual search
sub_file_handle = open(Vars['self_sub_file_name'],'w')
sub_file_handle.write('universe = standard \n')
sub_file_handle.write('executable = ' + Vars['code'] +'\n')
sub_file_handle.write('output = node_' + Vars['run_name'] + '_A.out.$(JobID)\n')
sub_file_handle.write('error = node_' + Vars['run_name'] + '_A.err.$(JobID)\n')
sub_file_handle.write('log = ' + Vars['submit_file_log_dir'].rstrip('/') + '/' + Vars['run_name'] + '.log\n' )
sub_file_handle.write('\n')
sub_file_handle.write('arguments = $(arglist)\n')
sub_file_handle.write('queue\n')

sub_file_handle.close()

#This sub file handles examination of the output
sub_file_examine_handle = open(Vars['self_sub_examine_file_name'],'w')
sub_file_examine_handle.write('universe = vanilla \n')
sub_file_examine_handle.write('executable = ' + Vars['python_dir'] + '/FstatFollowUpExamine.py \n')
sub_file_examine_handle.write('output = node_examine_' + Vars['run_name'] + '.out.$(JobID)\n')
sub_file_examine_handle.write('error = node_examine_' + Vars['run_name'] + '.err.$(JobID)\n')
sub_file_examine_handle.write('log = ' + Vars['submit_file_log_dir'].rstrip('/') + '/' + Vars['run_name'] + '.log\n')
sub_file_examine_handle.write('\n')
sub_file_examine_handle.write('arguments = $(arglist)\n')
sub_file_examine_handle.write('queue\n')

sub_file_examine_handle.close()

#This sub file handles the starting of the next run in the scheduler universe
#This subfile handles the scripts run in the scheduler universe so more dags can be run


sub_file_resubmit = open(Vars['self_resubmit_file_name'],'w')
sub_file_resubmit.write('universe= scheduler \n')
sub_file_resubmit.write('executable = ' + Vars['python_dir'].rstrip('/') + '/' + 'FstatFollowUpSubmit.py\n')
sub_file_resubmit.write('output = node_submit_' + Vars['run_name'] + '_C.out.$(JobID)\n')
sub_file_resubmit.write('error = node_submit_' + Vars['run_name'] + '_C.err.$(JobID)\n')
sub_file_resubmit.write('log = ' + Vars['submit_file_log_dir'].rstrip('/') + '/' + Vars['run_name'] + '.log\n')
sub_file_resubmit.write('\n')
sub_file_resubmit.write('arguments = $(arglist)\n')
sub_file_resubmit.write('queue\n')

sub_file_resubmit.close()



######################################
#This step submits the job to Condor

#os.chdir(Vars['current_directory'])

#condor_command = 'condor_submit_dag -outfile_dir ' + log_file_directory + ' ' + os.path.basename(Vars['self_dag_file_name'])
#os.system(condor_command)

sys.exit(0)
