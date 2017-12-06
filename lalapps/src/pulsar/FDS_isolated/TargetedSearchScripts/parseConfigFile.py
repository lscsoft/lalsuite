import os
import sys


def parse_config_file(config_file_name,Vars):

    ####################################
    #Parses the config file, if any
    #The config file is overriden by command line parameters

    try:
        config_file = os.path.abspath(config_file_name)
        sys.path.append(os.path.dirname(config_file_name))
        config = __import__(os.path.basename(config_file_name).rstrip('.py'));
    except:
        print "Could not find or load config file: " + config_file_name
        sys.exit(1)

    #Failed flag
    failed = False

    ### Self reference to config file
    Vars['config_file'] = os.path.abspath(config_file_name)

    ### Search name
    try:
        Vars['search_name'] = config.search_name
    except:
        print "search_name cannot be read"
        failed = True

    ### Work directory
    try:
        Vars['work_dir'] = os.path.abspath(config.work_dir)
        if not os.path.exists(os.path.dirname(Vars['work_dir'].rstrip('/'))):
            print "The preceding directories for the work directory cannot be found: " + Vars['work_dir']
            failed = True        
    except:
        print "work_dir cannot be read"
        failed = True

    ### Submit file log file directory
    try:
        Vars['submit_file_log_dir'] = os.path.abspath(config.submit_file_log_dir)

        if not os.path.exists(Vars['submit_file_log_dir']):
            print "The log directory to be used by the submit files cannot be found: " + Vars['submit_file_log_dir']
            failed = True
    except:
        print "submit_file_log_dir cannot be read"
        failed = True

    ### Code to be used for search
    try:
        Vars['code'] = os.path.abspath(config.code)
        if not os.path.exists(Vars['code']):
            print "The code listed in the config file cannot be found: " + Vars['code']
            failed = True
    except:
        print "code (the configuration variable) cannot be read"
        failed = True

    ### Directory containing the python scripts
    try:
        Vars['python_dir'] = os.path.abspath(config.python_dir)
        if not os.path.exists(Vars['code']):
            print "The python directory listed in the config file cannot be found: " + Vars['python_dir']
            failed = True
    except:
        print "python_dir cannot be read"
        failed = True

    ### Locations of the SFTs on the cluster or system
    try:
        Vars['data_location_file'] = config.data_location_file.split(',')
        for file_location in Vars['data_location_file']:
            if not os.path.exists(file_location):
                print "A data file listed in the config file cannot be found: " + file_location
                failed = True
    except:
        print "data_location_file cannot be read"
        failed = True

    ### Location of the file containing the parameters (outliers, sources, etc)
    try:
        Vars['parameter_source_file'] = config.parameter_source_file
        if not os.path.exists(Vars['parameter_source_file']):
            print "The file containing the parameters of the search cannot be found: " + Vars['parameter_source_file']
            failed = True
    except:
        print "parameter_source_file cannot be read"
        failed = True
        
    ### Earth and Sun Ephemeris directory
    try:
        Vars['ephem_dir'] = config.ephem_dir
    except:
        print "ephem_dir cannot be read"
        failed = True
      
    ### Earth and Sun Ephemeris year to use
    try:
        Vars['ephem_year'] = config.ephem_year
    except:
        print "ephem_year cannot be read"
        failed = True

    ### Band in Right Ascension (radians)
    try:
        if not ("alpha_band" in Vars):
            Vars['alpha_band'] = config.alpha_band
    except:
        print "alpha_band cannot be read"
        failed = True
  
    ### Band in Declination (radians)
    try:
        if not ("delta_band" in Vars):
            Vars['delta_band'] = config.delta_band
    except:
        print "delta_band cannot be read"
        failed = True
    
    ### Band in Frequency (Hz)
    try:
        if not ("f0_band" in Vars):
            Vars['f0_band'] = config.f0_band
    except:
        print "f0_band cannot be read"
        failed = True
      
    ### Band in 1st Spindown (f dot) (Hz/s)
    try:
        if not ("f1_band" in Vars):
            Vars['f1_band'] = config.f1_band
    except:
        print "f1_band cannot be read"
        failed = True
      
    ### Band in 2nd Spindown (f dot dot) (Hz/s^2)
    try:
        if not ("f2_band" in Vars):
            Vars['f2_band'] = config.f2_band
    except:
        print "f2_band cannot be read"
        failed = True

    ### Band in 3rd Spindown (f dot dot dot) (Hz/s^3)
    try:
        if not ("f3_band" in Vars):
            Vars['f3_band'] = config.f3_band
    except:
        print "f3_band cannot be read"
        failed = True

    ### Initial coherence time
    try:
        if not ("start_coherence_time" in Vars):
            Vars['coherence_time'] = config.start_coherence_time
    except:
        print "start_coherence_time cannot be read"
        failed = True

    ### First GPS time to be considered for data selection
    try:
        Vars['start_time'] = config.start_time
    except:
        print "start_time cannot be read"
        failed = True
    
    ### Last GPS time to be considered for data selection
    try:
        Vars['end_time'] = config.end_time
    except:
        print "end_time cannot be read"
        failed = True
      
    ### Time from which the frequency and spindown values come
    try:
        if not ("param_time" in Vars):
            Vars['param_time'] = config.param_time
    except:
        print "param_time cannot be read"
        failed = True

    ### IFOs to be used
    try:
        Vars['IFOs'] = config.IFOs
    except:
        print "IFOs cannot be read"
        failed = True 
    
    ### Maximum 1-D mismatch between template and possible signal
    try:
        Vars['mismatch'] = config.mismatch
    except:
        print "mismatch cannot be read"
        failed = True

    ### Number of candidates to keep per job 
    try:
        Vars['num_candidates_to_keep'] = config.num_candidates_to_keep
    except:
        print "num_candidates_to_keep cannot be read"
        failed = True

    ### Minimum Fstat values to keep
    try:
        Vars['min_Fstat_to_keep'] = config.min_Fstat_to_keep
    except:
        print "min_Fstat_to_keep cannot be read"
        failed = True

    ### Percentage of the largest 2F value to consider when averaging results to find the best spot to followup
    try:
        Vars['percent_largest_to_average'] = config.percent_largest_to_average
    except:
        print "percent_largest_to_average cannot be read"
        failed = True

    ### Grid type used by ComputeFStatistic
    try:
        Vars['grid_type'] = config.grid_type
    except:
        print "grid_type cannot be read"
        failed = True

    ### Metric type used by ComputeFStatistic
    try:
        Vars['metric_type'] = config.metric_type
    except:
        print "metric_type cannot be read"
        failed = True    

    ### Max number of jobs to run per event per step
    try:
        Vars['max_number_of_jobs'] = config.max_number_of_jobs
    except:
        print "max_number_of_jobs cannot be read"
        failed = True

    ### Length of SFTs to be used in seconds
    try:
        Vars['SFT_length'] = config.SFT_length
    except:
        print "SFT_length cannot be read"
        failed = True

    ### Frequency band wings to be added to the frequency band narrow passed both above and below
    try:
        Vars['frequency_band_wings'] = config.frequency_band_wings
    except:
        print "frequency_band_wings cannot be read"
        failed = True

    ### Iteration to start with
    try:
        if not ('iteration' in Vars):
            Vars['iteration'] = config.iteration
    except:
        print "iteration cannot be read"
        failed = True

    ### Multiplication factor by which to increase the coherence time each run if a significant event is found
    try:
        Vars['coherence_time_multiplier_per_iteration'] = config.coherence_time_multiplier_per_iteration
    except:
        print "coherence_time_multiplier_per_iteration cannot be read"
        failed = True

    ### Parameter space multiplier - multiplier times resolution of the previous iteration to determine new parameter space size to search

    try:
        Vars['parameter_space_multiplier'] = config.parameter_space_multiplier
    except:
        print "parameter_space_multiplier cannot be read"
        failed = True
        
### Parameter space multiplier - multiplier times resolution of the previous iteration to determine new parameter space size to search

    try:
        Vars['false_alarm_cutoff'] = config.false_alarm_cutoff
    except:
        print "false_alarm_cutoff cannot be read"
        failed = True

    ### Determines whether to leave intermediate files for debugging purposes or to clean up after itself
    try:
        Vars['messy'] = config.messy
    except:
        print "messy cannot be read"
        failed = True

    ### Create only soft links to the sfts on the nodes (True) or band pass and create concatenated files on the local storage drive (False)
    try:
        Vars['soft_link_sfts_only'] = config.soft_link_sfts_only
    except:
        print "soft_link_sfts cannot be read"
        failed = True

    if failed == True:
        sys.exit(1)
    else:
        return Vars


#### End parse_config_File

def generate_directory_names(identifier,Vars):

    Vars['events_directory'] = ''.join([Vars['work_dir'].rstrip('/'),'/events_',Vars['search_name'],'/'])

    Vars['false_alarm_results_directory'] = ''.join([Vars['work_dir'].rstrip('/'),'/false_alarm_results_',Vars['search_name'],'/'])

    Vars['current_directory'] = ''.join([Vars['work_dir'].rstrip('/'),'/',identifier])

    Vars['run_name'] = ''.join([identifier,'_',str(Vars['iteration'])])

    Vars['run_directory'] = ''.join([Vars['current_directory'].rstrip('/'),'/',Vars['run_name'],'_run'])

    Vars['final_sft_run_directory'] = ''.join([Vars['run_directory'].rstrip('/'),'/','final_sft_',Vars['run_name'],'/'])

    Vars['output_run_directory'] = ''.join([Vars['run_directory'].rstrip('/'),'/','output_',Vars['run_name'],'/'])

    Vars['false_alarm_file'] = ''.join([Vars['false_alarm_results_directory'].rstrip('/'),'/',identifier,'_false_alarm.txt'])

    return Vars

#### End generate_directory_names

def generate_file_names(identifier,Vars):
    Vars['self_log_file_name'] = ''.join([Vars['current_directory'].rstrip('/'),'/',identifier,'.log'])

    Vars['self_dag_file_name'] = ''.join([Vars['current_directory'].rstrip('/'),'/',identifier,'.dag'])

    Vars['self_sub_file_name'] = ''.join([Vars['current_directory'].rstrip('/'),'/',identifier,'.sub'])
    
    Vars['self_resubmit_file_name'] = ''.join([Vars['current_directory'].rstrip('/'),'/',identifier,'_resubmit.sub'])

    Vars['self_sub_log_file_name'] = ''.join([Vars['submit_file_log_dir'].rstrip('/'),'/',identifier,'.log'])

    Vars['self_sub_examine_file_name'] = ''.join([Vars['current_directory'].rstrip('/'),'/',identifier,'_examine.sub'])

    Vars['loudest_result_file'] = ''.join([Vars['current_directory'].rstrip('/'),'/loudest_results_',identifier,'.txt'])

    Vars['loudest_base_prefix'] = ''.join([Vars['run_name'],'_loudest_'])

    Vars['results_base_prefix'] = ''.join([Vars['run_name'],'_result_'])

    Vars['event_result_file'] = ''.join([Vars['events_directory'].rstrip('/'),'/',Vars['run_name'],'_result.txt'])

                                        
    
    return Vars

#### End generate_file_names
