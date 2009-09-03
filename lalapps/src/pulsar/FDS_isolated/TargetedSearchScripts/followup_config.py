# Configuration Script for the coherent followup pipeline.

##########
#Identifier for the search
search_name = 'pf_S5_outliers'

#######
#Directories for various files

work_dir = '/archive/home/josephb/PowerFluxFollowUp/'

submit_file_log_dir = '/usr1/josephb'

code = '/archive/home/josephb/PowerFluxFollowUp/lalapps_ComputeFStatistic_v2'
python_dir = '/archive/home/josephb/PowerFluxFollowUp'

data_location_file = '/archive/home/josephb/followScript/TotalListOnNodes_H1,/archive/home/josephb/followScript/TotalListOnNodes_H2,/archive/home/josephb/followScript/TotalListOnNodes_L1'

######
#Outlier or source parameter file
parameter_source_file = '/archive/home/josephb/PowerFluxFollowUp/short_list.csv'

#######
#Ephemeris information

ephem_dir = '/archive/home/josephb/head/opt/lscsoft/lal/share/lal'

ephem_year = '05-09'

########
#Parameter space information

alpha_band = 0.1
delta_band = 0.1

f0_band = 1e-2
f1_band = 2e-10
f2_band = 0
f3_band = 0

########
#IFOs, Time stamps and duration

start_coherence_time = 250*60*60 #240 hours converted to seconds

param_time = 846885755

start_time = 815155213
end_time = 877824097

IFOs = 'H1L1'

###########
#Other search options

iteration = 0
coherence_time_multiplier_per_iteration = 4.0
parameter_space_multiplier = 10

false_alarm_cutoff = 0.01

num_candidates_to_keep = 1000
min_Fstat_to_keep = 0
percent_largest_to_average = 0.90

grid_type = 2
metric_type = 1
mismatch = 0.15

max_number_of_jobs = 40
SFT_length = 1800

frequency_band_wings = 1.0 # in Hz, applied both above and below

soft_link_sfts_only = True

###########
#Clean up options - set true to leave intermediate files around
#set false to delete intermediate data products

messy = True


