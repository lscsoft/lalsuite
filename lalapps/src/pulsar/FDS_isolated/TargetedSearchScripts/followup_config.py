# Configuration Script for the coherent followup pipeline.

alpha_band = 0.1
delta_band = 0.1

f0_band = 1e-2
f1_band = 2e-10
f2_band = 0
f3_band = 0

start_coherence_time = 240*60*60 #240 hours converted to seconds

param_time = 846885755

start_time = 815155213
end_time = 877824097

IFOs = 'H1L1'

code = '/archive/home/josephb/PowerFluxFollowUp/lalapps_ComputeFStatistic_v2'
python_dir = '/archive/home/josephb/PowerFluxFollowUp'

data_location_file = '/archive/home/josephb/followScript/TotalListOnNodes_H1.txt,/archive/home/josephb/followScript/TotalListOnNodes_H2.txt,/archive/home/josephb/followScript/TotalListOnNodes_L1.txt'

ephem_dir = '/archive/home/josephb/head/opt/lscsoft/lal/share/lal'
ephem_year = '05-09'

num_candidates_to_keep = 1000
percent_largest_to_average = 0.90

grid_type = 2
metric_type = 1
mismatch = 0.15

max_number_of_jobs = 100





