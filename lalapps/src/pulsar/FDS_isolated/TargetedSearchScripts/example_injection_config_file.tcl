

# Configuration Script for the coherent followup pipeline.


foreach {var value} {
	do_injection 1
	injection_gps_time 846885755

	search_name  "pf_test_quiet_nsd_good_sft"
	work_dir  "/home/josephb/PowerFluxFollowUp/$search_name"
	scratch_dir "/local/user/josephb"	
	fstat_binary  "/home/josephb/PowerFluxFollowUp/lalapps_ComputeFStatistic_resamp_orig"

	parameter_source_file  "/home/josephb/PowerFluxFollowUp/primary_matches_quiet_no_NA.csv"

	ra_band  0.05
	dec_band  0.05
	
	f0_band  1e-2
	f1_band  0 
	f2_band  0
	f3_band  0
	
	f2 0
	f3 0
	
	coherence_time  "[expr 14*24*60*60]"
	param_time  846885755
	start_time  815155213
	end_time  875232014
	sft_prefix *
	individual_IFOs  False
	iteration  0
	data_iteration 0
	coherence_time_multiplier_per_iteration  4.0
	parameter_space_multiplier  10
	
	false_alarm_cutoff  0.10
	
	num_candidates_to_keep  1000
	min_Fstat_to_keep  0
	percent_largest_to_average  0.90
	
	grid_type  2
	metric_type  1
	mismatch  0.15
	
	max_number_of_jobs  1
	SFT_length  1800
	
	frequency_band_wings  2.0
	psd_estimation_wings 0.2
	
	
	messy  True
	
	ephem_dir  "/home/josephb/head/opt/lscsoft/lal/share/lal"
	ephem_year  "05-09"
	
	sft_location_files  {
		H1 "/home/josephb/PowerFluxFollowUp/TotalListOnNodes_H1.txt" 
		L1  "/home/josephb/PowerFluxFollowUp/TotalListOnNodes_L1.txt"
		}
	
	log_file_dir "/local/user/josephb"
	false_alarm_results_directory "$work_dir/false_alarm"
	events_directory "$work_dir/events"
	script_dir  "/home/josephb/PowerFluxFollowUp"
	psd_dir "/home/josephb/PowerFluxFollowUp/psd_freq_band"

	main_dag_file_name "$work_dir/${search_name}_${iteration}.dag"
	main_condorA_sub "$work_dir/${search_name}_prepare_${iteration}.sub"
	main_condorB_sub "$work_dir/${search_name}_submit_${iteration}.sub"
	
	veto_segment_file_name unknown
	segment_file_name unknown
	
	
	} {
	global $var
	set $var [subst $value]
	#puts "$var [set $var]"
	}


set main_condorA_sub_text {
universe= vanilla 
executable = ${script_dir}/FstatFollowUp_Inj.tcl
output = $work_dir/std_out_err/node_${search_name}_A.out.\$(JobID)
error = $work_dir/std_out_err/node_${search_name}_A.err.\$(JobID)
log = ${log_file_dir}/${search_name}.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
requirements = Memory > 500
should_transfer_files = NO
queue
}

set main_condorB_sub_text {
universe= scheduler
getenv= true
executable = ${script_dir}/MainSubmit_Inj.tcl
output = $work_dir/std_out_err/node_${search_name}_B.out.\$(JobID)
error = $work_dir/std_out_err/node_${search_name}_B.err.\$(JobID)
log = ${log_file_dir}/${search_name}.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
requirements = Memory > 500
should_transfer_files = NO
queue
}
