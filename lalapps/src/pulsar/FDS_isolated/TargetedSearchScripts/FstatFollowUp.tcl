#!/usr/bin/env tclsh

#FstatFollowUp.tcl - Uses the coherent ComputeFStatistic_v2 code to follow up on a given parameter space.

foreach {var value} {
	outlier_index unknown
	config_file unknown
	iteration unknown
	ra unknown
	dec unknown
	ra_band unknown
	dec_band unknown
	f0 unknown
	f0_band unknown
	param_time unknown
	fdot unknown
	f1_band unknown
	coherence_time unknown
	veto_segment_file_name unknown
	segment_file_name unknown
	veto_segment_list {}
	segment_list {}

	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}
	
source config_file

proc clear_gps {IFO gps} {
	foreach {start stop ifo} $veto_segment_list {
		if {($gps >= $start) & ($gps < $stop) & ($ifo == $IFO) } {
			return 0
			}
		}

	if {$segment_list == {}} { return 1 }
	foreach {start stop ifo} $segment_list {
		if {($gps >= $start) & ($gps < $stop) & ($ifo == $IFO)} {
			return 1
			}
		}
	return 0
	}

#If you create a segment list, it must have data for all ifos.

if {$segment_file_name != unknown} {
	set FILE [open $segment_file_name "r"]
	while { ![eof $FILE] } {
		gets $FILE line
		if {$line == ""} { continue }
		if {[regexp {^#} $line ]} { continue }
		lappend segment_list [lindex $line 0] [lindex $line 1] [lindex $line 2]
		}
	close $FILE
	}

if {$veto_segment_file_name != unknown} {
	set FILE [open $veto_segment_file_name "r"]
	while { ![eof $FILE] } {
		gets $FILE line
		if {$line == ""} { continue }
		if {[regexp {^#} $line ]} { continue }
		lappend veto_segment_list [lindex $line 0] [lindex $line 1] [lindex $line 2]
		}
	close $FILE
	}

set outlier_dir "$work_dir/idx$outlier_index"

file mkdir $outlier_dir

cd $outlier_dir

set LOG_FILE [open outlier_${outlier_index}_${iteration}.log "w"]

puts $LOG_FILE "$argv"

#############################################

set total_time  [expr {$end_time - $start_time}]

if {$total_time < 0} {
	puts stderr "No valid times: check your start and end times"
	exit -1
	}

if {$coherence_time > $total_time} {
	set coherence_time $total_time
	}

set sft_dir "$outlier_dir/sft_$iteration"
file mkdir $sft_dir

set results_dir "$outlier_dir/results_$iteration"
file mkdir $results_dir

set std_out_err_dir "$outlier_dir/std_out_err"
file mkdir $std_out_err_dir

##########################################################
#In this step, we loop through each IFO requested,
#examine the data available, use noise weighting to figure out
#which data to use, find the data on the cluster,
#strip out everything but a band around the source, and concatenate 
#the data into a single file per IFO (to ease I/O on the cluster) 
set start_freq [expr ($f0 - $f0_band/2.0 - $psd_estimation_wings/2.0)]
set freq_band_pass [expr ($start_freq + $psd_estimation_wings)]

set psd_file_name "$outlier_dir/Full_psd.txt"

if { ![file exists $psd_file_name] } {
	foreach {IFO sft_location} $sft_location_files { 

		set SFT_FILE [open $sft_location "r"]
		set node_directory {}
		set directory_count 0
		
		while { ![eof $SFT_FILE] } {
			gets $SFT_FILE line
			if {$line == ""} { continue }
			set new_node_directory [ file dirname $line ]
			if { $new_node_directory != $node_directory } {
				set node_directory $new_node_directory
				puts "${script_dir}/lalapps_FreqAverager_v2 -i ${node_directory}/*$IFO* -f $start_freq -b $freq_band_pass -o ${outlier_dir}/psd_${IFO}_${directory_count}.txt"
				exec ${script_dir}/lalapps_FreqAverager_v2 -i ${node_directory}/*$IFO* -f $start_freq -b $freq_band_pass -o ${outlier_dir}/psd_${IFO}_${directory_count}.txt 2>@stdout
				incr directory_count
				}
			}
		close SFT_FILE
		}
	exec cat $outlier_dir/psd_*_*.txt >> $psd_file_name
	exec rm psd_${IFO}_*_*.txt
	}

set PSD_FILE [open $psd_file_name "r"]
set L {}
while { ![eof $PSD_FILE] } {	
	gets $PSD_FILE line
	if {$line == ""} {continue}
	if {![clear_gps [lindex $line 1] [lindex $line 0] } { continue }
	lappend L $line
	}
close $PSD_FILE

set L [lsort -real -increasing -index 0 $L]

set cumultive_inverse_psd_sum {0}
set first_entry 0
set first_entry_GPS [lindex $L $first_entry 0]
set N [llength $L]
set cum_inverse_weight 0.0
set best_sum -1.0
set GPS_start -1
set GPS_end -1


for {set i 0} {$i < $N} {incr i} {
	#puts "$i [lindex $L $i]"
	set cum_inverse_weight [expr $cum_inverse_weight + 1.0/[lindex $L $i 2]]
	set GPS [lindex $L $i 0]
	lappend cumultive_inverse_psd_sum $cum_inverse_weight
	
	while {$GPS > $first_entry_GPS + $coherence_time} {
		incr first_entry 
		set first_entry_GPS [lindex $L $first_entry 0]
		}

	set sum [expr $cum_inverse_weight - [lindex $cumultive_inverse_psd_sum $first_entry] ] 
	if {$sum > $best_sum} {
		set best_sum $sum
		set GPS_start $first_entry_GPS
		set GPS_end $GPS
		}	
	}

#Determine what the frequency we should be band passing around is
set parameter_time_difference [expr $GPS_start - $param_time]
set freq_current [expr $f0 + ($fdot * $parameter_time_difference) + ($f2 * $parameter_time_difference * $parameter_time_difference) + ($f3 * $parameter_time_difference * $parameter_time_difference * $parameter_time_difference)]
set fdot_current [expr $fdot + ($f2 *  $parameter_time_difference) + ($f3 * $parameter_time_difference * $parameter_time_difference)]

set fdotdot_current [expr $f2 + ($f3 * $parameter_time_difference)]


set param_time $GPS_start
set f0  $freq_current
set f1  $fdot_current
set f2  $fdotdot_current

set lowest_freq  [expr $f0 - $f0_band/2.0 - $frequency_band_wings]
set highest_freq [expr $f0 + $f0_band/2.0 + $frequency_band_wings] 


foreach {IFO sft_location} $sft_location_files {
	
	set band_passed_directory "$outlier_dir/temp_sft_${iteration}_$IFO/"
	file mkdir $band_passed_directory
	
	set SFT_FILE [open $sft_location "r"]

	while { ![eof $SFT_FILE] } {	
		gets $SFT_FILE line
		if {$line == ""} {continue}
		set sft_gps [lindex [split $line "-"] end-1 ]
		if {![clear_gps $IFO $sft_gps] } { continue }
		if {($sft_gps > $GPS_start) & ($sft_gps < ($GPS_start + $coherence_time))} {
			#puts $line

			exec $script_dir/lalapps_ConvertToSFTv2 -i $line -I $IFO -o $band_passed_directory -f $lowest_freq -F $highest_freq 2>@stdout
			}
		}
	close $SFT_FILE
	
	exec find $band_passed_directory -name \*.sft | xargs cat >> $outlier_dir/sft_${iteration}/${IFO}_${GPS_start}-${coherence_time}.sft

	#exec rm -rf $band_passed_directory

	}

#The following are derived from "Parameter space metric for combined diurnal and
# orbital motion" by Ian Jones, Ben Owen, and David Whitbeck
#which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls5directed&action=view&page=5
#and also from "Monte-Carlo tests of ComputeFStatistic for 4D parameter-spaces" 
#by Reinhard Prix which can be found at: 
#http://www.lsc-group.phys.uwm.edu/cgi-bin/enote.pl?nb=puls3knownpulsardemod&action=view&page=2

set freq_resolution [expr 2 * ( sqrt(3*$mismatch/(pow(acos( -1.0),2)*pow($coherence_time,2)))) ]
set dfreq_resolution [expr 2 * ( sqrt( 4*5*9*$mismatch/(pow(acos( -1.0),2)*pow($coherence_time,4))))]

#Determine spacing in right ascension (ra), declination (dec), 
#using the template grid of the actual code




set freq_step_size  $f0_band/$max_number_of_jobs


set dag_name "$outlier_dir/run_${iteration}.dag"
set DAG_FILE [open $dag_name "w"]

set outlier_condorA_sub "$outlier_dir/run_${iteration}_A.sub"
set outlier_condorB_sub "$outlier_dir/run_${iteration}_B.sub"

for {set job 0} {$job < $max_number_of_jobs} {incr job} {

	puts $DAG_FILE "JOB A$job $outlier_condorA_sub"
	puts $DAG_FILE "VARS A$job JobID=\"$job\" argList=\" --Alpha [expr $ra - $ra_band/2.0 ] --Delta [expr $dec - $dec_band/2.0]  --AlphaBand $ra_band --DeltaBand $dec_band --Freq [expr $f0 - $f0_band/2.0 + $freq_step_size*$job] --FreqBand $freq_step_size --dFreq $freq_resolution --f1dot [expr $f1 - $f1_band/2.0] --f1dotBand $f1_band --df1dot $dfreq_resolution --DataFiles $sft_dir/* --ephemDir $ephem_dir,--ephemYear $ephem_year --NumCandidatesToKeep $num_candidates_to_keep --refTime $param_time --gridType $grid_type --metricType $metric_type --metricMismatch $mismatch --TwoFthreshold $min_Fstat_to_keep --outputFstat $results_dir/result_$job --outputLoudest $results_dir/loudest_$job \""
	puts $DAG_FILE ""

	}

puts $DAG_FILE "JOB B$job $outlier_condorB_sub"
	puts $DAG_FILE "VARS B$job JobID=\"$job\" outlier_index $outlier_index config_file $config_file work_dir \"$work_dir\""
	puts $DAG_FILE ""


puts -nonewline $DAG_FILE  "PARENT"
for {set job 0} {$job < $max_number_of_jobs} {incr job} {
	puts -nonewline $DAG_FILE  " A$job"
	}
puts $DAG_FILE "CHILD B0"

close $DAG_FILE

set outlier_condorA_sub_text {
universe = standard 
executable = $fstat_binary
output = $std_out_err_dir/node_${iteration}_A.out.\$(JobID)
error = $std_out_err_dir/node_${iteration}_A.err.\$(JobID)
log = ${log_file_dir}/${search_name}_${outlier_index}_${iteration}.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
queue
}

set FILE [open $outlier_condorA_sub "w"]
puts $FILE [subst $outlier_condorA_sub_text]
close $FILE

set outlier_condorB_sub_text {
universe = vanilla 
executable = $script_dir/FstatFollowUpExamine.tcl
output = $std_out_err_dir/node_${iteration}_0.out.\$(JobID)
error = $std_out_err_dir/node_${iteration}.err.\$(JobID)
log = ${log_file_dir}/${search_name}_${outlier_index}_${iteration}.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
queue
}

set FILE [open $outlier_condorB_sub "w"]
puts $FILE [subst $outlier_condorB_sub_text]
close $FILE
