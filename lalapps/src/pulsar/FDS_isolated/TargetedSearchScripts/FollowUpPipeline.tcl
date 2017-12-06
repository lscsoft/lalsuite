#!/usr/bin/env tclsh

#FollowUpPipeline.tcl - Main script which calls sub-scripts to perform a coherent follow up to continuous gravitational wave sources
foreach {var value} {
	help 0
	config_file unknown
	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	if {$var == "help"} {
		puts {
Usage: FollowUpPipeline.tcl [var value]

  help               display this message
  config_file         File which has default parameters for the follow-ups, i.e. ./followup_config.ini
}
		exit 0
		}
	set $var $value
	}

source $config_file

if { [file exists $work_dir] } { 
	puts "Working directory \"$work_dir\" already exists"
	}
	
file mkdir $work_dir

file mkdir $false_alarm_results_directory

file mkdir $events_directory

file copy -force $config_file "$work_dir/config_file_${iteration}"

file mkdir $work_dir/std_out_err

set DAG_FILE [open $main_dag_file_name "w"]

array set var_tags {dec dec ra ra delta dec alpha ra f0 f0 fdot f1 f1 f1 f2 f2 f3 f3 f4 f4 f5 f5 line_id line_id param_time param_time}

set TARGETS_FILE [open $parameter_source_file "r"]
gets $TARGETS_FILE header
set idx 0
foreach col $header {
	if { [
		catch {
			set column_$var_tags([string tolower $col]) $idx
			} ] } {
			set column_[string tolower $col] $idx
			}
	incr idx
	}
unset idx

set job 0

while { ![eof $TARGETS_FILE] } {
	gets $TARGETS_FILE line
	if {$line == ""} { continue }
	
	set values [ split $line "\t"]
	
	foreach var {ra dec f0 f1 line_id min_gps max_gps param_time} {
		if {[info exists column_$var]} {
			set [string toupper $var] [lindex $values [set column_$var]]
			} else {
			set [string toupper $var] NA
			}
		}
		
	if {$do_injection} {
		foreach var {line_id ra_orig dec_orig psi_orig phi_orig iota_orig f0_orig spindown_orig aplus across param_time} {
			if {[info exists column_$var]} {
				set [string toupper $var] [string trim [lindex $values [set column_$var]] {"}]
					#"
				} else {
				set [string toupper $var] NA
				}
			}
		puts $DAG_FILE "JOB A$job $main_condorA_sub"
		puts $DAG_FILE "VARS A$job JobID=\"$job\" outlier_index=\"$LINE_ID\" argList=\"outlier_index $LINE_ID config_file $config_file iteration $iteration ra $RA dec $DEC ra_band $ra_band dec_band  $dec_band f0 $F0 f0_band $f0_band param_time  $PARAM_TIME fdot $F1 f1_band  $f1_band coherence_time $coherence_time line_id $LINE_ID ra_orig $RA_ORIG dec_orig $DEC_ORIG psi_orig $PSI_ORIG phi_orig $PHI_ORIG iota_orig $IOTA_ORIG f0_orig $F0_ORIG spindown_orig $SPINDOWN_ORIG aplus $APLUS across $ACROSS min_gps $MIN_GPS max_gps $MAX_GPS \""
		puts $DAG_FILE ""

		incr job
		}	
		
	if {!$do_injection} {
		foreach var {line_id} {
			set [string toupper $var] [string trim [lindex $values [set column_$var]] {"}]
				#"
			}
		
		puts $DAG_FILE "JOB A$job $main_condorA_sub"
		puts $DAG_FILE "VARS A$job JobID=\"$job\" outlier_index=\"$LINE_ID\" argList=\"outlier_index $LINE_ID config_file $config_file iteration $iteration ra $RA dec $DEC ra_band $ra_band dec_band  $dec_band f0 $F0 f0_band $f0_band param_time  $PARAM_TIME fdot $F1 f1_band  $f1_band coherence_time $coherence_time min_gps $MIN_GPS max_gps $MAX_GPS \""
		puts $DAG_FILE ""
		
		incr job
		}
	}


puts $DAG_FILE "JOB B0 $main_condorB_sub"
puts $DAG_FILE "VARS B0 JobID=\"0\" argList=\"work_dir $work_dir iteration $iteration \""
puts $DAG_FILE ""

puts -nonewline $DAG_FILE "PARENT "
for {set j 0} {$j < $job} { incr j} {
	puts -nonewline $DAG_FILE "A$j "
	}
puts $DAG_FILE "CHILD B0"
puts $DAG_FILE ""
	
close $DAG_FILE
close $TARGETS_FILE

set FILE [open $main_condorA_sub "w"]
puts $FILE [subst $main_condorA_sub_text]
close $FILE

set FILE [open $main_condorB_sub "w"]
puts $FILE [subst $main_condorB_sub_text]
close $FILE

puts "Created a dag for $job outliers in $work_dir"