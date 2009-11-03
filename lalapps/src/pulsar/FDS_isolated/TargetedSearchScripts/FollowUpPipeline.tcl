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
Usage: FollowUpPipeline.py [var value]

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

file copy -force $config_file "$work_dir/config_file"

file mkdir $work_dir/std_out_err

set DAG_FILE [open $main_dag_file_name "w"]

set TARGETS_FILE [open $parameter_source_file "r"]
gets $TARGETS_FILE header
set idx 0
foreach col $header {
	set column_[string tolower $col] $idx
	incr idx
	}
unset idx

set job 0

while { ![eof $TARGETS_FILE] } {
	gets $TARGETS_FILE line
	if {$line == ""} { continue }
	
	set values [ split $line "\t"]
	
	foreach var {ra dec f0 fdot idx} {
		set [string toupper $var] [lindex $values [set column_$var]]
		}
		
	puts $DAG_FILE "JOB A$job $main_condorA_sub"
	puts $DAG_FILE "VARS A$job JobID=\"$job\" outlier_index=\"$IDX\" argList=\"outlier_index $IDX config_file $config_file iteration 0 ra $RA dec $DEC ra_band $ra_band dec_band  $dec_band f0 $F0 f0_band $f0_band param_time  $param_time fdot $FDOT f1_band  $f1_band coherence_time $coherence_time\""
	puts $DAG_FILE ""

	puts $DAG_FILE "JOB B$job $main_condorB_sub"
	puts $DAG_FILE "VARS B$job JobID=\"$job\" argList=\"outlier_index $IDX work_dir $work_dir \""
	puts $DAG_FILE ""
	
	puts $DAG_FILE "PARENT A$job CHILD B$job"
	puts $DAG_FILE ""
	incr job
	}
	
close $DAG_FILE
close $TARGETS_FILE

set FILE [open $main_condorA_sub "w"]
puts $FILE [subst $main_condorA_sub_text]
close $FILE

set FILE [open $main_condorB_sub "w"]
puts $FILE [subst $main_condorB_sub_text]
close $FILE

puts "Created a dag for $job outliers"