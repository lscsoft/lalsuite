#!/usr/bin/env tclsh

#CreatePSDEstimates.tcl - Script which produces averages of SFTs around outliers
foreach {var value} {
	help 0
	config_file unknown
	batch_size 10
	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	if {$var == "help"} {
		puts {
Usage: CreatePSDEstimates.tcl [var value]

  help               display this message
  config_file         File which has default parameters for the follow-ups, i.e. ./followup_config.ini
}
		exit 0
		}
	set $var $value
	}

source $config_file

set psd_estimate_sub "$work_dir/psd.sub"

file mkdir "$work_dir"
file mkdir "$work_dir/psd"
file mkdir "$work_dir/std_out_err"

set file_count 0
set file_list {}
puts $sft_location_files
foreach {IFO sft_location} $sft_location_files {
	if {$sft_location == ""} {continue}
	set line_count 0
	puts $IFO
	puts $sft_location
	set SFT_FILE [open $sft_location "r"]
	set name "$work_dir/${IFO}_sub_location_${file_count}.txt"
	set SFT_SUB_FILE [open $name "w"]
	lappend file_list $name
	while {![eof $SFT_FILE]} {
		if {$line_count > 999} {
			incr file_count
			close $SFT_SUB_FILE
			set name "$work_dir/${IFO}_sub_location_${file_count}.txt"
			set SFT_SUB_FILE [open $name "w"]
			lappend file_list $name
			set line_count 0
			}
		gets $SFT_FILE line
		puts $SFT_SUB_FILE $line
		incr line_count
		}
	close $SFT_SUB_FILE
	incr file_count
	close $SFT_FILE
	}
	
set DAG_FILE [open "$work_dir/psd.dag" "w"]

set job 0
foreach  {sub_sft_file} $file_list {
	puts $DAG_FILE "JOB A$job $psd_estimate_sub"
	puts $DAG_FILE "VARS A$job JobID=\"$job\" argList=\" --no-save  -f ${script_dir}/fast_averager.R --slave --args $sub_sft_file\""
	puts $DAG_FILE ""
	incr job
	}
	

close $DAG_FILE

set psd_estimate_sub_text {
universe= vanilla
executable = /usr/bin/R
output = $work_dir/std_out_err/psd_A.out.\$(JobID)
error = $work_dir/std_out_err/psd_A.err.\$(JobID)
log = ${log_file_dir}/psd_all.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
requirements = Memory > 500
queue
}

set FILE [open $psd_estimate_sub "w"]
puts $FILE [subst $psd_estimate_sub_text]
close $FILE

