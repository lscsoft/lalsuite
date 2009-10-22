#!/usr/bin/env tclsh

#CreatePSDEstimates.tcl - Script which produces averages of SFTs around outliers
foreach {var value} {
	help 0
	config_file unknown
	max_psd_jobs 50
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

set node_directory_list {}


foreach {IFO sft_location} $sft_location_files { 
	set SFT_FILE [open $sft_location "r"]
	set L {}
	while { ![eof $SFT_FILE] } {
		gets $SFT_FILE line
		if {$line == ""} { continue }
		lappend L [ file dirname $line ]
		}
	close $SFT_FILE
	lappend node_directory_list $IFO [lsort -unique -dictionary $L]
	}

set TARGETS_FILE [open $parameter_source_file "r"]
gets $TARGETS_FILE header
set idx 0
foreach col $header {
	set column_[string tolower $col] $idx
	incr idx
	}
unset idx

file mkdir "$work_dir"
file mkdir "$work_dir/psd"
file mkdir "$work_dir/std_out_err"

set DAG_FILE [open "$work_dir/psd.dag" "w"]

while { ![eof $TARGETS_FILE] } {
	gets $TARGETS_FILE line
	if {$line == ""} { continue }
	
	set values [ split $line "\t"]
	
	foreach var { f0 idx} {
		set [string toupper $var] [lindex $values [set column_$var]]
		}

	set start_freq [expr ($F0 - $f0_band/2.0 - $psd_estimation_wings/2.0)]
	set freq_band_pass [expr ($start_freq + $psd_estimation_wings)]



	set job 0
	foreach  {IFO directory_list} $node_directory_list {
		foreach {directory} $directory_list {
			puts $DAG_FILE "JOB A${IDX}_$job $psd_estimate_sub"
			puts $DAG_FILE "VARS A${IDX}_$job JobID=\"${IDX}_$job\" outlier_index=\"$IDX\" argList=\" -i ${directory}/*${IFO}*.sft -f $start_freq -b $freq_band_pass -o ${work_dir}/psd/psd_${IDX}_${job}_${IFO}.txt\""
			puts $DAG_FILE ""
	
			incr job
			}
		}
	}

set psd_estimate_sub_text {
universe= standard
executable = ${script_dir}/lalapps_FreqAverager_v2
output = $work_dir/std_out_err/psd_A.out.\$(JobID)
error = $work_dir/std_out_err/psd_A.err.\$(JobID)
log = ${log_file_dir}/${search_name}_psd.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
requirements = Memory > 500
queue
}

set FILE [open $psd_estimate_sub "w"]
puts $FILE [subst $psd_estimate_sub_text]
close $FILE

puts "PSD estimation dag and sub files created"