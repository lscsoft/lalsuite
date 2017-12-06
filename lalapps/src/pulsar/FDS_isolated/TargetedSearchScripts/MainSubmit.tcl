#!/usr/bin/env tclsh

#MainSubmit.tcl - Submits the new dag to be started in the scheduler universe

puts "[info hostname]"

foreach {var value} {
	help 0
	work_dir unknown
	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}

after 30000

cd "$work_dir/"

source config_file_${iteration}

set output_dag_name "${work_dir}/full_${search_name}_${iteration}.dag"
set OUTPUT_DAG [open $output_dag_name "w"]
set output_failed_name "${work_dir}/failed_MainSubmit_${iteration}.txt"
set OUTPUT_FAILED [open $output_failed_name "w"]

set A_jobs {}

set dir_count 0

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


while { ![eof $TARGETS_FILE] } {
	gets $TARGETS_FILE line
	if {$line == ""} { continue }
	
	set values [ split $line "\t"]
	
	set directory "${work_dir}/idx[string trim [lindex $values $column_line_id] {"}]"
		
	#"
	if { [file exists "${directory}/run_${iteration}.dag"] } {
		set TEMP_DAG [open "${directory}/run_${iteration}.dag" "r"]
		while { ![eof $TEMP_DAG] } {	
			gets $TEMP_DAG line
			if {[regexp {^JOB ([^ ]*) } $line {} id]} {
				lappend A_jobs $id
				}
			puts $OUTPUT_DAG $line
			}
		close $TEMP_DAG
		incr dir_count
		} else {
		puts $OUTPUT_FAILED "${directory}/run_${iteration}.dag" 
		}
		
	} 
	
close $OUTPUT_FAILED
	
set outlier_condorB_sub "$work_dir/run_${iteration}_examine.sub"

puts $OUTPUT_DAG "JOB B0 $outlier_condorB_sub"
puts $OUTPUT_DAG "VARS B0 JobID=\"0\" argList=\"work_dir $work_dir iteration $iteration \""
puts $OUTPUT_DAG ""


puts -nonewline $OUTPUT_DAG "PARENT "
puts -nonewline $OUTPUT_DAG "$A_jobs"
puts $OUTPUT_DAG " CHILD B0"
	
close $OUTPUT_DAG

set std_out_err_dir "$work_dir/std_out_err"

set outlier_condorB_sub_text {
universe = vanilla 
executable = $script_dir/FstatFollowUpExamine_Inj.tcl
output = $std_out_err_dir/node_${iteration}_examine.out.\$(JobID)
error = $std_out_err_dir/node_${iteration}_examine.err.\$(JobID)
log = ${log_file_dir}/${search_name}_examine_${iteration}.log
notify_user=jcbetzwieser@gmail.com
arguments = \$(arglist)
queue
}

set FILE [open $outlier_condorB_sub "w"]
puts $FILE [subst $outlier_condorB_sub_text]
close $FILE


exec condor_submit_dag $output_dag_name

