#!/usr/bin/env tclsh


#ExamineFollowUpResults.tcl - Examines the results from a FStatFollowUp.tcl run, finds the largest events, makes a few plots, and calls FStatFollowUp.tcl if another iteration is still possible. 

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


set param_vec {ra dec freq f1 twoF ref_time}

set FINAL_FA_FILE [open "$false_alarm_results_directory/FA_${search_name}_${iteration}.csv" "w"]
set FINAL_EVENTS_FILE [open "$events_directory/${search_name}_${iteration}" "w"]

puts $FINAL_FA_FILE "line_id\tfound_ra\tfound_dec\tfound_freq\tfound_f1\tparam_time\tfound_twoF\tra\tdec\tf0\tf1\ttotal_templates\tfa_prob"

puts $FINAL_EVENTS_FILE "line_id\tfound_ra\tfound_dec\tfound_freq\tfound_f1\tparam_time\tfound_twoF\tra\tdec\tf0\tf1\ttotal_templates\tfa_prob"

foreach {directory} [glob -directory $work_dir "idx*"] {
	regexp {.*/idx([^/]*)/?$} $directory {} line_id
	
	foreach var $param_vec {
	set temp_$var 0
	set best_$var 0
	}
	set number_of_templates 0
	
	foreach loudest_file [glob -nocomplain -directory ${directory}/results_${iteration}/ "loudest_*"] {
		
		set LOUDEST_FILE [open $loudest_file "r"]
		while { ![eof $LOUDEST_FILE] } {
			gets $LOUDEST_FILE line
			if {[regexp  {^%} $line]} {continue}
			set ident [string trim [lindex [split $line "="] 0]]
			if {$ident == "Alpha"} { set temp_ra [string trim [lindex [split $line "="] 1] ";"] }
			if {$ident == "Delta"} { set temp_dec [string trim [lindex [split $line "="] 1] ";"] }
			if {$ident == "Freq"} { set temp_freq [string trim [lindex [split $line "="] 1] ";"] }
			if {$ident == "f1dot"} { set temp_f1 [string trim [lindex [split $line "="] 1] ";"] }
			if {$ident == "twoF"} { set temp_twoF [string trim [lindex [split $line "="] 1] ";"] }
			if {$ident == "refTime"} {set temp_ref_time [string trim [lindex [split $line "="] 1] ";"] }
			}
		close $LOUDEST_FILE
		if {$temp_twoF > $best_twoF} {
			foreach var $param_vec {
				set best_$var [set [subst {temp_$var}]]
				}
			}
		}
	
	foreach histogram_file [glob -nocomplain -directory ${directory}/results_${iteration}/ "hist_*"] {
		set HIST_FILE [open $histogram_file "r"]
		while { ![eof $HIST_FILE] } {
			gets $HIST_FILE line
			if {[regexp  {^%} $line]} {continue}
			if {$line == ""} {continue}
			set number_of_templates [expr $number_of_templates + [lindex $line 2]]
			}
		close $HIST_FILE
		}
	
	set F_cutoff [expr $best_twoF/2.0]
	set alpha [expr ((1 + $F_cutoff)/2.0)*exp(-$F_cutoff)]
	if { $alpha < 1e-3 } {
		set log1ma -$alpha 
		} else {
		set log1ma [expr log(1-$alpha)]
		}
	set False_alarm_prob [expr 1 - exp($number_of_templates*$log1ma)]

	
	set weight 0
	set weighted_ra 0
	set weighted_dec 0
	set weighted_f0 0
	set weighted_f1 0
	set weighted_x_vec 0
	set weighted_y_vec 0
	set weighted_z_vec 0
	foreach results_file [glob -nocomplain -directory ${directory}/results_${iteration}/ "result_*"] {
		set RESULTS_FILE [open $results_file "r"]
		while { ![eof $RESULTS_FILE] } {
			gets $RESULTS_FILE line
			if {[regexp  {^%} $line]} {continue}
			if {$line == ""} {continue}
			set twoF [lindex [split $line " "] 6]
			set weight_factor [expr $twoF]
			if { ($twoF > [expr 0.75 * $best_twoF]) } {
				set template_ra [lindex [split $line " "] 1]
				set template_dec [lindex [split $line " "] 2]
				set weight [expr $weight + $weight_factor]
				set weighted_f0 [expr $weighted_f0 + [lindex [split $line " "] 0] * $weight_factor]
				set weighted_f1 [expr $weighted_f1 + [lindex [split $line " "] 3] * $weight_factor]
				set weighted_x_vec [expr $weighted_x_vec + [expr cos($template_ra) * cos($template_dec)] * $weight_factor]
				set weighted_y_vec [expr $weighted_y_vec + [expr sin($template_ra) * cos($template_dec)] * $weight_factor]
				set weighted_z_vec [expr $weighted_z_vec + [expr sin($template_dec)] * $weight_factor]
				}
			}
		close $RESULTS_FILE
		}
	if {$weight != 0} {
	
		set weighted_ra [expr atan2($weighted_y_vec /$weight,$weighted_x_vec/$weight)]
		if {$weighted_y_vec == 0.0 & $weighted_x_vec == 0.0} {
			set weighted_ra 0
			}
		set weighted_dec [expr asin($weighted_z_vec / $weight)]
		set weighted_f0 [expr $weighted_f0 / $weight]
		set weighted_f1 [expr $weighted_f1 / $weight]
		} else {
		set weighted_ra NA
		set weighted_dec NA
		set weighted_f0 NA
		set weighted_f1 NA
		}
	
		
	puts $FINAL_FA_FILE "${line_id}\t${best_ra}\t${best_dec}\t${best_freq}\t${best_f1}\t${best_ref_time}\t${best_twoF}\t${weighted_ra}\t${weighted_dec}\t${weighted_f0}\t${weighted_f1}\t${number_of_templates}\t${False_alarm_prob}"
	if {$False_alarm_prob < $false_alarm_cutoff} {
		puts $FINAL_EVENTS_FILE "${line_id}\t${best_ra}\t${best_dec}\t${best_freq}\t${best_f1}\t${best_ref_time}\t${best_twoF}\t${weighted_ra}\t${weighted_dec}\t${weighted_f0}\t${weighted_f1}\t${number_of_templates}\t${False_alarm_prob}"
		}
	}
	
	
close $FINAL_FA_FILE 
close $FINAL_EVENTS_FILE
	
	


