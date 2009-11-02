#!/usr/bin/env tclsh

foreach {var value} {
	source_dir unknown
	dest_dir unknown
	} {
	global $var
	set $var [subst $value]
	}

foreach {var value} $argv {
	global $var
	set $var $value
	}

puts $dest_dir
puts $source_dir

file mkdir "$dest_dir"


foreach filename [lsort [glob -directory "$source_dir" *.sft.txt]] {
	set gps_time [lindex [split $filename "-"] end-1]
	set IFO [lindex [split $filename "_"] 2]
	set DATA_FILE [open $filename "r"]
	gets $DATA_FILE line

	while {![eof $DATA_FILE]} {
		gets $DATA_FILE line
		set split_line [split $line "\t"]
		set OUTPUT_FILE [open "${dest_dir}/psd_[lindex $split_line 0]Hz.txt" "a"]
		puts $OUTPUT_FILE "$gps_time\t$IFO\t[ lindex $split_line end-1]\t[ lindex $split_line end]"
		close $OUTPUT_FILE
		flush stdout
		}
	close $DATA_FILE
	}