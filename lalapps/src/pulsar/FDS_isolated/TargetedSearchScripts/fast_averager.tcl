#!/usr/bin/env tclsh


while { ![eof stdin] } {
	gets stdin filename
	if { $filename == "" } { continue }
	if { [regexp {^#} $filename ] } { continue }

	set FILE [open $filename "r"]
	fconfigure $FILE -translation binary -encoding binary

	set header [read $FILE 48]
	binary scan $header "d1i1i1d1i1i1i2a2c2i1" key gps_sec gps_nsec tbase first_frequency_index nsamples crc detector padding1 comment_length

	foreach var {key gps_sec gps_nsec tbase first_frequency_index nsamples crc detector padding1 comment_length} {
		puts "$var=[set $var]"
		}

	set comment [read $FILE $comment_length]
	puts "comment=\"$comment\""

	set raw_data [read $FILE [expr $nsamples*8]]
	close $FILE


	binary scan $raw_data "f[expr $nsamples*2]" data

	puts [lrange $data 0 100]
	puts [llength $data]

	set f0 [expr $first_frequency_index/$tbase]
	set step [expr round($tbase)]
	for { set i 0 } { $i < $nsamples } { incr i $step } {
		set chunk {}
		set sum 0.0

		foreach {x y} [lrange $data [expr 2*$i] [expr 2*($i+$step)-1]] {
			set p [expr $x*$x+$y*$y]
			lappend chunk $p
			set sum [expr $sum +$p]
			}
		set chunk_sorted [lsort -real $chunk]
		set median [lindex $chunk_sorted [expr $step/4]]
		puts "[expr $f0+$i/$tbase] $median [expr $sum*2.0/$step]"
		}


	}
