#!/usr/bin/env tclsh

# 1st argument is output
# other arguments are eff files to merge

if { [ llength $argv ] == 1 } {
    set argv [split $argv " "]
}

set ofile [ lindex $argv 0 ]

set f1 [ lindex $argv 1 ]

set fid [ open $f1 "r" ]

while { [ gets $fid line ] } {

    set lin [ split $line ',' ]

    if { [ llength $lin ] < 1 } {
	if { [ eof $fid ] } {
	    break
	}
	continue
    }

    set ll [ string length $line ]
    set le [ expr $ll - 1 ]
    for { set nc 0 ; set i $le } { $i > 0 && $nc < 2 } { incr i -1 } {
	if { [ string index $line $i ] == "," } {
	    incr nc 1
	}
    }

    set str [ string range $line 0 $i ]
    set tocs($str) [ lrange $lin end-1 end ]

}

close $fid


for { set ifile 2 } { $ifile < $argc } { incr ifile 1 } {

    set f1 [ lindex $argv $ifile ]

    puts $f1

    set fid [ open $f1 "r" ]

    while { [ gets $fid line ] } {
	
	set lin [ split $line ',' ]

	if { [ llength $lin ] < 1 } {
	    if { [ eof $fid ] } {
		break
	    }
	    continue
	}

	set ll [ string length $line ]
	set le [ expr $ll - 1 ]
	for { set nc 0 ; set i $le } { $i > 0 && $nc < 2 } { incr i -1 } {
	    if { [ string index $line $i ] == "," } {
		incr nc 1
	    }
	}

	set toc [ string range $line 0 $i ]
	regsub -all {\\} $toc "\\\\\\" toc
	set val [ lrange $lin end-1 end ]
	set rep [ array get tocs $toc ]

	set to [ lindex $rep 0 ]
	set va [ lindex $rep 1 ]

	set va [ lreplace $va 0 0 [ expr [lindex $va 0] + [lindex $val 0]] ]
	set va [ lreplace $va 1 1 [ expr [lindex $va 1] + [lindex $val 1]] ]

	set tocs($to) $va
    }

    close $fid

}

set fid [ open $ofile w ]

foreach { toc val } [array get tocs] {
    set v1 [ lindex $val 0 ]
    set v2 [ lindex $val 1 ]
    puts $fid "$toc,$v1,$v2"
}

close $fid
