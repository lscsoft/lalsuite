#!/usr/bin/env tclshexe

set file [ lindex $argv 0 ]

set fid [ open $file "r" ]

set dat [ list ]

gets $fid line

while { [ eof $fid ] == 0 } {

    if { [ regexp {(/data/node[0-9]+/jsylvest/[^/]+/job[LH][12]\.jsylvest\.)([0-9]+)(\.[0-9]+\.bin) ([0-9]+)} $line junk p1 T p2 Tf ] } {
	lappend dat [ list $p1 $T $p2 $Tf ]
    }
    
    gets $fid line

}

close $fid

set dat [ lsort -index 1 -integer $dat ]

set N [ expr [llength $dat] - 1]

set i 0

while { $i < $N } {
    for { set j [expr $i + 1] } { $j < [llength $dat] } { incr j 1 } {

	if { [ lindex [ lindex $dat $i ] 1 ] != [ lindex [ lindex $dat $j ] 1 ] } {
	    set i $j
	    break
	}


	if { [ lindex [ lindex $dat $i ] 3 ] < [ lindex [ lindex $dat $j ] 3 ] } {

	    puts [ lindex [ lindex $dat $i ] 0 ][ lindex [ lindex $dat $i ] 1 ][ lindex [ lindex $dat $i ] 2 ]

	} else {

	    puts [ lindex [ lindex $dat $j ] 0 ][ lindex [ lindex $dat $j ] 1 ][ lindex [ lindex $dat $j ] 2 ]

	}

    }
}