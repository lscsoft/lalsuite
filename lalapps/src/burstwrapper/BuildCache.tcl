#!/usr/bin/env tclshexe
package require segments

set segFile [ lindex $argv 0 ]
set pat [ lindex $argv 1 ]
set server [ lindex $argv 2 ]
set observatory [ lindex $argv 3 ]
set type [ lindex $argv 4 ]
set fil [ lindex $argv 5 ]

set seglist [SegRead $segFile]

set segs [ SegSort $seglist ]

set T0 [ lindex [ lindex $segs 0 ] 0 ]
set T1 [ lindex [ lindex $segs end ] 1 ]

puts "type: $type, IFO: $observatory, start: $T0, end: $T1"

if { [ catch { exec /bin/bash -c "source $pat/setup.sh; $pat/ldg-client/bin/LSCdataFind --server $server --observatory $observatory --type $type --gps-start-time $T0 --gps-end-time $T1 --url-type file --lal-cache >> $fil" } cout ] } {
    puts $cout
    exit 1
}
