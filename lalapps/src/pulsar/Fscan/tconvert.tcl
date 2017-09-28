#-- "tconvert.tcl"
#-- GPS-UTC time conversion tools
#-- Extracted from 'tconvert' script by Peter Shawhan, January 2002

package provide tconvert 6.0

#==============================================================================
namespace eval tconvert {

    #-- Initialize namespace variables
    variable gotLeapSeconds 0
    variable gpsEpoch 0
    variable gpsLeaps {}
    variable sysLeaps {}

    #-- Commands to be imported with "namespace import tconvert::*"
    namespace export tconvert SysToGPS GPSToSys
}


##=========================================================================
## Name: tconvert
##

proc tconvert::tconvert { args } {
    variable gotLeapSeconds
    variable gpsLeaps
    variable sysLeaps

    ;##- Set defaults
    set gmtflag 1
    set format "%b %d %Y %T %Z"
    set timespec ""
    set debug 0

    ;##- Parse arguments
    set option ""
    foreach arg $args {
	if { [regexp {^-[A-Za-z]} $arg] } {
	    set option $arg
	    set val ""
	} else {
	    set val $arg
	}

	if { $option == "" } {
	    ;##- Append this to the time specification string
	    if { $timespec == "" } {
		set timespec $arg
	    } else {
		append timespec " " $arg
	    }

	} else {
	    ;##- Handle this option
	    switch -- $option {

		-local -
		-l {
		    set gmtflag 0
		    set option ""   ;##- The '-l' flag takes no argument
		}

		-gmt {
		    if { $val != "" } {
			set gmtflag $val
			set option ""  ;##- This flag takes a one-word argument
		    }
		}

		-debug -
		-d {
		    set debug 1
		    set option ""   ;##- The '-d' takes no argument
		}

		-format -
		-f {
		    if { $val != "" } {
			set format $val
			set option ""  ;##- This flag takes a one-word argument
		    }
		}

		default {
		    return -code error "tconvert ERROR: Invalid option $option"
		}

	    }
	}
    }


    ;##- If we haven't yet gotten leap-second information, do it now
    if { $gotLeapSeconds == 0 } {
	if { [catch {GetLeapSeconds $debug} msg] } {
	    return -code error "Error getting leap-second information. \
		    Message: $msg"
	}
	set gotLeapSeconds 1
	if {$debug} DumpLeapSeconds
    }

    ;##- Figure out whether we should interpret the timespec as GPS seconds
    ;##- or as a date-time string
    set timespec [string trim $timespec]
    if { [string is space $timespec] } { return "" }
    if { [regexp {[^0-9\.+\-*/()\s]} $timespec] } {
	;##- Contains a chararacter which is not a digit, decimal,
	;##- arithmetic symbol, parenthesis, or space
	set type "datetime"
    } elseif { [regexp {^-?[\d\.]+$} $timespec] } {
	;##- Consists entirely of digits and/or decimal points
	set gpstime $timespec
	set type "GPS"
    } elseif { [regexp {^\d{9,10}} $timespec] \
	    && ! [catch {expr $timespec} gpstime] \
	    && [string is double $gpstime] } {
	;##- Begins with 9 or 10 digits, and is an arithmetic expression
	set type "GPS"
    } else {
	set type "datetime"
    }

    if {$debug} { puts "Interpreting timespec ($timespec) as $type" }

    ;##- If the input was a GPS time, translate it now
    if { $type == "GPS" } {

	;##- Check for fractional seconds
	if { [regexp {^(-?\d+)(\.\d*)} $gpstime match intsec fracsec] } {
	    if {$debug} { puts "Found fractional GPS seconds: $fracsec" }
	    set gpstime $intsec
	    ;##- Modify output format to display the fractional part at the
	    ;##- appropriate place
	    regsub -all "%r" $format "%I:%M:%S %p" format
	    regsub -all "%T" $format "%H:%M:%S" format
	    regsub -all "%S" $format "%S$fracsec" format
	}

	if { ! [string is integer $gpstime] \
		|| $gpstime < 0 || $gpstime > 2147483647 } {
	    return -code error "Invalid GPS time"
	}
	if { $gpstime > 0 && $gpstime < 200000000 } {
	    return -code error "Unreasonable GPS time (prior to 1986)"
	}
	if { $gpstime > 1830000000 } {
	    return -code error "Unreasonable GPS time (after 2037)"
	}
	set systime [GPSToSys $gpstime $debug]

	;##- If this is an exact leap second, then the output should be a time
	;##- like "23:59:60".  Modify the format string to accomplish this.
	if { [lsearch -exact $gpsLeaps $gpstime] != -1 } {
	    if {$debug} { puts "This is an exact leap second" }
	    regsub -all "%r" $format "%I:%M:%S %p" format
	    regsub -all "%T" $format "%H:%M:%S" format
	    regsub -all "%S" $format "60" format
	}

	set outstring [clock format $systime -format $format -gmt $gmtflag]

	;##- Replace "GMT" with "UTC"
	regsub -all "GMT" $outstring "UTC" outstring

	return $outstring
    }

    ;##- If we get to this point, the input is a date/time string
    set leapcorr 0
    set fracsec ""
    set delsign 1

    ;##- Modify string to avoid a counterintuitive case: Tcl considers
    ;##- "jan 23 2003" to be "jan 23 20:03", but I want it to be
    ;##- interpreted as "jan 23, 2003".
    regsub -nocase {([a-z]{3,}\s+\d{1,2})\s+((?:19|20)\d\d)} $timespec \
	    {\1, \2} timespec

    ;##- Another special case: allow the form "jan 23, 20:03"
    regsub -nocase {([a-z]{3,}\s+\d{1,2})\s*,\s+(\d{1,2}:\d{1,2})} $timespec \
	    {\1 \2} timespec

    ;##- If input string begins with a date/time and ends with a relative
    ;##- number, assume it is a number of seconds
    if { [regexp {[^\d\s\+\-\*/].*[+-]\s*[\d\*]+$} $timespec] } {
	append timespec "sec"
    }

    ;##- If it just ends with "s", also assume this is seconds
    if { [regexp {[^\d\s\+\-\*/].*[+-]\s*[\d\*]+\s*s\s*$} $timespec] } {
	set timespec [string trim $timespec]
	append timespec "ec"
    }

    if { [catch {clock scan $timespec -gmt $gmtflag} systime] } {
	;##- 'clock scan' fails!  There are two special cases to check.

	if {$debug} {
	    puts "Initial 'clock scan' fails.  Checking special cases ..."
	}

	;##- First, see if this seems to be a leap second, during which
	;##- a UTC clock should properly read "23:59:60"
	if {$debug} {
	    puts -nonewline "  Checking if this is has the format of an\
		    exact leap second ..."
	}
	if { [regsub {:59:60} $timespec {:59:59} timespec2] \
		|| [regsub {5960} $timespec {5959} timespec2] } {
	    if {$debug} { puts " yes" }
	    set leapcorr -1
	} else {
	    if {$debug} { puts " no" }
	    set timespec2 $timespec
	}

	;##- If 'clock scan' still fails, check for fractional seconds
	;##- and remove if present
	if {$debug} { puts -nonewline "  Checking for fractional seconds ..." }
	if { [catch {clock scan $timespec2 -gmt $gmtflag} systime] } {
	    if { ! [regexp {^([^.]*)(\d)(\.\d*)([^.]*)$} $timespec2 \
		    match part1 digit fracsec part2] } {
		if {$debug} { puts " no decimal point found" }
		return -code error "Unable to parse date/time string"
	    }
	    ;##- See if removing the fractional part makes the string parsable
	    if { [catch {clock scan "${part1}${digit}${part2}" -gmt $gmtflag} \
		    systime] } {
		if {$debug} {
		    puts " found a decimal fraction, but removing it doesn't\
			    make string parsable"
		}
		return -code error "Unable to parse date/time string"
	    }

	    ;##- The string is parsable!
	    ;##- Make sure the fractional part was really on the seconds, by
	    ;##- modifying the least-significant-digit on the seconds and
	    ;##- verifying that the parsed time changes appropriately
	    if { $digit != 9 } {
		set newdigit [expr {$digit+1}]
		set delexpect 1
	    } else {
		set newdigit 8
		set delexpect -1
	    }
	    set systime2 \
		    [clock scan "${part1}${newdigit}${part2}" -gmt $gmtflag]
	    if { [expr {abs($systime2-$systime)}] != 1 } {
		if {$debug} {
		    puts " found a decimal fraction, but it is not on the\
			    number of seconds"
		}
		return -code error "Unable to parse date/time string"
	    }
	    set delsign [expr {($systime2-$systime)*$delexpect}]

	    ;##- The fractional part really WAS on the seconds!
	    ;##- Now we can go on.
	    if {$debug} {
		puts " successfully stripped out the fractional part\
			($fracsec)"
	    }

	}

	;##- If we get to this point, 'clock scan' has succeeded

	;##- If this looks like a leap second, make sure it corresponds to an
	;##- actual leap second
	if { $leapcorr == -1 } {
	    incr systime
	    if { [lsearch -exact $sysLeaps $systime] != -1 } {
		;##- This is an actual leap second!  We'll need to correct the
		;##- GPS time we come up with
		if {$debug} {
		    puts "Verified that this is an actual leap second"
		}
		set leapcorr -1
	    } else {
		if {$debug} { puts "This is NOT an actual leap second" }
		return -code error "That is not an actual leap second"
	    }
	}
    }

    ;##- If a U.S. time zone was specified, check daylight vs. standard time
    ;##- (only on Unix)
    if { [regexp -nocase {(?:\A|[^y])([ecmp][sd]t)} $timespec match inzone] \
	    && ! [regexp -nocase {(windows|macos|darwin)} \
	    $::tcl_platform(os)] } {
	set inzone [string toupper $inzone]
	set inzone1 [string index $inzone 0]
	if {$debug} {
	    puts -nonewline "Checking whether specified time zone ($inzone) is\
		    valid for this time ..."
	}

	switch -exact -- $inzone1 {
	    "E" { set checkTZ "EST5EDT" }
	    "C" { set checkTZ "CST6CDT" }
	    "M" { set checkTZ "MST7MDT" }
	    "P" { set checkTZ "PST8PDT" }
	}

	if { [info exists ::env(TZ)] } {
	    set saveZone $::env(TZ)
	} else {
	    set saveZone ""
	}

	set tzChanged 0
	if { ! [string equal $saveZone $checkTZ] } {
	    set ::env(TZ) $checkTZ
	    set tzChanged 1
	}

	set outzone [string toupper [clock format $systime -format %Z]]

	if { $tzChanged } {
	    ;##- Restore the original time zone
	    if { $saveZone != "" } {
		set ::env(TZ) $saveZone
	    } else {
		unset ::env(TZ)
	    }
	}

	if { ! [string equal $inzone $outzone] } {
	    if {$debug} { puts " no" }
	    return -code error "That time does not occur during $inzone"
	} else {
	    if {$debug} { puts " yes" }
	}

    }

    ;##- OK, convert to GPS time and print it out
    set gpstime [expr {[SysToGPS $systime $debug]+$leapcorr}]
    if { $gpstime < 0 } {
	return -code error "Time is prior to GPS 0 (Jan 6, 1980)"
    }
    if { $fracsec != "" } {
	if { $delsign == 1 || $fracsec == "." } {
	    set gpstime "$gpstime$fracsec"
	} else {
	    set length [expr {[string length $fracsec]-1}]
	    set modfrac [format %.${length}f [expr {1.0-$fracsec}]]
	    regsub {^0+\.} $modfrac {.} modfrac
	    incr gpstime -1
	    set gpstime "$gpstime$modfrac"
	}
    }

    return $gpstime
}


##=========================================================================
## Name: GPSToSys
##
## Convert a GPS time to system time

proc tconvert::GPSToSys { gpstime {debug 0} } {
    variable gpsEpoch
    variable gotLeapSeconds
    variable gpsLeaps

    ;##- If we haven't yet gotten leap-second information, do it now
    if { $gotLeapSeconds == 0 } {
	if { [catch {GetLeapSeconds $debug} msg] } {
	    return -code error "tconvert: Error getting leap-second\
		    information.  Message: $msg"
	}
	set gotLeapSeconds 1
	if {$debug} DumpLeapSeconds
    }

    set nleaps 0
    foreach gpsleap $gpsLeaps {
	if { $gpstime >= $gpsleap } { incr nleaps }
    }

    set systime [expr {$gpstime+$gpsEpoch-$nleaps}]

    return $systime
}


##=========================================================================
## Name: SysToGPS
##
## Convert system time to GPS time

proc tconvert::SysToGPS { systime {debug 0} } {
    variable gpsEpoch
    variable gotLeapSeconds
    variable sysLeaps

    ;##- If we haven't yet gotten leap-second information, do it now
    if { $gotLeapSeconds == 0 } {
	if { [catch {GetLeapSeconds $debug} msg] } {
	    return -code error "tconvert: Error getting leap-second\
		    information.  Message: $msg"
	}
	set gotLeapSeconds 1
	if {$debug} DumpLeapSeconds
    }

    set nleaps 0
    foreach sysleap $sysLeaps {
	if { $systime >= $sysleap } { incr nleaps }
    }

    set gpstime [expr {$systime-$gpsEpoch+$nleaps}]

    return $gpstime
}


##=========================================================================
## Name: GetLeapSeconds
##

proc tconvert::GetLeapSeconds { {debug 0} } {
    variable gpsEpoch
    variable gpsLeaps
    variable sysLeaps

    set gpsEpoch [clock scan "jan 6, 1980" -gmt 1]
    set timenow [clock seconds]

    set configdir ""
    set sysLeaps {}

    ;##- Figure out where the leap-seconds file lives.
    ;##- First, see if the TCLEAPSDIR environment variable is set
    #if { [info exists ::env(TCLEAPSDIR)] }
    if { 1 } {
	#if {$debug} {
	#    puts "Using TCLEAPSDIR environment variable for location\
	#	    of config dir"
	#}
	#set configdir $::env(TCLEAPSDIR)
        set configdir [pwd]
    } elseif { [info exists ::env(LIGOTOOLS)] } {
	;##- Use the LIGOTOOLS environment variable to locate the config dir;
	;##- then tcleaps.txt is in $LIGOTOOLS/config/public
	if {$debug} {
	    puts "Using LIGOTOOLS environment variable to determine location\
		    of config dir"
	}
	set configdir $::env(LIGOTOOLS)/config/public

	#-- If tcleaps.txt exists here but it's not writable, use home
	#-- directory instead
	if { [file exists $configdir/tcleaps.txt] && \
		 ! [file writable $configdir/tcleaps.txt] && \
		 [info exists ::env(HOME)] } {
	    if {$debug} { puts "$configdir/tcleaps.txt is\
                    not writable, so using home directory" }
	    set configdir $::env(HOME)
	}

    } else {
        ;##- Next, try to do so relative to this script's location
	set script [ScriptLocation]
	if {$debug} {puts "Script location is $script"}
	if { [regexp {^(.+)/packages/dataflow/[^/]+/} $script match ltdir] } {
	    set configdir "${ltdir}/config/public"
	    if {$debug} {
		puts "Using script's location to determine location\
			of config dir"
	    }
	} elseif { [info exists ::env(HOME)] } {
	    ;##- As a last resort, use the user's home directory
	    set configdir $::env(HOME)
	    if {$debug} { puts "Using home directory as config dir" }
	}
    }

    if {$debug} {
	if { $configdir != "" } {
	    puts "config dir is $configdir"
	    if { [file isdirectory $configdir] } {
		puts "config dir exists"
	    } else {
		puts "config dir does not exist"
	    }
	    puts "leap-seconds file is $configdir/tcleaps.txt"
	} else {
	    puts "Cannot determine location of config dir"
	}
    }

    ;##- Check whether the leap-seconds file exists and is readable
    if { $configdir != "" } {
	set leapfile $configdir/tcleaps.txt
    } else {
	set leapfile ""
    }
    if { $leapfile != "" && [file readable $leapfile] } {
	if {$debug} { puts "tcleaps.txt exists and is readable" }
	set expiration [ReadLeapSecFile $leapfile $debug]

	;##- Check whether the leap-second info has expired
	if { $timenow > $expiration } {
	    ;##- The leap-seconds file has expired, so we will have to try to
	    ;##- update it using information from the web.  But first save what
	    ;##- we have, in case we're unable to get any updated information.
	    if {$debug} {puts "  ... has EXPIRED"}

	    set sysLeapsExpired $sysLeaps
	    set gpsLeapsExpired $gpsLeaps
	    set sysLeaps {}
	    set gpsLeaps {}
	} else {
	    if {$debug} {puts "  ... is still valid"}
	    return
	}

    } else {
	if {$debug} { puts "tcleaps.txt does not exist or is not readable" }
    }

    ;##- If we get to this point, then the leap-seconds file either does not
    ;##- exist, is unreadable for some reason, or has expired.
    ;##- So read leap-second data from the web.

    if { [catch {package require http} msg] } {
	;##- We can't read from the web.  If we were able to read from the disk
	;##- file but that info had "expired", then print a warning message
	;##- and use it anyway.
	if { [info exists sysLeapsExpired] } {
	    #puts stderr "tconvert WARNING: Leap-second info in $leapfile is no\
	#	    longer certain to be valid, and there was an error loading\
	#	    the Tcl http package to get updated info from the web. \
	#	    Continuing with possibly-outdated info."
	    set sysLeaps $sysLeapsExpired
	    set gpsLeaps $gpsLeapsExpired
	    return
	} else {
	    return -code error "tconvert ERROR: There is no leap-second info\
		    cached on disk, and there was an error loading the Tcl\
		    http package to get info from the web."
	}
    }

    set contents ""

    ;##- Try to get the leapseconds file from various LDAS web servers
    foreach host [list www.ldas.ligo-wa.caltech.edu \
	    www.ldas.ligo-la.caltech.edu \
	    www.ldas-dev.ligo.caltech.edu \
	    www.ldas.ligo.caltech.edu \
	    ldas.mit.edu \
	    www.ldas-test.ligo.caltech.edu ] {
	if {$debug} {
	    puts -nonewline "Trying to get leapseconds file from $host ..."
	    flush stdout
	}

	set url "http://$host/ldas_outgoing/jobs/leapseconds"
	if { [catch {GetUrlContents $url 5000} contents] } {
	    if {$debug} {puts " FAILED"}
	    set contents ""
	    continue
	} else {
	    if {$debug} { puts " success"}
	    break
	}

    }

    if { $contents == "" } {
	;##- We can't read from the web.  If we were able to read from the disk
	;##- file but that info had "expired", then print a warning message
	;##- and use it anyway.
	if { [info exists sysLeapsExpired] } {
	    #puts stderr "tconvert WARNING: Leap-second info in $leapfile is no\
	#	    longer certain to be valid, and we were unable to get\
	#	    updated info from any LDAS web server. \
	#	    Continuing with possibly-outdated info."
	    set sysLeaps $sysLeapsExpired
	    set gpsLeaps $gpsLeapsExpired
	    return
	} else {
	    return -code error "tconvert ERROR: There is no leap-second info\
		    cached on disk, and we were unable to get updated info\
		    from any LDAS web server."
	}
    }

    ;##- Parse the contents to construct a list of leap seconds since 1980

    set pat {^ *((?:19|20)\d\d) +([A-Za-z]{3}) +(\d{1,2})}
    foreach {match year month mdate} \
	    [regexp -all -inline -line $pat $contents] {
	set time [clock scan "$month $mdate, $year" -gmt 1]
	if { $time > $gpsEpoch } {
	    lappend gpsLeaps [expr {$time-$gpsEpoch+[llength $sysLeaps]}]
	    lappend sysLeaps $time
	}
    }
    if {$debug} {
	puts "Found [llength $gpsLeaps] leap seconds since 6 Jan 1980"
    }

    ;##- Set a default expiration time just one day from now.  Hopefully this
    ;##- will be updated below; otherwise, we'll just check again tomorrow.
    set expiration [expr {$timenow+24*3600}]

    ;##- Finally, check the latest "Bulletin C" from the IERS to see how far
    ;##- in the future we can be sure that this list is valid.

    if {$debug} {
	puts -nonewline "Retrieving latest Bulletin C ..."
	flush stdout
    }
    set url "http://hpiers.obspm.fr/eoppc/bul/bulc/bulletinc.dat"
    if { [catch {GetUrlContents $url 5000} contents] } {
	if {$debug} {puts " FAILED"}
	;##- Failed to get Bulletin C.  But we know that leap seconds are
	;##- always at least 6 months apart, so set the expiration date to be
	;##- 5 months after the last known leap-second, if this is later than
	;##- the default we set above.
	set expir1 [expr {[lindex $sysLeaps end]+5*30*24*3600}]
	if { $expir1 > $expiration } { set expiration $expir1 }

    } else {
	if {$debug} {puts " success"}
	;##- Parse the bulletin.  ParseBulletinC returns an expiration date
	set expir1 [ParseBulletinC $contents $debug]
	if { $expir1 > $expiration } { set expiration $expir1 }
    }

    if {$debug} {
	puts "New expiration date is [clock format $expiration -gmt 1]"
    }

    ;##- If possible, update the leap-seconds file on local disk
    if { $leapfile != "" } {

	if { [catch \
		{WriteLeapSecFile $configdir $leapfile $expiration $debug} \
		msg] } {
	    ;##- We were unable to write the updated disk file.  If the
	    ;##- current disk file has expired, print a warning message.
	    if { [info exists sysLeapsExpired] } {
		puts stderr "tconvert NOTICE: Leap-second info in $leapfile\
			is no longer certain to be valid; We got valid\
			information from the web, but were unable to update\
			the local cache file: $msg"
	    }
	}

    }

    return
}


##=========================================================================
## Name: DumpLeapSeconds
##

proc tconvert::DumpLeapSeconds {} {
    variable gpsLeaps
    variable sysLeaps

    puts "Leap seconds:"
    foreach sys $sysLeaps gps $gpsLeaps {
	puts "[format %12d $sys]  \
		[clock format $sys -format "%d %b %Y" -gmt 1]  \
		GPS=[format %10d $gps]"
    }

    return
}


##=========================================================================
## Name: ReadLeapSecFile
##
## Sets tconvert::sysLeaps and tconvert::gpsLeaps
## Returns expiration time if file was read successfully, or 0 if there
## was some problem reading it.

proc tconvert::ReadLeapSecFile { leapfile {debug 0} } {
    variable gpsEpoch
    variable gpsLeaps
    variable sysLeaps

    if { [catch {open $leapfile r} fid] } {
	;##- Error opening file
	if {$debug} {puts "Error opening tcleaps.txt"}
	return 0
    }

    if { [catch {read $fid} contents] } {
	;##- Error reading file
	if {$debug} {puts "Error reading tcleaps.txt"}
	catch { close $fid }
	return 0
    }

    ;##- If we get here, file was read successfully
    if {$debug} {puts "Success reading tcleaps.txt"}
    close $fid

    ;##- Do some sanity checks on the file contents.  The first
    ;##- line should indicate the format, the second
    ;##- line should indicate the time through which this file is
    ;##- known to be valid, and one of the lines should be the
    ;##- system time corresponding to Jan 1, 1999

    set t1jan99 [clock scan "jan 1, 1999" -gmt 1]

    if { ! ( [regexp -nocase {^format (\d+)} $contents match format] \
	    && [regexp -nocase -lineanchor {^valid through +(-?\d+)} \
	    $contents match expiration] \
	    && [regexp -nocase -lineanchor {^gps epoch +(-?\d+)} \
	    $contents match epoch] \
	    && $epoch == $gpsEpoch
	    && [regexp -lineanchor "^$t1jan99" $contents] ) } {
	;##- File seems corrupted!
	if {$debug} {puts "  ... FAILS sanity checks"}
	return 0
    }

    ;##- If we get here, then the file looks good
    if {$debug} {
	puts "  ... passes sanity checks"
	puts "  ... has format $format"
	puts "  ... expires at [clock format $expiration -gmt 1]"
    }

    ;##- Get leap-second times
    foreach time [regexp -all -inline -lineanchor {^-?\d+} $contents] {
	if { $time > $gpsEpoch } {
	    lappend gpsLeaps [expr {$time-$gpsEpoch+[llength $sysLeaps]}]
	    lappend sysLeaps $time
	}
    }
    if {$debug} { puts "  ... contains [llength $gpsLeaps] leap seconds" }

    return $expiration
}


##=========================================================================
## Name: WriteLeapSecFile
##

proc tconvert::WriteLeapSecFile { configdir leapfile expiration {debug 0} } {
    variable gpsEpoch
    variable sysLeaps

    ;##- If the config directory does not exist, try to create it
    if { ! [file isdirectory $configdir] \
	    && [file writable [file dirname $configdir]] \
	    && ( ! [info exists $::env(HOME)] || $configdir != $::env(HOME) ) \
	} {
	if {$debug} {
	    puts -nonewline "Attempting to create config directory ..."
	    flush stdout
	}
	if { [catch {file mkdir $configdir}] } {
	    if {$debug} {puts " FAILED to make config dir"}
	    return -code error "Failed to create directory $configdir"
	}
	if { [catch {file attributes $configdir -permissions 0777}] } {
	    if {$debug} {
		puts " made config dir, but failed to set permissions"
	    }
	    ;##- Continue despite this error
	} else {
	    if {$debug} {puts " success"}
	}
    }

    ;##- Make sure directory exists and is writable now
    if { ! [file writable $configdir] } {
	if {$debug} { puts "config directory is not writable" }
	return -code error "No permission to write in directory $configdir"
    }
    if {$debug} { puts "config directory is writable" }

    ;##- If file already exists but we do not have write permission for it,
    ;##- then delete it first before rewriting it
    if { [file exists $leapfile] && ! [file writable $leapfile] } {
	if { [catch {file delete $leapfile}] } {
	    if {$debug} { puts "No permission to rewrite tcleaps.txt" }
	    return -code error "No permission to rewrite $leapfile"
	}
    }

    ;##- Open the file for writing
    if { [catch {open $leapfile w} fid] } {
	if {$debug} { puts "Error opening tcleaps.txt for writing" }
	return -code error "Error writing $leapfile"
    }

    ;##- Write the info to the file, then close it
    puts $fid "format 1"
    puts $fid "valid through $expiration"
    puts $fid "GPS epoch $gpsEpoch"
    foreach time $sysLeaps {
	puts $fid $time
    }
    close $fid
    if {$debug} { puts "Successfully wrote tcleaps.txt" }

    ;##- Try to make the file writable by all
    if { [ catch {file attributes $leapfile -permissions 0666} ] } {
	if {$debug} { puts "Error making tcleaps.txt writable by all" }
	;##- Continue despite this error
    }

    return
}


##=========================================================================
## Name: ScriptLocation
##
## Description:
##   Returns the full path to this script, following symbolic links if needed.

proc tconvert::ScriptLocation {} {

    set script $::argv0

    while { [file type $script] == "link" \
	    || [file pathtype $script] == "relative" } {
	if { [file type $script] == "link" } {
	    set linkto [file readlink $script]
	    if { [file pathtype $linkto] == "absolute" } {
		set script $linkto
	    } else {
		set script [file dirname $script]/$linkto
	    }
	} else {
	    set script "[pwd]/$script"
	}
    }

    ;##- Now rationalize the full path

    ;##- Condense adjacent slashes
    regsub -all {/{2,}} $script {/} script

    ;##- Remove redundant directories
    while { [regexp {/\./} $script] } {
	regsub {/\./} $script {/} script
    }

## This is incorrect in some obscure cases
##    ;##- Un-nest parent-directory references
##    while { [regexp {./\.\./} $script] } {
##	regsub {/[^/]+/\.\./} $script {/} script
##    }

    return $script
}


##=========================================================================
## Name: GetUrlContents
##

proc tconvert::GetUrlContents { url {timeout 0} } {

    if { [catch {http::geturl $url -timeout $timeout} httpvar] } {
	return -code error "GetUrlContents: http::geturl failed"
    }

    ;##- Parse the http response string to see if transfer succeeded
    upvar #0 $httpvar httpstate
    if { [info exists httpstate(http)] } {
	regexp -- {^([^\s]+)\s([^\s]+)\s(.*)$} $httpstate(http) \
		match httpversion httpcode status
	if { [info exists httpcode] && $httpcode == "200" } {
	    ;##- Success!
	    set contents $httpstate(body)
	    http::cleanup $httpvar
	    return $contents
	}
    }
    http::cleanup $httpvar

    return -code error "GetUrlContents: http::geturl failed"
}


##=========================================================================
## Name: ParseBulletinC
##
## This isn't very robust, so I hope they keep the same general formatting
## in the future!  Why can't some authority provide a nice, reliable,
## machine-readable URL which indicates how long the current set of leap
## seconds is guaranteed to be valid?

proc tconvert::ParseBulletinC { contents {debug 0} } {
    variable gpsEpoch
    variable gpsLeaps
    variable sysLeaps

    set expiration 0

    ;##- Get the date of the bulletin
    set bultime 0
    if { [regexp -nocase {paris, +(\d{1,2} +\w+ +20\d\d)} $contents \
	    match date] \
	    || [regexp -nocase {paris, (20\d\d +\w+ +\d{1,2})} $contents \
	    match date] } {

	if {$debug} {puts "Bulletin C date string is $date"}
	if { [catch {ScanEnglishOrFrenchDate $date} bultime] } {
	    if {$debug} {puts "  Error parsing date string"}
	    ;##- Error parsing date of bulletin.  We won't be able to set
	    ;##- the expiration date any farther into the future.
	} else {
	    if {$debug} {
		puts "Bulletin C date is [clock format $bultime -gmt 1]"
	    }
	}

    } else {
	if {$debug} {puts "Unable to identify Bulletin C date"}
    }

    ;##- Now, check the last UTC-TAI difference indicated in the Bulletin.
    ;##- This could be more up-to-date than the LDAS file we got.
    if { [regexp -nocase -linestop \
	    {from ([^,]+).* until further notice.*utc *- *tai *= *-} \
	    $contents match date] } {

	if {$debug} {puts "Bulletin C leap date is $date"}
	if { [catch {ScanEnglishOrFrenchDate $date} time] } {
	    if {$debug} {puts "  Error parsing leap date"}
	    ;##- We were unable to determine the time of the leap second
	    ;##- from Bulletin C, but we know that if a new leap second is
	    ;##- being announced, it will be ~6 months after the date of
	    ;##- the bulletin.  And it should show up in the LDAS
	    ;##- leapseconds file before that.  So set the expiration
	    ;##- time for our leap-second info to be 4 months after the
	    ;##- date of the bulletin, OR 5 months after the latest leap
	    ;##- second we know about, if either of these is later than
	    ;##- our default expiration time.
	    if { $bultime > 0 } {
		set expir1 [expr {$bultime+4*30*24*3600}]
		if { $expir1 > $expiration } { set expiration $expir1 }
	    }
	    set expir1 [expr {[lindex $sysLeaps end]+5*30*24*3600}]
	    if { $expir1 > $expiration } { set expiration $expir1 }

	} else {
	    ;##- See whether this leap second is already in our list.
	    ;##- If not, append it.
	    if { $time > [lindex $sysLeaps end] } {
		if {$debug} { puts "  Appending this leap date to our list" }
		lappend gpsLeaps [expr {$time-$gpsEpoch+[llength $sysLeaps]}]
		lappend sysLeaps $time
	    } else {
		if {$debug} { puts "  This leap date is already in our list" }
	    }

	    ;##- We were able to parse Bulletin C, so set the expiration
	    ;##- time for our leap-second info to be just after the next
	    ;##- Bulletin C is scheduled to come out (~6 months from now),
	    ;##- OR just after the final leap second takes effect, OR one
	    ;##- week from now, whichever is latest.
	    if { $bultime > 0 } {
		set expir1 [expr {$bultime+7*30*24*3600}]
		if { $expir1 > $expiration } { set expiration $expir1 }
	    }
	    set expir1 [expr {$time+7*30*24*3600}]
	    if { $expir1 > $expiration } { set expiration $expir1 }
	    set expir1 [expr {[clock seconds]+7*24*3600}]
	    if { $expir1 > $expiration } { set expiration $expir1 }

	}

    } else {
	if {$debug} {puts "Unable to find leap date"}
    }

    return $expiration
}


##=========================================================================
## Name: ScanEnglishOrFrenchDate
##

proc tconvert::ScanEnglishOrFrenchDate { date } {

    if { ! ( [regexp {(\d{1,2})\s+(\w+)\s+((?:19|20)\d\d)} $date \
	    match mdate month year] \
	    || [regexp {((?:19|20)\d\d)\s+(\w+)\s+(\d{1,2})} $date \
	    match year month mdate] ) } {
	return -code error "Cannot parse components of date"
    }

    ;##- Handle either English or French month names
    switch -glob [string tolower $month] {
	jan* { set month jan }
	fe* {set month feb}
	mar* {set month mar}
	ap* - av* {set month apr}
	may - mai {set month may}
	jun* - juin {set month jun}
	jul* - juil* {set month jul}
	aug* - ao* {set month aug}
	sep* {set month sep}
	oct* {set month oct}
	nov* {set month nov}
	d* {set month d}
	default {
	    return -code error \
		    "Cannot understand month as either English or French"
	}
    }

    if { [catch {clock scan "$month $mdate, $year" -gmt 1} systime] } {
	return -code error "Error while scanning modified date string"
    }

    return $systime
}


#==============================================================================
# The following commands are executed immediately when this file is sourced,
# i.e. when you do 'package require tconvert'

namespace import tconvert::*

