#!/bin/sh
#
# showFscansOverTime.tcl, started by Greg Mendell Feburary 2012
#
# Don't change these lines - they are executed only by bash, not Tcl \
if [ -z "${LIGOTOOLS}" ]; then \
  echo "Shell variable LIGOTOOLS undefined - please \"eval \`<ligotools path>/bin/use_ligotools\`\" and try again"; \
  exit 1; \
fi; \
exec tclshexe "$0" "$@"

proc SendMail { subject msg } {
  if {$::notifyList > ""} {
    if { [ exec uname ] == "SunOS" } {
      set mailer /usr/ucb/mail
      exec $mailer -s $subject $::notifyList << $msg
    } else {
      set mailer mail
      exec echo $msg | $mailer -s $subject $::notifyList
    }
  }
}

proc PrintHelpAndExit {} {
   global argv0;
   puts "";
   puts "Usage:";
   puts "";
   puts "$argv0 <fscan_location> <channel_name> <start_time> <end_time>";
   puts "";
   puts "Finds fscan plots under the given location, for the given channel, and puts them onto an html page, showing the fscans over time.";
   puts "";
   exit;
}


##### SET DEFAULT VALUES AND SOURCE FILES #####

package require tconvert
#source tconvert.tcl;

set startingPWD [exec pwd];

# Set up the numbers of the months
set monthNums(Jan) "01";set monthNums(Feb) "02";set monthNums(Mar) "03";set monthNums(Apr) "04";
set monthNums(May) "05";set monthNums(Jun) "06";set monthNums(Jul) "07";set monthNums(Aug) "08";
set monthNums(Sep) "09";set monthNums(Oct) "10";set monthNums(Nov) "11";set monthNums(Dec) "12";

# Set month names
set ::monthNames(01) "Jan";set ::monthNames(02) "Feb";set ::monthNames(03) "Mar";
set ::monthNames(04) "Apr";set ::monthNames(05) "May";set ::monthNames(06) "Jun";
set ::monthNames(07) "Jul";set ::monthNames(08) "Aug";set ::monthNames(09) "Sep";
set ::monthNames(10) "Oct";set ::monthNames(11) "Nov";set ::monthNames(12) "Dec";

set startTime "none";
set endTime "none";

##### MAIN CODE STARTS HERE #####

# Parse command line arguements

if {$argv == "" || $argv == "-h" || $argv == "?" || $argv == "--help"} {
   PrintHelpAndExit;
}

if {$argc < 4} {
   PrintHelpAndExit;
}

# Get the command line options

set fscanLocation [lindex $argv 0];
set inputChanName [lindex $argv 1];
set chanName [join [split $inputChanName ":"] "_"]; # replace : with _ in the channel name.
set startTime [lindex $argv 2];
set endTime [lindex $argv 3];

if {$startTime == "none"} {
   puts "Could not determine start time from $resourceFile or $startTimeFromFile.";
   PrintHelpAndExit;
}

if {$endTime == "none"} {
   puts "Could not determine end time from $resourceFile.";
   PrintHelpAndExit;
} elseif {$endTime == "now"} {
  set endTime [ tconvert now ];
  set endTime [ expr $endTime - $timeLag];
} elseif { [string index $endTime 0] == "+"} {
  set endTime [string range $endTime 1 end];
  set endTime [ expr $startTime + $endTime];
}

set duration [expr $endTime - $startTime];

set startTime [lindex $argv 2];
set endTime [lindex $argv 3];

set outputFile "fscans_"
append outputFile $chanName;
append outputFile "_";
append outputFile $startTime
append outputFile "_to_";
append outputFile "$endTime.html";

# Format start and end time
set formatStartTime [tconvert -l -f "%a %b %d %I:%M:%S %p %Z %Y" $startTime]
set formatEndTime [tconvert -l -f "%a %b %d %I:%M:%S %p %Z %Y" $endTime]

#set formatTimeNow $formatEndTime;
#set timeYear [lindex $formatTimeNow 6];
#set timeMonthChar [lindex $formatTimeNow 1];
#set timeMonthNum $monthNums($timeMonthChar);
#set timeDayNum [lindex $formatTimeNow 2];
#set timeDayChar [lindex $formatTimeNow 0];
#set timeAMPM [lindex $formatTimeNow 4];
#set timeTime [split [lindex $formatTimeNow 3] ":"]
## Convert to 24 hr time
#if {$timeAMPM == "PM"} {
#   set timeTimeHr [string trimleft [lindex $timeTime 0] "0"];
#   if {$timeTimeHr < 12} {
#      set timeTimeHr [expr $timeTimeHr + 12]; # Only add 12 for 1 to 11 PM, not for 12 PM.
#   }
#   set timeTime [list $timeTimeHr [lindex $timeTime 1] [lindex $timeTime 2]];
#}
#set timeTime [join $timeTime "_"];
#set timeZone [lindex $formatTimeNow 5];
#set joinFormatTimeNow [join [list $timeYear $timeMonthNum $timeDayNum $timeTime $timeZone $timeDayChar] "_"]

############# Begin finding png files under the fscanLocation between the startTime and endTime ######################
set count 0;

# Get a list of directories that contain SFTs
set fullPathDirList [glob -nocomplain "$fscanLocation/*"];
set fullPathDirList [lsort $fullPathDirList];
#puts $fullPathDirList;

set bandList [];

# Check if this location has png files between the start time and end time
foreach dirPath $fullPathDirList {
        set pathList [split $dirPath "/"];
        set thisDateAndTime [lindex $pathList end];
        #puts $thisDateAndTime
        set thisDateAndTimeList [split $thisDateAndTime "_"];
        set thisYear [lindex $thisDateAndTimeList 1];      # Get the yyyy year
        set thisMonthNum [lindex $thisDateAndTimeList 2];  # Get the numerical month (01, 02, 03, ....12)
        set thisMonthName $::monthNames($thisMonthNum);    # Get month name (Jan, Feb, Mar,...)
        set thisDayNum [lindex $thisDateAndTimeList 3];    # Get the numerical day (01, 02, 03, ....31)
        set thisTime [lrange $thisDateAndTimeList 4 6];    # Get time of day
        set thisHour [lindex $thisTime 0];                 # Get hour of the day
        if {$thisHour > 23} {
           # Handle bug where hour can be 24.
           set thisHour [expr $thisHour - 12];
           set thisMin [lindex $thisTime 1];
           set thisSec [lindex $thisTime 2];               # (format hh:mm:ss)
           set thisTime "$thisHour:$thisMin:$thisSec";
        } else {
           set thisTime [join $thisTime ":"];              # (format hh:mm:ss)
        }
        set thisTimeZone [lindex $thisDateAndTimeList 7];  # Get the time zone (PDT, CDT, ...)
        set thisGPSTime [tconvert "$thisMonthName $thisDayNum $thisYear $thisTime $thisTimeZone"]
        #puts "This date and time = $thisMonthName $thisDayNum $thisYear $thisTime $thisTimeZone";
        #puts "This GPS Time = $thisGPSTime";
        # Note using greater than $startTime and not greater than or equal $startTime
        # because directory give names give the end time of the last SFTs in that directory.
        if {($thisGPSTime > $startTime) && ($thisGPSTime <= $endTime)} {
           set pngList [glob -nocomplain "$dirPath/$chanName/*.png"];
           foreach pngFileWithPath $pngList {
                   set pngFile [file tail $pngFileWithPath];
                   if {[string match "*2.png*" $pngFile]} {
                      set band [lindex [split $pngFile "_"] 1]; # get the start frequency of the band
                      if {$count == 0} {
                        lappend bandList $band;
                        set pngSpectraOutList($band) [];
                      }
                      lappend pngSpectraOutList($band) $pngFileWithPath; 
                   } else {
                      set band [lindex [split $pngFile "_"] 1]; # get the start frequency of the band
                      if {$count == 0} {
                        set pngSpectrogramOutList($band) [];
                      }
                      lappend pngSpectrogramOutList($band) $pngFileWithPath;
                   }
           }
           incr count;
        }
}
############# End finding png files under the fscanLocation between the startTime and endTime ######################

set bandList [lsort -real $bandList]

############# Begin output the html file with the spectrograms and spectra #################

set fid [open "$outputFile" "w"];

puts $fid "<html>"
puts $fid "<head>"
puts $fid "<meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">"
puts $fid "<title>FSCANPLOTSOVERTIME</title>"
puts $fid "</head>"
puts $fid "<body>"
puts $fid "<div style=\"text-align: center;\">"
puts $fid "<h1>Fscan plots for $inputChanName from $formatStartTime to $formatEndTime</h1>"
puts $fid "</div>"
puts $fid  "<br>"
puts $fid  "<div style=\"text-align: center;\">"
puts $fid  "</div>"
puts $fid  "<br>" 
puts $fid  "<table style=\"width: 100%; text-align: left;\" border=\"1\" cellpadding=\"2\" cellspacing=\"2\">"
puts $fid  "<tbody>"

foreach band $bandList {

    # For each band output a rwo of spectrgrams and spectra
    set ind 0;
    puts $fid "  <tr>"
    foreach pngSpectrogram $pngSpectrogramOutList($band) {
      set pngSpectra [lindex $pngSpectraOutList($band) $ind];
      incr ind;
      puts $fid "  <td style=\"vertical-align: top;\">"
      #puts $fid "    <a href=\"$pngSpectrogram\"><img alt=\"\" src=\"$pngSpectrogram\" style=\"border: 0px solid ; width: 576px; height: 432px;\"></a><br>"
      #puts $fid "    <a href=\"$pngSpectra\"><img alt=\"\" src=\"$pngSpectra\" style=\"border: 0px solid ; width: 576px; height: 432px;\"></a><br>"
      puts $fid "    <a href=\"$pngSpectrogram\"><img alt=\"\" src=\"$pngSpectrogram\" style=\"border: 0px solid ; width: 150; height: 100;\"></a><br>"
      puts $fid "    <a href=\"$pngSpectra\"><img alt=\"\" src=\"$pngSpectra\" style=\"border: 0px solid ; width: 150; height: 100;\"></a><br>"
      puts $fid "  </td>"
    }
    puts $fid "  </tr>"

}

close $fid
############# End output the html file with the spectrograms and spectra #################

exit;
#### MAIN CODE ENDS HERE #####
