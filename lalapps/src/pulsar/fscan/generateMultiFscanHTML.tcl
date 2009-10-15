#!/bin/sh
#
# generateMultiFscanHTML.tcl, started by Greg Mendell April 09,2009
#
# Generates html frame view to navigate fscan html pages by date and channel name.
#
# \
exec tclsh "$0" "$@"
# Don't change these lines - they are executed only by bash, not Tcl \
if [ -z "${LIGOTOOLS}" ]; then \
  echo "Shell variable LIGOTOOLS undefined - please \"eval \`<ligotools path>/bin/use_ligotools\`\" and try again"; \
  exit 1; \
fi; \
exec tclshexe "$0" "$@"


# Base names of html output files (written as .tmp and moved to .html).
set ::outputFscanNavigationHTMLFile "fscanNavigation";
set ::outputFscanDateHTMLFile "fscanDates";
set ::outputFscanTimeHTMLFile "fscanTimes";
set ::outputFscanChannelHTMLFile "fscanChannels";
set ::outputFscanOutputHTMLFile "fscanOutput";

proc PrintHelpAndExit {} {
   global argv0;
   puts "";
   puts "Usage:";
   puts "";
   puts "$argv0 <fscanDirectory> \[channelSortingList\]";
   puts "";
   puts "Runs glob */*/fscan*.html from the fscanDirectory and generates html frame view to navigate fscan html pages by date and channel name."
   puts "";
   puts "The fscanDirectory should be a relative path from the working directory where this script is run.";
   puts "";
   puts "The channelSortingList is an optional comma delimited list of strings used to sort channels names into catagories, e.g., \"DARM,MAG,MIC,V\"."
   puts "";
   exit;
}

# Set month names and max days each month (will correct Feb if leap year below).
set ::monthNamesAndMaxDays(01) "January 31";set ::monthNamesAndMaxDays(02) "February 28";set ::monthNamesAndMaxDays(03) "March 31";
set ::monthNamesAndMaxDays(04) "April 30";set ::monthNamesAndMaxDays(05) "May 31";set ::monthNamesAndMaxDays(06) "June 30";
set ::monthNamesAndMaxDays(07) "July 31";set ::monthNamesAndMaxDays(08) "August 31";set ::monthNamesAndMaxDays(09) "September 30";
set ::monthNamesAndMaxDays(10) "October 31";set ::monthNamesAndMaxDays(11) "November 30";set ::monthNamesAndMaxDays(12) "December 31";

# Procedure for writing a calendar for a given year and month, with links to fscans that occur during that month.
proc makeCalender {thisFid thisYearAndMonth} {

  # thisFid is the file id of the file in which to write
  # thisYearAndMonth is in the format fscans_yyyy_mm

  # Find the year, month number, month name, and max days this month
  set thisYearAndMonthList [split $thisYearAndMonth "_"];
  set thisYear [lindex $thisYearAndMonthList 1];
  set thisMonthNum [lindex $thisYearAndMonthList 2];
  set thisMonthNamesAndMaxDaysList [split $::monthNamesAndMaxDays($thisMonthNum)];
  set thisMonthName [lindex $thisMonthNamesAndMaxDaysList 0];
  set maxDays [lindex $thisMonthNamesAndMaxDaysList 1];
  if {($thisMonthNum == "02") && ([expr thisYear % 4] == 0)} {
     set maxDays 29; # Its Februrary of a leap year.
  }

  # Find the column on the calendar on which the 1st occurs:
  set dayCharList [list "Sun" "Mon" "Tue" "Wed" "Thu" "Fri" "Sat"];
  set thisDayNum [lindex $::dayNumList($thisYearAndMonth) 0];
  set thisDayChar [lindex $::dayCharList($thisYearAndMonth) 0];
  set thisDayIndex [lsearch -exact $dayCharList $thisDayChar];
  set thisDayNumTmp [string trimleft $thisDayNum "0"];
  set thisDayOffset [expr $thisDayNumTmp % 7];
  set thisDayOffset [expr $thisDayOffset - 1];
  set firstDayColOffset [expr $thisDayIndex - $thisDayOffset];
  if {$firstDayColOffset < 0} {
     set firstDayColOffset [expr 7 + $firstDayColOffset]; 
  }

  # label the top of the calendar and the days of the week
  puts $thisFid "$thisMonthName $thisYear"
  puts $thisFid "<table style=\"text-align: left;\" border=\"1\" cellpadding=\"1\" cellspacing=\"1\">"
  puts $thisFid "<tbody>"
  puts $thisFid "    <tr>"
  puts $thisFid "      <th>Su</th>"
  puts $thisFid "      <th>Mo</th>"
  puts $thisFid "      <th>Tu</th>"
  puts $thisFid "      <th>We</th>"
  puts $thisFid "      <th>Th</th>"
  puts $thisFid "      <th>Fr</th>"
  puts $thisFid "      <th>Sa</th>"
  puts $thisFid "    </tr>"

  # Start the first row of numbers:
  puts $thisFid "    <tr>"
  set thisColOffset 0;
  while {$thisColOffset < $firstDayColOffset} {
    puts $thisFid "      <td></td>"
    incr thisColOffset
  }

  # Loop through the days of the month
  for {set dayNum 1} {$dayNum <= $maxDays} {incr dayNum} {

    # Get thisDate in the format fscans_yyyy_mm_dd
    # Get thisDayNum in the format dd
    set thisDate $thisYearAndMonth;
    append thisDate "_";
    if {$dayNum < 10} {
       append thisDate "0";
       append thisDate $dayNum;
       set thisDayNum "0";
       append thisDayNum $dayNum;
    } else {
       append thisDate $dayNum;
       set thisDayNum $dayNum;
    }

    # Check if fscans exist for this day in this year and month:
    set thisDayIndex [lsearch -exact $::dayNumList($thisYearAndMonth) $thisDayNum];

    # If fscans exist for this date, then link to the times or channels for this date, else just give the day number.
    if {$thisDayIndex > -1} {
      set thisDateAndTime [lindex $::dirList($thisYearAndMonth) $thisDayIndex]
      if {[llength $::timeList($thisDate)] > 1} {
        puts $thisFid "      <td bgcolor=lightgreen><a href=\"$::fscanDirectory/$thisDateAndTime/$::outputFscanTimeHTMLFile.html\" target=\"fscanChannels\">$dayNum</a></td>";
      } else {
        puts $thisFid "      <td bgcolor=lightgreen><a href=\"$::fscanDirectory/$thisDateAndTime/$::outputFscanChannelHTMLFile.html\" target=\"fscanChannels\">$dayNum</a></td>";
      }
    } else {
       puts $thisFid "      <td>$dayNum</td>";
    }
    incr thisColOffset
    # Check if we need to start a new column
    if {$thisColOffset == 7} {
        puts $thisFid "    </tr>"
        # Start the first row of numbers:
        puts $thisFid "    <tr>"
        set thisColOffset 0;
    }
  }

  # Finish off last row and the table
  while {$thisColOffset < 7} {
    puts $thisFid "      <td></td>"
    incr thisColOffset
  }
  puts $thisFid "    </tr>"
  puts $thisFid "</tbody>"
  puts $thisFid "</table>"
  puts $thisFid " <br>"

}

#### MAIN code starts here

# Set default values
set ::fscanDirectory "none";

# Comma delimited strings used to sort channels names into catagories, e.g., "DARM,MAG,MIC,V"
set ::channelSortingList "none";

# save the working directory
set workingDir [pwd];
#puts $workingDir;

if {$argc < 1} {
   PrintHelpAndExit;
}

# Get the fscan parent directory
set ::fscanDirectory [lindex $argv 0];
if {$::fscanDirectory == ""} {
   puts "Error: fscanDirectory not given.";
   PrintHelpAndExit;
}


if {$argc >= 2} {
   set ::channelSortingList [lindex $argv 1];
   set ::channelSortingList [split $::channelSortingList ","];
}

# Change to fscanDirectory 
cd $::fscanDirectory;

# Get a list of date/channel/fscan.html files
set fullPathDirList [glob "*/*/fscan*.html"];
set fullPathDirList [lsort $fullPathDirList];
#puts $fullPathDirList;

# For each of full paths, make a list of dates and times that fscan exist (each entry has the format fscans_yyyy_mm_dd_hh_mm_ss_tmz_day)
# For each YearAndMonth, construct list of numeric and charactor days that fscans exist and set a flag if fscans exist for more than one time on that date.
# For each Date construct a list paths to fscans that exist for that date.
# For each DateAndTime construct a list of channels and fscan html files.
# We end up with these globally available lists:
# ::dateAndTimeList
# ::chanList($thisDateAndTime)
# ::fscanHTMLList($thisDateAndTime)
# ::timeList($thisDate)
# ::dayNumList($thisYearAndMonth)
# ::dayCharList($thisYearAndMonth)
# ::dirList($thisYearAndMonth)
set ::dateAndTimeList "";  
set lastYearAndMonth "";
set lastDateAndTime "";
set lastDate "";
foreach dirPath $fullPathDirList {
        set pathList [split $dirPath "/"];
        set thisDateAndTime [lindex $pathList 0];
        set thisDateAndTimeList [split $thisDateAndTime "_"];
        set thisYearAndMonth [lrange $thisDateAndTimeList 0 2]; # Get the year and month and day (format fscans_yyyy_mm)
        set thisYearAndMonth [join $thisYearAndMonth "_"];
        set thisDate [lrange $thisDateAndTimeList 0 3];         # Get the year month and day (format fscans_yyyy_mm_dd)
        set thisDate [join $thisDate "_"];
        set thisTime [lrange $thisDateAndTimeList 4 7];         # Get time of day (format hh_mm_ss_tmz)
        set thisTime [join $thisTime "_"];
        set thisDayNum [lindex $thisDateAndTimeList 3];         # Get the numerical day (01, 02, 03, ....31)
        set thisDayChar [lindex $thisDateAndTimeList 8];        # Get the character day (Sun, Mon, Tue, Wed, Thu, Fri, Sat)
        set thisChannel [lindex $pathList 1];
        set thisFscanHTML [lindex $pathList 2];
        if {$thisDateAndTime != $lastDateAndTime} {
           lappend ::dateAndTimeList $thisDateAndTime;
           # Starting a new list of channels and fscan html files for this date
           set ::chanList($thisDateAndTime) $thisChannel;
           set ::fscanHTMLList($thisDateAndTime) $thisFscanHTML;
           # Update list of times that fscans occur for a given date.
           if {$thisDate != $lastDate} {
              # Start a new list of paths for each time fscans exist for this date
              set ::timeList($thisDate) $thisTime;
              set ::dirList($thisDate) $thisDateAndTime;
              if {$thisYearAndMonth != $lastYearAndMonth} {
                # Start new lists of days and directories for which fscans exist this year and month
                set ::dayNumList($thisYearAndMonth) $thisDayNum;
                set ::dayCharList($thisYearAndMonth) $thisDayChar;
                set ::dirList($thisYearAndMonth) $thisDateAndTime;
              } else {
                # Append to the existing lists
                lappend ::dayNumList($thisYearAndMonth) $thisDayNum;
                lappend ::dayCharList($thisYearAndMonth) $thisDayChar;
                lappend ::dirList($thisYearAndMonth) $thisDateAndTime;
              }
           } else {
              # Append to the existing lists
              lappend ::timeList($thisDate) $thisTime;
              lappend ::dirList($thisDate) $thisDateAndTime;
           }
        } else {
           # Append to the existing lists
           lappend ::chanList($thisDateAndTime) $thisChannel;
           lappend ::fscanHTMLList($thisDateAndTime) $thisFscanHTML;
        }
        set lastYearAndMonth $thisYearAndMonth;
        set lastDateAndTime $thisDateAndTime;
        set lastDate $thisDate;
}

cd $workingDir

# Open files
set fidNav [open "$::outputFscanNavigationHTMLFile.tmp" "w"]
set fidDate [open "$::outputFscanDateHTMLFile.tmp" "w"]
set fidOut [open "$::outputFscanOutputHTMLFile.tmp" "w"]

set thisDateAndTime [lindex $::dateAndTimeList 0]

puts $fidNav "<html>"
puts $fidNav "<head>"
puts $fidNav "  <meta content=\"text/html; charset=ISO-8859-1\" http-equiv=\"content-type\">"
puts $fidNav "  <title>Fscans ($::fscanDirectory)</title>"
puts $fidNav "</head>"
puts $fidNav "<FRAMESET COLS=\"23%,77%\" rows=\"100%\">"
puts $fidNav "  <FRAME src=\"fscanDates.html\" name=\"fscanDates\">"
puts $fidNav "  <FRAMESET COLS=\"23%,77%\" rows=\"100%\">"
puts $fidNav "    <FRAME src=\"$::fscanDirectory/$thisDateAndTime/$::outputFscanChannelHTMLFile.html\" name=\"fscanChannels\">"
puts $fidNav "    <FRAME src=\"fscanOutput.html\" name=\"fscanOutput\">"
puts $fidNav "  </FRAMESET>"
puts $fidNav "</FRAMESET>"
puts $fidNav "<body>"
puts $fidNav "</body>"
puts $fidNav "</html>"
close $fidNav

puts $fidDate "<center>"
puts $fidDate "<h1>Fscans ($::fscanDirectory)</h1><br>"
puts $fidDate "<h2>Select a date below and then a channel (or time and then channel) from the next column.<h2><br>"

puts $fidOut "Fscans"
close $fidOut;

set lastYearAndMonth "";
set lastDateAndTime "";
set lastDate "";

#puts $::dateAndTimeList
foreach thisDateAndTime $::dateAndTimeList {

        #puts $thisDateAndTime

        set thisDateAndTimeList [split $thisDateAndTime "_"];

        set thisYearAndMonth [lrange $thisDateAndTimeList 0 2]; # Get the year and month and day (format fscans_yyyy_mm)
        set thisYearAndMonth [join $thisYearAndMonth "_"];
        set thisDate [lrange $thisDateAndTimeList 0 3];         # Get the year month and day (format fscans_yyyy_mm_dd)
        set thisDate [join $thisDate "_"];
        set thisTime [lrange $thisDateAndTimeList 4 7];         # Get time of day (format hh_mm_ss_tmz)
        set thisTime [join $thisTime "_"];

        # If we have not yet made the calendar for this month, then make it:
        if {$thisYearAndMonth != $lastYearAndMonth} {
           makeCalender $fidDate $thisYearAndMonth;
        }

        # Change to the directory for this date and time
        cd $::fscanDirectory/$thisDateAndTime

        # If more than one set of fscans exist for this date, then make links to the channel lists for each time.
        if {[llength $::timeList($thisDate)] > 1} {
           set fidTime [open "$::outputFscanTimeHTMLFile.tmp" "w"]
           puts $fidTime "<center>"
           puts $fidTime "<h2>Times</h2><br>"
           puts $fidTime "$thisDateAndTime <br>"
           puts $fidTime " <br>"
           set timeInd 0;
           foreach thisTime $::timeList($thisDate) {
               set thisDateAndTimeDir [lindex $::dirList($thisDate) $timeInd]
               puts $fidTime "<a href=\"../$thisDateAndTimeDir/$::outputFscanChannelHTMLFile.html\" target=\"fscanChannels\">$thisTime</a> <br>";
               puts $fidTime " <br>"
               incr timeInd
           }
           puts $fidTime "</center>"
           close $fidTime;
           exec mv "$::outputFscanTimeHTMLFile.tmp" "$::outputFscanTimeHTMLFile.html"
        }

        set fidChan [open "$::outputFscanChannelHTMLFile.tmp" "w"]
        puts $fidChan "<center>"
        puts $fidChan "<h2>Channels</h2><br>"
        puts $fidChan "$thisDateAndTime <br>"
        puts $fidChan " <br>"
# Sort channels by comma delimited strings given in $::channelSortingList.
if {$::channelSortingList != "none"} {

    # Initialize a set of temporary lists for each channel type
    foreach chanSortString $::channelSortingList {
       set newTmpChanList($chanSortString)   ""
       set newFscanHTMLList($chanSortString) ""
    }
    set newTmpChanList(Other)   ""
    set newFscanHTMLList(Other) ""

    # Go through the list of channels, adding them to the appropriate list based on the string that matches
    set chanInd 0;
    foreach thisChan $::chanList($thisDateAndTime) {
        set thisFscanHTML [lindex $::fscanHTMLList($thisDateAndTime) $chanInd];
        set noMatch 1;
        foreach chanSortString $::channelSortingList {
           if {[string match "*$chanSortString*" $thisChan]} {
             lappend newTmpChanList($chanSortString) $thisChan;
             lappend newFscanHTMLList($chanSortString) $thisFscanHTML;
             set noMatch 0;
             break;
           }
        }
        if {$noMatch} {
             lappend newTmpChanList(Other) $thisChan;
             lappend newFscanHTMLList(Other) $thisFscanHTML;
        }
        incr chanInd;
    }

    # Finally, output the channels
    foreach chanSortString $::channelSortingList {
        if {[llength $newTmpChanList($chanSortString)] > 0} {
          puts $fidChan " <br>"
          puts $fidChan "$chanSortString Channels: <br>"
          #puts $fidChan " <br>"
          set chanInd 0;
          foreach thisChan $newTmpChanList($chanSortString) {
              set thisFscanHTML [lindex $newFscanHTMLList($chanSortString) $chanInd];
              puts $fidChan "<a href=\"$thisChan/$thisFscanHTML\" target=\"fscanOutput\">$thisChan</a> <br>";
              puts $fidChan " <br>"
              incr chanInd;
          }
        }
    }
    if {[llength $newTmpChanList(Other)] > 0} {
      puts $fidChan " <br>"
      puts $fidChan "Other Channels: <br>"
      #puts $fidChan " <br>"
      set chanInd 0;
      foreach thisChan $newTmpChanList(Other) {
          set thisFscanHTML [lindex $newFscanHTMLList(Other) $chanInd];
          puts $fidChan "<a href=\"$thisChan/$thisFscanHTML\" target=\"fscanOutput\">$thisChan</a> <br>";
          puts $fidChan " <br>"
          incr chanInd;
      }
    }
} else {

        set chanInd 0;
        foreach thisChan $::chanList($thisDateAndTime) {
            set thisFscanHTML [lindex $::fscanHTMLList($thisDateAndTime) $chanInd];
            puts $fidChan "<a href=\"$thisChan/$thisFscanHTML\" target=\"fscanOutput\">$thisChan</a> <br>";
            puts $fidChan " <br>"
            incr chanInd;
        }
}
        puts $fidChan "</center>"
        close $fidChan;
        exec mv "$::outputFscanChannelHTMLFile.tmp" "$::outputFscanChannelHTMLFile.html"
        cd $workingDir
        #puts $::chanList($thisDateAndTime)
        #puts $::fscanHTMLList($thisDateAndTime)
        #puts ""

        set lastYearAndMonth $thisYearAndMonth;
        set lastDateAndTime $thisDateAndTime;
        set lastDate $thisDate;
}

# Close outputFscanDateHTMLFile file:
puts $fidDate "</center>"
close $fidDate

exec mv "$::outputFscanNavigationHTMLFile.tmp" "$::outputFscanNavigationHTMLFile.html"
exec mv "$::outputFscanDateHTMLFile.tmp" "$::outputFscanDateHTMLFile.html"
exec mv "$::outputFscanOutputHTMLFile.tmp" "$::outputFscanOutputHTMLFile.html"

puts "Program $argv0 has finished.";

exit 0
