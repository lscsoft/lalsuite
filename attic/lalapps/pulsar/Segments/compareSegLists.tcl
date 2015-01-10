#!/bin/sh
#
#
# Started 08/24/06 by Greg Mendell
#
# Computes the intersection, difference, or percent overlap per segment between
# the contents of two files with a list of GPS segments or timestamps
#
# \
exec tclsh "$0" "$@"

proc Exit {code} {
     global argv0;
     if {$code > 1} {
           puts "The $argv0 script finished with an error."
        } elseif {$code == 1} {
           # Displayed help.  Just quit without a message.
        } else {
           # Just exit
        }
        exit $code;
}

proc PrintHelp {} {
    global argv0;
    puts " ";
    puts "Computes the intersection, difference, or percent overlap per segment between";
    puts "the contents of two files with a list of GPS segments or timestamps";    
    puts " ";
    puts "The files need to have a list of either GPS segments of the form \[start_time,end_time) with";
    puts "start_time in the first column and end_time in the second column (all other columns are ignored)";
    puts "or GPS timestamps (i.e., start_times) in the first column (all other columns are ignored). In the";
    puts "latter case a duration must be given such that the end_time = start_time + duration for each timestamp."
    puts " ";
    puts "Usage:";
    puts " ";
    puts "$argv0 -f1 <filename1> -f2 <filename2> \[-t1 <duration1>\] \[-t2 <duration2>\] <-i | -d | -o | -u> \[-s\] \[-T\]";
    puts " ";
    puts " -f1 <filename1>     Required; name of first file with list of segments or timestamps.";
    puts " -f2 <filename2>     Required; name of second file with list of segments or timestamps.";
    puts " -t1 <duration1>     Optional; if given then filename1 contains timestamps as defined above.";
    puts " -t2 <duration2>     Optional; if given then filename2 contains timestamps as defined above.";    
    puts " ";
    puts " -i | -d | -o | -u   Specifies the operation between segments from filename1 and filename2."
    puts "                     One of these is required:";
    puts "                     -i:  find the intersection between the segment lists.";
    puts "                     -d:  find the difference between the segment lists.";
    puts "                     -o:  find the percent overlap on a per segment basis."
    puts "                     -u:  find the union between the segment lists.";
    puts " "
    puts " -s                  Optional; give size of each interval";
    puts " -T                  Optional; give totals";
    puts " ";    
    puts "The output is printed to stdout.";
    puts " ";
    puts "Examples: "
    puts " ";
    puts "Find the intersection between times S4 H1 HW injections were running and S4 H1 SFTs exists:";
    puts " ";
    puts "$argv0 -f1 S4/S4H1_HWInjRunning_try2.txt -f2 S4/SFT_TStamps_S4_LHO.txt -t2 1800 -i > S4H1_Intersect_HWInjRuntry2andSFTs.txt";
    puts " ";
    puts "Find the overlap between each S4 H1 Einstein@Home segment and times S4 H1 HW injections were running and SFTs exist:";
    puts " ";
    puts "$argv0 -f1 S4/S4H1_EatH_segs.txt -f2 S4/S4H1_Intersect_HWInjRuntry2andSFTs.txt -o -s -T";
    puts " ";
    Exit 1;
}

#######################
# MAIN CODE STARTS HERE 
#######################

if { $argc < 1 || $argv == "-h" || $argv == "--help"} {
   PrintHelp;
}

source sets.tcl; # source the procedures for handling segment lists

# get the command line arguments
set operation "unknown";
set filename1 "unknown";
set filename2 "unknown";
set duration1 0;
set duration2 0;
set optList [split "-f1:-f2:-t1:-t2:-i:-d:-o:-u:-s:-T" ":"];
set opt "unknown";
set giveSizes 0;
set giveTotals 0;
foreach element $argv {
        set optInd [lsearch -exact $optList $element];
        if {$optInd > -1} {
                if {$element == "-i"} {
                   set operation "intersection";
                } elseif {$element == "-d"} {
                   set operation "difference";
                } elseif {$element == "-o"} {
                   set operation "overlap";
                } elseif {$element == "-u"} {
                   set operation "union";
                } elseif {$element == "-s"} {
                   set giveSizes 1;
                } elseif {$element == "-T"} {
                   set giveTotals 1;
                } else {
                   set opt $element; # Option used with next element below
                }
        } else {
                if {$opt == "unknown"} {
                   puts "Unknown option: $element";
                   Exit 2;
                } elseif {$opt == "-f1"} {
                        set filename1 $element;
                } elseif {$opt == "-f2"} {
                        set filename2 $element;
                } elseif {$opt == "-t1"} {
                        set duration1 $element;
                } elseif {$opt == "-t2"} {
                        set duration2 $element;
                }
                set opt "unknown";
        }
}

# Validate the arguments and set up related values
if {$operation == "unknown"} {
      puts "Missing -i, -d, -o, or -u option.";
      Exit 2;
}

# Get or generate first list of segments from filename1
if {$filename1 == "unknown"} {
      puts "Missing -f1 <filename1> option.";
      Exit 2;
} else {
   if {[file exists $filename1]} {
       # Read the file into buffer.
       set fid [open $filename1 "r"];
       set buffer [read $fid];
       close $fid;
     } else {
       puts "Error: filename1 = $filename1 does not exist.";
       Exit 2;
     }
     set fileLines1 [list ];   
     foreach line [split $buffer "\n"] {
          set line [string trimleft $line]; #trim leading spaces.
          if {[string index $line 0] == "#" || $line == ""} {
                    # line is blank or a comment, skip
          } else {
            lappend fileLines1 $line;
          }
     }
     set segList1 {};     
     foreach interval $fileLines1 {
          set splitInterval [split $interval];
          set start_time [lindex $splitInterval 0];
          if {$duration1 == 0} {
              set end_time [lindex $splitInterval 1];
          } else {
              set end_time [expr $start_time + $duration1];
          }
          #puts "Interval start_time1 = $start_time";
          #puts "Interval end_time1 = $end_time\n"; 
          set interval [list $start_time $end_time];
          set segList1 [ setIntervalUnion $segList1 $interval ];
     }
} 

# Get or generate second list of segments from filename2
if {$filename2 == "unknown"} {
      puts "Missing -f2 <filename1> option.";
      Exit 2;
} else {
   if {[file exists $filename2]} {
       # Read the file into buffer.
       set fid [open $filename2 "r"];
       set buffer [read $fid];
       close $fid;
     } else {
       puts "Error: filename2 = $filename2 does not exist.";
       Exit 2;
     }
     set fileLines2 [list ];   
     foreach line [split $buffer "\n"] {
          set line [string trimleft $line]; #trim leading spaces.
          if {[string index $line 0] == "#" || $line == ""} {
                    # line is blank or a comment, skip
          } else {
            lappend fileLines2 $line;
          }
     }
     set segList2 {};     
     foreach interval $fileLines2 {
          set splitInterval [split $interval];
          set start_time [lindex $splitInterval 0];
          if {$duration2 == 0} {
              set end_time [lindex $splitInterval 1];
          } else {
              set end_time [expr $start_time + $duration2];
          }
          #puts "Interval start_time2 = $start_time";
          #puts "Interval end_time2 = $end_time\n";
          set interval [list $start_time $end_time];
          set segList2 [ setIntervalUnion $segList2 $interval ];
     }      
} 

set outputSegs {};
if {$operation == "intersection"} {
   set outputSegs [ setIntersection $segList1 $segList2 ]
} elseif {$operation == "difference"} {
   set outputSegs [ setDifference $segList1 $segList2 ]
} elseif {$operation == "union"} {
   set outputSegs [ setUnion $segList1 $segList2 ]
} elseif {$operation == "overlap"} {
  if {$giveSizes > 0} {
     puts "Segments From File 1   Duration   Overlap Duration     Fractional Overlap";     
  } else {
     puts "Segments From File 1      Fractional Overlap With Segments From File 2";
  }
  set totalSegList1 0;
  set totalOverlap 0;
  foreach interval $segList1 {
          #puts "interval = {$interval}";
          set intervalSize [ setSize "{$interval}" ];
          set totalSegList1 [ expr $totalSegList1 + $intervalSize ];
          set overlapSegs [ setIntervalIntersection $segList2 $interval];
          set overlapSize [ setSize $overlapSegs ];
          set totalOverlap [ expr $totalOverlap + $overlapSize ];
          set fracOverlap [ expr 1.0*$overlapSize / (1.0*$intervalSize) ]
          if {$giveSizes > 0} {
             puts "$interval     $intervalSize       $overlapSize                 $fracOverlap";
          } else {
             puts "$interval        $fracOverlap";
          }
  }
  if {$giveTotals > 0} {
      set totalFracOverlap [ expr 1.0*$totalOverlap / (1.0*$totalSegList1) ];
      puts "Total Duration, File 1 Segs = $totalSegList1 s. Total overlap = $totalOverlap s. Total fractional overlap = $totalFracOverlap.";
  }  
}

#puts "segList1 = $segList1";
#puts "segList2 = $segList2";
foreach interval $outputSegs {
   set start_time [lindex $interval 0];
   set end_time [lindex $interval 1];
   if {$giveSizes > 0} {
      set intervalSize [ setSize "{$interval}" ];
      puts "$start_time $end_time $intervalSize";
   } else {
      puts "$start_time $end_time";
   }
} 
if {[llength $outputSegs] > 0} {
   if {$giveTotals > 0} {
      set totalSize [ setSize $outputSegs ];
      puts "Total duration = $totalSize seconds.";
   }
}

Exit 0;
