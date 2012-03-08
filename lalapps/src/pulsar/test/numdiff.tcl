#!/bin/sh
#
#
# Started 04/23/08 by Greg Mendell
#
# Compares two ascii files line by line, splitting each line and
# comparing each item on a line.  If items are numeric, then the
# difference within a specified fractional or absolute tolerance
# is checked.
#
# \
exec tclsh "$0" "$@"

set FALSE 0;
set TRUE  1;

proc Exit {code} {
     exit $code;
}

proc printDifference {lineItems1 lineItems2 lineNumber1 lineNumber2 filename1 filename2 lineCount} {
     puts stderr "Line $lineNumber1 of $filename1 and line $lineNumber2 of $filename2 differ:"
     puts stderr "< $lineItems1";
     puts stderr "> $lineItems2";
     puts stderr " ";
}

proc PrintHelp {} {
    global argv0;
    puts " ";
    puts "Compares two ascii files line by line, splitting each line and"
    puts "comparing each item on a line. If items are numeric, then the"
    puts "difference within a specified fractional or absolute tolerance"
    puts "is checked."
    puts " ";
    puts "Usage:";
    puts " ";
    puts "$argv0 <filename1> <filename2> \[-iL <numLines>\] \[-iC <charString>\] \[-d <delimiter>\} \[-a|-f|-n <tolerance>\]";
    puts " ";
    puts " <filename1>           Required; name of first file.";
    puts " <filename2>           Required; name of second file.";
    puts " -iL <ignoreNumLines>  Optional; ignore first numLines. (Default: 0.)";
    puts " -iC <ignoreString>    Optional; ignore lines starting with this string. (Default: do not ignore any lines.)";
    puts " -d <delimiter>        Optional; use this delimiter to split each line. (Default: split line using spaces.";
    puts " -a|-f|-n <tolerance>  Optional; give absolute (-a) or fractional (-f) or number of decimal places (-n) toleranance to use comparing numbers.";
    puts "                       (Default: Numbers must agree exactly.)";
    puts "                       For example, if x1 is number if filename1 and x2 is the corresponding number in filename2, then";
    puts " ";
    puts "                       -a 0.001 means abs(x1 - x2) < 0.001";
    puts "                       -f 0.001 means 2.0*abs(x1 - x2)/abs(x1 + x2) < 0.001";
    puts "                       -n 1 means compare x1 and x2 digit by digit and only the last digits can differ.";
    puts " ";
    Exit 1;
}

#######################
# MAIN CODE STARTS HERE 
#######################

if { $argc < 2 || $argv == "-h" || $argv == "--help"} {
   PrintHelp;
}

# get the command line arguments
set filename1 [lindex $argv 0];
set filename2 [lindex $argv 1];
set ignoreNumLines 0;
set ignoreString "";
set delimiter " ";
set tolerance 0;
set numericComparisonType "";

set optList [split "-iL:-iC:-d:-a:-f:-n" ":"];
set opt "unknown";
foreach element [lrange $argv 2 end] {
        set optInd [lsearch -exact $optList $element];
        if {$optInd > -1} {
                if {$element == "-iL"} {
                   set opt "-iL";
                } elseif {$element == "-iC"} {
                   set opt "-iC";
                } elseif {$element == "-d"} {
                   set opt "-d";
                } elseif {($element == "-a") || ($element == "-f") || ($element == "-n")} {
                   set opt "tolerance";
                   set numericComparisonType $element;
                } else {
                   set opt $element; # Option used with next element below
                }
        } else {
                if {$opt == "unknown"} {
                   puts stderr "Unknown option: $element";
                   Exit 2;
                } elseif {$opt == "-iL"} {
                        set ignoreNumLines $element;
                } elseif {$opt == "-iC"} {
                        set ignoreString $element;
                } elseif {$opt == "-d"} {
                        set delimiter $element;
                } elseif {$opt == "tolerance"} {
                        set tolerance $element;
                }
                set opt "unknown";
        }
}

# Uncomment to DEBUG:
#puts " ";
#puts "filename1 = $filename1";
#puts "filename2 = $filename2";
#puts "ignoreNumLines = $ignoreNumLines";
#puts "ignoreString = $ignoreString";
#puts "delimiter = $delimiter";
#puts "tolerance = $tolerance";
#puts " ";

if {$ignoreString > ""} {
  set ignoreStringLen [string length $ignoreString];
  set ignoreStringLen [expr $ignoreStringLen - 1];
} else {
  set ignoreStringLen -1;
}

# Read in contents of the first file.
if {[file exists $filename1]} {
    # Read the file into buffer.
    set fid [open $filename1 "r"];
    set buffer [read $fid];
    close $fid;
} else {
    puts stderr "Error: filename1 = $filename1 does not exist.";
    Exit 2;
}
set lineList1 [split $buffer "\n"];

# Read in contents of the second file.
if {[file exists $filename2]} {
    # Read the file into buffer.
    set fid [open $filename2 "r"];
    set buffer [read $fid];
    close $fid;
} else {
    puts stderr "Error: filename2 = $filename2 does not exist.";
    Exit 2;
}
set lineList2 [split $buffer "\n"];

# Make a list of lines to compare for each file:
set fileLines1 [list ];
set lineNumberList1 [list ];
set lineNumber1 0;
foreach line [lrange $lineList1 $ignoreNumLines end] {
      incr lineNumber1;
      if {$ignoreString > ""} {
         if {[string range $line 0 $ignoreStringLen] == "$ignoreString"} {
                 # skip
         } else {
           lappend fileLines1 $line;
           lappend lineNumberList1 $lineNumber1;
         }
      } else {
         lappend fileLines1 $line;
         lappend lineNumberList1 $lineNumber1;
      }

}
set fileLines2 [list ];
set lineNumberList2 [list ];
set lineNumber2 0;
foreach line [lrange $lineList2 $ignoreNumLines end] {
      incr lineNumber2;
      if {$ignoreString > ""} {
         if {[string range $line 0 $ignoreStringLen] == "$ignoreString"} {
                 # skip
         } else {
           lappend fileLines2 $line;
           lappend lineNumberList2 $lineNumber2;
         }
      } else {
         lappend fileLines2 $line;
         lappend lineNumberList2 $lineNumber2;
      }
}

if {[llength $fileLines1] != [llength $fileLines2]} {
   puts stderr "Files $filename1 and $filename2 have a different number of lines to compare! Exiting..."
   Exit 3;
}

set exitCode 0;
set lineCount 0;
# Go through the files line by line
foreach line1 $fileLines1 {
   set line2 [lindex $fileLines2 $lineCount]; # Get the corresponding line from the second file.
   set lineNumber1 [lindex $lineNumberList1 $lineCount]; # line number from original file
   set lineNumber2 [lindex $lineNumberList2 $lineCount]; # line number from original file
   incr lineCount;  # Gives which line we are on
   # Split the lines into items to compare
   set lineItems1 [split $line1 $delimiter];
   set lineItems2 [split $line2 $delimiter];
   # Uncomment to DEBUG:
   #puts "Comparing line $lineCount:";
   #puts "< $lineItems1";
   #puts "> $lineItems2";
   if {[llength $lineItems1] != [llength $lineItems2]} {
      printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
   }  else {
      set itemCount 0;
      foreach item1 $lineItems1 {
           set item2 [lindex $lineItems2 $itemCount];
           # Uncomment to DEBUG:
           #puts "item1 = $item1";
           #puts "item2 = $item2";
           #puts " ";
           incr itemCount;
           if {$tolerance > 0} {
             # Check if data is numerical or string
             if {([string is double $item1]) && ([string trimleft $item1] > "")} {
                set item1IsNumber $TRUE;
             } else {
                set item1IsNumber $FALSE;
             }
             if {([string is double $item2]) && ([string trimleft $item2] > "")} {
                set item2IsNumber $TRUE;
             } else {
                set item2IsNumber $FALSE;
             }
             # Uncomment to DEBUG:
             #puts "item1IsNumber = $item1IsNumber";
             #puts "item2IsNumber = $item2IsNumber";
             #puts " ";
             # Check that items are of the same type
             if {$item1IsNumber != $item2IsNumber} {
                printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                break;
             } else {
               # Handle items depending on whether numerical or not
               if {$item1IsNumber} {
                   if {$numericComparisonType == "-a"} {
                     # Compare absolute difference
                     if { [expr {abs($item1 - $item2) > $tolerance}] } {
                        printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                        break;
                     }
                   } elseif {$numericComparisonType == "-f"} {
                     # Compare fractional difference
                     set meanValue [expr {0.5*(abs($item1) + abs($item2))}];
                     # if meanValue is 0 then both numbers are zero and equal to each other.
                     if {$meanValue > 0} {
                       if { [expr {abs($item1 - $item2)/$meanValue > $tolerance}] } {
                          printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                          break;
                       }
                     } 
                   } elseif {$numericComparisonType == "-n"} {
                     # Compare, ignoring last $tolerance digits.
                     set len1 [string length $item1];
                     set len1 [expr $len1 - 1 - $tolerance];
                     set len2 [string length $item2];
                     set len2 [expr $len2 - 1 - $tolerance];
                     if {$len1 != $len2} {
                        printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                        break;
                     } else {
                       if {$len1 > -1} {
                          set item1 [string range $item1 0 $len1];
                          set item2 [string range $item2 0 $len2];
                          if {$item1 != $item2} {
                             printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                             break;
                          }
                       }
                     }
                   } else {
                     # Should not be able to happen, but just for safety:
                     puts stderr "ERROR: tolerance given but with invalid option: $numericComparisonType. This must be -a, -f, or -n. Exiting...";
                     Exit 3;
                   }
                   # END if {$numericComparisonType == "-a"} elseif...
               } else {
                 if {$item1 != $item2} {
                   printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                   break;
                 }
               }
               # END if {$item1IsNumber} else
             }
             # END if {$item1IsNumber != $item2IsNumber} else
           } else {
             if {$item1 != $item2} {
                printDifference $lineItems1 $lineItems2 $lineNumber1 $lineNumber2 $filename1 $filename2 $lineCount; set exitCode 3;
                break;
             }
           }
           # END if {$tolerance > 0} else
      }
      # END foreach item1 $lineItems1
   }
   # END if {[llength $lineItems1] != [llength $lineItems2]}
}
# END foreach line1 $fileLines1

Exit $exitCode;
