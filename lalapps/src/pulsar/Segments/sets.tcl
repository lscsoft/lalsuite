#!/bin/sh
#
#
# First written by Philip Charlton, 2003.
#
# Implementation of various functions relating to sets and intervals on
# real line, used in creatRDS and elsewhere.
#
# An "interval" is a Tcl list { a b } reprenting the semi-open
# interval [a, b).
#
# A "set" is a list of "intervals". Note that a set need not be normalised
# ie. it's possible that some intervals overlap, or intervals are not ordered
# etc.
#
# \
exec tclsh "$0" "$@"

#
# Max function
#
proc max { a b } {
  expr $a > $b ? $a : $b
}

#
# Min function
#
proc min { a b } {
  expr $a < $b ? $a : $b
}

#
# Return the union of intervals s and t
#
proc intervalUnion { s t } {
  set result {}

  foreach { a b } $s { continue }
  foreach { x y } $t { continue }
     
  if { $y < $a } {
    lappend result [ list $x $y ]
    lappend result [ list $a $b ]
  } elseif { $x > $b } {
    lappend result [ list $a $b ]
    lappend result [ list $x $y ]
  } else {
    set first [ min $a $x ]
    set last [ max $b $y ]
    lappend result [ list $first $last ]
  }

  return $result
}

#
# Return the intersection of intervals s and t
#
proc intervalIntersection { s t } {
  foreach { a b } $s { continue }
  foreach { x y } $t { continue }
     
  set first [ max $a $x ]
  set last [ min $b $y ]

  if { $first < $last } {
    set result [ list $first $last ]
  } else {
    set result {}
  }

  return $result
}

#
# Return the difference of intervals s and t ie. s\t
#
# Can return the empty set, a single interval or two intervals
#
# I don't really know an efficient way to do this...
#
proc intervalDifference { s t } {
  foreach { a b } $s { continue }
  foreach { x y } $t { continue }
  
  if { $x <= $a } {
    if { $y <= $a } {
      set result [ list [ list $a $b ] ]
    } elseif { $y < $b } {
      set result [ list [ list $y $b ] ]
    } else {
      # y >= b
      set result {}
    }
  } elseif { $x < $b } {
    set result [ list [ list $a $x ] ]
    if { $y < $b } {
      lappend result [ list $y $b ]
    }
  } else {
    # $x >= $b
    set result [ list [ list $a $b ] ]
  }

  return $result
}

#
# Return the union of set S and and interval t
#
# The algorithm is: for each element of s of S, take the
# intersection of s and t. If this is empty, add s to the output.
# If it is non-empty, replace t by the union of s and t and continue.
# After S is exhausted, add the last t to the output and re-sort.
# 
proc setIntervalUnion { S t } {
  set result {}

  foreach s $S {
    set sut [ intervalUnion $s $t ]
    if { [ llength $sut ] == 2 } {
      foreach { a b } $s { continue }
      lappend result [ list $a $b ]
    } else {
      set t [ lindex $sut 0 ]
    }
  }

  foreach { a b } $t { continue }
  lappend result [ list $a $b ]
  set result [ lsort -integer -index 0 $result ]

  return $result
}

proc setUnion { S T } {
  set result $S

  foreach t $T {
    set result [ setIntervalUnion $result $t ]
  }

  return $result
}

#
# Return the intersection of a set S and interval t
#
proc setIntervalIntersection { S t } {
  set result {}

  foreach s $S {
    set st [ intervalIntersection $s $t ]
    if { [ llength $st ] != 0 } {
      lappend result $st
    }
  }

  return $result
}

#
# Return the intersection of sets S and T
#
proc setIntersection { S T } {
  set result {}

  foreach s $S {
    foreach t $T {
      set st [ intervalIntersection $s $t ]
      if { [ llength $st ] != 0 } {
        lappend result $st
      }
    }
  }

  return $result
}

#
# Return the difference of a set S and and interval t ie. S\t
#
proc setIntervalDifference { S t } {

  set result {}
  foreach s $S {
    set st [ intervalDifference $s $t ]
    foreach interval $st {
      lappend result $interval
    }
  }

  return $result
}

#
# Return the difference of sets S and T ie. S\T
#
proc setDifference { S T } {
  set result $S

  foreach t $T {
    set result [ setIntervalDifference $result $t ]
  }

  return $result
}

#
# Return the total size of a set, in seconds
#
proc setSize { S } {
  set size 0

  foreach s $S {
    foreach { a b } $s { continue }
    set size [ expr $size + $b - $a ]
  }

  return $size
}
