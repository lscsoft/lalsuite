# Remove initial copyright block
1 {

  # Start loop
  : loop1

  # Skip commented lines
  \|^//| {
    n
    b loop1
  }

}

# For any non-comment lines
\|^ *[^ /]| {

  # Start Doxygen \verbatim block
  i /// \\code

  # Start loop
  : loop2

  # Skip empty lines
  \|^$| {
    n
    b loop2
  }

  # End Doxygen \code block at next comment line
  \|^ *//| {
    i /// \\endcode
    i ///
    b end2
  }

  # Comment out and print line in \code block
  s|^|/// |
  p

  # Read next line and continue loop
  n
  b loop2

  # End loop
  : end2
}

# Remove Emacs local variables
\|^// Local Variables:$| {

  # Start loop
  : loop3

  # Skip rest of file
  n
  b loop3

}

# Remove any leading spaces from comments
s|^ *//|//|

# Escape percent signs
s|%|\\%|g

# Comment any empty lines
s|^$|///|

# Ensure line is a Doxygen comment
s|^///*|///|

# Print line
p
