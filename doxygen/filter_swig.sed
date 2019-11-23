# For any non-comment lines
\|^ *[^ /]| {

  # Start Doxygen \verbatim block
  i /// \\code

  # Start loop
  : loop

  # Skip empty lines
  \|^$| {
    n
    b loop
  }

  # End Doxygen \code block at next comment line
  \|^ *//| {
    i /// \\endcode
    i ///
    b end
  }

  # Comment out and print line in \code block
  s|^|/// |
  p

  # Read next line and continue loop
  n
  b loop

  # End loop
  : end
}

# Remove any leading spaces from comments
s|^ *//|//|

# Escape percent signs
s|%|\\%|g

# Comment any empty lines
s|^$|///|

# Print line
p
