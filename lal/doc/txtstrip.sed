# strip blank lines
/^$/d

# get rid of the #! line
/^#!/d

# meta-comment: don't include in documentation
/^###/d

# verbatim text:
/^#verbatim/,/^#\/verbatim/{
/^#verbatim/s/.*//
/^#\/verbatim/s/.*//
s/^ *\(.*\)/	\1/
s/|| fail//
}

# ignore text:
/^#ignore/,/^#\/ignore/d

# irrelevant markups
s/^##verse//
s/^##\/verse//

# comment lines are normal text:
/^#/s/^##* *\(.*\)/\1/
