# end of document
$a\
\\endinput

# get rid of the #! line
/^#!/d

# meta-comment: don't include in documentation
/^###/d

# tex line:
/^##/{
/^##verse/,/^##\/verse/{
/^##verse/c\
\\begin{verse}
/^##\/verse$/c\
\\end{verse}
s/##*\(.*\)/\1\\\\ /
}
s/##/\\noindent/
s/"\([^"]*\)"/{\\small\\color{blue}\\verb|\1|}/g
}

# verbatim text:
/^#verbatim/,/^#\/verbatim/{
/^#verbatim/c\
\\begin{color}{blue}\
\\small\
\\begin{alltt}
/^#\/verbatim/c\
\\end{alltt}\
\\end{color}
s/|| fail//
# because we're using alltt, need to escape \'s 
s/\\$/\\verb{\\{/
s/#.*/  {\\normalfont\\color{black}\\itshape{}&}/
}

# ignore text:
/^#ignore/,/^#\/ignore/d

# comment line:
/^#/{
s/^#/\\noindent/
s/"\([^"]*\)"/{\\small\\color{blue}\\verb|\1|}/g
s/.*/{\\itshape &}/
}
