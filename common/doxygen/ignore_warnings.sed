# remove blank lines
/^\$$/d

# -- ignore some Doxygen warnings due to Doxygen bugs

# https://bugzilla.gnome.org/show_bug.cgi?id=742151
/^citelist/d

# https://bugzilla.gnome.org/show_bug.cgi?id=743604
/warning: Duplicate anchor/d
/warning: multiple use of section label/d

# https://bugzilla.gnome.org/show_bug.cgi?id=743605
/warning: explicit link request/d

# -- ignore problematic lines

\|lalinference/bayespputils.py:7015: warning: documented symbol 'ligolw LIGOLWContentHandler lalinference::bayespputils::endElement' was not declared or defined|d
\|lalpulsar/NstarTools.py:143: warning: '\\code' command is not allowed in section title, ending section title|d
