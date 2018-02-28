#!/bin/sh

BASE="ligolw_sqlite_test"

#
# can ligolw_sqlite detect a missing input file?
#

echo
echo "ligolw_sqlite test 1:  detect missing input"
echo "--------------------------------------------------------------------"
rm -vf ${BASE}.sqlite
ligolw_sqlite --verbose --extract /dev/null --database ${BASE}.sqlite && exit 1
[ -f ${BASE}.sqlite ] && exit 1
echo
echo "ligolw_sqlite test 1:  success"
echo "ligolw_sqlite detected missing input and did not create an empty database"

#
# do ligolw_add and ligolw_sqlite produce identical output when used to
# merge two .xml files into a single .xml document?  confirm that URLs can
# be used to identify local source files
#

echo
echo "ligolw_sqlite test 2:  merge .xml.gz files and compare to ligolw_add"
echo "--------------------------------------------------------------------"
rm -vf ${BASE}_ref.xml ${BASE}.sqlite ${BASE}_output.xml
ligolw_add --verbose --output ${BASE}_ref.xml ${BASE}_input.xml.gz ${BASE}_input.xml.gz
ligolw_sqlite --verbose --replace --database ${BASE}.sqlite --extract ${BASE}_output.xml file://${PWD}/${BASE}_input.xml.gz ${BASE}_input.xml.gz
cmp ${BASE}_ref.xml ${BASE}_output.xml || exit
rm -vf ${BASE}.sqlite ${BASE}_output.xml
echo
echo "ligolw_sqlite test 2:  success"
echo "ligolw_add and ligolw_sqlite produced identical merged documents"

#
# does ligolw_sqlite produce the same document as above if the input
# documents are .sqlite files (converted from the original .xml input).
# confirm that URLs can be used to identify local source files
#

echo
echo "ligolw_sqlite test 3:  merge .sqlite files and compare to ligolw_add"
echo "--------------------------------------------------------------------"
ligolw_sqlite --verbose --preserve-ids --replace --database ${BASE}_input.sqlite ${BASE}_input.xml.gz
ligolw_sqlite --verbose --replace --database ${BASE}.sqlite --extract ${BASE}_output.xml file://${PWD}/${BASE}_input.sqlite ${BASE}_input.sqlite
cmp ${BASE}_ref.xml ${BASE}_output.xml || exit
rm -vf ${BASE}_ref.xml ${BASE}_input.sqlite ${BASE}.sqlite ${BASE}_output.xml
echo
echo "ligolw_sqlite test 3:  success"
echo "ligolw_add and ligolw_sqlite produced identical merged documents"
