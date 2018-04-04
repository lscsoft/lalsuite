#!/bin/sh

if [[ $# -ne 2 ]]; then
    echo usage: $0 git_id_a git_id_b
    exit 1
fi

GIT_ID_A=$1
GIT_ID_B=$2

# this script additionally expects:
#
# 1 or more pairs of identically-named xml data products exist in
#   ./$GIT_ID_A/ and ./$GIT_ID_B/
# lalsuite glue is installed in default OS directories

# wrap this whole script in a block combining stdout & stderr
{

# exit immediately if any command exits with a non-zero status.
#set -e
# treat unset variables as an error when performing parameter expansion.
set -u
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

# for debugging
ls -1 */*.xml*

for f in */*.xml*; do
    gunzip -v $f 2> /dev/null
done

RETVAL=0
for f in $GIT_ID_A/*.xml; do
    basef=$(basename $f)
    ./ligolw_diff --exclude-tables=processgroup:process:table --exclude-columns=search_summarygroup:search_summary:lal_cvs_tag --columns $f $GIT_ID_B/$basef > $basef.diff 2>&1
    if [ $? -ne 0 ]; then
	RETVAL=1
    fi
    echo; echo "===== BEGIN XML DIFF OUTPUT of $f ====="; echo
    cat $basef.diff
    echo; echo "===== END XML DIFF OUTPUT of $f ====="; echo
    echo; echo "ligolw_diff returned $RETVAL"

    echo; echo "As a sanity-check, let's try a simple diff too..."
    echo; echo "===== BEGIN SIMPLE DIFF OUTPUT of $f ====="; echo
    diff $f $GIT_ID_B/$basef | tee $basef.simplediff 2>&1
    echo; echo "===== END SIMPLE DIFF OUTPUT of $f ====="; echo
done

# put any files we want Metronome to keep for us into results.tar.gz
#find . -type f \( -name '*.xml*' -or -name \*.txt \) -print0 | xargs -0 tar zcvf results.tar.gz *diff
find . -name head -prune -or \( -name \*.xml -or -name \*.txt \) -print0 | xargs -0 tar zcvf results.tar.gz *diff

# end of stdout/stderr-combining block
} 2>&1

exit $RETVAL
