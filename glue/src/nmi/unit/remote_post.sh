#!/bin/sh

# wrap this whole script in a block combining stdout & stderr
{

# exit immediately if any command exits with a non-zero status.
set -e
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

# if anything failed, pack up everything for debugging
if [[ "$_NMI_STEP_FAILED" != "" ]]; then
    DEBUG=.
fi

# put any files we want Metronome to keep for us into results.tar.gz
find . -type f \( -name '*.xml*' -or -name \*.txt \) -print0 | xargs -0 tar zcvf results.tar.gz $DEBUG

# end of stdout/stderr-combining block
} 2>&1
