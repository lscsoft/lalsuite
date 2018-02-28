#!/bin/sh
# pack up output for return to metronome

# wrap this whole script in a block combining stdout & stderr
{

# exit immediately if any command exits with a non-zero status.
set -e
# print (unexpanded) shell input lines as they are read
set -v
# print (expanded) commands before they're executed
set -x

# if anything failed, pack up the src/ dir for debugging
if [[ "$_NMI_STEP_FAILED" != "" ]]; then
    DEBUG=src
fi

tar zcvf results.tar.gz opt $DEBUG

# end of stdout/stderr-combining block
} 2>&1
