#!/bin/sh

# This test takes a while... disable it if this is a problem
THETEST=./CLRTest
#THETEST=./NoTest

if test -r CLRindata.asc
then
  $THETEST
  exitcode=$?
else
  cp $srcdir/CLRindata.asc .
  $THETEST
  exitcode=$?
  rm CLRindata.asc
fi
exit $exitcode
