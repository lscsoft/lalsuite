#!/bin/sh
if test -r a.pgm
then
  TrackSearchTest
  exitcode=$?
else
  cp $srcdir/a.pgm .
  TrackSearchTest
  exitcode=$?
  rm a.pgm
fi
exit $exitcode
