#!/bin/sh

fail() {
  echo "FAIL: $0"
  exit 1
}

ln -sf ${1} || fail

dirname=`pwd`
destdir=`pwd`
frfiles=`ls $dirname/*.gwf`
test -n "$frfiles" || fail

rm -f calibration.cache || fail
for file in $frfiles; do
  basename=`basename $file .gwf`
  IFS_save="$IFS"
  IFS=-
  command='printf "%s %s %s %s "'
  eval $command $basename $out >> calibration.cache || fail
  IFS="$IFS_save"
  command='echo "file://localhost$destdir/$basename.gwf"'
  eval $command $out >> calibration.cache || fail
done
