#!/bin/sh
ephemdir=./ephems
sftdir=./somedata

params=" -I 2 -y 03 -E $ephemdir -D $sftdir -f 300.0 -b 0.1 -a 2.2 -z 0.003 -d 0.8 -c 0.003 "

if [ x$1 = x ]; then
    prog="./lalapps_ComputeFStatistic";
else
    prog=$1;
fi

echo "Starting run 1 ..."
$prog $params
echo "run 1 complete"
echo `diff --brief Fstats Fstats.ref`

echo "Starting run 2 ..."
$prog @testrun.cfg
echo "run 2 complete"
echo `diff --brief Fstats Fstats.ref`

