#!/bin/sh

prog="./lalapps_ring"
args="\
  --verbose \
  --debug-level=3 \
  --data-cache=H-600000000-600000180.cache \
  --calibration-cache=H1-CAL-600000000-600006000.cache \
  --gps-start-time=600000000 --gps-start-time-ns=0 \
  --gps-end-time=600000070 --gps-end-time-ns=0 \
  --channel=H1:LSC-AS_Q \
  --random-seed=11 \
  --sample-rate=1024 \
  --segment-duration=4 \
  --cutoff-frequency=100 \
  --bank-min-frequency=150 \
  --bank-max-frequency=400 \
  --bank-min-quality=10 \
  --bank-max-quality=20 \
  --bank-max-mismatch=0.03 \
  --threshold=5 \
  --write-raw-data \
  --write-data \
  --write-response \
  --write-spectrum \
  --write-inv-spectrum \
  --write-bank \
  --write-segment \
  --only-segment-numbers=18\
  --write-filter-output \
  --only-template-numbers=0 \
  --inject-file=HL-INJECTIONS_1-600000000-180.xml \
  "

if test ${1:-"none"} = "-n" ; then 
  echo $prog $args
  exit 0
fi

echo $prog $args
$prog $args
