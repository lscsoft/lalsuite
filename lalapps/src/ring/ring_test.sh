#!/bin/sh
./lalapps_ring -f "*.F" -i filterpar.in -o events.out \
  || exit $?
sed -n '/#/p' events.out \
    > maxevents.out \
  || exit $?
sed '/#/d' events.out \
    | sort -r -k2 \
    | sort -k1,1 -t. -s \
    | uniq -w9 \
    >> maxevents.out \
  || exit $?
