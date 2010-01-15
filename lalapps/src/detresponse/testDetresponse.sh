#!/bin/sh
./lalapps_detresponse -S -NTEST -r0. -d1.57 -o1.023 -Dtest -s2001 -u2048 -i600. -v0 -e0 || exit 1
./lalapps_detresponse -W -NTEST -o1.023 -Dtest -s2001 -u2048 -i600. -v0 -e0 || exit 1

## clean up files
rm -f TEST_*.txt

## tests passed
exit 0
