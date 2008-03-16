#!/bin/sh
./lalapps_detresponse -S -NTEST -r0. -d1.57 -o1.023 -Dtest -s2001 -u2048 -i600. -v4 -e7
./lalapps_detresponse -W -NTEST -o1.023 -Dtest -s2001 -u2048 -i600. -v4 -e7
