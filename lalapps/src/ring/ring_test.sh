#!/bin/sh
d=${srcdir:-.}
cmd="./lalapps_ring --verbose --debug-level=3 --frame-path=$d --frame-files=\*.F --response-file=$d/response.asc --sample-rate=1024 --filter-params=$d/filterpar.in"
echo $cmd
eval $cmd || exit $?
#sed -n '/#/p' events.out \
#    > maxevents.out \
#  || exit $?
#sed '/#/d' events.out \
#    | sort -r -k2 \
#    | sort -k1,1 -t. -s \
#    | uniq -w9 \
#    >> maxevents.out \
#  || exit $?
