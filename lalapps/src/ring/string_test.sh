#!/bin/sh
d=${srcdir:-.}
cmd="./stringfind -v -d 1 -f $d/\*.F -r $d/response.asc -i $d/filterpar.in -o events.out"
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
