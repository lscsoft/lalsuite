#!/bin/sh
d=${srcdir:-.}
cmd="./lalapps_ring_pipe -f S2_H1_dag.ini -l . -d -r"
echo $cmd
eval $cmd || exit $?

exit 0
