#!/bin/sh
d=${srcdir:-.}
cmd="./lalapps_ring --verbose --debug-level=1 --sample-rate=4096 --filter-segsz=8192 --filter-speclen=0 --filter-flow=0.001 --filter-fhighpass=-1 --filter-fmin=150 --filter-fmax=151 --filter-qmin=10 --filter-qmax=11 --filter-maxmm=0.1 --filter-thresh=4 --bank-end-template=0 --output-format=ascii --test-gaussian-data --test-white-spectrum --test-zero-data --test-inject=3.0,150,10,8,0 --strain-data --write-data --write-filter-output --write-format=ascii"
echo $cmd
eval $cmd || exit $?

exit 0
