#!/bin/sh
d=${srcdir:-.}
cmd="./lalapps_ring --verbose --debug-level=1 --frame-path=. --frame-files=\*.F --response-file=./response.asc --sample-rate=1024 --filter-segsz=8192 --filter-speclen=1024 --filter-flow=100 --filter-fhighpass=90 --filter-fmin=150 --filter-fmax=151 --filter-qmin=10 --filter-qmax=11 --filter-maxmm=0.1 --filter-thresh=4 --filter-scale=1e20 --bank-end-template=0 --output-format=ascii --inject-file=HL-INJECTIONS_1-600000000-180.xml --write-data --write-filter-output --write-format=ascii"
echo $cmd
eval $cmd || exit $?

exit 0

# Note: generate the injection file with (some variant of) this command:
../power/lalapps_binj --gps-start-time=600000000 --gps-end-time=600000180 --time-step=120 --coordinates=ZENITH --quality=10 --hpeak=3e-13 --waveform=Ringdown
