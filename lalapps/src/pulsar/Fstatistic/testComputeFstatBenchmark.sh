## functions

function assert {
    if eval test "$@"; then
        echo "OK: 'test $@' passed"
    else
        echo "ERROR: 'test $@' failed"
        exit 1
    fi
}

## common arguments to lalapps_ComputeFstatBenchmark

Nseg=50
common_args="--numSegments=${Nseg} --Tseg=86400 --numFreqBins=1000"

## run lalapps_ComputeFstatBenchmark

cmd="lalapps_ComputeFstatBenchmark ${common_args} --FstatMethod=DemodBest --outputInfo=demod.txt"
echo "=== $cmd ==="
eval $cmd
echo "--- $cmd ---"
echo

cmd="lalapps_ComputeFstatBenchmark ${common_args} --FstatMethod=ResampBest --outputInfo=resamp.txt"
echo "=== $cmd ==="
eval $cmd
echo "--- $cmd ---"
echo

## check for expected output

for file in demod.txt resamp.txt; do

    lines=`sed -n '/^%/d;/./p' "${file}" | wc -l | sed 's|[^0-9]||g'`
    assert "X${lines}" = "X${Nseg}"

    assert -f "${file}.pars"

done
