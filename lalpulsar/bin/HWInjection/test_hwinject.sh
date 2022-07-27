if test "${LALFRAME_ENABLED}" = false; then
    echo "Skipping test: requires LALFrame"
    exit 77
fi

# create injection input files
mfd_args="--hardwareTDD --actuation=./H1PCALXactuationfunction_withDelay.txt --actuationScale=1.0 --Tsft=20 --duration=5"
for inj in 0 1 2; do
    echo "lalpulsar_Makefakedata_v4 @Pulsar${inj}_StrainAmp.cfg --logfile=pulsar${inj}.log ${mfd_args}" > "in.${inj}"
done
echo "===== $0: injection input files ====="
grep . in.*
echo "----- $0: injection input files -----"
echo

# hardware injection base command
hwinj_base_cmd="lalpulsar_hwinject -I H1 -n 3 -G 1341748395"

# generate hardware injection text output
hwinj_cmd="${hwinj_base_cmd} -T -X"
echo "===== $0: run ${hwinj_cmd} ====="
eval "${hwinj_cmd}" > out || true   # lalpulsar_hwinject may not exit nicely
outlines=`cat out | wc -l`
if test ${outlines} -lt 1000; then
    echo "ERROR: got only ${outlines} lines of output"
    exit 1
else
    echo "OK: got ${outlines} lines of output"
fi
echo "----- $0: run ${hwinj_cmd} -----"
echo

# output MFDv4 logs
echo "===== $0: output MFDv4 logs ====="
grep . pulsar*.log
echo "----- $0: output MFDv4 logs -----"
echo

# compare text output to reference result
echo "===== $0: compare text output to reference result ====="
paste out out.ref | awk '
BEGIN {
   exitcode = 0
   maxrelerr = 0
}
FNR == 1 {
   if ($1 != 1234.5 || $2 != 1234.5) {
      printf "ERROR: line %5i: incorrect magic numbers %g, %g\n", FNR, $1, $2
      exitcode = 1
   }
}
NF == 4 {
   if ($1 != $3) {
      printf "ERROR: line %5i: mismatched GPS times %g, %g\n", FNR, $1, $3
      exitcode = 1
   }
   relerr = sqrt( ($2 - $4)^2 )
   if (relerr > 2e-4) {
      printf "ERROR: line %5i: %0.8f does not compare to reference %0.8f (relerr=%g)\n", FNR, $2, $4, relerr
      exitcode = 1
   }
   if (maxrelerr < relerr) {
      maxrelerr = relerr
   }
}
END {
   if (exitcode != 0) {
      printf "ERROR: compare text output to reference result (maxrelerr=%g)\n", maxrelerr
   } else {
      print "OK: compare text output to reference result"
   }
   exit(exitcode)
}'
echo "----- $0: compare text output to reference result -----"
echo

# generate hardware injection frame output
hwinj_cmd="${hwinj_base_cmd} -F 0"
echo "===== $0: run ${hwinj_cmd} ====="
eval "${hwinj_cmd}" >/dev/null || true   # lalpulsar_hwinject may not exit nicely
gwf_files=
for n in `seq 395 1 400`; do
    gwf_file="CW_Injection-1341748${n}-1.gwf"
    if test -f "${gwf_file}"; then
        gwf_files="${gwf_files} ${gwf_file}"
    fi
done
if test "x${gwf_files}" = x; then
    echo "ERROR: did not get any frame files"
    exit 1
else
    echo "OK: got frame files:${gwf_files}"
fi
( echo 1234.5; lalfr-print ${gwf_files} | grep '^1341748' | sort -n -k1,1 ) > out
outlines=`cat out | wc -l`
if test ${outlines} -lt 1000; then
    echo "ERROR: got only ${outlines} lines of output"
    exit 1
else
    echo "OK: got ${outlines} lines of output"
fi
echo "----- $0: run ${hwinj_cmd} -----"
echo

# compare frame output to reference result
echo "===== $0: compare frame output to reference result ====="
paste out out.ref | awk '
BEGIN {
   exitcode = 0
   maxrelerr = 0
}
FNR == 1 {
   if ($1 != 1234.5 || $2 != 1234.5) {
      printf "ERROR: line %5i: incorrect magic numbers %g, %g\n", FNR, $1, $2
      exitcode = 1
   }
}
NF == 4 {
   if ($1 != $3) {
      printf "ERROR: line %5i: mismatched GPS times %g, %g\n", FNR, $1, $3
      exitcode = 1
   }
   relerr = sqrt( ($2 - $4)^2 )
   if (relerr > 2e-4) {
      printf "ERROR: line %5i: %0.8f does not compare to reference %0.8f (relerr=%g)\n", FNR, $2, $4, relerr
      exitcode = 1
   }
   if (maxrelerr < relerr) {
      maxrelerr = relerr
   }
}
END {
   if (exitcode != 0) {
      printf "ERROR: compare frame output to reference result (maxrelerr=%g)\n", maxrelerr
   } else {
      print "OK: compare frame output to reference result"
   }
   exit(exitcode)
}'
echo "----- $0: compare frame output to reference result -----"
