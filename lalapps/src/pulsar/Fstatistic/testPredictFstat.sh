##---------- names of codes and input/output files
mfd_code="lalapps_Makefakedata_v4"
saf_code="SemiAnalyticF"
pfs_code="lalapps_PredictFstat"

Ftolerance=0.05

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours
maxStartTime=$(echo $startTime $duration | awk '{printf "%.1f", $1 + $2}'); ## convention is [minStartTime,maxStartTime) so maxStartTime==endTime will yield the full set of SFTs
mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
noiseSqrtSh=0.5

cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

IFO=H1

# input parameters for timestamps test
## ---------- data parameters ----------
timestamps="./H1-timestamps.dat"
echo "${startTime}" > ${timestamps}
for i in `seq 1 $(( ${duration} / ${Tsft} - 1 ))`
do
  echo "$((${startTime} + $i * ${Tsft}))" >> ${timestamps}
done

echo
echo "----- Test parameters: ----------"
echo "startTime=${startTime}; duration=${duration}; Tsft=${Tsft}"
echo "h0=${h0}; cosi=${cosi}; psi=${psi}; phi0=${phi0}"
echo "Alpha=${Alpha}; Delta=${Delta}; Freq=${Freq}"
echo "sqrtSn=${noiseSqrtSh}; IFO=${IFO}; mfd_FreqBand=${mfd_FreqBand}"
echo "--------------------------------"

# this part of the command-line is compatible with SemiAnalyticF:
saf_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=. --noiseSqrtSh=$noiseSqrtSh --randSeed=1"

## ---------- Run MFDv4 ----------
cmdline="$mfd_code $mfd_CL";
echo $cmdline
echo -n "Running ${mfd_code} ... "
if ! eval "$cmdline"; then
    echo "FAILED:"
    echo $cmdline
    exit 1
else
    echo "OK."
fi

## ---------- Run SemiAnalyticF ----------
cmdline="$saf_code $saf_CL  --sqrtSh=$noiseSqrtSh"
echo $cmdline
echo -n "Running ${saf_code} ... "
if ! tmp=`eval $cmdline 2> /dev/null`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resSAF=`echo $tmp | awk '{printf "%g", 2.0 * $1}'`

pfs_CL_common=" --Alpha=$Alpha --Delta=$Delta --cosi=$cosi --h0=$h0 --psi=$psi"
## ---------- Run PredictFstat{NoiseWeights} ----------
outfile_pfs1="__tmp_PFS1.dat";
pfs_CL="${pfs_CL_common} --DataFiles=\"./*.sft\" --Freq=$Freq --outputFstat=$outfile_pfs1"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
echo -n "Running ${pfs_code}{NoiseWeights} ... "
if ! tmp=`eval ${cmdline}`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS1=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Run PredictFstat{NoiseWeights,maxStartTime} (testing one-sided constraints) ----------
outfile_pfs1a="__tmp_PFS1a.dat";
pfs_CL="${pfs_CL_common} --DataFiles=\"./*.sft\" --Freq=$Freq --outputFstat=$outfile_pfs1a --maxStartTime=$maxStartTime"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
echo -n "Running ${pfs_code}{NoiseWeights,maxStartTime} ... "
if ! tmp=`eval ${cmdline}`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS1a=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Run PredictFstat{assumeSqrtSX} ----------
outfile_pfs0="__tmp_PFS0.dat";
pfs_CL="${pfs_CL_common} --DataFiles=\"./*.sft\" --outputFstat=$outfile_pfs0 --assumeSqrtSX=${noiseSqrtSh}"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
echo -n "Running ${pfs_code}{assumeSqrtSX} ... "
if ! tmp=`eval $cmdline`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS0=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Run PredictFstat{timestamps,assumeSqrtSX} ----------
outfile_pfs2="__tmp_PFS2.dat";
pfs_CL="${pfs_CL_common} --timestampsFiles=${timestamps} --Tsft=$Tsft --outputFstat=$outfile_pfs2 --assumeSqrtSX=${noiseSqrtSh} --IFOs=$IFO"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
echo -n "Running ${pfs_code}{timestamps,assumeSqrtSX} ... "
if ! tmp=`eval $cmdline`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS2=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Run PredictFstat{minStartTime+duration,assumeSqrtSX} ----------
outfile_pfs3="__tmp_PFS3.dat";
pfs_CL="${pfs_CL_common} --minStartTime=$startTime --duration=$duration --Tsft=$Tsft --outputFstat=$outfile_pfs3 --assumeSqrtSX=${noiseSqrtSh} --IFOs=$IFO"
cmdline="$pfs_code $pfs_CL"
echo $cmdline
echo -n "Running ${pfs_code}{minStartTime+duration,assumeSqrtSX} ... "
if ! tmp=`eval $cmdline`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS3=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Comparing results ----------
echo
echo "SemiAnalyticF:              2F_SA  = $resSAF"
echo "PredictFstat{assumeSqrtSX}: 2F_PF0 = $resPFS0"
echo "PredictFstat{NoiseWeights}: 2F_PF1 = $resPFS1"
echo "PredictFstat{NoiseWeights,maxStartTime}: 2F_PF1a = $resPFS1a"
echo "PredictFstat{timestamps,assumeSqrtSX}: 2F_PF2 = $resPFS2"
echo "PredictFstat{minStartTime+duration,assumeSqrtSX}: 2F_PF2 = $resPFS3"
echo

awk_reldevPercent='{printf "%.4f", 100.0 * sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'
tolerance0=1	## percent
tolerance1=20	## percent
tolerance2=0.0001	## percent

eps0=$(echo $resSAF $resPFS0 | awk "$awk_reldevPercent");
eps1=$(echo $resPFS0 $resPFS1 | awk "$awk_reldevPercent");
eps1a=$(echo $resPFS1 $resPFS1a | awk "$awk_reldevPercent");
eps2=$(echo $resPFS0 $resPFS2 | awk "$awk_reldevPercent");
eps3=$(echo $resPFS0 $resPFS3 | awk "$awk_reldevPercent");

res=0;
fail0=$(echo $eps0 $tolerance0 | awk "$awk_isgtr")
fail1=$(echo $eps1 $tolerance1 | awk "$awk_isgtr")
fail1a=$(echo $eps1a $tolerance2 | awk "$awk_isgtr")
fail2=$(echo $eps2 $tolerance2 | awk "$awk_isgtr")
fail3=$(echo $eps3 $tolerance2 | awk "$awk_isgtr")

echo -n "Relative deviation 2F_PF{assumeSqrtSX} wrt 2F_SA = ${eps0}% (tolerance = ${tolerance0}%)"
if [ "$fail0" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi

echo -n "Relative deviation 2F_PF{NoiseWeights} wrt 2F_PF{assumeSqrtSX} = ${eps1}% (tolerance = ${tolerance1}%)"
if [ "$fail1" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi

echo -n "Relative deviation 2F_PF{NoiseWeights,maxStartTime} wrt 2F_PF{NoiseWeights} = ${eps1a}% (tolerance = ${tolerance2}%)"
if [ "$fail1a" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi

echo -n "Relative deviation 2F_PF{timestamps,assumeSqrtSX} wrt 2F_PF{assumeSqrtSX} = ${eps2}% (tolerance = ${tolerance2}%)"
if [ "$fail2" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi

echo -n "Relative deviation 2F_PF{minStartTime+duration,assumeSqrtSX} wrt 2F_PF{assumeSqrtSX} = ${eps3}% (tolerance = ${tolerance2}%)"
if [ "$fail3" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi
echo

exit $res;
