#!/bin/sh

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## test PredictFStat by comparison with SemiAnalyticF
extra_args="$@"

builddir="./";
injectdir="../Injections/"

##---------- names of codes and input/output files
mfd_code="lalapps_Makefakedata_v4"
saf_code="lalapps_SemiAnalyticF"
pfs_code="lalapps_PredictFStat"

mfd_path="${injectdir}${mfd_code}"
saf_path="${builddir}${saf_code}"
pfs_path="${builddir}${pfs_code}"

SFTdir="./testPredictFStat_sfts"

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

Ftolerance=0.05
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
noiseSqrtSh=0.5
h0Double=2  ## same SNR for noiseSqrtSh=1 as assumed in SignalOnly
cosi=-0.3

psi=0.6
phi0=1.5

Freq=100.12345
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');



IFO=H1

echo
echo "----- Test parameters: ----------"
echo "startTime=${startTime}; duration=${duration}; Tsft=${Tsft}"
echo "h0=${h0}; cosi=${cosi}; psi=${psi}; phi0=${phi0}"
echo "Alpha=${Alpha}; Delta=${Delta}; Freq=${Freq}"
echo "sqrtSn=${noiseSqrtSh}; IFO=${IFO}; mfd_FreqBand=${mfd_FreqBand}"
echo "--------------------------------"

if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir
else
    rm -f ${SFTdir}/*
fi

# this part of the command-line is compatible with SemiAnalyticF:
saf_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
# concatenate this with the mfd-specific switches:
mfd_CL="${saf_CL} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --noiseSqrtSh=$noiseSqrtSh --randSeed=1"

## ---------- Run MFDv4 ----------
cmdline="$mfd_path $mfd_CL";
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${mfd_code} ... "
if ! eval "$cmdline 2> /dev/null"; then
    echo "FAILED:"
    echo $cmdline
    exit 1
else
    echo "OK."
fi

## ---------- Run SemiAnalyticF ----------
cmdline="$saf_path $saf_CL  --sqrtSh=$noiseSqrtSh"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${saf_code} ... "
if ! tmp=`eval $cmdline 2> /dev/null`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resSAF=`echo $tmp | awk '{printf "%g", 2.0 * $1}'`

pfs_CL_common=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --cosi=$cosi --psi=$psi --Freq=$Freq --DataFiles=\"${SFTdir}/*.sft\""
## ---------- Run PredictFStat{NoiseWeights} ----------
outfile_pfs1="__tmp_PFS1.dat";
pfs_CL="${pfs_CL_common} --h0=$h0 --outputFstat=$outfile_pfs1"
cmdline="$pfs_path $pfs_CL"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${pfs_code}{NoiseWeights} ... "
if ! tmp=`eval ${cmdline}`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS1=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Run PredictFStat{SignalOnly} ----------
outfile_pfs0="__tmp_PFS0.dat";
pfs_CL="${pfs_CL_common} --h0=$h0Double --outputFstat=$outfile_pfs0 --SignalOnly"
cmdline="$pfs_path $pfs_CL"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${pfs_code}{SignalOnly} ... "
if ! tmp=`eval $cmdline`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi
resPFS0=`echo $tmp | awk '{printf "%g", $1}'`

## ---------- Comparing results ----------
echo
echo "SemiAnalyticF:              2F_SA  = $resSAF"
echo "PredictFStat{SignalOnly}:   2F_PF0 = $resPFS0"
echo "PredictFStat{NoiseWeights}: 2F_PF1 = $resPFS1"
echo

awk_reldevPercent='{printf "%.1f", 100.0 * sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'
tolerance0=1	## percent
tolerance1=20	## percent

eps0=$(echo $resSAF $resPFS0 | awk "$awk_reldevPercent");
eps1=$(echo $resSAF $resPFS1 | awk "$awk_reldevPercent");

res=0;
fail0=$(echo $eps0 $tolerance0 | awk "$awk_isgtr")
fail1=$(echo $eps1 $tolerance1 | awk "$awk_isgtr")

echo -n "Relative deviation 2F_PF{SignalOnly}   wrt 2F_SA = ${eps0}% (tolerance = ${tolerance0}%)"
if [ "$fail0" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi

echo -n "Relative deviation 2F_PF{NoiseWeights} wrt 2F_SA = ${eps1}% (tolerance = ${tolerance1}%)"
if [ "$fail1" ]; then
    echo " ==> FAILED."
    res=1;
else
    echo " ==> OK."
fi
echo

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_pfs0 $outfile_pfs1
fi

exit $res;
