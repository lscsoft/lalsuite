#!/bin/bash

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

builddir="./";
injectdir="../Injections/"

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
cfs_code="${builddir}lalapps_ComputeFStatistic_v2"
pfs_code="${builddir}lalapps_PredictFStat"
SFTdir="./testBinarySFTs"

# ---------- fixed parameter of our test-signal
Tsft=20
startTime=850468893
duration=86400

h0=0.5
cosi=0.5
psi=0.3
phi0=1.3

Alpha=6.275
Delta=1.25
refTime=731163327.0
Freq=200.0
f1dot=1e-09

IFO=H1

## fixed binary parameters
##orbitasini=2.94
orbitasini=4.94036957341473
dopplermax=0.02
##orbitArgp=5.2
orbitArgp=2.88
TpSec=22665
TpNS=1234567
if [ -z "$orbitEcc" ]; then
    orbitEcc=0.35
fi
if [ -z "$orbitPeriod" ]; then
    orbitPeriod=5981.1234
fi

## mfd-specific bands
mfd_FreqBand=20;
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir;
else
    rm -f $SFTdir/*;
fi

echo "--------------------------------------------------"
echo "Injection parameters:"
echo "Alpha=$Alpha, Delta=$Delta"
echo "h0=$h0, cosi=$cosi, psi=$psi, phi0=$phi0"
echo "Freq=$Freq, f1dot=$f1dot"
echo "orbitasini=$orbitasini, orbitPeriod=$orbitPeriod, orbitEcc=$orbitEcc, orbitArgp=$orbitArgp, orbitTpSSBsec=$TpSec, orbitTpSSBnan=$TpNS"
echo "--------------------------------------------------"
echo

## ---------- generate data with binary injection
outfileCFS="CFS.dat";
mfd_CL1="--refTime=${refTime} --Alpha=$Alpha --Delta=$Delta --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --f1dot=$f1dot --orbitasini=$orbitasini --orbitPeriod=$orbitPeriod --orbitEcc=$orbitEcc --orbitArgp=$orbitArgp --orbitTpSSBsec=$TpSec --orbitTpSSBnan=$TpNS"
mfd_CL="${mfd_CL1} --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/catsft.sft --outSingleSFT --IFO=$IFO"
cmdline="$mfd_code $mfd_CL"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running $cmdline ... "
if ! eval "$cmdline &> /dev/null"; then
    echo "lalapps_Makefakedata_v4 failed"
    exit;
fi
echo "done."

## ---------- predict F-stat value 'PFS'
cmdline="${pfs_code} --Alpha=$Alpha --Delta=$Delta --Freq=$Freq --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --IFO=$IFO --DataFiles='$SFTdir/*.sft' --SignalOnly"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running $cmdline ... "
resPFS0=`$cmdline`
resPFS=`echo $resPFS0 | awk '{printf "%g", $1}'`
echo "done:         2F_PFS = $resPFS"

## ---------- targeted CFSv2 search
cfs_CL=" --refTime=${refTime} --Alpha=$Alpha --Delta=$Delta  --Freq=$Freq --f1dot=$f1dot --orbitasini=$orbitasini --orbitPeriod=$orbitPeriod --orbitEcc=$orbitEcc --orbitArgp=$orbitArgp --orbitTpSSBsec=$TpSec --orbitTpSSBnan=$TpNS --DataFiles='$SFTdir/*.sft' --SignalOnly --dopplermax=$dopplermax"
cmdline="$cfs_code $cfs_CL  --outputLoudest=$outfileCFS"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running $cmdline ... "
if ! eval "$cmdline &> /dev/null"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi
resCFS=$(cat ${outfileCFS} | awk '{if($1=="twoF") {printf "%.11g", $3}}')
res_h0=$(cat ${outfileCFS} | awk '{if($1=="h0") {printf "%.11g", $3}}')
res_cosi=$(cat ${outfileCFS} | awk '{if($1=="cosi") {printf "%.11g", $3}}')
res_psi=$(cat ${outfileCFS} | awk '{if($1=="psi") {printf "%.11g", $3}}')
res_phi0=$(cat ${outfileCFS} | awk '{if($1=="phi0") {printf "%.11g", $3}}')
echo "done: 2F_CFS = $resCFS"
echo "Estimated parameters: h0=${res_h0}, cosi=${res_cosi}, psi=${res_psi}, phi0=${res_phi0}"

## ---------- check agreement between 'PFS' and 'CFS' ----------
tolerancePercent=5	## percent
awk_reldevPercent='{printf "%.1f", 100.0 * sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
awk_isgtr='{if($1>$2) {print "1"}}'

reldevPercent=$(echo $resCFS $resPFS | awk "$awk_reldevPercent")
echo -n "Relative 2F deviation = ${reldevPercent}%  (tolerance = ${tolerancePercent}%)"
fail=$(echo $reldevPercent $tolerancePercent | awk "$awk_isgtr")
if [ "$fail" ]; then
    echo " ==> FAILED"
    retstatus=1
else
    echo " ==> OK"
    retstatus=0
fi
echo

## -------------------------------------------
## clean up files
## -------------------------------------------
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfileCFS
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old

exit $retstatus