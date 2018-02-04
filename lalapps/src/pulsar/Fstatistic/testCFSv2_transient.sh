#!/bin/bash

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

builddir="./";
injectdir="../Injections/"

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

## ----- allow user-control of hotloop variant to use
if [ -n "$FSTAT_METHOD" ]; then
    FstatMethod="--FstatMethod=${FSTAT_METHOD}"
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
saf_code="${builddir}lalapps_SemiAnalyticF"
## allow user to specify a different CFSv2 version to test by passing as cmdline-argument
if test $# -eq 0 ; then
    cfsv2_code="${builddir}lalapps_ComputeFstatistic_v2"
else
    cfsv2_code="$@"
fi

Dterms=8
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
Tdata=144000		## total data+observation span: 40 hours
tauInj=72000		## injected signal duration: 20 hours
tauInjDays=$(echo $tauInj | awk '{printf "%g", $1/(24*3600) }');
t0Inj=$(echo $startTime $Tdata $tauInj | awk '{printf "%d", $1 + 0.5*$2 - 0.5*$3 }');

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5

h0=1
cosi=-0.3
psi=0.6
phi0=1.5

Freq=100.12345
f1dot=-1e-10;

## mfd-specific bands
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

## cfs search bands
NFreq=100;
cfs_FreqBand=$(echo $Tdata | awk '{printf "%.16g", 1.0 / $1 }');	## fix band to 1/T so we're close to signal peak always
cfs_Freq=$(echo $Freq $cfs_FreqBand | awk '{printf "%.16g", $1 - $2 / 2.0}');
cfs_dFreq=$(echo $cfs_FreqBand $NFreq | awk '{printf "%.16g", $1 / $2 }');
cfs_nCands=$NFreq	## toplist length: keep all cands

cfs_f1dotBand=0;
cfs_f1dot=$(echo $f1dot $cfs_f1dotBand | awk '{printf "%.16g", $1 - $2 / 2.0}');
##Nf1dot=10
cfs_df1dot=1 ##$(echo $cfs_f1dotBand $Nf1dot | awk '{printf "%g", $1 / $2}');

noiseSqrtSh=5

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi

IFO=H1

## ----- define output directory and files
testDir=testCFSv2transient.d
rm -rf $testDir
mkdir -p $testDir
SFTdir=${testDir}

outfile_Fstat1=${testDir}/testCFSv2_run1.dat
outfile_Loudest1=${testDir}/Fstat_loudest_run1.dat
outfile_transient1=${testDir}/testCFSv2_tCW_run1.dat
outfile_transientMap1=${testDir}/testCFSv2_tCW_Fstatmap_run1.dat
outfile_Fstat2=${testDir}/testCFSv2_run2.dat
outfile_Loudest2=${testDir}/Fstat_loudest_run2.dat
outfile_transient2=${testDir}/testCFSv2_tCW_run2.dat
outfile_transientMap2=${testDir}/testCFSv2_tCW_Fstatmap_run2.dat
outfile_transientMap3=${testDir}/testCFSv2_tCW_Fstatmap_run3.dat
outfile_atoms=${testDir}/testCFSv2_atoms

## awk commands for results comparisons
awk_iseq='{if($1==$2) {print "1"}}'
awk_isgtr='{if($1>$2) {print "1"}}'
awk_absdev='{printf "%.2e", sqrt(($1-$2)^2) }'
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'

## Tolerance of comparison
Tolerance=5e-2	## 5%
ToleranceT=$( echo $Tsft | awk '{printf "%g", 1.5*$1 }')
Tolerance2Fratio=0.15

##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Signal"
echo "----------------------------------------------------------------------"
echo

# this part of the command-line is compatible with SemiAnalyticF:
base_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --Tsft=$Tsft --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0"
# concatenate this with the SAF-specific switches:
saf_CL="${base_CL} --startTime=$t0Inj --duration=$tauInj"
# concatenate this with the mfd-specific switches:
mfd_CL="${base_CL} --startTime=$startTime --duration=$Tdata --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir/$IFO-sfts.sft --f1dot=$f1dot --outSingleSFT  --refTime=$startTime --transientWindowType=rect --transientStartTime=$t0Inj --transientTauDays=$tauInjDays"
if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL --randSeed=1"
echo $cmdline;
if ! eval "$cmdline 2> /dev/null"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi


echo
echo "----------------------------------------------------------------------"
echo "STEP 2: Comparing CW and tCW results over a freq grid at (t0,tau)=(T0,Tdata): "
echo "----------------------------------------------------------------------"
echo

cfs_CL_base="--IFO=$IFO --Alpha=$Alpha --Delta=$Delta --DataFiles='$SFTdir/*.sft' --Dterms=${Dterms} ${FstatMethod} --refTime=$startTime --TwoFthreshold=0" #  --NumCandidatesToKeep=${cfs_nCands}
if [ "$haveNoise" != "true" ]; then
    cfs_CL_base="$cfs_CL_base --SignalOnly"
fi

cfs_searchBand="--Freq=$cfs_Freq --FreqBand=$cfs_FreqBand --dFreq=$cfs_dFreq --f1dot=$cfs_f1dot --f1dotBand=$cfs_f1dotBand --df1dot=$cfs_df1dot"

cfs_CL_run1="${cfs_CL_base} ${cfs_searchBand} --outputFstat=${outfile_Fstat1} --outputLoudest=${outfile_Loudest1} --outputTransientStats=${outfile_transient1} --outputTransientStatsAll=${outfile_transientMap1}" # --outputFstatAtoms=${outfile_atoms}"

cmdline="$cfsv2_code $cfs_CL_run1"
echo $cmdline;
if ! eval "$cmdline 2> /dev/null"; then
    echo "Error.. something failed when running '$cfsv2_code' ..."
    exit 1;
fi
echo

## work around toplist-sorting bugs in CFSv2: manually sort before comparing
sort -o ${outfile_Fstat1} ${outfile_Fstat1}

## check for matching filelengths of standard-CW, tCW and tCW-Fstatmap output
## (should be over same frequency vector)
echo "--------- Comparing file lengths between F-stat and tCW output files ---------"
Nlines_Fstat=$(grep -v ^% ${outfile_Fstat1} | wc -l)
Nlines_trans=$(grep -v ^% ${outfile_transient1} | wc -l)
Nlines_transMap=$(grep -v ^% ${outfile_transientMap1} | wc -l)
if [[ $Nlines_Fstat -ne $Nlines_trans || $Nlines_Fstat -ne $Nlines_transMap ]]; then
    echo "==> ERROR: file lengths differ: len(outputFstat)=${Nlines_Fstat}, len(--outputTransientStats)=${Nlines_trans}, len(--outputTransientStatsAll=${Nlines_transMap}."
    exit 2
else
    echo "==> OK."
fi
echo

## check for matching 2F<->max2F when using full-length transient window

topline_Fout=$(grep -v ^% $outfile_Fstat1 | sort -nr -k7,7  | head -1)
# echo 'best Fstat result:'
# echo $topline_Fout
topline_tCWout=$(grep -v ^% $outfile_transient1 | sort -nr -k9,9 | head -1)
# echo 'best tCW result by 2Fmax:'
# echo $topline_tCWout
topline_tCWout_Bstat=$(grep -v ^% $outfile_transient1 | sort -nr -k10,10 | head -1)
# echo 'best tCW result by Bstat:'
# echo $topline_tCWout_Bstat
topline_tCWmap=$(grep -v ^% $outfile_transientMap1 | sort -nr -k9,9 | head -1)

Fout_top2F_freq=$(echo $topline_Fout | awk '{print $1}')
Fout_top2F_2F=$(echo $topline_Fout | awk '{print $7}')
tCWout_top2F_freq=$(echo $topline_tCWout | awk '{print $1}')
tCWout_top2F_max2F=$(echo $topline_tCWout | awk '{print $9}')
tCWout_top2F_Bstat=$(echo $topline_tCWout | awk '{print $10}')
tCWout_topB_freq=$(echo $topline_tCWout_Bstat | awk '{print $1}')
tCWout_topB_max2F=$(echo $topline_tCWout_Bstat | awk '{print $9}')
tCWout_topB_Bstat=$(echo $topline_tCWout_Bstat | awk '{print $10}')
tCWmap_top2F_freq=$(echo $topline_tCWmap | awk '{print $1}')
tCWmap_top2F_max2F=$(echo $topline_tCWmap | awk '{print $9}')

reldev_freq_Fout_tCWout_top2F=$(echo $Fout_top2F_freq $tCWout_top2F_freq | awk "$awk_reldev")
reldev_freq_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_freq $tCWout_topB_freq | awk "$awk_reldev")
reldev_2F_Fout_tCWout_top2F=$(echo $Fout_top2F_2F $tCWout_top2F_max2F | awk "$awk_reldev")
reldev_2F_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_max2F $tCWout_topB_max2F | awk "$awk_reldev")
reldev_Bstat_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_Bstat $tCWout_topB_Bstat | awk "$awk_reldev")
reldev_freq_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_freq $tCWmap_top2F_freq | awk "$awk_reldev")
reldev_2F_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_max2F $tCWmap_top2F_max2F | awk "$awk_reldev")

fail1=$(echo $reldev_freq_Fout_tCWout_top2F $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_freq_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
fail3=$(echo $reldev_2F_Fout_tCWout_top2F $Tolerance | awk "$awk_isgtr")
fail4=$(echo $reldev_2F_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
fail5=$(echo $reldev_Bstat_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
fail6=$(echo $reldev_freq_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")
fail7=$(echo $reldev_2F_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")

echo "--------- Comparing results for (t0,tau)=(T0,Tdata): [Tolerance = ${Tolerance}] ---------"
echo "                                        freq                 2F         Bstat"
echo "==>  outputFstat:                       $Fout_top2F_freq    $Fout_top2F_2F ---"
echo "==>  --outputTransientStats(top2F):     $tCWout_top2F_freq $tCWout_top2F_max2F  $tCWout_top2F_Bstat"
echo "==>  --outputTransientStats(topB):      $tCWout_topB_freq $tCWout_topB_max2F  $tCWout_topB_Bstat"
echo "==>  --outputTransientStatsAll(top2F):  $tCWmap_top2F_freq $tCWmap_top2F_max2F ---"

if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" -o "$fail5" -o "$fail6" -o "$fail7" ]; then
    echo "==> *FAILED*"
    exit 2
else
    echo "==> OK"
fi


echo
echo "----------------------------------------------------------------------"
echo " STEP 3: Rerunning with full (t0,tau) search grid: "
echo "----------------------------------------------------------------------"
echo

# currently not using outputTransientStatsAll here,
# as it would be a ~80MB file
# which is not prohibitive, but does slow things down a bit,
# and the above targeted test case should suffice for now

taumin=$(echo $Tsft | awk '{printf "%d", 2*$1 }');
#tauBand=$(echo $Tdata $taumin | awk '{printf "%d", $1-$2 }');
tauBand=$Tdata
t0min=$startTime
t0Band=$(echo $Tdata $taumin | awk '{printf "%d", $1-$2 }');

cfs_CL_run2="${cfs_CL_base} ${cfs_searchBand} --outputFstat=${outfile_Fstat2} --outputLoudest=${outfile_Loudest2} --outputTransientStats=${outfile_transient2} --transient-WindowType=rect --transient-t0Epoch=$t0min --transient-t0Band=$t0Band --transient-dt0=$Tsft --transient-tau=$taumin --transient-tauBand=$tauBand --transient-dtau=$Tsft"
cmdline="$cfsv2_code $cfs_CL_run2"
echo $cmdline;
if ! eval "$cmdline 2> /dev/null"; then
    echo "Error.. something failed when running '$cfsv2_code' ..."
    exit 1;
fi
echo

## work around toplist-sorting bugs in CFSv2: manually sort before comparing
sort -o ${outfile_Fstat2} ${outfile_Fstat2}

echo "--------- Checking that --outputFstat and --outputLoudest produced unchanged results: ---------"
if ! eval "diff -I '[%][%].*' ${outfile_Loudest1} ${outfile_Loudest2}"; then
    echo "Error: --outputLoudest produced different results between runs."
    exit 1
fi
if ! eval "diff -I '[%][%].*' ${outfile_Fstat1} ${outfile_Fstat2}"; then
    echo "Error: --outputFstat produced different results between runs."
    exit 1
fi
echo "==> OK."
echo

## check for matching filelengths again
## (we're only outputting the maxlikelihood (t0,tau) results at each Doppler point)
echo "--------- Comparing file lengths between F-stat and tCW output files ---------"
Nlines_Fstat=$(grep -v ^% ${outfile_Fstat2} | wc -l)
Nlines_trans=$(grep -v ^% ${outfile_transient2} | wc -l)
if [ $Nlines_Fstat -ne $Nlines_trans ]; then
    echo "==> ERROR: file lengths differ: len(outputFstat)=${Nlines_Fstat}, len(--outputTransientStats)=${Nlines_trans}."
    exit 2
else
    echo "==> OK."
fi
echo

## check that max(max2F) is achieved at true (t0,tau) parameters

topline_Fout=$(grep -v ^% $outfile_Fstat2 | sort -nr -k7,7  | head -1)
# echo 'best Fstat result:'
# echo $topline_Fout
topline_tCWout=$(grep -v ^% $outfile_transient2 | sort -nr -k9,9 | head -1)
# echo 'best tCW result by 2Fmax:'
# echo $topline_tCWout
topline_tCWout_Bstat=$(grep -v ^% $outfile_transient2 | sort -nr -k10,10 | head -1)
# echo 'best tCW result by Bstat:'
# echo $topline_tCWout_Bstat

Fout_top2F_freq=$(echo $topline_Fout | awk '{print $1}')
Fout_top2F_2F=$(echo $topline_Fout | awk '{print $7}')
tCWout_top2F_freq=$(echo $topline_tCWout | awk '{print $1}')
tCWout_top2F_max2F=$(echo $topline_tCWout | awk '{print $9}')
tCWout_top2F_Bstat=$(echo $topline_tCWout | awk '{print $10}')
tCWout_topB_freq=$(echo $topline_tCWout_Bstat | awk '{print $1}')
tCWout_topB_max2F=$(echo $topline_tCWout_Bstat | awk '{print $9}')
tCWout_topB_Bstat=$(echo $topline_tCWout_Bstat | awk '{print $10}')
tCWout_top2F_t0ML=$(echo $topline_tCWout | awk '{print $7}')
tCWout_top2F_tauML=$(echo $topline_tCWout | awk '{print $8}')
tCWout_topB_t0MP=$(echo $topline_tCWout_Bstat | awk '{print $11}')
tCWout_topB_tauMP=$(echo $topline_tCWout_Bstat | awk '{print $12}')

reldev_freq_Fout_tCWout_top2F=$(echo $Fout_top2F_freq $tCWout_top2F_freq | awk "$awk_reldev")
reldev_freq_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_freq $tCWout_topB_freq | awk "$awk_reldev")
reldev_2F_Fout_tCWout_top2F=$(echo $Fout_top2F_2F $tCWout_top2F_max2F | awk "$awk_reldev")
reldev_2F_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_max2F $tCWout_topB_max2F | awk "$awk_reldev")
reldev_Bstat_tCWout_top2F_tCWout_topB=$(echo $tCWout_top2F_Bstat $tCWout_topB_Bstat | awk "$awk_reldev")
absdev_t0inj_t0ML=$(echo $t0Inj $tCWout_top2F_t0ML | awk "$awk_absdev")
absdev_t0inj_t0MP=$(echo $t0Inj $tCWout_topB_t0MP | awk "$awk_absdev")
absdev_tauinj_tauML=$(echo $tauInj $tCWout_top2F_tauML | awk "$awk_absdev")
absdev_tauinj_tauMP=$(echo $tauInj $tCWout_topB_tauMP | awk "$awk_absdev")

# fail1=$(echo $reldev_freq_Fout_tCWout_top2F $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_freq_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
#fail3=$(echo $reldev_2F_Fout_tCWout_top2F $Tolerance | awk "$awk_isgtr")
fail4=$(echo $reldev_2F_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
fail5=$(echo $reldev_Bstat_tCWout_top2F_tCWout_topB $Tolerance | awk "$awk_isgtr")
fail6=$(echo $absdev_t0inj_t0ML $ToleranceT | awk "$awk_isgtr")
fail7=$(echo $absdev_t0inj_t0MP $ToleranceT | awk "$awk_isgtr")
fail8=$(echo $absdev_tauinj_tauML $ToleranceT | awk "$awk_isgtr")
fail9=$(echo $absdev_tauinj_tauMP $ToleranceT | awk "$awk_isgtr")

echo "--------- Comparing results for (t0,tau) search grid: [Tolerance = ${Tolerance} resp. ${ToleranceT}] ---------"
echo "                                     freq                 2F         Bstat     t0[ML/MP]      tau[ML/MP]"
echo "==>  outputFstat:                    $Fout_top2F_freq    $Fout_top2F_2F ---       ---            ---"
echo "==>  --outputTransientStats(top2F):  $tCWout_top2F_freq $tCWout_top2F_max2F   $tCWout_top2F_Bstat $tCWout_top2F_t0ML      $tCWout_top2F_tauML"
echo "==>  --outputTransientStats(topB):   $tCWout_topB_freq $tCWout_topB_max2F   $tCWout_topB_Bstat $tCWout_topB_t0MP $tCWout_topB_tauMP"
echo "(injection: t0=$t0Inj, tau=$tauInj)"
if [ "$fail2" -o "$fail4" -o "$fail5" -o "$fail6" -o "$fail7" -o "$fail8" -o "$fail9" ]; then
    echo "==> *FAILED*"
    exit 2
else
    echo "==> OK"
fi
echo

## check for expected 2F->max2F improvement at injection point with true (t0,tau) parameters


echo "--------- Checking if transient max2F is:
                consistent with SemiAnalyticF
                and about a factor 2 higher than full 2F ---------"

echo -n "Running '$saf_code' ... "
echo
cmdline="$saf_code $saf_CL --sqrtSh=$sqrtSh"
echo $cmdline
if ! resF=`eval "$cmdline  2> /dev/null"`; then
    echo "Error ... something failed running '$saf_code' ..."
    exit 1;
fi
saf2F=$(echo $resF | awk '{printf "%g", 2.0 * $1}')

ratio_max2F_2F=$(echo $tCWout_top2F_max2F $Fout_top2F_2F | awk '{print $1/$2}')
ratio_max2F_saf2F=$(echo $tCWout_top2F_max2F $saf2F | awk '{print $1/$2}')
reldev_ratio_max2F_2F=$(echo $ratio_max2F_2F 2 | awk "$awk_reldev")
reldev_ratio_max2F_saf2F=$(echo $ratio_max2F_saf2F 1 | awk "$awk_reldev")
fail1=$(echo $reldev_ratio_max2F_2F $Tolerance2Fratio | awk "$awk_isgtr")
fail2=$(echo $reldev_ratio_max2F_saf2F $Tolerance2Fratio | awk "$awk_isgtr")

echo "transient max2F: $tCWout_top2F_max2F"
echo "full 2F:         $Fout_top2F_2F ($ratio_max2F_2F lower instead of expected factor 2.0)"
echo "SemiAnalyticF:   $saf2F ($ratio_max2F_saf2F different instead of expected equality)"

if [ "$fail1" -o "$fail2" ]; then
    echo "==> *FAILED*"
    exit 2
else
    echo "==> OK"
fi


echo
echo "----------------------------------------------------------------------"
echo " STEP 4: Rerunning with full (t0,tau) search grid at a single freq: "
echo "----------------------------------------------------------------------"
echo

cfs_CL_run3="${cfs_CL_base} --Freq=$Freq --FreqBand=0 --dFreq=0 --f1dot=$f1dot --f1dotBand=0 --df1dot=0 --outputTransientStatsAll=${outfile_transientMap3} --transient-WindowType=rect --transient-t0Epoch=$t0min --transient-t0Band=$t0Band --transient-dt0=$Tsft --transient-tau=$taumin --transient-tauBand=$tauBand --transient-dtau=$Tsft"
cmdline="$cfsv2_code $cfs_CL_run3"
echo $cmdline;
if ! eval "$cmdline 2> /dev/null"; then
    echo "Error.. something failed when running '$cfsv2_code' ..."
    exit 1;
fi
echo

topline_tCWmap=$(grep -v ^% $outfile_transientMap3 | sort -nr -k9,9 | head -1)
tCWmap_top2F_freq=$(echo $topline_tCWmap | awk '{print $1}')
tCWmap_top2F_t0=$(echo $topline_tCWmap | awk '{print $7}')
tCWmap_top2F_tau=$(echo $topline_tCWmap | awk '{print $8}')
tCWmap_top2F_max2F=$(echo $topline_tCWmap | awk '{print $9}')

reldev_freq_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_freq $tCWmap_top2F_freq | awk "$awk_reldev")
reldev_t0_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_t0ML $tCWmap_top2F_t0 | awk "$awk_reldev")
reldev_tau_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_tauML $tCWmap_top2F_tau | awk "$awk_reldev")
reldev_2F_tCWout_top2F_tCWmap=$(echo $tCWout_top2F_max2F $tCWmap_top2F_max2F | awk "$awk_reldev")
fail1=$(echo $reldev_freq_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_t0_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")
fail3=$(echo $reldev_tau_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")
fail4=$(echo $reldev_2F_tCWout_top2F_tCWmap $Tolerance | awk "$awk_isgtr")

echo "--------- Comparing results for (t0,tau)=(T0,Tdata): [Tolerance = ${Tolerance}] ---------"
echo "                                        freq                 2F        t0        tau"
echo "==>  --outputTransientStats(top2F):     $tCWout_top2F_freq $tCWout_top2F_max2F  $tCWout_top2F_t0ML $tCWout_top2F_tauML"
echo "     (from step 3, full grid)"
echo "==>  --outputTransientStatsAll(top2F):  $tCWmap_top2F_freq $tCWmap_top2F_max2F $tCWmap_top2F_t0 $tCWmap_top2F_tau"

if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo "==> *FAILED*"
    exit 2
else
    echo "==> OK"
fi


## -------------------------------------------
## clean up files
## -------------------------------------------
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDir
fi
echo

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old
