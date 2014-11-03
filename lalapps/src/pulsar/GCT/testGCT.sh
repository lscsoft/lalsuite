#!/bin/bash

## set LAL debug level
echo "Setting LAL_DEBUG_LEVEL=${LAL_DEBUG_LEVEL:-msglvl1,memdbg}"
export LAL_DEBUG_LEVEL

#NOCLEANUP="1"
#DEBUG=1

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

# The only thing where 'dirsep' can and should be used is in paths of the SFT files,
# as in fact SFTfileIO is the only code that requires it to be set properly. Other
# file references should be handled by the shell (or wine) and converted if necessary.
dirsep="/"

builddir="./";
injectdir="../Injections/"
fdsdir="../FDS_isolated/"


if [ "`echo $1 | sed 's%.*/%%'`" = "wine" ]; then
    builddir="./";
    injectdir="$1 ./"
    fdsdir="$1 ./"
    dirsep='\'
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v5"
cfs_code="${fdsdir}lalapps_ComputeFStatistic_v2"
if test $# -eq 0 ; then
    gct_code="${builddir}lalapps_HierarchSearchGCT"
else
    gct_code="$@"
fi

testDirBase="testGCT.d"
testDir="./${testDirBase}";
if [ -d "$testDir" ]; then
    rm -rf $testDir
fi
mkdir -p "$testDir"

SFTdir="${testDirBase}"
SFTfiles="$SFTdir${dirsep}*.sft"
SFTfiles_H1="$SFTdir${dirsep}H1-*.sft"
SFTfiles_L1="$SFTdir${dirsep}L1-*.sft"

## Tolerance of comparison
Tolerance=5e-2	## 5%

## ---------- fixed parameter of our test-signal -------------
Alpha=3.1
Delta=-0.5
h0=1.0
cosi=-0.3
psi=0.6
phi0=1.5
Freq=100.123456789
f1dot=-1e-9
f2dot=1e-18

## perfectly targeted search in sky
AlphaSearch=$Alpha
DeltaSearch=$Delta

## Produce skygrid file for the search
skygridfile="${testDir}/tmpskygridfile.dat"
echo "$AlphaSearch $DeltaSearch" > $skygridfile

mfd_FreqBand=0.20;
mfd_fmin=100;
numFreqBands=2;	## produce 'frequency-split' SFTs used in E@H

Dterms=8
RngMedWindow=50

gct_FreqBand=0.01
gct_F1dotBand=2.0e-10
gct_F2dotBand=0
gct_dFreq=0.000002 #"2.0e-6"
gct_dF1dot=1.0e-10
gct_dF2dot=0
gct_nCands=10

sqrtSh=1

## --------- Generate fake data set time stamps -------------
echo "----------------------------------------------------------------------"
echo " STEP 0: Generating time-stamps and segments files "
echo "----------------------------------------------------------------------"
echo

Tsft=1800
startTime=852443819
refTime=862999869
Tsegment=90000
if [ -n "$NSEGMENTS" ]; then
    Nsegments=${NSEGMENTS}
else
    Nsegments=3
fi


seggap=$(echo ${Tsegment} | awk '{printf "%.0f", $1 * 1.12345}')

tsFile_H1="${testDir}/timestampsH1.dat"  # for makefakedata
tsFile_L1="${testDir}/timestampsL1.dat"  # for makefakedata
segFile="${testDir}/segments.dat"

tmpTime=$startTime
iSeg=1
while [ $iSeg -le $Nsegments ]; do
    t0=$tmpTime
    t1=$(($t0 + $Tsegment))
    TspanHours=`echo $Tsegment | awk '{printf "%.7f", $1 / 3600.0 }'`
        ## first and last segment will be single-IFO only
    if [ $iSeg -eq 1 -o $iSeg -eq $Nsegments ]; then
        NSFT=`echo $Tsegment $Tsft |  awk '{print int(1.0 * $1 / $2 + 0.5) }'`
    else	## while all other segments are 2-IFO
        NSFT=`echo $Tsegment $Tsft |  awk '{print int(2.0 * $1 / $2 + 0.5) }'`
    fi
    echo "$t0 $t1 $TspanHours $NSFT" >> $segFile

    segs[$iSeg]=$tmpTime # save seg's beginning for later use
    echo "Segment: $iSeg of $Nsegments	GPS start time: ${segs[$iSeg]}"

    Tseg=$Tsft
    while [ $Tseg -le $Tsegment ]; do
        if [ $iSeg -ne 1 ]; then
	    echo "${tmpTime} 0" >> $tsFile_H1
        fi
        if [ $iSeg -ne $Nsegments ]; then	            ## we skip segment N for L1
	    echo "${tmpTime} 0" >> $tsFile_L1
        fi
	tmpTime=$(($tmpTime + $Tsft))
	Tseg=$(($Tseg + $Tsft))
    done

    tmpTime=$(($tmpTime + $seggap))
    iSeg=$(($iSeg + 1))
done

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Signal"
echo "----------------------------------------------------------------------"
echo

FreqStep=`echo $mfd_FreqBand $numFreqBands |  awk '{print $1 / $2}'`
mfd_fBand=`echo $FreqStep $Tsft |  awk '{print ($1 - 1.5 / $2)}'`	## reduce by 1/2 a bin to avoid including last freq-bins

# construct common MFD cmd
mfd_CL_common="--Band=${mfd_fBand} --injectionSources=\"{Freq=$Freq; f1dot=$f1dot; f2dot=$f2dot; Alpha=$Alpha; Delta=$Delta; psi=$psi; phi0=$phi0; h0=$h0; cosi=$cosi; refTime=$refTime}\" --Tsft=$Tsft --randSeed=1000 --outSingleSFT --IFOs=H1,L1 --timestampsFiles=${tsFile_H1},${tsFile_L1}"

if [ "$sqrtSh" != "0" ]; then
    mfd_CL_common="$mfd_CL_common --sqrtSX=${sqrtSh},${sqrtSh}";
fi

iFreq=1
while [ $iFreq -le $numFreqBands ]; do
    mfd_fi=`echo $mfd_fmin $iFreq $FreqStep | awk '{print $1 + ($2 - 1) * $3}'`

    cmdline="$mfd_code $mfd_CL_common --fmin=$mfd_fi --outSFTdir=${SFTdir} --outLabel=freqBand$iFreq"
    if [ -n "$DEBUG" ]; then
        cmdline="$cmdline"
    else
        cmdline="$cmdline &> /dev/null"
    fi
    echo "$cmdline";
    if ! eval "$cmdline"; then
        echo "Error.. something failed when running '$mfd_code' ..."
        exit 1
    fi

    iFreq=$(( $iFreq + 1 ))

done


echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run CFSv2 over segments"
echo "----------------------------------------------------------------------"
echo
outfile_cfs="${testDir}/CFS.dat";

if [ ! -r "$outfile_cfs" ]; then
    tmpfile_cfs="${testDir}/__tmp_CFS.dat";
    cfs_CL_common=" --Alpha=$Alpha --Delta=$Delta --Freq=$Freq --f1dot=$f1dot --f2dot=$f2dot --outputLoudest=$tmpfile_cfs --refTime=$refTime --Dterms=$Dterms --RngMedWindow=$RngMedWindow --outputSingleFstats"
    if [ "$sqrtSh" = "0" ]; then
        cfs_CL_common="$cfs_CL_common --SignalOnly";
    fi

    TwoFsum=0
    TwoFsum_L1=0
    TwoFsum_H1=0

    for ((iSeg=1; iSeg <= $Nsegments; iSeg++)); do
        minStartGPS=${segs[$iSeg]}
        maxStartGPS=$(($minStartGPS + $Tsegment))
        cfs_CL="$cfs_code $cfs_CL_common --minStartTime=$minStartGPS --maxStartTime=$maxStartGPS"

        # ----- get multi-IFO + single-IFO F-stat values
        cmdline="$cfs_CL --DataFiles='$SFTfiles'"
        if [ -n "$DEBUG" ]; then
            cmdline="$cmdline"
        else
            cmdline="$cmdline &> /dev/null"
        fi
        echo "$cmdline"
        if ! eval "$cmdline"; then
	    echo "Error.. something failed when running '$cfs_code' ..."
	    exit 1
        fi

        resCFS=$(cat ${tmpfile_cfs} | awk '{if($1=="twoF") {printf "%.11g", $3}}')
        TwoFsum=$(echo $TwoFsum $resCFS | awk '{printf "%.11g", $1 + $2}')

        if [ $iSeg -eq 1 ]; then	## segment 1 has no H1 SFTs
            resCFS_L1=$(cat ${tmpfile_cfs} | awk '{if($1=="twoF0") {printf "%.11g", $3}}')
            TwoFsum_L1=$(echo $TwoFsum_L1 $resCFS_L1 | awk '{printf "%.11g", $1 + $2}')
        elif [ $iSeg -eq $Nsegments ]; then	## segment N has no L1 SFTs
            resCFS_H1=$(cat ${tmpfile_cfs} | awk '{if($1=="twoF0") {printf "%.11g", $3}}')	## therefore 'H1' is the first and only detector
            TwoFsum_H1=$(echo $TwoFsum_H1 $resCFS_H1 | awk '{printf "%.11g", $1 + $2}')
        else
            resCFS_H1=$(cat ${tmpfile_cfs} | awk '{if($1=="twoF0") {printf "%.11g", $3}}')	## 'H1' is first
            TwoFsum_H1=$(echo $TwoFsum_H1 $resCFS_H1 | awk '{printf "%.11g", $1 + $2}')

            resCFS_L1=$(cat ${tmpfile_cfs} | awk '{if($1=="twoF1") {printf "%.11g", $3}}')	## 'L1' second
            TwoFsum_L1=$(echo $TwoFsum_L1 $resCFS_L1 | awk '{printf "%.11g", $1 + $2}')
        fi
    done

    TwoFAvg=$(echo    $TwoFsum    $Nsegments | awk '{printf "%-12.11g", $1 / ($2)}')
    TwoFAvg_H1=$(echo $TwoFsum_H1 $Nsegments | awk '{printf "%-12.11g", $1 / ($2-1)}')	## H1 has one segment less (the first one)
    TwoFAvg_L1=$(echo $TwoFsum_L1 $Nsegments | awk '{printf "%-12.11g", $1 / ($2-1)}')	## L1 also one segment less (the last one)
    echo "$TwoFAvg	$TwoFAvg_H1	$TwoFAvg_L1" > $outfile_cfs
else
    echo "CFS result file '$outfile_cfs' exists already ... reusing it"
    cfs_res=$(cat $outfile_cfs)
    TwoFAvg=$(echo $cfs_res | awk '{print $1}')
    TwoFAvg_H1=$(echo $cfs_res | awk '{print $2}')
    TwoFAvg_L1=$(echo $cfs_res | awk '{print $3}')
fi

echo
echo "==>   Average <2F_multi>=$TwoFAvg, <2F_H1>=$TwoFAvg_H1, <2F_L1>=$TwoFAvg_L1"

## ---------- run GCT code on this data ----------------------------------------

gct_CL_common="--gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles' --skyGridFile='$skygridfile' --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand  --df2dot=$gct_dF2dot --f2dot=$f2dot --f2dotBand=$gct_F2dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segFile --Dterms=$Dterms --blocksRngMed=$RngMedWindow"
if [ "$sqrtSh" = "0" ]; then
    gct_CL_common="$gct_CL_common --SignalOnly";
fi

BSGL_flags="--computeBSGL --Fstar0=10 --oLGX='0.5,0.5' --recalcToplistStats"

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 3: run HierarchSearchGCT using Resampling (perfect match) and segment-list file and --recalcToplistStats"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT_RS="${testDir}/GCT_RS.dat"
timingsfile_RS="${testDir}/timing_RS.dat"

cmdline="$gct_code $gct_CL_common --FstatMethod=ResampGeneric --fnameout='$outfile_GCT_RS' --outputTiming='$timingsfile_RS' --recalcToplistStats ${BSGL_flags}"
if [ -n "$DEBUG" ]; then
    cmdline="$cmdline"
else
    cmdline="$cmdline &> /dev/null"
fi
echo "$cmdline"
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi
topline=$(sort -nr -k7,7 $outfile_GCT_RS | head -1)
freqGCT_RS=$(echo $topline | awk '{print $1}')
resGCT_RS=$(echo $topline | awk '{print $7}')
resGCT_RS_H1=$(echo $topline | awk '{print $9}')
resGCT_RS_L1=$(echo $topline | awk '{print $10}')
resGCT_RSr=$(echo $topline  | awk '{print $11}')
resGCT_RSr_H1=$(echo $topline  | awk '{print $13}')
resGCT_RSr_L1=$(echo $topline  | awk '{print $14}')

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 4: run HierarchSearchGCT using LALDemod (perfect match) and --tStack and --nStacksMax and --recalcToplistStats"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT_DM="${testDir}/GCT_DM.dat"
timingsfile_DM="${testDir}/timing_DM.dat"

cmdline="$gct_code $gct_CL_common --FstatMethod=DemodOptC --fnameout='$outfile_GCT_DM' --outputTiming='$timingsfile_DM' ${BSGL_flags}"
if [ -n "$DEBUG" ]; then
    cmdline="$cmdline"
else
    cmdline="$cmdline &> /dev/null"
fi

echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

topline=$(sort -nr -k7,7 $outfile_GCT_DM | head -1)
resGCT_DM=$(echo $topline  | awk '{print $7}')
resGCT_DM_H1=$(echo $topline  | awk '{print $9}')
resGCT_DM_L1=$(echo $topline  | awk '{print $10}')
freqGCT_DM=$(echo $topline | awk '{print $1}')
resGCT_DMr=$(echo $topline  | awk '{print $11}')
resGCT_DMr_H1=$(echo $topline  | awk '{print $13}')
resGCT_DMr_L1=$(echo $topline  | awk '{print $14}')

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 5: run HierarchSearchGCT using LALDemod (perfect match) and --tStack and --nStacksMax and --computeBSGL"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT_DM_BSGL="${testDir}/GCT_DM_BSGL.dat"
timingsfile_DM_BSGL="${testDir}/timing_DM_BSGL.dat"

cmdline="$gct_code $gct_CL_common --FstatMethod=DemodOptC ${BSGL_flags} --SortToplist=2 --fnameout='$outfile_GCT_DM_BSGL' --outputTiming='$timingsfile_DM_BSGL'"
if [ -n "$DEBUG" ]; then
    cmdline="$cmdline"
else
    cmdline="$cmdline &> /dev/null"
fi

echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

topline=$(sort -nr -k7,7 $outfile_GCT_DM_BSGL | head -1)
resGCT_DM_BSGL=$(echo $topline  | awk '{print $7}')
resGCT_DM_H1_BSGL=$(echo $topline  | awk '{print $9}')
resGCT_DM_L1_BSGL=$(echo $topline  | awk '{print $10}')
freqGCT_DM_BSGL=$(echo $topline | awk '{print $1}')

echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 6: run HierarchSearchGCT using LALDemod (perfect match) with 'dual' toplist: 1st=F, 2nd=BSGL"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_GCT_DM_DUAL="${testDir}/GCT_DM_DUAL.dat"
timingsfile_DM_DUAL="${testDir}/timing_DM_DUAL.dat"

cmdline="$gct_code $gct_CL_common --FstatMethod=DemodOptC --SortToplist=3 ${BSGL_flags} --fnameout='$outfile_GCT_DM_DUAL' --outputTiming='$timingsfile_DM_DUAL'"
if [ -n "$DEBUG" ]; then
    cmdline="$cmdline"
else
    cmdline="$cmdline &> /dev/null"
fi

echo $cmdline
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

diff_F=`diff -I "[%][%].*" ${outfile_GCT_DM} ${outfile_GCT_DM_DUAL}`
diff_BSGL=`diff -I "[%][%].*" ${outfile_GCT_DM_BSGL} ${outfile_GCT_DM_DUAL}-BSGL`
if [ -n "${diff_F}" -o -n "${diff_BSGL}" ]; then
    echo "Error: 'dual' toplist handling seems to differ from indidual 'F' or 'BSGL'-sorted toplists"
    exit 1
fi

## ---------- compute relative differences and check against tolerance --------------------
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'

freqreldev_RS=$(echo $Freq $freqGCT_RS | awk "$awk_reldev")
reldev_RS=$(echo $TwoFAvg $resGCT_RS | awk "$awk_reldev")
reldev_RS_H1=$(echo $TwoFAvg_H1 $resGCT_RS_H1 | awk "$awk_reldev")
reldev_RS_L1=$(echo $TwoFAvg_L1 $resGCT_RS_L1 | awk "$awk_reldev")

reldev_RSr=$(echo $TwoFAvg $resGCT_RSr | awk "$awk_reldev")
reldev_RSr_H1=$(echo $TwoFAvg_H1 $resGCT_RSr_H1 | awk "$awk_reldev")
reldev_RSr_L1=$(echo $TwoFAvg_L1 $resGCT_RSr_L1 | awk "$awk_reldev")

reldev_DM=$(echo $TwoFAvg $resGCT_DM | awk "$awk_reldev")
reldev_DM_H1=$(echo $TwoFAvg_H1 $resGCT_DM_H1 | awk "$awk_reldev")
reldev_DM_L1=$(echo $TwoFAvg_L1 $resGCT_DM_L1 | awk "$awk_reldev")

reldev_DMr=$(echo $TwoFAvg $resGCT_DMr | awk "$awk_reldev")
reldev_DMr_H1=$(echo $TwoFAvg_H1 $resGCT_DMr_H1 | awk "$awk_reldev")
reldev_DMr_L1=$(echo $TwoFAvg_L1 $resGCT_DMr_L1 | awk "$awk_reldev")

freqreldev_DM=$(echo $Freq $freqGCT_DM | awk "$awk_reldev")
reldev_DM_BSGL=$(echo $TwoFAvg $resGCT_DM_BSGL | awk "$awk_reldev")
reldev_DM_H1_BSGL=$(echo $TwoFAvg_H1 $resGCT_DM_H1_BSGL | awk "$awk_reldev")
reldev_DM_L1_BSGL=$(echo $TwoFAvg_L1 $resGCT_DM_L1_BSGL | awk "$awk_reldev")
freqreldev_DM_BSGL=$(echo $Freq $freqGCT_DM_BSGL | awk "$awk_reldev")

# ---------- Check relative deviations against tolerance, report results ----------
retstatus=0
awk_isgtr='{if($1>$2) {print "1"}}'

echo
echo "--------- Timings ------------------------------------------------------------------------------------------------"
awk_timing='{printf "c0ic = %-6.1e s, c1co = %-6.1e s, c0Demod = %-6.1e s,  (%s)", $8, $9, $10, $11, $12}'
timing_DM=$(sed '/^%.*/d' $timingsfile_DM | awk "$awk_timing")
timing_DM_BSGL=$(sed '/^%.*/d' $timingsfile_DM_BSGL | awk "$awk_timing")
timing_RS=$(sed '/^%.*/d' $timingsfile_RS | awk "$awk_timing")
echo " GCT-LALDemod:      $timing_DM"
echo " GCT-LALDemod-BSGL:  $timing_DM_BSGL"
echo " GCT-Resamp:        $timing_RS"

echo
echo "--------- Compare results ----------------------------------------------------------------------------------------"
echo "                     	<2F_multi>	<2F_H1>  	<2F_L1>  	 @ Freq [Hz]     	(reldev, reldev_H1, reldev_L1, reldev_Freq)"
echo    "==>  CFSv2:         	$TwoFAvg 	$TwoFAvg_H1   	$TwoFAvg_L1   	 @ $Freq 	[Tolerance = ${Tolerance}]"

echo -n "==>  GCT-DM: 		$resGCT_DM	$resGCT_DM_H1	$resGCT_DM_L1  	 @ $freqGCT_DM 	($reldev_DM, $reldev_DM_H1, $reldev_DM_L1, $freqreldev_DM)"
fail1=$(echo $freqreldev_DM $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_DM $Tolerance     | awk "$awk_isgtr")
fail3=$(echo $reldev_DM_H1 $Tolerance  | awk "$awk_isgtr")
fail4=$(echo $reldev_DM_L1 $Tolerance  | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo " ==> *FAILED*"
    retstatus=1
else
    echo " ==> OK"
fi

echo -n "==>  GCT-DM-recalc:	$resGCT_DMr	$resGCT_DMr_H1	$resGCT_DMr_L1	 @ $freqGCT_DM	($reldev_DMr, $reldev_DMr_H1, $reldev_DMr_L1, $freqreldev_DM)"
fail2r=$(echo $reldev_DMr $Tolerance     | awk "$awk_isgtr")
fail3r=$(echo $reldev_DMr_H1 $Tolerance  | awk "$awk_isgtr")
fail4r=$(echo $reldev_DMr_L1 $Tolerance  | awk "$awk_isgtr")
if [ "$fail2r" -o "$fail3r" -o "$fail4r" ]; then
    echo " ==> *FAILED*"
    retstatus=1
else
    echo " ==> OK"
fi

echo -n "==>  GCT-LALDemodBSGL:	$resGCT_DM_BSGL 	$resGCT_DM_H1_BSGL	$resGCT_DM_L1_BSGL	 @ $freqGCT_DM_BSGL 	($reldev_DM_BSGL, $reldev_DM_H1_BSGL, $reldev_DM_L1_BSGL, $freqreldev_DM_BSGL)"
fail1=$(echo $freqreldev_DM_BSGL $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_DM_BSGL $Tolerance     | awk "$awk_isgtr")
fail3=$(echo $reldev_DM_H1_BSGL $Tolerance  | awk "$awk_isgtr")
fail4=$(echo $reldev_DM_L1_BSGL $Tolerance  | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo " ==> *FAILED*"
    retstatus=1
else
    echo " ==> OK"
fi

echo -n "==>  GCT-Resamp: 	$resGCT_RS 	$resGCT_RS_H1 	$resGCT_RS_L1  	 @ $freqGCT_RS 	($reldev_RS, $reldev_RS_H1, $reldev_RS_L1, $freqreldev_RS)"
fail1=$(echo $freqreldev_RS $Tolerance | awk "$awk_isgtr")
fail2=$(echo $reldev_RS     $Tolerance | awk "$awk_isgtr")
fail3=$(echo $reldev_RS_H1  $Tolerance | awk "$awk_isgtr")
fail4=$(echo $reldev_RS_L1  $Tolerance | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" -o "$fail3" -o "$fail4" ]; then
    echo " ==> *FAILED*"
    retstatus=1
else
    echo " ==> OK"
fi

echo -n "==>  GCT-RS-recalc: 	$resGCT_RSr 	$resGCT_RSr_H1 	$resGCT_RSr_L1  	 @ $freqGCT_RS 	($reldev_RSr, $reldev_RSr_H1, $reldev_RSr_L1, $freqreldev_RS)"
fail2r=$(echo $reldev_RSr     $Tolerance | awk "$awk_isgtr")
fail3r=$(echo $reldev_RSr_H1  $Tolerance | awk "$awk_isgtr")
fail4r=$(echo $reldev_RSr_L1  $Tolerance | awk "$awk_isgtr")
if [ "$fail2r" -o "$fail3r" -o "$fail4r" ]; then
    echo " ==> *FAILED*"
    retstatus=1
else
    echo " ==> OK"
fi


echo "----------------------------------------------------------------------"

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDir
    echo "Cleaned up."
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old

exit $retstatus
