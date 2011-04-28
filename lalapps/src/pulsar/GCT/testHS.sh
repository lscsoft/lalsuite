#!/bin/bash

NORESAMP="1"
#NOCLEANUP="1"
#SEPIFOVETO="1"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
    fdsdir="../FDS_isolated/"
else
    srcdir=.
fi
dirsep=/
if [ "`echo $1 | sed 's%.*/%%'`" = "wine" ]; then
    builddir="./";
    injectdir="$1 ./"
    fdsdir="$1 ./"
    dirsep='\'
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
pfs_code="${fdsdir}lalapps_PredictFStat"
if test $# -eq 0 ; then
    gct_code="${builddir}lalapps_HierarchSearchGCT"
else
    gct_code="$@"
fi

SFTdir="TestSFTs"
SFTfiles="$SFTdir${dirsep}*.sft"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    if [ -n "$LALPULSAR_PREFIX" ]; then
	export LAL_DATA_PATH="${LALPULSAR_PREFIX}/share/lalpulsar";
    else
	echo
	echo "Need environment-variable LALPULSAR_PREFIX, or LAL_DATA_PATH to be set"
	echo "to your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
	echo "This might indicate an incomplete LAL+LALPULSAR installation"
	echo
	exit 1
    fi
fi


## Tolerance of comparison
Tolerance="0.05"

## ---------- fixed parameter of our test-signal -------------
Alpha="3.1"
Delta="-0.5"
h0="1.0"
cosi="-0.3"
psi="0.6"
phi0="1.5"
Freq="100.123456789"
f1dot="-1e-9"

AlphaSearch=$Alpha
DeltaSearch=$Delta

## Produce skygrid file for the search
skygridfile="tmpskygridfile.dat"
echo $AlphaSearch" "$DeltaSearch > $skygridfile

mfd_FreqBand=0.20;
mfd_fmin=100;
numFreqBands=4;	## produce 'frequency-split' SFTs used in E@H

gct_FreqBand="0.01"
gct_F1dotBand="2.0e-10"
gct_dFreq="0.000002" #"2.0e-6"
gct_dF1dot="1.0e-10"
gct_nCands="1000"

noiseSqrtSh="0"

## --------- Generate fake data set time stamps -------------
echo "----------------------------------------------------------------------"
echo " STEP 0: Generating fake data time-stamps file "
echo "----------------------------------------------------------------------"
echo

Tsft="1800"
startTime="852443819"
refTime="862999869"
Tsegment="90000"
Nsegments="14"
seggap=$(echo "scale=0; ${Tsegment} * 1.12345" | bc)
tsfile="timestampsTEST.txt"
segFile="segments.txt"
rm -rf $tsfile $segFile
tmpTime=$startTime
ic1="1"
while [ "$ic1" -le "$Nsegments" ];
do
    t0=$tmpTime
    t1=`echo $t0 $Tsegment | LC_ALL=C awk '{print $1 + $2}'`
    TspanHours=`echo $Tsegment | LC_ALL=C awk '{printf "%.7f", $1 / 3600.0 }'`
    NSFT=`echo $Tsegment $Tsft | LC_ALL=C awk '{print int(2.0 * $1 / $2 + 0.5) }'`
    echo "$t0 $t1 $TspanHours $NSFT" >> $segFile
    segs[${ic1}]=$tmpTime # save seg's beginning for later use
    echo "Segment: "$ic1" of "$Nsegments"   GPS start time: "${segs[${ic1}]}

    ic2=$Tsft
    while [ "$ic2" -le "$Tsegment" ];
    do
	echo ${tmpTime}" 0" >> $tsfile
	tmpTime=$(echo "scale=0; ${tmpTime} + ${Tsft}" | bc)
	ic2=$(echo "scale=0; ${ic2} + ${Tsft}" | bc)
    done

    tmpTime=$(echo "scale=0; ${tmpTime} + ${seggap}" | bc | LC_ALL=C awk '{printf "%.0f",$1}')
    ic1=$(echo "scale=0; ${ic1} + 1" | bc)
done

## ------------------------------------------------------------

if [ "$noiseSqrtSh" != 0 ]; then
    sqrtSh=$noiseSqrtSh
    haveNoise=true;
else
    sqrtSh=1;	## for SemiAnalyticF signal-only case
    haveNoise=false;
fi


##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Signal"
echo "----------------------------------------------------------------------"
echo
if [ ! -d "$SFTdir" ]; then
    mkdir -p $SFTdir;
else
    rm -f $SFTdir/*;
fi

FreqStep=`echo $mfd_FreqBand $numFreqBands | LC_ALL=C awk '{print $1 / $2}'`
mfd_fBand=`echo $FreqStep $Tsft | LC_ALL=C awk '{print ($1 - 1.5 / $2)}'`	## reduce by 1/2 a bin to avoid including last freq-bins
iFreq=1
while [ "$iFreq" -le "$numFreqBands" ]; do
    mfd_fi=`echo $mfd_fmin $iFreq $FreqStep | LC_ALL=C awk '{print $1 + ($2 - 1) * $3}'`

    # construct common MFD cmd
    mfd_CL_common=" --fmin=$mfd_fi --Band=${mfd_fBand} --Freq=$Freq --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --h0=$h0 --cosi=$cosi --ephemYear=05-09 --generationMode=1 --timestampsFile=$tsfile --refTime=$refTime --Tsft=$Tsft --randSeed=1000 --outSingleSFT"

    if [ "$haveNoise" = true ]; then
        mfd_CL_common="$mfd_CL_common --noiseSqrtSh=$sqrtSh";
    fi

    # for H1:
    cmdline="$mfd_code $mfd_CL_common --IFO=H1 --outSFTbname=${SFTdir}\\${dirsep}H1-${mfd_fi}_${FreqStep}.sft";
    echo "$cmdline";
    if ! eval "$cmdline"; then
        echo "Error.. something failed when running '$mfd_code' ..."
        exit 1
    fi

    # for L1:
    cmdline="$mfd_code $mfd_CL_common --IFO=L1 --outSFTbname=${SFTdir}\\${dirsep}L1-${mfd_fi}_${FreqStep}.sft";
    echo "$cmdline";
    if ! eval "$cmdline"; then
        echo "Error.. something failed when running '$mfd_code' ..."
        exit 1
    fi

    iFreq=$(( $iFreq + 1 ))

done


echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run PredictFstat (perfect match) "
echo "----------------------------------------------------------------------"
echo


TwoFsum="0"
TwoFsum1="0"
TwoFsum2="0"

for ((x=1; x <= $Nsegments; x++))
  do
    outfile_pfs="__tmp_PFS.dat";

    startGPS=${segs[${x}]}
    endGPS=$(echo "scale=0; ${startGPS} + ${Tsegment}" | bc | awk '{printf "%.0f",$1}')
    #echo "Segment: "$x"  "$startGPS" "$endGPS
    # construct pfs cmd
    pfs_CL=" --Alpha=$Alpha --Delta=$Delta --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='$SFTfiles' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS"
    if [ "$haveNoise" = false ]; then
        pfs_CL="$pfs_CL --SignalOnly";
    fi

    cmdline="$pfs_code $pfs_CL"
    echo "  "$cmdline
    if ! tmp=`eval $cmdline`; then
	echo "Error.. something failed when running '$pfs_code' ..."
	exit 1
    fi
    resPFS=$(cat ${outfile_pfs} | grep 'twoF_expected' | awk -F';' '{print $1}' | awk '{print $3}')
    #resPFS=`echo $tmp | awk '{printf "%g", $1}'`
    TwoFsum=$(echo "scale=6; ${TwoFsum} + ${resPFS}" | bc);

    if [ -n "$SEPIFOVETO" ]; then
	# H1 only
	pfs_CL=" --Alpha=$Alpha --Delta=$Delta --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='${SFTfiles}H1*' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS"
        if [ "$haveNoise" = false ]; then
            pfs_CL="$pfs_CL --SignalOnly";
        fi
	cmdline="$pfs_code $pfs_CL"
	echo "  "$cmdline
	if ! tmp=`eval $cmdline`; then
            echo "Error.. something failed when running '$pfs_code' ..."
            exit 1
	fi
	resPFS1=$(cat ${outfile_pfs} | grep 'twoF_expected' | awk -F';' '{print $1}' | awk '{print $3}')
	TwoFsum1=$(echo "scale=6; ${TwoFsum1} + ${resPFS1}" | bc);

	# L1 only
	pfs_CL=" --Alpha=$Alpha --Delta=$Delta --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='${SFTfiles}L1*' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS"
        if [ "$haveNoise" = false ]; then
            pfs_CL="$pfs_CL --SignalOnly";
        fi
        cmdline="$pfs_code $pfs_CL"
        echo "  "$cmdline
        if ! tmp=`eval $cmdline`; then
            echo "Error.. something failed when running '$pfs_code' ..."
	    exit 1
        fi
        resPFS2=$(cat ${outfile_pfs} | grep 'twoF_expected' | awk -F';' '{print $1}' | awk '{print $3}')
        TwoFsum2=$(echo "scale=6; ${TwoFsum2} + ${resPFS2}" | bc);

	echo "Segment: "$x"   2F: "$resPFS"    (H1 only: "$resPFS1"  L1 only: "$resPFS2")"
    else
	echo "Segment: "$x"   2F: "$resPFS
    fi
    echo
done
TwoFsum=$(echo "scale=6; ${TwoFsum} / ${Nsegments}" | bc);
if [ -n "$SEPIFOVETO" ]; then
    TwoFsum1=$(echo "scale=6; ${TwoFsum1} / ${Nsegments}" | bc);
    TwoFsum2=$(echo "scale=6; ${TwoFsum2} / ${Nsegments}" | bc);
    echo
    echo "==>   Average 2F: "$TwoFsum"  in H1: "$TwoFsum1"  in L1: "$TwoFsum2
else
    echo
    echo "==>   Average 2F: "$TwoFsum
fi

edat="earth05-09.dat"
sdat="sun05-09.dat"

echo
echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run HierarchSearchGCT using Resampling (perfect match)"
echo "----------------------------------------------------------------------"
echo

if [ -e "checkpoint.cpt" ]; then
    rm checkpoint.cpt # delete checkpoint to start correctly
fi

outfile_gct1="__tmp_GCT1.dat"

gct_CL=" --useResamp --fnameout=$outfile_gct1 --gridType1=3 --tStack=$Tsegment --nCand1=$gct_nCands --nStacksMax=$Nsegments --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles'  --ephemE=$edat --ephemS=$sdat --skyGridFile='./$skygridfile' --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime "
if [ "$haveNoise" = false ]; then
    gct_CL="$gct_CL --SignalOnly";
fi

cmdline="$gct_code $gct_CL"
echo $cmdline
if [ -z "$NORESAMP" ]; then
    if ! tmp=`eval $cmdline`; then
	echo "Error.. something failed when running '$gct_code' ..."
	exit 1
    fi
    resGCT1=$(cat $outfile_gct1 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $6}')
    freqGCT1=$(cat $outfile_gct1 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $1}')
else
    echo
    echo "Not run with resampling."
fi


echo
echo "----------------------------------------------------------------------"
echo " STEP 4: run HierarchSearchGCT without Resampling (perfect match)"
echo "----------------------------------------------------------------------"
echo

if [ -e "checkpoint.cpt" ]; then
    rm checkpoint.cpt # delete checkpoint to start correctly
fi

outfile_gct2="__tmp_GCT2.dat"

gct_CL=" --fnameout=$outfile_gct2 --gridType1=3 --tStack=$Tsegment --nCand1=$gct_nCands --nStacksMax=$Nsegments --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles'  --ephemE=$edat --ephemS=$sdat --skyGridFile='./$skygridfile'  --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime "
if [ "$haveNoise" = false ]; then
    gct_CL="$gct_CL --SignalOnly";
fi

if [ -n "$SEPIFOVETO" ]; then
    gct_CL=$gct_CL" --SepDetVeto"
fi

cmdline="$gct_code $gct_CL > >(tee stdout.log) 2> >(tee stderr.log >&2)"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi
resGCT2=$(cat $outfile_gct2 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $6}')
freqGCT2=$(cat $outfile_gct2 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $1}')

if [ -n "$SEPIFOVETO" ]; then
    resGCT2H1=$(cat stderr.log | grep "average2F\[1\]=" | awk '{print $4}')
    resGCT2L1=$(cat stderr.log | grep "average2F\[2\]=" | awk '{print $4}')
    oZ1=$(cat stderr.log | grep "average2F\[1\]=" | awk '{print $6}')
    oZ2=$(cat stderr.log | grep "average2F\[2\]=" | awk '{print $6}')
fi


if [ -z "$NORESAMP" ]; then
reldev1=$(echo "scale=5; ($TwoFsum - $resGCT1)/(0.5 * ($TwoFsum + $resGCT1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
freqreldev1=$(echo "scale=13; (($Freq - $freqGCT1)/$Freq) " | bc | awk '{ if($1>=0) {printf "%.6f",$1} else {printf "%.6f",$1*(-1)}}')
reldev3=$(echo "scale=5; ($resGCT1 - $resGCT2)/(0.5 * ($resGCT2 + $resGCT1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
fi

reldev2=$(echo "scale=5; ($TwoFsum - $resGCT2)/(0.5 * ($TwoFsum + $resGCT2))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
freqreldev2=$(echo "scale=13; (($Freq - $freqGCT2)/$Freq) " | bc | awk '{ if($1>=0) {printf "%.13f",$1} else {printf "%.12f",$1*(-1)}}')

if [ -n "$SEPIFOVETO" ]; then
    reldev2H1=$(echo "scale=5; ($TwoFsum1 - $resGCT2H1)/(0.5 * ($TwoFsum1 + $resGCT2H1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
    reldev2L1=$(echo "scale=5; ($TwoFsum2 - $resGCT2L1)/(0.5 * ($TwoFsum2 + $resGCT2L1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
fi
freqreldev2B=$(echo "scale=13; (($Freq - $freqGCT2)/${gct_dFreq})" | bc | awk '{ if($1>=0) {printf "%.12f",$1} else {printf "%.12f",$1*(-1)}}')


echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 5: re-run HierarchSearchGCT with a segment list file instead of --tStack and --nStacksMax"
echo "----------------------------------------------------------------------------------------------------"
echo
if [ -e "checkpoint.cpt" ]; then
    rm checkpoint.cpt # delete checkpoint to start correctly
fi

outfile_gct5="__tmp_GCT5.dat"

gct_CL=" --SignalOnly --fnameout=$outfile_gct5 --gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles'  --ephemE=$edat --ephemS=$sdat --skyGridFile='./$skygridfile'  --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segFile -d1"

cmdline="$gct_code $gct_CL > >(tee stdout.log) 2> >(tee stderr.log >&2)"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi
resGCT5=$(cat $outfile_gct5 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $6}')
freqGCT5=$(cat $outfile_gct5 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $1}')
reldev5=$(echo "scale=5; ($resGCT2 - $resGCT5)/(0.5 * ($resGCT2 + $resGCT5))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')

echo
echo "----------------------------------------------------------------------"
echo "==>  Predicted:      "$TwoFsum


# Check predicted 2F against search code output

if [ -z "$NORESAMP" ]; then
if [ `echo $reldev1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, Resamp:    "$resGCT1"  ("$reldev1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, Resamp:    "$resGCT1"  ("$reldev1")     OK."
fi
fi

if [ `echo $reldev2" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, no Resamp: "$resGCT2"  ("$reldev2")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, no Resamp: "$resGCT2"  ("$reldev2")     OK."
fi

if [ `echo $reldev5" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, seg-list: "$resGCT5"  ("$reldev5")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, seg-list: "$resGCT5"  ("$reldev5")     OK."
fi


if [ -z "$NORESAMP" ]; then
if [ `echo $reldev3" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, Resamp vs. no-Resamp:     "$reldev3
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, Resamp vs. no-Resamp:    "$reldev3"      OK."
fi
fi

# Check relative error in frequency
if [ -z "$NORESAMP" ]; then
echo
echo "==>  Signal frequency: "$Freq"  Found at: "$freqGCT1
if [ `echo $freqreldev1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, Resamp.     Rel. dev. in frequency: "$freqreldev1
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, Resamp.     Rel. dev. in frequency: "$freqreldev1"   OK."
fi
fi

echo
echo "==>  Signal frequency: "$Freq"  Found at: "$freqGCT2
if [ `echo $freqreldev2" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, no Resamp.  Rel. dev. in frequency: "$freqreldev2
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, no Resamp.  Rel. dev. in frequency: "$freqreldev2"   OK."
    echo "     offset as fraction of frequency bins: "$freqreldev2B
fi


if [ -n "$SEPIFOVETO" ]; then
    echo

    if [ `echo $reldev2H1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
	echo "==>  GCT, H1 only deviation:     "$reldev2H1
	echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
	exit 2
    else
	echo "==>  GCT, H1 only deviation:    "$reldev2H1"      OK."
    fi

    if [ `echo $reldev2L1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
        echo "==>  GCT, L1 only deviation:     "$reldev2L1
        echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
        exit 2
    else
        echo "==>  GCT, L1 only deviation:    "$reldev2L1"      OK."
    fi


    Z1=$(echo "scale=6; ($TwoFsum / $TwoFsum1)" | bc )
    Z2=$(echo "scale=6; ($TwoFsum / $TwoFsum2)" | bc )
    echo "( Predicted: Z = <F>_H1L1 / <F>_H1 = "$Z1"   Obtained: "$oZ1" )"
    echo "( Predicted: Z = <F>_H1L1 / <F>_L1 = "$Z2"   Obtained: "$oZ2" )"
fi


echo "----------------------------------------------------------------------"


## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $skygridfile $tsfile $outfile_pfs $outfile_gct1 $outfile_gct2 $outfile_gct5 checkpoint.cpt stderr.log stdout.log $segFile
    echo "Cleaned up."
fi
