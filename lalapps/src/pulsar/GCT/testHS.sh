#!/bin/bash

#NORESAMP="1"
#NOCLEANUP="1"

## make sure we work in 'C' locale here to avoid awk sillyness
LC_ALL_old=$LC_ALL
export LC_ALL=C

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

AlphaSearch=$Alpha
DeltaSearch=$Delta

## Produce skygrid file for the search
skygridfile="tmpskygridfile.dat"
echo $AlphaSearch" "$DeltaSearch > $skygridfile

mfd_FreqBand=0.20;
mfd_fmin=100;
numFreqBands=4;	## produce 'frequency-split' SFTs used in E@H

gct_FreqBand=0.01
gct_F1dotBand=2.0e-10
gct_dFreq=0.000002 #"2.0e-6"
gct_dF1dot=1.0e-10
gct_nCands=1000

sqrtSh=0

## --------- Generate fake data set time stamps -------------
echo "----------------------------------------------------------------------"
echo " STEP 0: Generating fake data time-stamps file "
echo "----------------------------------------------------------------------"
echo

Tsft=1800
startTime=852443819
refTime=862999869
Tsegment=90000
Nsegments=14
seggap=$(echo ${Tsegment} | awk '{printf "%.0f", $1 * 1.12345}')

tsFile="timestampsTEST.txt"
if [ -r $tsFile ]; then
    have_tsFile=true;
    echo "Timestamps file '$tsFile' already exists ... reusing it."
fi
segFile="segments.txt"
if [ -r $segFile ]; then
    have_segFile=true;
    echo "Segments file '$segFile' already exists ... reusing it."
fi
echo

tmpTime=$startTime
iSeg=1
while [ $iSeg -le $Nsegments ]; do
    t0=$tmpTime
    ## only write segment-file if it is not found already
    if [ "$have_segFile" != "true" ]; then
        t1=$(($t0 + $Tsegment))
        TspanHours=`echo $Tsegment | awk '{printf "%.7f", $1 / 3600.0 }'`
        NSFT=`echo $Tsegment $Tsft |  awk '{print int(2.0 * $1 / $2 + 0.5) }'`
        echo "$t0 $t1 $TspanHours $NSFT" >> $segFile
    fi
    segs[$iSeg]=$tmpTime # save seg's beginning for later use
    echo "Segment: $iSeg of $Nsegments	GPS start time: ${segs[$iSeg]}"

    Tseg=$Tsft
    while [ $Tseg -le $Tsegment ]; do
        ## only write timestamps-file if it is not found already
        if [ "$have_tsFile" != "true" ]; then
	    echo ${tmpTime}" 0" >> $tsFile
        fi
	tmpTime=$(($tmpTime + $Tsft))
	Tseg=$(($Tseg + $Tsft))
    done

    tmpTime=$(($tmpTime + $seggap))
    iSeg=$(($iSeg + 1))
done

## ------------------------------------------------------------
echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Signal"
echo "----------------------------------------------------------------------"
echo
if [ ! -d "$SFTdir" ]; then
    mkdir -p $SFTdir;
fi

FreqStep=`echo $mfd_FreqBand $numFreqBands |  awk '{print $1 / $2}'`
mfd_fBand=`echo $FreqStep $Tsft |  awk '{print ($1 - 1.5 / $2)}'`	## reduce by 1/2 a bin to avoid including last freq-bins

# construct common MFD cmd
mfd_CL_common="--Band=${mfd_fBand} --Freq=$Freq --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --h0=$h0 --cosi=$cosi --ephemYear=05-09 --generationMode=1 --timestampsFile=$tsFile --refTime=$refTime --Tsft=$Tsft --randSeed=1000 --outSingleSFT"

if [ "$sqrtSh" != "0" ]; then
    mfd_CL_common="$mfd_CL_common --noiseSqrtSh=$sqrtSh";
fi

iFreq=1
while [ $iFreq -le $numFreqBands ]; do
    mfd_fi=`echo $mfd_fmin $iFreq $FreqStep | awk '{print $1 + ($2 - 1) * $3}'`

    # for H1:
    SFTname="${SFTdir}${dirsep}H1-${mfd_fi}_${FreqStep}.sft"
    if [ ! -r $SFTname ]; then
        cmdline="$mfd_code $mfd_CL_common --fmin=$mfd_fi --IFO=H1 --outSFTbname=$SFTname"
        echo "$cmdline";
        if ! eval "$cmdline >& /dev/null"; then
            echo "Error.. something failed when running '$mfd_code' ..."
            exit 1
        fi
    else
        echo "SFT '$SFTname' exists already ... reusing it"
    fi

    # for L1:
    SFTname="${SFTdir}${dirsep}L1-${mfd_fi}_${FreqStep}.sft"
    if [ ! -r $SFTname ]; then
        cmdline="$mfd_code $mfd_CL_common --fmin=$mfd_fi --IFO=L1 --outSFTbname=$SFTname";
        echo "$cmdline";
        if ! eval "$cmdline >& /dev/null"; then
            echo "Error.. something failed when running '$mfd_code' ..."
            exit 1
        fi
    else
        echo "SFT '$SFTname' exists already ... reusing it"
    fi

    iFreq=$(( $iFreq + 1 ))

done


echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run PredictFstat (perfect match) "
echo "----------------------------------------------------------------------"
echo
outfile_pfs="avgPFS.dat";

if [ -r "$outfile_pfs" ]; then
    echo "PFS result file '$outfile_pfs' exists already ... reusing it"
    TwoFAvg=$(cat $outfile_pfs)
else
    tmpfile_pfs="__tmp_PFS.dat";
    pfs_CL_common=" --Alpha=$Alpha --Delta=$Delta --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='$SFTfiles' --outputFstat=$tmpfile_pfs --ephemYear=05-09"
    if [ "$sqrtSh" = "0" ]; then
        pfs_CL_common="$pfs_CL_common --SignalOnly";
    fi

    TwoFsum=0
    for ((x=1; x <= $Nsegments; x++)); do
        startGPS=${segs[${x}]}
        endGPS=$(($startGPS + $Tsegment))
        # construct pfs cmdline
        pfs_CL="$pfs_code $pfs_CL_common --minStartTime=$startGPS --maxEndTime=$endGPS"
        echo "$pfs_CL"
        if ! tmp=`eval $pfs_CL`; then
	    echo "Error.. something failed when running '$pfs_code' ..."
	    exit 1
        fi
        resPFS=$(cat ${tmpfile_pfs} | grep 'twoF_expected' | awk '{printf "%g\n", $3}')
        TwoFsum=$(echo $TwoFsum $resPFS | awk '{printf "%.6g", $1 + $2}')
    done

    TwoFAvg=$(echo $TwoFsum $Nsegments | awk '{printf "%.11g", $1 / $2}')
    echo $TwoFAvg > $outfile_pfs
fi

echo
echo "==>   Average 2F: $TwoFAvg"



edat="earth05-09.dat"
sdat="sun05-09.dat"


gct_CL_common="--gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles' --skyGridFile='./$skygridfile' --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segFile --ephemE=$edat --ephemS=$sdat"
if [ "$sqrtSh" = "0" ]; then
    gct_CL_common="$gct_CL_common --SignalOnly";
fi

echo
echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 3: run HierarchSearchGCT using Resampling (perfect match) and segment-list file"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly
outfile_gct1="__tmp_GCT1.dat"

if [ -z "$NORESAMP" ]; then
    cmdline="$gct_code $gct_CL_common --fnameout=$outfile_gct1 --useResamp"
    echo $cmdline
    if ! tmp=`eval $cmdline >& /dev/null`; then
	echo "Error.. something failed when running '$gct_code' ..."
	exit 1
    fi
    topline=$(sort -nr -k6,6 $outfile_gct1 | head -1)
    resGCT1=$(echo $topline | awk '{print $6}')
    freqGCT1=$(echo $topline | awk '{print $1}')
else
    echo
    echo "Not run with resampling."
fi


echo
echo "----------------------------------------------------------------------------------------------------"
echo " STEP 4: run HierarchSearchGCT without Resampling (perfect match) and --tStack and --nStacksMax"
echo "----------------------------------------------------------------------------------------------------"
echo

rm -f checkpoint.cpt # delete checkpoint to start correctly

outfile_gct2="__tmp_GCT2.dat"
cmdline="$gct_code $gct_CL_common --fnameout=$outfile_gct2"
echo $cmdline
if ! eval "$cmdline >& /dev/null"; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

topline=$(sort -nr -k6,6 $outfile_gct2 | head -1)
resGCT2=$(echo $topline  | awk '{print $6}')
freqGCT2=$(echo $topline | awk '{print $1}')

## ---------- compute relative differences and check against tolerance --------------------
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*($1+$2)) }'
echo "awk_reldev = $awk_reldev"

if [ -z "$NORESAMP" ]; then
    reldev1=$(echo $TwoFAvg $resGCT1 | awk "$awk_reldev")
    freqreldev1=$(echo $Freq $freqGCT1 | awk "$awk_reldev")
fi

reldev2=$(echo $TwoFAvg $resGCT2 | awk "$awk_reldev")
freqreldev2=$(echo $Freq $freqGCT2 | awk "$awk_reldev")

# ---------- Check relative deviations against tolerance, report results ----------
retstatus=0
awk_isgtr='{if($1>$2) {print "1"}}'

echo
echo "----------------------------------------------------------------------"
echo "==>  Predicted:		$TwoFAvg	@ $Freq Hz	(Tolerance = $Tolerance)"

if [ -z "$NORESAMP" ]; then
    echo -n "==>  GCT-Resamp:	$resGCT1	@ $freqGCT1 Hz	($reldev1, $freqreldev1"
    fail1=$(echo $reldev1 $Tolerance | awk "$awk_isgtr")
    fail2=$(echo $freqreldev1 $Tolerance | awk "$awk_isgtr")
    if [ "$fail1" -o "$fail2" ]; then
        echo " ... FAILED)"
        retstatus=1
    else
        echo " ... OK)"
    fi
fi

echo -n "==>  GCT-NO-Resamp:	$resGCT2	@ $freqGCT2 Hz 	($reldev2, $freqreldev2"
fail1=$(echo $reldev2 $Tolerance | awk "$awk_isgtr")
fail2=$(echo $freqreldev2 $Tolerance | awk "$awk_isgtr")
if [ "$fail1" -o "$fail2" ]; then
    echo " ... FAILED)"
    retstatus=1
else
    echo " ... OK)"
fi

echo "----------------------------------------------------------------------"

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $skygridfile $outfile_pfs $tmpfile_pfs $outfile_gct1 $outfile_gct2 checkpoint.cpt stderr.log stdout.log $segFile $tsFile
    echo "Cleaned up."
fi

## restore original locale, just in case someone source'd this file
export LC_ALL=$LC_ALL_old

exit $retstatus
