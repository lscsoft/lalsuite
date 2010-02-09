#!/bin/sh

RESAMP="0"

## allow 'make test' to work from builddir != srcdir
if [ -n "${srcdir}" ]; then
    builddir="./";
    injectdir="../Injections/"
    fdsdir="../FDS_isolated/"
else
    srcdir=.
fi

##---------- names of codes and input/output files
mfd_code="${injectdir}lalapps_Makefakedata_v4"
pfs_code="${fdsdir}lalapps_PredictFStat"
gct_code="${builddir}lalapps_HierarchSearchGCT"

SFTdir="./TestSFTs"

# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    if [ -n "$LALPULSAR_PREFIX" ]; then
	export LAL_DATA_PATH="${LALPULSAR_PREFIX}/share/lalpulsar";
    else
	echo
	echo "Need environment-variable LAL_PREFIX, or LAL_DATA_PATH to be set"
	echo "to your ephemeris-directory (e.g. /usr/local/share/lal)"
	echo "This might indicate an incomplete LAL installation"
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
Freq="100.12345"
f1dot="-1e-9"

AlphaSearch=$Alpha
DeltaSearch=$Delta

## Produce skygrid file for the search
skygridfile="./tmpskygridfile.dat"
echo $AlphaSearch" "$DeltaSearch > $skygridfile

mfd_FreqBand="2.0"
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

gct_FreqBand="0.02"
gct_F1dotBand="2.0e-10"
gct_dFreq="2.0e-6"
gct_dF1dot="1.0e-10"
gct_nCands="100000"

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
Nsegments="131"
seggap=$(echo "scale=0; ${Tsegment} * 1.12345" | bc) 
tsfile="timestampsTEST.txt"
rm -rf $tsfile
tmpTime=$startTime
ic1="1"
while [ "$ic1" -le "$Nsegments" ];
do
    segs[${ic1}]=$tmpTime # save seg's beginning for later use
    echo "Segment: "$ic1" of "$Nsegments"   GPS start time: "${segs[${ic1}]}
    
    ic2=$Tsft
    while [ "$ic2" -le "$Tsegment" ];
    do
	echo ${tmpTime}" 0" >> $tsfile
	tmpTime=$(echo "scale=0; ${tmpTime} + ${Tsft}" | bc)
	ic2=$(echo "scale=0; ${ic2} + ${Tsft}" | bc)
    done
    
    tmpTime=$(echo "scale=0; ${tmpTime} + ${seggap}" | bc | awk '{printf "%.0f",$1}')
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

IFO=LHO

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

# construct MFD cmd:
mfd_CL=" --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --h0=$h0 --cosi=$cosi --ephemYear=05-09 --generationMode=1 --timestampsFile=$tsfile --IFO=$IFO --refTime=$refTime --Tsft=$Tsft "

if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run PredictFstat (perfect match) "
echo "----------------------------------------------------------------------"
echo

TwoFsum="0"
for ((x=1; x <= $Nsegments; x++))
  do
    outfile_pfs="__tmp_PFS.dat";
    
    startGPS=${segs[${x}]}
    endGPS=$(echo "scale=0; ${startGPS} + ${Tsegment}" | bc | awk '{printf "%.0f",$1}')
    echo "Segment: "$x"  "$startGPS" "$endGPS
    # construct pfs cmd
    pfs_CL=" --Alpha=$Alpha --Delta=$Delta --IFO=$IFO --h0=$h0 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='$SFTdir/*' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS --SignalOnly"

    cmdline="$pfs_code $pfs_CL"
    echo $cmdline
    if ! tmp=`eval $cmdline`; then
	echo "Error.. something failed when running '$pfs_code' ..."
	exit 1
    fi
    resPFS=`echo $tmp | awk '{printf "%g", $1}'`
    TwoFsum=$(echo "scale=6; ${TwoFsum} + ${resPFS}" | bc);
    #echo " Segment: "$x"   2F: "$resPFS"   Sum of 2F: "$TwoFsum
    echo
done
TwoFsum=$(echo "scale=6; ${TwoFsum} / ${Nsegments}" | bc);
echo
echo "==>   Average 2F: "$TwoFsum


edat=${LAL_DATA_PATH}"/earth05-09.dat"
sdat=${LAL_DATA_PATH}"/sun05-09.dat"

echo
echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run HierarchSearchGCT using Resampling (perfect match)"
echo "----------------------------------------------------------------------"
echo

outfile_gct1="__tmp_GCT1.dat"
                                                                                           
gct_CL=" --useResamp --SignalOnly --fnameout=$outfile_gct1 --gridType=3 --tStack=$Tsegment --nCand1=$gct_nCands --nStacksMax=$Nsegments --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTdir/*'  --ephemE=$edat --ephemS=$sdat --skyGridFile='$skygridfile' --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime "

cmdline="$gct_code $gct_CL"
echo $cmdline
if [ -z "$RESAMP" ]; then
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

outfile_gct2="__tmp_GCT2.dat"

gct_CL=" --SignalOnly --fnameout=$outfile_gct2 --gridType=3 --tStack=$Tsegment --nCand1=$gct_nCands --nStacksMax=$Nsegments --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTdir/*'  --ephemE=$edat --ephemS=$sdat --skyGridFile='$skygridfile'  --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime "

cmdline="$gct_code $gct_CL"
echo $cmdline
if ! tmp=`eval $cmdline`; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi
resGCT2=$(cat $outfile_gct2 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $6}')
freqGCT2=$(cat $outfile_gct2 | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $1}')

if [ -z "$RESAMP" ]; then
reldev1=$(echo "scale=5; ($TwoFsum - $resGCT1)/(0.5 * ($TwoFsum + $resGCT1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
freqreldev1=$(echo "scale=8; (($Freq - $freqGCT1)/$Freq) " | bc | awk '{ if($1>=0) {printf "%.6f",$1} else {printf "%.6f",$1*(-1)}}')
reldev3=$(echo "scale=5; ($resGCT1 - $resGCT2)/(0.5 * ($resGCT2 + $resGCT1))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
fi

reldev2=$(echo "scale=5; ($TwoFsum - $resGCT2)/(0.5 * ($TwoFsum + $resGCT2))" | bc | awk '{ if($1>=0) {printf "%.4f",$1} else {printf "%.4f",$1*(-1)}}')
freqreldev2=$(echo "scale=8; (($Freq - $freqGCT2)/$Freq) " | bc | awk '{ if($1>=0) {printf "%.6f",$1} else {printf "%.6f",$1*(-1)}}')


echo
echo "----------------------------------------------------------------------"
echo "==>  Predicted:      "$TwoFsum

# Check predicted 2F against search code output

if [ -z "$RESAMP" ]; then
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

if [ -z "$RESAMP" ]; then
if [ `echo $reldev3" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT, Resamp vs. no-Resamp:     "$reldev3
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, Resamp vs. no-Resamp:    "$reldev3"      OK." 
fi
fi

echo
echo "==>  Signal frequency: "$Freq"  Found at: "$freqGCT1  

# Check relative error in frequency
if [ -z "$RESAMP" ]; then
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
fi

echo "----------------------------------------------------------------------"


## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $skygridfile $tsfile $outfile_pfs $outfile_gct1 $outfile_gct2
    echo "Cleaned up."
fi

