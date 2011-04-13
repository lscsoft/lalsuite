#!/bin/bash


## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
# if [ -n "${srcdir}" ]; then
#     builddir="./";
#     injectdir="../Injections/"
#     fdsdir="../FDS_isolated/"
# else
#     srcdir=.
# fi
builddir="./";
injectdir="../Injections/"
fdsdir="../FDS_isolated/"

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
cfs_code="${fdsdir}lalapps_ComputeFStatistic_v2"
gct_code="${builddir}lalapps_HierarchSearchGCT"

SFTdir="TestSFTs"
SFTfiles="$SFTdir${dirsep}*"


# test if LAL_DATA_PATH has been set ... needed to locate ephemeris-files
if [ -z "$LAL_DATA_PATH" ]; then
    if [ -n "$LALPULSAR_PREFIX" ]; then
	export LAL_DATA_PATH=".:${LALPULSAR_PREFIX}/share/lalpulsar";
    else
	echo
	echo "Need environment-variable LALPULSAR_PREFIX, or LAL_DATA_PATH to be set"
	echo "to your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
	echo "This might indicate an incomplete LAL+LALPULSAR installation"
	echo
	exit 1
    fi
fi

## ---------- fixed parameter of our test-signal -------------
Alpha="3.1"
Delta="-0.5"
h0="1.0"
h0H1=$h0
h0L1=$h0
# h0L1="0.1"
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

mfd_FreqBand="2.0"
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

gct_FreqBand="0.01"
gct_F1dotBand="2.0e-10"
gct_dFreq="0.000002" #"2.0e-6"
gct_dF1dot="1.0e-10"
gct_nCands="100"

edat="earth05-09.dat"
sdat="sun05-09.dat"

noiseSqrtSh="1"

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
midTime="853731034.5"
tsfile="./timestamps.txt"
segFile="./segments.txt"
rm -rf $tsfile $segFile
tmpTime=$startTime
ic1="1"
while [ "$ic1" -le "$Nsegments" ];
do
    t0=$tmpTime
    t1=`echo $t0 $Tsegment | awk '{print $1 + $2}'`
    TspanHours=`echo $Tsegment | awk '{printf "%.7f", $1 / 3600.0 }'`
    NSFT=`echo $Tsegment $Tsft | awk '{print int(2.0 * $1 / $2 + 0.5) }'`
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

# construct MFD cmd for H1:
mfd_CL=" --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --h0=$h0H1 --cosi=$cosi --ephemYear=05-09 --generationMode=1 --timestampsFile=$tsfile --IFO=H1 --refTime=$refTime --Tsft=$Tsft --randSeed=10400"

if [ "$haveNoise" = true ]; then
    mfd_CL="$mfd_CL --noiseSqrtSh=$sqrtSh";
fi

cmdline="$mfd_code $mfd_CL";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

# construct MFD cmd for L1:
mfd_CL=" --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --h0=$h0L1 --cosi=$cosi --ephemYear=05-09 --generationMode=1 --timestampsFile=$tsfile --IFO=L1 --refTime=$refTime --Tsft=$Tsft --randSeed=7001"

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

    # H1 only
    pfs_CL=" --Alpha=$Alpha --Delta=$Delta --h0=$h0H1 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='${SFTfiles}H1*' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS"
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
    pfs_CL=" --Alpha=$Alpha --Delta=$Delta --h0=$h0L1 --cosi=$cosi --psi=$psi --phi0=$phi0 --Freq=$Freq --DataFiles='${SFTfiles}L1*' --outputFstat=$outfile_pfs --ephemYear=05-09 --minStartTime=$startGPS --maxEndTime=$endGPS"
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
    echo
done
TwoFsum=$(echo "scale=6; ${TwoFsum} / ${Nsegments}" | bc);

TwoFsum1=$(echo "scale=6; ${TwoFsum1} / ${Nsegments}" | bc);
TwoFsum2=$(echo "scale=6; ${TwoFsum2} / ${Nsegments}" | bc);
echo
echo "==>   Average 2F: "$TwoFsum"  in H1: "$TwoFsum1"  in L1: "$TwoFsum2


echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run HierarchSearchGCT without Resampling (perfect match), with segment list file"
echo "----------------------------------------------------------------------"
echo

if [ -e "checkpoint.cpt" ]; then
    rm checkpoint.cpt # delete checkpoint to start correctly
fi

outfile_gct="HS_GCT_LV.dat"

gct_CL=" -d1 --fnameout=$outfile_gct --gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles'  --ephemE=$edat --ephemS=$sdat --skyGridFile='./$skygridfile'  --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segFile --outputFX"
if [ "$haveNoise" = false ]; then
    gct_CL="$gct_CL --SignalOnly";
fi

cmdline="$gct_code $gct_CL > >(tee stdout.log) 2> >(tee stderr.log >&2)"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 4: run External ComputeFStat"
echo "----------------------------------------------------------------------"
echo

outfile_cfs="fstat_loudest.dat"
cfs_CL="--DataFiles='$SFTfiles' --outputLoudest='$outfile_cfs' --TwoFthreshold=0.0 --ephemYear=05-09 --refTime=$refTime --Freq=$Freq --f1dot=$f1dot --Alpha=$AlphaSearch --Delta=$DeltaSearch"

cmdline="$cfs_code -v0 $cfs_CL"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

twoFcfs=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
twoFcfs=$(echo "scale=11; ($twoFcfs/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')

cmdline="$cfs_code -v0 $cfs_CL --IFO='H1'"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

twoFcfsX1=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
twoFcfsX1=$(echo "scale=11; ($twoFcfsX1/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')

cmdline="$cfs_code -v0 $cfs_CL --IFO='L1'"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

twoFcfsX2=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
twoFcfsX2=$(echo "scale=11; ($twoFcfsX2/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')

echo "==>  CFS at f="$Freq": F ="$twoFcfs" F1="$twoFcfsX1" F2="$twoFcfsX2

freqGCT=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{printf "%.10f",$1}')
f1dotGCT=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{printf "%.10e",$4}')

twoFcfs_gct="0"
twoFcfs_gctX1="0"
twoFcfs_gctX2="0"
startGPS="0"
endGPS="0"

for ((x=1; x <= $Nsegments; x++))
  do

    startGPS=${segs[${x}]}
    endGPS=$(echo "scale=0; ${startGPS} + ${Tsegment}" | bc | awk '{printf "%.0f",$1}')
    echo "Segment: "$x"  "$startGPS" "$endGPS

    cfs_CL="--DataFiles='$SFTfiles' --outputLoudest='$outfile_cfs' --TwoFthreshold=0.0 --ephemYear=05-09 --refTime=$refTime --Freq=$freqGCT --f1dot=$f1dotGCT --Alpha=$AlphaSearch --Delta=$DeltaSearch --minStartTime=$startGPS --maxEndTime=$endGPS"

    cmdline="$cfs_code -v0 $cfs_CL"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi

    twoFcfs_gct_seg=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_gct=$(echo "scale=6; ${twoFcfs_gct} + ${twoFcfs_gct_seg}" | bc);

    cmdline="$cfs_code -v0 $cfs_CL --IFO='H1'"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi

    twoFcfs_gctX1_seg=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_gctX1=$(echo "scale=6; ${twoFcfs_gctX1} + ${twoFcfs_gctX1_seg}" | bc);

    cmdline="$cfs_code -v0 $cfs_CL --IFO='L1'"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi

    twoFcfs_gctX2_seg=$(echo | sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_gctX2=$(echo "scale=6; ${twoFcfs_gctX2} + ${twoFcfs_gctX2_seg}" | bc);

done

twoFcfs_gct=$(echo "scale=11; ($twoFcfs_gct/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')
twoFcfs_gctX1=$(echo "scale=11; ($twoFcfs_gctX1/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')
twoFcfs_gctX2=$(echo "scale=11; ($twoFcfs_gctX2/$Nsegments)" | bc | awk '{ printf "%.6f",$1}')

echo "==>  CFS at f="$freqGCT": F ="$twoFcfs_gct" F1="$twoFcfs_gctX1" F2="$twoFcfs_gctX2

echo
echo "----------------------------------------------------------------------"
echo " STEP 5: Comparing results"
echo "----------------------------------------------------------------------"
echo

# broad tolerance needed for comparison prediction, since data contains noise . Also do not exit for these, just warn. (if in doubt, run testHS.sh for plain GCT checks)
Tolerance_pred=1e-1
# but different implementations should be very similar
Tolerance_lv=1e-5


# Check relative error in frequency
reldev_pred_freq=$(echo "scale=13; (($Freq - $freqGCT)/$Freq) " | bc | awk '{ if($1>=0) {printf "%.12f",$1} else {printf "%.12f",$1*(-1)}}')
reldev_pred_freq_bins=$(echo "scale=13; (($Freq - $freqGCT)/${gct_dFreq})" | bc | awk '{ if($1>=0) {printf "%.12f",$1} else {printf "%.12f",$1*(-1)}}')

echo "==>  Signal frequency: "$Freq"  Found at: "$freqGCT
if [ `echo $reldev_pred_freq" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT: Rel. dev. in frequency: "$reldev_pred_freq
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT, no Resamp.  Rel. dev. in frequency: "$reldev_pred_freq"   OK."
    echo "     offset as fraction of frequency bins: "$reldev_pred_freq_bins
fi
echo


# Check predicted/recomputed 2F against plain GCT search code output
echo "==>  Predicted              : F ="$TwoFsum" F1="$TwoFsum1" F2="$TwoFsum2
reldev_pred_cfs=$(echo "scale=11; ($TwoFsum - $twoFcfs)/(0.5 * ($TwoFsum + $twoFcfs))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_cfsX1=$(echo "scale=11; ($TwoFsum1 - $twoFcfsX1)/(0.5 * ($TwoFsum1 + $twoFcfsX1))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_cfsX2=$(echo "scale=11; ($TwoFsum2 - $twoFcfsX2)/(0.5 * ($TwoFsum2 + $twoFcfsX2))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_cfsgct=$(echo "scale=11; ($TwoFsum - $twoFcfs_gct)/(0.5 * ($TwoFsum + $twoFcfs_gct))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_cfsgctX1=$(echo "scale=11; ($TwoFsum1 - $twoFcfs_gctX1)/(0.5 * ($TwoFsum1 + $twoFcfs_gctX1))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_cfsgctX2=$(echo "scale=11; ($TwoFsum2 - $twoFcfs_gctX2)/(0.5 * ($TwoFsum2 + $twoFcfs_gctX2))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_cfs_cfsgct=$(echo "scale=11; ($twoFcfs - $twoFcfs_gct)/(0.5 * ($twoFcfs + $twoFcfs_gct))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_cfs_cfsgctX1=$(echo "scale=11; ($twoFcfsX1 - $twoFcfs_gctX1)/(0.5 * ($twoFcfsX1 + $twoFcfs_gctX1))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_cfs_cfsgctX2=$(echo "scale=11; ($twoFcfsX2 - $twoFcfs_gctX2)/(0.5 * ($twoFcfsX2 + $twoFcfs_gctX2))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
echo "==>  CFS at f="$Freq" : F ="$twoFcfs" F1="$twoFcfsX1" F2="$twoFcfsX2
if [ `echo $reldev_pred_cfs" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
    echo "                diff vs pred: "$reldev_pred_cfs",   "$reldev_pred_cfsX1",   "$reldev_pred_cfsX2
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
    echo "                diff vs pred: "$reldev_pred_cfs",   "$reldev_pred_cfsX1",   "$reldev_pred_cfsX2"     OK."
fi
echo "==>  CFS at f="$freqGCT": F ="$twoFcfs_gct" F1="$twoFcfs_gctX1" F2="$twoFcfs_gctX2
if [ `echo $reldev_cfs_cfsgct" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
echo "                diff vs prev: "$reldev_cfs_cfsgct",   "$reldev_cfs_cfsgctX1",   "$reldev_cfs_cfsgctX2
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
echo "                diff vs cfs1: "$reldev_cfs_cfsgct",   "$reldev_cfs_cfsgctX1",   "$reldev_cfs_cfsgctX2"     OK."
fi

twoFGCT=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $6}')
reldev_pred_gct=$(echo "scale=11; ($TwoFsum - $twoFGCT)/(0.5 * ($TwoFsum + $twoFGCT))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
if [ `echo $reldev_pred_gct" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT plain  : F ="$twoFGCT"  (dev vs pred: "$reldev_pred_gct")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
    echo "==>  GCT plain  : F ="$twoFGCT"  (dev vs pred: "$reldev_pred_gct")     OK."
fi

# Check recomputed 2F against previous value and single-detector F against prediction
twoFGCTLV=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $7}')
twoFGCTLVX1=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $8}')
twoFGCTLVX2=$(cat $outfile_gct | sed -e '/%/d;' | sort -nr -k6,6 | head -1 | awk '{print $9}')
reldev_gct_LV=$(echo "scale=11; ($twoFGCT - $twoFGCTLV)/(0.5 * ($twoFGCT + $twoFGCTLV))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_LVX1=$(echo "scale=11; ($TwoFsum1 - $twoFGCTLVX1)/(0.5 * ($TwoFsum1 + $twoFGCTLVX1))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_pred_LVX2=$(echo "scale=11; ($TwoFsum2 - $twoFGCTLVX2)/(0.5 * ($TwoFsum2 + $twoFGCTLVX2))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
maxdev_gct_LV=$(echo | sed '/^ *%%/d;s/%%.*//' $outfile_gct | awk '{if(max==""){max=($6-$7)/($6+$7)}; if($1!="100.123554") { if(($6-$7)/($6+$7)>max) {max=($6-$7)/($6+$7)}; } } END {printf "%.10f",max}' )
maxdevfreq_gct_LV=$(echo | sed '/^ *%%/d;s/%%.*//' $outfile_gct | awk '{if(max==""){max=($6-$7)/($6+$7); maxfreq=$1}; if(($6-$7)/($6+$7)>max) {max=($6-$7)/($6+$7); maxfreq=$1}; } END {printf "%.10f",maxfreq}' )
if [ `echo $reldev_gct_LV" "$Tolerance_lv | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : F ="$twoFGCTLV"  (dev vs plain: "$reldev_gct_LV")"
    echo "                  (max dev: "$maxdev_gct_LV" at freq="$maxdevfreq_gct_LV")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT LV     : F ="$twoFGCTLV"  (dev vs plain: "$reldev_gct_LV")     OK."
    echo "                  (max dev: "$maxdev_gct_LV" at freq="$maxdevfreq_gct_LV")"
fi

if [ `echo $reldev_pred_LVX1" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : F1="$twoFGCTLVX1"  (dev vs pred: "$reldev_pred_LVX1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
    echo "==>  GCT LV     : F1="$twoFGCTLVX1"  (dev vs pred: "$reldev_pred_LVX1")     OK."
fi

if [ `echo $reldev_pred_LVX2" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : F2="$twoFGCTLVX2"  (dev vs pred: "$reldev_pred_LVX2")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
    echo "==>  GCT LV     : F2="$twoFGCTLVX2"  (dev vs pred: "$reldev_pred_LVX2")     OK."
fi

reldev_cfs_LV=$(echo "scale=11; ($twoFcfs - $twoFGCTLV)/(0.5 * ($twoFcfs + $twoFGCTLV))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_cfs_LVX1=$(echo "scale=11; ($twoFcfsX1 - $twoFGCTLVX1)/(0.5 * ($twoFcfsX1 + $twoFGCTLVX1))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
reldev_cfs_LVX2=$(echo "scale=11; ($twoFcfsX2 - $twoFGCTLVX2)/(0.5 * ($twoFcfsX2 + $twoFGCTLVX2))" | bc | awk '{ if($1>=0) {printf "%.10f",$1} else {printf "%.10f",$1*(-1)}}')
if [ `echo $reldev_cfs_LV" "$Tolerance_pred | awk '{if($1>$2) {print "1"}}'` ];then
echo "                diff vs cfs1: "$reldev_cfs_LV",   "$reldev_cfs_LVX1",   "$reldev_cfs_LVX2
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
#     exit 2
else
echo "                diff vs cfs1: "$reldev_cfs_LV",   "$reldev_cfs_LVX1",   "$reldev_cfs_LVX2"     OK."
fi