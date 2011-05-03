#!/bin/bash

#NOCLEANUP="1"

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

## make targeted search in sky coordinates
AlphaSearch=$Alpha
DeltaSearch=$Delta

## generate data with noise
noiseSqrtSh="1.0"

## Produce skygrid file for the search
skygridfile="tmpskygridfile.dat"
echo $AlphaSearch" "$DeltaSearch > $skygridfile

## frequency bands for makefakedata, HSGCT
mfd_FreqBand="2.0"
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

gct_FreqBand="0.01"
gct_F1dotBand="2.0e-10"
gct_dFreq="0.000002" #"2.0e-6"
gct_dF1dot="1.0e-10"
gct_nCands="1000"  # number of candidates in output

edat="earth09-11.dat"
sdat="sun09-11.dat"
RngMedWindow=101 # running median window, needs to be equal for HSGCT and CFS

## --------- Generate fake data set time stamps -------------
echo "----------------------------------------------------------------------"
echo " STEP 0: Generating fake data time-stamps file "
echo "----------------------------------------------------------------------"
echo

Tsft="1800"
startTime="952443819"
refTime="962999869"
Tsegment="90000"
Nsegments="14"
## insert gaps between segments
seggap=$(echo $Tsegment | awk '{printf "%f", $1 * 1.12345}');
# midTime="853731034.5"
tsfileH1="./timestampsH1.dat"  # for makefakedata
tsfileL1="./timestampsL1.dat"  # for makefakedata
segfile="./segments.dat"   # for HSGCT
rm -rf $tsfileH1 $tsfileL1 $segfile
tmpTime=$startTime
ic1="1"
while [ "$ic1" -le "$Nsegments" ];
do
    t0=$tmpTime
    t1=`echo $t0 $Tsegment | awk '{print $1 + $2}'`
    TspanHours=`echo $Tsegment | awk '{printf "%.7f", $1 / 3600.0 }'`
    NSFT=`echo $Tsegment $Tsft | awk '{print int(2.0 * $1 / $2 + 0.5) }'`
    if [ `echo $ic1 | awk '{if($1==1) {print "1"}}'` ];then
       NSFT=`echo $Tsegment $Tsft | awk '{print int(1.0 * $1 / $2 + 0.5) }'`
    fi
    if [ `echo $ic1" "$Nsegments | awk '{if($1==$2) {print "1"}}'` ];then
       NSFT=`echo $Tsegment $Tsft | awk '{print int(1.0 * $1 / $2 + 0.5) }'`
    fi
    echo "$t0 $t1 $TspanHours $NSFT" >> $segfile
    segs[${ic1}]=$tmpTime # save seg's beginning for later use
    echo "Segment: "$ic1" of "$Nsegments"   GPS start time: "${segs[${ic1}]}

    ic2=$Tsft
    while [ "$ic2" -le "$Tsegment" ];
    do
        if [ `echo $ic1 | awk '{if($1>1) {print "1"}}'` ];then
            echo ${tmpTime}" 0" >> $tsfileH1
        fi
        if [ `echo $ic1" "$Nsegments | awk '{if($1<$2) {print "1"}}'` ];then
            echo ${tmpTime}" 0" >> $tsfileL1
        fi
        tmpTime=$(echo $tmpTime $Tsft | awk '{printf "%.0f", $1 + $2}');
        ic2=$(echo $ic2 $Tsft | awk '{printf "%.0f", $1 + $2}');
    done

    tmpTime=$(echo $tmpTime $seggap | awk '{printf "%.0f", $1 + $2}');
    ic1=$(echo $ic1 | awk '{printf "%.0f", $1 + 1}');
done


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

## construct MFD cmdline
mfd_CL=" --fmin=$mfd_fmin --Band=$mfd_FreqBand --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --cosi=$cosi --ephemYear=09-11 --generationMode=1 --refTime=$refTime --Tsft=$Tsft --noiseSqrtSh=$noiseSqrtSh"

## detector H1
cmdline="$mfd_code $mfd_CL --IFO=H1 --h0=$h0H1  --timestampsFile=$tsfileH1  --randSeed=1000";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## detector L1:
cmdline="$mfd_code $mfd_CL --IFO=L1 --h0=$h0L1 --timestampsFile=$tsfileL1  --randSeed=1001";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi


echo
echo "----------------------------------------------------------------------"
echo " STEP 2: run HierarchSearchGCT"
echo "----------------------------------------------------------------------"
echo

if [ -e "checkpoint.cpt" ]; then
    rm checkpoint.cpt # delete checkpoint to start correctly
fi

outfile_gct="HS_GCT_LV.dat"

gct_CL=" -d1 --fnameout=$outfile_gct --gridType1=3 --nCand1=$gct_nCands --skyRegion='allsky' --Freq=$Freq --DataFiles='$SFTfiles'  --ephemE=$edat --ephemS=$sdat --skyGridFile='./$skygridfile'  --printCand1 --semiCohToplist --df1dot=$gct_dF1dot --f1dot=$f1dot --f1dotBand=$gct_F1dotBand --dFreq=$gct_dFreq --FreqBand=$gct_FreqBand --refTime=$refTime --segmentList=$segfile --outputFX --blocksRngMed=$RngMedWindow --outputTiming='./timing.dat'"

cmdline="$gct_code $gct_CL > >(tee stdout.log) 2> >(tee stderr.log >&2)"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$gct_code' ..."
    exit 1
fi

## get recovered signal frequency from GCT
freqGCT=$(sed -e '/%/d;' $outfile_gct | sort -nr -k6,6 | head -1 | awk '{printf "%.16f",$1}')
f1dotGCT=$(sed -e '/%/d;' $outfile_gct | sort -nr -k6,6 | head -1 | awk '{printf "%.10e",$4}')

echo
echo "----------------------------------------------------------------------"
echo " STEP 3: run External ComputeFStat "
echo "----------------------------------------------------------------------"
echo

outfile_cfs="fstat_loudest.dat"
## initialise summation variables
twoFcfs="0"
twoFcfsH1="0"
twoFcfsL1="0"

for ((x=1; x <= $Nsegments; x++))
  do

    ## find GPS times of segment
    startGPS=${segs[${x}]}
    endGPS=$(echo $startGPS $Tsegment | awk '{printf "%.0f", $1 + $2}');
    echo "Segment: "$x"  "$startGPS" "$endGPS

    ## construct ComputeFStatistic command lines
    cfs_CL="--DataFiles='$SFTfiles' --outputLoudest='$outfile_cfs' --TwoFthreshold=0.0 --ephemYear=09-11 --refTime=$refTime --Alpha=$AlphaSearch --Delta=$DeltaSearch --Freq=$freqGCT --f1dot=$f1dotGCT --minStartTime=$startGPS --maxEndTime=$endGPS --RngMedWindow=$RngMedWindow"

    ## multi-IFO
    cmdline="$cfs_code -v0 $cfs_CL"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi
    twoFcfs_seg=$(sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs=$(echo $twoFcfs $twoFcfs_seg | awk '{printf "%.6f", $1 + $2}');

    ## detector H1
    if [ `echo $x | awk '{if($1>1) {print "1"}}'` ];then
      cmdline="$cfs_code -v0 $cfs_CL --IFO='H1'"
      echo $cmdline
      if ! eval $cmdline; then
        echo "Error.. something failed when running '$cfs_code' ..."
        exit 1
      fi
      twoFcfs_seg=$(sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
      twoFcfsH1=$(echo $twoFcfsH1 $twoFcfs_seg | awk '{printf "%.6f", $1 + $2}');
   fi

    ## detector L1
    if [ `echo $x" "$Nsegments | awk '{if($1<$2) {print "1"}}'` ];then
      cmdline="$cfs_code -v0 $cfs_CL --IFO='L1'"
      echo $cmdline
      if ! eval $cmdline; then
        echo "Error.. something failed when running '$cfs_code' ..."
        exit 1
      fi
      twoFcfs_seg=$(sed 's/\;//' $outfile_cfs | awk '{if($1=="twoF"){printf "%.6f",$3}}')
      twoFcfsL1=$(echo $twoFcfsL1 $twoFcfs_seg | awk '{printf "%.6f", $1 + $2}');
    fi

done

## get averages
twoFcfs=$(echo $twoFcfs $Nsegments             | awk '{printf "%.6f", $1/$2}');
twoFcfsH1=$(echo $twoFcfsH1 $Nsegments         | awk '{printf "%.6f", $1/($2-1)}');
twoFcfsL1=$(echo $twoFcfsL1 $Nsegments         | awk '{printf "%.6f", $1/($2-1)}');

echo
echo "----------------------------------------------------------------------"
echo " STEP 4: Comparing results"
echo "----------------------------------------------------------------------"
echo

Tolerance=1e-3


## check relative error in frequency (injection vs. GCT)
reldev_freq=$(echo $Freq $freqGCT                 | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/$2} else {printf "%.12f", (-1)*($1-$2)/$2}}');
reldev_freq_bins=$(echo $Freq $freqGCT $gct_dFreq | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/$3} else {printf "%.12f", (-1)*($1-$2)/$3}}');

echo "==>  Signal frequency: "$Freq"  Found at: "$freqGCT
if [ `echo $reldev_freq" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT: Rel. dev. in frequency: "$reldev_freq
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT: Rel. dev. in frequency:          "$reldev_freq"   OK."
    echo "     offset as fraction of frequency bins: "$reldev_freq_bins
fi
echo


## get Fstats values from GCT output file
twoFGCT=$(sed -e '/%/d;'  $outfile_gct | sort -nr -k6,6 | head -1 | awk '{print $6}')
twoFGCTLV=$(sed -e '/%/d;'  $outfile_gct | sort -nr -k6,6 | head -1 | awk '{print $7}')
twoFGCTLVH1=$(sed -e '/%/d;'  $outfile_gct | sort -nr -k6,6 | head -1 | awk '{print $9}')
twoFGCTLVL1=$(sed -e '/%/d;'  $outfile_gct | sort -nr -k6,6 | head -1 | awk '{print $8}')
## compute relative errors of Fstats from CFS, GTC plain and GTCLV
reldev_cfs_gct=$(echo $twoFcfs $twoFGCT        | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
reldev_gct_LV=$(echo $twoFGCT $twoFGCTLV       | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
reldev_cfs_LV=$(echo $twoFcfs $twoFGCTLV       | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
reldev_cfs_LVH1=$(echo $twoFcfsH1 $twoFGCTLVH1 | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
reldev_cfs_LVL1=$(echo $twoFcfsL1 $twoFGCTLVL1 | awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
## get maximum deviations (over all candidates) of freq, Fstat between plain GCT and LV
maxdev_gct_LV=$(sed '/^%.*/d' $outfile_gct | awk 'BEGIN { maxDiff = 0; } { relDiff=($6-$7)/($6+$7); if( relDiff^2 > maxDiff^2 ) {maxDiff=relDiff}; } END {printf "%.12f", maxDiff}' )
maxdevfreq_gct_LV=$(sed '/^%.*/d' $outfile_gct | awk 'BEGIN { maxDiff = 0; maxFreq = 0;} { relDiff=($6-$7)/($6+$7); if( relDiff^2 > maxDiff^2 ) {maxDiff=relDiff; maxFreq=$1}; } END {printf "%.12f", maxFreq}' )


echo "==>  CFS at f="$freqGCT" : F="$twoFcfs" FH1="$twoFcfsH1" FL1="$twoFcfsL1

## check plain GCT search code output 2F against externally computed 2F
if [ `echo $reldev_cfs_gct" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT plain  : F  ="$twoFGCT"  (diff vs cfs:   "$reldev_cfs_gct")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT plain  : F  ="$twoFGCT"  (diff vs cfs:   "$reldev_cfs_gct")     OK."
fi

## check LV-recomputed 2F against plain GCT value
if [ `echo $reldev_gct_LV" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : F  ="$twoFGCTLV"  (diff vs plain: "$reldev_gct_LV")"
    echo "                                   (maximum diff:  "$maxdev_gct_LV" at freq="$maxdevfreq_gct_LV")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT LV     : F  ="$twoFGCTLV"  (diff vs plain: "$reldev_gct_LV")     OK."
    echo "                                    (maximum diff:  "$maxdev_gct_LV" at freq="$maxdevfreq_gct_LV")"
fi

## check recomputed 2F against external CFS
if [ `echo $reldev_cfs_LV" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
echo "                                    (diff vs cfs:   "$reldev_cfs_LV")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
echo "                                    (diff vs cfs:   "$reldev_cfs_LV")     OK."
fi

## check single-detector F1 against external CFS (at GCT freq)
if [ `echo $reldev_cfs_LVH1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : FH1="$twoFGCTLVH1"   (diff vs cfs:   "$reldev_cfs_LVH1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT LV     : FH1="$twoFGCTLVH1"   (diff vs cfs:   "$reldev_cfs_LVH1")     OK."
fi

## check single-detector F2 against external CFS (at GCT freq)
if [ `echo $reldev_cfs_LVL1" "$Tolerance | awk '{if($1>$2) {print "1"}}'` ];then
    echo "==>  GCT LV     : FL1="$twoFGCTLVL1"   (diff vs cfs:   "$reldev_cfs_LVL1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "==>  GCT LV     : FL1="$twoFGCTLVL1"   (diff vs cfs:   "$reldev_cfs_LVL1")     OK."
fi


echo "----------------------------------------------------------------------"

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $skygridfile $tsfileH1 $tsfileL1 $segfile $outfile_gct $outfile_cfs checkpoint.cpt stderr.log stdout.log timing.dat
    echo "Cleaned up."
fi
