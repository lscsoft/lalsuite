#!/bin/bash

#NOCLEANUP="1"

## take user-arguments
extra_args="$@"

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
cfs_code="${fdsdir}lalapps_ComputeFStatistic_v2"

SFTdir="testCFSv2_singleF_sfts"
SFTfiles="$SFTdir${dirsep}*"

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

## ---------- fixed parameter of our test-signal -------------
Alpha="1.42"
Delta="-0.5"
h0="1.0"
cosi="-0.3"
psi="0.6"
phi0="1.5"
Freq="100.5"
f1dot="0.0"

## generate data with noise
noiseSqrtSh="1.0"
Tsft="1800"
startTime="952443819"
duration="90000"
refTime="962999869"

## frequency band for makefakedata
mfd_FreqBand="0.5"
mfd_fmin=$(echo $Freq $mfd_FreqBand | LC_ALL=C awk '{printf "%g", $1 - $2 / 2.0}');

edat="earth09-11.dat"
sdat="sun09-11.dat"

cfs_FreqBand="0.1";
cfs_fmin=$(echo $Freq $cfs_FreqBand | LC_ALL=C awk '{printf "%8f", $1 - $2 / 2.0}');
cfs_toplist_cands="1000"

SFTdir="TestSFTs"
SFTfiles="$SFTdir${dirsep}*"

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
mfd_CL=" --Tsft=$Tsft --startTime=$startTime --duration=$duration --fmin=$mfd_fmin --Band=$mfd_FreqBand --h0=$h0 --Freq=$Freq --outSFTbname=$SFTdir --f1dot=$f1dot --Alpha=$Alpha --Delta=$Delta --psi=$psi --phi0=$phi0 --cosi=$cosi --generationMode=1 --refTime=$refTime --noiseSqrtSh=$noiseSqrtSh"

## detector H1
cmdline="$mfd_code $mfd_CL --IFO=H1 --randSeed=1000";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

## detector L1:
cmdline="$mfd_code $mfd_CL --IFO=L1 --randSeed=1001";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi


timing_awk='BEGIN { timingsum = 0; counter=0; } { timingsum=timingsum+$9; counter=counter+1; } END {printf "%.3g", timingsum/counter}'

echo
echo "----------------------------------------------------------------------"
echo " STEP 2a: run standard ComputeFStatistic_v2"
echo "----------------------------------------------------------------------"
echo

outfile_cfs_loudest="fstat_loudest.dat"
outfile_cfs_all="fstat_all.dat"
timingsfile="cfs_timing.dat"

    ## construct ComputeFStatistic command lines
    cfs_CL=" --DataFiles='$SFTfiles' --TwoFthreshold=0.0 --Alpha=$Alpha --Delta=$Delta --Freq=$cfs_fmin --FreqBand=$cfs_FreqBand --clusterOnScanline=2"

    ## multi-IFO
    cmdline="$cfs_code $cfs_CL --outputFstat='$outfile_cfs_all' --outputLoudest='$outfile_cfs_loudest' --outputTiming='$timingsfile'"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi
    twoFcfs_multi=$(sed 's/\;//' $outfile_cfs_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_multi_all=$(sed -e '/%/d;'  $outfile_cfs_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')

    ## detector H1
    cmdline="$cfs_code  $cfs_CL --outputFstat='$outfile_cfs_all' --outputLoudest='$outfile_cfs_loudest' --IFO='H1'"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi
    twoFcfs_H1=$(sed 's/\;//' $outfile_cfs_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_H1_all=$(sed -e '/%/d;'  $outfile_cfs_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')

    ## detector L1
    cmdline="$cfs_code $cfs_CL --outputFstat='$outfile_cfs_all' --outputLoudest='$outfile_cfs_loudest' --IFO='L1'"
    echo $cmdline
    if ! eval $cmdline; then
      echo "Error.. something failed when running '$cfs_code' ..."
      exit 1
    fi
    twoFcfs_L1=$(sed 's/\;//' $outfile_cfs_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
    twoFcfs_L1_all=$(sed -e '/%/d;'  $outfile_cfs_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')
    timing_plain=$(sed '/^%.*/d' $timingsfile | LC_ALL=C awk "$timing_awk")


echo
echo "----------------------------------------------------------------------"
echo " STEP 2b: run standard ComputeFStatistic_v2 with toplist"
echo "----------------------------------------------------------------------"
echo

outfile_cfs_toplist_loudest="fstat_toplist_loudest.dat"
outfile_cfs_toplist_all="fstat_toplist_all.dat"
timingsfile_toplist="cfs_timing_toplist.dat"

     cmdline="$cfs_code $cfs_CL --outputFstat='$outfile_cfs_toplist_all' --outputLoudest='$outfile_cfs_toplist_loudest'  --outputTiming='$timingsfile_toplist' --NumCandidatesToKeep=$cfs_toplist_cands"
     echo $cmdline
     if ! eval $cmdline; then
       echo "Error.. something failed when running '$cfs_code' ..."
       exit 1
     fi
     twoFcfs_toplist_multi=$(sed 's/\;//' $outfile_cfs_toplist_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
     twoFcfs_toplist_multi_all=$(sed -e '/%/d;'  $outfile_cfs_toplist_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')
     timing_toplist=$(sed '/^%.*/d' $timingsfile_toplist | LC_ALL=C awk "$timing_awk" )

echo
echo "----------------------------------------------------------------------"
echo " STEP 3a: run ComputeFStatistic_v2 with single-IFO F-stats"
echo "----------------------------------------------------------------------"
echo

outfile_cfs_singleF_loudest="fstat_singleF_loudest.dat"
outfile_cfs_singleF_all="fstat_singleF_all.dat"
timingsfile_singleF="cfs_timing_singleF.dat"

     cmdline="$cfs_code $cfs_CL --outputSingleF --outputFstat='$outfile_cfs_singleF_all' --outputLoudest='$outfile_cfs_singleF_loudest' --outputTiming='$timingsfile_singleF'"
     echo $cmdline
     if ! eval $cmdline; then
       echo "Error.. something failed when running '$cfs_code' ..."
       exit 1
     fi
     twoFcfs_singleF_multi=$(sed 's/\;//' $outfile_cfs_singleF_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
     twoFcfs_singleF_H1=$(sed 's/\;//' $outfile_cfs_singleF_loudest | LC_ALL=C awk '{if($1=="twoF0"){printf "%.6f",$3}}')
     twoFcfs_singleF_L1=$(sed 's/\;//' $outfile_cfs_singleF_loudest | LC_ALL=C awk '{if($1=="twoF1"){printf "%.6f",$3}}')
     twoFcfs_singleF_multi_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')
     twoFcfs_singleF_H1_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$8}')
     twoFcfs_singleF_L1_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$9}')
     timing_singleF=$(sed '/^%.*/d' $timingsfile_singleF | LC_ALL=C awk "$timing_awk" )

echo
echo "----------------------------------------------------------------------"
echo " STEP 3b: run ComputeFStatistic_v2 with single-IFO F-stats and toplist"
echo "----------------------------------------------------------------------"
echo

outfile_cfs_singleF_toplist_loudest="fstat_singleF_toplist_loudest.dat"
outfile_cfs_singleF_toplist_all="fstat_singleF_toplist_all.dat"
timingsfile_singleF_toplist="cfs_timing_singleF_toplist.dat"

     cmdline="$cfs_code  $cfs_CL --outputSingleF --outputFstat='$outfile_cfs_singleF_toplist_all' --outputLoudest='$outfile_cfs_singleF_toplist_loudest'  --outputTiming='$timingsfile_singleF_toplist' --NumCandidatesToKeep=$cfs_toplist_cands"
     echo $cmdline
     if ! eval $cmdline; then
       echo "Error.. something failed when running '$cfs_code' ..."
       exit 1
     fi
     twoFcfs_singleF_toplist_multi=$(sed 's/\;//' $outfile_cfs_singleF_toplist_loudest | LC_ALL=C awk '{if($1=="twoF"){printf "%.6f",$3}}')
     twoFcfs_singleF_toplist_H1=$(sed 's/\;//' $outfile_cfs_singleF_toplist_loudest | LC_ALL=C awk '{if($1=="twoF0"){printf "%.6f",$3}}')
     twoFcfs_singleF_toplist_L1=$(sed 's/\;//' $outfile_cfs_singleF_toplist_loudest | LC_ALL=C awk '{if($1=="twoF1"){printf "%.6f",$3}}')
     twoFcfs_singleF_toplist_multi_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_toplist_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$7}')
     twoFcfs_singleF_toplist_H1_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_toplist_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$8}')
     twoFcfs_singleF_toplist_L1_all=$(sed -e '/%/d;'  $outfile_cfs_singleF_toplist_all | sort -nr -k7,7 | head -1 | LC_ALL=C awk '{printf "%6f",$9}')
     timing_singleF_toplist=$(sed '/^%.*/d' $timingsfile_singleF_toplist | LC_ALL=C awk "$timing_awk" )


echo
echo "----------------------------------------------------------------------"
echo " STEP 4: Comparing results"
echo "----------------------------------------------------------------------"
echo

Tolerance=1e-5
echo "standard CFS_v2 :"
echo "Loudest candidate        : 2F_multi = "$twoFcfs_multi" , 2F_H1 = "$twoFcfs_H1" , 2F_L1 = "$twoFcfs_L1
echo "Highest signal from list : 2F_multi = "$twoFcfs_multi_all" , 2F_H1 = "$twoFcfs_H1_all" , 2F_L1 = "$twoFcfs_L1_all
echo "With toplist: loudest cand 2F_multi = "$twoFcfs_toplist_multi", highest from list : 2F_multi = "$twoFcfs_toplist_multi_all

echo "Rerun with --singleF=TRUE :"
reldev_F_multi=$(echo $twoFcfs_multi $twoFcfs_singleF_multi | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_F_multi" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_multi = "$twoFcfs_singleF_multi"  (diff vs plain: "$reldev_F_multi")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_multi = "$twoFcfs_singleF_multi"  (diff vs plain: "$reldev_F_multi")     OK."
fi
reldev_F_H1=$(echo $twoFcfs_H1 $twoFcfs_singleF_H1 | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_F_H1" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_H1    = "$twoFcfs_singleF_H1"  (diff vs plain: "$reldev_F_H1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_H1    = "$twoFcfs_singleF_H1"  (diff vs plain: "$reldev_F_H1")     OK."
fi
reldev_F_L1=$(echo $twoFcfs_L1 $twoFcfs_singleF_L1 | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_F_L1" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_L1    = "$twoFcfs_singleF_L1"  (diff vs plain: "$reldev_F_L1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_L1    = "$twoFcfs_singleF_L1"  (diff vs plain: "$reldev_F_L1")     OK."
fi
echo "Highest signal from list : 2F_multi = "$twoFcfs_singleF_multi_all" , 2F_H1 = "$twoFcfs_singleF_H1_all" , 2F_L1 = "$twoFcfs_singleF_L1_all

echo "Second rerun with --singleF=TRUE and toplist:"
reldev_toplist_F_multi=$(echo $twoFcfs_multi $twoFcfs_singleF_toplist_multi | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_toplist_F_multi" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_multi = "$twoFcfs_singleF_toplist_multi"  (diff vs plain: "$reldev_toplist_F_multi")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_multi = "$twoFcfs_singleF_toplist_multi"  (diff vs plain: "$reldev_toplist_F_multi")     OK."
fi
reldev_toplist_F_H1=$(echo $twoFcfs_H1 $twoFcfs_singleF_toplist_H1 | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_toplist_F_H1" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_H1    = "$twoFcfs_singleF_toplist_H1"  (diff vs plain: "$reldev_toplist_F_H1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_H1    = "$twoFcfs_singleF_toplist_H1"  (diff vs plain: "$reldev_toplist_F_H1")     OK."
fi
reldev_toplist_F_L1=$(echo $twoFcfs_L1 $twoFcfs_singleF_toplist_L1 | LC_ALL=C awk '{ if(($1-$2)>=0) {printf "%.12f", ($1-$2)/(0.5*($1+$2))} else {printf "%.12f", (-1)*($1-$2)/(0.5*($1+$2))}}');
if [ `echo $reldev_toplist_F_L1" "$Tolerance | LC_ALL=C awk '{if($1>$2) {print "1"}}'` ];then
    echo "     2F_L1    = "$twoFcfs_singleF_toplist_L1"  (diff vs plain: "$reldev_toplist_F_L1")"
    echo "OUCH... results differ by more than tolerance limit. Something might be wrong..."
    exit 2
else
    echo "     2F_L1    = "$twoFcfs_singleF_toplist_L1"  (diff vs plain: "$reldev_toplist_F_L1")     OK."
fi
echo "Highest signal from list : 2F_multi = "$twoFcfs_singleF_toplist_multi_all" , 2F_H1 = "$twoFcfs_singleF_toplist_H1_all" , 2F_L1 = "$twoFcfs_singleF_toplist_L1_all


echo
echo "----------------------------------------------------------------------"
echo " STEP 5: Timings"
echo "----------------------------------------------------------------------"
echo

echo "Timings for search targeted in alpha, delta, f1dot and with freqband="$cfs_FreqBand" , with "$cfs_toplist_cands" candidates in the toplists"
    echo "total time CFS_v2 plain           : "$timing_plain"s"
    echo "total time CFS_v2 toplist         : "$timing_toplist"s"
    echo "total time CFS_v2 singleF         : "$timing_singleF"s"
    echo "total time CFS_v2 singleF toplist : "$timing_singleF_toplist"s"

echo "----------------------------------------------------------------------"

# clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir $outfile_cfs_loudest $outfile_cfs_all $outfile_cfs_toplist_loudest $outfile_cfs_toplist_all $outfile_cfs_singleF_loudest $outfile_cfs_singleF_all $outfile_cfs_singleF_toplist_loudest $outfile_cfs_singleF_toplist_all $timingsfile $timingsfile_toplist $timingsfile_singleF $timingsfile_singleF_toplist checkpoint.cpt stderr.log stdout.log timing.dat
    echo "Cleaned up."
fi
