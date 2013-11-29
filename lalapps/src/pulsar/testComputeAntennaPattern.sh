#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## take user-arguments
extra_args="$@"

## allow 'make test' to work from builddir != srcdir
if [ -z "${srcdir}" ]; then
    srcdir=`dirname $0`
fi

builddir="./";
injectdir="./Injections/"
fdsdir="./FDS_isolated/"

## check for LALPulsar data directory
if [ -z "${LALPULSAR_DATADIR}" ]; then
    echo "Need environment-variable LALPULSAR_DATADIR to be set to"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    exit 1
fi

##---------- names of codes
cap_code="${builddir}lalapps_ComputeAntennaPattern"
pds_code="${builddir}lalapps_PrintDetectorState"
mfd_code="${injectdir}lalapps_Makefakedata_v4 -E ${LALPULSAR_DATADIR}"
pfs_code="${fdsdir}lalapps_PredictFStat"

tolerance=1e-3
Tsft=1800

## awk commands needed for testing
awk_reldev='{printf "%.2e", sqrt(($1-$2)*($1-$2))/(0.5*sqrt(($1+$2)*($1+$2))) }'
awk_isgtr='{if($1>$2) {print "1"}}'

# ---------- common test parameters
IFO=H1
timestamp1=852443819
timestamp2=852445619
timestamp3=852447419
alpha=4.649850989853494
delta=-0.506281802989210

outCAP=antenna_pattern_test.dat
outPDS=detector_state_test.dat


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test0: internal consistency of different computeAM implementations in PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

pds_cmdline="${pds_code} -I $IFO -a $alpha -d $delta -t $timestamp1 >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi

pds_out_lalcomp=$(grep LALComputeAM ${outPDS} | tr ',' ' ')
a_lalcomp=$( echo $pds_out_lalcomp | awk '{print $3}' )
b_lalcomp=$( echo $pds_out_lalcomp | awk '{print $4}' )
echo "==> LALComputeAM:                    a=$a_lalcomp, b=$b_lalcomp"
pds_out_lalget=$(grep LALGetAMCoeffs ${outPDS} | tr ',' ' ')
a_lalget=$( echo $pds_out_lalget | awk '{print $3}' )
b_lalget=$( echo $pds_out_lalget | awk '{print $4}' )
echo "    LALGetAMCoeffs:                  a=$a_lalget, b=$b_lalget"
pds_out_lalnew=$(grep LALNewGetAMCoeffs ${outPDS} | tr ',' ' ')
a_lalnew=$( echo $pds_out_lalnew | awk '{print $3}' )
b_lalnew=$( echo $pds_out_lalnew | awk '{print $4}' )
echo "    LALNewGetAMCoeffs:               a=$a_lalnew, b=$b_lalnew"
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a_xlal=$( echo $pds_out_xlal | awk '{print $3}' )
b_xlal=$( echo $pds_out_xlal | awk '{print $4}' )
echo "    XLALComputeAntennaPatternCoeffs: a=$a_xlal, b=$b_xlal"

reldev_a_lalget=$(echo $a_lalcomp $a_lalget | awk "$awk_reldev")
reldev_b_lalget=$(echo $b_lalcomp $b_lalget | awk "$awk_reldev")
reldev_a_lalnew=$(echo $a_lalcomp $a_lalnew | awk "$awk_reldev")
reldev_b_lalnew=$(echo $b_lalcomp $b_lalnew | awk "$awk_reldev")
reldev_a_xlal=$(echo $a_lalcomp $a_xlal | awk "$awk_reldev")
reldev_b_xlal=$(echo $b_lalcomp $b_xlal | awk "$awk_reldev")

fail_a_get=$(echo $reldev_a_lalget $tolerance | awk "$awk_isgtr")
fail_b_get=$(echo $reldev_b_lalget $tolerance | awk "$awk_isgtr")
fail_a_new=$(echo $reldev_a_lalnew $tolerance | awk "$awk_isgtr")
fail_b_new=$(echo $reldev_b_lalnew $tolerance | awk "$awk_isgtr")
fail_a_xlal=$(echo $reldev_a_xlal $tolerance | awk "$awk_isgtr")
fail_b_xlal=$(echo $reldev_b_xlal $tolerance | awk "$awk_isgtr")

if [ "$fail_a_get" -o "$fail_b_get" -o "$fail_a_new" -o "$fail_b_new" -o "$fail_a_xlal" -o "$fail_b_xlal" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test1: single-sky-point, single-timestamp comparison to PrintDetectorState";
echo "----------------------------------------------------------------------------------------------------"

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --detector=$IFO --timeGPS=$timestamp1 --outputFile=$outCAP --Alpha=$alpha --Delta=$delta"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

eval $(sed -i '/^\%\%/d' $outCAP)
a_cap=$(awk '{print $4}' $outCAP)
b_cap=$(awk '{print $5}' $outCAP)
echo "==> lalapps_ComputeAntennaPattern:   a=$a_cap, b=$b_cap"

reldev_a_cap_xlal=$(echo $a_xlal $a_cap | awk "$awk_reldev")
reldev_b_cap_xlal=$(echo $b_xlal $b_cap | awk "$awk_reldev")
fail_a_cap=$(echo $reldev_a_cap_xlal $tolerance | awk "$awk_isgtr")
fail_b_cap=$(echo $reldev_b_cap_xlal $tolerance | awk "$awk_isgtr")

if [ "$fail_a_cap" -o "$fail_b_cap" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test2: sky grid";
echo "----------------------------------------------------------------------------------------------------"

skygrid=./skygrid_test.dat
alpha1=0.0
delta1=0.0
alpha2=0.0
delta2=0.5
alpha3=3.0
delta3=-0.5
echo "$alpha1 $delta1\n$alpha2 $delta2\n$alpha3 $delta3" >> $skygrid

## ----- run PrintDetectorState 3 times
rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha1 -d $delta1 -t $timestamp1 >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a1_pds=$( echo $pds_out_xlal | awk '{print $3}')
b1_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha2 -d $delta2 -t $timestamp1 >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a2_pds=$( echo $pds_out_xlal | awk '{print $3}')
b2_pds=$( echo $pds_out_xlal | awk '{print $4}')

rm $outPDS
pds_cmdline="${pds_code} -I $IFO -a $alpha3 -d $delta3 -t $timestamp1 >> $outPDS"
echo $pds_cmdline;
if ! eval $pds_cmdline; then
    echo "Error.. something failed when running '$pds_code' ..."
    exit 1
fi
pds_out_xlal=$(grep XLALComputeAntennaPatternCoeffs ${outPDS} | tr ',' ' ')
a3_pds=$( echo $pds_out_xlal | awk '{print $3}')
b3_pds=$( echo $pds_out_xlal | awk '{print $4}')

## ----- run ComputeAntennaPattern
cap_cmdline="${cap_code} --detector=$IFO --timeGPS=$timestamp1 --outputFile=$outCAP --skyGridFile=$skygrid"

echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi

eval $(sed -i '/^\%\%/d' $outCAP)
a1_cap=$( awk 'NR==1 {print $4}' $outCAP)
b1_cap=$( awk 'NR==1 {print $5}' $outCAP)
a2_cap=$( awk 'NR==2 {print $4}' $outCAP)
b2_cap=$( awk 'NR==2 {print $5}' $outCAP)
a3_cap=$( awk 'NR==3 {print $4}' $outCAP)
b3_cap=$( awk 'NR==3 {print $5}' $outCAP)
echo "==> alpha=0.0 delta= 0.0: lalapps_ComputeAntennaPattern: a=$a1_cap, b=$b1_cap / lalapps_PrintDetectorState: a=$a1_pds, b=$b1_pds"
echo "    alpha=0.0 delta= 0.5: lalapps_ComputeAntennaPattern: a=$a2_cap, b=$b2_cap / lalapps_PrintDetectorState: a=$a2_pds, b=$b2_pds"
echo "    alpha=3.0 delta=-0.5: lalapps_ComputeAntennaPattern: a=$a3_cap, b=$b3_cap / lalapps_PrintDetectorState: a=$a3_pds, b=$b3_pds"

reldev_a1=$(echo $a1_pds $a1_cap | awk "$awk_reldev")
reldev_b1=$(echo $b1_pds $b1_cap | awk "$awk_reldev")
fail_a1=$(echo $reldev_a1 $tolerance | awk "$awk_isgtr")
fail_b1=$(echo $reldev_b1 $tolerance | awk "$awk_isgtr")
reldev_a2=$(echo $a2_pds $a2_cap | awk "$awk_reldev")
reldev_b2=$(echo $b2_pds $b2_cap | awk "$awk_reldev")
fail_a2=$(echo $reldev_a2 $tolerance | awk "$awk_isgtr")
fail_b2=$(echo $reldev_b2 $tolerance | awk "$awk_isgtr")
reldev_a3=$(echo $a3_pds $a3_cap | awk "$awk_reldev")
reldev_b3=$(echo $b3_pds $b3_cap | awk "$awk_reldev")
fail_a3=$(echo $reldev_a3 $tolerance | awk "$awk_isgtr")
fail_b3=$(echo $reldev_b3 $tolerance | awk "$awk_isgtr")

if [ "$fail_a1" -o "$fail_b1" -o "$fail_a2" -o "$fail_b2" -o "$fail_a3" -o "$fail_b3" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
fi


echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test3: summing / averaging over timestamps from file";
echo "----------------------------------------------------------------------------------------------------"

## ----- run ComputeAntennaPattern with single-timestamp input, output
cap_cmdline="${cap_code} --detector=$IFO --timeGPS=$timestamp1,$timestamp2,$timestamp3 --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --mthopTimeStamps=0"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi
eval $(sed -i '/^\%\%/d' $outCAP)
a1_cap=$( awk 'NR==1 {print $4}' $outCAP)
b1_cap=$( awk 'NR==1 {print $5}' $outCAP)
a2_cap=$( awk 'NR==2 {print $4}' $outCAP)
b2_cap=$( awk 'NR==2 {print $5}' $outCAP)
a3_cap=$( awk 'NR==3 {print $4}' $outCAP)
b3_cap=$( awk 'NR==3 {print $5}' $outCAP)

## ----- externally compute sum, mean for test
asum=$( echo $a1_cap $a2_cap $a3_cap | awk '{print $1+$2+$3}' )
bsum=$( echo $b1_cap $b2_cap $b3_cap | awk '{print $1+$2+$3}' )
amean=$( echo $asum | awk '{print $1/3}' )
bmean=$( echo $bsum | awk '{print $1/3}' )

## ----- make timestampsfile
timestampsfile=./timestamps_test.dat
echo "$timestamp1 0\n$timestamp2 0\n$timestamp3 0" >> $timestampsfile

## ----- run ComputeAntennaPattern with timestampsfile input, summed output
cap_cmdline="${cap_code} --detector=$IFO --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --mthopTimeStamps=1"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi
eval $(sed -i '/^\%\%/d' $outCAP)
a_cap_sum=$(awk '{print $3}' $outCAP)
b_cap_sum=$(awk '{print $4}' $outCAP)
reldev_asum=$(echo $a_cap_sum $asum | awk "$awk_reldev")
reldev_bsum=$(echo $b_cap_sum $bsum | awk "$awk_reldev")
fail_asum=$(echo $reldev_asum $tolerance | awk "$awk_isgtr")
fail_bsum=$(echo $reldev_bsum $tolerance | awk "$awk_isgtr")

## ----- run ComputeAntennaPattern with timestampsfile input, mean output
cap_cmdline="${cap_code} --detector=$IFO --timeStampsFile=$timestampsfile --outputFile=$outCAP --Alpha=$alpha --Delta=$delta --mthopTimeStamps=2"
echo $cap_cmdline;
if ! eval $cap_cmdline; then
    echo "Error.. something failed when running '$cap_code' ..."
    exit 1
fi
eval $(sed -i '/^\%\%/d' $outCAP)
a_cap_mean=$(awk '{print $3}' $outCAP)
b_cap_mean=$(awk '{print $4}' $outCAP)
reldev_amean=$(echo $a_cap_mean $amean | awk "$awk_reldev")
reldev_bmean=$(echo $b_cap_mean $bmean | awk "$awk_reldev")
fail_amean=$(echo $reldev_amean $tolerance | awk "$awk_isgtr")
fail_bmean=$(echo $reldev_bmean $tolerance | awk "$awk_isgtr")

echo "==> --mthopTimeStamps=0, recomputed sum:  asum=$asum, bsum=$bsum"
echo "    --mthopTimeStamps=1,            sum:  asum=$a_cap_sum, bsum=$b_cap_sum"
echo "    --mthopTimeStamps=0, recomputed mean: <a> =$amean, <b> =$bmean"
echo "    --mthopTimeStamps=2,            mean: <a> =$a_cap_mean, <b> =$b_cap_mean"

if [ "$fail_asum" -o "$fail_bsum" -o "$fail_amean" -o "$fail_bmean" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
fi

echo "----------------------------------------------------------------------------------------------------"
echo "ComputeAntennaPattern Test4: comparing A,B,C,D with PredictFStat";
echo "----------------------------------------------------------------------------------------------------"

A_cap=$(awk '{print $5}' $outCAP)
B_cap=$(awk '{print $6}' $outCAP)
C_cap=$(awk '{print $7}' $outCAP)
D_cap=$(awk '{print $8}' $outCAP)

sftfile=./H1_test.sft
outPFS=./pfs_test.dat

mfd_cmdline="${mfd_code} --IFO=$IFO --outSingleSFT --outSFTbname=$sftfile --fmin=59.95 --Band=0.1 --timestampsFile=$timestampsfile"
echo $mfd_cmdline;
if ! eval $mfd_cmdline; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

pfs_cmdline="${pfs_code} --IFO=$IFO --h0=1e-24 --cosi=0 --psi=0 --phi0=0 --Freq=60 --Alpha=$alpha --Delta=$delta --DataFiles=$sftfile --outputFstat=$outPFS --SignalOnly --printFstat=0"
echo $pfs_cmdline;
if ! eval $pfs_cmdline; then
    echo "Error.. something failed when running '$pfs_code' ..."
    exit 1
fi

A_pfs=$(grep 'A =' ${outPFS} | tr -d 'A =;')
B_pfs=$(grep 'B =' ${outPFS} | tr -d 'B =;')
C_pfs=$(grep 'C =' ${outPFS} | tr -d 'C =;')
D_pfs=$(grep 'D =' ${outPFS} | tr -d 'D =;')

echo "==> lalapps_ComputeAntennaPattern:   A=$A_cap, B=$B_cap, C=$C_cap, D=$D_cap"
echo "    lalapps_PredictFStat:            A=$A_pfs,   B=$B_pfs,   C=$C_pfs,   D=$D_pfs"

reldev_A=$(echo $A_cap $A_pfs | awk "$awk_reldev")
fail_A=$(echo $reldev_A $tolerance | awk "$awk_isgtr")
reldev_B=$(echo $B_cap $B_pfs | awk "$awk_reldev")
fail_B=$(echo $reldev_B $tolerance | awk "$awk_isgtr")
reldev_C=$(echo $C_cap $C_pfs | awk "$awk_reldev")
fail_C=$(echo $reldev_C $tolerance | awk "$awk_isgtr")
reldev_D=$(echo $D_cap $D_pfs | awk "$awk_reldev")
fail_D=$(echo $reldev_D $tolerance | awk "$awk_isgtr")

if [ "$fail_A" -o "$fail_B" -o "$fail_C" -o "$fail_D" ]; then
    echo "==> FAILED at tolerance=$tolerance"
    retstatus=1
else
    echo "==> OK at tolerance=$tolerance"
    echo
    echo "========== OK. All ComputeAntennaPattern tests PASSED. =========="
    echo
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm $outCAP
    rm $outPDS
    rm $skygrid
    rm $timestampsfile
    rm $sftfile
    rm $outPFS
    echo "Cleaned up."
fi

exit $retstatus
