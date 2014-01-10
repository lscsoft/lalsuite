#!/bin/sh

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## test of lalapps_pulsar_crosscorr_v2; currently just makes sure it runs without errors
extra_args="$@"

builddir="./";
injectdir="../Injections/"

##---------- names of codes and input/output files
mfd_code="lalapps_Makefakedata_v4"
pcc_code="lalapps_pulsar_crosscorr_v2"

mfd_path="${injectdir}${mfd_code}"
pcc_path="${builddir}${pcc_code}"

SFTdir="./testPulsarCrossCorr_v2_sfts"

if [ -n "${LALPULSAR_DATADIR}" ]; then
    mfd_code="${mfd_code} -E ${LALPULSAR_DATADIR}"
fi

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

Ftolerance=0.05
# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=827884814
duration=7200
endTime=827892014
refTime=827884814

# Sky coordinates of Sco X-1
alphaRad=4.2756992385
deltaRad=-0.272973858335

mfd_fmin=149.8
mfd_Band=0.4
mfd_h0=3e-22
mfd_cosi=0
mfd_psi=0
mfd_phi0=0
mfd_Freq=150.0
mfd_noiseSqrtSh=3e-23
mfd_seed1=201401090
mfd_seed2=201401091
mfd_ifo1=H1
mfd_ifo2=L1

pcc_fStart=149.9995
pcc_fBand=.001
pcc_maxLag=1800
pcc_orbitAsiniSec=1.40
pcc_orbitAsiniSecBand=0.10
pcc_orbitPSec=68023.7136
pcc_orbitTimeAsc=1245967374
pcc_orbitTimeAscBand=20
pcc_numBins=1

if [ ! -d "$SFTdir" ]; then
    mkdir $SFTdir
else
    rm -f ${SFTdir}/*
fi

mfd_CL="--fmin=$mfd_fmin --Band=$mfd_Band --Freq=$mfd_Freq --outSFTbname=$SFTdir --noiseSqrtSh=$mfd_noiseSqrtSh --Alpha=$alphaRad --Delta=$deltaRad --Tsft=$Tsft --startTime=$startTime --duration=$duration --h0=$mfd_h0 --cosi=$mfd_cosi --psi=$mfd_psi --phi0=$mfd_phi0"
mfd_CL1="${mfd_CL} --IFO=$mfd_ifo1 --randSeed=$mfd_seed1"
mfd_CL2="${mfd_CL} --IFO=$mfd_ifo2 --randSeed=$mfd_seed2"

pcc_CL="--startTime=$startTime --endTime=$endTime --sftLocation='$SFTdir/*.sft' --fStart=$pcc_fStart --fBand=$pcc_fBand --alphaRad=$alphaRad --deltaRad=$deltaRad --maxLag=$pcc_maxLag --orbitAsiniSec=$pcc_orbitAsiniSec --orbitAsiniSecBand=$pcc_orbitAsiniSecBand --orbitPSec=$pcc_orbitPSec --orbitTimeAsc=$pcc_orbitTimeAsc --orbitTimeAscBand=$pcc_orbitTimeAscBand --numBins=$pcc_numBins"

## ---------- Run MFDv4 ----------
cmdline="$mfd_path $mfd_CL1";
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${mfd_code} ... "
if ! eval "$cmdline 2> /dev/null"; then
    echo "FAILED:"
    echo $cmdline
    exit 1
else
    echo "OK."
fi

cmdline="$mfd_path $mfd_CL2";
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${mfd_code} ... "
if ! eval "$cmdline 2> /dev/null"; then
    echo "FAILED:"
    echo $cmdline
    exit 1
else
    echo "OK."
fi

## ---------- Run PulsarCrossCorr_v2 ----------
cmdline="$pcc_path $pcc_CL"
if [ "$DEBUG" ]; then echo $cmdline; fi
echo -n "Running ${pcc_code} ... "
if ! tmp=`eval $cmdline 2> /dev/null`; then
    echo "FAILED:"
    echo $cmdline
    exit 1;
else
    echo "OK."
fi

## clean up files
if [ -z "$NOCLEANUP" ]; then
    rm -rf $SFTdir
fi

exit $res;
