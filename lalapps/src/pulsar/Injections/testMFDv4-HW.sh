#!/bin/bash

##
## special test-script to test hardware-injection feature in mfd_v4.
## This test is similar to testMakefakedata.sh, except that we compare
## the binary timeseries output by mfd_v4 versus that of mfd_v2.
##
## We use S5 pulsar0 parameters for the test.
##

## run all LALApps programs with memory debugging
export LAL_DEBUG_LEVEL="${LAL_DEBUG_LEVEL},memdbg"

## allow 'make test' to work from builddir != srcdir
builddir="./";
testDIR="./mfd_TEST-HW"

v2_code="${builddir}lalapps_makefakedata_v2"
v4_code="${builddir}lalapps_Makefakedata_v4"
comp_code="${builddir}lalapps_compareTS"

v2_strain="${testDIR}/v2-strain.dat"
v4_strain="${testDIR}/v4-strain.dat"

tol=1e-4	## error tolerance for v2-v4 strain comparison

## ----- user-controlled level of debug-output detail
if [ -n "$DEBUG" ]; then
    debug=${DEBUG}
else
    debug=0	## default=quiet
fi

if [ -z "${LAL_DATA_PATH}" ]; then
    echo
    echo "Need environment-variable LAL_DATA_PATH to be set to include"
    echo "your ephemeris-directory (e.g. /usr/local/share/lalpulsar)"
    echo "This might indicate an incomplete LAL+LALPULSAR installation"
    echo
    exit 1
fi

if [ -n "${LALPULSAR_DATADIR}" ]; then
    v2_code="${v2_code} -E ${LALPULSAR_DATADIR}"
fi

#prepare test subdirectory
if [ ! -d "$testDIR" ]; then
    mkdir $testDIR
else
## cleanup: remove previous timeseries output
    rm -f $testDIR/* || true
fi

# input parameters
Tsft=20
duration=100       ## 100 seconds
nTsft=5  	   ## duration / Tsft
IFO=LHO
fmin=0
Band=8192
startTime=705283213

## ----- S5 pulsar0 parameters: http://blue.ligo-wa.caltech.edu/scirun/S5/HardwareInjection/Details/pulsar/test2_H1/
refTime=751680013       ## pulsar reference time in SSB frame
aPlus=1.6100e-21        ## plus-polarization signal amplitude
aCross=1.5685e-21       ## cross-polarization signal amplitude
psi=0.770087086         ## polarization angle psi
phi0=2.66               ## phase at tRef
Freq=265.5771052   	## GW frequency at tRef
Delta=-0.981180225      ## latitude (delta,declination) in radians
Alpha=1.248816734       ## longitude (alpha, right ascension) in radians

f1dot=-4.15E-12         ## spindown parameters d/dt f0
f2dot=0.0
f3dot=0.0
## --------------------------------------------------pp

## now we generate an equivalent mfd_v2 input file:
v2_cfg=In.data-v2
v4_cfg=In.data-v4
v4_log=v4.log

v2_CL="-i $v2_cfg -I $IFO -S $refTime -G $startTime -b"
v4_CL="-I $IFO @${v4_cfg} --Tsft=$Tsft --startTime=$startTime --duration=$duration -l $v4_log --generationMode=1 -b"

## produce In.data-v2 file for makefakedata_v2
echo "
$Tsft   %Tsft_in_sec
$nTsft  %nTsft
$fmin   %first_SFT_frequency_in_Hz
$Band   %SFT_freq_band_in_Hz
0.0     %sigma
$aPlus  %Aplus
$aCross %Across
$psi    %psi
$phi0   %phi0
$Freq   %f0
$Delta  %latitude_in_radians
$Alpha  %longitude_in_radians
3       %max_spin-down_param_order
$f1dot  %value of first spindown
$f2dot  %value of second spindown
$f3dot  %value of third spindown
" > $v2_cfg

## produce In.data-v4 file for makefakedata_v4
echo "
## ----- S5 pulsar0 parameters: http://blue.ligo-wa.caltech.edu/scirun/S5/HardwareInjection/Details/pulsar/test2_H1/
refTime         = $refTime              ## pulsar reference time in SSB frame

aPlus           = $aPlus                ## plus-polarization signal amplitude
aCross          = $aCross               ## cross-polarization signal amplitude
psi             = $psi                  ## polarization angle psi
phi0            = $phi0                 ## phase at tRef
Freq            = $Freq                 ## GW frequency at tRef
Delta           = $Delta                ## latitude (delta,declination) in radians
Alpha           = $Alpha                ## longitude (alpha, right ascension) in radians

f1dot           = $f1dot                ## spindown parameters d/dt f0
f2dot           = $f2dot
f3dot           = $f3dot
## --------------------------------------------------
" > $v4_cfg

echo "--------------------------------------------------------------------------------"
echo "Testing binary strain output of pulsar-signal for hardware injections ... "
echo "--------------------------------------------------------------------------------"
echo
echo "1) Running mfd_v2 reference-code:"
cmdline="$v2_code $v2_CL > ${v2_strain}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$v2_code' ..."
    exit 1
fi

echo "2) Running mfd_v4 injection code:"
cmdline="$v4_code $v4_CL > ${v4_strain}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$v4_code' ..."
    exit 1
fi


echo
echo "3) Comparison of resulting binary strains ..."
echo
cmdline="$comp_code -1 ${v2_strain} -2 ${v4_strain} --relErrorMax=$tol"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. comparison failed using '$comp_code' ..."
    exit 1
else
    echo "OK. Test passed within tolerance of $tol."
    echo
fi


## clean up files [allow turning off via 'NOCLEANUP' environment variable
if [ -z "$NOCLEANUP" ]; then
    rm -rf $testDIR $v2_cfg $v4_cfg $v4_log
fi


exit 0
