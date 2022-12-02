if test "${LALFRAME_ENABLED}" = false; then
    echo "Skipping test: requires LALFrame"
    exit 77
fi

mfdv4_CODE="lalpulsar_Makefakedata_v4"
mfdv5_CODE="lalpulsar_Makefakedata_v5"
cmp_CODE="lalpulsar_compareSFTs"
dump_code="lalpulsar_dumpSFT"
split_code="lalpulsar_splitSFTs"

testDIR="."
testDIRnb="./nb"

# prepare test subdirectories
mkdir -p $testDIRnb

tol=1e-3;	## tolerance on relative difference between SFTs in comparison
# input parameters
## ---------- data parameters ----------
Tsft=1800
nTsft=20
timestamps=./testT8_1800.txt

## excercise non-integer cycle gaps in heterodyned timeseries
fmin=299
Band=10
fmax=$(echo $fmin $Band | LC_ALL=C awk '{printf "%.7g", $1 + $2}');

IFO1=H1
IFO2=L1
sqrtSn1=0.5;
sqrtSn2=1.2;	## for comparison with 2 calls to mfdv4 and fixed see, the 2nd IFO noise must be 0

timestamps1="${testDIR}/H1-timestamps.dat"
timestamps2="${testDIR}/L1-timestamps.dat"

cp ${timestamps} ${timestamps1}
echo "701210229 0" > ${timestamps2}
echo "701230229 0" >> ${timestamps2}
echo "701240229 0" >> ${timestamps2}

## ---------- signal parameters ----------
## ----- signal 1
s1_h0=0.73
s1_cosi=0.1
s1_psi=0.5
s1_phi0=0.9
s1_refTime=701210229
s1_Alpha=1.7
s1_Delta=0.9
s1_Freq=299.12
s1_f1dot=0
s1_f2dot=0
## ----- signal 2
s2_h0=2.5
s2_cosi=-0.5
s2_aPlus=`echo "scale = 6; 0.5 * ${s2_h0} * (1.0 + ${s2_cosi} * ${s2_cosi})" | bc`
s2_aCross=`echo "scale = 6; ${s2_h0} * ${s2_cosi}" | bc`
s2_psi=1.2
s2_phi0=1.5
s2_refTime=711210229
s2_Alpha=3.7
s2_Delta=-0.5
s2_Freq=300.12
s2_f1dot=-1e-10
s2_f2dot=0
## ----- signal 3
s3_h0=3.1
s3_cosi=0.5
s3_psi=-1.2
s3_phi0=2.5
s3_refTime=721210229
s3_Alpha=0.5
s3_Delta=1.2
s3_Freq=300.00
s3_f1dot=-1e-9
s3_f2dot=-2e-19
s3_transientWindowType=rect
s3_transientStartTime=701230229
s3_transientTauDays=0.5
s3_orbitasini=1.44
s3_orbitPeriod=68400
s3_orbitTp=57058.434479MJD
s3_orbitArgp=1
s3_orbitEcc=0.1
# --------------------

injString="{Alpha=${s1_Alpha};Delta=${s1_Delta};refTime=${s1_refTime};Freq=${s1_Freq};f1dot=${s1_f1dot};f2dot=${s1_f2dot};h0=${s1_h0};cosi=${s1_cosi};psi=${s1_psi};phi0=${s1_phi0};}"

## ---------- signal file 1 ----------
injFile1=${testDIR}/injectionS1.dat
echo "[Pulsar 1]" >> ${injFile1}
echo "Alpha = ${s1_Alpha}" >> ${injFile1}
echo "Delta = ${s1_Delta}" >> ${injFile1}
echo "refTime = ${s1_refTime}" >> ${injFile1}
echo "Freq = ${s1_Freq}" >> ${injFile1}
echo "f1dot = ${s1_f1dot}" >> ${injFile1}
echo "f2dot = ${s1_f2dot}" >> ${injFile1}
echo "h0 = ${s1_h0}" >> ${injFile1}
echo "cosi = ${s1_cosi}" >> ${injFile1}
echo "psi = ${s1_psi}" >> ${injFile1}
echo "phi0 = ${s1_phi0}" >> ${injFile1}
echo >> ${injFile1}
## ---------- signal file 2 ----------
injFile2=${testDIR}/injectionS2.dat
echo "Alpha = ${s2_Alpha}" >> ${injFile2}
echo "Delta = ${s2_Delta}" >> ${injFile2}
echo "refTime = ${s2_refTime}" >> ${injFile2}
echo "Freq = ${s2_Freq}" >> ${injFile2}
echo "f1dot = ${s2_f1dot}" >> ${injFile2}
echo "f2dot = ${s2_f2dot}" >> ${injFile2}
echo "aPlus = ${s2_aPlus}" >> ${injFile2}
echo "aCross = ${s2_aCross}" >> ${injFile2}
echo "psi = ${s2_psi}" >> ${injFile2}
echo "phi0 = ${s2_phi0}" >> ${injFile2}
echo >> ${injFile2}
## ---------- add section for Pulsar 3 into signal-file 2 ----------
echo "[Pulsar 3]" >> ${injFile2}
echo "Alpha = ${s3_Alpha}" >> ${injFile2}
echo "Delta = ${s3_Delta}" >> ${injFile2}
echo "refTime = ${s3_refTime}" >> ${injFile2}
echo "Freq = ${s3_Freq}" >> ${injFile2}
echo "f1dot = ${s3_f1dot}" >> ${injFile2}
echo "f2dot = ${s3_f2dot}" >> ${injFile2}
echo "h0 = ${s3_h0}" >> ${injFile2}
echo "cosi = ${s3_cosi}" >> ${injFile2}
echo "psi = ${s3_psi}" >> ${injFile2}
echo "phi0 = ${s3_phi0}" >> ${injFile2}
echo "transientWindowType = ${s3_transientWindowType}"  >> ${injFile2}
echo "transientStartTime = ${s3_transientStartTime}"  >> ${injFile2}
echo "transientTauDays = ${s3_transientTauDays}"  >> ${injFile2}
echo "orbitasini = ${s3_orbitasini}" >> ${injFile2}
echo "orbitPeriod = ${s3_orbitPeriod}" >> ${injFile2}
echo "orbitTp = ${s3_orbitTp}" >> ${injFile2}
echo "orbitArgp = ${s3_orbitArgp}" >> ${injFile2}
echo "orbitEcc = ${s3_orbitEcc}" >> ${injFile2}

echo >> ${injFile2}

## ---------- output parameters ----------
sftsv4_1=${testDIR}/${IFO1}-sftsv4.sft
sftsv4_2=${testDIR}/${IFO2}-sftsv4.sft
sftsv5_1_meth1=${testDIR}/H-*_mfdv5-*.sft
sftsv5_2_meth1=${testDIR}/L-*_mfdv5-*.sft

sftsv5_1_meth2=${testDIR}/H-*_mfdv5meth2-*.sft
sftsv5_2_meth2=${testDIR}/L-*_mfdv5meth2-*.sft
sftsv5_meth2=${testDIR}/*_mfdv5meth2-*.sft

## ----------
## produce SFTs for 2 detectors, containing Gaussian noise + N signals, compare between mfdv4 and mfdv5
## ----------
window="tukey"
windowParam="0.5"

echo
echo "========== MFDv4 =========="
echo
mfdv4_CL="$mfdv4_CODE ${mfdv4_extra} --fmin=$fmin --Band=$Band --generationMode=0 --outSingleSFT"
echo "---------- mfdv4: inject first signal ----------"
sig1="--refTime=${s1_refTime} --h0=${s1_h0} --cosi=${s1_cosi} --psi=${s1_psi} --phi0=${s1_phi0} --Freq=${s1_Freq} --Alpha=${s1_Alpha} --Delta=${s1_Delta} --f1dot=${s1_f1dot} --f2dot=${s1_f2dot}"
##----- first IFO
out_IFO1="--IFO=${IFO1} --timestampsFile=${timestamps1} --window=${window} --windowParam=${windowParam} --outSFTbname=${sftsv4_1} --noiseSqrtSh=${sqrtSn1} --randSeed=1"
cmdline="$mfdv4_CL ${sig1} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --timestampsFile=${timestamps2} --window=${window} --windowParam=${windowParam} --outSFTbname=${sftsv4_2} --noiseSqrtSh=${sqrtSn2} --randSeed=2"
cmdline="$mfdv4_CL  ${sig1} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "---------- mfdv4: inject second signal on top ----------"
sig2="--refTime=${s2_refTime} --h0=${s2_h0} --cosi=${s2_cosi} --psi=${s2_psi} --phi0=${s2_phi0} --Freq=${s2_Freq} --Alpha=${s2_Alpha} --Delta=${s2_Delta} --f1dot=${s2_f1dot} --f2dot=${s2_f2dot}"
##----- first IFO
out_IFO1="--IFO=${IFO1} --noiseSFTs=${sftsv4_1} --outSFTbname=${sftsv4_1}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --noiseSFTs=${sftsv4_2} --outSFTbname=${sftsv4_2}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
echo "---------- mfdv4: inject third signal on top ----------"
sig2="--refTime=${s3_refTime} --h0=${s3_h0} --cosi=${s3_cosi} --psi=${s3_psi} --phi0=${s3_phi0} --Freq=${s3_Freq} --Alpha=${s3_Alpha} --Delta=${s3_Delta} --f1dot=${s3_f1dot} --f2dot=${s3_f2dot} --transientWindowType=${s3_transientWindowType} --transientStartTime=${s3_transientStartTime} --transientTauDays=${s3_transientTauDays} --orbitasini=${s3_orbitasini} --orbitPeriod=${s3_orbitPeriod} --orbitTp=${s3_orbitTp} --orbitArgp=${s3_orbitArgp} --orbitEcc=${s3_orbitEcc}"

##----- first IFO
out_IFO1="--IFO=${IFO1} --noiseSFTs=${sftsv4_1} --outSFTbname=${sftsv4_1}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO1}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi
##----- second IFO
out_IFO2="--IFO=${IFO2} --noiseSFTs=${sftsv4_2} --outSFTbname=${sftsv4_2}"
cmdline="$mfdv4_CL ${sig2} ${out_IFO2}"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv4_CODE' ..."
    exit 1
fi

echo
echo "========== MFDv5 =========="
echo
mfdv5_CL="$mfdv5_CODE ${mfdv5_extra} --outSingleSFT --outSFTdir=${testDIR}"

echo "----- Method 1: single multi-IFO, multi-signal call"
outIFOs="--IFOs=${IFO1},${IFO2} --timestampsFiles=${timestamps1},${timestamps2} --sqrtSX=${sqrtSn1},${sqrtSn2} --randSeed=1"
sig13="--injectionSources='${injFile1},${injFile2}'"
cmdline="$mfdv5_CL ${outIFOs} ${sig13} --fmin=$fmin --Band=$Band --SFTWindowType=${window} --SFTWindowParam=${windowParam}"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi

echo
echo "----- Method 2: and again the same, using different input methods"
outIFOs="--IFOs=${IFO1},${IFO2} --timestampsFiles=${timestamps1},${timestamps2} --sqrtSX=${sqrtSn1},${sqrtSn2} --randSeed=1"
sig1="--injectionSources='${injString}'"
sig23="--injectionSources='${injFile2}'"
cmdline1="$mfdv5_CL ${outIFOs} ${sig1} --SFTWindowType=${window} --SFTWindowParam=${windowParam} --outLabel='mfdv5meth2' --fmin=$fmin --Band=$Band"
echo $cmdline1;
if ! eval $cmdline1; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi
cmdline2="$mfdv5_CL ${sig23} --noiseSFTs='${sftsv5_meth2}' --outLabel='mfdv5meth2'"
echo $cmdline2;
if ! eval $cmdline2; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi

echo
echo "--------------------------------------------------"
echo "Comparison of resulting (concatenated) SFTs:"
echo "--------------------------------------------------"

echo "---------- compare mfdv5 Method 1 SFTs ----------"
cmdline="$cmp_CODE -V -e ${tol} -1 ${sftsv4_1} -2 '${sftsv5_1_meth1}'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

cmdline="$cmp_CODE -V -e ${tol} -1 ${sftsv4_2} -2 '${sftsv5_2_meth1}'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

echo
echo "---------- compare mfdv5 Method 2  SFTs ----------"
cmdline="$cmp_CODE -V -e ${tol} -1 ${sftsv4_1} -2 '${sftsv5_1_meth2}'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

cmdline="$cmp_CODE -V -e ${tol} -1 ${sftsv4_2} -2 '${sftsv5_2_meth2}'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "Failed. SFTs produced by makefakedata_v4 and makefakedata_v5 differ by more than ${tol}!"
    exit 2
else
    echo "OK."
fi

echo
echo "--------------------------------------------------"
echo "Test frame in/output and sub-band selection"
echo "--------------------------------------------------"

startTime=1238112018

mfdv5_cmd_Base="$mfdv5_CODE --outSFTdir=$testDIR --outSingleSFT=False"

cmdline="$mfdv5_cmd_Base --IFOs=H1 --fmin=0.0 --Band=512.0 --sqrtSX=$sqrtSn1 --startTime=$startTime --duration=$Tsft --outFrameDir=$testDIR"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi

cache="${testDIR}/H1_frame_cache.lcf"
pwd=$("pwd")
echo "H H1_mfdv5 $startTime $Tsft file://localhost/${pwd}/${testDIR}/H-H1_mfdv5-${startTime}-${Tsft}.gwf" > $cache

cmdline="$mfdv5_cmd_Base --inFrames=$cache --inFrChannels=H1:mfdv5 --fmin=$fmin --Band=$Band --outLabel=mfdv5fromframes"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdv5_CODE' ..."
    exit 1
fi

SFTs_orig="H-1_H1_${Tsft}SFT_mfdv5-${startTime}-${Tsft}.sft"
SFTs_split="H-1_H1_${Tsft}SFT_mfdv5_NBF0${fmin}Hz0W00${Band}Hz0-${startTime}-${Tsft}.sft"
SFTs_fromframes="H-1_H1_${Tsft}SFT_mfdv5fromframes-${startTime}-${Tsft}.sft"
dumps_split="${testDIR}/dump_${SFTs_split}.txt"
dumps_fromframes="${testDIR}/dump_${SFTs_fromframes}.txt"
SFTs_orig="${testDIR}/${SFTs_orig}"
SFTs_split="${testDIRnb}/${SFTs_split}"
SFTs_fromframes="${testDIR}/${SFTs_fromframes}"

fmax=$(echo "$fmin + $Band" | bc)
cmdline="$split_code -fs $fmin -fb $Band -fe $fmax -n $testDIRnb -- $SFTs_orig"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$split_code' ..."
    exit 1
fi

cmdline="$dump_code --dataOnly --SFTfiles=$SFTs_split > $dumps_split"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dump_code' ..."
    exit 1
fi

cmdline="$dump_code --dataOnly --SFTfiles=$SFTs_fromframes > $dumps_fromframes"
echo $cmdline
if ! eval $cmdline; then
    echo "Error.. something failed when running '$dump_code' ..."
    exit 1
fi

diff_cmd="diff -s $dumps_split $dumps_fromframes"
echo $diff_cmd
if ! eval $diff_cmd; then
    echo "error: $dumps_split and $dumps_fromframes should be equal."
    exit 1
fi

echo
echo "-------------------------------------------------------"
echo " Test proper failure whenever windows are misspecified "
echo "-------------------------------------------------------"
mkdir -p testwinpar/
base_call="$mfdv5_CODE --fmin 10 --Band 1 --IFOs H1 --startTime 100000000 --duration=1800 --outSFTdir=testwinpar/"
base_sft="./testwinpar/H-1_H1_1800SFT_mfdv5-100000000-1800.sft" # Produced by `base_call`
echo "Creating ${base_sft}..."
echo ${base_call}
if ! eval ${base_call}; then
    echo "Should not have failed!"
    exit 1
fi

no_param_call="$mfdv5_CODE --noiseSFTs=./testwinpar/H-1_H1_1800SFT_mfdv5-100000000-1800.sft --SFTWindowType=tukey"
echo "Calling with --noiseSFTs and --SFTWindowType=tukey without specifying window param, should fail:"
echo ${no_param_call}
if ! eval ${no_param_call}; then
    echo "Failed, as expected."
else
    echo "Did not fail, but it should have!"
    exit 1
fi

echo
echo "-------------------------------------------------------"
echo " Test creation of SFTs with public filenames "
echo "-------------------------------------------------------"
mkdir -p pubSFT/
cmdline="$mfdv5_CODE --fmin 10 --Band 1 --IFOs H1 --startTime 200000000 --duration=1800 -O 100 -R 1 --outLabel=testPubSFT --outSFTdir=pubSFT/"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$mfdv5_CODE' ..."
    exit 1
fi
for pubSFTname in "./pubSFT/H-1_H1_1800SFT_O100SIM+R1+CtestPubSFT+WRECT-200000000-1800.sft"; do
    if ! test -f "${pubSFTname}"; then
        echo "Error: '$mfdv5_CODE' failed to create public SFT '${pubSFTname}' as expected"
        exit 1
    fi
done
cmdline="$mfdv5_CODE --fmin 10 --Band 1 --IFOs 'H1,L1' --startTime 300000000 --duration=1800 --SFTWindowType=hann -O 100 -R 1 --outFrChannels='H1:test123,L1:test456' --outSFTdir=pubSFT/"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error: something failed when running '$mfdv5_CODE' ..."
    exit 1
fi
for pubSFTname in "./pubSFT/H-1_H1_1800SFT_O100SIM+R1+Ctest123+WHANN-300000000-1800.sft" "./pubSFT/L-1_L1_1800SFT_O100SIM+R1+Ctest456+WHANN-300000000-1800.sft"; do
    if ! test -f "${pubSFTname}"; then
        echo "Error: '$mfdv5_CODE' failed to create public SFT '${pubSFTname}' as expected"
        exit 1
    fi
done
