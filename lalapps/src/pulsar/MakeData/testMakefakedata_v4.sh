mfdCODE="lalapps_Makefakedata_v4"
cmpCODE="lalapps_compareSFTs"
extractCODE="lalapps_splitSFTs"

testDIR1="./mfdv4_TEST1"
testDIR2="./mfdv4_TEST2"
testDIR3="./mfdv4_TEST3"

# prepare test subdirectory
mkdir $testDIR1
mkdir $testDIR2
mkdir $testDIR3

tol="1e-4";	## tolerance on relative difference between SFTs in comparison
# input parameters
## FIXED
Tsft=1800
nTsft=20
timestamps="testT8_1800.txt"
refTime=701210229

## excercise non-integer cycle gaps in heterodyned timeseries
fmin=299.1001
Band=9.9998
fmax=$(echo $fmin $Band | awk '{printf "%.7g", $1 + $2}');
fUpper=311

## VARY
IFO=H1
aPlus=1.5
aCross=0.7
psi=0.5
phi0=0.9
Freq=300.2
alpha=1.7
delta=0.9

f1dot="-1.e-9"
f2dot="1e-14"

echo "------------------------------------------------------------"
echo " SIGNAL-ONLY - compare heterodyned SFTs with 'exact' ones"
echo "------------------------------------------------------------"

echo
echo "mfd_v4: producing SFTs via heterodyned timeseries (generationMode=0 [ALL_AT_ONCE] )..."
echo

mfdCL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$alpha --Delta=$delta --IFO=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=0"
cmdline="$mfdCODE $mfdCL --outSFTbname=$testDIR1";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi
# generate concatenated SFT
mfdCL="${mfdCL} --outSFTbname=${testDIR1}.sft --outSingleSFT"
cmdline="$mfdCODE $mfdCL";
echo "$cmdline (concatenated SFT version)";
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi

echo
echo "mfd_v4: producing SFTs via heterodyned timeseries (generationMode=1 [PER_SFT] )..."
echo

mfdCL="--Tsft=$Tsft --fmin=$fmin --Band=$Band --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$alpha --Delta=$delta --IFO=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=1"
cmdline="$mfdCODE $mfdCL --outSFTbname=$testDIR2";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi
# generate concatenated SFT
mfdCL="${mfdCL} --outSFTbname=${testDIR2}.sft --outSingleSFT"
cmdline="$mfdCODE $mfdCL";
echo "$cmdline (concatenated SFT version)";
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi

echo
echo "mfd_v4: producing SFTs via 'exact' timeseries (non-heterodyned)..."
echo

mfdCL="--Tsft=$Tsft --fmin=0 --Band=$fUpper --aPlus=$aPlus --aCross=$aCross --psi=$psi --phi0=$phi0 --Freq=${Freq} --Alpha=$alpha --Delta=$delta --IFO=$IFO --timestampsFile=$timestamps --refTime=$refTime --f1dot=$f1dot --f2dot=$f2dot --generationMode=0"
cmdline="$mfdCODE $mfdCL --outSFTbname=$testDIR3";
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$mfdCODE' ..."
    exit 1
fi

echo "... and extracting the relevant frequency band ..."
echo

## extract relevant frequency-band
fmax2=`echo "${fmax} + 0.0006" | bc`   ## add 0.0006 to bandwidth to get same number of SFT bins as MDFv4
Band2=`echo "${Band} + 0.0006" | bc`   ## add 0.0006 to bandwidth to get same number of SFT bins as MDFv4
extractCL="-fs $fmin -fe $fmax2 -fb $Band2 -n ${testDIR3} -- $testDIR3/*.sft"
cmdline="$extractCODE $extractCL"
echo $cmdline;
if ! eval $cmdline; then
    echo "Error.. something failed when running '$extractCODE' ..."
    exit 1
fi

echo
echo "comparison of resulting SFTs:"

cmdline="$cmpCODE -V -e $tol -1 '${testDIR1}/*.sft' -2 '${testDIR3}/*SFT_NB*.sft'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... SFTs differ by more than $tol. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
cmdline="$cmpCODE -V -e $tol -1 '${testDIR2}/*.sft' -2 '${testDIR3}/*SFT_NB*.sft'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... SFTs differ by more than $tol. Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
echo "comparison of concatenating SFTs:"

cmdline="$cmpCODE -V -e 1e-10 -1 '${testDIR1}/*.sft' -2 '${testDIR1}.sft'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... concatenated SFTs differ! Something might be wrong..."
    exit 2
else
    echo "OK."
fi

echo
cmdline="$cmpCODE -V -e 1e-10 -1 '${testDIR2}/*.sft' -2 '${testDIR2}.sft'"
echo ${cmdline}
if ! eval $cmdline; then
    echo "OUCH... concatenated SFTs differ! Something might be wrong..."
    exit 2
else
    echo "OK."
fi
