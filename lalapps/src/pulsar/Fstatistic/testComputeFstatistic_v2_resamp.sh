##---------- names of codes and input/output files
mfd_code="lalapps_Makefakedata_v5"
cmp_code="lalapps_compareFstats"
cfs_code="lalapps_ComputeFstatistic_v2"

# ---------- fixed parameter of our test-signal
Tsft=1800;
startTime=711595934
duration=144000		## 40 hours
refTime=611595934	## ~3.5 years prior to startTime

mfd_FreqBand=2.0;

Alpha=2.0
Delta=-0.5
AlphaBand=0.1
DeltaBand=0.1
dAlpha=0.05
dDelta=0.05

h0=1
cosi=-0.3
psi=0.6
phi0=1.5

Freq=100.12345
f1dot=-1e-10;

## mfd-specific bands
mfd_fmin=$(echo $Freq $mfd_FreqBand | awk '{printf "%g", $1 - $2 / 2.0}');

## cfs search bands
Nfreq=500
cfs_dFreq=$(echo $duration | awk '{printf "%.16g", 1.0 / ( 2.0 * $1 ) }');		## 1/(4T) frequency resolution
cfs_FreqBand=$(echo $Nfreq $cfs_dFreq | awk '{printf "%.16g", $1 * $2 }');		## band corresponding to fixed number of frequency bins
cfs_Freq=$(echo $Freq $Nfreq $cfs_dFreq | awk '{printf "%.16g", $1 - $2/2 * $3}');	## center search bin on exact signal frequency (for parameter estimation)

Nf1dot=10
cfs_df1dot=$(echo $duration | awk '{printf "%g", 1.0 / ($1 * $1) }');			## 1/(T^2) resolution in f1dot
cfs_f1dotBand=$(echo $Nf1dot $cfs_df1dot | awk '{printf "%.16g", $1 * $2 }');		## band corresponding to fixed number of f1dot bins
cfs_f1dot=$(echo $f1dot $Nf1dot $cfs_cf1dot | awk '{printf "%.16g", $1 - $2/2 * $3 }');	## center f1dot band bin on exact signal f1dot (for parameter estimation)

cfs_nCands=100000	## toplist length: keep N cands

sqrtSX="5,4"
haveNoise=true;

## ----- define output directory and files
SFTdir=.
SFTdir2=./injSourcesComp
mkdir -p $SFTdir2

outfile_Comp=./Fstat_Comp.dat
loudest_Comp=./Loudest_Comp.dat

outfile_Demod=./Fstat_Demod.dat
loudest_Demod=./Loudest_Demod.dat

outfile_Resamp=./Fstat_Resamp.dat
loudest_Resamp=./Loudest_Resamp.dat

outfile_Resamp_otfn=./Fstat_Resamp_otfn.dat
loudest_Resamp_otfn=./Loudest_Resamp_otfn.dat

timefile_Comp=./timing_Comp.dat
timefile_Demod=./timing_Demod.dat
timefile_Resamp=./timing_Resamp.dat
timefile_Resamp_otfn=./timing_Resamp_otfn.dat

##--------------------------------------------------
## test starts here
##--------------------------------------------------

echo
echo "----------------------------------------------------------------------"
echo " STEP 1: Generate Fake Data"
echo "----------------------------------------------------------------------"
echo
injectionSources="{refTime=${refTime}; Freq=$Freq; f1dot=$f1dot; Alpha=$Alpha; Delta=$Delta; h0=$h0; cosi=$cosi; psi=$psi; phi0=$phi0;}"
dataSpec="--Tsft=$Tsft --startTime=$startTime --duration=$duration --sqrtSX=${sqrtSX} --fmin=$mfd_fmin --Band=$mfd_FreqBand"
mfd_CL_comp="${dataSpec} --injectionSources='${injectionSources}' --outSingleSFT --outSFTdir=${SFTdir2} --randSeed=1 --IFOs=H1,L1"
cmdline="$mfd_code $mfd_CL_comp "
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

echo
mfd_CL="${dataSpec} --outSingleSFT --outSFTdir=${SFTdir} --randSeed=1 --IFOs=H1,L1"
cmdline="$mfd_code $mfd_CL "
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$mfd_code' ..."
    exit 1
fi

cfs_CL=" --refTime=${refTime} --Dterms=8 --Alpha=$Alpha --Delta=$Delta  --AlphaBand=$AlphaBand --DeltaBand=$DeltaBand --dAlpha=$dAlpha --dDelta=$dDelta --Freq=$cfs_Freq --FreqBand=$cfs_FreqBand --dFreq=$cfs_dFreq --f1dot=$cfs_f1dot --f1dotBand=$cfs_f1dotBand --df1dot=${cfs_df1dot} --NumCandidatesToKeep=${cfs_nCands}"
if [ "$haveNoise" != "true" ]; then
    cfs_CL="$cfs_CL --SignalOnly"
fi
echo
echo "----------------------------------------------------------------------"
echo "STEP 2: run directed CFS_v2 with LALDemod with injectionSources in MFDv5"
echo "----------------------------------------------------------------------"
echo
cmdline="$cfs_code $cfs_CL --FstatMethod=DemodBest --DataFiles='$SFTdir2/*.sft' --outputFstat=$outfile_Comp --outputTiming=$timefile_Comp --outputLoudest=${loudest_Comp} "
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo "STEP 3: run directed CFS_v2 with LALDemod method"
echo "----------------------------------------------------------------------"
echo
cmdline="$cfs_code $cfs_CL --injectionSources='${injectionSources}' --FstatMethod=DemodBest --DataFiles='$SFTdir/*.sft' --outputFstat=$outfile_Demod --outputTiming=$timefile_Demod --outputLoudest=${loudest_Demod} "
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 4: run directed CFS_v2 with resampling"
echo "----------------------------------------------------------------------"
echo
cmdline="$cfs_code $cfs_CL --injectionSources='${injectionSources}' --FstatMethod=ResampBest --DataFiles='$SFTdir/*.sft' --outputFstat=$outfile_Resamp --outputTiming=$timefile_Resamp  --outputLoudest=${loudest_Resamp}"
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 5: run directed CFS_v2 with resampling and on-the-fly noise"
echo "----------------------------------------------------------------------"
echo
cmdline="$cfs_code $cfs_CL --injectionSources='${injectionSources}' --injectSqrtSX=${sqrtSX} --IFOs='H1,L1' --Tsft=${Tsft} --randSeed=1 --timestampsFiles='./H1_test_timestamps.txt,./L1_test_timestamps.txt' --FstatMethod=ResampBest  --outputFstat=$outfile_Resamp_otfn --outputTiming=$timefile_Resamp_otfn  --outputLoudest=${loudest_Resamp_otfn}"
echo $cmdline;
if ! eval "$cmdline"; then
    echo "Error.. something failed when running '$cfs_code' ..."
    exit 1;
fi

echo "----------------------------------------"
echo " STEP 6: Comparing results: "
echo "----------------------------------------"

sort -o $outfile_Comp $outfile_Comp
sort -o $outfile_Demod $outfile_Demod
sort -o $outfile_Resamp $outfile_Resamp
sort -o $outfile_Resamp_otfn $outfile_Resamp_otfn

## compare absolute differences instead of relative, allow deviations of up to sigma=sqrt(8)~2.8
echo
cmdline="$cmp_code -1 ./${outfile_Demod} -2 ./${outfile_Comp} --tol-L1=8e-2 --tol-L2=8e-2 --tol-angle=0.08 --tol-atMax=8e-2"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong with injectionSources..."
    exit 2
else
    echo "	==> OK."
fi

echo
cmdline="$cmp_code -1 ./${outfile_Demod} -2 ./${outfile_Resamp} --tol-L1=8e-2 --tol-L2=8e-2 --tol-angle=0.08 --tol-atMax=8e-2"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong with Resamp/Demod..."
    exit 2
else
    echo "	==> OK."
fi

#laxer tolerances as we have different noiuse realizations in this testcase. Incuded just for smoke test
echo
cmdline="$cmp_code -1 ./${outfile_Resamp} -2 ./${outfile_Resamp_otfn} --tol-L1=7e-1 --tol-L2=7e-1 --tol-angle=0.2 --tol-atMax=2e-1"
echo -n $cmdline
if ! eval $cmdline; then
    echo "==> OUCH... files differ. Something might be wrong with --injectSqrtSX et al."
    exit 2
else
    echo "      ==> OK."
fi

echo
echo "----------------------------------------------------------------------"
echo " STEP 7: Sanity-check Resampling parameter estimation: "
echo "----------------------------------------------------------------------"
echo

if grep -q 'nan;' ${loudest_Resamp}; then
    echo "ERROR: ${loudest_Resamp} contains NaNs!"
    exit 2
fi

esth0=$(grep '^h0' ${loudest_Resamp} | awk -F '[ ;]*' '{print $3}')
estdh0=$(grep '^dh0' ${loudest_Resamp} | awk -F '[ ;]*' '{print $3}')

estPhi0=$(grep '^phi0' ${loudest_Resamp} | awk -F '[ ;]*' '{print $3}')
estdPhi0=$(grep '^dphi0' ${loudest_Resamp} | awk -F '[ ;]*' '{print $3}')

echo "Estimated h0 = $esth0 +/- $estdh0: Injected h0 = $h0"
h0inrange=$(echo $h0 $esth0 $estdh0 | awk '{printf "%i\n", (($1 - $2)^2 < 3 * $3^2)}')
if test x$h0inrange != x1; then
    echo "ERROR: estimated h0 was not within 3 x sigma of injected h0!"
    echo
    exit 2
else
    echo "OK: Estimated h0 is within 3 x sigma of injected h0"
    echo
fi

echo "Estimated phi0 = $estPhi0 +/- $estdPhi0: Injected phi0 = $phi0"
phi0inrange=$(echo $phi0 $estPhi0 $estdPhi0 | awk '{printf "%i\n", (($1 - $2)^2 < 3 * $3^2)}')
if test x$phi0inrange != x1; then
    echo "ERROR: estimated phi0 was not within 3 x sigma of injected phi0!"
    echo
    exit 2
else
    echo "OK: Estimated phi0 is within 3 x sigma of injected phi0"
    echo
fi
