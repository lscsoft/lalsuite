# Test the generation of all available waveform approximants by using
# the GenerateSimulation test code

m1=20.
m2=10.
fmin=10.
s1x=0.1
s2x=0.2
s1y=0.4
s2y=0.3
s1z=0.5
s2z=0.5
lam1=500.
lam2=400.
srate=16384.

# Test TD non-spinning, non-tidal approximants
dom="TD"
approx="TaylorEt EOBNRv2 EOBNRv2HM IMRPhenomA"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done


# Test TD aligned-spin, non-tidal approximants
approx="SEOBNRv1 IMRPhenomB"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1z $s1z --spin2z $s2z --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test TD non-spinning, tidal approximants
approx="TaylorT1 TaylorT2 TaylorT3 TaylorT4"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --tidal-lambda1 $lam1 --tidal-lambda2 $lam2 --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test TD precessing, tidal approximants
approx="SpinTaylorT2 SpinTaylorT4"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1x $s1x --spin1y $s1y --spin1z $s1z --spin2x $s2x --spin2y $s2y --spin2z $s2z --tidal-lambda1 $lam1 --tidal-lambda2 $lam2 --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test TD precessing, non-tidal approximants
# FIXME: Disable until PhenSpin is fixed
#approx="PhenSpinTaylor PhenSpinTaylorRD"

#for i in $approx;
#do
#    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1x $s1x --spin1y $s1y --spin1z $s1z --spin2x $s2x --spin2y $s2y --spin2z $s2z --f-min $fmin --phase-order -1 --amp-order -1"
#    echo Running: $cmd
#    $cmd
#    if [ $? = 0 ]; then
#        echo Success.
#    else
#        echo Error occurred. Return code: $?
#        exit 1
#    fi
#done

# Test FD non-spinning, non-tidal approximants
dom="FD"
approx="IMRPhenomA"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test FD aligned-spin, non-tidal approximants
approx="IMRPhenomB IMRPhenomC TaylorF2RedSpin"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1z $s1z --spin2z $s2z --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test FD aligned-spin, tidal approximants
approx="TaylorF2 TaylorF2RedSpinTidal"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1z $s1z --spin2z $s2z --tidal-lambda1 $lam1 --tidal-lambda2 $lam2 --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done



# Test FD precessing, non-tidal approximants
approx="IMRPhenomP"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1x $s1x --spin1y $s1y --spin1z $s1z --spin2x $s2x --spin2y $s2y --spin2z $s2z --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done

# Test FD single-spin precessing, non-tidal approximants
approx="SpinTaylorF2"

for i in $approx;
do
    cmd="./GenerateSimulation --domain $dom --approximant $i --outname simulation-$dom-$i.dat --m1 $m1 --m2 $m2 --spin1x $s1x --spin1y $s1y --spin1z $s1z --f-min $fmin --phase-order -1 --amp-order -1"
    echo Running: $cmd
    $cmd
    if [ $? = 0 ]; then
        echo Success.
    else
        echo Error occurred. Return code: $?
        exit 1
    fi
done
