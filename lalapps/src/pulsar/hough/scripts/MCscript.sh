#!/bin/bash

# example script for Condor usage on Merlin cluster
# based on scripts written by Marialessandra Papa, Badri Krishnan and Steffen Grunewald
# last modifications: 2003-11-13

# define some useful inline functions to do floating point arithmetics
# beware that awk is slow!
# first two are to be run like "if testxy ; then" : 0 is OK, 1 is bad
# next are elementary float point maths
# if you've got more complex formulae please use awk directly as shown below
testlt(){
return `echo $1 $2 | awk '{if($1<$2){print 0}else{print 1}}'`
}
testgt(){
return `echo $1 $2 | awk '{if($1>$2){print 0}else{print 1}}'`
}
add(){
echo $1 $2 | awk '{print $1+$2}'
}
sub (){
echo $1 $2 | awk '{print $1-$2}'
}
mult(){
echo $1 $2 | awk '{print $1*$2}'
}
div(){
echo $1 $2 | awk '{print $1/$2}'
}

process=$$

# save initial working directory
workdir=`pwd`

iam=badkri

# tell us where you are working
echo job $1 running as process $process owned by $iam on `hostname`

# choose detector 
# L is L1 and H is H1
det=L1
#det=H1
#det=H2

# the detector option in the hough driver is a number
# 2 is L1 and 3 is H1
if [ $det = L1 ]; then
    detnum=2
fi
if [ $det = H1 ]; then
    detnum=3
fi
if [ $det = H2 ]; then
    detnum=3
fi

mkdir -p $workdir/$det
outdir=$workdir/$det/MC_${det}_$1
sftdir=/sft/S2-LIGO/S2_${det}_Funky-v3Cal30MinSFTs

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf /scratch/tmp/$iam/$1.run; exit 0" 0 1 2 9 15

mkdir -p /scratch/tmp/$iam/$1.run
cd /scratch/tmp/$iam/$1.run

# ----------------------------------
# copy files from some central place
# ----------------------------------

cp -f $workdir/sun00-04.dat .
cp -f $workdir/earth00-04.dat .
cp -f $workdir/MCInjectComputeHough .
cp -f $workdir/${det}_freqbands.run .
cp -f $workdir/${det}_h0range.run .
echo finished copying stuff


# choose freq band
freqindex=`echo $1 1 | awk '{print $1+$2}'`
freqbandinfo=`./${det}_freqbands.run $freqindex`
startfreq=`echo $freqbandinfo | awk '{print $1}'`
endfreq=`echo $freqbandinfo | awk '{print $2}'`
band=`echo $endfreq $startfreq | awk '{print $1-$2}'`
echo frequency range is from $startfreq Hz with bandwidth $band Hz

# choose range for h0
h0rangeinfo=`./${det}_h0range.run $freqindex`
lowh0=`echo $h0rangeinfo | awk '{print $1}'`
highh0=`echo $h0rangeinfo | awk '{print $2}'`
echo range of h0 is from $lowh0 to $highh0


# now run MCInject
./MCInjectComputeHough -d 0 -i ${detnum} -E ./earth00-04.dat -S ./sun00-04.dat -D $sftdir -o $outdir -f $startfreq -b $band -H 10 $lowh0 $highh0 -L 10000

echo finished running driver

# cleanup will be done by trap!
