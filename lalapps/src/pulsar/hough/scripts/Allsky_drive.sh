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
#det=L1
#det=H1
det=H2

mkdir -p $workdir/$det
sftdir=/sft/S2-LIGO/S2_${det}_Funky-v3Cal30MinSFTs
# the detector option in the hough driver is a number
# 2 is L1 and 3 is H1 or H2
if [ $det = L1 ]; then
    detnum=2
fi
if [ $det = H1 ]; then
    detnum=3
fi
if [ $det = H2 ]; then
    detnum=3
fi

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf /scratch/tmp/$iam/$1.run; exit 0" 0 1 2 9 15

# make the temporary directory where we work 
mkdir -p /scratch/tmp/$iam/$1.run
cd /scratch/tmp/$iam/$1.run

# -----------------------------------
# copy files from some central place
# -----------------------------------

cp -f $workdir/sun00-04.dat .
cp -f $workdir/earth00-04.dat .
cp -f $workdir/DriveHoughColor_velo .
cp -f $workdir/skypatches.run .
echo finished copying stuff

# choose start of freq band to be analyzed
startfreq=200.0

# the 1Hz frequency band analyzed by the i^{th} node is "startfreq+i"
freq=`add $1 $startfreq` 
echo frequency analysed is 1Hz band starting from $freq


# loop over the skypatches
j=1
while [ $j -lt 24 ] ; do

    # set skypatch parameters    
    skypatchinfo=`./skypatches.run $j`
    alpha=`echo $skypatchinfo | awk '{print $1}'`
    delta=`echo $skypatchinfo | awk '{print $2}'`
    sizealpha=`echo $skypatchinfo | awk '{print $3}'`
    sizedelta=`echo $skypatchinfo | awk '{print $4}'`

    # make the output directory to store results
    outdir=$workdir/$det/skypatch_$j
    mkdir -p $outdir

    # now run the driver
    ./DriveHoughColor_velo -d 0 -i ${detnum} -w 25 -E ./earth00-04.dat -S ./sun00-04.dat -D $sftdir -o $outdir/freq_$1_ -p $alpha $delta -s $sizedelta $sizealpha -f $freq -b 1.0
    echo finished running driver for skypatch $j
    
    let j+=1
done 
