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

# needs to be changed
#sftdir=/sft/S4-LIGO/sft_1800.20060227.S4_v4/S4-${det}.1800-sft
sftdir=/home/badkri/S4/SFTv2/

# temporary scratch directory where we work in
tempoutdir=/scratch/tmp/$iam/$1.run

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf ${tempoutdir}; exit 0" 0 1 2 9 15

# make the temporary directory where we work 
mkdir -p $tempoutdir

# make the temporary output directory where results are stored
mkdir -p $tempoutdir/$det

# change to temporary working dir
cd $tempoutdir

# -----------------------------------
# copy files from some central place
# -----------------------------------

rsync -a $workdir/sun05-09.dat .
rsync -a $workdir/earth05-09.dat .
rsync -a $workdir/DriveHough_v3 .
rsync -a $workdir/skypatchfile .
echo finished copying stuff

# choose start of freq range to be analyzed
startfreq=200.0

# the 1Hz frequency band analyzed by the i^{th} node is "startfreq+i"
freq=`add $1 $startfreq` 
echo frequency analysed is 1Hz band starting from $freq

# now run the driver
./DriveHough_v3 -d 0 --ifo=$det --f0=$freq --fSearchBand=1.0 --skyfile=./skypatchfile --weighAM --weighNoise --earthEphemeris=./earth05-09.dat --sunEphemeris=./sun05-09.dat --sftDir=$sftdir --dirnameOut=$tempoutdir/$det --fbasenameOut=freq_${freq}_  --printStats --printSigma

echo finished running driver

# now copy the results to the home directory
mkdir -p $workdir/$det
rsync -a $tempoutdir/$det $workdir

echo done!  

