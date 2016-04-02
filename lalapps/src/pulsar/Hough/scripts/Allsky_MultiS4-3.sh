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

# choose start of freq range to be analyzed
startfreq=850.0
freqband=0.5

process=$$

# save initial working directory
workdir=`pwd`

iam=badkri

# tell us where you are working
echo job $1 running as process $process owned by $iam on `hostname`

# sft filename pattern
sftdir="/scratch/S4-SFTv2/*.sft"

# temporary scratch directory where we work in
tempoutdir=/scratch/tmp/$iam/$1.Multi-3.run

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf ${tempoutdir}; exit 0" 0 1 2 9 15

# make the temporary directory where we work 
mkdir -p $tempoutdir

# make the temporary output directory where results are stored
mkdir -p ${tempoutdir}/Multi

# change to temporary working dir
cd $tempoutdir

# -----------------------------------
# copy files from some central place
# -----------------------------------

rsync -a $workdir/sun05-09.dat .
rsync -a $workdir/earth05-09.dat .
rsync -a $workdir/DriveHoughMulti .
rsync -a $workdir/skyfileS4c .
rsync -a $workdir/S4lines_H1_xavi.txt.v2 .
rsync -a $workdir/S4lines_L1_xavi.txt.v2 .
rsync -a $workdir/S4lines_H2.txt .

echo finished copying stuff

# the 1Hz frequency band analyzed by the i^{th} node is "startfreq+i"
temp=`mult $1 $freqband`
freq=`add $temp $startfreq` 
echo frequency analysed is 1Hz band starting from $freq

# now run the driver
./DriveHoughMulti -d 0 --f0=$freq --fSearchBand=$freqband --skyfile=./skyfileS4c --weighAM=1 --weighNoise=1 --earthEphemeris=./earth05-09.dat --sunEphemeris=./sun05-09.dat --sftDir=$sftdir --dirnameOut=${tempoutdir}/Multi --fbasenameOut=freq_${freq}_  --linefiles=./S4lines_H1_xavi.txt.v2,./S4lines_L1_xavi.txt.v2,./S4lines_H2.txt --printLog

echo finished running driver

# now copy the results to the home directory
mkdir -p ${workdir}/Multi
rsync -a ${tempoutdir}/Multi $workdir

echo done!  

