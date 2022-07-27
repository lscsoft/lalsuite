#!/bin/bash

# example script for Condor usage on Merlin cluster

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
startfreq=50.0
freqband=0.25

process=$$

# save initial working directory
workdir=`pwd`

#iam=badkri
iam=sintes

# tell us where you are working
echo job $1 running as process $process owned by $iam on `hostname`

# sft filename pattern
#sftdir="/scratch/S4-SFTv2/*SFT*.*"
sftdir="/home/badkri/S4/S4-SFTv2-48-110/*.sft"

# temporary scratch directory where we work in
tempoutdir=/scratch/tmp/$iam/$1.MultiMC-1.run

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf ${tempoutdir}; exit 0" 0 1 2 9 15


# make the temporary directory where we work 
mkdir -p $tempoutdir

# make the temporary output directory where results are stored
mkdir -p ${tempoutdir}/MultiMC

# change to temporary working dir
cd $tempoutdir

# ----------------------------------
# copy files from some central place
# ----------------------------------

rsync -a $workdir/sun05-09.dat .
rsync -a $workdir/earth05-09.dat .
rsync -a $workdir/skyfileS4c .
rsync -a $workdir/S4lines_H1_xavi.txt.v2 .
rsync -a $workdir/S4lines_L1_xavi.txt.v2 .
rsync -a $workdir/S4lines_H2.txt .
rsync -a $workdir/MCInjectHoughMulti .
rsync -a $workdir/Multi_W50_100predictedUL .

echo finished copying stuff

# the frequency band analyzed by the i^{th} node is "startfreq+freqband*i"
temp=`mult $1 $freqband`
freq=`add $temp $startfreq` 
echo frequency range analysed is from $freq Hz with bandwidth $freqband Hz

#input file with predicted UL
kk=./Multi_W50_100predictedUL
firstline=1
temp=`add $1 $firstline`

# choose range for h0
h0rangeinfo=`tail --lines=+$temp ${kk} | head --lines=1`
#freq=`echo $h0rangeinfo | awk '{print $1}'`
h0ini=`echo $h0rangeinfo | awk '{print $4}'`
lowh0=`mult $h0ini 0.9`
highh0=`mult $h0ini 1.1`

echo range of h0 is from $lowh0 to $highh0

sleeptime=`mult $1 20`
sleep $sleeptime


# now run the MC
./MCInjectHoughMulti -d 0 --f0=$freq --fSearchBand=$freqband  --skyfile=./skyfileS4c --weighAM=1 --weighNoise=1 --earthEphemeris=./earth05-09.dat --sunEphemeris=./sun05-09.dat --sftDir=${sftdir} --dirnameOut=${tempoutdir}/MultiMC --fnameout=/MCfreq_${freq} --linefiles=./S4lines_H1_xavi.txt.v2,./S4lines_L1_xavi.txt.v2,./S4lines_H2.txt --printLog --nMCloop=500 --h0Min=$lowh0 --h0Max=$highh0 --nh0=10 --fast=1

echo finished running injections

# now copy the results to the home directory
mkdir -p ${workdir}/MultiMC
rsync -a ${tempoutdir}/MultiMC $workdir

echo done!
