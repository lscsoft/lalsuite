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

# choose start of freq range to be analyzed 193.95-0.05
startfreq=193.90
freqband=0.1

process=$$

# save initial working directory
workdir=`pwd`

#iam=badkri
iam=sintes

# tell us where you are working
echo job $1 running as process $process owned by $iam on `hostname`

# sft filename pattern
#sftdir="/scratch/S4-SFTv2/H-1_H1*.sft"
sftdir="/scratch/S4-SFTv2/L-1_L1*.sft"
#sftdir="/home/badkri/S4/S4-SFTv2-48-110/L-1_L1*.sft"
#sftdir="/scratch/S4-SFTv2/*SFT*.*"

# temporary scratch directory where we work in
tempoutdir=/scratch/tmp/$iam/$1.L1INJECT.run

# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf ${tempoutdir}; exit 0" 0 1 2 9 15

# make the temporary directory where we work 
mkdir -p $tempoutdir

# make the temporary output directory where results are stored
mkdir -p ${tempoutdir}/L1.P8

# change to temporary working dir
cd $tempoutdir

# -----------------------------------
# copy files from some central place
# -----------------------------------

rsync -a $workdir/sun05-09.dat .
rsync -a $workdir/earth05-09.dat .
rsync -a $workdir/DriveHoughMulti .
rsync -a $workdir/S4lines_H1_xavi.txt.v2 .
rsync -a $workdir/S4lines_L1_xavi.txt.v2 .
rsync -a $workdir/S4lines_H2.txt .
rsync -a $workdir/ts-inject-S4_L1.txt .
rsync -a $workdir/ts-inject-S4_H1.txt .
rsync -a $workdir/ts-inject-S4_H2.txt .
rsync -a $workdir/ts-inject-S4_Multi.txt .
rsync -a $workdir/skypulsar.8 .

echo finished copying stuff


#temp=`mult $1 $freqband`
#freq=`add $temp $startfreq` 
freq=$startfreq
echo frequency analysed is $freqband Hz band starting from $freq

#sleeptime=`mult $1 50`
#sleep $sleeptime

# now run the driver
# with AM and noise weights
./DriveHoughMulti -d 0 --f0=$freq --fSearchBand=$freqband --nfSizeCylinder=45 --timeStampsFile=./ts-inject-S4_L1.txt  --skyfile=./skypulsar.8 --weighAM=1 --weighNoise=1 --printMaps=1 --printStats=1 --earthEphemeris=./earth05-09.dat --sunEphemeris=./sun05-09.dat --sftDir=$sftdir --dirnameOut=${tempoutdir}/L1.P8 --fbasenameOut=freq_${freq}_  --linefiles=./S4lines_H1_xavi.txt.v2,./S4lines_L1_xavi.txt.v2,./S4lines_H2.txt --printLog

echo finished running driver

# now copy the results to the home directory
mkdir -p ${workdir}/L1.P8
rsync -a ${tempoutdir}/L1.P8 $workdir

echo done!  

