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
startdir=`pwd`

iam=badkri

# tell us where you are working
echo job $1 running as process $process owned by $iam on `hostname`

# choose detector 
# L is L1 and H is H1
det=L1
#det=H1
#det=H2

# the detector option in the hough driver is a number and choose harmonicsfile
# 2 is L1 and 3 is H1
if [ $det = L1 ]; then
    detnum=2
    harmonicsfile=harmonicsS2LLO4K_200_400.txt
fi
if [ $det = H1 ]; then
    detnum=3
    harmonicsfile=harmonicsS2LHO4K_200_400.txt
fi
if [ $det = H2 ]; then
    detnum=3
    harmonicsfile=harmonicsS2LHO2K_200_400.txt
fi


# the directory where the sfts are located
#sftdir=/sft/S2-LIGO/S2_${det}_Funky-v3Cal30MinSFTs
sftdir=/scratch/tmp/badkri/${det}sfts


# set a trap so if the job ends some cleanup is done
trap "/bin/rm -rf /scratch/tmp/$iam/$1.run; exit 0" 0 1 2 9 15

# create temporary scratch directory where we work in 
tempworkdir=/scratch/tmp/${iam}/$1.run
mkdir -p $tempworkdir

# change to temporary working dir
cd $tempworkdir

# create temporary directory to write the output
tempoutdir=${tempworkdir}/${det}
mkdir -p $tempoutdir

# ----------------------------------
# copy files from some central place
# ----------------------------------

rsync -a $startdir/sun00-04.dat .
rsync -a $startdir/earth00-04.dat .
rsync -a $startdir/MCInjectS2 .
rsync -a $startdir/${det}_h0range.run .
rsync -a ${startdir}/${harmonicsfile} .
echo finished copying stuff

# choose freq band
startfreq=200.0
freq=`add $1 $startfreq` 
band=1.0
echo frequency range is from $freq Hz with bandwidth $band Hz

# choose range for h0
h0rangeinfo=`./${det}_h0range.run $freq`
lowh0=`echo $h0rangeinfo | awk '{print $1}'`
highh0=`echo $h0rangeinfo | awk '{print $2}'`
echo range of h0 is from $lowh0 to $highh0

# now run MCInject
./MCInjectS2 -d 0 -i ${detnum} -f $freq -b 1.0 -E ./earth00-04.dat -S ./sun00-04.dat -D $sftdir -o $tempoutdir/MC_$det_$freq -N 1000 -m $lowh0 -M $highh0 -n 10 -H $harmonicsfile

echo finished running injections

# create the output directory
mkdir -p $startdir/$det
rsync -a $tempoutdir $startdir

echo done!
