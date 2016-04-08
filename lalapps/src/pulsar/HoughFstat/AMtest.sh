#!/bin/bash

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

#process=$$

# save initial working directory
startdir=`pwd`

h0=3.5e-24
mismatch=0.4
sftdir=/local_data/badkri/fakesfts/
#th=4.0
outfile=AM_4_test_$h0
#echo $outfile
#cosiota=1.0

#sleep 10

# choose sky location
#alpha=0.0
#delta=-1.570796327
#pi2=1.570796327
# number of times signal is injected

top=1000
j=0

# loop $top times over code
while [ $j -lt $top ] ; do

        echo $j
        #------------------------------------------------------
        # create random input parameters within chosen ranges
        #------------------------------------------------------

        # output of allangles are seven random numbers in these ranges:
        # (0,2*pi) (-pi/2,pi/2) (0,2*pi) (-1,1) (0,2*pi) (0,1) (0,1)
        # they are used respectively to find random values of
        # alpha, delta, phi0, cosiota, psi, freq, spindown
        allangles=`./makeangles`

        # first create psi, cosiota and phi0
        alpha=`echo $allangles | awk '{print $1}'`
        delta=`echo $allangles | awk '{print $2}'`
	cosiota=`echo $allangles | awk '{print $4}'`
	psi=`echo $allangles | awk '{print $5}'`
	# if we want cosiota random
	./HoughValidateAM -E ./earth00-19-DE405.dat -S ./sun00-19-DE405.dat -D $sftdir -r $alpha -l $delta -m $h0 -c $cosiota -M $mismatch >> $outfile

	# if we want cosiota at a fixed value
	#./HoughValidateAM -E ./earth00-19-DE405.dat -S ./sun00-19-DE405.dat -D $sftdir -r $alpha -l $delta -m $h0 -c 0.5 -p $psi -M $mismatch >> $outfile

	#sleep 30
        let j+=1
done
