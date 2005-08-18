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

h0=5.0e-24
sftdir=/local_data/sintes/fakesfts/
th=1.6
outfile=AMtest

# choose sky location
alpha=0.0
delta=-1.570796327
pi2=1.570796327
# number of times signal is injected

top=5
j=0

# loop $top times over code
while [ $j -lt $top ] ; do
        #------------------------------------------------------
        # create random input parameters within chosen ranges
        #------------------------------------------------------

        # output of allangles are seven random numbers in these ranges:
        # (0,2*pi) (-1,1) (0,2*pi) (-1,1) (-1,1) (0,1) (0,1)
        # they are used respectively to find random values of
        # psi, cosiota, phi0, alpha, delta, frequency, spindown
        allangles=`./makeangles`

        # first create psi, cosiota and phi0
        alpha=`echo $allangles | awk '{print $1}'`
        delta=`echo $allangles | awk '{print $2}'`
	./HoughValidateAM -E ./earth00-04.dat -S ./sun00-04.dat -D $sftdir -r $alpha -l $delta -t $th -m $h0 >> $outfile

        let j+=1
done
