#!/bin/bash

# example script for to generate the bandlist into an ascii file

# -----------------------------------
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

# -----------------------------------


nbands=172
let nbands+=1
j=1

while [ $j -lt $nbands ] ; do
   # choose freq band
  ./L1_freqbands.run $j >> L1bandlist
  let j+=1
done
