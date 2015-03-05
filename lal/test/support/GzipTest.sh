#!/bin/bash

builddir="./";

testexec=${builddir}GzipTest

filetext="This is some test text"
fileout=testzip.txt

# create a text file
echo $filetext > $fileout

if [ $? -ne "0" ]; then
  echo "Could not create file $fileout"
  exit 2
fi

# copy the file to a temporary place for comparison at the end
fileoutcpy="${fileout}_copy"
cp $fileout $fileoutcpy

# try zipping the file
$testexec -f $fileout -g

if [ $? -ne "0" ]; then
  echo "Error gzipping file!"
  rm $fileout
  rm $fileoutcpy
  exit 2
fi

# try unzipping the file
fileoutzipped="${fileout}.gz"
$testexec -f $fileoutzipped -u

if [ $? -ne "0" ]; then
  echo "Error guzipping file!"
  rm $fileoutzipped
  rm $fileoutcpy
  exit 2
fi

# test that the zipped and unzipped file is the same as the input
diffcheck=`diff -s $fileout $fileoutcpy`

if [ "$diffcheck" != "Files $fileout and $fileoutcpy are identical" ]; then
  echo "Error, unzipped file is not the same as input"
  rm $fileout
  rm $fileoutcpy
  exit 2
fi

rm $fileout
rm $fileoutcpy

echo "Gzipping and unzipping of a text file is working correctly"

exit 0
