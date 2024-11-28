#!/bin/bash

for i in `seq 50 0.25 99.75`;
  do
  echo $i
  cat /home/badkri/S4/H1-50-100/nstarfiles/freq_${i}_nstar >> nstarH1_50_100.txt
done    
