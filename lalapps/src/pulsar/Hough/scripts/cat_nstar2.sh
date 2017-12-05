#!/bin/bash

for i in `seq 800 0.25 899.75`;
  do
  echo $i
  cat /home/badkri/S4/Multi/nstarfiles/freq_${i}_nstar >>  nstarH1_800_900W.txt
done    
