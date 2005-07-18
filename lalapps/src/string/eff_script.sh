#!/bin/sh

#must add a check here for the command line argument

echo L1ctriggers.xml > L1.txt
echo H1ctriggers.xml > H1.txt
echo H2ctriggers.xml > H2.txt
ls -1 HL-* > inj.txt

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig H1.txt --input-burstinj inj.txt \
--output-trig H1cout.xml --output-inj-made H1cinj.xml --output-inj-found H1cinjfound.xml \
--best-peaktime --compare-choice comparebytime  --noplayground --gps-start-time  795169179 \
--gps-end-time 795211620 --printresult --max-confidence $1  >& /dev/null

mv BinjFindResults.dat BinjFindResultsH1triple.dat

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig H2.txt --input-burstinj inj.txt \
--output-trig H2cout.xml --output-inj-made H2cinj.xml --output-inj-found H2cinjfound.xml \
--best-peaktime --compare-choice comparebytime --noplayground --gps-start-time  795169179 \
--gps-end-time 795211620 --printresult --max-confidence $1  >& /dev/null

mv BinjFindResults.dat BinjFindResultsH2triple.dat

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig L1.txt --input-burstinj inj.txt \
--output-trig L1cout.xml --output-inj-made L1cinj.xml --output-inj-found L1cinjfound.xml \
--best-peaktime --compare-choice comparebytime --noplayground --gps-start-time 795169179 \
--gps-end-time 795211620 --printresult --max-confidence $1  >& /dev/null

mv BinjFindResults.dat BinjFindResultsL1triple.dat

echo Triple coincident efficiencies at SNR threshold $1 
echo --------------------------------------------------
echo
echo H1:
cat *H1triple.dat
echo
echo H2:
cat *H2triple.dat
echo
echo L1:
cat *L1triple.dat
echo
