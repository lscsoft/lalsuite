#!/bin/sh

#must add a check here for the command line argument

ls -1 L1triggers/* > L1.txt
ls -1 H1triggers/* > H1.txt
ls -1 H2triggers/* > H2.txt

ls -1 HL-* > inj.txt

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig H1.txt --input-burstinj inj.txt \
--output-trig H1out.xml --output-inj-made H1inj.xml --output-inj-found H1injfound.xml \
--best-peaktime --compare-choice comparebytime  --noplayground --gps-start-time 795169179  \
--gps-end-time 795211620 --printresult --max-confidence $1 >& /dev/null

mv BinjFindResults.dat BinjFindResultsH1single.dat

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig H2.txt --input-burstinj inj.txt \
--output-trig H2out.xml --output-inj-made H2inj.xml --output-inj-found H2injfound.xml \
--best-peaktime --compare-choice comparebytime --noplayground --gps-start-time  795169179 \
--gps-end-time 795211620 --printresult --max-confidence $1  >& /dev/null

mv BinjFindResults.dat BinjFindResultsH2single.dat

~/lscsoft/lalapps/src/power/lalapps_binj_find --input-trig L1.txt --input-burstinj inj.txt \
--output-trig L1out.xml --output-inj-made L1inj.xml --output-inj-found L1injfound.xml \
--best-peaktime --compare-choice comparebytime --noplayground --gps-start-time 795169179 \
--gps-end-time 795211620 --printresult --max-confidence $1  >& /dev/null

mv BinjFindResults.dat BinjFindResultsL1single.dat

echo Single IFO efficiencies at SNR threshold $1
echo --------------------------------------------
echo
echo H1:
cat *H1single.dat
echo
echo H2:
cat *H2single.dat
echo
echo H1:
cat *L1single.dat
echo

ls -1 L1ctriggers/* > L1.txt
ls -1 H1ctriggers/* > H1.txt
ls -1 H2ctriggers/* > H2.txt

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
echo H1:
cat *L1triple.dat
echo
