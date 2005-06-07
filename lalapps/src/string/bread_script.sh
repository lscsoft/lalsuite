#!/bin/sh

mkdir H1triggers/ >& /dev/null
cp triggers/H1* H1triggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output H1triggers.xml --globtrig H1triggers 

echo H1 single IFO triggers:
lwtprint H1triggers.xml -t sngl_burst |wc -l

mkdir H2triggers/ >& /dev/null
cp triggers/H2* H2triggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output H2triggers.xml --globtrig H2triggers

echo H2 single IFO triggers:
lwtprint H2triggers.xml -t sngl_burst |wc -l

mkdir L1triggers/ >& /dev/null
cp triggers/L1* L1triggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output L1triggers.xml --globtrig L1triggers

echo L1 single IFO triggers:
lwtprint L1triggers.xml -t sngl_burst |wc -l

mkdir H1ctriggers/ >& /dev/null
cp coincidences/H1*H1H2* H1ctriggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output H1ctriggers.xml --globtrig H1ctriggers

echo H1 triple coincidence survivors:
lwtprint H1ctriggers.xml -t sngl_burst |wc -l

mkdir H2ctriggers/ >& /dev/null
cp coincidences/H2*H2L1* H2ctriggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output H2ctriggers.xml --globtrig H2ctriggers

echo H2 triple coincidence survivors:
lwtprint H2ctriggers.xml -t sngl_burst |wc -l

mkdir L1ctriggers/ >& /dev/null
cp coincidences/L1*L1H1* L1ctriggers/
~/lscsoft/lalapps/src/power/lalapps_bread  --output L1ctriggers.xml --globtrig L1ctriggers

echo L1 triple coincidence survivors:
lwtprint L1ctriggers.xml -t sngl_burst |wc -l
