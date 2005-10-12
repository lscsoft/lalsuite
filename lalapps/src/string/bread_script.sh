#!/bin/sh

rm -rf H1triggers/;mkdir H1triggers/ >& /dev/null

for file in triggers/H1* ; do cp $file H1triggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H1triggers.xml --globtrig H1triggers --cluster stringcluster
echo H1 single IFO triggers:
~/ligotools/bin/lwtprint H1triggers.xml -t sngl_burst |wc -l

rm -rf H2triggers/;mkdir H2triggers/ >& /dev/null

for file in triggers/H2* ; do cp $file H2triggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H2triggers.xml --globtrig H2triggers  --cluster stringcluster
echo H2 single IFO triggers:
~/ligotools/bin/lwtprint H2triggers.xml -t sngl_burst |wc -l

rm -rf L1triggers/;mkdir L1triggers/ >& /dev/null

for file in triggers/L1* ; do cp $file L1triggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output L1triggers.xml --globtrig L1triggers --cluster stringcluster 
echo L1 single IFO triggers:
~/ligotools/bin/lwtprint L1triggers.xml -t sngl_burst |wc -l

rm -rf H1ctriggers/;mkdir H1ctriggers/ >& /dev/null

for file in coincidences/H1*H1L1*M* ; do cp $file H1ctriggers ; done
for file in coincidences/H1*H1L1*P* ; do cp $file H1ctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H1ctriggers.xml --globtrig H1ctriggers  --cluster stringcluster
echo H1 triple coincidence survivors:
~/ligotools/bin/lwtprint H1ctriggers.xml -t sngl_burst |wc -l

rm -rf H2ctriggers/;mkdir H2ctriggers/ >& /dev/null

for file in coincidences/H2*H2L1*M* ; do cp $file H2ctriggers ; done
for file in coincidences/H2*H2L1*P* ; do cp $file H2ctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H2ctriggers.xml --globtrig H2ctriggers  --cluster stringcluster
echo H2 triple coincidence survivors:
~/ligotools/bin/lwtprint H2ctriggers.xml -t sngl_burst |wc -l

rm -rf L1ctriggers/;mkdir L1ctriggers/ >& /dev/null

for file in coincidences/L1*L1H1*M* ; do cp $file L1ctriggers ; done
for file in coincidences/L1*L1H1*P* ; do cp $file L1ctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output L1ctriggers.xml --globtrig L1ctriggers  --cluster stringcluster
echo L1 triple coincidence survivors:
~/ligotools/bin/lwtprint L1ctriggers.xml -t sngl_burst |wc -l


rm -rf H1dctriggers/;mkdir H1dctriggers/ >& /dev/null

for file in doublecoincidencesH1L1/H1*H1L1*M* ; do cp $file H1dctriggers ; done
for file in doublecoincidencesH1L1/H1*H1L1*P* ; do cp $file H1dctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H1dctriggers.xml --globtrig H1dctriggers  --cluster stringcluster
echo H1-L1 double coincidence survivors:
~/ligotools/bin/lwtprint H1dctriggers.xml -t sngl_burst |wc -l

rm -rf L1dctriggers/;mkdir L1dctriggers/ >& /dev/null

for file in doublecoincidencesH1L1/L1*L1H1*M* ; do cp $file L1dctriggers ; done
for file in doublecoincidencesH1L1/L1*L1H1*P* ; do cp $file L1dctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output L1dctriggers.xml --globtrig L1dctriggers  --cluster stringcluster
echo L1-H1 double coincidence survivors:
~/ligotools/bin/lwtprint L1dctriggers.xml -t sngl_burst |wc -l

rm -rf H1H2dctriggers/;mkdir H1H2dctriggers/ >& /dev/null

for file in coincidences/H1*H1H2*M* ; do cp $file H1H2dctriggers ; done
for file in coincidences/H1*H1H2*P* ; do cp $file H1H2dctriggers ; done

~/lscsoft/lalapps/src/power/lalapps_bread  --output H1H2dctriggers.xml --globtrig H1H2dctriggers  --cluster stringcluster
echo H1-H2 double coincidence survivors:
~/ligotools/bin/lwtprint H1H2dctriggers.xml -t sngl_burst |wc -l

