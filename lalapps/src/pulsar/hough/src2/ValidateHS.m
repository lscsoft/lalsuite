## Octave script for comparing the 1st stage Fstat calculation of 
## HierarchicalSearch.c with ComputeFStatistic_v2. We run
## HierarchicalSearch using 1 stack and compare the outputs on the same
## parameter space.  The
## reference time is allowed to be different from the start time of the
## first sft. 


## Set the parameters

fStart=310;
fBand=0.1;
fdot=0;
fdotBand=1.0e-9;
nStacks1=1;
maxEndTime=820631477; 
minStartTime=0;
refTime=600000000;
DataFiles="/local_data/badkri/fakesfts/H-1_H1*.sft";
skyGridFile="./skygrid";
dFreq=0.000055555555555555555;
df1dot=1.5432e-09;

## run HierarchicalSearch.c
cmdline = sprintf("HierarchicalSearch --sftData=%s --followUp=0 \
--printFstat1=1 \
    --ephemE=/local_data/badkri/lscsoft/share/lal/earth05-09.dat \
    --ephemS=/local_data/badkri/lscsoft/share/lal/sun05-09.dat \
    --fStart=%.12g --fBand=%.12g --fdot=%.12g --fdotBand=%.12g \
--nStacks1=%d --skyGridFile=%s --gridType=3 --blocksRngMed=50 \
--Dterms=16 --refTime=%d --minStartTime=%d --maxEndTime=%d", \
		  DataFiles, fStart, fBand, fdot, fdotBand, nStacks1, \
		  skyGridFile, refTime, minStartTime, maxEndTime); 

[output,status] = system(cmdline)



## run ComputeFStatistic_v2
cmdline = sprintf("ComputeFStatistic_v2 --Freq=%.12g --f1dot=%.12g \
--FreqBand=%.12g --f1dotBand=%.12g --dFreq=%.12g --df1dot=%.12g \
--minStartTime=%d --maxEndTime=%d --refTime=%d \
--DataFiles=%s --skyGridFile=%s \
--ephemDir=/local_data/badkri/lscsoft/share/lal/ --ephemYear=05-09 \
--TwoFthreshold=0 --gridType=3 --outputLabel=CFSv2 --outputFstat=CFSv2 \
", fStart, fdot, fBand, fdotBand, dFreq, df1dot, minStartTime, \
		  maxEndTime, refTime, DataFiles, skyGridFile);


load out/HS_fstatVec1.txt

load CFSv2

F1 = CFSv2(:,7);
F2 = 2*HS_fstatVec1(:,5);
diff = 100*(F1 - F2)./F1;

mean(diff)
std(diff)
max(diff)
min(diff)

plot(diff)

clear


