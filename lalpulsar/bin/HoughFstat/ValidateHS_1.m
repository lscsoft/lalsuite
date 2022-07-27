#! /usr/bin/octave -q
## Another script for comparing HierarchicalSearch.c with CFSv2.
## This is more general that ValidateHS.c because it goes through 
## the results produved by HS 1-by-1.  Thus it can deal with the case
## when we use multiple stacks

fStart=310;
fBand=0.1;
fdot=0;
fdotBand=1.0e-9;
nStacks1=10;
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

#[output,status] = system(cmdline)


## load output of HS
load out/HS_fstatVec1.txt
F_HS = 2*HS_fstatVec1(:,5);
freq_HS = HS_fstatVec1(:,1);
alpha_HS = HS_fstatVec1(:,2);
delta_HS = HS_fstatVec1(:,3);
fdot_HS = HS_fstatVec1(:,4);



cmdline = sprintf("ComputeFStatistic_v2 --Freq=%.12g --f1dot=%.12g \
--FreqBand=%.12g --f1dotBand=%.12g --dFreq=%.12g --df1dot=%.12g \
--minStartTime=%d --maxEndTime=%d --refTime=%d \
--DataFiles=%s --skyGridFile=%s \
--ephemDir=/local_data/badkri/lscsoft/share/lal/ --ephemYear=05-09 \
    --TwoFthreshold=0 --gridType=3 --outputLabel=CFSv2 --outputFstat=CFSv2 \
    ", fStart, fdot, fBand, fdotBand, dFreq, df1dot, minStartTime, \
 		  maxEndTime, refTime, DataFiles, skyGridFile)


#for index=1:length(Freq) 
for index=1:length(F_HS)
 
  cmdline = sprintf("ComputeFStatistic_v2 --Freq=%.12g --f1dot=%.12g \
  --FreqBand=%.12g --f1dotBand=%.12g --dFreq=%.12g --df1dot=%.12g \
  --dAlpha=%.12g  --dDelta=%.12g --minStartTime=%d --maxEndTime=%d --refTime=%d \
  --DataFiles=%s --Alpha=%.12g --Delta=%.12g --AlphaBand=0 --DeltaBand=0 \
  --ephemDir=/local_data/badkri/lscsoft/share/lal/ --ephemYear=05-09 \
      --TwoFthreshold=0 --gridType=0 --outputLabel=CFSv2 \
      --outputFstat=CFSv2 ", freq_HS(index), fdot_HS(index), 0.1, 0, dFreq, df1dot, \
		    0.2, 0.2, minStartTime, maxEndTime, refTime, DataFiles, \
		    alpha_HS(index), delta_HS(index));
  
  [output, status] = system(cmdline);

  index

  load CFSv2
    
  F_CFSv2(index) = CFSv2(1 ,7);

  clear CFSv2 

endfor

diff=F_CFSv2 - F_HS;
