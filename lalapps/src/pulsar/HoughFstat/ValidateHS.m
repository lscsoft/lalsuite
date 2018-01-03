#!/usr/bin/octave -q 
## Octave script for comparing the 1st stage Fstat calculation of 
## HierarchicalSearch.c with ComputeFStatistic_v2. We run
## HierarchicalSearch using 1 stack and compare the outputs on the same
## parameter space.
## Will be generalized to allow for multiple stacks and further validations

## miscellaneous parameters 
nStacks1 = 1;
startTime = 833786244;
duration = 200*1800;
#refTime = 600000000;
refTime = startTime;
##ephemDir = '/local_data/badkri/lscsoft/share/lal/';

## doppler parameter space
fStart = 310;
fBand = 0.01;
fdot = 0;
fdotBand = 0;
#fdotBand=1.0e-9;

## frequency and spindown resolution
dFreq = 1/duration;
df1dot = 1/duration^2

## frequency band for fake sfts 
## -- should be larger than the search band 
mfdfmin = fStart - 0.5;
mfdband = fBand + 1.0;

## injected signal params for makefakedata
signalAlpha = 0;
signalDelta = 0;
signalh0 = 0;
signalcosi = 0;
signalpsi = 0;
signalphi0 = 0;
signalFreq = fStart + 0.25;
signalF1dot = 0;

## cleanup from previous runs
system("rm -rf ./fakesfts");
unlink ("CFSv2");
unlink ("outHS_fstatVec1.dat");

## run makefakedata
system("mkdir -p ./fakesfts");
cmd = sprintf("lalapps_Makefakedata --outSFTbname=./fakesfts/ \
    --IFO=H1 --noiseSqrtSh=1.0e-23\
    --ephemYear=05-09 --fmin=%.12g --Band=%.12g --Alpha=%.12g \
    --Delta=%.12g --h0=%.5e --cosi=%.12g --phi0=%.12g --psi=%.12g \
    --Freq=%.12g --f1dot=%.5e --startTime=%.12g --duration=%.12g", \ 
	      mfdfmin, mfdband, signalAlpha, signalDelta, \
	      signalh0, signalcosi, signalphi0, signalpsi, \
	      signalFreq, signalF1dot, startTime, duration)

[output,status] = system(cmd);
if ( status != 0 )
  error ("Failed to create SFTs! output = '%s'!", output );
endif

## the fake sfts are created 
## -- now run HierarchicalSearch.c and CFSv2

## fake sfts created previously
DataFiles = "./fakesfts/*.sft";

## sky grid file
skyGridFile = "./skygrid";


## run ComputeFStatistic_v2
cmd = sprintf("ComputeFStatistic_v2 --Freq=%.12g --f1dot=%.12g \
--FreqBand=%.12g --f1dotBand=%.12g --dFreq=%.12g --df1dot=%.12g \
--Dterms=8 --refTime=%d --DataFiles='%s' --gridFile=%s \
--ephemYear=05-09 --TwoFthreshold=0 --gridType=3 \
--outputLabel=CFSv2 --outputFstat=CFSv2 --RngMedWindow=101",\ 
	      fStart, fdot, fBand, fdotBand, dFreq, df1dot, \
 	      refTime, DataFiles, skyGridFile )

[output,status] = system(cmd);
if ( status != 0 )
  error("CFSv2 failed: output = '%s'", output );
endif

## now HierarchicalSearch.c without hough or followup stages
cmd = sprintf("HierarchicalSearch --followUp=0 --DataFiles1='%s' \
--Freq=%.12g --FreqBand=%.12g --f1dot=%.12g --f1dotBand=%.12g \
--method=-1 \
--nStacks1=%.12g --gridType1=3 --skyGridFile=%s --refTime=%.12g \
--printFstat1=1 --fnameout=./outHS", DataFiles, \
	      fStart, fBand, fdot, fdotBand, nStacks1, skyGridFile, \
	      refTime)
	      
[output,status] = system(cmd);
if ( status != 0 )
  error("HierarchicalSearch failed! out = '%s'", output );
endif

## compare outputs from the two codes
load "CFSv2";
load "outHS_fstatVec1.dat";

F1 = CFSv2(:,7);
F2 = outHS_fstatVec1(:,5);
diff = F1 - F2;

avg = mean(diff)
stdev = std(diff)
maximumdiff = max(diff)
minimumdiff = min(diff)

plot(diff)




