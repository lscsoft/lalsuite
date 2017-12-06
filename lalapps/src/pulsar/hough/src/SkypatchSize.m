#!/usr/bin/octave

%% Set the parameters for makefakedata
%% Signal parameters
FreqSignal = 310;
fdotSignal = 0;
AlphaSignal = 0.0;
DeltaSignal = 0.0;
h0 = 5.0e-25;
cosiota = 0.5;
psiSignal = 0;
phi0Signal = 0;

%% other info
IFO="H1";
ephemDir = "/home/badkri/lscsoft/share/lal/"
ephemYear="05-09"
outputDir = "./tempsfts/"

%% sft info
fmin = FreqSignal - 1.0;
Band = 2.0;
startTime = 833752044;
duration = 3600000;
Sh = 1.0e-23;

cmdline = sprintf("rm -rf %s", outputDir);
system(cmdline);

cmdline = sprintf("mkdir -p %s", outputDir);
system(cmdline);

cmdline = sprintf("lalapps_Makefakedata --outSFTbname=%s --IFO=%s \
--ephemDir=%s  --ephemYear=%s --startTime=%.12g --duration=%.12g \
--fmin=%.12g --Band=%.12g --Alpha=%.12g --Delta=%.12g --h0=%.12g \
--cosi=%.12g --psi=%.12g --phi0=%.12g --Freq=%.12g --f1dot=%.12g \
--noiseSqrtSh=%.12g", outputDir, IFO, ephemDir, ephemYear, startTime, \
		  duration, fmin, Band, AlphaSignal, DeltaSignal, h0, \
		  cosiota, psiSignal, phi0Signal, FreqSignal, \
		  fdotSignal, Sh);


[output,stat] = system(cmdline)



%% parameters for ValidateHoughMulti
fStart = FreqSignal - 0.5;
fBand = 1.0;
AlphaSearch = AlphaSignal;
DeltaSearch = DeltaSignal;
FreqSearch = FreqSignal;
fdotSearch = fdotSignal;

%% search with weights 

mismatchWeight = 0:0.01:0.6;

for index = 1:length(mismatchWeight)

  index

  AlphaWeight = AlphaSignal + mismatchWeight(index);
  DeltaWeight = DeltaSignal + mismatchWeight(index);


  cmdline = sprintf("./ValidateHoughMulti --sftDir=%s%s \
  --fStart=%.12g --fSearchBand=%.12g  --Alpha=%.12g --Delta=%.12g \
  --Freq=%.12g --fdot=%.12g --AlphaWeight=%.12g --DeltaWeight=%.12g \
  --weighAM=1 --weighNoise=1", outputDir, "*.sft", fStart, fBand, \
		    AlphaSearch, DeltaSearch, FreqSearch, fdotSearch, \
		    AlphaWeight, DeltaWeight);
  

  [output,status] = system(cmdline);
  
  load tempout;
  
  numberCount = tempout(1)
  mean = tempout(2)
  sigma = tempout(3)
  
  sigWeight(index) = (numberCount - mean)/sigma

endfor

%% search without weights 
cmdline = sprintf("./ValidateHoughMulti --sftDir=%s%s \
--fStart=%.12g --fSearchBand=%.12g  --Alpha=%.12g --Delta=%.12g \
--Freq=%.12g --fdot=%.12g --weighAM=0 --weighNoise=0", outputDir,
		  "*.sft", fStart, fBand, \ 
		  AlphaSearch, DeltaSearch, FreqSearch, fdotSearch);


[output,status] = system(cmdline);

load tempout;

numberCount = tempout(1)
mean = tempout(2)
sigma = tempout(3)

sigNoWeight = (numberCount - mean)/sigma


for index = 1:length(mismatchWeight)
  sigNoWeightVec(index) = sigNoWeight;
endfor

plot(mismatchWeight, sigWeight,";With Weights;",  mismatchWeight, sigNoWeightVec, "3;Without Weights;" )
xlabel("Sky Position Mismatch (radians)")
ylabel("Significance")


##mismatchPosition(0.1,0,1,1)

