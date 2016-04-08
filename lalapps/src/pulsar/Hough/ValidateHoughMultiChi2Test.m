#! /usr/bin/octave

## Reads a list of templates and number counts produced 
## by the Hough driver (DriveHoughMultiChi2Test.c) and validates 
## the number counts by doing the Hough search for each 
## template in the list


## Set the parameters

fStart=92.1;
fBand=0.8;
DataFiles="/home/llucia/S4-SFTv2-48-110/*.sft";
##TimeStampsFile="/home/llucia/chi2/ts-inject-S4_Multi.txt";
##skyfile="./skypatchfile";

system("rm -rf ./outMultiChi2Test");
system("mkdir -p ./outMultiChi2Test");

for index = 1:1000

Freq = fStart+fBand*rand; 
fdot =-(1e-09)*rand;
Alpha = rand*2*pi;
Delta = acos(2*rand-1)-(pi/2);

AlphaWeight = Alpha;
DeltaWeight = Delta;

cmdline = sprintf("/home/llucia/CVSDIR/lalapps/src/pulsar/hough/src/ValidateChi2Test --sftDir=%s \
--fStart=%.12g --fSearchBand=%.12g  --Alpha=%.12g --Delta=%.12g \
--Freq=%.12g --fdot=%.12g --AlphaWeight=%.12g --DeltaWeight=%.12g -d0",\
		  DataFiles, fStart, fBand, Alpha, \
		  Delta, Freq, \
		  fdot, AlphaWeight, DeltaWeight); 

system ( cmdline );
cat = "cat ./tempout >> ./outMultiChi2Test/Chi2_out";
system (cat);
system("rm -f ./tempout");

endfor 


