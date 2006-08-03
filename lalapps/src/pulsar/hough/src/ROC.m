

f0 = 78.5;
fBand = 0.25;

signalFreq = 78.6;

FAthr = 1.0e-10;

alpha = 1.0;
delta = 1.0;

h0 = 1.0e-23;
phi0 = 0;
psi = 0;
cosi = 0.5;
f1dot = 0;

mfdfmin = f0 - 0.2;
mfdband = fBand + 0.4;

cmdline = sprintf("lalapps_Makefakedata --outSFTbname=./ROC-SFT/ \
    --IFO=H1 --ephemDir=/home/badkri/lscsoft/share/lal/ \
    --ephemYear=05-09 --fmin=%.12g --Band=%.12g --Alpha=%.12g \
--Delta=%.12g --h0=%.5e --cosi=%.12g --phi0=%.12g --psi=%.12g \
--Freq=%.12g --f1dot=%.5e --noiseSFTs='/home/badkri/ROC-78/*.sft'", \ 
       mfdfmin, mfdband, alpha, delta, h0, cosi, phi0, psi, \
		  signalFreq, f1dot)


[output, status] = system(cmdline);


x = sin(delta)*cos(alpha);
y = sin(delta)*sin(alpha);
z = cos(delta);

load skyfileS4c

alphaVec = skyfileS4c(:,1);
deltaVec = skyfileS4c(:,2);

sizeAlphaVec = skyfileS4c(:,3);
sizeDeltaVec = skyfileS4c(:,4);

for i = 1:length(alphaVec)

  thisAlpha = alphaVec(i);
  thisDelta = deltaVec(i);

  thisX = sin(thisDelta) * cos(thisAlpha);
  thisY = sin(thisDelta) * sin(thisAlpha);
  thisZ = cos(thisDelta);

  xprodX(i) = x*thisX + y*thisY + z*thisZ;

endfor


[maxval,maxind] = max(xprodX);

skyPar(1) = alphaVec(maxind);
skyPar(2) = deltaVec(maxind);
skyPar(3) = sizeAlphaVec(maxind);
skyPar(4) = sizeDeltaVec(maxind);

save -text tempskyfile1 skyPar

system("grep -v '#' tempskyfile1 > tempskyfile");

system("rm tempskyfile1");

system("rm -rf ./outMulti/*");

cmdline = sprintf("./DriveHoughMulti --f0=%.12g --fSearchBand=%.12g --skyfile=./tempskyfile --houghFalseAlarm=%.5e --sftDir='./ROC-SFT/*.sft' --printEvents=1 --linefiles=../../LineInfo/S4lines_H1_xavi.txt.v2.thinlines,../../LineInfo/S4lines_L1_xavi.txt.v2.thinlines",f0,fBand, FAthr)

[output, status] = system(cmdline)

##system("rm tempskyfile");

clear
