#/usr/bin/octave

## start frequency of band to be analyzed 
##f0 = 140.5;
f0 = 78.5;

## length of band to be analyzed
fBand = 0.25;

## injected signal amplitude
h0Start = 9.0e-24;
steph0 = 1.0e-24;

## Threshold on significance for detection
sigmaThr = 7;
FAthr = 4;

## number of software injections
numInjections = 2;
numH0 = 2;

## get skypatch info 
load skyfileS4c

alphaVec = skyfileS4c(:,1);
deltaVec = skyfileS4c(:,2);
sizeAlphaVec = skyfileS4c(:,3);
sizeDeltaVec = skyfileS4c(:,4);


for h0Index = 1:numH0

  h0 = h0Start + (h0Index - 1)*steph0;

  ## initialize numDetections
  numDetections = 0;

  for loopIndex = 1:numInjections
    
    loopIndex
    
    ## signal params
    
    ## generate 7 random numbers in [0,1)
    randvals = rand(1,7);
    
    ## set frequency etc. 
    signalFreq = f0 + fBand * randvals(1);
    signalAlpha = 2*pi*randvals(2);
    signalDelta = pi*randvals(3) - pi/2;
    phi0 =  2*pi*randvals(4);
    psi =  2*pi*randvals(5);
    cosi =  2*randvals(6) - 1;
    signalF1dot = 0;
    
    ## frequency band to be used by makefakedata -- freq. band to be
    ## analyzed plus some wings
    mfdfmin = f0 - 0.1;
    mfdband = fBand + 0.2;
    
    ## run makefakedata
    cmdline = sprintf("lalapps_Makefakedata --outSFTbname=./ROC-SFT/ \
	--IFO=H1 --ephemDir=/home/badkri/lscsoft/share/lal/ \
      --ephemYear=05-09 --fmin=%.12g --Band=%.12g --Alpha=%.12g \
    --Delta=%.12g --h0=%.5e --cosi=%.12g --phi0=%.12g --psi=%.12g \
    --Freq=%.12g --f1dot=%.5e --noiseSFTs='/home/badkri/ROC-78/*.sft'", \ 
		      mfdfmin, mfdband, signalAlpha, signalDelta, h0, cosi, phi0, psi, \
		      signalFreq, signalF1dot);
    
    [output, status] = system(cmdline);
  
    ## now find skypatch closest to injected signal sky position
    ## this is given by the skypatch for which X1.X2 is maximum with X1
    ## being the signal position and X2 the locations of the centers of the
    ## various skypatches 
    
    ## cartesian components of signal location on unit sphere
    x = cos(signalDelta)*cos(signalAlpha);
    y = cos(signalDelta)*sin(signalAlpha);
    z = sin(signalDelta);
            
    ## cartesian components of skypatch centers
    thisX = cos(deltaVec) .* cos(alphaVec);
    thisY = cos(deltaVec) .* sin(alphaVec);
    thisZ = sin(deltaVec);
    
    ## calculate inner products and maximize
    xprodX = x.*thisX + y.*thisY + z.*thisZ;
    [maxval,maxind] = max(xprodX);
    
    ## now we know which skypatch hough driver should analyze
    skyPar(1,1) = alphaVec(maxind);
    skyPar(1,2) = deltaVec(maxind);
    skyPar(1,3) = sizeAlphaVec(maxind);
    skyPar(1,4) = sizeDeltaVec(maxind);
    
    ## write skypatch info to file which hough driver can use
    save -text tempskyfile1 skyPar
    system("grep -v '#' tempskyfile1 > tempskyfile");
    system("rm tempskyfile1");
    
    ## now run hough driver
    system("rm -rf ./outMulti/*");
    cmdline = sprintf("./DriveHoughMulti --f0=%.12g --fSearchBand=%.12g --skyfile=./tempskyfile \
    --houghThreshold=%.5e --sftDir='./ROC-SFT/*.sft' --printEvents=1 \
    --earthEphemeris=/home/badkri/lscsoft/share/lal/earth05-09.dat \
    --sunEphemeris=/home/badkri/lscsoft/share/lal/sun05-09.dat \
    --linefiles=../../LineInfo/S4lines_H1_xavi.txt.v2.thinlines,../../LineInfo/S4lines_L1_xavi.txt.v2.thinlines",f0,fBand, \
		      FAthr);
    
    [output, status] = system(cmdline);
    
    ## some cleanup
    system("rm tempskyfile");
    system("rm ROC-SFT/*");
    
    ## now analyze output to see whether signal was detected
    
    ##load output of hough driver
    load outMulti/skypatch_1/HMevents
    
    ## find templates closest in frequency
    closeFreqInd = find( abs(1800*(HMevents(:,4) - signalFreq)) < 0.5);
    
    closeAlpha = HMevents(closeFreqInd, 2);
    closeDelta = HMevents(closeFreqInd, 3);
    closeSig = HMevents(closeFreqInd, 1);
    closeFdot = HMevents(closeFreqInd, 5);
    
    ## among the templates closest in frequency, find the one closest in
    ## fdot
    
    
    ## among the templates close in freq, find the one closest in
    ## skyposition 
    closeX = cos(closeDelta) .* cos(closeAlpha);
    closeY = cos(closeDelta) .* sin(closeAlpha);
    closeZ = sin(closeDelta);
    
    xprodCloseX = x*closeX + y*closeY + z*closeZ;
    
    [dump, closestIndex] = max(xprodCloseX);
    
    ## now we have the closest template -- do we have a detection? 
    
    if closeSig(closestIndex) > sigmaThr    
      numDetections += 1;
    endif

    clear HMevents    
    
    ##freq = HMnstar(:,3);
    ##nstar = HMnstar(:,2);
    ##plot(freq,nstar);
    
    ##[maxval, maxind] = max(nstar);
    
    ## check if signal was detected
    ##if abs( freq(maxind) - signalFreq) < 2/1800
    ##  detection = 1
    ##else 
    ##  detection = 0
    ##endif
    
  endfor ## end of injections

  numDetections
  numInjections
  efficiency(h0Index) = numDetections/numInjections
  h0Vector(h0Index) = h0

endfor ## end of all h0 values

save -text ROC.out efficiency h0Vector