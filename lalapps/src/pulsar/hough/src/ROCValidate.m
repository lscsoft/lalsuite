#/usr/bin/octave

## produce the ROC curve for the Hough search for
## a given set of input SFTs


## start frequency of band to be analyzed 
##f0 = 140.5;
f0 = 78.5;

## length of band to be analyzed
fBand = 0.25;

## injected signal amplitude
h0Start = 1.0e-24;
steph0 = 1.0e-24;

## Threshold on significance for detection
sigmaThr = 7;
FAthr = 4;

## number of software injections
numInjections = 50;
numH0 = 10;

pixelSize = 1.0/(1.0e-4 * f0 * 2 * 1800)

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
    signalF1dot = -randvals(6)*2.2e-9;
    
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

    ## displace the template
    randvals1 = rand(1,4);
    displaceAlpha = signalAlpha + (randvals1(1) - 0.5)*pixelSize;
    displaceDelta = signalDelta + (randvals1(2) - 0.5)*pixelSize;
    displaceFreq = signalFreq + (randvals1(3) - 0.5)/1800;
    displaceF1dot = signalF1dot + (randvals1(4) -0.5)*2.2e-10;

    ## get the hough result    
    cmdline = sprintf("./ValidateHoughMulti --fStart=%.12g --fSearchBand=%.12g \ 
    --sftDir='./ROC-SFT/*.sft' --Alpha=%.12g --Delta=%.12g --Freq=%.12g \
    --fdot=%.12g --AlphaWeight=%.12g --DeltaWeight=%.12g", f0, fBand, \
		      displaceAlpha, displaceDelta, displaceFreq, \
		      displaceF1dot, alphaVec(maxind), deltaVec(maxind));

    [output, status] = system(cmdline);
        
    ## some cleanup
    system("rm ROC-SFT/*");
    
    ## now analyze output to see whether signal was detected
    
    ##load output of hough driver
    load tempout
    sigma = (tempout(1) - tempout(2))/tempout(3)

    %% save sigma value 
    saveSigma(loopIndex, h0Index) = sigma;

    if sigma > sigmaThr    
      numDetections += 1;
    endif

  endfor ## end of injections

  numDetections
  numInjections
  efficiency(h0Index) = numDetections/numInjections
  h0Vector(h0Index) = h0

endfor ## end of all h0 values

save -text ROC.out efficiency h0Vector