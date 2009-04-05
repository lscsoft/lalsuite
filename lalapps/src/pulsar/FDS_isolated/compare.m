ts = load("ts.dat");
tddSFT = load("tddSFT.dat");
Tsft = 1800;
numTS = 1;

tdd0 = load("tdd.dat.00");
tddEnd = load( sprintf("tdd.dat.0%d", numTS - 1) );
t0 = tdd0(1,1);
tEnd = tddEnd(end,1);
dt = tdd0(2,1) - tdd0(1,1);

Tspan = tEnd + dt - t0;
numTimeSamples = round(Tspan / dt);

tdd = zeros(numTimeSamples,2);
tdd(:,1) = (0:numTimeSamples-1) * dt;
for i = 1:numTS
  thisTdd = load( sprintf("tdd.dat.0%d", i - 1) );
  t0_n = thisTdd(1,1) - t0;
  i0_n = find ( abs(tdd(:,1) - t0_n) < 1e-1 );
  i1_n = i0_n + length(thisTdd) - 1;
  tdd(i0_n:i1_n, 2) = thisTdd(:,2);
endfor

%% nock out timesamples falling into gaps, if any
binMask = ones(1,numTimeSamples);	%% start with unity mask
gaps = diff ( ts(:,1) ) - Tsft;
indGaps = find ( gaps != 0 )

tddGaps = tdd;
for i = indGaps
  ii = round ( i * Tsft / dt );		%% start-index of this gap
  Ni = round(gaps(i) / dt);		%% number of gap-bins for this gap
  tddGaps ( ii:(ii+Ni), 2) = 0;		%% set this gap-bin to 0
endfor


lft0 = fft ( tddGaps(:,2) );

numBins0 = numTimeSamples /2 + 1;
lftPos0 = dt * lft0(1:numBins0);
df = 1 / Tspan;

fMax0 = (numBins0 - 1) * df;
fMin0 = 0;
freqBins0 = df * (0:numBins0-1);

lft1 = load("lft.dat");  ## computed by computeLFTfromSFTs
fmin1 = lft1(1,1);
fmax1 = lft1(end,1);

imin0 = find ( abs(freqBins0 - fmin1) < df/2)

clf;
hold on;
plot ( freqBins0(imin0:end), abs(lftPos0(imin0:end)).^2, "b+;FFT full timeseries;" );
plot ( lft1(:,1), lft1(:,2).^2 + lft1(:,3).^2, "r-x;LFTfromSFTs_Band;" );
axis([fmin1  fmax1]);
hold off;


