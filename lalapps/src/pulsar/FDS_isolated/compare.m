ts = load("ts.dat");
tddSFT = load("tddSFT.dat");
Tsft = 1800;
numTS = 1;

tdd0 = load("tdd.dat.00");
tddEnd = load( sprintf("tdd.dat.0%d", numTS - 1) );
t0 = tdd0(1,1);
tEnd = tddEnd(end,1);

numBinsSFT = 3602;
dt = tdd0(2,1) - tdd0(1,1);

Tspan = tEnd + dt - t0
numTimeSamples = round(Tspan / dt)

tdd = zeros(numTimeSamples,2);
tdd(:,1) = (0:numTimeSamples-1) * dt;
for i = 1:numTS
  thisTdd = load( sprintf("tdd.dat.0%d", i - 1) );
  t0_n = thisTdd(1,1) - t0;
  i0_n = find ( abs(tdd(:,1) - t0_n) < 1e-1 );
  i1_n = i0_n + length(thisTdd) - 1;
  tdd(i0_n:i1_n, 2) = thisTdd(:,2);
endfor


# %% nock out timesamples falling into gaps, if any
# binMask = ones(1,numTimeSamples);	%% start with unity mask
# gaps = diff ( ts(:,1) ) - Tsft;
# indGaps = find ( gaps != 0 )

tddGaps = tdd;
# for i = indGaps
#   ii = round ( i * Tsft / dt );		%% start-index of this gap
#   Ni = round(gaps(i) / dt);		%% number of gap-bins for this gap
#   tddGaps ( ii:(ii+Ni), 2) = 0;		%% set this gap-bin to 0
# endfor


lft = fft ( tddGaps(:,2) );


numBins = numTimeSamples /2 + 1;

Tspan = numTimeSamples * dt;
df = 1 / Tspan;

fMax = (numBins - 1) * df;
f0 = 100.1234;
freqBins = 0:df:fMax;
freqBins += f0;


lftPos = dt * lft(1:numBins);

lft0 = load("lft.dat");

clg;
hold on;
plot ( freqBins, abs(lftPos).^2, "b+;FFT timeseries;" );
plot ( lft0(:,1), lft0(:,2).^2 + lft0(:,3).^2, "r-x;LFTfromSFTs_Band;" );
%%axis([100.12  100.13]);
hold off;


