ts = load("ts.dat");
tddSFT = load("tddSFT.dat");
Tsft = 1800;

tddOrig = load("tdd.dat.00");

t0 = tddOrig(1,1);
tEnd = tddOrig(end,1);
Norig = length(tddOrig);
dt0 = (tEnd - t0)/(Norig - 1);
Tspan = Norig * dt0;
Tspan0 = Tsft * round ( Tspan / Tsft );		%% make Tspan' integer multiple of Tsft

numTimeSamples0 = round(Tspan0 / dt0);
tdd0 = zeros(numTimeSamples0,2);
tdd0(:,1) = tddOrig(1:numTimeSamples0,1);	%% truncate original timeseries to Tspan'
tdd0(:,2) = tddOrig(1:numTimeSamples0,2);

%% nock out timesamples falling into gaps, if any
gaps = diff ( ts(:,1) ) - Tsft;
indGaps = find ( gaps != 0 )

tddGaps0 = tdd0;
for i = indGaps
  ii = round ( i * Tsft / dt0 );	%% start-index of this gap
  Ni = round(gaps(i) / dt0);		%% number of gap-bins for this gap
  tddGaps0 ( ii:(ii+Ni), 2) = 0;		%% set this gap-bin to 0
endfor

lftGaps0 = fft ( tddGaps0(:,2) );
lft0 = fft ( tdd0(:,2) );

numBins0 = numTimeSamples0 /2 + 1;
lftGaps0Pos = dt0 * lftGaps0(1:numBins0);
lft0Pos = dt0 * lft0(1:numBins0);

df0 = 1 / Tspan0;

fMax0 = (numBins0 - 1) * df0;
fMin0 = 0;
freqBins0 = df0 * (0:numBins0-1);

lft1 = load("lft.dat");  ## computed by computeLFTfromSFTs
fmin1 = lft1(1,1);
fmax1 = lft1(end,1);
numBins1 = length(lft1);
df1 = (fmax1 - fmin1)/(numBins1);

imin0 = find ( abs(freqBins0 - fmin1) < df0/2)

clf;
hold on;
plot ( freqBins0(imin0:end), real(lftGaps0Pos(imin0:end)), "b+;FFT TS w gaps:Re;" );
plot ( freqBins0(imin0:end), imag(lftGaps0Pos(imin0:end)), "bo;FFT TS w gaps:Im;" );

plot ( freqBins0(imin0:end), real(lft0Pos(imin0:end)), "g+-;FFT TS NO gaps:Re;" );
plot ( freqBins0(imin0:end), imag(lft0Pos(imin0:end)), "go-;FFT TS NO gaps:Im;" );

plot ( lft1(:,1), lft1(:,2), "r-x;LFTfromSFTs_Band:Re;" );
plot ( lft1(:,1), lft1(:,3), "r-*;LFTfromSFTs_Band:Im;" );
%% axis([101.12 101.13]); %% axis([fmin1  fmax1]);
axis([12.123 12.124]);
hold off;


