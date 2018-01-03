
close all
clear all

ifo = 'H1';
time = '815155213'
time2 = '822760000'

%ifo = 'H1';
%time = '822760000'
%time2 = '824695694'

%ifo = 'H1';
%time = '824695694'
%time2 = '824862720'

%ifo = 'H1';
%time = '824862720'
%time2 = '835044014'

%ifo = 'H1';
%time = '835044014'
%time2 = '843942254'

%ifo = 'H2';
%time = '815155213'
%time2 = '822760000'

%ifo = 'H2';
%time = '822760000'
%time2 = '824862720'

%ifo = 'H2';
%time = '824862720'
%time2 = '824949188'

%ifo = 'H2';
%time = '824949188'
%time2 = '846138794'

%ifo = 'L1';
%time = '816019213'
%time2 = '822760000'

%ifo = 'L1';
%time = '822760000'
%time2 = '824497827'

%ifo = 'L1';
%time = '824497827'
%time2 = '824862720'

%ifo = 'L1';
%time = '824862720'
%time2 = '825465258'

%ifo = 'L1';
%time = '825465258'
%time2 = '999999999'

filename = sprintf('/home/josephb/XaviCal/cluster/L1_825_phase_outFinal.txt',ifo,time,time2)
data = load(filename);

freq = 59.55;% load('/home/josephb/XaviCal/cluster/FreqPhase.txt');
st = 1;
%for st = 1:12
mask = [1:1:length(data)];

reHt = data([st:12:end],1);
imHt = data([st:12:end],2);
reHf = data([st:12:end],3);
imHf = data([st:12:end],4);

magHt = sqrt(reHt.^2 + imHt.^2);
phaseHt = atan2(imHt,reHt);

magHf = sqrt(reHf.^2 + imHf.^2);
phaseHf = atan2(imHf,reHf);

medianMagHt = median(magHt)
stdMagHt = std(magHt)

medianMagHf = median(magHf)
stdMagHf = std(magHf)

A = (reHt + i*imHt)./(reHf + i*imHf); %Complex transfer function between Ht and Hf
magA = abs(A);
realA = real(A);
imagA = imag(A);
phaseA = phase(A);
medianA = median(real(A))+i*median(imag(A));%median(A);%mean(A);
medianMagA = median(magA);
medianPhaseA = median(mod(phaseHt - phaseHf+pi,2*pi)-pi);

N = reHt + i*imHt - (reHf +i*imHf);
absN = abs(N);

AHf = medianA*exp(-i*median(phaseA))*(reHf + i*imHf);

sortedMagHtHf = sort(magHt./magHf);
rng = floor(length(sortedMagHtHf)*[ 1-0.9973 1-0.9545 1-0.683 ]/2);
for ii = 1:length(rng)
  if rng(ii) > 0
    bands(2*ii-1,:) = sortedMagHtHf(rng(ii), :);
    bands(2*ii,:) = sortedMagHtHf(length(sortedMagHtHf)-rng(ii), :);
  else
    bands(2*ii-1,:) = sortedMagHtHf(1,:);
    bands(2*ii,:) = sortedMagHtHf(length(sortedMagHtHf)-1, :);
  end
end

sortedPhaseHtHf = sort(mod(phaseHt - phaseHf+pi,2*pi));
for ii = 1:length(rng)
  if rng(ii) > 0
    pbands(2*ii-1,:) = sortedPhaseHtHf(rng(ii), :);
    pbands(2*ii,:) = sortedPhaseHtHf(length(sortedPhaseHtHf)-rng(ii), :);
  else
    pbands(2*ii-1,:) = sortedPhaseHtHf(1,:);
    pbands(2*ii,:) = sortedPhaseHtHf(length(sortedPhaseHtHf)-1, :);
  end
end
pbands = pbands -pi;

sortedRealA = sort(realA);
rng = floor(length(sortedRealA)*[ 1-0.9973 1-0.9545 1-0.683 ]/2);
for ii = 1:length(rng)
  if rng(ii) > 0
    rAbands(2*ii-1,:) = sortedRealA(rng(ii), :);
    rAbands(2*ii,:) = sortedRealA(length(sortedRealA)-rng(ii), :);
  else
    rAbands(2*ii-1,:) = sortedRealA(1,:);
    rAbands(2*ii,:) = sortedRealA(length(sortedRealA)-1, :);
  end
end

sortedImagA = sort(imagA);
rng = floor(length(sortedImagA)*[ 1-0.9973 1-0.9545 1-0.683 ]/2);
for ii = 1:length(rng)
  if rng(ii) > 0
    iAbands(2*ii-1,:) = sortedImagA(rng(ii), :);
    iAbands(2*ii,:) = sortedImagA(length(sortedImagA)-rng(ii), :);
  else
    iAbands(2*ii-1,:) = sortedImagA(1,:);
    iAbands(2*ii,:) = sortedImagA(length(sortedImagA)-1, :);
  end
end

meanRBands = mean(abs(rAbands([5:6])-median(real(A))));%median(real(A))));
meanIBands = mean(abs(iAbands([5:6])-median(imag(A))));%median(imag(A))));
A0m(st) = medianA - 2*(meanRBands+i*meanIBands)/sqrt(length(A));
A0p(st) = medianA + 2*(meanRBands+i*meanIBands)/sqrt(length(A));
mA(st) = medianA;
mrb(st) = meanRBands;
mib(st)= meanIBands;

mean(A)
median(A)
if st == 13
figure(st+10)
hist(real(A),100)
title('Real A');
figure(st+30)
hist(imag(A),100)
title('Imag A');
storeMean = mean(A);
storeMedian = median(A);
storeNewMedian = medianA;
storeMeanphase = mean(mod(phase(A)+pi,2*pi)-pi);
storeMedianPhase = median(mod(phase(A)+pi,2*pi)-pi);
storePhaseMean = phase(mean(A));

figure(st+50)
scatter(real(A),imag(A))
hold on
axis([-2 2 -2 2]);
scatter(real(storeMean),imag(storeMean),'r','MarkerFaceColor','r')
hold on
scatter(real(storeMedian),imag(storeMedian),'k','MarkerFaceColor','k')
legend('Complex A','Mean','Median')



end
%end
absA0m = abs(A0m)
absA0p = abs(A0p)
phaseA0m = phase(A0m)
phaseA0p = phase(A0p)
absmA = abs(mA)
phasemA = phase(mA)

lab = sprintf('%s using 1800 seconds,%s-%s',ifo,time,time2)

figure(1)
subplot(2,2,1)
plot(freq,absA0m,freq,absA0p,freq,absmA)
legend('A0 minus','A0 plus','A1')
xlabel('Freq in Hz')
ylabel('|A0|')
title(sprintf('A0 and A1 magnitude, %s',lab))

%figure(2)
subplot(2,2,2)
plot(freq,phaseA0m*180/pi,freq,phaseA0p*180/pi,freq,phasemA*180/pi)
legend('A0 minus','A0 plus','A1')
xlabel('Freq in Hz')
ylabel('phase(A0) in degrees')
title(sprintf('A0 and A1 phase, %s',lab))

%figure(3)
subplot(2,2,3)
plot(freq,real(A0m),freq,real(A0p),freq,real(mA))
legend('A0 minus','A0 plus','A1')
xlabel('Freq in Hz')
ylabel('real(A)')
title(sprintf('A0 and A1 real, %s',lab))

%figure(4)
subplot(2,2,4)
plot(freq,imag(A0m),freq,imag(A0p),freq,imag(mA))
legend('A0 minus','A0 plus','A1')
xlabel('Freq in Hz')
ylabel('imag(A)')
title(sprintf('A0 and A1 imag, %s',lab))
