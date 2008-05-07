
close all
clear all

%ifo = 'H1';
%time = '815155213'
%time2 = '822760000'

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

ifo = 'L1';
time = '824862720'
time2 = '825465258'

%ifo = 'L1';
%time = '825465258'
%time2 = '999999999'

filename = sprintf('/home/josephb/XaviCal/cluster/L1_8248_phase_outFinal.txt');
data = load(filename);

%mask = [1:12:length(data)];

%st = 1;
reHt = data(:,1);%[st:12:end],1);
imHt = data(:,2);%[st:12:end],2);
reHf = data(:,3);%[st:12:end],3);
imHf = data(:,4);%[st:12:end],4);

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
medianA = median(A);
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

A0m = medianA + 2*(rAbands(5)+i*iAbands(5))/sqrt(length(A));
A0p = medianA + 2*(rAbands(6)+i*iAbands(6))/sqrt(length(A));


bins = [0:0.01:5];
mask_last = [1:(length(bins)-1)];

figure(2)
subplot(2,2,1)
%hist(absN,bins*median(absN));
[N,B] = hist(absN/median(absN),bins);
bar(B(mask_last)*median(absN),N(mask_last))
hold on


x0 = min(B(find(N == max(N))));
xdata = B(mask_last);
ydata = N(mask_last);
options = optimset('MaxFunEvals',1000,'MaxIter',1000);
LSQ_sig = lsqcurvefit(@(x,xdata) (xdata/x^2).*exp(-(xdata).^2/(2*x^2))*sum(ydata)/sum((xdata/x^2).*exp(-(xdata).^2/(2*x^2))),x0,xdata,ydata);

dist = [B]*median(absN);
Bsig = LSQ_sig*median(absN);
rayl = raylpdf(dist,Bsig);
rayl2 = sum(N)*rayl/sum(rayl);
plot(dist,rayl2,'r')
temp = sprintf('Rayleigh distribution with sig %1.3d',Bsig);
legend('Histogram of |n|',temp)

xlabel('|N| = |Ht-Hf|')
ylabel('Number of occurances')
title(sprintf('Absolute value of N, |N| = |Ht-Hf|, IFO %s, epoch %s-%s',ifo,time,time2))
%axis([min(bins*median(absN)) max(bins*median(absN)) 0 max(N)]) 

figure(2)
subplot(2,2,3)
%hist(magHt,bins*median(magHt));
[N,B] = hist(magHt/median(magHt),bins);
bar(B(mask_last)*median(magHt),N(mask_last))
hold on


x0 = min(B(find(N == max(N))));
xdata = B(mask_last);
ydata = N(mask_last);
options = optimset('MaxFunEvals',1000,'MaxIter',1000);
LSQ_sig = lsqcurvefit(@(x,xdata) (xdata/x^2).*exp(-(xdata).^2/(2*x^2))*sum(ydata)/sum((xdata/x^2).*exp(-(xdata).^2/(2*x^2))),x0,xdata,ydata);

dist = [B]*median(magHt);
Bsig = LSQ_sig*median(magHt);
rayl = raylpdf(dist,Bsig);
rayl2 = sum(N)*rayl/sum(rayl);
plot(dist,rayl2,'r')
temp = sprintf('Rayleigh distribution with sig %1.3d',Bsig);
legend('Histogram of |Ht|',temp)

xlabel('Magnitude of Ht')
ylabel('Number of occurances')
title(sprintf('Ht, Crab band, 0.0006 bins every 0.006 Hz, IFO %s, epoch %s-%s',ifo,time,time2))
%axis([min(bins*median(absN)) max(bins*median(absN)) 0 max(N)]) 

figure(2)
subplot(2,2,4)
%hist(magHt,bins*median(magHt));
[N,B] = hist(magHf/median(magHf),bins);
bar(B(mask_last)*median(magHf),N(mask_last))
hold on


x0 = min(B(find(N == max(N))));
xdata = B(mask_last);
ydata = N(mask_last);
options = optimset('MaxFunEvals',1000,'MaxIter',1000);
LSQ_sig = lsqcurvefit(@(x,xdata) (xdata/x^2).*exp(-(xdata).^2/(2*x^2))*sum(ydata)/sum((xdata/x^2).*exp(-(xdata).^2/(2*x^2))),x0,xdata,ydata);

dist = [B]*median(magHf);
Bsig = LSQ_sig*median(magHf);
rayl = raylpdf(dist,Bsig);
rayl2 = sum(N)*rayl/sum(rayl);
plot(dist,rayl2,'r')
temp = sprintf('Rayleigh distribution with sig %1.3d',Bsig);
legend('Histogram of |Ht|',temp)

xlabel('Magnitude of Ht')
ylabel('Number of occurances')
title(sprintf('Hf, Crab band, 0.0006 bins every 0.006 Hz, IFO %s, epoch %s-%s',ifo,time,time2))
%axis([min(bins*median(absN)) max(bins*median(absN)) 0 max(N)]) 


figure(2)
subplot(2,2,2)
pdiff = mod(phaseHt - phaseHf+pi,2*pi)-pi;
hist(pdiff,100)

xlabel('Phase(Ht) - Phase(Hf), radians')
ylabel('Number of occurances')
title(sprintf('Difference in phase between Ht and Hf, IFO %s, epoch %s-%s',ifo,time,time2))
legend(sprintf('Median phase of phase(Ht) - phase(Hf) is %1.3d radians',medianPhaseA))