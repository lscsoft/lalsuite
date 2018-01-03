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

ifo = 'H1';
time = '835044014'
time2 = '843942254'

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

filename = sprintf('/archive/home/josephb/XaviCal/S5%sCrab_1800_%s-%s/results/outTotal.txt',ifo,time,time2)
data2 = load(filename);
data(:,[1 2 3 4]) = data2(:,[1 2 3 4]); %Only a single data point, displayed 3 times for ease of reading
data(:,[5 6 7]) = data2(:,[2 3 4]);
data(:,[8,9,10]) = data2(:,[2 3 4]);

times = data(:,1);

freqdata= data(:,[4:3:end]);
freqdataHderr = data(:,[3:3:end]);
freqdataHoft = data(:,[2:3:end]);
%frequencies = load('./FreqsCrab_591.txt'); 
frequencies = [59.51 59.545 59.58];

figure(11);
freqstd = std(freqdata,0,2);
freqmean = median(freqdata);
sortedfreq = sort(freqdata, 1);
imax = floor(0.0005*length(freqdata));
rng = floor(length(freqdata)*[ 1-0.9973 1-0.9545 1-0.683 ]/2);
for ii = 1:length(rng)
    if rng(ii) > 0
    bands(2*ii-1,:) = sortedfreq(rng(ii), :);
    bands(2*ii,:) = sortedfreq(length(freqdata)-rng(ii), :);
    else
    bands(2*ii-1,:) = sortedfreq(1, :);
    bands(2*ii,:) = sortedfreq(length(freqdata)-1, :);  
    1
    end
end

set(gca, 'Layer', 'top');
a = plot(frequencies, freqmean, '.', ...
             frequencies, bands, '-');

colors = [0.7 0.7 0.7;
          0.5 0.5 0.5;
          0.3 0.3 0.3 ];
set(gca, 'Layer', 'bottom');
clear b;
for ii = 1:length(rng)
  xfill = [frequencies(1:end) frequencies(end:-1:1)]; 
  yfill = [bands(2*ii-1,:)'; bands(2*ii, end:-1:1)'];
  b(ii) = patch(xfill, yfill, colors(ii,:));
end
  
hold on;
a = semilogx(frequencies, freqmean, '.');
set(a, 'Color', [0.8 0.1 0.1], 'MarkerSize', 15);
hold off;


legend([a b], 'Median', '3 \sigma', '2 \sigma', '1 \sigma','Location','SouthEast');

xlabel('Frequency [Hz]');
ylabel('FFT(h(t)) / h(f)');
title(sprintf('%s comparison using 1800 seconds,%s-%s', ifo,time,time2));
axis([59.50 59.59 0.8 1.2]);


figure(22)
%mask = find(freqdata(:,1) ~= max(freqdata(:,1)));
%quick = freqdata(mask);
%clear freqdata
%freqdata = quick;

totalbins = floor(sqrt(length(data)));
[freqN,freqDist] = hist(freqdata(:,1),totalbins);
hist(freqdata(:,1),totalbins);
hold all;

sig = max(abs(freqmean(1) - bands(5,1)),abs(freqmean(1) - bands(6,1)));
mn = freqmean;%freqDist(find(freqN == max(freqN)));
x1 = freqDist;
y1 = gaussmf(x1,[sig mn]);
y2 = y1*sum(freqN)/sum(y1);
y3 = y1*max(freqN);
plot(x1,y2,'r')
temp = sprintf('Gaussian with mean %1.3f and std %1.3f',mn,sig);
legend('Histogram of FFT(H(t))/H(f) ratios',temp)
ylabel('Number of SFTs')
xlabel('Ratio of FFT(H(t))/H(f)')
title(sprintf('%s comparison using 1800 seconds,%s-%s', ifo,time,time2));

rmsHoft = sqrt(mean(freqdataHoft(:,1).^2))
meanHoft = mean(freqdataHoft(:,1))
rmsHderr = sqrt(mean(freqdataHderr(:,1).^2))
meanHderr = mean(freqdataHderr(:,1))
n = freqdataHoft(:,1) - freqdataHderr(:,1);
rmsn = sqrt(mean(n.^2))
meann = mean(n)
stdHoft = sqrt(rmsHoft^2 - meanHoft^2)
stdHderr = sqrt(rmsHderr^2 - meanHderr^2)
stdn = sqrt(rmsn^2 - meann^2)
finalHoft = sqrt(rmsHderr^2 + rmsn^2)
finalHderr = sqrt(rmsHoft^2 + rmsn^2)

hoftdiff = (rmsHoft - finalHoft)/rmsHoft
hderrdiff = (rmsHderr - finalHderr)/rmsHderr

figure(33)
newn = n;
for z = 1:10
    clear oldn
    oldn = newn;
    clear newn
    mask = find(abs(oldn) ~= max(abs(oldn)));
    newn = oldn(mask);
end
n2 = newn;
totalbins = floor(sqrt(length(n2)));
[freqN2,freqDist2] = hist(abs(n2),totalbins);
hist(abs(n2),totalbins);
hold all;

dist = [freqDist2];
Bsig = freqDist2(find(freqN2 == max(freqN2)));
rayl = raylpdf(dist,Bsig);
rayl2 = sum(freqN2)*rayl/sum(rayl);
plot(dist,rayl2,'r')

temp = sprintf('Rayleigh distribution normalized to the area of the histogram');
legend('Histogram of FFT(H(t))-H(f) differences',temp)
xlabel('Number of SFTs')
ylabel('Ratio of FFT(H(t))/H(f)')
title(sprintf('%s comparison using 1800 seconds,%s-%s', ifo,time,time2));

figure(44)
newn = n;
for z = 1:10
    clear oldn
    oldn = newn;
    clear newn
    mask = find(abs(oldn) ~= max(abs(oldn)));
    newn = oldn(mask);
end
n2 = newn;
totalbins = floor(sqrt(length(n2)));
[freqN2,freqDist2] = hist((n2),totalbins);
hist((n2),totalbins);
hold all;

sig = std(n2);
mn = mean(n2);%freqDist(find(freqN == max(freqN)));
x1 = freqDist2;
y1 = gaussmf(x1,[sig mn]);
y2 = y1*sum(freqN)/sum(y1);
y3 = y1*max(freqN);
plot(x1,y2,'r')
temp = sprintf('Gaussian with mean %1.3d and std %1.3d',mn,sig);
legend('Histogram of FFT(H(t)) -H(f) differences',temp)
xlabel('Difference: FFT(H(t)) - H(f)')
ylabel('Number')
title(sprintf('%s comparison using 1800 seconds,%s-%s', ifo,time,time2));

figure(66)
hist(freqdataHderr,sqrt(length(freqdataHderr)))

figure(77)
hist(freqdataHoft,sqrt(length(freqdataHoft)))


