cd /home/siemens/revtex4/csLIGOprd/matlabfiles
fstart=40;
fend=160;%89.8750;
df=1/8;

%frequency=fstart:df:fend;

% seismicReference = 3.3e-21 * (30 ./ frequency).^14;
% thermalReference = 3.8e-23 * (100 ./ frequency).^2;
% shotReference = 1.13e-23 * sqrt(1 + (frequency / 90).^2);
% referenceSpectrum = sqrt(shotReference.^2 + ...
%                          thermalReference.^2 + ...
%                          seismicReference.^2);
                     
load lho4k_060102_strain.txt

frequency=lho4k_060102_strain(fstart/df:fend/df,1);
referenceSpectrum=lho4k_060102_strain(fstart/df:fend/df,2);

sigmasq=4*df*sum(frequency.^(-8/3)./referenceSpectrum.^2);

sigma=sqrt(sigmasq)

rhothreshold=4;
A=rhothreshold/sigma * sqrt(5)
