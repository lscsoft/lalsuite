fstart=40;
fend=150;
df=0.01;

frequency=fstart:df:fend;

seismicReference = 3.3e-21 * (30 ./ frequency).^14;
thermalReference = 3.8e-23 * (100 ./ frequency).^2;
shotReference = 1.13e-23 * sqrt(1 + (frequency / 90).^2);
referenceSpectrum = sqrt(shotReference.^2 + ...
                         thermalReference.^2 + ...
                         seismicReference.^2);
                     
sigmasq=4*df*sum(frequency.^(-8/3)./referenceSpectrum.^2);

sigma=sqrt(sigmasq)

rhothreshold=6;
A=rhothreshold/sigma

