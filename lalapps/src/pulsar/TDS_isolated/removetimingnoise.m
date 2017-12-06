function [hetdata, tnoise] = removetimingnoise(ephemerisFile, dataFile, ...
    hetParams, outputFile)

% function hetdata = removetimingnoise(ephemerisFile, dataFile);
%
% This function will read in an pulsar ephemeris file (in the format of my
% Crab pulsar ephemeris e.g. TOA(MJD) TOA(sec) f0(Hz) f1(1e-15 Hz^2)), and 
% read in a heterodyned data file (in the format of the output of the fine 
% heterodyne output e.g time(GPS) real[data] imag[data]). It will also take
% in a structure called hetParams, which contains the pulsar rotation 
% frequency (f0), frequency derivative (f1), second frequency derivative 
% (f2) and frequency epoch (fepoch (MJD)) used in the fine heterodyne. It
% will then work out the phase difference between the ephemeris and the
% heterodyne phase (i.e. the timing noise) and remove this from the data,
% which will be output as a structure containing the time, real and 
% imaginary parts of the data. It will also be written out to the
% outputFile.

% load the ephemeris file
ephemeris = load(ephemerisFile);

% load the heterodyned data file
data = load(dataFile);

% convert the ephemeris MJD times to GPS
times = MJD_to_GPS(ephemeris(:,1)) + ephemeris(:,2);

% convert frequency epoch to GPS
epoch = MJD_to_GPS(hetParams.fepoch);

% find the ephemeris point before the start of the data and after the end
% of the data
vals = find(times < data(1,1));
starttime = vals(end);

vals = find(times > data(end,1));
endtime = vals(1);

% calculate the phases at each time
phase(1) = 0;
phasehet(1) = 2*(hetParams.f0*(times(starttime)-epoch) + ...
    0.5*hetParams.f1*(times(starttime)-epoch)^2 + ...
    (1/6)*hetParams.f2*(times(starttime)-epoch)^3);

% update frequency values to the ephemeris time
f0 = hetParams.f0 + hetParams.f1*(times(starttime)-epoch) + ...
    0.5*hetParams.f2*(times(starttime)-epoch)^2;
f1 = hetParams.f1 + hetParams.f2*(times(starttime)-epoch);
f2 = hetParams.f2;

% assume that f and fdot are good enough give phase to nearest integer,
% which is a very good assumption over a month
for j=2:(endtime-starttime)+1
    % 2 times for gw timing noise
    phase(j) = phase(j-1) + 2*round(ephemeris(starttime+j-2,3)*...
        (times(starttime+j-1)-times(starttime+j-2)) + ...
        0.5*ephemeris(starttime+j-2,4)*1e-15*(times(starttime+j-1) ...
        -times(starttime+j-2))^2);
    
    % calculate the phase of the heterodyne parameters at each time
    phasehet(j) = phasehet(j-1) + 2*(f0*(times(starttime+j-1)...
        -times(starttime+j-2)) + 0.5*f1*(times(starttime+j-1)-...
        times(starttime+j-2))^2 + (1/6)*f2*...
        (times(starttime+j-1)-times(starttime+j-2))^3);
    
    % update frequency values to the current ephemeris time
    f0 = f0 + f1*(times(starttime+j-1)-times(starttime+j-2)) + ...
    0.5*f2*(times(starttime+j-1)-times(starttime+j-2))^2;
    f1 = f1 + f2*(times(starttime+j-1)-times(starttime+j-2));
end

% subtract the heterodyne phase from the ephemeris phase to get the timing
% noise
tnoise.times = times(starttime:endtime);
tnoise.dphase = phase - phasehet;

% get phases for each point within the heterodyned data file
dphase = interp1(tnoise.times, tnoise.dphase, data(:,1), 'spline');

plot(data(:,1), dphase);

tnoise.times = data(:,1);
%tnoise.dphase = dphase;
tnoise.dphase = 2*pi*mod(dphase, 1);
dphase = tnoise.dphase;

hetdata.times = data(:,1);
hetdata.real = data(:,2).*cos(-dphase) - data(:,3).*sin(-dphase);
hetdata.imag = data(:,2).*sin(-dphase) + data(:,3).*cos(-dphase);

fp = fopen(outputFile, 'w');

fprintf(fp, '%f\t%e\t%e\n', [hetdata.times'; hetdata.real'; hetdata.imag']);

fclose(fp);
