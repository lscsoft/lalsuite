function [residual, times, fit] = timingnoise(ephemfile, starttime, ...
    endtime, type)

% function [residual, times, fit] = timingnoise(ephemfile, starttime,...
%   endtime, type)
%
% this function reads in the ephemeris data from the ephemeris file (given
% as MJD, first TOA aon MJD (secs), spin freq, spin-down (1e-15)). It finds
% the data between starttime and endtime. It will calculate the residual in
% either 'phase' or 'frequency' given in type.

% load Crab ephemeris
ephemeris = load(ephemfile);

% find data between start and end times
vals = find(ephemeris(:,1) >= starttime & ephemeris(:,1) <= endtime);

times = (ephemeris(vals,1)-44244)*86400 - 51.184 + ephemeris(vals,2);
%times = MJD_to_GPS(ephemeris(vals,1)) + ephemeris(vals,2);
times = times - times(1); % set time epoch

if strcmp(type,'phase')
    % calculate phases at each time
    phase(1) = 0;

    % assume that f and fdot are good enough give phase to nearest integer,
    % which is a very good assumption over a month
    for j=2:length(vals)
        % 2 times for gw timing noise
        phase(j) = phase(j-1) + 2*round(ephemeris(vals(j-1),3)*(times(j) ...
            -times(j-1)) + 0.5*ephemeris(vals(j-1),4)*1e-15*(times(j) ...
            -times(j-1))^2);
    end

    % fit 3rd order polynomial to phase (f, fdot, and fdotdot)
    p = polyfit(times, phase', 3);
    %p = polyfit(times, phase', 4);
    %p = polyfit(times, phase', 5);

    % calculate phases at each time using polyfit values
    phasefit = p(4) + p(3)*times + p(2)*times.^2 + p(1)*times.^3;
    %phasefit = p(5) + p(4)*times + p(3)*times.^2 + p(2)*times.^3 + p(1)*times.^4;
    %phasefit = p(6) + p(5)*times + p(4)*times.^2 + p(3)*times.^3 + p(2)*times.^4 + p(1)*times.^5;

    fit = p;

    % calculate the residual
    residual = phase' - phasefit;

elseif strcmp(type,'frequency')
    % calculate a 2nd order fit to the frequency (i.e. include fdotdot)
    % 2 times for GW timing noise
    p = polyfit(times, 2*ephemeris(vals,3), 2);
    %p = polyfit(times, ephemeris(vals,3), 2);

    freqfit = p(3) + p(2)*times + p(1)*times.^2;
    fit = p;

    % calculate the residual
    residual = 2*ephemeris(vals,3) - freqfit;
    %residual = ephemeris(vals,3) - freqfit;
else
    error('type must be either frequency or phase');
end