function GPS = MJD_to_GPS(MJD)

% function GPS = MJD_to_GPS(MJD)
%
% This function takes in an MJD times and converts it to GPS via the
% equation given at http://www.stjarnhimlen.se/comp/time.html

Tdiff = MJD + (2400000.5-2451545.0);
meanAnomaly = 357.53 + 0.98560028*Tdiff; % mean anomaly in degrees
meanAnomaly = meanAnomaly*pi/180; % mean anomaly in rads

% time diff from TDB to TT in secs
TDBtoTT = 0.001658*sin(meanAnomaly) + 0.000014*sin(2*meanAnomaly);

GPS = (MJD-44244)*86400 - 51.184 - TDBtoTT;