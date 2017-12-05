#! /usr/bin/octave -q

%%
%% Copyright (C) 2006 Reinhard Prix
%%
%%  This program is free software; you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation; either version 2 of the License, or
%%  (at your option) any later version.
%%
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%
%%  You should have received a copy of the GNU General Public License
%%  along with with program; see the file COPYING. If not, write to the
%%  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%%  MA  02111-1307  USA
%%

%% script to generate ephemeris-file for LISA to be read by LALIntitBarycenter()
%% for now, simply implement Eq.(2.1) in MLDC Challenge1 Document  [Draft v1.0]
%% which is the one used in MLDC1
%% NOTE however that we use the GUIDING-CENTER as ephemeris, i.e. the center-of-mass from Eq.(2.1),
%% which simplifies the expression quite a bit...

%% some useful constants
YRSID_SI	= 31558149.8;		%% Sidereal year (1994), s
AU_SI		= 1.4959787066e11;	%% Astronomical unit, m
C_SI 		= 299792458; 					%% Speed of light in vacuo, m s^-1

## MLDC time-origin "t=0" (MUST be equal to setting in lal/packages/pulsar/include/LISAspecifics.h)
LISA_TIME_ORIGIN = 700000000;

if ( nargin != 2 )
  error ("\nNeed exactly two input arguments: start-GPS and end-GPS!");
endif

startGPS = str2num ( argv{1} );
endGPS   = str2num ( argv{2} );

duration = endGPS - startGPS;

%% fields to write in the first line of ephem-file
gpsYr = ceil(startGPS); %% 'gpsYr' entry, not really used anywhere so we simply set it to start-time
tStep = 14400;	%% hardcoded default for now
nSteps = ceil( duration / tStep - 1e-6) + 1;

lastGPS = startGPS + (nSteps-1) * tStep;

%% timesteps to compute LISA ephemeris for
tiGPS = [startGPS:tStep:lastGPS]' ;

%% ----- implement LISA orbit Eq.(2.1) -----
%% default at 't=0'
kappa = 0;

%% orbital frequency
OmOrb = 2 * pi / YRSID_SI;

alpha = OmOrb * (tiGPS - LISA_TIME_ORIGIN) + kappa;

%% center-of-mass motion of Eq.(2.1):
rOrb = AU_SI / C_SI;	%% ephemeris-files are in units of r/c

%% orbital positions
xG = rOrb * cos ( alpha );
yG = rOrb * sin ( alpha );
zG = 0 * xG;
%% velocities
x1G = - OmOrb * yG;
y1G =   OmOrb * xG;
z1G =   0 * xG;
%% accelerations
x2G = - OmOrb^2 * xG;
y2G = - OmOrb^2 * yG;
z2G = 0 * xG;

fname = "ephemMLDC.dat";
[fid, msg] = fopen (fname, "wb");
if ( fid == -1 )
  error("Failed to open ephemeris-file '%s'  for writing: error = %s", fname, msg );
endif

%% write header-line to file:
fprintf (fid, "%d   %f   %d\n",  gpsYr, tStep, nSteps );

%% write the time-series of [tGPS, pos, vel, acc ]:
ephem = [ tiGPS,  xG, yG, zG,  x1G, y1G, z1G,   x2G, y2G, z2G ];

fprintf (fid, "%.16g    %.16g  %.16g  %.16g    %.16g  %.16g  %.16g    %.16g %.16g %.16g\n", ephem' );

fclose ( fid );
