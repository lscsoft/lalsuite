#! /usr/bin/octave -q 
##
## Copyright (C) 2006 Reinhard Prix
##
## factorial-function fact(n) taken from octave-forge
## Copyright (C) 2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

1;

## factorial(n)
##   ==prod(1:n)
function x = fact(n)
  if (any (n (:) < 0))
    error ("factorial: n be be a scalar or matrix of positive integers");
  endif
  if (isscalar (n))
    x = prod(2:n);
  else
    n (find (n < 1)) = 1;
    m = max(n(:));
    c = cumprod(1:m);
    x = c(floor(n));
  endif
endfunction

function ret = ExtrapolatePulsarSpinRange ( range0, epoch1 )
  ## extrapolate the pulsar-spin range 'range0' from epoch0 to epoch1

  dtau = epoch1 - range0.epoch;
  
  numSpins = length( range0.fkdot );
  if ( numSpins != length ( range0.fkdotBand ) )
    error ("Error: range0.fkdotBand has to have same length as range0.fkdot!\n");
  endif

  fkdotIntervals0 = [ range0.fkdot; range0.fkdot .+ range0.fkdotBand ];
  fkdotIntervals1 = zeros ( 2, numSpins);	## hold result-vector

  s = numSpins - 1;
  ks = [ 0 : s ];
  dtaukfact = (dtau .^ ks) ./ fact(ks);

  for l = 0 : s
    for k = 0 : (s - l)
      fkdotIntervals1(1, l+1) +=  min( fkdotIntervals0(:, k+l+1) * dtaukfact(k+1) );
      fkdotIntervals1(2, l+1) +=  max( fkdotIntervals0(:, k+l+1) * dtaukfact(k+1) );
    endfor
  endfor

  ret.epoch = epoch1;
  ret.fkdot = fkdotIntervals1(1,:);
  ret.fkdotBand = fkdotIntervals1(2,:) .- fkdotIntervals1(1,:);
  
endfunction

function ret = extrapolateSpins ( fkdotRef, tauRef, tauX )
  ## ret = extrapolateSpins ( fkdotRef, tauRef, tauX )
  ## propagate spins 'fkdotRef' from referenct-time 'tauRef'
  ## to new time(s) tauX, returns 'ret.Freq(tauX)' 'ret.f1dot(tauX)'
  FreqRef = fkdotRef.Freq;
  f1dotRef = f2dotRef = f3dotRef = 0;
  if ( isfield(fkdotRef, "f1dot" )) f1dotRef = fkdotRef.f1dot; endif
  if ( isfield(fkdotRef, "f2dot" )) f2dotRef = fkdotRef.f2dot; endif
  if ( isfield(fkdotRef, "f3dot" )) f3dotRef = fkdotRef.f3dot; endif
  if ( isfield(fkdotRef, "f4dot" )) error("Sorry, only 3 spindowns supported right now"); endif

  ## f(t) = f(t0 + (t-t0) ) = f(t0) + f1dot(t0) (t - t0) + 1/2 * f2dot (t-t0)^2 + ...
  dTau = tauX - tauRef;

  ret.Freq  = FreqRef + f1dotRef * dTau + (1/2) * f2dotRef * dTau.^2 + (1/6) * f3dotRef * dTau.^3;
  ret.f1dot = f1dotRef + f2dotRef * dTau + (1/2) * f3dotRef * dTau.^2;
  ret.f2dot = f2dotRef + f3dotRef * dTau;
  ret.f3dot = f3dotRef;

endfunction

## extrapolate phase phi0 @ epoch0 "backward" to phi1 @ epoch1, given fkdot1
function phi1 = extrapolatePhase ( fkdot1, epoch1, phi0, epoch0 )

  Freq1 = fkdot1.Freq;
  f1dot1 = f2dot1 = f3dot1 = 0;
  if ( isfield(fkdot1, "f1dot" )) f1dot1 = fkdot1.f1dot; endif
  if ( isfield(fkdot1, "f2dot" )) f2dot1 = fkdot1.f2dot; endif
  if ( isfield(fkdot1, "f3dot" )) f3dot1 = fkdot1.f3dot; endif
  if ( isfield(fkdot1, "f4dot" )) error("Sorry, only 3 spindowns supported right now"); endif

  dTau = epoch0 - epoch1;

  phi1 = phi0 - 2*pi* ( 
		       mod( Freq1 * dTau , 1 ) + 
		       mod( (1/2) * f1dot1 * dTau^2, 1) + 
		       mod( (1/6) * f2dot1 * dTau^3, 1) + 
		       mod( (1/24)* f3dot1 * dTau^4, 1) );

  phi1 = mod ( phi1, 2*pi);
  if ( phi1 < 0 )
    phi1 += 2 * pi;
  endif

endfunction


## ---------- extrapolate phases/spins for comparison with ExtrapolatePulsarSpinsTest ----- 

epoch0 = 714180733;
epoch1 = 714180733 + 94608000;
epoch2 = 714180733 - 94608000;

range0.epoch = epoch0;
range0.fkdot = [ 300,  -1.0000e-07, 1.0000e-15, -1.0000e-22 ];
range0.fkdotBand = [0, -1.0e-7,  1.0e-15, -1.0e-22];

save_prec = output_precision;
output_precision = 16;
rangeResult = ExtrapolatePulsarSpinRange( range0, epoch2 )

## test simple spin-extrapolation
fkdot0.Freq  =  300;
fkdot0.f1dot = -1.0000e-07;
fkdot0.f2dot =  1.0000e-15;
fkdot0.f3dot = -1.0000e-22;

fkdot1 = extrapolateSpins(fkdot0, epoch0, epoch1 )

phi0 = 1;
phi1 = extrapolatePhase ( fkdot1, epoch1, phi0, epoch0 )

output_precision = save_prec;