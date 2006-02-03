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



range0.epoch = 714180733;
range0.fkdot = [ 300,  -1.0000e-07, 1.0000e-15, -1.0000e-22, 2.0000e-29 ];
range0.fkdotBand = [0, -1.0e-7,  1.0e-15, -1.0e-22, -1.0e-30 ];

epoch1 = 714180733 + 94608000;

save_prec = output_precision;
output_precision = 16;
rangeResult = ExtrapolatePulsarSpinRange( range0, epoch1 )

output_precision = save_prec;