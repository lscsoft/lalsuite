function dt = time_min(fa,fb,Qa,Qb)

%
% time_min is called by snr_vs_ds2 in order to calculate the time that
% minimises ds^2(ifo1,ifo2).
%
% Numerical inputs:
%
% fa = f in ifo1
% fb = f in ifo2
%
% Qa = Q in ifo1
% Qb = Q in ifo2
%
% Usage Example:
%
% dt = time_min(f_H1,f_L1,Q_H1,Q_L1);
%
% Author: Fiona Speirits (based on Neil Cornish's lalapps function
% XLAL3DRingTimeMinimum)

f = (fa+fb)/2;
Q = (Qa+Qb)/2;
Q2 = Q*Q;

gtt = ( pi*pi * f*f ) * ( 1 + 4*Q2 ) / ( Q2 );
gtf = - ( pi * Q ) * ( 1 + 4*Q2) / ( 1 + 2*Q2 );
gtQ = ( pi * f ) * ( 1 - 2*Q2 ) / ( ( 1 + 2*Q2 )*( 1 + 2*Q2 ) );
  
df = fb - fa;
dQ = Qb - Qa;
  
dt = -(gtf * df + gtQ * dQ)/gtt;