function ds2 = metric3d(fa,fb,Qa,Qb,dt)

%
% metric3d is called by snr_vs_ds2 in order to calculate ds^2(ifo1,ifo2)
% based on the 3D(f,Q,t) coincidence test.
%
% Numerical inputs:
%
% fa = f in ifo1
% fb = f in ifo2
%
% Qa = Q in ifo1
% Qb = Q in ifo2
%
% dt = t_b-t_a
%
% Usage Example:
%
% ds2_H1L1 = metric3d(f_H1,f_L1,Q_H1,Q_L1,dt);
%
% Author: Fiona Speirits

f = (fa+fb)/2;
Q = (Qa+Qb)/2;
Q2 = Q^2;

gQQ= ( 1 + 28*Q2*Q2 + 128*Q2*Q2*Q2 + 64*Q2*Q2*Q2*Q2) / ( 4 * Q2 * ( 1 + 6*Q2 + 8*Q2*Q2 )*( 1 + 6*Q2 + 8*Q2*Q2 ) );

gff = ( 1 + 6*Q2 + 16*Q2*Q2) / ( 4 *f*f * ( 1 + 2*Q2 ) );
gtt = ( pi^2 * f*f ) * ( 1 + 4*Q2 ) / ( Q2 );
gQf = - ( 1 + 2*Q2 + 8*Q2*Q2 ) / ( 4*Q*f * ( 1 + 6*Q2 + 8*Q2*Q2 ) );
gtf = - ( pi * Q ) * ( 1 + 4*Q2) / ( 1 + 2*Q2 );
gtQ = ( pi * f ) * ( 1 - 2*Q2 ) / ( ( 1 + 2*Q2 )*( 1 + 2*Q2 ) );

df = fb - fa;
dQ = Qb - Qa;

ds2 = ( gQQ * dQ*dQ + gff * df*df + gtt * dt*dt + gQf * 2*dQ*df + gtf * 2*dt*df + gtQ * 2*dt*dQ );



