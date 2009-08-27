
function ds_sq = get_ds_sq (f1, f2, Q1, Q2, t1, t2)

%
% function ds_sq = get_ds_sq ( f1, f2, Q1, Q2, t1, t2 )
%
% Compute 3D metric distance in (f,Q,t) space between two ringdown templates.
% If t1 and t2 are not specified, assume dt = 0.
% 

  if exist('t1','var') && exist('t2','var')
    dt = t2 - t1;
  else
    dt = 0;
  end

  f = (f1+f2)./2;
  Q = (Q1+Q2)./2;

  gQQ = ( 1 + 28*Q.^4 + 128*Q.^6 + 64*Q.^8 ) ./ ( 4 * Q.^2 .* ( 1 + 6*Q.^2 + 8*Q.^4 ).*( 1 + 6*Q.^2 + 8*Q.^4 ) );
  gff = ( 1 + 6*Q.^2 + 16*Q.^4) ./ ( 4 * f.*f .* ( 1 + 2*Q.^2 ) );
  gtt = ( pi*pi * f.*f ) .* ( 1 + 4*Q.^2 ) ./ ( Q.^2 );
  gQf = - ( 1 + 2*Q.^2 + 8*Q.^4 ) ./ ( 4*Q.*f .* ( 1 + 6*Q.^2 + 8*Q.^4 ) );
  gtf = - ( pi * Q ) .* ( 1 + 4*Q.^2) ./ ( 1 + 2*Q.^2 );
  gtQ = ( pi * f ) .* ( 1 - 2*Q.^2 ) ./ ( ( 1 + 2*Q.^2 ).*( 1 + 2*Q.^2 ) );

  df = f2 - f1;
  dQ = Q2 - Q1;

  ds_sq = ( gQQ .* dQ.*dQ + gff .* df.*df + gtt .* dt.*dt + gQf .* 2.*dQ.*df + gtf .* 2.*dt.*df + gtQ .* 2.*dt.*dQ );

