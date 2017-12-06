#! /usr/bin/octave -q

## mini-Monte Carlo to test different antenna-pattern functions against each other, 
## namely '0' = LALComputeAM(), '1' = LALGetAMCoeffs(), '2' = LALNewGetAMCoeffs() 

numTrials = 2000;

cmd = "./NewGetAMCoeffsTest >> amCoeffsTest.dat";

printf ("Running %d trials: ", numTrials );
for i=1:numTrials
  printf (" %d", i );
  [ out, status ] = system ( cmd );
  if ( status != 0 )
    error ("Something failed running '%s' .. \n", cmd );
  endif
endfor
printf (". done.\n");



## load results and histogram them

res = load("amCoeffsTest.dat");

hold off;
clg;
## avgerr01
hist ( res(:,4), 100 );		
title ("Average error LALComputeAM() - LALGetAMCoeffs()");
plot2eps ( "avgerr_01.eps");

## avgerr02
hist ( res(:,5), 100 );		
title ("Average error LALComputeAM() -  LALNewGetAMCoeffs()");
plot2eps ( "avgerr_02.eps");

## avgerr12
hist ( res(:,6), 100 );		
title ("Average error LALGetAMCoeffs() -  LALNewGetAMCoeffs()");
plot2eps ( "avgerr_12.eps");

## maxerr01
hist ( res(:,7), 100 );		
title ("Maximal error LALComputeAM() - LALGetAMCoeffs()");
plot2eps ( "maxerr_01.eps");

## maxerr02
hist ( res(:,8), 100 );		
title ("Maximal error LALComputeAM() -  LALNewGetAMCoeffs()");
plot2eps ( "maxerr_02.eps" );

## maxerr12
hist ( res(:,9), 100 );		
title ("Maximal error LALGetAMCoeffs() -  LALNewGetAMCoeffs()");
plot2eps ( "maxerr_12.eps" );



