/* Author: Duncan Brown */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConstants.h>

int main ( int argc, char *argv[] )
{
  float m1,m2,M,eta,fmin,m,c0,c2,c3,c4,x,x2,x3,x4,x8,chirpTime,fmax;

  if ( argc != 4 )
  {
    fprintf( stderr, "usage: %s m1 m2 flow\n", argv[0] );
    exit( 1 );
  }

  m1 = atof( argv[1] );
  m2 = atof( argv[2] );
  fmin = atof( argv[3] );
  fprintf( stdout, "m1 = %f\tm2 = %f\tfLow = %f\n", m1, m2, fmin );

  M = m1 + m2;
  eta = ( m1 * m2 ) / ( M * M );
  m = 2 * ( m1 > m2 ? m2 : m1 );
  fprintf( stdout, "eta = %f\tm = %f\n", eta, M );

  fmax = 1.0 / (6.0 * sqrt(6.0) * LAL_PI * m * LAL_MTSUN_SI);
  fprintf( stdout, "isco freq = %f Hz", fmax );

  c0 = 5*M*LAL_MTSUN_SI/(256*eta);
  c2 = 743.0/252.0 + eta*11.0/3.0;
  c3 = -32*LAL_PI/3;
  c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
  x  = pow(LAL_PI*M*LAL_MTSUN_SI*fmin, 1.0/3.0);
  x2 = x*x;
  x3 = x*x2;
  x4 = x2*x2;
  x8 = x4*x4;
  chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

  fprintf( stdout, "length = %f seconds\n", chirpTime );

  return 0;
}
