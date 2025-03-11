/*
*  Copyright (C) 2007 Duncan Brown
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include <stdio.h>
#include <stdlib.h>

extern void srand48 (long int);
extern double drand48 (void);

int main ( int argc, char *argv[] )
{
  int    i, ntmplts;
  double mmin, mmax;
  FILE   *fp = NULL;
  FILE   *bankfp = NULL;
  FILE   *txtfp = NULL;
  long int seed = 0;

  if ( argc != 4 )
  {
    fprintf( stderr, "useage: %s [mimumum mass] [maximum mass] [number of "
        "templates]\n", argv[0] );
    exit( 1 );
  }

  if ( (mmin = atof( argv[1] )) < 0 || (mmax = atof( argv[2] )) < 0 )
  {
    fprintf( stderr, "template masses must be positive.\n" );
    exit( 1 );
  }
  if ( (ntmplts = atoi( argv[3] )) < 1 )
  {
    fprintf( stderr, "number of templates must be greater than zero.\n" );
    exit( 1 );
  }

  /*
   * Note: /dev/random can be slow after the first few accesses, which is why we're using urandom instead.
   * [Cryptographic safety isn't a concern here at all]
   */
  if ( ! (fp = fopen( "/dev/urandom", "r" )) )
  {
    perror( "could not open /dev/urandom" );
    exit ( 1 );
  }

  for ( i = 0; i < 8; ++i )
  {
    long int rbyte = (long int) fgetc( fp );
    seed += rbyte << ( i * 8 );
  }

  fclose( fp );

  srand48( seed );

  if ( ! (bankfp = fopen( "bank.ilwd", "w" )) )
  {
    perror( "could not open bank.ilwd" );
    exit ( 1 );
  }

  if ( ! (txtfp = fopen( "bank.txt", "w" )) )
  {
    perror( "could not open bank.txt" );
    exit ( 1 );
  }

  fprintf( bankfp, "<?ilwd?>\n"
      "<ilwd name='tmpltBank::sequence' comment='seed:%ld' size='7'>\n"
      "<lstring name='real:domain' size='4'>TIME</lstring>\n"
      "<int_4u name='gps_sec:start_time' units='sec'>0</int_4u>\n"
      "<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n"
      "<int_4u name='gps_sec:stop_time' units='sec'>0</int_4u>\n"
      "<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n"
      "<real_8 name='time:step_size' units='sec'>1.0e+00</real_8>\n"
      "<real_8 dims='2,%d' name='data' ndim='2' units='mass,mass'>",
      seed, ntmplts );
  fprintf( txtfp, "#seed = %ld\n", seed );

  for ( i = 0; i < 2 * ntmplts; ++i )
  {
    double mass = mmin + ( mmax - mmin ) * drand48();
    fprintf( bankfp, "%e", mass );
    fprintf( txtfp, "%e", mass );
    if ( i % 2 )
    {
      fprintf( txtfp, "\n" );
    }
    else
    {
      fprintf( txtfp, " " );
    }
    if ( i != (2 * ntmplts - 1) )
    {
      fprintf( bankfp, " " );
    }
  }

  fprintf( bankfp, "</real_8>\n</ilwd>" );

  fclose( bankfp );
  fclose( txtfp );

  exit( 0 );

}
