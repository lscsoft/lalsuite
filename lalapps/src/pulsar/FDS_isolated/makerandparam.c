/*
*  Copyright (C) 2007 Iraj Gholami
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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Iraj Gholami
 * \brief
 * This code prints five + alpha random numbers to stdout,
 * where alpha specifies the number of additional random numbers to be
 * printed out, and can be specified by -n option.
 * Those random numbers are uniformaly distributed over the ranges
 * [0, 2pi) [-pi/2, pi/2) [0, 2pi) [0, 2pi) [-1, 1) and [0,1) ...
 */

/*
 *
 * gcc -static -Wall -g -o makerandparam makerandparam.c  -I/afs/aeiw/grawave/Linux/lal/include/ 
 *
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define LAL_PI        3.1415926535897932384626433832795029L  /* pi */

int main(int argc, char *argv[])
{
  int option;
  int count;
  int vrbflg = 0;
  int prmsize = 0;
  int i;
  double randval;
  unsigned int randomSeed;
  enum {krandom, user} randSeedType;

  randSeedType = krandom;

  while ( (option = getopt(argc,argv,"hvs:n:")) != -1 ) {
    switch (option) {
    case 's':
      randSeedType = user;
      randomSeed = atoi(optarg);
      break;
    case 'v':
      vrbflg = 1;
      break;
    case 'n':
      prmsize = atoi(optarg);
      break;
    case 'h':
      fprintf(stderr,"\t   This code prints 1 + N random numbers to stdout,\n");
      fprintf(stderr,"\t   The first output lies in [-pi/2, pi/2] and is uniform in sin(theta), the remaining N in [0:1]\n");
      fprintf(stderr,"Usage:\n");
      fprintf(stderr,"./makerand [-s <UINT>] [-n <UINT>] [-v] \n");
      fprintf(stderr,"-s <UINT> seed\n");
      fprintf(stderr,"-n <UINT> number of additional N random numbers in the range [0:1]\n");
      fprintf(stderr,"-v verbose output\n");
      exit(0);
    }
  }

  if ( randSeedType == krandom )
    {
      FILE   *fpRand = NULL;

      randomSeed = 0;
      /* using /dev/random is slower than /dev/urandom. */
      /*  fpRand = fopen( "/dev/random", "r" ); */
      fpRand = fopen( "/dev/urandom", "r" );
      if ( fpRand )
      {
	/* read one int */
	count=fread(&randomSeed, sizeof(randomSeed), 1, fpRand);
	if (count != 1) {
	  perror( "error: unable to read a seed from /dev/random" );
	  exit( 1 );
	}
	//	randomSeed = abs(randomSeed);
      }
      else
      {
        perror( "error: unable to open /dev/random" );
        exit( 1 );
      }
      fclose( fpRand );
    }

  srand( randomSeed );

  if ( vrbflg ) fprintf( stderr, "Random-seed used: %u\n", randomSeed );


  // get angle delta such that cos(delta) is uniform in (-1, 1)
  randval = (double)rand()/(RAND_MAX+1.0);
  randval = 2.0* randval -1;
  randval = acos(randval) - LAL_PI / 2.0;
  printf("%10.6f ", randval);

  // get N more random-numbers in the range [0,1)
  for (i=0; i < prmsize; i++) 
    {
      randval = (double) rand()/(RAND_MAX+1.0);
      printf("%10.6f ", randval);
    }

    return 0;

} /* main() */
