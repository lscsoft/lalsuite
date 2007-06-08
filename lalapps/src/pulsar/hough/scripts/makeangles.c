/*
*  Copyright (C) 2007 Alicia Sintes Olives
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

// This code writes seven numbers to stdout
// They come from uniform distributions over the respective ranges
// [0, 2pi), [-1, 1), [0, 2pi) [-1,1), [-1,1), [0,1) and [0,1)

// Compile code like this:
// cc -o makeangles makeangles.c -lm

#include <stdlib.h>
#include <math.h>
#include <stdio.h>


main (){
  FILE *fp;
  long int seed;
  double randval;
  int count;

  // open random device to get seed value
  fp=fopen("/dev/urandom", "r");

  // check that open worked
  if (!fp){
    fprintf(stderr, "Unable to open /dev/urandom\n");
    return 1;
  }

  // read one long int
  count=fread(&seed, sizeof(seed), 1, fp);
  if (count != 1){
    fprintf(stderr, "Unable to read seed from /dev/urandom\n");
    return 1;
  }
  
  // set seed
  srand48(seed);

  // get random value, print phi0
  randval=drand48();
  printf("%f ", 2.0*M_PI*randval);

  // get random value, print cos iota
  randval=drand48();
  printf("%f ", 2.0*randval-1.0);

  // get random value, print psi
  randval=drand48();
  printf("%f ", 2.0*M_PI*randval);

  // get random value, print alpha
  randval=drand48();
  printf("%f ", 2.0*randval-1.0);

  // get random value, print delta
  randval=drand48();
  printf("%f ", 2.0*randval-1.0);

  // get random value, for freq
  randval=drand48();
  printf("%f ", randval);

  // get random value, for spndown
  randval=drand48();
  printf("%f\n", randval);

  return 0;
}
