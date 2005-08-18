// This code writes seven numbers to stdout
// They come from uniform distributions over the respective ranges
// alpha delta and the rest
// [0, 2pi), [-1, 1), [0, 2pi) [-1,1), [-1,1), [0,1) and [0,1)

// Compile code like this:
// cc -o makeangles makeNewAngles.c -lm

#include <stdlib.h>
#include <math.h>
#include <stdio.h>


main (){
  FILE *fp;
  long int seed;
  double randval,kkcos;
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


  // get random value, print alpha [0,2pi)
  randval=drand48();
  printf("%f ", 2.0*M_PI*randval);
  
  // get random value, print delta (uniform in sky)
  randval=drand48();
  kkcos = 2.0* randval -1.0;
  printf("%f ", acos(kkcos) -0.5*M_PI);
  
    

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
