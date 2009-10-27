
  Spinning MCMC code
  
  
  To compile the MCMC code, you'll need	the mcmc_*.c, the header files in ./include/,
  "standard" libraries  libm, libgsl, libgslcblas and libfftw3.  You'll also need to
  get libFrame to read data files, e.g. from the ligotools package, and link it in 
  your makefile.
  
  To compile the waveform code (to produce just signal output for a certain detector), 
  you'll need the waveform_*.c, the header files in ./include/ and libfftw3.
  
  
  The library libinspiral is no longer needed, this code is included in remez.h and 
  mcmc_3rdparty.c.
  
  
  This directory (./doc/) contains an example Makefile and an example input file.
  These are supposed to be changed as the code changes, you can use them as reference
  to maintain your own files.
  
