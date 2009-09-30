#ifndef _CUDAPLAN_H
#define _CUDAPLAN_H

/* suppress warnings from cuda headers */
#pragma GCC system_header

#include <lal/LALDatatypes.h>
#include <fftw3.h>
#include <cufft.h>

/* This is only included by CudaComplexFFT.c
 * or any other file that needs to access the
 * FFT plan structure.
 */

struct tagCOMPLEX8FFTPlan
{
  INT4       sign;
  UINT4      size;
  cufftHandle plan;
  COMPLEX8  *d_input;
  COMPLEX8  *d_output;
};

struct tagCOMPLEX16FFTPlan
{
  INT4       sign;
  UINT4      size;
  fftw_plan  plan;
};

#endif	/* _CUDAPLAN_H */
