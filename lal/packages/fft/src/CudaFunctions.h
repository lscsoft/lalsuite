/* supress warnings from cuda headers */
#pragma GCC system_header

#include <lal/LALDatatypes.h>
#include <cufft.h>

REAL4 *XLALCudaMallocReal(UINT4);
COMPLEX8 *XLALCudaMallocComplex(UINT4);
void XLALCudaFree(void *);


/*
 * NOTE:  the "C" and "C++" prototypes must be kept synchronized!  The
 * compiler will not check this for you!
 */


#ifdef __cplusplus
extern "C" int cudafft_execute_r2c(cufftHandle plan,
    cufftComplex *output, const cufftReal *input,
    cufftComplex *d_output, cufftReal *d_input, UINT4 size);

extern "C" int cudafft_execute_c2r(cufftHandle plan,
    cufftReal *output, const cufftComplex *input,
    cufftReal *d_output, cufftComplex *d_input, UINT4 size);

extern "C" int cudafft_execute_c2c(cufftHandle plan,
  cufftComplex *output, const cufftComplex *input,
  cufftComplex *d_output, cufftComplex *d_input,
  INT4 direction, UINT4 size);
#else
int cudafft_execute_r2c(cufftHandle plan,
    cufftComplex *output, const cufftReal *input,
    cufftComplex *d_output, cufftReal *d_input, UINT4 size);

int cudafft_execute_c2r(cufftHandle plan,
    cufftReal *output, const cufftComplex *input,
    cufftReal *d_output, cufftComplex *d_input, UINT4 size);

int cudafft_execute_c2c(cufftHandle plan,
  cufftComplex *output, const cufftComplex *input,
  cufftComplex *d_output, cufftComplex *d_input,
  INT4 direction, UINT4 size);
#endif
