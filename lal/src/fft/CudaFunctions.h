/* supress warnings from cuda headers */
#pragma GCC system_header

#include <stdlib.h>
#include <stdio.h>

#include <lal/LALDatatypes.h>
#include <cufft.h>

REAL4 *XLALCudaMallocReal(UINT4);
COMPLEX8 *XLALCudaMallocComplex(UINT4);
void XLALCudaFree(void *);


#define XLALCUDACHECK(e) (XLALCudaError(e, __FILE__, __LINE__))
#define XLALCUDAFFTCHECK(e) (XLALCudaFFTError(e, __FILE__, __LINE__))

#ifdef __cplusplus
extern "C" {
#endif

void XLALCudaError(cudaError_t error, const char *file, int line);

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

#ifdef __cplusplus
}
#endif
