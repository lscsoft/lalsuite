#include <lal/LALDatatypes.h>
#include "CudaFunctions.h"

int cudafft_execute_r2c(cufftHandle plan,
    cufftComplex *output, const cufftReal *input,
    cufftComplex *d_output, cufftReal *d_input,UINT4 size)
{
    UINT4 inputBytes = size * sizeof(cufftReal);
    UINT4 outputBytes = (size/2 + 1) * sizeof(cufftComplex);

    cudaMemcpy( d_input, input, inputBytes, cudaMemcpyHostToDevice );

    cufftExecR2C(plan, d_input, d_output);

    cudaMemcpy( output, d_output, outputBytes, cudaMemcpyDeviceToHost );

    return 0;
}

int cudafft_execute_c2r(cufftHandle plan,
    cufftReal *output, const cufftComplex *input,
    cufftReal *d_output, cufftComplex *d_input, UINT4 size)
{
    UINT4 inputBytes = (size/2 + 1) * sizeof(cufftComplex);
    UINT4 outputBytes = size * sizeof(cufftReal);

    cudaMemcpy( d_input, input, inputBytes, cudaMemcpyHostToDevice );

    cufftExecC2R(plan, d_input, d_output);

    cudaMemcpy( output, d_output, outputBytes, cudaMemcpyDeviceToHost );

    return 0;
}

int cudafft_execute_c2c(cufftHandle plan,
    cufftComplex *output, const cufftComplex *input,
    cufftComplex *d_output, cufftComplex *d_input,
    INT4 direction, UINT4 size)
{
    UINT4 nBytes = size * sizeof(cufftComplex);

    cudaMemcpy( d_input, input, nBytes, cudaMemcpyHostToDevice );

    cufftExecC2C(plan, d_input, d_output, direction);

    cudaMemcpy( output, d_output, nBytes, cudaMemcpyDeviceToHost );

    return 0;
}
