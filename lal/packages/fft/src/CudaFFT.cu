#include <lal/LALDatatypes.h>
#include "CudaFunctions.h"

void XLALCudaError(cudaError_t error, const char *file, int line)
{
    if(error != cudaSuccess)	
    {	   
        fprintf( stderr, "%s:%d %s\n", file, line, cudaGetErrorString(error));
        exit(1);
    }
}

void XLALCudaFFTError(cufftResult_t error, const char *file, int line)
{
    if(error != CUFFT_SUCCESS) 
    {
	/* As there are no GetErrorString function available for CUDA FFT, 
	 * the error messages had to be hard-coded, 
	 * and needs to be updated with new CUDA releases. 
	 */
	switch( error )
	{
	    case CUFFT_INVALID_PLAN:
	      fprintf( stderr, "%s:%d The plan handle is invalid\n", file, line );
	      break;
	
	    case CUFFT_INVALID_VALUE:
	      fprintf( stderr, "%s:%d The input data and/or output data is not valid\n", file, line );
	      break;
	      
	    case CUFFT_INTERNAL_ERROR:
	      fprintf( stderr, "%s:%d Internal driver error is detected\n", file, line );
	      break;

	    case CUFFT_EXEC_FAILED:
	      fprintf( stderr, "%s:%d CUFFT failed to execute the transform on GPU\n", file, line );
	      break;

	    case CUFFT_SETUP_FAILED:
	      fprintf( stderr, "%s:%d CUFFT library failed to initialize\n", file, line );
	      break;

	    default:
		fprintf( stderr, "%s:%d Cuda FFT Error: %d\n", file, line, error);
	}
	exit(1);
    }
}

int cudafft_execute_r2c(cufftHandle plan,
    cufftComplex *output, const cufftReal *input,
    cufftComplex *d_output, cufftReal *d_input,UINT4 size)
{
    UINT4 inputBytes = size * sizeof(cufftReal);
    UINT4 outputBytes = (size/2 + 1) * sizeof(cufftComplex);

    XLALCUDACHECK(cudaMemcpy( d_input, input, inputBytes, cudaMemcpyHostToDevice ));

    XLALCUDAFFTCHECK(cufftExecR2C(plan, d_input, d_output));

    XLALCUDACHECK(cudaMemcpy( output, d_output, outputBytes, cudaMemcpyDeviceToHost ));

    return 0;
}

int cudafft_execute_c2r(cufftHandle plan,
    cufftReal *output, const cufftComplex *input,
    cufftReal *d_output, cufftComplex *d_input, UINT4 size)
{
    UINT4 inputBytes = (size/2 + 1) * sizeof(cufftComplex);
    UINT4 outputBytes = size * sizeof(cufftReal);

    XLALCUDACHECK(cudaMemcpy( d_input, input, inputBytes, cudaMemcpyHostToDevice ));

    XLALCUDAFFTCHECK(cufftExecC2R(plan, d_input, d_output));

    XLALCUDACHECK(cudaMemcpy( output, d_output, outputBytes, cudaMemcpyDeviceToHost ));

    return 0;
}

int cudafft_execute_c2c(cufftHandle plan,
    cufftComplex *output, const cufftComplex *input,
    cufftComplex *d_output, cufftComplex *d_input,
    INT4 direction, UINT4 size)
{
    UINT4 nBytes = size * sizeof(cufftComplex);

    XLALCUDACHECK(cudaMemcpy( d_input, input, nBytes, cudaMemcpyHostToDevice ));

    XLALCUDAFFTCHECK(cufftExecC2C(plan, d_input, d_output, direction));

    XLALCUDACHECK(cudaMemcpy( output, d_output, nBytes, cudaMemcpyDeviceToHost ));

    return 0;
}
