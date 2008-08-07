#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/CudaPlan.h>
#include <cufft.h>

#define	NUMTHREADSX  512

extern "C" void Cuda_CalculateChisq(UINT4 *chisqBin, COMPLEX8 *qtilde, ComplexFFTPlan *plan, COMPLEX8Vector **qBinVecPtr, UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, REAL4 *chisq, COMPLEX8 *q);

/* The kernel function, calculations that will be executed in GPU. */
__global__ void Cuda_CalculateChisq_kernel(UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, REAL4 *d_chisq, COMPLEX8 *d_data, COMPLEX8 *d_q)
{
    int i;
    UINT4 tx = threadIdx.x;
    UINT4 bx = blockIdx.x;
    UINT4 numThreadsX = blockDim.x;
    UINT4 j = tx / 2;
    UINT4 index = numThreadsX / 2;
    UINT4 startIndex, endIndex;

    __shared__ REAL4 tempChisq[NUMTHREADSX];

    switch( tx % 2 )
    {
	case 0:
	startIndex = 0;
	endIndex = numChisqBins/2;
	break;

	case 1:
	startIndex = numChisqBins/2;
	endIndex = numChisqBins;
    }

    tempChisq[tx] = 0.0;
    for( i = startIndex; i < endIndex; i++ )
    {
	REAL4 Xl = d_data[i*numPoints + bx*index + j].re;
	REAL4 Yl = d_data[i*numPoints + bx*index + j].im;
	REAL4 deltaXl = chisqNorm * Xl -
	    (chisqNorm * d_q[bx*index + j].re / (REAL4) numChisqBins);
	REAL4 deltaYl = chisqNorm * Yl - 
	    (chisqNorm * d_q[bx*index + j].im / (REAL4) numChisqBins);

	tempChisq[tx] += deltaXl * deltaXl + deltaYl * deltaYl;
    }

    __syncthreads();

    if( tx % 2 == 0 )
    {
	d_chisq[bx*index + j] = 0.0;
	for( i = tx; i < tx+2; i++ )
	{
	    d_chisq[bx*index + j] += tempChisq[i];
	}
    }
}

/* 
 * The host code for calling the GPU functions.
 * A few low level functions are used in order to access their 
 * intermediate results without the need of copying them back 
 * to the host memory. This should improve their efficiency.
 */
void Cuda_CalculateChisq(UINT4 *chisqBin, COMPLEX8 *qtilde, ComplexFFTPlan *plan, COMPLEX8Vector **qBinVecPtr, UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm, REAL4 *chisq, COMPLEX8 *q)
{
    int i;
    unsigned int numBlocksX;
    unsigned int numThreadsX;

    /* FFT plan for batch executions. */
    cufftHandle batchPlan;
    COMPLEX8 *d_qtildeBin;

    UINT4 qtildeBinBytes = sizeof(COMPLEX8) * numPoints * numChisqBins;
    UINT4 chisqBytes = sizeof(REAL4) * numPoints;
    UINT4 dataBytes = sizeof(COMPLEX8) * numPoints * numChisqBins;
    UINT4 qBytes = sizeof(COMPLEX8) * numPoints;
    COMPLEX8 *d_data, *d_q;
    REAL4 *d_chisq;

    numThreadsX = (unsigned int) NUMTHREADSX;
    /* numPoints must be divisible by numThreadsX */
    numBlocksX = numPoints * 2/ numThreadsX;

    /* Special FFT plan for batch executions. */
    cufftPlan1d( &batchPlan, numPoints, CUFFT_C2C, numChisqBins );
    cudaMalloc( (void **)&d_qtildeBin, qtildeBinBytes );

    /* 
     * d_data is the output from the batch FFT, and input 
     * for the Chisq calculations.
     */
    cudaMalloc( (void **)&d_data,  dataBytes );
    /* d_q is an input for Chisq calculations. */
    cudaMalloc( (void **)&d_q, qBytes );
    /* d_chisq is the final output for this function. */
    cudaMalloc( (void **)&d_chisq, chisqBytes );

    cudaMemcpy( d_q, q, qBytes, cudaMemcpyHostToDevice );

    cudaMemset( d_qtildeBin, 0, qtildeBinBytes );
    
    for( i = 0; i < numChisqBins; i++ )
    {
	cudaMemcpy( &d_qtildeBin[i*numPoints] + chisqBin[i], 
		qtilde + chisqBin[i], 
		(chisqBin[i+1] - chisqBin[i]) * sizeof(COMPLEX8), 
		cudaMemcpyHostToDevice );
    }
    
    /* Actual batch executions of the FFT. */
    cufftExecC2C( batchPlan, 
	    (cufftComplex *)d_qtildeBin, (cufftComplex *)d_data, 
	    (int)plan->sign );

    dim3 grid( numBlocksX, 1, 1 );
    dim3 block( numThreadsX, 1, 1 );

    cudaThreadSynchronize();
    /* Kernel function call, for executing the floating operation loops. */
    Cuda_CalculateChisq_kernel<<< grid, block >>>( numPoints, numChisqBins, chisqNorm, d_chisq, d_data, d_q);

    cudaMemcpy( chisq, d_chisq, chisqBytes, cudaMemcpyDeviceToHost );

    cufftDestroy(batchPlan);
    cudaFree(d_data);
    cudaFree(d_qtildeBin);
    cudaFree(d_q);
    cudaFree(d_chisq);
}
