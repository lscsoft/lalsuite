/*
*  Copyright (C) 2010 Karsten Wiesner
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

/*-----------------------------------------------------------------------
 *
 * File Name: Chisq_GPU.cu
 *
 * Author: Wiesner, K.
 *
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <cufft.h>
#include <lal/LALAtomicDatatypes.h>

static void CudaError(cudaError_t error, const char *file, int line)
{
        if(error != cudaSuccess)
        {
                printf("%s:%d %s\n", file, line, cudaGetErrorString(error));
                exit(-1);
        }
}

#define CUDACHECK(e) (CudaError(e, __FILE__, __LINE__))

#define NUM_THREADS_MAX 512

////////////////////////////////////////////////////////////////////////////////
//! Chi Squared test kernel
//! @param g_chisq    output data in global memory
//! @param g_q        input data in global memory
//! @param g_data     time series (output of IFFT) in global memory
//! @param chisqNorm  normalization constant
////////////////////////////////////////////////////////////////////////////////
__global__ void
chisqKernel( REAL4* g_chisq, COMPLEX8* g_q, COMPLEX8 *g_data, 
	     UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm) 
{
  for (unsigned l=0; l < numChisqBins; l++)
    {
      unsigned j= blockIdx.x * blockDim.x + threadIdx.x;
      
      REAL4 Xl = crealf(g_data[l*numPoints + j]);
      REAL4 Yl = cimagf(g_data[l*numPoints + j]);
      
      REAL4 deltaXl = chisqNorm * Xl -
	(chisqNorm * crealf(g_q[j]) / (REAL4) (numChisqBins));
      REAL4 deltaYl = chisqNorm * Yl -
	(chisqNorm * cimagf(g_q[j]) / (REAL4) (numChisqBins));
      
      g_chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
    }
}

////////////////////////////////////////////////////////////////////////////////
//! Chi Squared test lalinspiral interface
////////////////////////////////////////////////////////////////////////////////
extern "C"
void Chisq_GPU (REAL4* chisq, COMPLEX8* q, COMPLEX8* qtilde, UINT4* chisqBin,
		UINT4 numPoints, UINT4 numChisqBins, REAL4 chisqNorm )
{
  
  // cudaSetDevice( 0 ) already done by cufft.

  COMPLEX8* d_q;        // input snr timeseries
  CUDACHECK( cudaMalloc( (void**) &d_q, numPoints * sizeof(COMPLEX8)) );
  //printf("Allocated device d_q\n");
  CUDACHECK( cudaMemcpy( d_q, q, numPoints * sizeof(COMPLEX8), cudaMemcpyHostToDevice) );
  //printf("Memcopy q to device d_q done\n");

  REAL4*    d_chisq;    // output chisq calculation
  CUDACHECK( cudaMalloc( (void**) &d_chisq, numPoints * sizeof(REAL4)) );
  CUDACHECK( cudaMemset( d_chisq, 0, numPoints  * sizeof(REAL4) ));
  //printf("Allocated and zero initialized device d_chisq\n");

  // d_data is the output from the batch FFT, and input
  // for the Chisq calculations
  COMPLEX8* d_data;
  CUDACHECK( cudaMalloc( (void **)&d_data, numPoints * numChisqBins * sizeof(COMPLEX8)) );
  CUDACHECK( cudaMemset( d_data, 0, numPoints * numChisqBins * sizeof(COMPLEX8) ));  
  //printf("Allocated and zero initialized device d_data\n");

  COMPLEX8 *d_qtildeBin;	
  CUDACHECK(cudaMalloc( (void **)&d_qtildeBin, numPoints * numChisqBins * sizeof(COMPLEX8)) );
  CUDACHECK( cudaMemset( d_qtildeBin, 0, numPoints * numChisqBins * sizeof(COMPLEX8) ));
  //printf("Allocated and zero initialized device q_tildeBin\n");

  //printf("\nCopy portions of qtilde to bins of zeroinitialized d_qtildeBin device memory:\n");
  //printf("d_qtildeBin (dest) points to =%p\n", d_qtildeBin);
  //printf("qtilde (src) points to = %p\n", qtilde);
  
  for( unsigned i = 0; i < numChisqBins; i++ )
    {
      //printf ("Copy bin #%d\n", i);
      //printf("dest=%p\n", (d_qtildeBin + i * numPoints + chisqBin[i]));
      //printf("src=%p   %p    0x%lx      \n", (qtilde + chisqBin[i]), qtilde, chisqBin[i] * sizeof(COMPLEX8));
      //printf("num=0x%lx\n", (chisqBin[i+1] - chisqBin[i]) * sizeof(COMPLEX8));
      
      CUDACHECK( cudaMemcpy( &d_qtildeBin[i*numPoints] + chisqBin[i],
			     (qtilde + chisqBin[i]),
			     (chisqBin[i+1] - chisqBin[i]) * sizeof(COMPLEX8),      
			     cudaMemcpyHostToDevice) );
    }
 
  //printf("\nStarting batch executions of %d cuda IFFTs -----------------------------\n", numChisqBins);
  
  cufftHandle batchPlan;
  cufftPlan1d( &batchPlan, numPoints, CUFFT_C2C, numChisqBins );
  
  cudaEvent_t start, stop;
  CUDACHECK(cudaEventCreate(&start));
  CUDACHECK(cudaEventCreate(&stop));
  CUDACHECK(cudaEventRecord(start, 0));
	
  cufftExecC2C( batchPlan, 
		(cufftComplex *)d_qtildeBin, (cufftComplex *)d_data, 
		CUFFT_INVERSE );
		
  cufftDestroy(batchPlan);	

  // GPU chisq calculation ----------------------------------------------------------
  unsigned numThreadsX = (unsigned) NUM_THREADS_MAX;
  unsigned numBlocksX = numPoints / numThreadsX;

  dim3  grid( numBlocksX, 1, 1);
  dim3  threads( numThreadsX, 1, 1);

  //printf("\nStarting grid of kernels: %d blocks with %d threads ----------------\n", numBlocksX, numThreadsX);	
  
  chisqKernel<<< grid, threads >>>( d_chisq, d_q, d_data, numPoints, numChisqBins, chisqNorm );
  
  cudaThreadSynchronize(); // implicit in cudaPrintfDisplay
  // cudaPrintfDisplay(stdout, true);
  
  // copy result from device to host ----------------------------------------
  CUDACHECK( cudaMemcpy( chisq, d_chisq, numPoints * sizeof(REAL4), cudaMemcpyDeviceToHost) );

  // cleanup device memory
  CUDACHECK(cudaFree(d_q));
  CUDACHECK(cudaFree(d_chisq));
  CUDACHECK(cudaFree(d_data));
  CUDACHECK(cudaFree(d_qtildeBin));

}
