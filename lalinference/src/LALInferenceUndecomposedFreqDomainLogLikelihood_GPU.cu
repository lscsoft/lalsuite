/*
*  Copyright (C) 2011 Alex Ayerdi
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
 * File Name: LALInferenceUndecomposedFreqDomainLogLikelihood_GPU.cu
 *
 * Author: Ayerdi, A.
 *
 *
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <cufft.h>
#include <cuda.h>

#include "LALInferenceUndecomposedFreqDomainLogLikelihood_GPU.h"

#define MAX_THREADS 512

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
    }                         
}

__global__ void chisquared_LogLikelihood_Kernel(REAL8 *d_sum, int lower, int dataSize,
						COMPLEX16 *freqModelhPlus_Data,
						COMPLEX16 *freqModelhCross_Data,
						COMPLEX16 *freqData_Data,
						REAL8 *oneSidedNoisePowerSpectrum_Data,
						double FplusScaled,
						double FcrossScaled,
						double deltaF,
						double twopit,
						double deltaT,
						double TwoDeltaToverN)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < dataSize)
	{
		
		__shared__ REAL8 ssum[MAX_THREADS];

		idx += lower; //accounts for the shift that was made in the original loop

		memset(ssum, 0, MAX_THREADS * sizeof(*ssum));

		int tid = threadIdx.x;
		int bid = blockIdx.x;

		REAL8 plainTemplateReal = FplusScaled * freqModelhPlus_Data[idx].re  
                          	          +  FcrossScaled * freqModelhCross_Data[idx].re;
		REAL8 plainTemplateImag = FplusScaled * freqModelhPlus_Data[idx].im
					  +  FcrossScaled * freqModelhCross_Data[idx].im;
		
                
		/* do time-shifting...             */
		/* (also un-do 1/deltaT scaling): */
		double f = ((double) idx) * deltaF;

		/* real & imag parts of  exp(-2*pi*i*f*deltaT): */
		double re = cos(twopit * f);
		double im = - sin(twopit * f);

		REAL8 templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
		REAL8 templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
		double dataReal     = freqData_Data[idx].re / deltaT;
		double dataImag     = freqData_Data[idx].im / deltaT;
		
		/* compute squared difference & 'chi-squared': */
		double diffRe       = dataReal - templateReal;         // Difference in real parts...
		double diffIm       = dataImag - templateImag;         // ...and imaginary parts, and...
				
		double diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
		
		ssum[tid] = ((TwoDeltaToverN * diffSquared) / oneSidedNoisePowerSpectrum_Data[idx]);

		/*****   REDUCTION    *****/
		
		__syncthreads(); //all the temps should have data before we add them up

		for (int i = blockDim.x / 2; i > 0; i >>= 1) { /* per block */
			if (tid < i)
			   ssum[tid] += ssum[tid + i];

			__syncthreads();
		}

		d_sum[bid] = ssum[0];
		
	}
}

__global__ void chisquared_LogLikelihood_Kernel_2(double *d_re, double *d_im, double deltaF, double twopit, int lower, int shift, double dataSize)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx < dataSize)
	{
		idx += lower;

		double f = ((double) idx) * deltaF;

		d_re[idx - lower + shift] = cos(twopit * f);
		d_im[idx - lower + shift] = - sin(twopit * f);
	}
}

__global__ void reduction_Kernel(REAL8 *d_temp, REAL8 *d_sum, int dataSize)
{
	__shared__ REAL8 ssum[MAX_THREADS];

	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int tid = threadIdx.x;
	int bid = blockIdx.x;
	
	ssum[tid] = (idx < dataSize) ? d_temp[idx] : 0;

	__syncthreads();

	for (int i = blockDim.x / 2; i > 0; i >>= 1) { /* per block */
		if (tid < i)
		   ssum[tid] += ssum[tid + i];

		__syncthreads();
	}

	if (tid == 0) d_sum[bid] = ssum[0];
}

REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood_GPU (LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
                              LALInferenceTemplateFunction *_template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  //static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit;
  double timeTmp;
  double mc;
  int different;
  UINT4 logDistFlag=0;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceVariables intrinsicParams;

  logDistFlag=LALInferenceCheckVariable(currentParams, "logdistance");
  if(LALInferenceCheckVariable(currentParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
    LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
	if(logDistFlag)
		 distMpc = exp(*(REAL8*)LALInferenceGetVariable(currentParams,"logdistance"));
	else
		 distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  //XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  //UandA.units    = MST_RAD;
  //UandA.accuracy = LALLEAPSEC_LOOSE;
  //LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
	if(logDistFlag)
			LALInferenceRemoveVariable(&intrinsicParams, "logdistance");
	else
			LALInferenceRemoveVariable(&intrinsicParams, "distance");
  // TODO: add pointer to template function here.
  // (otherwise same parameters but different template will lead to no re-computation!!)

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;

  //float totalTime;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      _template(dataPtr);
      if(XLALGetBaseErrno()==XLAL_FAILURE) /* Template generation failed in a known way, set -Inf likelihood */
          return(-DBL_MAX);

      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
//	if (!timeDomainWarning) {
//	  timeDomainWarning = 1;
//	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
//	}
        LALInferenceExecuteFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

    //FILE *testout=fopen("test_likeliLAL.txt","w");
    //fprintf(testout, "f PSD dataRe dataIm signalRe signalIm\n");
    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    // printf("deltaF %g, Nt %d, deltaT %g\n", deltaF, dataPtr->timeData->data->length, dataPtr->timeData->deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);

    /*****    CUDA  SUBSTITUTION  HERE    *****/

    int dataSize = upper-lower + 1;

    int numThreadsPerBlock = 0;

    int remainder[9]; //hold the remainders of the modulus operations
    int j = 0;
    //calculate the number of threads needed for this particular dataSize
    for (i = MAX_THREADS; i > 0; i >>= 1)
    {
        if (dataSize % i == 0) //split up the number of threads evenly across the dataset
        {
            numThreadsPerBlock = i;  
            break;
        }

        remainder[j] = dataSize % i; //store the remainder
        j++;

    }

    int numBlocks = dataSize / numThreadsPerBlock;

    //recalculate a better number of threads per block if the number of blocks is exceeding hardware limits
    if (numBlocks > 65000 || (numBlocks > 65000 && numThreadsPerBlock == 1))
    {
        int tMin = remainder[0]; //initialize the minimum remainder
        int tracker = 0; //tracker for which index had the minimum remainder
        for (j = 1; j < 9; j++)
        {
            if (remainder[j] < tMin && remainder[j] > 0)
            {
                tMin = remainder[j];
                tracker = j;
            }
        }
	
        switch (tracker)
        {
	    //choose the number of threads per block based on the least remainder
            case 0: numThreadsPerBlock = 512; break;
            case 1: numThreadsPerBlock = 256; break;
            case 2: numThreadsPerBlock = 128; break;
            case 3: numThreadsPerBlock = 64; break;
            case 4: numThreadsPerBlock = 32; break;
            case 5: numThreadsPerBlock = 16; break;
            case 6: numThreadsPerBlock = 8; break;
            case 7: numThreadsPerBlock = 4; break;
            default: numThreadsPerBlock = 2; break;
        }

        numBlocks = dataSize / numThreadsPerBlock + 1; //we need to give it an extra so it accounts for the data not evenly distributed
    }

    dim3 dimGrid(numBlocks); 
    dim3 dimBlock(numThreadsPerBlock);

    size_t memSize = dataSize * sizeof(double);

    //for concurrent streams
    //cudaMallocHost of host data
    //cudaMalloc of device data
    //make two streams cudaStream_t stream1 cudaStream_t stream2
    //create the streams cudaStreamCreate(&stream1)  cudaStreamCreate(&stream2)
    //create event handlers  cudaEvent_t start_event, stop_event;
    //create the event cudaEventCreate(&start_event) cudaEventCreate(&stop_event);
    //start recording the event cudaEventRecord(start_event, 0);
    //launch kernel to streams kernel<<<dimGrid,dimBlocks,0,stream1>>>(args);  kernel<<<dimGrid,dimBlocks,0,stream1>>>(args);
    //optional: after each launch you need an event cudaEventRecord(kernelEvent[i], streams[i])
    //copy the data cudaMemcpyAsync(host, device, sizeof(memorySize), cudaMemcpyDeviceToHost, stream1/2)
    //record the event cudaEventRecord(stop_event, 0)
    //cudaEventSynchronize(stop_event)
    //cudaStreamDestroy(stream1/2);
    //cudaEventDestroy(start_event);
    //cudaEventDestroy(stop_event);
    //cudaFreeHost(host);
    //cudaFree(device);

    double *d_re;
    double *d_im;

    double *h_re;
    double *h_im;

    //h_re = (double *)malloc(memSize);
    //h_im = (double *)malloc(memSize);

    cudaMallocHost((void **)&h_re, memSize);
    cudaMallocHost((void **)&h_im, memSize);

    cudaMalloc((void **)&d_re, memSize);
    cudaMalloc((void **)&d_im, memSize);

    cudaStream_t stream1;
    cudaStream_t stream2;
    cudaStreamCreate(&stream1);
    cudaStreamCreate(&stream2);

    cudaEvent_t start_event, stop_event;
    cudaEventCreate(&start_event);
    cudaEventCreate(&stop_event);

    cudaEventRecord(start_event, 0);

    chisquared_LogLikelihood_Kernel_2<<<dimGrid, dimBlock, 0, stream1>>>(d_re, d_im, deltaF, twopit, lower, 0, dataSize / 2);

    chisquared_LogLikelihood_Kernel_2<<<dimGrid, dimBlock, 0, stream2>>>(d_re, d_im, deltaF, twopit, lower + (dataSize / 2), dataSize / 2, dataSize / 2);

    // block until the device has completed
    //cudaThreadSynchronize();

    //checkCUDAError("kernel execution");

    cudaMemcpyAsync(h_re, d_re, memSize, cudaMemcpyDeviceToHost, stream1);
    cudaMemcpyAsync(h_im, d_im, memSize, cudaMemcpyDeviceToHost, stream2);

    cudaEventRecord(stop_event, 0);
    cudaEventSynchronize(stop_event);

    cudaStreamDestroy(stream1);
    cudaStreamDestroy(stream2);
    cudaEventDestroy(start_event);
    cudaEventDestroy(stop_event);

    //checkCUDAError("cudaMemcpy");

    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      templateReal = (plainTemplateReal*h_re[i - lower] - plainTemplateImag*h_im[i - lower]) / deltaT;
      templateImag = (plainTemplateReal*h_im[i - lower] + plainTemplateImag*h_re[i - lower]) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         // Difference in real parts...
      diffIm       = dataImag - templateImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  // ...squared difference of the 2 complex figures.
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      
      chisquared  += temp;
      dataPtr->loglikelihood -= temp;
    }
    
    cudaFree(d_re);
    cudaFree(d_im);
    cudaFreeHost(h_re);
    cudaFreeHost(h_im);

    //free(h_re);
    //free(h_im);

    //size_t sumMemSize = numBlocks * sizeof(REAL8);
    //size_t complexMemSize = (upper + 1) * sizeof(COMPLEX16); 
    //size_t realMemSize = (upper + 1) * sizeof(REAL8);

    //the host pointer to be used for the overall summation in host
    //REAL8 *h_sum;   

    //h_sum = (REAL8 *)malloc(sumMemSize);
    //memset(h_sum, 0, sumMemSize);

    //the device pointer for the per-block reduction summation
    //REAL8 *d_sum;   

    //device pointers for data from dataPtr
    //COMPLEX16 *d_freqModelhPlus_data;
    //COMPLEX16 *d_freqModelhCross_data;
    //COMPLEX16 *d_freqData_data;
    //REAL8 *d_oneSidedNoisePowerSpectrum_data;

    //device memory allocation
    //cudaMalloc((void **)&d_sum, sumMemSize);
   
    //cudaMalloc((void **)&d_freqModelhPlus_data, complexMemSize);
    //cudaMalloc((void **)&d_freqModelhCross_data, complexMemSize);
    //cudaMalloc((void **)&d_freqData_data, complexMemSize);
    //cudaMalloc((void **)&d_oneSidedNoisePowerSpectrum_data, realMemSize);

    //cudaMemcpy(d_freqModelhPlus_data, dataPtr->freqModelhPlus->data->data, complexMemSize, cudaMemcpyHostToDevice);
    //cudaMemcpy(d_freqModelhCross_data, dataPtr->freqModelhCross->data->data, complexMemSize, cudaMemcpyHostToDevice);
    //cudaMemcpy(d_freqData_data, dataPtr->freqData->data->data, complexMemSize, cudaMemcpyHostToDevice);
    //cudaMemcpy(d_oneSidedNoisePowerSpectrum_data, dataPtr->oneSidedNoisePowerSpectrum->data->data, realMemSize, cudaMemcpyHostToDevice);
    
    //cudaEventRecord(start, 0);

    //launch the kernel
    //chisquared_LogLikelihood_Kernel<<<dimGrid, dimBlock>>>(d_sum, lower, dataSize,
							//d_freqModelhPlus_data,
							//d_freqModelhCross_data,
							//d_freqData_data,
							//d_oneSidedNoisePowerSpectrum_data,
	    				   		//FplusScaled,
            				   		//FcrossScaled,
            				   		//deltaF,
            				   		//twopit,
            				   		//deltaT,
            				   		//TwoDeltaToverN);

    //chisquared_LogLikelihood_Kernel<<<dimGrid, dimBlock>>>(d_re, d_im);

    // block until the device has completed
    //cudaThreadSynchronize();

    // check if kernel execution generated an error
    //checkCUDAError("kernel execution");

    //cudaEventRecord(stop, 0);
    //cudaEventSynchronize(stop);

    //cudaMemcpy(h_sum, d_sum, sumMemSize, cudaMemcpyDeviceToHost);

    //checkCUDAError("cudaMemcpy");

    //for (i = 0; i < numBlocks; i++)
    //{
    //    chisquared += h_sum[i];
    //    dataPtr->loglikelihood -= h_sum[i];
    //}
    
    //cudaFree(d_sum);
    //cudaFree(d_freqModelhPlus_data);
    //cudaFree(d_freqModelhCross_data);
    //cudaFree(d_freqData_data);
    //cudaFree(d_oneSidedNoisePowerSpectrum_data);
    //free(h_sum);

    //cudaEventElapsedTime(&time, start, stop);
    //printf ("Time for the kernel: %f ms\n", time);
    //totalTime += time;

    dataPtr = dataPtr->next;
	//fclose(testout);
  }
  //printf ("Total Time for kernel and loop: %f ms\n", totalTime);
  loglikeli = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!
  LALInferenceDestroyVariables(&intrinsicParams);
  return(loglikeli);
}

