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
#include <math.h>

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

__global__ void chisquared_LogLikelihood_Kernel(double *d_re, double *d_im, double deltaF, double twopit, int lower, double dataSize)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	if (idx < dataSize)
	{
		idx += lower;
		double f = ((double) idx) * deltaF;
		float cosValue;
		float sinValue;
		sincosf((float)(twopit * f), &sinValue, &cosValue);

		d_re[idx - lower] = (double)cosValue;
		d_im[idx - lower] = (double)(- sinValue);
	}
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
  REAL8 loglikeli;
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
    size_t realMemSize = dataSize * sizeof(REAL8);

    double *d_re;
    double *d_im;

    double *h_re;
    double *h_im;

    h_re = (double *)malloc(memSize);
    h_im = (double *)malloc(memSize);

    cudaMalloc((void **)&d_re, memSize);
    cudaMalloc((void **)&d_im, memSize);

    REAL8 *h_plainTemplateReal;
    h_plainTemplateReal = (REAL8 *)malloc(realMemSize);

    REAL8 *h_plainTemplateImag;
    h_plainTemplateImag = (REAL8 *)malloc(realMemSize);

    REAL8 *h_dataReal;
    h_dataReal = (REAL8 *)malloc(realMemSize);

    REAL8 *h_dataImag;
    h_dataImag = (REAL8 *)malloc(realMemSize);

    chisquared_LogLikelihood_Kernel<<<dimGrid, dimBlock>>>(d_re, d_im, deltaF, twopit, lower, dataSize);

    for (i=lower; i<=upper; ++i){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      h_plainTemplateReal[i - lower] = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      h_plainTemplateImag[i - lower] = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;
      h_dataReal[i - lower]     = dataPtr->freqData->data->data[i].re / deltaT;
      h_dataImag[i - lower]     = dataPtr->freqData->data->data[i].im / deltaT;
    }
    // block until the device has completed
    cudaThreadSynchronize();
    
    cudaMemcpy(h_re, d_re, memSize, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_im, d_im, memSize, cudaMemcpyDeviceToHost);
    
    for (i=lower; i<=upper; ++i){
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      templateReal = (h_plainTemplateReal[i - lower]*h_re[i - lower] - h_plainTemplateImag[i - lower]*h_im[i - lower]) / deltaT;
      templateImag = (h_plainTemplateReal[i - lower]*h_im[i - lower] + h_plainTemplateImag[i - lower]*h_re[i - lower]) / deltaT;

      /* compute squared difference & 'chi-squared': */
      diffRe       = h_dataReal[i - lower] - templateReal;         // Difference in real parts...
      diffIm       = h_dataImag[i - lower] - templateImag;         // ...and imaginary parts, and...
      diffSquared  = diffRe*diffRe + diffIm*diffIm;  // ...squared difference of the 2 complex figures.
      REAL8 temp = ((TwoDeltaToverN * diffSquared) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
      
      chisquared  += temp;
      dataPtr->loglikelihood -= temp;
    }
    
    cudaFree(d_re);
    cudaFree(d_im);
    free(h_re);
    free(h_im);
    free(h_plainTemplateReal);
    free(h_plainTemplateImag);
    free(h_dataReal);
    free(h_dataImag);

    dataPtr = dataPtr->next;
	//fclose(testout);
  }
  //printf ("Total Time for kernel and loop: %f ms\n", totalTime);
  loglikeli = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!
  LALInferenceDestroyVariables(&intrinsicParams);
  return(loglikeli);
}

