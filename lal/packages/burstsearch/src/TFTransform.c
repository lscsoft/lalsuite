/*----------------------------------------------------------------------- 
 * 
 * File Name: TFTransform.c 
 * 
 * Author: Eanna Flanagan
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * TFTransform
 * 
 * SYNOPSIS 
 *  Routines to compute time-frequency planes from either time domain
 *  data or frequency domain data.
 * 
 * DIAGNOSTICS
 * 
 * CALLS
 *
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include "LALRCSID.h"


NRCSID (TFTRANSFORMC, "$Id$");



#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "LALStdlib.h"
#include "SeqFactories.h"
#include "RealFFT.h"
#include "ComplexFFT.h"
#include "TFTransform.h"





void
CreateRealDFTParams ( 
                     Status                         *status, 
                     RealDFTParams                  **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
		     )
{
  INITSTATUS (status, "CreateRealDFTParams", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  /* 
   * Check return structure: dftParams should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP); 
  ASSERT (*dftParams == NULL, status, TFTRANSFORM_EALLOCP, 
          TFTRANSFORM_MSGEALLOCP);


  ASSERT (params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);  

  ASSERT (params->length > 0, status, TFTRANSFORM_EPOSARG, 
          TFTRANSFORM_MSGEPOSARG);

  ASSERT( (sign==1) || (sign==-1), status, TFTRANSFORM_EINCOMP,
          TFTRANSFORM_MSGEINCOMP);

  /*  Assign memory for *dftParams   */
  *dftParams = (RealDFTParams *) LALMalloc(sizeof(RealDFTParams));
  
  /*  Make sure that the allocation was succesful */
  ASSERT (*dftParams, status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);


  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      EstimateFwdRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
    }
  else
    {
      EstimateInvRealFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
    }
  CHECKSTATUSPTR (status);

  SCreateVector (status->statusPtr, &((*dftParams)->window), params->length);
  CHECKSTATUSPTR (status);

  LALWindow (status->statusPtr, ((*dftParams)->window), params);
  CHECKSTATUSPTR (status);

  (*dftParams)->sumofsquares = params->sumofsquares;
  (*dftParams)->windowType = params->type;
  
  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
CreateComplexDFTParams ( 
                     Status                         *status, 
                     ComplexDFTParams               **dftParams, 
                     LALWindowParams                   *params,
                     INT2                           sign
		     )
{


  INITSTATUS (status, "CreateComplexDFTParams", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  /* 
   * Check return structure: dftParams should point to a valid pointer
   * which should not yet point to anything.
   *
   */
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP); 
  ASSERT (*dftParams == NULL, status, TFTRANSFORM_EALLOCP, 
          TFTRANSFORM_MSGEALLOCP);

  ASSERT (params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);  

  ASSERT (params->length > 0, status, TFTRANSFORM_EPOSARG, 
          TFTRANSFORM_MSGEPOSARG);

  ASSERT( (sign==1) || (sign==-1), status, TFTRANSFORM_EINCOMP,
          TFTRANSFORM_MSGEINCOMP);


  /*  Assign memory for *dftParams   */
  *dftParams = (ComplexDFTParams *) LALMalloc(sizeof(ComplexDFTParams));

  /*  Make sure that the allocation was succesful */
  ASSERT (*dftParams, status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);


  /* fill in some values */
  (*dftParams)->window = NULL;
  (*dftParams)->plan = NULL;

  if(sign==1)
    {
      EstimateFwdComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
    }
  else
    {
      EstimateInvComplexFFTPlan (status->statusPtr, &((*dftParams)->plan), 
                              params->length);
    }
  CHECKSTATUSPTR (status);

  SCreateVector (status->statusPtr, &((*dftParams)->window), params->length);
  CHECKSTATUSPTR (status);

  LALWindow (status->statusPtr, ((*dftParams)->window), params);
  CHECKSTATUSPTR (status);

  (*dftParams)->sumofsquares = params->sumofsquares;
  (*dftParams)->windowType = params->type;
  
  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





void
DestroyRealDFTParams (
		      Status                 *status, 
		      RealDFTParams          **dftParams
		      )
{
  INITSTATUS (status, "DestroyRealDFTParams", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (*dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that data pointed to is non-null */
  ASSERT ((*dftParams)->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP); 
  ASSERT ((*dftParams)->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP); 

  /* Ok, now let's free allocated storage */
  SDestroyVector (status->statusPtr, &((*dftParams)->window));
  CHECKSTATUSPTR (status);
  DestroyRealFFTPlan (status->statusPtr, &((*dftParams)->plan));
  CHECKSTATUSPTR (status);
  LALFree ( *dftParams );      /* free DFT parameters structure itself */

  *dftParams = NULL;	       /* make sure we don't point to freed struct */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




void
DestroyComplexDFTParams (
		         Status                 *status, 
		         ComplexDFTParams       **dftParams
		        )
{
  INITSTATUS (status, "DestroyComplexDFTParams", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not null */
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (*dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that data pointed to is non-null */
  ASSERT ((*dftParams)->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP); 
  ASSERT ((*dftParams)->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP); 

  /* Ok, now let's free allocated storage */
  SDestroyVector (status->statusPtr, &((*dftParams)->window));
  CHECKSTATUSPTR (status);
  DestroyComplexFFTPlan (status->statusPtr, &((*dftParams)->plan));
  CHECKSTATUSPTR (status);
  LALFree ( *dftParams );      /* free DFT parameters structure itself */

  *dftParams = NULL;	       /* make sure we don't point to freed struct */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}







void
ComputeFrequencySeries (
		    Status                   *status,
		    COMPLEX8FrequencySeries  *freqSeries,
		    REAL4TimeSeries          *timeSeries,
		    RealDFTParams            *dftParams
		    )
{
  REAL4Vector *tmp = NULL;
  REAL4        fac;
  INT4         i;
  INT4         n;

  INITSTATUS (status, "ComputeFrequencySeries", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (timeSeries, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (timeSeries->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (dftParams->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (dftParams->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);


  /* make sure sizes are reasonable and agree */
  n = dftParams->plan->size;
  ASSERT (n > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);
  ASSERT (timeSeries->data->length == n, status,
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  ASSERT (freqSeries->data->length == n/2 + 1, status,
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  ASSERT (dftParams->window->length == n, status,
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);


  /* copy over data into frequency series structure */
  freqSeries->epoch = timeSeries->epoch;
  freqSeries->f0    = timeSeries->f0;
  ASSERT(timeSeries->deltaT > 0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);
  freqSeries->deltaF = 1/((REAL8)(timeSeries->data->length)*timeSeries->deltaT); 
  freqSeries->name = timeSeries->name;
  freqSeries->sampleUnits = timeSeries->sampleUnits; /* fix this later */


  /* compute normalization factor */
  ASSERT(dftParams->sumofsquares > 0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);
  fac = 1 / sqrt( dftParams->sumofsquares);


  /* create temporary vector */
  CreateVector (status->statusPtr, &tmp, n);
  CHECKSTATUSPTR (status);


  /* compute windowed version of time series data */
  for (i = 0; i < n; i++)
  {
    tmp->data[i] = fac*timeSeries->data->data[i]*dftParams->window->data[i];
  };

  /* compute the DFT */
  FwdRealFFT (status->statusPtr, freqSeries->data, tmp, dftParams->plan);
  CHECKSTATUSPTR (status);

  DestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




void
CreateTFPlane (
	       Status                               *status,
	       COMPLEX8TimeFrequencyPlane           **tfp,
	       TFPlaneParams                        *input
	       )
{
  INITSTATUS (status, "CreateTFPlane", TFTRANSFORMC);

  /* Check input structure: report if NULL */
  ASSERT (input, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
      
  /* Make sure that input parameters are reasonable */
  ASSERT (input->timeBins > 0, status, 
          TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);
  ASSERT (input->freqBins > 0, status,
          TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);
  ASSERT (input->deltaT > 0.0, status,
          TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  /* 
   * Check return structure: tfp should point to a valid pointer
   * which should not yet point to anything.
   *
   */

  ASSERT (tfp != NULL, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (*tfp == NULL, status, TFTRANSFORM_EALLOCP, TFTRANSFORM_MSGEALLOCP);


  /*  Assign memory for *tfp   */
  *tfp = (COMPLEX8TimeFrequencyPlane *) LALMalloc(sizeof(COMPLEX8TimeFrequencyPlane));
  
  /*  Make sure that the allocation was succesful */
  ASSERT (*tfp, status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);

  /* assign memory for params field */
  (*tfp)->params = (TFPlaneParams *) LALMalloc(sizeof (TFPlaneParams));

  /*  Make sure that the allocation was succesful */
  ASSERT ((*tfp)->params, status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);
  


  /* 
   *  Fill some of the fields with nominal values pending the 
   *  allocation of correct values for these fields.
   */
  (*tfp)->name = NULL;
  (*tfp)->sampleUnits=NULL;
  (*tfp)->epoch.gpsSeconds=0;
  (*tfp)->epoch.gpsNanoSeconds=0;
  (*tfp)->data = NULL;    /* until allocated below */
  

  /* 
   * Allocate storage 
   */

  {
    size_t tlength;
    tlength = input->timeBins * input->freqBins * sizeof(COMPLEX8);
    (*tfp)->data = (COMPLEX8 *) LALMalloc (tlength);
  }
  if (NULL == (*tfp)->data)
  {
    /* Must free storage pointed to by *tfp */
    LALFree ((void *) *tfp);
    ABORT (status, TFTRANSFORM_EMALLOC, TFTRANSFORM_MSGEMALLOC);
    return;
  }
 

  /* 
   *  Set timeBins, freqBins etc.
   *  by copying the values from the input structure 
   */
  *((*tfp)->params) = *input;

  /* Normal exit */
  RETURN (status);
}



void
DestroyTFPlane (
	       Status                               *status,
	       COMPLEX8TimeFrequencyPlane           **tfp
	       )
{
  INITSTATUS (status, "DestroyTFPlane", TFTRANSFORMC);

  /* make sure that arguments are not null */
  ASSERT (tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (*tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that data pointed to is non-null */
  ASSERT ((*tfp)->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);

  /* make sure that pointer to parameter structure is non-null */
  ASSERT ((*tfp)->params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP); 

  /* Ok, now let's free allocated storage */
  LALFree ( (*tfp)->params );      /* free parameter structure         */
  LALFree ( (*tfp)->data );        /* free allocated data              */
  LALFree ( *tfp );	           /* free TF plane struct itself      */

  *tfp = NULL;	     	      /* make sure we don't point to freed struct */

  /* Normal exit */
  RETURN (status);
}



void
TimeSeriesToTFPlane (
            Status                               *status,
	    COMPLEX8TimeFrequencyPlane           *tfp,
	    REAL4TimeSeries                      *timeSeries,
	    VerticalTFTransformIn                *input
	    )

{
  REAL4Vector        *tmp  = NULL;
  COMPLEX8Vector     *tmp1 = NULL;
  REAL4              fac;
  INT4               i;
  INT4               j;
  INT4               nt;
  INT4               nf;
  INT4               nforig;
  INT4               ntotal;

  INT4               flow1;
  INT4               fhigh1;
  INT4               tseglength;

  INITSTATUS (status, "TimeSeriesToTFPlane", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  

  /* make sure that arguments are not NULL */
  ASSERT (timeSeries, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (timeSeries->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (timeSeries->data->data, status, TFTRANSFORM_ENULLP,
          TFTRANSFORM_MSGENULLP);

  ASSERT (input, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window->data, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);



  /* make sure that output structure is not NULL */
  ASSERT (tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);




  /*
   *
   *
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   *
   *
   */



  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT (nt > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  /* 
   * Next compute nforig = total number of bins in frequncy domain
   * before chopping  (original number of bins in frequency domain)
   * = (size of chunks of data in time domain) / 2, 
   *
   * Note tfp->params->deltaT is required sampling time of TFplane
   * while timeSeries->deltaT is sampling time of input time series.
   *
   */

  ASSERT( timeSeries->deltaT>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  
  ASSERT( tfp->params->deltaT>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  

  nforig = (INT4)( (tfp->params->deltaT) / (2.0*(timeSeries->deltaT)) );
  tseglength = 2 * nforig;

  ASSERT( nforig>0, status, TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);  
  ASSERT( tseglength == input->dftParams->plan->size, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  ASSERT( tseglength == input->dftParams->window->length, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* Supplied FFT plan must be in forward direction */
  ASSERT( input->dftParams->plan->sign==1, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  
  /* Input hetrydyne frequency must be non-negative */
  ASSERT(timeSeries->f0 >= 0.0, status, TFTRANSFORM_EPOSARG,
         TFTRANSFORM_MSGEPOSARG);

  /* sumofsquares parameter must be positive */
  ASSERT(input->dftParams->sumofsquares>0.0, status, TFTRANSFORM_EPOSARG, 
	 TFTRANSFORM_MSGEPOSARG);  

  /* compute total length of data to be used to construct TF plane */
  ntotal = 2 * nt * nforig;
  ASSERT(input->startT + ntotal <= timeSeries->data->length, status, 
         TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* 
   * Actual number of number of frequency bins to be used,
   * after bandpass filtering (chopping)
   *
   */
  nf = tfp->params->freqBins;   
  ASSERT( nf>0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  /*
   * Dealing with hetrodyned time Series (i.e. timeSeries->f0 > 0 )
   * 
   * Let f = original real frequency in time series.
   * Let fbar = rescaled frequency = f - f0.
   *
   * When we compute DFT of segment in the time domain, we get nforig
   * frequency bins separated by deltaF = 1 / tfp->params->deltaT
   * with rescaled frequencies satisfying
   *   0 \le fbar \le nforig * deltaF
   * 
   * We want to retain only the portion of the DFT with true frequencies
   * f with flow \le f \le flow + nf * deltaF, which is equivalent to 
   * retaining rescaled frequencies fbar with
   *   flow-f0 \le fbar \le flow-f0 + nf * deltaF.
   *
   * So, we just use flow - f0 in place of flow.
   *
   */

  /* low frequency cutoff in units of frequency resolution of TF plane */
  flow1 = (INT4)( (tfp->params->flow - timeSeries->f0)* tfp->params->deltaT);

  /* high frequency cutoff in units of frequency resolution of TF plane */
  fhigh1 = flow1 + nf;


  /* 
   * check is specified number of frequency bins ok
   * note that nforig = (nyquist frequency in units of frequency
   * resolution of TF plane). 
   *
   * use nforig+1 in ASSERT below rather than nforig since 
   * the DFT of a sequence of 2*n
   * real numbers is a sequence of n+1 complex numbers (with the first
   * and last being real) rather than a sequence of n complex numbers
   *
   */
  ASSERT( fhigh1 <= nforig+1, status, TFTRANSFORM_EINCOMP, 
	  TFTRANSFORM_MSGEINCOMP);  



  /* 
   * copy some of the information from the input timeseries
   * to the TFplane structure
   *
   */
  tfp->sampleUnits = timeSeries->sampleUnits; /* is this correct?
						 check later */ 


  /* 
   *  
   *  Compute the starting epoch tfp->epoch.  If we are to use data displaced
   *  from the start of the timeseries, then we need to convert
   *  from LIGOTimeGPS to REAL8 and back again to do the addition.
   *
   */
  if(input->startT)
    {
      LIGOTimeGPS t1,t2;
      REAL8 t;
      t1 = timeSeries->epoch;
      t = (REAL8)(t1.gpsSeconds) + (REAL8)(t1.gpsNanoSeconds)/1000000000.0 + (REAL8)(input->startT ) * timeSeries->deltaT;
      t2.gpsSeconds = (INT4)(t);
      t = t - (REAL8)(t2.gpsSeconds);
      t2.gpsNanoSeconds = (INT4)( t * 1000000000.0);
      tfp->epoch = t2;
    }
  else
    {
      tfp->epoch = timeSeries->epoch;
    };

  /* This TF plane is a vertical type */
  tfp->planeType = verticalPlane;


  /*
   *
   *
   *  Now start the computation of the TF transform.
   *
   *
   */

  /* create temporary vectors */
  SCreateVector (status->statusPtr, &tmp, tseglength);
  CHECKSTATUSPTR (status);
  CCreateVector (status->statusPtr, &tmp1, tseglength/2+1);
  CHECKSTATUSPTR (status);


  fac = 1/sqrt(input->dftParams->sumofsquares);

  /* loop over successive data segments in the time domain */
  for(i=0; i< nt; i++)
    {
      INT4 offset = input->startT + i * tseglength;

      /* Get the segment of time domain data and window it */
      for(j=0; j< tseglength; j++)
	{
	  tmp->data[j] = fac * timeSeries->data->data[offset+j] 
	    * input->dftParams->window->data[j];
	}
      
      /* Do the FFT */
      FwdRealFFT (status->statusPtr, tmp1, tmp, input->dftParams->plan);
      CHECKSTATUSPTR (status);
      
      /* Copy the result into appropriate spot in output structure */
      offset = i * nf;
      for(j=0; j< nf; j++)
	{
	  tfp->data[offset+j] = tmp1->data[flow1+j];
	}
       
    }

  DestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  CDestroyVector (status->statusPtr, &tmp1);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





void
FreqSeriesToTFPlane (
	       Status                               *status,
	       COMPLEX8TimeFrequencyPlane           *tfp,
	       COMPLEX8FrequencySeries              *freqSeries,
	       HorizontalTFTransformIn              *input
		    )
{
  COMPLEX8Vector     *tmp  = NULL;
  COMPLEX8Vector     *tmp1 = NULL;
  REAL4              fac;
  INT4               i;
  INT4               j;
  INT4               nt;
  INT4               nf;
  INT4               ntotal;

  INT4               flow1;
  INT4               fseglength;

  INITSTATUS (status, "FreqSeriesToTFPlane", TFTRANSFORMC);
  ATTATCHSTATUSPTR (status);

  

  /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (freqSeries->data->data, status, TFTRANSFORM_ENULLP,
          TFTRANSFORM_MSGENULLP);

  ASSERT (input, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->plan->plan, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);
  ASSERT (input->dftParams->window->data, status, TFTRANSFORM_ENULLP, 
          TFTRANSFORM_MSGENULLP);



  /* make sure that output structure is not NULL */
  ASSERT (tfp, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->params, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);
  ASSERT (tfp->data, status, TFTRANSFORM_ENULLP, TFTRANSFORM_MSGENULLP);




  /*
   *
   *
   *  make sure input parameters are reasonable, compatible with
   *  each other, etc.
   *
   *
   */



  nt = tfp->params->timeBins;   /* Number of time bins */
  ASSERT (nt > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  nf = tfp->params->freqBins;   /* Number of frequency bins */
  ASSERT (nf > 0, status, TFTRANSFORM_EPOSARG, TFTRANSFORM_MSGEPOSARG);

  /* 
   * Next compute fseglength = size of segments in freq domain
   * Round off to nearest integer rather than rounding down.
   */

  ASSERT( freqSeries->deltaF>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  
  ASSERT( tfp->params->deltaT>0.0, status, TFTRANSFORM_EPOSARG, 
         TFTRANSFORM_MSGEPOSARG);  
  fseglength = (INT4)( 0.5+ 1.0/(tfp->params->deltaT * freqSeries->deltaF));  

  /*
   * Length of segments in freq domain must exceed number of time bins
   * in final TF plane
   */
  ASSERT( fseglength >= nt, status, TFTRANSFORM_EINCOMP, 
          TFTRANSFORM_MSGEINCOMP);

  ASSERT( fseglength == input->dftParams->plan->size, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);
  ASSERT( fseglength == input->dftParams->window->length, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* Supplied FFT plan must be in inverse direction */
  ASSERT( input->dftParams->plan->sign==-1, status, 
          TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* Input lowest frequency must be non-negative */
  ASSERT(freqSeries->f0 >= 0.0, status, TFTRANSFORM_EPOSARG,
         TFTRANSFORM_MSGEPOSARG);

  /* Lowest freq of time freq plane >= lowest freq of freq series */
  ASSERT(tfp->params->flow >= freqSeries->f0, status, 
         TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* low frequency cutoff */
  flow1 = (INT4)( (tfp->params->flow - freqSeries->f0)/ freqSeries->deltaF);

  /* compute total number of data points to be used to construct TF plane */
  ntotal = nf * fseglength;

  /* make sure have enough data points in freq series */
  ASSERT(ntotal + flow1<= freqSeries->data->length, status, 
         TFTRANSFORM_EINCOMP, TFTRANSFORM_MSGEINCOMP);

  /* sumofsquares parameter must be positive */
  ASSERT( input->dftParams->sumofsquares>0.0, status, TFTRANSFORM_EPOSARG, 
          TFTRANSFORM_MSGEPOSARG);


  /* 
   * copy some of the information from the input frequency series
   * to the TFplane structure
   *
   */
  tfp->sampleUnits = freqSeries->sampleUnits; /* is this correct?
						 check later */ 
  tfp->epoch = freqSeries->epoch;

  /* the parameter input->startT is not used by this function */



  /* This TF plane is a horizontal type */
  tfp->planeType = horizontalPlane;

  /*
   *
   *
   *  Now start the computation of the TF transform.
   *
   *
   */

  /* create temporary vectors */
  CCreateVector (status->statusPtr, &tmp, fseglength);
  CHECKSTATUSPTR (status);
  CCreateVector (status->statusPtr, &tmp1, fseglength);
  CHECKSTATUSPTR (status);

  fac = 1/sqrt(input->dftParams->sumofsquares);

  /* loop over successive data segments in the freq domain */
  for(i=0; i< nf; i++)
    {
      INT4 offset = flow1 + i * fseglength;

      /* Get the segment of freq domain data and window it */
      for(j=0; j< fseglength; j++)
	{
	  tmp->data[j].re = fac * freqSeries->data->data[offset+j].re * input->dftParams->window->data[j]; 
	  tmp->data[j].im = fac * freqSeries->data->data[offset+j].im * input->dftParams->window->data[j]; 
	}

      COMPLEX8VectorFFT (status->statusPtr, tmp1, tmp, input->dftParams->plan);
      CHECKSTATUSPTR (status);
            
      /* Copy the result into appropriate spot in output structure */
      for(j=0; j< nt; j++)
	{
	  tfp->data[j*nf+i] = tmp1->data[j];
	}
       
    }

  CDestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  CDestroyVector (status->statusPtr, &tmp1);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



