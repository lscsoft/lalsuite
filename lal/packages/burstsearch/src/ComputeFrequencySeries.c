/******** <lalVerbatim file="ComputeFrequencySeriesCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (COMPUTEFREQUENCYSERIESC, "$Id$");


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="ComputeFrequencySeriesCP"> ********/
void
LALComputeFrequencySeries (
		    LALStatus                   *status,
		    COMPLEX8FrequencySeries  *freqSeries,
		    REAL4TimeSeries          *timeSeries,
		    RealDFTParams            *dftParams
		    )
/******** </lalVerbatim> ********/
{
  REAL4Vector *tmp = NULL;
  REAL4        fac;
  INT4         i;
  INT4         n;

  INITSTATUS (status, "LALComputeFrequencySeries", COMPUTEFREQUENCYSERIESC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (freqSeries->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (timeSeries, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (timeSeries->data, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (dftParams, status, TFTRANSFORMH_ENULLP, TFTRANSFORMH_MSGENULLP);
  ASSERT (dftParams->plan, status, TFTRANSFORMH_ENULLP, 
          TFTRANSFORMH_MSGENULLP);
  ASSERT (dftParams->window, status, TFTRANSFORMH_ENULLP, 
          TFTRANSFORMH_MSGENULLP);


  /* make sure sizes are reasonable and agree */
  /* MODIFIED -- JC
   * n = dftParams->plan->size;
   */
  n = timeSeries->data->length;
  ASSERT (n > 0, status, TFTRANSFORMH_EPOSARG, TFTRANSFORMH_MSGEPOSARG);
  /*
  ASSERT ((INT4)timeSeries->data->length == n, status,
          TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);
  */
  ASSERT ((INT4)freqSeries->data->length == n/2 + 1, status,
          TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);
  ASSERT ((INT4)dftParams->window->length == n, status,
          TFTRANSFORMH_EINCOMP, TFTRANSFORMH_MSGEINCOMP);


  /* copy over data into frequency series structure */
  freqSeries->epoch = timeSeries->epoch;
  freqSeries->f0    = timeSeries->f0;
  ASSERT(timeSeries->deltaT > 0.0, status, TFTRANSFORMH_EPOSARG, 
         TFTRANSFORMH_MSGEPOSARG);
  freqSeries->deltaF = 1/((REAL8)(timeSeries->data->length)*timeSeries->deltaT); 
  /* 
   * OMITTED
   *
  freqSeries->name = timeSeries->name;
  freqSeries->sampleUnits = timeSeries->sampleUnits;
   */


  /* compute normalization factor */
  ASSERT(dftParams->sumofsquares > 0.0, status, TFTRANSFORMH_EPOSARG, 
         TFTRANSFORMH_MSGEPOSARG);
  fac = 1 / sqrt( dftParams->sumofsquares);


  /* create temporary vector */
  LALCreateVector (status->statusPtr, &tmp, n);
  CHECKSTATUSPTR (status);


  /* compute windowed version of time series data */
  for (i = 0; i < n; i++)
  {
    tmp->data[i] = fac*timeSeries->data->data[i]*dftParams->window->data[i];
  };

  /* compute the DFT */
  LALForwardRealFFT (status->statusPtr, freqSeries->data, tmp, dftParams->plan);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &tmp);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



