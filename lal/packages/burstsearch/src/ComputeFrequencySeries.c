/******** <lalVerbatim file="ComputeFrequencySeriesCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (COMPUTEFREQUENCYSERIESC, "$Id$");


#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <lal/ComplexFFT.h>
#include <lal/LALStdlib.h>
#include <lal/LALErrno.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
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
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(timeSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(timeSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(dftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(dftParams->plan, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(dftParams->window, status, LAL_NULL_ERR, LAL_NULL_MSG);


  /* make sure sizes are reasonable and agree */
  ASSERT(timeSeries->data->length > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  n = timeSeries->data->length;

  ASSERT(freqSeries->data->length == n/2 + 1, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  ASSERT(dftParams->window->length == n, status, LAL_RANGE_ERR, LAL_RANGE_MSG);


  /* copy over data into frequency series structure */
  freqSeries->epoch = timeSeries->epoch;
  freqSeries->f0    = timeSeries->f0;
  ASSERT(timeSeries->deltaT > 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  freqSeries->deltaF = 1/((REAL8)(timeSeries->data->length)*timeSeries->deltaT); 
  /* compute normalization factor */
  ASSERT(dftParams->sumofsquares > 0.0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);
  fac = sqrt(n)*timeSeries->deltaT / sqrt( dftParams->sumofsquares);

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



