/******** <lalVerbatim file="ComputeTFPlanesCV">
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/
 
#include <lal/LALRCSID.h>


NRCSID (COMPUTETFPLANESC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/ExcessPower.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

/******** <lalVerbatim file="ComputeTFPlanesCP"> ********/
void
LALComputeTFPlanes (
		 LALStatus                             *status,
		 TFTiling                           *tfTiling,
		 COMPLEX8FrequencySeries            *freqSeries,
		 UINT4                              windowShift
		 )
/******** </lalVerbatim> ********/
{
  INT4               i;

  INITSTATUS (status, "LALComputeTFPlanes", COMPUTETFPLANESC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data->data, status, LAL_NULL_ERR, LAL_NULL_MSG);

  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->firstTile, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->numPlanes > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  /* compute the ith TFPlane from input frequency series */
  for(i=0; i<tfTiling->numPlanes; i++)
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;
      HorizontalTFTransformIn       transformparams;

      thisPlane = tfTiling->tfp + i;
      ASSERT(thisPlane, status, LAL_NULL_ERR, LAL_NULL_MSG);
      ASSERT(*thisPlane, status, LAL_NULL_ERR, LAL_NULL_MSG);

      thisDftParams = tfTiling->dftParams + i;
      ASSERT(thisDftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
      ASSERT(*thisDftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);


      /* setup input structure for computing TF transform */
      transformparams.startT=0;  /* not used for horizontal transforms */
      transformparams.dftParams=*thisDftParams;
      transformparams.windowShift = windowShift;

      /* Compute TF transform */
      LALInfo(status->statusPtr, "Converting Frequency series to TFPlane");
      CHECKSTATUSPTR (status);
      LALModFreqSeriesToTFPlane( status->statusPtr, *thisPlane, freqSeries, 
                           &transformparams); 
      CHECKSTATUSPTR (status);
    }

  /* set flags saying TF planes have been computed, but not EP or sorted */
  tfTiling->planesComputed=TRUE;
  tfTiling->excessPowerComputed=FALSE;
  tfTiling->tilesSorted=FALSE;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALModComputeTFPlanes (
		 LALStatus                             *status,
		 TFTiling                           *tfTiling,
		 COMPLEX8FrequencySeries            *freqSeries,
		 UINT4                              windowShift
		 )
/******** </lalVerbatim> ********/
{
 
  INITSTATUS (status, "LALComputeTFPlanes", COMPUTETFPLANESC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT(freqSeries, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(freqSeries->data->data, status, LAL_NULL_ERR, LAL_NULL_MSG);

  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->firstTile, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->numPlanes > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);


  /* compute the TFPlane from input frequency series */
  
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;
      HorizontalTFTransformIn       transformparams;

      thisPlane = tfTiling->tfp;
      ASSERT(thisPlane, status, LAL_NULL_ERR, LAL_NULL_MSG);
      ASSERT(*thisPlane, status, LAL_NULL_ERR, LAL_NULL_MSG);

      thisDftParams = tfTiling->dftParams;
      ASSERT(thisDftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);
      ASSERT(*thisDftParams, status, LAL_NULL_ERR, LAL_NULL_MSG);


      /* setup input structure for computing TF transform */
      transformparams.startT=0;  /* not used for horizontal transforms */
      transformparams.dftParams=*thisDftParams;
      transformparams.windowShift = windowShift;
      /* Compute TF transform */
      LALInfo(status->statusPtr, "Converting Frequency series to TFPlane");
      CHECKSTATUSPTR (status);
      LALModModFreqSeriesToTFPlane( status->statusPtr, *thisPlane, freqSeries, 
                           &transformparams); 
      CHECKSTATUSPTR (status);
    }

  /* set flags saying TF planes have been computed, but not EP or sorted */
  tfTiling->planesComputed=TRUE;
  tfTiling->excessPowerComputed=FALSE;
  tfTiling->tilesSorted=FALSE;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

