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

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/Thresholds.h>
#include <lal/ExcessPower.h>
#include <lal/Random.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

/******** <lalVerbatim file="ComputeTFPlanesCP"> ********/
void
LALComputeTFPlanes (
		 LALStatus                             *status,
		 TFTiling                           *tfTiling,
		 COMPLEX8FrequencySeries            *freqSeries
		 )
/******** </lalVerbatim> ********/
{
  INT4               i;

  INITSTATUS (status, "LALComputeTFPlanes", COMPUTETFPLANESC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT (freqSeries, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (freqSeries->data, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (freqSeries->data->data, status, EXCESSPOWERH_ENULLP,
          EXCESSPOWERH_MSGENULLP);

  ASSERT (tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->numPlanes>0, status, EXCESSPOWERH_EPOSARG,
          EXCESSPOWERH_MSGEPOSARG);

  /* compute the ith TFPlane from input frequency series */
  for(i=0; i<tfTiling->numPlanes; i++)
    {
      COMPLEX8TimeFrequencyPlane    **thisPlane;
      ComplexDFTParams              **thisDftParams;
      HorizontalTFTransformIn       transformparams;

      thisPlane = tfTiling->tfp + i;
      ASSERT(thisPlane, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
      ASSERT(*thisPlane, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);

      thisDftParams = tfTiling->dftParams + i;
      ASSERT(thisDftParams, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
      ASSERT(*thisDftParams, status, EXCESSPOWERH_ENULLP, 
              EXCESSPOWERH_MSGENULLP);


      /* setup input structure for computing TF transform */
      transformparams.startT=0;  /* not used for horizontal transforms */
      transformparams.dftParams=*thisDftParams;

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

