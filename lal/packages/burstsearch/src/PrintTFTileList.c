/********* <lalVerbatim file="PrintTFTileListCV"> ********
Author: Eanna Flanagan
$Id$
********* </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (PRINTTFTILELISTC, "$Id$");


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
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

static void 
PrintTFTile (
	     LALStatus                                 *status,
	     FILE                                   *fp,
	     TFTile                                 *tfTile,
	     TFTiling                               *tfTiling
	     )
{
  /*
   *  this function prints out information contained in one time-frequency
   *  tile.
   */

  INT4 t1;
  INT4 t2;
  INT4 f1;
  INT4 f2;
  INT4 dof;
  COMPLEX8TimeFrequencyPlane *thisPlane=NULL;
  REAL8 flow;
  REAL8 epoch;
  REAL8 deltaT;
  

  INITSTATUS (status, "PrintTFTile", PRINTTFTILELISTC);
  ATTATCHSTATUSPTR (status);

  ASSERT(tfTile, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->tfp, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(fp, status, LAL_NULL_ERR, LAL_NULL_MSG);

  t1 = tfTile->tstart;
  t2 = tfTile->tend;
  f1 = tfTile->fstart;
  f2 = tfTile->fend;
  dof = 2*(t2-t1+1)*(f2-f1+1);
  thisPlane = *(tfTiling->tfp + tfTile->whichPlane);
  flow = thisPlane->params->flow;
  epoch = (REAL8)(thisPlane->epoch.gpsSeconds) + (REAL8)(thisPlane->epoch.gpsNanoSeconds)/1000000000.0;
  deltaT = thisPlane->params->deltaT;

  fprintf(fp,"\n");
  fprintf(fp," Time frequency tile in TF Plane number %d\n", 
          tfTile->whichPlane);
  fprintf(fp," Frequency interval: %f Hz to %f Hz, i.e.,  bins %d to %d of %d\n",
	 flow + (REAL8)(f1)/deltaT, flow + (REAL8)(f2+1)/deltaT, f1, f2,
	 thisPlane->params->freqBins);
  fprintf(fp," Time interval    :  %f s to %f s, i.e.,  bins %d to %d of %d\n",
	 epoch + (REAL8)(t1)*deltaT, epoch + (REAL8)(t2+1)*deltaT, t1, t2,
	 thisPlane->params->timeBins);
  fprintf(fp," Total number of degrees of freedom:  %d\n",dof);
  fprintf(fp," Excess power:  %f,   1 / alpha    :  %f\n", tfTile->excessPower,
         1/tfTile->alpha);
  {
    /* print out effective number of sigma */
    REAL8 sigma2;
    Chi2ThresholdIn input;

    input.dof = 1.0;
    input.falseAlarm = tfTile->alpha;
    LALChi2Threshold (status->statusPtr, &sigma2, &input);
    CHECKSTATUSPTR (status);
    fprintf(fp," Effective number of sigma:  %f\n", sqrt(sigma2));
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/******** <lalVerbatim file="PrintTFTileListCP"> ********/
void 
LALPrintTFTileList (
		 LALStatus                                 *status,
		 FILE                                   *fp,
		 TFTiling                               *tfTiling,
		 INT4                                   maxTiles
		 )
/******** </lalVerbatim> ********/
{
  TFTile *thisTile;
  INT4   tileCount=0;

  INITSTATUS (status, "LALPrintTFTileList", PRINTTFTILELISTC);
  ATTATCHSTATUSPTR (status);

  ASSERT(fp, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG); 
  ASSERT(tfTiling->firstTile, status, LAL_NULL_ERR, LAL_NULL_MSG); 
  ASSERT(tfTiling->excessPowerComputed, status, EXCESSPOWERH_EORDER, EXCESSPOWERH_MSGEORDER);


  thisTile = tfTiling->firstTile;
  while ( (thisTile != NULL) && (tileCount < maxTiles))
    {
      PrintTFTile (status->statusPtr, fp, thisTile, tfTiling);
      CHECKSTATUSPTR(status);
      tileCount++;
      thisTile = thisTile->nextTile;
    }

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




















