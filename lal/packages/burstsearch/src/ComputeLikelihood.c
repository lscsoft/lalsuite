/******** <lalVerbatim file="ComputeLikelihoodCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>


NRCSID (COMPUTELIKELIHOODC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

/******** <lalVerbatim file="ComputeLikelihoodCP"> ********/
void
LALComputeLikelihood (
		   LALStatus                             *status,
		   REAL8                              *lambda,
		   TFTiling                           *tfTiling
		   )
/******** </lalVerbatim> ********/
{
  REAL8              avglambda=0.0;
  TFTile             *thisTile;


  INITSTATUS (status, "LALComputeLikelihood", COMPUTELIKELIHOODC);
  ATTATCHSTATUSPTR (status);


    /* make sure that arguments are not NULL */
  ASSERT(tfTiling, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(tfTiling->firstTile, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(lambda, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* make sure LALComputeExcessPower() has been already called */
  ASSERT(tfTiling->excessPowerComputed, status, EXCESSPOWERH_EORDER, EXCESSPOWERH_MSGEORDER);

  thisTile = tfTiling->firstTile;
  while (thisTile != NULL)
    {
      if(thisTile->firstCutFlag)
	{
	  INT4 t1;
	  INT4 t2;
	  INT4 f1;
	  INT4 f2;
	  REAL8 dof;
	  REAL8 loclambda;
	  REAL8 rho4;

	  t1=thisTile->tstart;
	  t2=thisTile->tend;
	  f1=thisTile->fstart;
	  f2=thisTile->fend;
	  dof = (REAL8)(2*(t2-t1+1)*(f2-f1+1));
	  rho4 = thisTile->excessPower * thisTile->excessPower;

	  loclambda = dof / (rho4 * thisTile->alpha);

	  avglambda += loclambda * thisTile->weight;
	}

      /* go onto next tile */
      thisTile = thisTile->nextTile;
    }
  
  /* compute the likelihood averaged over TF tiles */
  avglambda /= (REAL8)(tfTiling->numTiles);

  /* return value of statistic */
  *lambda=avglambda;
  

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



