/******** <lalVerbatim file="AddWhiteNoiseCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>


NRCSID (ADDWHITENOISEC, "$Id$");


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <lal/ExcessPower.h>
#include <lal/LALErrno.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/Thresholds.h>


#define TRUE 1
#define FALSE 0


extern INT4 lalDebugLevel;

/******** <lalVerbatim file="AddWhiteNoiseCP"> ********/
void
LALAddWhiteNoise (
	       LALStatus        *status,
	       COMPLEX8Vector   *v,
	       REAL8             noiseLevel
	       )
/******** </lalVerbatim> ********/
{
  /*
   *
   *  Add white noise to complex vector
   *
   */

  RandomParams           *params=NULL;
  REAL4Vector            *vr=NULL;
  REAL4Vector            *vi=NULL;
  INT4                   i;

  INITSTATUS (status, "LALAddWhiteNoise", ADDWHITENOISEC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT(v, status, LAL_NULL_ERR, LAL_NULL_MSG);
  ASSERT(v->data, status, LAL_NULL_ERR, LAL_NULL_MSG);

  /* make sure length of series is nonzero */
  ASSERT(v->length > 0, status, LAL_RANGE_ERR, LAL_RANGE_MSG);

  
  /* Seed Random Number Generator with current time for seed */
  LALCreateRandomParams (status->statusPtr, &params, 0);  
  CHECKSTATUSPTR (status);

  /* create temporary vectors */
  LALSCreateVector (status->statusPtr, &vr, v->length);
  CHECKSTATUSPTR (status);
  LALSCreateVector (status->statusPtr, &vi, v->length);
  CHECKSTATUSPTR (status);
  
  /* Fill temporary vectors with Gaussian deviates */
  LALNormalDeviates (status->statusPtr, vr, params);
  CHECKSTATUSPTR (status);
  LALNormalDeviates (status->statusPtr, vi, params);
  CHECKSTATUSPTR (status);

  for(i=0;i<(INT4)v->length;i++) 
    {
      v->data[i].re += noiseLevel * vr->data[i];
      v->data[i].im += noiseLevel * vi->data[i];
    }
    
  LALSDestroyVector (status->statusPtr, &vr);
  CHECKSTATUSPTR (status);
  
  LALSDestroyVector (status->statusPtr, &vi);
  CHECKSTATUSPTR (status);
  
  LALDestroyRandomParams (status->statusPtr, &params);
  CHECKSTATUSPTR (status);
  
  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}    


