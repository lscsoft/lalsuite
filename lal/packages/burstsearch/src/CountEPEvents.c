/********* <lalVerbatim file="CountEPEventsCV"> ********
Author: Flanagan, E
$Id$
********** </lalVerbatim> ********/


#include <lal/LALRCSID.h>


NRCSID (COUNTEPEVENTSC, "$Id$");


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

/******** <lalVerbatim file="CountEPEventsCP"> ********/
void
LALCountEPEvents (
               LALStatus                               *status,
               INT4                                 *numEvents,
               TFTiling                             *tfTiling,
               REAL8                                alphaThreshold
               )
/******** </lalVerbatim> ********/
{
  INT4               tileCount;
  TFTile             *thisTile;

  INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
  ATTATCHSTATUSPTR (status);


  /* make sure that arguments are not NULL */
  ASSERT (tfTiling, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (tfTiling->firstTile, status, EXCESSPOWERH_ENULLP, 
          EXCESSPOWERH_MSGENULLP);
  ASSERT (numEvents, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (alphaThreshold > 0.0, status, EXCESSPOWERH_EPOSARG,
          EXCESSPOWERH_MSGEPOSARG);

  /* 
   *  Should call this routine after calling SortTFTiling, so tiles
   *  are already in order 
   *
   */

  /* check that already sorted */
  ASSERT( tfTiling->tilesSorted, status, EXCESSPOWERH_EORDER,
          EXCESSPOWERH_MSGEORDER);

  thisTile = tfTiling->firstTile;
  tileCount=0;

  while ( (thisTile != NULL) && (thisTile->alpha <= alphaThreshold))
    {
      tileCount++;
      thisTile = thisTile->nextTile;
    }
  
  /* return number of events */
  *numEvents = tileCount;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



