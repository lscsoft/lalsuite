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
#include <lal/BurstSearch.h>


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


/******** <lalVerbatim file="TFTilesToBurstEventsCP"> ********/
void
LALTFTileToBurstEvent (
               LALStatus                            *status,
               BurstEvent                           *burstEvent,
               TFTile                               *event,
               INT8                                  tstart,
               REAL4                                 flow
               )
/******** </lalVerbatim> ********/
{
  INT8 dummyNS;

  INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
  ATTATCHSTATUSPTR(status);

  /* make sure that arguments are not NULL */
  ASSERT (event, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (burstEvent, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);

  dummyNS = tstart + (INT8) (1e9 * event->tstart * event->deltaT);
  burstEvent->startTime        = (INT4) (dummyNS/1000000000L);
  burstEvent->startTimeNS      = (INT4) (dummyNS%1000000000L);
  burstEvent->duration         = (REAL4) (event->tend - event->tstart + 1) * 
    event->deltaT;
  burstEvent->centralFrequency = flow + 
    (REAL4) (event->fstart + event->fend + 1) / (2.0 * event->deltaT);
  burstEvent->bandwidth        = (REAL4) (event->fend - event->fstart + 1) / 
    (event->deltaT);
  burstEvent->amplitude        = event->excessPower;
  burstEvent->excessPower      = event->excessPower;
  burstEvent->confidence       = event->alpha;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


