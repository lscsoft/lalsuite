/********* <lalVerbatim file="CountEPEventsCV"> ********
Author: Flanagan, E and Brady, P R
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
#include <lal/EPSearch.h>
#include <lal/LIGOMetadataTables.h>

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
               SnglBurstTable                       *burstEvent,
               TFTile                               *event,
               INT8                                  tstart,
               EPSearchParams                       *params  
               )
/******** </lalVerbatim> ********/
{
  INT8 dummyNS, peakNS;
  REAL4 flow;

  INITSTATUS (status, "LALCountEPEvents", COUNTEPEVENTSC);
  ATTATCHSTATUSPTR(status);

  /* make sure that arguments are not NULL */
  ASSERT (event, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);
  ASSERT (burstEvent, status, EXCESSPOWERH_ENULLP, EXCESSPOWERH_MSGENULLP);

  flow=(REAL4)params->tfTilingInput->flow;
  dummyNS = tstart + (INT8) (1e9 * event->tstart * event->deltaT);
  burstEvent->start_time.gpsSeconds     = (INT4) (dummyNS/1000000000L);
  burstEvent->start_time.gpsNanoSeconds = (INT4) (dummyNS%1000000000L);
  burstEvent->duration         = (REAL4) (event->tend - event->tstart + 1) * 
    event->deltaT;
  peakNS = dummyNS + (0.5 * burstEvent->duration * 1e9);
  burstEvent->peak_time.gpsSeconds     = (INT4) (peakNS/1000000000L);
  burstEvent->peak_time.gpsNanoSeconds = (INT4) (peakNS%1000000000L);
  burstEvent->central_freq = flow + 
    (REAL4) (event->fstart + event->fend + 1) / (2.0 * event->deltaT);
  burstEvent->bandwidth        = (REAL4) (event->fend - event->fstart + 1) / 
    (event->deltaT);
  burstEvent->amplitude        = event->excessPower;
  burstEvent->snr              = event->excessPower;
  burstEvent->next             = NULL;

  /* THIS NEEDS TO FIXED */
  sprintf( burstEvent->ifo ,"XX");
  strncpy( burstEvent->ifo,  params->channelName, 2);
  sprintf( burstEvent->search, "power");
  sprintf( burstEvent->channel, params->channelName);

  /* report log(confidence) to avoid problems with REAL4's */
  /* changement by Stefan Ballmer and Erik Katsavounidis, 9/16/2002 */
  if(event->alpha<0) {
      burstEvent->confidence     =   LAL_REAL4_MAX;
  } else {
      if(event->alpha==0) {
          burstEvent->confidence   = - LAL_REAL4_MAX;
      } else {
          burstEvent->confidence   =   log(event->alpha);
      }
  }
  /* end of Stefan's changement */

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


