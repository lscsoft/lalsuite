/**** <lalVerbatim file="BurstSearchHV"> ***********
Author:Brady, P.
$Id$
***** </lalVerbatim> ***********************************/

#ifndef _BURSTSEARCH_H
#define _BURSTSEARCH_H

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (BURSTSEARCHH, "$Id$");


typedef struct
tagBurstEvent
{
  INT4                                  startTime;
  INT4                                  startTimeNS;
  REAL4                                 duration;
  REAL4                                 centralFrequency;
  REAL4                                 bandwidth;
  REAL4                                 amplitude;
  REAL4                                 excessPower;
  REAL4                                 confidence;
  struct tagBurstEvent                 *nextEvent;
}
BurstEvent;


typedef struct
tagBurstEventList
{
  BurstEvent                            *burstEvent;
  CHAR                                  *ifo;
  CHAR                                  *search;
  INT4                                   numEvents;
}
BurstEventList;

void
LALTFTileToBurstEvent (
               LALStatus                            *status,
               BurstEvent                           *burstEvent,
               TFTile                               *event,
               INT8                                  tstart,
               REAL4                                 flow
               );

#ifdef  __cplusplus
}
#endif

#endif /* _BURSTSEARCH_H */
