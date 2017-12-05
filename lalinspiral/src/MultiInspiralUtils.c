/*-----------------------------------------------------------------------
 *
 * File Name: MultiInspiralUtils.c
 *
 * Author: Bose, S.
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

/**
 * \author Bose, S.
 * \file
 *
 * \brief Provides a set of utilities for manipulating \c multiInspiralTables.
 *
 */

/*
 * A few quickies for convenience.
 */

static INT8 end_time(const MultiInspiralTable *x)
{
	return(XLALGPSToINT8NS(&x->end_time));
}

#if 0
/* functions currently unused */
static INT4 end_time_sec(const MultiInspiralTable *x)
{
	return(x->end_time.gpsSeconds);
}

static INT4 end_time_nsec(const MultiInspiralTable *x)
{
	return(x->end_time.gpsNanoSeconds);
}
#endif


void
LALFreeMultiInspiral (
    LALStatus          *status,
    MultiInspiralTable **eventHead
    )

{
  INITSTATUS(status);
  XLALFreeMultiInspiral( eventHead );
  RETURN( status );
}


int
XLALFreeMultiInspiral (
    MultiInspiralTable **eventHead
    )

{
  EventIDColumn        *eventId;

  while ( (eventId = (*eventHead)->event_id) )
    {
      /* free any associated event_id's */
      (*eventHead)->event_id = (*eventHead)->event_id->next;
      LALFree( eventId );
    }
  LALFree( *eventHead );

  return (0);
}


MultiInspiralTable *
XLALSortMultiInspiral (
    MultiInspiralTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    )

{
  INT4                  i;
  INT4                  numEvents = 0;
  MultiInspiralTable    *thisEvent = NULL;
  MultiInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    XLALPrintInfo(
      "XLALSortMultiInspiral: Empty MultiInspiralTable passed as input\n" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (MultiInspiralTable **)
    LALCalloc( numEvents, sizeof(MultiInspiralTable *) );
  for ( i = 0, thisEvent = eventHead; i < numEvents;
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  /* qsort the array using the specified function */
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisEvent = eventHead = eventHandle[0];
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = thisEvent->next = eventHandle[i];
  }
  thisEvent->next = NULL;

  /* free the internal memory */
  LALFree( eventHandle );

  return( eventHead );
}


REAL4
XLALMultiInspiralStat(
    MultiInspiralTable         *multiInspiral,
    MultiInspiralClusterChoice  multiStat
    )
{
  REAL4 statValue = 0;

  if ( multiStat == cohsnr )
  {
    statValue = multiInspiral->snr;
  }
  else if ( multiStat == effCohSnr )
  {
    /* CHECK: Replace tau5 with a proper column for storing effCohSnr*/
    statValue = multiInspiral->tau5;
  }
  else if ( multiStat == nullstat )
  {
    /* CHECK: Invert it since clustering chooses the trigger with the
       MINIMUM null-stat value; add 0.001 to keep statValue from blowing up*/
    statValue = 1.0 / (multiInspiral->null_statistic + 1.0e-6);
  }
  else if ( multiStat == snrByNullstat )
  {
    REAL4 tmp_cohsnr = multiInspiral->snr;
    REAL4 tmp_nullstat = multiInspiral->null_statistic;
    if ( tmp_nullstat > 0.0 )
    {
      statValue = tmp_cohsnr / tmp_nullstat;
    }
    else
    {
      statValue = 1.0e6 * tmp_cohsnr;
    }
  }
  else if ( multiStat == autoCorrCohSqByNullstat )
  {
    REAL4 tmp_cohsnr = multiInspiral->autoCorrCohSq;
    REAL4 tmp_nullstat = multiInspiral->null_statistic;
    statValue = tmp_cohsnr / tmp_nullstat;
  }
  else if ( multiStat == crossCorrCohSqByNullstat  )
  {
    REAL4 tmp_cohsnr = multiInspiral->crossCorrCohSq;
    REAL4 tmp_nullstat = multiInspiral->null_statistic;
    statValue = tmp_cohsnr / tmp_nullstat;
  }
  else if ( multiStat == autoCorrNullSqByNullstat )
  {
    REAL4 tmp_cohsnr = multiInspiral->autoCorrNullSq;
    REAL4 tmp_nullstat = multiInspiral->null_statistic;
    statValue = tmp_cohsnr / tmp_nullstat;
  }
  else if ( multiStat == crossCorrNullSqByNullstat )
  {
    REAL4 tmp_cohsnr = multiInspiral->crossCorrNullSq;
    REAL4 tmp_nullstat = multiInspiral->null_statistic;
    statValue = tmp_cohsnr / tmp_nullstat;
  }
  else
  {
    statValue = 0;
  }
  return( statValue );
}



int
XLALClusterMultiInspiralTable (
    MultiInspiralTable         **inspiralList,
    INT8                         dtimeNS,
    MultiInspiralClusterChoice   clusterchoice
    )

{
  MultiInspiralTable     *thisEvent=NULL;
  MultiInspiralTable     *prevEvent=NULL;
  MultiInspiralTable     *nextEvent=NULL;
  int                    numMultiClust = 0;

  if ( !inspiralList )
  {
    XLAL_ERROR(XLAL_EIO);
  }

  if ( ! *inspiralList )
  {
    XLALPrintInfo(
      "XLALClusterMultiInspiralTable: Empty coincList passed as input\n" );
    return( 0 );
  }

  thisEvent = *inspiralList;
  nextEvent = (*inspiralList)->next;
  *inspiralList = NULL;

  while ( nextEvent )
  {
    INT8 thisTime = XLALGPSToINT8NS( &(thisEvent->end_time) );
    INT8 nextTime = XLALGPSToINT8NS( &(nextEvent->end_time) );;

    /* find events within the cluster window */
    if ( (nextTime - thisTime) < dtimeNS )
    {
      REAL4 thisStat = XLALMultiInspiralStat( thisEvent, clusterchoice );
      REAL4 nextStat = XLALMultiInspiralStat( nextEvent, clusterchoice );

      if ( nextStat > thisStat )
      {
        /* displace previous event in cluster */
        if( prevEvent )
        {
          prevEvent->next = nextEvent;
        }
        XLALFreeMultiInspiral( &thisEvent );
        thisEvent = nextEvent;
        nextEvent = thisEvent->next;
      }
      else
      {
        /* otherwise just dump next event from cluster */
        thisEvent->next = nextEvent->next;
        XLALFreeMultiInspiral ( &nextEvent );
        nextEvent = thisEvent->next;
      }
    }
    else
    {
      /* otherwise we keep this unique event trigger */
      if ( ! *inspiralList )
      {
        *inspiralList = thisEvent;
      }
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
      nextEvent = thisEvent->next;
      ++numMultiClust;
    }
  }

    /* store the last event */
  if ( ! (*inspiralList) )
  {
    *inspiralList = thisEvent;
  }
  ++numMultiClust;

  return(numMultiClust);
}


MultiInspiralTable *
XLALTimeCutMultiInspiral(
    MultiInspiralTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )

{
  MultiInspiralTable    *inspiralEventList = NULL;
  MultiInspiralTable    *thisEvent = NULL;
  MultiInspiralTable    *prevEvent = NULL;
  INT8                  startTimeNS = XLALGPSToINT8NS( startTime );
  INT8                  endTimeNS = XLALGPSToINT8NS( endTime );


  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = eventHead;

  while ( thisEvent )
  {
    MultiInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( end_time(tmpEvent) >= startTimeNS &&
        end_time(tmpEvent) < endTimeNS )
    {
      /* keep this template */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
    }
    else
    {
      /* discard this template */
      XLALFreeMultiInspiral ( &tmpEvent );
    }
  }
  eventHead = inspiralEventList;

  return (eventHead);
}



MultiInspiralTable *
XLALSNRCutMultiInspiral (
    MultiInspiralTable          *eventHead,
    REAL4                       snrCut
    )

{
  MultiInspiralTable    *thisEvent = NULL;
  MultiInspiralTable    *prevEvent = NULL;

  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    MultiInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->snr >= snrCut )
    {
      /* keep this template */
      if ( ! eventHead  )
      {
        eventHead = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
    }
    else
    {
      /* discard this template */
      XLALFreeMultiInspiral ( &tmpEvent );
    }
  }
  return( eventHead );
}



MultiInspiralTable *
XLALPlayTestMultiInspiral(
    MultiInspiralTable          *eventHead,
    LALPlaygroundDataMask       *dataType
    )

{
  MultiInspiralTable    *inspiralEventList = NULL;
  MultiInspiralTable    *thisEvent = NULL;
  MultiInspiralTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      MultiInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALGPSToINT8NS( &(tmpEvent->end_time) );
      isPlay = XLALINT8NanoSecIsPlayground( triggerTime );

      if ( ( (*dataType == playground_only)  && isPlay ) ||
          ( (*dataType == exclude_play) && ! isPlay) )
      {
        /* keep this trigger */
        if ( ! inspiralEventList  )
        {
          inspiralEventList = tmpEvent;
        }
        else
        {
          prevEvent->next = tmpEvent;
        }
        tmpEvent->next = NULL;
        prevEvent = tmpEvent;
        ++numTriggers;
      }
      else
      {
        /* discard this template */
        XLALFreeMultiInspiral ( &tmpEvent );
      }
    }
    eventHead = inspiralEventList;
    if ( *dataType == playground_only )
    {
      XLALPrintInfo( "Kept %d playground triggers \n", numTriggers );
    }
    else if ( *dataType == exclude_play )
    {
      XLALPrintInfo( "Kept %d non-playground triggers \n", numTriggers );
    }
  }
  else if ( *dataType == all_data )
  {
    XLALPrintInfo( "Keeping all triggers since all_data specified\n" );
  }
  else
  {
    XLALPrintInfo( "Unknown data type, returning no triggers\n" );
    eventHead = NULL;
  }

  return(eventHead);
}

/*CHECK: This function is not needed since XLALCountMultiInspiralTable
  in LIGOMetadataUtils in support already does that job */

INT4 XLALCountMultiInspiral( MultiInspiralTable *head )

{
  INT4 length;
  MultiInspiralTable *event;

  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}


int
LALCompareMultiInspiralByTime (
    const void *a,
    const void *b
    )

{
  LALStatus     status;
  const MultiInspiralTable *aPtr = *((const MultiInspiralTable * const *)a);
  const MultiInspiralTable *bPtr = *((const MultiInspiralTable * const *)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  ta = XLALGPSToINT8NS( &(aPtr->end_time) );
  tb = XLALGPSToINT8NS( &(bPtr->end_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


int
XLALMultiInspiralIfos (
    MultiInspiralTable  *multiInspiral,
    char                *ifos
    )

{
  int                   ifosMatch  = 0;

  if ( !multiInspiral )
  {
    return ( 0 );
  }

  /* check that the multi inspiral trigger is of the correct type */
  if ( strcmp(multiInspiral->ifos, ifos) == 0 )
  {
    ifosMatch = 1;
  }
  return( ifosMatch );
}


int
XLALMultiInspiralIfosCut(
    MultiInspiralTable **multiHead,
    char                *ifos
    )

{
  MultiInspiralTable    *prevTrig = NULL;
  MultiInspiralTable    *thisTrig = NULL;
  int                    numTrig = 0;

  thisTrig = *multiHead;
  *multiHead = NULL;

  while ( thisTrig )
  {
    MultiInspiralTable *tmpTrig = thisTrig;
    thisTrig = thisTrig->next;

    if ( XLALMultiInspiralIfos( tmpTrig, ifos ) )
    {
      /* ifos match so keep tmpTrig */
      if ( ! *multiHead  )
      {
        *multiHead = tmpTrig;
      }
      else
      {
        prevTrig->next = tmpTrig;
      }
      tmpTrig->next = NULL;
      prevTrig = tmpTrig;
      ++numTrig;
    }
    else
    {
      /* discard tmpCoinc */
      XLALFreeMultiInspiral( &tmpTrig );
    }
  }

  return( numTrig );
}


MultiInspiralTable *
XLALMultiInspiralSlideCut(
    MultiInspiralTable **multiHead,
    int                  slideNum
    )
{
  MultiInspiralTable    *prevMulti      = NULL;
  MultiInspiralTable    *thisMulti      = NULL;
  MultiInspiralTable    *slideHead      = NULL;
  MultiInspiralTable    *thisSlideMulti = NULL;

  UINT8 idNumber = 0;

  if( slideNum < 0 )
  {
    slideNum = 5000 - slideNum;
  }

  thisMulti = *multiHead;
  *multiHead = NULL;

  while ( thisMulti )
  {
    idNumber = thisMulti->event_id->id;

    if ( (int) ((idNumber % 1000000000) / 100000) == slideNum )
    {
      /* add thisMulti to the slideMulti list */
      if ( slideHead )
      {
        thisSlideMulti = thisSlideMulti->next = thisMulti;
      }
      else
      {
        slideHead = thisSlideMulti = thisMulti;
      }

      /* remove from multiHead list */
      if ( prevMulti )
      {
        prevMulti->next = thisMulti->next;
      }

      thisMulti = thisMulti->next;
      thisSlideMulti->next = NULL;
    }
    else
    {
      /* move along the list */
      if( ! *multiHead )
      {
        *multiHead = thisMulti;
      }

      prevMulti = thisMulti;
      thisMulti = thisMulti->next;
    }
  }
  return( slideHead );
}
