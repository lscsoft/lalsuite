/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglBurstUtils.c
 *
 * Author: Brown, D. A.  and Brady, P. R.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglBurstUtilsCV">
Author: Brown, D. A. and Brady, P. R.
$Id$
</lalVerbatim> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/BurstSearch.h>

NRCSID( SNGLBURSTUTILSC, "$Id$" );

#define NANOSEC  (1000000000LL)

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglBurstUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglBurstUtilsCP}
\idx{LALSortSnglBurst()}
\idx{LALCompareSnglBurstByMass()}
\idx{LALCompareSnglBurstByTime()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.
 
\vfill{\footnotesize\input{SnglBurstUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALSortSnglBurst(
    LALStatus          *status,
    SnglBurstTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglBurstTable    *thisEvent = NULL;
  SnglBurstTable   **eventHandle = NULL;

  INITSTATUS( status, "LALSortSnglBurst", SNGLBURSTUTILSC );

  /* count the number of events in the linked list */
  for ( thisEvent = *eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    LALWarning( status, "No events in list to sort" );
    RETURN( status );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (SnglBurstTable **) 
    LALCalloc( numEvents, sizeof(SnglBurstTable *) );
  for ( i = 0, thisEvent = *eventHead; i < numEvents; 
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  /* qsort the array using the specified function */
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), comparfunc );
  
  /* re-link the linked list in the right order */
  thisEvent = *eventHead = eventHandle[0];
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = thisEvent->next = eventHandle[i];
  }
  thisEvent->next = NULL;

  /* free the internal memory */
  LALFree( eventHandle );

  RETURN( status );
}



/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
LALCompareSnglBurstByTime(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  SnglBurstTable *aPtr = *((SnglBurstTable **)a);
  SnglBurstTable *bPtr = *((SnglBurstTable **)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

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


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
LALCompareSnglBurstByTimeAndFreq(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  SnglBurstTable *aPtr = *((SnglBurstTable **)a);
  SnglBurstTable *bPtr = *((SnglBurstTable **)b);
  INT8 ta, tb;
  REAL4 flowa, flowb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

  flowa = aPtr->central_freq - 0.5 * aPtr->bandwidth;
  flowb = bPtr->central_freq - 0.5 * bPtr->bandwidth;

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta == tb && flowa > flowb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else if ( ta == tb && flowa < flowb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSnglBurst(
    LALStatus             *status,
    SnglBurstTable        *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta, tb;
  REAL4 dm1, dm2;
  REAL4 sigmaRatio;

  INITSTATUS (status, "LALCompareSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  params->match=0;

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->start_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->start_time) );

  if( labs(ta-tb) < params->dtime ){
    params->match = 1;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALClusterSnglBurstTable (
	      LALStatus        *status,
              SnglBurstTable   *burstEvent
	      )
/* </lalVerbatim> */
{
  SnglBurstTable     *thisEvent=NULL,*prevEvent=NULL;

  INITSTATUS (status, "LALClusterSnglBurstTable", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  thisEvent = burstEvent->next;
  prevEvent = burstEvent;
  while (thisEvent != NULL)
  {
    /* cornerTime[0][0] = thisEvent->start_time,  etc */
    INT8    cornerTime[2][2];
    REAL4   cornerFreq[2][2];

    /* compute the time in nanosec for each event trigger */
    LALGPStoINT8(status->statusPtr, &(cornerTime[1][0]), &(thisEvent->start_time));
    CHECKSTATUSPTR(status);
    cornerTime[1][1] = cornerTime[1][0] + ( NANOSEC * thisEvent->duration );

    LALGPStoINT8(status->statusPtr, &(cornerTime[0][0]), &(prevEvent->start_time));
    CHECKSTATUSPTR(status);
    cornerTime[0][1] = cornerTime[0][0] + ( NANOSEC * prevEvent->duration );

    /* compute the start and stop frequencies */
    cornerFreq[1][0] = thisEvent->central_freq - 0.5 * thisEvent->bandwidth;
    cornerFreq[1][1] = cornerFreq[1][0] + thisEvent->bandwidth;
    cornerFreq[0][0] = prevEvent->central_freq - 0.5 * prevEvent->bandwidth;
    cornerFreq[0][1] = cornerFreq[0][0] + prevEvent->bandwidth;

    /* find overlapping events */
    if ( ( (cornerTime[1][0]-cornerTime[0][0]) * 
          (cornerTime[1][0]-cornerTime[0][1]) <= 0 ) &&
        ( ((cornerFreq[1][0]-cornerFreq[0][0]) * 
           (cornerFreq[1][0]-cornerFreq[0][1]) <=0 ) ||
          ((cornerFreq[1][1]-cornerFreq[0][0]) * 
           (cornerFreq[1][1]-cornerFreq[0][1]) <=0 ) ) )
    {
      cornerFreq[0][0] = cornerFreq[1][0] < cornerFreq[0][0] ? cornerFreq[1][0] :
        cornerFreq[0][0];
      cornerFreq[0][1] = cornerFreq[1][1] > cornerFreq[0][1] ? cornerFreq[1][1] :
        cornerFreq[0][1];
      cornerTime[0][1] = cornerTime[1][1] > cornerTime[0][1] ? cornerTime[1][1] :
        cornerTime[0][1];
      prevEvent->central_freq = 0.5 * (cornerFreq[0][0]+cornerFreq[0][1]);
      prevEvent->bandwidth = (cornerFreq[0][1]-cornerFreq[0][0]);
      prevEvent->duration = (REAL4)(cornerTime[0][1]-cornerTime[0][0])/1.0e9;

      if ( prevEvent->confidence > thisEvent->confidence )
      {
        prevEvent->confidence = thisEvent->confidence;
      }

      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFree(thisEvent);
      thisEvent = prevEvent->next;
    }
    else
    {
      /* otherwise keep this as a unique trigger */
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

#undef NANOSEC
