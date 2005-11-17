/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglRingdownUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., Fairhurst, S. and Messaritaki, E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglRingdownUtilsCV">
Author: Brown, D. A., Fairhurst, S. and Messaritaki, E.
$Id$
</lalVerbatim> 
#endif


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

NRCSID( SNGLRINGDOWNUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglRingdownUtils.c}}

Provides a set of utilities for manipulating \texttt{snglRingdownTable}s.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglRingdownUtilsCP}
\idx{LALSortSnglRingdown()}
\idx{LALCompareSnglRingdownByMass()}
\idx{LALCompareSnglRingdownByPsi()}
\idx{LALCompareSnglRingdownByTime()}
\idx{LALCompareSnglRingdown()}
\idx{LALCompareRingdowns()}
\idx{LALClusterSnglRingdownTable()}
\idx{LALTimeCutSingleRingdown()}
\idx{LALalphaFCutSingleInspiral()}
\idx{LALIfoCutSingleInspiral()}
\idx{LALIfoCountSingleInspiral()}
\idx{LALTimeSlideSingleInspiral()} 
\idx{LALPlayTestSingleInspiral()}
\idx{LALCreateTrigBank()}
\idx{LALIncaCoincidenceTest()}
\idx{LALTamaCoincidenceTest()}


\subsubsection*{Description}

The function \texttt{LALFreeSnglInspiral()} frees the memory associated to a
single inspiral table.  The single inspiral table may point to a linked list
of EventIDColumns.  Thus, it is necessary to free all event ids associated
with the single inspiral.

The function \texttt{LALSortSnglInspiral()} sorts a list of single inspiral
tables.  The function simply calls qsort with the appropriate comparison
function, \texttt{comparfunc}.  It then ensures that the head of the sorted
list is returned.  There then follow several comparison functions for single
inspiral tables.  \texttt{LALCompareSnglInspiralByMass ()} first compares the
\texttt{mass1} entry of the two inspiral tables, returning 1 if the first mass
is larger and -1 if the second is larger.  In the case that the \texttt{mass1}
fields are equal, a similar comparsion is performed on \texttt{mass2}.  If
these also agree, 0 is returned.  \texttt{LALCompareSnglInspiralByPsi()}
compares the \texttt{Psi0} and \texttt{Psi3} fields in two single inspiral
tables.  The function is analogous to the mass comparison described above.
\texttt{LALCompareSnglInspiralByTime} compares the end times of two single
inspiral tables, returnng 1 if the first time is larger, 0 if equal and -1 if
the second time is larger.

\texttt{LALCompareSnglInspiral()} tests whether two single inspiral tables
pass a coincidence test.  The coincidence parameters are given by
\texttt{params} which is a \texttt{SnglInspiralAccuracy} structure.  It tests
first that the \texttt{ifo} fields are different.  If they are, it then tests
for time and mass coincidence, where mass coincidence may be any one of
\texttt{psi0\_and\_psi3}, \texttt{m1\_and\_m2}, \texttt{mchirp\_and\_eta}.
Finally, if the test is on \texttt{m1\_and\_m2}, consistency of effective
distances is also checked.  If the two single inspiral tables pass
coincidences the \texttt{params.match} is set to 1, otherwise it is set to
zero.

\texttt{LALClusterSnglInspiralTable ()} clusters single inspiral triggers
within a time window \texttt{dtimeNS}.  The triggers are compared either by
\texttt{snr}, \texttt{snr\_and\_chisq} or \texttt{snrsq\_over\_chisq}.  The
"loudest" trigger, as determined by the selected algorithm, within each time
window is returned.

\texttt{LALTimeCutSingleInspiral()} takes in a linked list of single inspiral
tables and returns only those which occur after the given \texttt{startTime}
and before the \texttt{endTime}.

\texttt{LALalphaFCutSingleInspiral()} takes in a linked list of single
inspiral tables and returns only those triggers which have alphaF values below
a specific alphaFcut. It is relevant for the BCV search only.

\texttt{LALIfoCutSingleInspiral()} scans through a linked list of single
inspiral tables and returns those which are from the requested \texttt{ifo}.
On input, \texttt{eventHead} is a pointer to the head of a linked list of
single inspiral tables.  On output, this list contains only single inspirals
from the requested \texttt{ifo}.

\texttt{LALIfoCountSingleInspiral()} scans through a linked list of single
inspiral tables and counts the number which are from the requested IFO.  
This count is returned as \texttt{numTrigs}.

\texttt{LALTimeSlideSingleInspiral()} performs a time slide on the triggers
contained in the \texttt{triggerList}.  The time slide for each instrument is
specified by \texttt{slideTimes[LAL\_NUM\_IFO]}.  If \texttt{startTime} and
\texttt{endTime} are specified, then the time slide is performed on a ring.  If
the slide takes any trigger outside of the window
\texttt{[startTime,endTime]}, then the trigger is wrapped to be in
this time window.

\texttt{LALPlayTestSingleInspiral()} tests whether single inspiral events
occured in playground or non-playground times.  It then returns the requested
subset of events which occurred in the times specified by \texttt{dataType}
which must be one of \texttt{playground\_only}, \texttt{exclude\_play} or
\texttt{all\_data}.  

\texttt{LALCreateTrigBank()} takes in a list of single inspiral tables and
returns a template bank.  The function tests whether a given template produced
multiple triggers.  If it did, only one copy of the template is retained.
Triggers are tested for coincidence in \texttt{m1\_and\_m2} or
\texttt{psi0\_and\_psi3}. 

\texttt{LALIncaCoincidenceTest()} performs a coincidence test between triggers
from two interferometers.  It tests pairs of events for both time and mass
coincidence and returns two equal length lists of coincident events.  Note
that if an event in one detector is coincident with several events in the
other detector, the output lists will contain several copies of this event.

\texttt{LALTamaCoincidenceTest()} also performs a coincidence test between
triggers from two interferometers, but with a slightly different coincidence
test.  First, it locates all triggers in the second instrument which are
coincident with triggers in the first instrument.  Then, it clusters these
triggers using the appropriate \texttt{clusterchioce}.  Finally, it tests for
mass coincidence between the first trigger and the clustered trigger from the
second instrument.


\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent LALCalloc, LALFree, LALGPStoINT8, LALINT8NanoSecIsPlayground.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SnglInspiralUtilsCV}}

</lalLaTeX>
#endif

/*
 * A few quickies for convenience.
 */

static INT8 start_time(const SnglRingdownTable *x)
{
	return(XLALGPStoINT8(&x->start_time));
}

static INT4 start_time_sec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsSeconds);
}

static INT4 start_time_nsec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsNanoSeconds);
}


/* <lalVerbatim file="SnglRingdownutilsCP"> */
void
LALFreeSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;

  INITSTATUS( status, "LALFreeSnglRingdown", SNGLRINGDOWNUTILSC );
  while ( (*eventHead)->event_id )
  {
    /* free any associated event_id's */
    eventId = (*eventHead)->event_id;
    (*eventHead)->event_id = (*eventHead)->event_id->next;
    LALFree( eventId );
  }
  LALFree( *eventHead );
  RETURN( status );
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
int
XLALFreeSnglRingdown (
    SnglRingdownTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;

  while ( (*eventHead)->event_id )
  {
    /* free any associated event_id's */
    eventId = (*eventHead)->event_id;
    (*eventHead)->event_id = (*eventHead)->event_id->next;
    LALFree( eventId );
  }
  LALFree( *eventHead );

  return (0);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALSortSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable   **eventHandle = NULL;

  INITSTATUS( status, "LALSortSnglRingdown", SNGLRINGDOWNUTILSC );

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
  eventHandle = (SnglRingdownTable **) 
    LALCalloc( numEvents, sizeof(SnglRingdownTable *) );
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


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
int
LALCompareSnglRingdownByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  const SnglRingdownTable *aPtr = *((const SnglRingdownTable * const *)a);
  const SnglRingdownTable *bPtr = *((const SnglRingdownTable * const *)b);
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



/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALClusterSnglRingdownTable (
    LALStatus                  *status,
    SnglRingdownTable          *ringdownEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  SnglRingdownTable     *thisEvent=NULL;
  SnglRingdownTable     *prevEvent=NULL;

  INITSTATUS( status, "LALClusterSnglRingdownTable", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( ringdownEvent, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  thisEvent = ringdownEvent->next;
  prevEvent = ringdownEvent;

  while ( thisEvent )
  {
    INT8 currTime;
    INT8 prevTime;

    /* compute the time in nanosec for each event trigger */
    LALGPStoINT8(status->statusPtr, &currTime, &(thisEvent->start_time));
    CHECKSTATUSPTR(status);

    LALGPStoINT8(status->statusPtr, &prevTime, &(prevEvent->start_time));
    CHECKSTATUSPTR(status);

    /* find events within the cluster window */
    if ( (currTime - prevTime) < dtimeNS )
    {
      /* displace previous event in cluster */
      if ( (clusterchoice == snr) && (thisEvent->snr > prevEvent->snr))
      {
        memcpy( prevEvent, thisEvent, sizeof(SnglRingdownTable) );
        thisEvent->event_id = NULL;
      }
      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFreeSnglRingdown ( status->statusPtr, &thisEvent );
      thisEvent = prevEvent->next;
    }
    else 
    {
      /* otherwise we keep this unique event trigger */
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALIfoCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    CHAR                       *ifo
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *eventList = NULL;
  SnglRingdownTable    *prevEvent = NULL;
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoScanSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that eventHead is non-null */
  ASSERT( eventHead, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_ringdown tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent = *eventHead;
  
  while ( thisEvent )
  {
    SnglRingdownTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( ! strcmp( tmpEvent->ifo, ifo ) )
    {
      /* ifos match so keep this event */
      if ( ! eventList  )
      {
        eventList = tmpEvent;
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
      LALFreeSnglRingdown ( status->statusPtr, &tmpEvent );
    }
  }
  *eventHead = eventList; 


  DETATCHSTATUSPTR (status);
  RETURN (status);
}  


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALIfoCountSingleRingdown(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglRingdownTable          *input,
    InterferometerNumber        ifoNumber 
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCountSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !(*numTrigs), status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_ringdown tables and return a
     pointer to the head of a linked list of tables for a specific IFO */
  for( thisEvent = input; thisEvent; thisEvent = thisEvent->next )
  {
    if ( ifoNumber == XLALIFONumber(thisEvent->ifo) )
    {
      /* IFOs match so count this trigger */
      ++(*numTrigs);
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}  



