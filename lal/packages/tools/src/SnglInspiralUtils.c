/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., Fairhurst, S. and Messaritaki, E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglInspiralUtilsCV">
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
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

NRCSID( SNGLINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglInspiralUtils.c}}

Provides a set of utilities for manipulating \texttt{snglInspiralTable}s.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglInspiralUtilsCP}
\idx{LALSortSnglInspiral()}
\idx{LALCompareSnglInspiralByMass()}
\idx{LALCompareSnglInspiralByPsi()}
\idx{LALCompareSnglInspiralByTime()}
\idx{LALCompareSnglInspiral()}
\idx{LALCompareInspirals()}
\idx{LALClusterSnglInspiralTable()}
\idx{LALTimeCutSingleInspiral()}
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

static INT8 end_time(const SnglInspiralTable *x)
{
	return(XLALGPStoINT8(&x->end_time));
}

static INT4 end_time_sec(const SnglInspiralTable *x)
{
	return(x->end_time.gpsSeconds);
}

static INT4 end_time_nsec(const SnglInspiralTable *x)
{
	return(x->end_time.gpsNanoSeconds);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALFreeSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;

  INITSTATUS( status, "LALFreeSnglInspiral", SNGLINSPIRALUTILSC );
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

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALFreeSnglInspiral (
    SnglInspiralTable **eventHead
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

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALSortSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;

  INITSTATUS( status, "LALSortSnglInspiral", SNGLINSPIRALUTILSC );

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
  eventHandle = (SnglInspiralTable **) 
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );
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


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByMass (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const SnglInspiralTable *aPtr = *((const SnglInspiralTable * const *)a);
  const SnglInspiralTable *bPtr = *((const SnglInspiralTable * const *)b);

  if ( aPtr->mass1 > bPtr->mass1 )
  {
    return 1;
  }
  else if ( aPtr->mass1 < bPtr->mass1 )
  {
    return -1;
  }
  else if ( aPtr->mass2 > bPtr->mass2 )
  {
    return 1;
  }
  else if ( aPtr->mass2 < bPtr->mass2 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByPsi (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const SnglInspiralTable *aPtr = *((const SnglInspiralTable * const *)a);
  const SnglInspiralTable *bPtr = *((const SnglInspiralTable * const *)b);

  if ( aPtr->psi0 > bPtr->psi0 )
  {
    return 1;
  }
  else if ( aPtr->psi0 < bPtr->psi0 )
  {
    return -1;
  }
  else if ( aPtr->psi3 > bPtr->psi3 )
  {
    return 1;
  }
  else if ( aPtr->psi3 < bPtr->psi3 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  const SnglInspiralTable *aPtr = *((const SnglInspiralTable * const *)a);
  const SnglInspiralTable *bPtr = *((const SnglInspiralTable * const *)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->end_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->end_time) );

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


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALCompareSnglInspiral (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    SnglInspiralAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta, tb;
  REAL4 dm1, dm2;
  REAL4 dmchirp, deta;
  REAL4 dpsi0, dpsi3;

  INITSTATUS( status, "LALCompareSnglInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  params->match = 1;

  /* check that triggers come from different IFOs */
  if( strcmp(aPtr->ifo, bPtr->ifo) )
  {
    LALInfo( status, "Triggers from different IFOs");
    params->match = 1;
  }
  else
  {
    LALInfo( status, "Triggers from same IFO");
    params->match = 0;
  }

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

  /* compare on trigger time coincidence */
  if ( labs( ta - tb ) < params->dt && params->match)
  {
    LALInfo( status, "Triggers pass time coincidence test");
    params->match = 1;
  }
  else if ( labs( ta - tb ) < params->dt && params->match)
  {
    LALInfo( status, "Triggers fail time coincidence test" );
    params->match = 0;
  }

  /* perform the mass parameter test */
  if( params->match )
  {
    /* compare psi0 and psi3 parameters */
    if ( params->test == psi0_and_psi3 )
    {
      dpsi0 = fabs( aPtr->psi0 - bPtr->psi0 );
      dpsi3 = fabs( aPtr->psi3 - bPtr->psi3 );

      if ( dpsi0 <= params->dpsi0 && dpsi3 <= params->dpsi3 )
      {
        LALInfo( status, "Triggers are coincident in psi0 and psi3" );
        params->match = 1;
      }
      else
      {
        LALInfo( status, "Triggers are not coincident in psi0 and psi3" );
        params->match = 0;
      }
    }
    else if ( params->test == m1_and_m2 )
    {  
      dm1 = fabs( aPtr->mass1 - bPtr->mass1 );
      dm2 = fabs( aPtr->mass2 - bPtr->mass2 );

      /* compare mass1 and mass2 parameters */
      if ( dm1 <= params->dm && dm2 <= params->dm )
      {
        LALInfo( status, "Triggers are coincident in mass1 and mass2" );
        params->match = 1;
      }
      else
      {
        LALInfo( status, "Triggers are not coincident in mass1 and mass2" );
        params->match = 0;
      }
    }
    else if ( params->test == mchirp_and_eta )
    {  
      dmchirp = fabs( aPtr->mchirp - bPtr->mchirp );
      deta = fabs( aPtr->eta - bPtr->eta );

      /* compare mchirp and eta parameters */
      if ( dmchirp <= params->dmchirp && deta <= params->deta )
      {
        LALInfo( status, "Triggers are coincident in mchirp and eta" );
        params->match = 1;
      }
      else
      {
        LALInfo( status, "Triggers fail mchirp, eta coincidence test" );
        params->match = 0;
      }
    }
    else
    {
      LALInfo( status, "error: unknown test\n" );
      params->match = 0;
    }
  }

  /* check for distance consistency */
  if ( params->match && params->test == m1_and_m2 )
  {
    if ( fabs( (aPtr->eff_distance - bPtr->eff_distance) / aPtr->eff_distance) 
        < params->epsilon / bPtr->snr + params->kappa )
    {
      LALInfo( status, "Triggers are coincident in eff_distance" );
      params->match = 1;
    }
    else
    {
      LALInfo( status, "Triggers fail eff_distance coincidence test" );
      params->match = 0;
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALCompareInspirals (
    LALStatus                *status,
    SnglInspiralTable        *aPtr,
    SnglInspiralTable        *bPtr,
    InspiralAccuracyList     *params
    )
/* </lalVerbatim> */
{
  INT8    ta,  tb;
  REAL4   dmass1, dmass2;
  REAL4   dmchirp, deta;
  REAL4   dpsi0, dpsi3;
  InterferometerNumber ifoaNum,  ifobNum;
  SnglInspiralAccuracy aAcc, bAcc;

  INITSTATUS( status, "LALCompareInspirals", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  params->match = 1;

  /* check that triggers come from different IFOs */
  if( strcmp(aPtr->ifo, bPtr->ifo) )
  {
    LALInfo( status, "Triggers from different IFOs");
    params->match = 1;
  }
  else
  {
    LALInfo( status, "Triggers from same IFO");
    params->match = 0;
    goto exit;
  }

  ifoaNum = XLALIFONumber( aPtr->ifo );
  ifobNum = XLALIFONumber( bPtr->ifo );
  
  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

  /* compare on trigger time coincidence */
  aAcc = params->ifoAccuracy[ifoaNum];
  bAcc = params->ifoAccuracy[ifobNum];

  

  /* XXX Need to add light travel time between sites XXX */
  if ( labs( ta - tb ) < (aAcc.dt + bAcc.dt)
      + params->lightTravelTime[ifoaNum][ifobNum])
  {
    LALInfo( status, "Triggers pass time coincidence test");
    params->match = 1;
  }
  else 
  {
    LALInfo( status, "Triggers fail time coincidence test" );
    params->match = 0;
    goto exit;
  }

  /* compare psi0 and psi3 parameters */
  if ( params->test == psi0_and_psi3 )
  {
    dpsi0 = fabs( aPtr->psi0 - bPtr->psi0 );
    dpsi3 = fabs( aPtr->psi3 - bPtr->psi3 );
    
    if ( ( dpsi0 <= (aAcc.dpsi0 + bAcc.dpsi0) )
        && ( dpsi3 <= (aAcc.dpsi3 + bAcc.dpsi3) ))
    {
      LALInfo( status, "Triggers are coincident in psi0 and psi3" );
      params->match = 1;
    }
    else
    {
      LALInfo( status, "Triggers are not coincident in psi0 and psi3" );
      params->match = 0;
      goto exit;
    }
  }
  else if ( params->test == m1_and_m2 )
  {  
    dmass1 = fabs( aPtr->mass1 - bPtr->mass1 );
    dmass2 = fabs( aPtr->mass2 - bPtr->mass2 );

    /* compare mass1 and mass2 parameters */
    if ( (dmass1 <= (aAcc.dm + bAcc.dm) )
      && (dmass2 <= (aAcc.dm + bAcc.dm) ))
    {
      LALInfo( status, "Triggers are coincident in mass1 and mass2" );
      params->match = 1;
    }
    else
    {
      LALInfo( status, "Triggers are not coincident in mass1 and mass2" );
      params->match = 0;
      goto exit;
    }
  }
  else if ( params->test == mchirp_and_eta )
  {  
    dmchirp = fabs( aPtr->mchirp - bPtr->mchirp );
    deta = fabs( aPtr->eta - bPtr->eta );

    /* compare mchirp and eta parameters */
    if ( (dmchirp <= (aAcc.dmchirp + bAcc.dmchirp))
          && (deta <= (aAcc.deta + bAcc.deta)) )
    {
      LALInfo( status, "Triggers are coincident in mchirp and eta" );
      params->match = 1;
    }
    else
    {
      LALInfo( status, "Triggers fail mchirp, eta coincidence test" );
      params->match = 0;
      goto exit;
    }
  }
  else
  {
    LALInfo( status, "error: unknown test\n" );
    params->match = 0;
    goto exit;
  }

exit:
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALClusterSnglInspiralTable (
    LALStatus                  *status,
    SnglInspiralTable          *inspiralEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  SnglInspiralTable     *thisEvent=NULL;
  SnglInspiralTable     *prevEvent=NULL;

  INITSTATUS( status, "LALClusterSnglInspiralTable", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( inspiralEvent, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  thisEvent = inspiralEvent->next;
  prevEvent = inspiralEvent;

  while ( thisEvent )
  {
    INT8 currTime;
    INT8 prevTime;

    /* compute the time in nanosec for each event trigger */
    LALGPStoINT8(status->statusPtr, &currTime, &(thisEvent->end_time));
    CHECKSTATUSPTR(status);

    LALGPStoINT8(status->statusPtr, &prevTime, &(prevEvent->end_time));
    CHECKSTATUSPTR(status);

    /* find events within the cluster window */
    if ( (currTime - prevTime) < dtimeNS )
    {
      /* displace previous event in cluster */
      if ( 
          (
           (clusterchoice == snr_and_chisq) && 
           (thisEvent->snr > prevEvent->snr) && 
           (thisEvent->chisq < prevEvent->chisq)
          ) || (
            (clusterchoice == snrsq_over_chisq) &&
            (thisEvent->snr)*(thisEvent->snr)/(thisEvent->chisq) > 
            (prevEvent->snr)*(prevEvent->snr)/(prevEvent->chisq)
            ) || (
              (clusterchoice == snr) && (thisEvent->snr > prevEvent->snr)
              )
         )
      {
        memcpy( prevEvent, thisEvent, sizeof(SnglInspiralTable) );
        thisEvent->event_id = NULL;
      }

      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFreeSnglInspiral ( status->statusPtr, &thisEvent );
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



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALTimeCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;


  INITSTATUS( status, "LALTimeCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );




  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->end_time.gpsSeconds >= startTime->gpsSeconds &&
        tmpEvent->end_time.gpsSeconds < endTime->gpsSeconds )
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
      LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
    }
  }
  *eventHead = inspiralEventList; 

  DETATCHSTATUSPTR (status);
  RETURN (status);

}  

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALTimeCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    INT8                        startTimeNS,
    INT8                        endTimeNS
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
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
      XLALFreeSnglInspiral ( &tmpEvent );
    }
  }
  *eventHead = inspiralEventList; 

  return (0);
}  


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALalphaFCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       alphaFcut
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;


  INITSTATUS( status, "LALalphaFCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );




  /* Remove all the triggers with alphaF > alphaFcut */

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( (tmpEvent->alpha * pow(tmpEvent->f_final,(2.0/3.0))) <= alphaFcut )
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
      LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
    }
  }
  *eventHead = inspiralEventList; 

  DETATCHSTATUSPTR (status);
  RETURN (status);

}  



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALIfoCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    CHAR                       *ifo
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *eventList = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  SnglInspiralTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoScanSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that eventHead is non-null */
  ASSERT( eventHead, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent = *eventHead;
  
  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
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
      LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
    }
  }
  *eventHead = eventList; 


  DETATCHSTATUSPTR (status);
  RETURN (status);
}  


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALIfoCountSingleInspiral(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglInspiralTable          *input,
    InterferometerNumber        ifoNumber 
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCountSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !(*numTrigs), status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_inspiral tables and return a
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

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALTimeSlideSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable          *triggerList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime,
    LIGOTimeGPS                 slideTimes[LAL_NUM_IFO]
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent   = NULL;
  INT8                  startTimeNS = 0;
  INT8                  endTimeNS   = 0;
  INT8                  slideNS     = 0;
  INT8                  trigTimeNS  = 0;
  INITSTATUS( status, "LALTimeSlideSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* time slide triggers by a time = slideTime, except those from the
   * instrument skipIfo which are left untouched. If you want to slide 
   * all triggers, simply set skipIfo = LAL_UNKNOWN_IFO */

  
  /* check that input non-null */
  ASSERT( triggerList, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  
  if ( startTime )
  {
    LALGPStoINT8( status->statusPtr, &startTimeNS, startTime );
  }

  if ( endTime )
  {
    LALGPStoINT8( status->statusPtr, &endTimeNS, endTime );
  }
  
  for( thisEvent = triggerList; thisEvent; thisEvent = thisEvent->next )
  {
    /* calculate the slide time in nanoseconds */
    LALGPStoINT8( status->statusPtr, &slideNS, 
        &(slideTimes[XLALIFONumber(thisEvent->ifo)]) );
    /* and trig time in nanoseconds */
    LALGPStoINT8( status->statusPtr, &trigTimeNS, &(thisEvent->end_time));
    trigTimeNS += slideNS;
    
    if ( startTimeNS && trigTimeNS < startTimeNS )
    {
      /* if before startTime, then wrap trigger time */
      trigTimeNS += endTimeNS - startTimeNS;
    }
    else if ( endTimeNS && trigTimeNS > endTimeNS )
    {
      /* if after endTime, then wrap trigger time */
      trigTimeNS -= endTimeNS - startTimeNS;
    }
   
    /* convert back to LIGOTimeGPS */
    LALINT8toGPS( status->statusPtr, &(thisEvent->end_time), &trigTimeNS );
  }         

  DETATCHSTATUSPTR (status);
  RETURN (status);
}  




/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALPlayTestSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  INITSTATUS( status, "LALPlayTestSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = *eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      SnglInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      LALGPStoINT8( status->statusPtr, &triggerTime, &(tmpEvent->end_time) );
      LALINT8NanoSecIsPlayground( status->statusPtr, &isPlay, &triggerTime );

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
        LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
      }
    }
    *eventHead = inspiralEventList; 
    if ( *dataType == playground_only )
    {
      /*LALInfo( status, "Kept %d playground triggers \n", numTriggers );*/
    }
    else if ( *dataType == exclude_play )
    {
      /*LALInfo( status, "Kept %d non-playground triggers \n", numTriggers );*/
    }
  }
  else if ( *dataType == all_data )
  {
    LALInfo( status, "Keeping all triggers since all_data specified\n" );
  }
  else
  {
    LALInfo( status, "Unknown data type, returning no triggers\n" );
    *eventHead = NULL;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}  


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALCreateTrigBank(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    SnglInspiralParameterTest  *test
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *trigBankList = NULL;
  SnglInspiralTable   **eventHandle = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  INT4 numTriggers = 0;
  INT4 numEvents = 0;
  INT4 i = 0;

  INITSTATUS( status, "LALCreateTrigBank", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  /* count the number of events */
  for ( thisEvent = *eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }

  eventHandle = (SnglInspiralTable **) 
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );

  for ( i = 0, thisEvent = *eventHead; i < numEvents; 
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  if ( *test == m1_and_m2 )
  {     
    LALInfo( status, "sorting events by mass... " ); 
    qsort( eventHandle, numEvents, sizeof(eventHandle[0]), 
        LALCompareSnglInspiralByMass );
    LALInfo( status, "done\n" ); 
  }
  else if ( *test == psi0_and_psi3 )
  { 
    LALInfo( status, "sorting events by psi... " );
    qsort( eventHandle, numEvents, sizeof(eventHandle[0]),
        LALCompareSnglInspiralByPsi );
    LALInfo( status, "done\n" );
  }
  else
  {
    ABORT( status, LIGOMETADATAUTILSH_ETEST, LIGOMETADATAUTILSH_MSGETEST );
  }

  /* create a linked list of sorted templates */
  LALInfo( status, "discarding template with duplicate masses: " );

  numTriggers = 0;
  trigBankList = prevEvent = eventHandle[0];
  if ( trigBankList ) numTriggers = 1;

  for ( i = 1; i < numEvents; ++i )
  {
    if ( *test == m1_and_m2 )
    {
      if ( (prevEvent->mass1 == eventHandle[i]->mass1)  &&
          (prevEvent->mass2 == eventHandle[i]->mass2) ) 
      {
        /* discard the event as it is a duplicate */
        LALFreeSnglInspiral( status->statusPtr, &(eventHandle[i]) );
        LALInfo( status, "-" );
      }
      else
      {
        /* add the event to the linked list */
        prevEvent = prevEvent->next = eventHandle[i];
        LALInfo( status, "+" );
      }
    }
    else if ( *test == psi0_and_psi3 )
    {
      if ( (prevEvent->psi0 == eventHandle[i]->psi0)  &&
          (prevEvent->psi3 == eventHandle[i]->psi3) )
      {
        /* discard the event as it is a duplicate */
        LALFreeSnglInspiral( status->statusPtr, &(eventHandle[i]) );
        LALInfo( status, "-" );
      }
      else
      {
        /* add the event to the linked list */
        prevEvent = prevEvent->next = eventHandle[i];
        LALInfo( status, "+" );
      }
    }
    else
    {
      ABORT( status, LIGOMETADATAUTILSH_ETEST, LIGOMETADATAUTILSH_MSGETEST );
    }
  }

  /* if the list is non-emnpty, make sure it is terminated */
  if ( prevEvent ) prevEvent->next = NULL;

  LALFree( eventHandle );

  /* return the head of the linked list in eventHead */

  *eventHead = trigBankList; 

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALIncaCoincidenceTest(
    LALStatus                  *status,
    SnglInspiralTable         **ifoAOutput,
    SnglInspiralTable         **ifoBOutput,
    SnglInspiralTable          *ifoAInput,
    SnglInspiralTable          *ifoBInput,
    SnglInspiralAccuracy       *errorParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *currentTrigger[2];
  SnglInspiralTable    *coincidentEvents[2];
  SnglInspiralTable    *outEvent[2];
  SnglInspiralTable    *currentEvent;

  INT8 ta,tb;
  INT4 j;

  INITSTATUS( status, "LALIncaCoincidenceTest", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  memset( currentTrigger, 0, 2 * sizeof(SnglInspiralTable *) );
  memset( coincidentEvents, 0, 2 * sizeof(SnglInspiralTable *) );
  memset( outEvent, 0, 2 * sizeof(SnglInspiralTable *) );


  if ( ! ifoAInput )
  {
    LALInfo( status, "No input triggers from IFO A, exiting");
  }

  if ( ! ifoBInput )
  {
    LALInfo( status, "No input triggers from IFO B, exiting");
  }

  currentTrigger[1] = ifoBInput;

  for( currentTrigger[0]=ifoAInput; currentTrigger[0]; 
      currentTrigger[0] = currentTrigger[0]->next  )
  {
    LALGPStoINT8( status->statusPtr, &ta, &(currentTrigger[0]->end_time) );

    /* spin ifo b until the current trigger is within the coinicdence */
    /* window of the current ifo a trigger                            */
    while ( currentTrigger[1] )
    {
      LALGPStoINT8( status->statusPtr, &tb, &(currentTrigger[1]->end_time) );

      if ( tb > ta - errorParams->dt )
      {
        /* we have reached the time coinicidence window */
        break;
      }
      currentTrigger[1] = currentTrigger[1]->next;
    }

    /* look for coincident events in B within the time window */
    currentEvent = currentTrigger[1];

    while ( currentTrigger[1] )
    {
      LALGPStoINT8( status->statusPtr, &tb, &(currentTrigger[1]->end_time) );

      if (tb > ta + errorParams->dt )
      {
        /* we are outside the time coincidence so move to next event */
        break;
      }
      else
      {
        /* call the LAL function which compares events parameters */
        LALCompareSnglInspiral( status->statusPtr, currentTrigger[0],
            currentTrigger[1], errorParams );
      }

      if ( errorParams->match )
      {
        /* store this event for output */
        LALInfo( status, "    >>> found coincidence <<<" );

        for ( j = 0; j < 2; ++j )
        {
          if ( ! coincidentEvents[j] )
          {
            coincidentEvents[j] = outEvent[j] = (SnglInspiralTable *) 
              LALCalloc( 1, sizeof(SnglInspiralTable) );
          }
          else
          {
            outEvent[j] = outEvent[j]->next = (SnglInspiralTable *) 
              LALCalloc( 1, sizeof(SnglInspiralTable) );
          }

          memcpy( outEvent[j], currentTrigger[j], sizeof(SnglInspiralTable) );
          outEvent[j]->next = NULL;
        }  
      }

      currentTrigger[1] = currentTrigger[1]->next;

    } /* end loop over current events */

    /* go back to saved current IFO B trigger */
    currentTrigger[1] = currentEvent;

  } /* end loop over ifo A events */

  *ifoAOutput = coincidentEvents[0];
  *ifoBOutput = coincidentEvents[1];

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALTamaCoincidenceTest(
    LALStatus                  *status,
    SnglInspiralTable         **ifoAOutput,
    SnglInspiralTable         **ifoBOutput,
    SnglInspiralTable          *ifoAInput,
    SnglInspiralTable          *ifoBInput,
    SnglInspiralAccuracy       *errorParams,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *currentTrigger[2];
  SnglInspiralTable    *coincidentEvents[2];
  SnglInspiralTable    *outEvent[2];
  SnglInspiralTable    *currentEvent = NULL;
  SnglInspiralTable    *timeCoincHead = NULL;
  SnglInspiralTable    *thisTimeCoinc = NULL;

  INT8 ta,tb;
  INT4 j;

  INITSTATUS( status, "LALIncaCoincidenceTest", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  memset( currentTrigger, 0, 2 * sizeof(SnglInspiralTable *) );
  memset( coincidentEvents, 0, 2 * sizeof(SnglInspiralTable *) );
  memset( outEvent, 0, 2 * sizeof(SnglInspiralTable *) );

  if ( ! ifoAInput )
  {
    LALInfo( status, "No input triggers from IFO A, exiting");
  }

  if ( ! ifoBInput )
  {
    LALInfo( status, "No input triggers from IFO B, exiting");
  }

  currentTrigger[1] = ifoBInput;

  for( currentTrigger[0]=ifoAInput; currentTrigger[0]; 
      currentTrigger[0] = currentTrigger[0]->next  )
  {
    LALGPStoINT8( status->statusPtr, &ta, &(currentTrigger[0]->end_time) );

    LALInfo( status, printf("  using IFO A trigger at %d + %10.10f\n",
          currentTrigger[0]->end_time.gpsSeconds, 
          ((REAL4) currentTrigger[0]->end_time.gpsNanoSeconds * 1e-9) ));

    /* spin ifo b until the current trigger is within the coinicdence */
    /* window of the current ifo a trigger                            */
    while ( currentTrigger[1] )
    {
      LALGPStoINT8( status->statusPtr, &tb, &(currentTrigger[1]->end_time) );

      if ( tb > ta - errorParams->dt )
      {
        /* we have reached the time coinicidence window */
        break;
      }
      currentTrigger[1] = currentTrigger[1]->next;
    }


    /* look for coincident events in B within the time window */
    currentEvent = currentTrigger[1];

    while ( currentTrigger[1] )
    {
      LALGPStoINT8( status->statusPtr, &tb, &(currentTrigger[1]->end_time) );

      if (tb > ta + errorParams->dt )
      {
        /* we are outside the time coincidence so move to next event */
        LALInfo( status, "outside the time coincidence window\n" );
        break;
      }
      else
      {
        /* store all time coincident triggers */
        if ( ! timeCoincHead )
        {
          timeCoincHead = thisTimeCoinc = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        else
        {
          thisTimeCoinc = thisTimeCoinc->next = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
        }

        memcpy( thisTimeCoinc, currentTrigger[1], 
            sizeof(SnglInspiralTable) );

        thisTimeCoinc->next = NULL;
      }
      currentTrigger[1] = currentTrigger[1]->next;


    }  /* end loop over current events */


    /* take the loudest time coincident trigger and compare other params */
    if ( timeCoincHead )
    {
      LALClusterSnglInspiralTable ( status->statusPtr, timeCoincHead, 
          2 * errorParams->dt, clusterchoice);

      currentTrigger[1] = timeCoincHead;


      /* call the LAL function which compares events parameters */
      LALCompareSnglInspiral( status->statusPtr, currentTrigger[0],
          currentTrigger[1], errorParams );

      if ( errorParams->match )
      {
        /* store this event for output */
        LALInfo( status, "    >>> found coincidence <<<\n" );

        for ( j = 0; j < 2; ++j )
        {
          if ( ! coincidentEvents[j] )
          {
            coincidentEvents[j] = outEvent[j] = (SnglInspiralTable *) 
              LALCalloc( 1, sizeof(SnglInspiralTable) );
          }
          else
          {
            outEvent[j] = outEvent[j]->next = (SnglInspiralTable *) 
              LALCalloc( 1, sizeof(SnglInspiralTable) );
          }

          memcpy( outEvent[j], currentTrigger[j], sizeof(SnglInspiralTable) );
          outEvent[j]->next = NULL;
        }  
      }

      /* reset the list of time coincident triggers to null */
      LALFreeSnglInspiral( status->statusPtr, &timeCoincHead );
      timeCoincHead = NULL;
    }
    /* go back to saved current IFO B trigger */
    currentTrigger[1] = currentEvent;

  } /* end loop over ifo A events */

  *ifoAOutput = coincidentEvents[0];
  *ifoBOutput = coincidentEvents[1];

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALMaxSnglInspiralOverIntervals(
    SnglInspiralTable         **eventHead,
    INT4                       deltaT
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *nextEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  /* if there are no events, then no-op */
  if ( ! *eventHead )
    return (0);

  inspiralEventList = *eventHead;
  thisEvent = *eventHead;
  nextEvent = thisEvent->next;

  while ( nextEvent )
  {
    if ( end_time_sec(nextEvent) == end_time_sec(thisEvent) &&
        end_time_nsec(nextEvent)/deltaT == end_time_nsec(thisEvent)/deltaT )
    {
      if ( nextEvent->snr > thisEvent->snr )
      {
        /* replace thisEvent with nextEvent */
        XLALFreeSnglInspiral ( &thisEvent );

        /* deal with start of the list */
        if (prevEvent)
          prevEvent->next = nextEvent;
        else
          inspiralEventList = nextEvent;

        /* standard stuff */
        thisEvent = nextEvent;
        nextEvent = thisEvent->next;
      }
      else
      {
        /* get rid of nextEvent */
        thisEvent->next = nextEvent->next;
        XLALFreeSnglInspiral ( &nextEvent );
        nextEvent = thisEvent->next;
      }
    }
    else
    {
      /* step to next set of events */
      prevEvent=thisEvent;
      thisEvent=nextEvent;
      nextEvent = thisEvent->next;
    }
  }

  *eventHead = inspiralEventList; 

  return (0);
}  

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
INT4 XLALCountSnglInspiral( SnglInspiralTable *head )
/* </lalVerbatim> */
{
  INT4 length;
  SnglInspiralTable *event;
  
  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}


