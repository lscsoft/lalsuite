/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglInspiralUtilsCV">
Author: Brown, D. A.
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

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglInspiralUtilsCP}
\idx{LALSortSnglInspiral()}
\idx{LALCompareSnglInspiralByMass()}
\idx{LALCompareSnglInspiralByPsi()}
\idx{LALCompareSnglInspiralByTime()}
\idx{LALCompareSnglInspiral()}
\idx{LALClusterSnglInspiralTable()}
\idx{LALTimeCutSingleInspiral()}
\idx{LALIfoScanSingleInspiral()}
\idx{LALPlayTestSingleInspiral()}
\idx{LALCreateTrigBank()}
\idx{LALIncaCoincidenceTest()}
\idx{LALTamaCoincidenceTest()}
                                                                                

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SnglInspiralUtilsCV}}

</lalLaTeX>
#endif

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
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

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
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

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
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);
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

  params->match = 0;

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

  /* compare on trigger time coincidence */
  if ( labs( ta - tb ) < params->dt )
  {
    if ( params->test == psi0_and_psi3 )
    {
      dpsi0 = fabs( aPtr->psi0 - bPtr->psi0 );
      dpsi3 = fabs( aPtr->psi3 - bPtr->psi3 );

      if ( dpsi0 <= params->dpsi0 && dpsi3 <= params->dpsi3 )
      {
        params->match = 1;
        LALInfo( status, "Triggers are coincident in psi0 and psi3" );
      }
      else
      {
        LALInfo( status, "Triggers are not coincident in psi0 and psi3" );
      }
    }
    else if ( params->test == m1_and_m2 )
    {  
      dm1 = fabs( aPtr->mass1 - bPtr->mass1 );
      dm2 = fabs( aPtr->mass2 - bPtr->mass2 );

      /* compare mass1 and mass2 parameters */
      if ( dm1 <= params->dm && dm2 <= params->dm )
      {
        if ( fabs( 
              (aPtr->eff_distance - bPtr->eff_distance) / aPtr->eff_distance
              ) < params->epsilon / bPtr->snr + params->kappa )
        {
          params->match = 1;
          LALInfo( status, "Triggers are coincident in eff_distance" );
        }
        else
        {
          LALInfo( status, "Triggers fail eff_distance coincidence test" );
        }
      }
      else
      {
        LALInfo( status, "Triggers fail mass coincidence test" );
      }
    }
    else if ( params->test == mchirp_and_eta )
    {  
      dmchirp = fabs( aPtr->mchirp - bPtr->mchirp );
      deta = fabs( aPtr->eta - bPtr->eta );

      /* compare mchirp and eta parameters */
      if ( dmchirp <= params->dmchirp && deta <= params->deta )
      {
        params->match = 1;
        LALInfo( status, "Triggers are coincident in mchirp and eta" );
      }
      else
      {
        LALInfo( status, "Triggers fail mchirp, eta coincidence test" );
      }
    }
    else
    {
      LALInfo( status, "error: unknown test\n" );
    }
  }
  else
  {
    LALInfo( status, "Triggers fail time coincidence test" );
  }

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
      }

      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFree( thisEvent );
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
      LALFree( tmpEvent );
    }
  }
  *eventHead = inspiralEventList; 

  DETATCHSTATUSPTR (status);
  RETURN (status);

}  


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALIfoScanSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **output,
    SnglInspiralTable          *input,
    CHAR                       *ifo
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !output, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  for( thisEvent = input; thisEvent; thisEvent = thisEvent->next );
  {
    SnglInspiralTable *keptEvent = NULL;

    if ( !strcmp(thisEvent->ifo, ifo) ) 
    {
      /* IFOs match so write this entry to the output table */
      if ( ! output  )
      {
        *output = keptEvent = (SnglInspiralTable *) 
          LALMalloc( sizeof(SnglInspiralTable) );
      }
      else
      {
        keptEvent = keptEvent->next = (SnglInspiralTable *) 
          LALMalloc( sizeof(SnglInspiralTable) );
      }
      keptEvent = thisEvent;
      keptEvent->next = NULL;
    }
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
        LALFree( tmpEvent );
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
        LALFree( eventHandle[i] );
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
        LALFree( eventHandle[i] );
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
              LALMalloc( sizeof(SnglInspiralTable) );
          }
          else
          {
            outEvent[j] = outEvent[j]->next = (SnglInspiralTable *) 
              LALMalloc( sizeof(SnglInspiralTable) );
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
  SnglInspiralTable    *currentEventHead, *currentEvent;

  INT8 ta,tb;
  INT4 j;

  INITSTATUS( status, "LALIncaCoincidenceTest", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  memset( currentTrigger, 0, 2 * sizeof(SnglInspiralTable *) );

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
    currentEventHead = currentEvent = currentTrigger[1];

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
        currentEvent = currentEvent->next = currentTrigger[1];
      }
    }
    currentEvent->next = NULL;

    LALClusterSnglInspiralTable ( status->statusPtr, currentEventHead, 
        2 * errorParams->dt, clusterchoice);

    currentTrigger[1] = currentEventHead;

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
            LALMalloc( sizeof(SnglInspiralTable) );
        }
        else
        {
          outEvent[j] = outEvent[j]->next = (SnglInspiralTable *) 
            LALMalloc( sizeof(SnglInspiralTable) );
        }

        memcpy( outEvent[j], currentTrigger[j], sizeof(SnglInspiralTable) );
        outEvent[j]->next = NULL;
      }  
    }

    /* go back to saved current IFO B trigger */
    currentTrigger[1] = currentEvent;

  } /* end loop over ifo A events */

  *ifoAOutput = coincidentEvents[0];
  *ifoBOutput = coincidentEvents[1];

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

