/*----------------------------------------------------------------------- 
 * 
 * File Name: CoincInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., and Fairhurst, S
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="CoincInspiralUtilsCV">
Author: Fairhurst, S.
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

NRCSID( COINCINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{CoincInspiralUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CoincInspiralUtilsCP}
\idx{LALCreateTwoIFOCoincList()}
\idx{LALAddSnglInspiralToCoinc()}
\idx{LALSnglInspiralCoincTest()}


\subsubsection*{Description}

\texttt{LALCreateTwoIFOCoincList()} takes in a linked list of single inspiral
tables and returns a list of two instrument coincidences.  The coincidence
requirements are given by the \texttt{accuracyParams}.  When single inspirals
from two different instruments are found to be coincident, the code creates a
new \texttt{coincInspiralTable} and uses \texttt{LALAddSnglInspiralToCoinc()}
to add the single inspirals to the coinc.  The function returns
\texttt{coincOutput} which is a pointer to the head of a linked list of
\texttt{CoincInspiralTable}s.

\texttt{LALCreateNIFOCoincList()} takes linked list of
\texttt{CoincInspiralTable}s, assumed to contain (N-1) ifo coincidences and
creates all N ifo coincidences.  Both the input and output list of
\texttt{CoincInspiralTable}s are passed as \texttt{coincHead}.  

\texttt{LALAddSnglInspiralToCoinc()} adds a pointer to a single inspiral table
to a coinc inspiral table.  Upon entry, if \texttt{coincPtr} points to a
\texttt{NULL} coinc inspiral table, the table is created before a pointer to
the single inspiral table is added.  Additionally, an \texttt{eventId} table is
created for the single inspiral table.  This points to both the single and
coinc inspirals.  If an \texttt{eventId} already exists for the single
inspiral, another eventId table is added to the linked list.  The linked list
of \texttt{eventId}s associated to a single inspiral table allow us to easily
determine which coincident events each single is a part of.

\texttt{LALSnglInspiralCoincTest()} tests for coincidence between a single
inspiral and a coinc inspiral.  It works by testing for coincidence between
each non-null entry in the coinc inspiral and the single.  This is done using
\texttt{LALCompareSnglInspiral()}.  If all members of the coinc are found to be
coincident with the single, the \texttt{accuracyParams.match} is set to 1,
otherwise to 0.

\texttt{LALExtractSnglInspiralFromCoinc()} extracts the information from a
linked list of \texttt{coincInspiralTable}s and returns it as a linked list of
\texttt{snglInspiralTable}s.  Thus, the output \texttt{snglPtr} is a pointer to
a linked list of single inspiral tables.  That list contains only single
inspirals which are found in coincidence.  In order to preserve the coincidence
information, we assign to each coincident event an integer value.  This is
stored in the \texttt{UINT4 id} field of the \texttt{eventIDColumn} of each
single inspiral which forms part of the coincidence.  The \texttt{id} is set
equal to $10^{9} \times$ \texttt{gpsStartTime} $+ 10^{5} \times$
\texttt{slideNum} $+$ event number. We do not assign multiple \texttt{id}
values to a given single inspiral table, but instead make multiple copies of
the table, each with a unique \texttt{id}.  

\texttt{LALCoincCutSnglInspiral()} extracts all single inspirals from a
specific ifo which are in coinc inspirals.  The output \texttt{snglPtr} is a
pointer to a linked list of single inspiral tables.  That list contains only
single inspirals from the specified \texttt{ifo} which are found in
coincidence. 


\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{CoincInspiralUtilsCV}}

</lalLaTeX>
#endif

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALCreateTwoIFOCoincList(
    LALStatus                  *status,
    CoincInspiralTable        **coincOutput,
    SnglInspiralTable          *snglInput,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable            *currentTrigger[2];
  INT8                          currentTriggerNS[2];
  CoincInspiralTable           *coincHead = NULL;
  CoincInspiralTable           *thisCoinc = NULL;
  INT4                          numEvents = 0;
  INT4                          ifoNumber;
  INT8                          maxTimeDiff = 0;

  INITSTATUS( status, "LALCreateTwoIFOCoincList", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincOutput, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );
  memset( currentTrigger, 0, 2 * sizeof(SnglInspiralTable *) );

  
  /* calculate the maximum time delay 
   * set it equal to 2 * worst IFO timing accuracy plus
   * light travel time for earth's diameter 
   * (detectors can't be further apart than this) */
  
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    maxTimeDiff = (maxTimeDiff > accuracyParams->ifoAccuracy[ifoNumber].dt) ?
      maxTimeDiff : accuracyParams->ifoAccuracy[ifoNumber].dt;
  }
  maxTimeDiff *= 2;    
  maxTimeDiff += (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );
  
  for ( currentTrigger[0] = snglInput; currentTrigger[0]->next;
      currentTrigger[0] = currentTrigger[0]->next)
  {

    /* calculate the time of the trigger */
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[0], 
        &(currentTrigger[0]->end_time) );

    /* set next trigger for comparison */
    currentTrigger[1] = currentTrigger[0]->next;
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
          &(currentTrigger[1]->end_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {
      /* check that triggers pass coincidence test */
      LALCompareInspirals( status->statusPtr, currentTrigger[0], 
          currentTrigger[1], accuracyParams );

      /* test whether we have coincidence */
      if ( accuracyParams->match )
      {
        LALInfo( status, "Found double coincident trigger,");
        /* create a 2 IFO coinc and store */
        if ( ! coincHead  )
        {
          coincHead = thisCoinc = (CoincInspiralTable *) 
            LALCalloc( 1, sizeof(CoincInspiralTable) );
        }
        else
        {
          thisCoinc = thisCoinc->next = (CoincInspiralTable *) 
            LALCalloc( 1, sizeof(CoincInspiralTable) );
        }

        /* Add the two triggers to the coinc */
        LALAddSnglInspiralToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[0] );
        LALAddSnglInspiralToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[1] );

        ++numEvents;

      }

      /* scroll on to the next sngl inspiral */
      
      if ( (currentTrigger[1] = currentTrigger[1]->next) )
      {
        LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
            &(currentTrigger[1]->end_time) );
      }
      else
      {
        LALInfo(status, "Second trigger has reached end of list");
        break;
      }
    }
  }

  *coincOutput = coincHead;

  DETATCHSTATUSPTR (status);
  RETURN (status);
} 


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALCreateNIFOCoincList(
    LALStatus                  *status,
    CoincInspiralTable        **coincHead,
    InspiralAccuracyList       *accuracyParams,
    INT4                        N
    )
/* </lalVerbatim> */
{
  INT4                          numEvents  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          ifoNum     = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincInspiralTable           *thisCoinc     = NULL;
  CoincInspiralTable           *otherCoinc    = NULL;
  CoincInspiralTable           *nIfoCoincHead = NULL;
  CoincInspiralTable           *thisNIfoCoinc = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALCreateThreeIFOCoincList", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  for( thisCoinc = *coincHead; thisCoinc; thisCoinc = thisCoinc->next)
  {
    /* check that this is an (N-1) coinc */
    if ( thisCoinc->numIfos == N - 1 )
    {
      /* look up the first single inspiral */
      for ( firstEntry = 0; firstEntry < LAL_NUM_IFO; firstEntry++)
      {
        if ( thisCoinc->snglInspiral[firstEntry] )
        {
          LALInfo( status, "Found the first entry in the coinc" );
          break;
        }
      }

      /* get the list of event IDs for this first entry */
      eventIDHead = thisCoinc->snglInspiral[firstEntry]->event_id;

      /* loop over the (N-1) ifo coincs that first entry is a member of 
       * try to find an N ifo coinc */
      for( thisID = eventIDHead; thisID; thisID = thisID->next )
      {
        otherCoinc = thisID->coincInspiralTable;

        if( otherCoinc->numIfos == N - 1 )
        {
          /* loop over all singles which are alphabetically before the 
           * first one in thisCoinc */
          for( ifoNumber = 0; ifoNumber < firstEntry; ifoNumber++ )
          {
            /* test whether we have an N ifo coincidence */
            accuracyParams->match = 0;

            if ( otherCoinc->snglInspiral[ifoNumber] )
            {
              LALSnglInspiralCoincTest( status->statusPtr, thisCoinc, 
                  otherCoinc->snglInspiral[ifoNumber], accuracyParams );
            }

            if ( accuracyParams->match )
            {
              LALInfo( status, "We have found an N ifo coinc, storing");
              ++numEvents;
              
              /* create a N IFO coinc and store */
              if ( ! nIfoCoincHead  )
              {
                nIfoCoincHead = thisNIfoCoinc = (CoincInspiralTable *) 
                  LALCalloc( 1, sizeof(CoincInspiralTable) );
              }
              else
              {
                thisNIfoCoinc = thisNIfoCoinc->next = (CoincInspiralTable *) 
                  LALCalloc( 1, sizeof(CoincInspiralTable) );
              }

              /* add the single to the new N coinc */
              LALAddSnglInspiralToCoinc( status->statusPtr, &thisNIfoCoinc,
                  otherCoinc->snglInspiral[ifoNumber] );

              /* add the triggers from the (N-1) coinc to the new N coinc */
              for( ifoNum = 0; ifoNum < LAL_NUM_IFO; ifoNum++ )
              {
                if( thisCoinc->snglInspiral[ifoNum] )
                {
                  LALAddSnglInspiralToCoinc( status->statusPtr, &thisNIfoCoinc,
                      thisCoinc->snglInspiral[ifoNum] );
                }
              }
            } /* closes: if ( accuracyParams->match ) */
          }
        }
      } /* closes: for( thisID = eventIDHead; thisID; thisID->next ) */
    }
  } /* closes: for( thisCoinc = coincHead; thisCoinc; 
     *              thisCoinc = thisCoinc->next) */

  /* append the N ifo coincs to the end of the linked list */
  thisCoinc->next = nIfoCoincHead;


  DETATCHSTATUSPTR (status);
  RETURN (status);

}



/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALAddSnglInspiralToCoinc(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable          *snglInspiral
    )
/* </lalVerbatim> */
{
  CoincInspiralTable  *coincInspiral = NULL;
  EventIDColumn       *eventId = NULL;
  
  INITSTATUS( status, "LALAddSnglInspiralToCoinc", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglInspiral, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  coincInspiral = *coincPtr;
 
  /* allocate memory for new coinc if it doesn't exist */
  if (! coincInspiral )
  {
    coincInspiral = (CoincInspiralTable *) 
      LALCalloc( 1, sizeof(CoincInspiralTable) );
  }
  
  switch ( (snglInspiral->ifo)[0] ) 
  {
    case 'G':
      coincInspiral->snglInspiral[LAL_IFO_G1] = snglInspiral;
      break;

    case 'H':
      if ( !strcmp( snglInspiral->ifo, "H1" ) )
      {
        coincInspiral->snglInspiral[LAL_IFO_H1] = snglInspiral;
      }
      else if (!strcmp( snglInspiral->ifo, "H2" ) )
      {
        coincInspiral->snglInspiral[LAL_IFO_H2] = snglInspiral;
      }
      else
      {
        /* Invalid Hanford Detector */
        ABORT( status, LIGOMETADATAUTILSH_EDET, LIGOMETADATAUTILSH_MSGEDET );
      } 
      break;

    case 'L':
      coincInspiral->snglInspiral[LAL_IFO_L1] = snglInspiral;
      break;

    case 'T':
      coincInspiral->snglInspiral[LAL_IFO_T1] = snglInspiral;
      break;

    case 'V':
      coincInspiral->snglInspiral[LAL_IFO_V1] = snglInspiral;
      break;

    default:
      /* Invalid Detector Site */
      ABORT( status, LIGOMETADATAUTILSH_EDET, LIGOMETADATAUTILSH_MSGEDET );
  }

  ++(coincInspiral->numIfos);

  /* create an eventId for the single, populate it with the single and coinc */
  if ( ! snglInspiral->event_id )
  {
    eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    snglInspiral->event_id = eventId;
  }
  else
  {
     for( eventId = snglInspiral->event_id; eventId->next; 
         eventId = eventId->next);
     eventId = eventId->next = (EventIDColumn *) 
         LALCalloc( 1, sizeof(EventIDColumn) );
  }
  eventId->snglInspiralTable = snglInspiral;
  eventId->coincInspiralTable = coincInspiral;

  *coincPtr = coincInspiral;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}





/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALSnglInspiralCoincTest(
    LALStatus                  *status,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;
  
  INITSTATUS( status, "LALSnglInspiralCoincTest", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  /* Loop over sngl_inspirals contained in coinc_inspiral */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    thisCoincEntry = coincInspiral->snglInspiral[ifoNumber];

    if ( thisCoincEntry )
    {
      /* snglInspiral entry exists for this IFO, perform coincidence test */
      if ( ifoNumber == XLALIFONumber(snglInspiral->ifo) )
      {
        LALInfo( status, "We already have a coinc from this IFO" );
        accuracyParams->match = 0;
      }

      else
      {
        LALCompareInspirals ( status->statusPtr, snglInspiral, 
            thisCoincEntry, accuracyParams );
      }
      /* set match to zero if no match.  Keep same if match */
      match *= accuracyParams->match;
    }
  }
  /* returm errorParams->match to be 1 if we match, zero otherwise */
  accuracyParams->match = match;
  if ( accuracyParams->match == 0 ) 
    LALInfo( status, "Coincidence test failed" );
  if ( accuracyParams->match == 1 ) 
    LALInfo( status, "Coincidence test passed" );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALExtractSnglInspiralFromCoinc(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr,
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )
/* </lalVerbatim> */
{
  SnglInspiralTable  *snglHead = NULL;
  SnglInspiralTable  *thisSngl = NULL;
  SnglInspiralTable  *thisCoincEntry = NULL;
  CoincInspiralTable *thisCoinc = NULL;
  EventIDColumn      *eventId = NULL;
  UINT4               eventNum = 1;
  INT4                j;
  
  INITSTATUS( status, "LALExtractCoincSngls", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( snglPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! (*snglPtr), status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  
  if ( !coincInspiral )
  {
    LALInfo( status, "Empty coincInspiral passed to LALExtractCoincSngls");
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }

  /* loop over the linked list of coinc inspirals */
  for( thisCoinc = coincInspiral; thisCoinc; thisCoinc = thisCoinc->next,
      ++eventNum)
  {
    /* loop over the interferometers */
    for ( j = 0; j < LAL_NUM_IFO; j++)
    {
      thisCoincEntry = thisCoinc->snglInspiral[j];

      if ( thisCoincEntry )
      {
        /* allocate memory for a new sngl inspiral */
        if ( !snglHead )
        {
          thisSngl = snglHead = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
        }

        /* copy thisCoincEntry into our list */
        memcpy( thisSngl, thisCoincEntry, sizeof(SnglInspiralTable) );
        thisSngl->next = NULL;
        
        /* create an eventId and populate the id */
        eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
        eventId->id = LAL_INT8_C(1000000000) * (INT8) gpsStartTime->gpsSeconds
          + (INT8) eventNum;
        
        if ( slideNum < 0 )
        {
          eventId->id += LAL_INT8_C(100000)* (-1 *slideNum + 5000);
        }
        else 
        {
          eventId->id += LAL_INT8_C(100000) * slideNum;
        }
        thisSngl->event_id = eventId;
        eventId->snglInspiralTable = thisSngl;
      }
    }
  }

  *snglPtr = snglHead;
        

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALCoincCutSnglInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead
    )
/* </lalVerbatim> */
{
  SnglInspiralTable  *eventList = NULL;
  SnglInspiralTable  *prevEvent = NULL;
  SnglInspiralTable  *thisEvent = NULL;

  INITSTATUS( status, "LALCoincCutSnglInspiral", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that eventHead is non-null */
  ASSERT( eventHead, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
 

  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent = *eventHead;
  
  while ( thisEvent )
  {
    INT4           keepEvent = 0;
    EventIDColumn *eventID   = NULL;

    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->event_id )
    {
      for( eventID = tmpEvent->event_id; tmpEvent; tmpEvent = tmpEvent->next )
      {
        /* there is an eventID, check if there's a coinc inspiral */
        eventID = tmpEvent->event_id;
      
        if ( eventID->coincInspiralTable )
        {
          keepEvent = 1;
          break;
        }
      
        eventID = eventID->next;
      }
    }
     
    if ( keepEvent )
    {
      /* event is in a coinc so keep it */
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
      /* discard this trigger since it's not in a coinc */
      LALFreeSnglInspiral ( status->statusPtr, &tmpEvent );
    }
  }

  *eventHead = eventList; 

  DETATCHSTATUSPTR (status);
  RETURN (status);
}  

