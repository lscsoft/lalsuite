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
\idx{LALSnglInspiralLookup()}
\idx{LALAddSnglInspiralToCoinc()}
\idx{LALCreateNewCoinc()}
\idx{LALSnglInspiralCoincTest()}

\subsubsection*{Description}

\noindent Blah.

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
    SnglInspiralAccuracy       *errorParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable            *currentTrigger[2];
  INT8                          currentTriggerNS[2];
  CoincInspiralTable           *coincHead = NULL;
  CoincInspiralTable           *thisCoinc = NULL;
  INT4                          maxTC = 0;
  INT4                          numEvents = 0;

  INITSTATUS( status, "LALCreateTwoIFOCoincList", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincOutput, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );
  memset( currentTrigger, 0, 2 * sizeof(SnglInspiralTable *) );

  maxTC = errorParams->dt + 44;
  
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

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTC )
    {
      /* check that triggers pass coincidence test */
      LALCompareSnglInspiral( status->statusPtr, currentTrigger[0], 
          currentTrigger[1], errorParams );

      /* test whether we have coincidence */
      if ( errorParams->match )
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
      currentTrigger[1] = currentTrigger[1]->next;
      LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
            &(currentTrigger[1]->end_time) );

    }
  }

  *coincOutput = coincHead;

  DETATCHSTATUSPTR (status);
  RETURN (status);
} 


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALSnglInspiralLookup(
    LALStatus          *status,
    SnglInspiralTable **snglInspiralPtr,
    CoincInspiralTable *coincInspiral,
    char               *ifo 
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *snglInspiralEntry;

  INITSTATUS( status, "LALSnglInspiralLookup", COINCINSPIRALUTILSC );

  switch ( ifo[0] ) 
  {
    case 'G':
      snglInspiralEntry = coincInspiral->G1Inspiral;
      break;

    case 'H':
      if ( !strcmp( ifo, "H1" ) )
      {
        snglInspiralEntry = coincInspiral->H1Inspiral;
      }
      else if (!strcmp( ifo, "H2" ) )
      {
        snglInspiralEntry = coincInspiral->H2Inspiral;
      }
      else
      {
        /* Invalid Hanford Detector */
        snglInspiralEntry = NULL;
      } 
      break;

    case 'L':
      snglInspiralEntry = coincInspiral->L1Inspiral;
      break;

    case 'T':
      snglInspiralEntry = coincInspiral->T1Inspiral;
      break;

    case 'V':
      snglInspiralEntry = coincInspiral->V1Inspiral;
      break;

    default:
      /* Invalid Detector Site */
      snglInspiralEntry = NULL;
  }
  *snglInspiralPtr = snglInspiralEntry;

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
      coincInspiral->G1Inspiral = snglInspiral;
      break;

    case 'H':
      if ( !strcmp( snglInspiral->ifo, "H1" ) )
      {
        coincInspiral->H1Inspiral = snglInspiral;
      }
      else if (!strcmp( snglInspiral->ifo, "H2" ) )
      {
        coincInspiral->H2Inspiral = snglInspiral;
      }
      else
      {
        /* Invalid Hanford Detector */
        ABORT( status, LIGOMETADATAUTILSH_EDET, LIGOMETADATAUTILSH_MSGEDET );
      } 
      break;

    case 'L':
      coincInspiral->L1Inspiral = snglInspiral;
      break;

    case 'T':
      coincInspiral->T1Inspiral = snglInspiral;
      break;

    case 'V':
      coincInspiral->V1Inspiral = snglInspiral;
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
LALCreateNewCoinc(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral
    )
/* </lalVerbatim> */
{
  CoincInspiralTable  *newCoincInspiral;
  INITSTATUS( status, "LALAddSnglInspiralToCoinc", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( *coincPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglInspiral, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  /* allocate memory for the new IFO coinc table */
  newCoincInspiral = (CoincInspiralTable *) 
    LALCalloc( 1, sizeof(CoincInspiralTable) );


  /* copy over the single IFO event */
  memcpy( newCoincInspiral, coincInspiral, sizeof(CoincInspiralTable) );


  /* add the additional sngl_inspiral to the new coinc */
  LALAddSnglInspiralToCoinc( status->statusPtr, &newCoincInspiral, 
      snglInspiral );

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALSnglInspiralCoincTest(
    LALStatus                  *status,
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    SnglInspiralAccuracy       *errorParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  j = 0;
  static char         *ifoList[] = {"Unknown IFO", "G1", "H1", "H2", 
    "L1", "T1", "V1"};


  INITSTATUS( status, "LALSnglInspiralCoincTest", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  /* Loop over sngl_inspirals contained in coinc_inspiral */
  for ( j = 1; j < 7; j++)
  {
    LALSnglInspiralLookup( status->statusPtr, &thisCoincEntry, coincInspiral, 
        ifoList[j] );

    if ( thisCoincEntry )
    {
      /* snglInspiral entry exists for this IFO, perform coincidence test */
      if ( !strcmp( snglInspiral->ifo, thisCoincEntry->ifo ) )
      {
        LALInfo( status, "We already have a coinc from this IFO" );
        errorParams->match = 0;
      }

      else
      {
        LALCompareSnglInspiral ( status->statusPtr, snglInspiral, 
            thisCoincEntry, errorParams );
      }
      /* set match to zero if no match.  Keep same if match */
      match *= errorParams->match;
    }
  }
  /* returm errorParams->match to be 1 if we match, zero otherwise */
  errorParams->match = match;
  if ( errorParams->match == 0 ) LALInfo( status, "Coincidence test failed" );
  if ( errorParams->match == 1 ) LALInfo( status, "Coincidence test passed" );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}



