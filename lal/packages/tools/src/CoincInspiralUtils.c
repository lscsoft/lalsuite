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
\idx{LALSnglInspiralLookup()}


\subsubsection*{Description}

\texttt{LALCreateTwoIFOCoincList()} takes in a linked list of single inspiral
tables and returns a list of two instrument coincidences.  The coincidence
requirements are given by the \texttt{errorParams}.  When single inspirals from
two different instruments are found to be coincident, the code creates a new
\texttt{coincInspiralTable} and uses \texttt{LALAddSnglInspiralToCoinc()} to
add the single inspirals to the coinc.  The function returns
\texttt{coincOutput} which is a pointer to the head of a linked list of
\texttt{CoincInspiralTable}s.

\texttt{LALAddSnglInspiralToCoinc()} adds a pointer to a single inspiral table
to a coinc inspiral table.  Upon entry, if \texttt{coincPtr} points to a
\texttt{NULL} coinc inspiral table, the table is created before a pointer to
the single inspiral table is added.  Additionally, an \texttt{eventId} table is
created for the single inspiral table.  This points to both the single and
coinc inspirals.  If an \texttt{eventId} already exists for the single
inspiral, another eventId table is added to the linked list.  The linked list
of \texttt{eventId}s associated to a single inspiral table allow us to easily
determine which coincident events each single is a part of.

The function \texttt{LALSnglInspiralLookup()} can be used to retrieve single
inspiral entries from a coinc inspiral table.  Given an \texttt{ifo} and a
\texttt{coincInspiral}, it returns a pointer to the relevant single inspiral
table.

\texttt{LALSnglInspiralCoincTest()} tests for coincidence between a single
inspiral and a coinc inspiral.  It works by testing for coincidence between
each non-null entry in the coinc inspiral and the single.  This is done using
\texttt{LALCompareSnglInspiral()}.  If all members of the coinc are found to be
coincident with the single, the \texttt{errorParams.match} is set to 1,
otherwise to 0.

\texttt{LALExtractCoincSngls()} extracts the information from a linked list of
\texttt{coincInspiralTable}s and returns it as a linked list of
\texttt{snglInspiralTable}s.  Thus, the output \texttt{snglPtr} is a pointer to
a linked list of single inspiral tables.  That list contains only single
inspirals which are found in coincidence.  In order to preserve the coincidence
information, we assign to each coincident event an integer value.  This is
stored in the \texttt{UINT4 id} field of the \texttt{eventIDColumn} of each
single inspiral which forms part of the coincidence.  We do not assign multiple
\texttt{id} values to a given single inspiral table, but instead make multiple
copies of the table, each with a unique \texttt{id}.  


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
LALSnglInspiralLookup(
    LALStatus            *status,
    SnglInspiralTable   **snglInspiralPtr,
    CoincInspiralTable   *coincInspiral,
    InterferometerLabel   ifo 
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *snglInspiralEntry;

  INITSTATUS( status, "LALSnglInspiralLookup", COINCINSPIRALUTILSC );

  switch ( ifo ) 
  {
    case g1:
      snglInspiralEntry = coincInspiral->G1Inspiral;
      break;

    case h1:
      snglInspiralEntry = coincInspiral->H1Inspiral;
      break;

    case h2:
      snglInspiralEntry = coincInspiral->H2Inspiral;
      break;

    case l1:
      snglInspiralEntry = coincInspiral->L1Inspiral;
      break;

    case t1:
      snglInspiralEntry = coincInspiral->T1Inspiral;
      break;

    case v1:
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
  
  INITSTATUS( status, "LALSnglInspiralCoincTest", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  /* Loop over sngl_inspirals contained in coinc_inspiral */
  for ( j = 1; j < 7; j++)
  {
    LALSnglInspiralLookup( status->statusPtr, &thisCoincEntry, coincInspiral, 
        j );

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



/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALExtractCoincSngls(
    LALStatus                  *status,
    SnglInspiralTable         **snglPtr,
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime
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
  ASSERT( coincInspiral, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  /* loop over the linked list of coinc inspirals */
  for( thisCoinc = coincInspiral; thisCoinc; thisCoinc = thisCoinc->next,
      ++eventNum)
  {
    /* loop over the interferometers */
    for ( j = 1; j < 7; j++)
    {
      LALSnglInspiralLookup( status->statusPtr, &thisCoincEntry, thisCoinc, 
        j );

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
        /* XXX Commented until id is changed from UINT4 to UINT8 
        eventId->id = LAL_INT8_C(1000000000) * (INT8) gpsStartTime->gpsSeconds 
          + (INT8) eventNum; XXX */
       
        thisSngl->event_id = eventId;
        eventId->snglInspiralTable = thisSngl;
      }
    }
  }

  *snglPtr = snglHead;
        

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
