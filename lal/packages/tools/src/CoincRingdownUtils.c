/*----------------------------------------------------------------------- 
 * 
 * File Name: CoincRingdownUtils.c
 *
 * Author: Goggin, L. M. based on CoincInspiralUtils.c by  Brady, P. R., Brown, D. A.,
 * and Fairhurst, S
 * 
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="CoincRingdownUtilsCV">
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
#include <lal/GenerateRing.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>

NRCSID( COINCRINGDOWNUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{CoincRingdownUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CoincRingdownUtilsCP}
\idx{LALCreateTwoIFOCoincList()}
\idx{LALCreateNIFORingdownCoincList()}
\idx{LALRemoveRepeatedRingdownCoincs()}
\idx{LALFreeCoincRingdown()}
\idx{XLALFreeCoincRingdown()}
\idx{LALAddSnglRingdownToCoinc()}
\idx{LALSnglInspiralCoincTest()}
\idx{LAExtractSnglInspiralFromCoinc()}
\idx{XLALRecreateCoincFromSngls()}
\idx{XLALGenerateCoherentBank()}
\idx{XLALInspiralDistanceCut()}
\idx{LALCoincCutSnglInspiral()}

\subsubsection*{Description}

\texttt{LALCreateTwoIFOCoincList()} takes in a linked list of single inspiral
tables and returns a list of two instrument coincidences.  The coincidence
requirements are given by the \texttt{accuracyParams}.  When single inspirals
from two different instruments are found to be coincident, the code creates a
new \texttt{coincInspiralTable} and uses \texttt{LALAddSnglInspiralToCoinc()}
to add the single inspirals to the coinc.  The function returns
\texttt{coincOutput} which is a pointer to the head of a linked list of
\texttt{CoincInspiralTable}s.

\texttt{LALCreateNIFORingdownCoincList()} takes linked list of
\texttt{CoincInspiralTable}s, assumed to contain (N-1) ifo coincidences and
creates all N ifo coincidences.  Both the input and output list of
\texttt{CoincInspiralTable}s are passed as \texttt{coincHead}.  

\texttt{LALRemoveRepeatedRingdownCoincs()} will remove any lower order coincidences 
if they are contained in a higher order coincidence.  For example, if an H1-L1
double coincident trigger is also part of an H1-H2-L1 triple coincident 
trigger, the double coincident trigger will be removed.  The head of the list 
of coincident triggers is passed and returned as \texttt{coincHead}.

\texttt{XLALFreeCoincRingdown()} \texttt{LALFreeCoincRingdown()} and  free the
memory associated to the \texttt{CoincInspiralTable} pointed to by
\texttt{coincPtr}.  This entails freeing the \texttt{CoincInspiralTable} as
well as any \texttt{eventId}s which point to the coinc.

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
stored in the \texttt{UINT8 id} field of the \texttt{eventIDColumn} of each
single inspiral which forms part of the coincidence.  The \texttt{id} is set
equal to $10^{9} \times$ \texttt{gpsStartTime} $+ 10^{5} \times$
\texttt{slideNum} $+$ event number. We do not assign multiple \texttt{id}
values to a given single inspiral table, but instead make multiple copies of
the table, each with a unique \texttt{id}.  

\texttt{XLALRecreateCoincFromSngls()} is used to recreate a list of coinc
inspirals from a list of \texttt{snglInspiralTable}s with populated 
\texttt{eventIDColumn}.  The code searches for entries in 
\texttt{snglInspiral} which have the same numerical value of the \texttt{id} 
field in the \texttt{eventIDColumn}.  

\texttt{XLALGenerateCoherentBank()} is used to generate a coherent bank from
a list of \texttt{coincInspiralTable}s.  The coherent bank has the same mass
parameters for each ifo.  These are currently chosen as the mass parameters
of the trigger in the coinc with the highest \texttt{snr}.  If the 
\texttt{ifos} field is not \texttt{NULL}, then a template is generated for
every ifo in \texttt{ifos}.  If it is \texttt{NULL} then templates are only
generated for those ifos which have triggers in the coinc.

\texttt{XLALInspiralDistanceCut()} is used to perform a distance cut between 
the triggers in a coincidence.  The distance cut uses the following algorithm:

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

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
LALCreateTwoIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincOutput,
    SnglRingdownTable          *snglInput,
    RingdownAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglRingdownTable            *currentTrigger[2];
  INT8                          currentTriggerNS[2];
  CoincRingdownTable           *coincHead = NULL;
  CoincRingdownTable           *thisCoinc = NULL;
  INT4                          numEvents = 0;
  INT4                          ifoNumber;
  INT8                          maxTimeDiff = 0;

  INITSTATUS( status, "LALCreateTwoIFOCoincList", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincOutput, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! *coincOutput, status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );

  memset( currentTriggerNS, 0, 2 * sizeof(INT8) );
  memset( currentTrigger, 0, 2 * sizeof(SnglRingdownTable *) );

  
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
        &(currentTrigger[0]->start_time) );

    /* set next trigger for comparison */
    currentTrigger[1] = currentTrigger[0]->next;
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
          &(currentTrigger[1]->start_time) );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {
      /* check that triggers pass coincidence test */
      LALCompareRingdowns( status->statusPtr, currentTrigger[0], 
          currentTrigger[1], accuracyParams );

      /* test whether we have coincidence */
      if ( accuracyParams->match )
      {
        LALInfo( status, "Found double coincident trigger,");
        /* create a 2 IFO coinc and store */
        if ( ! coincHead  )
        {
          coincHead = thisCoinc = (CoincRingdownTable *) 
            LALCalloc( 1, sizeof(CoincRingdownTable) );
        }
        else
        {
          thisCoinc = thisCoinc->next = (CoincRingdownTable *) 
            LALCalloc( 1, sizeof(CoincRingdownTable) );
        }

        /* Add the two triggers to the coinc */
        LALAddSnglRingdownToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[0] );
        LALAddSnglRingdownToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[1] );

        ++numEvents;

      }

      /* scroll on to the next sngl ringdown */
      
      if ( (currentTrigger[1] = currentTrigger[1]->next) )
      {
        LALGPStoINT8( status->statusPtr, &currentTriggerNS[1], 
            &(currentTrigger[1]->start_time) );
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

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
LALCreateNIFORingdownCoincList(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead,
    RingdownAccuracyList       *accuracyParams,
    INT4                        N
    )
/* </lalVerbatim> */
{
  INT4                          numEvents  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          ifoNum     = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincRingdownTable           *thisCoinc     = NULL;
  CoincRingdownTable           *lastCoinc     = NULL;
  CoincRingdownTable           *otherCoinc    = NULL;
  CoincRingdownTable           *nIfoCoincHead = NULL;
  CoincRingdownTable           *thisNIfoCoinc = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALCreateNIFORingdownCoincList", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  for( thisCoinc = *coincHead; thisCoinc; thisCoinc = thisCoinc->next)
  {
    lastCoinc = thisCoinc;

    /* check that this is an (N-1) coinc */
    if ( thisCoinc->numIfos == N - 1 )
    {
      /* look up the first single inspiral */
      for ( firstEntry = 0; firstEntry < LAL_NUM_IFO; firstEntry++)
      {
        if ( thisCoinc->snglRingdown[firstEntry] )
        {
          LALInfo( status, "Found the first entry in the coinc" );
          break;
        }
      }

      /* get the list of event IDs for this first entry */
      eventIDHead = thisCoinc->snglRingdown[firstEntry]->event_id;

      /* loop over the (N-1) ifo coincs that first entry is a member of 
       * and try to find an N ifo coinc */
      for( thisID = eventIDHead; thisID; thisID = thisID->next )
      {
        otherCoinc = thisID->coincRingdownTable;

        if( otherCoinc->numIfos == N - 1 )
        {
          /* loop over all singles which are alphabetically before the 
           * first one in thisCoinc */
          for( ifoNumber = 0; ifoNumber < firstEntry; ifoNumber++ )
          {
            /* test whether we have an N ifo coincidence */
            accuracyParams->match = 0;

            if ( otherCoinc->snglRingdown[ifoNumber] )
            {
              LALSnglRingdownCoincTest( status->statusPtr, thisCoinc, 
                  otherCoinc->snglRingdown[ifoNumber], accuracyParams );
            }

            if ( accuracyParams->match )
            {
              LALInfo( status, "We have found an N ifo coinc, storing");
              ++numEvents;
              
              /* create a N IFO coinc and store */
              if ( ! nIfoCoincHead  )
              {
                nIfoCoincHead = thisNIfoCoinc = (CoincRingdownTable *) 
                  LALCalloc( 1, sizeof(CoincRingdownTable) );
              }
              else
              {
                thisNIfoCoinc = thisNIfoCoinc->next = (CoincRingdownTable *) 
                  LALCalloc( 1, sizeof(CoincRingdownTable) );
              }

              /* add the single to the new N coinc */
              LALAddSnglRingdownToCoinc( status->statusPtr, &thisNIfoCoinc,
                  otherCoinc->snglRingdown[ifoNumber] );

              /* add the triggers from the (N-1) coinc to the new N coinc */
              for( ifoNum = 0; ifoNum < LAL_NUM_IFO; ifoNum++ )
              {
                if( thisCoinc->snglRingdown[ifoNum] )
                {
                  LALAddSnglRingdownToCoinc( status->statusPtr, &thisNIfoCoinc,
                      thisCoinc->snglRingdown[ifoNum] );
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
  if ( lastCoinc )
  {
    lastCoinc->next = nIfoCoincHead;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void
LALRemoveRepeatedRingdownCoincs(
    LALStatus                  *status,
    CoincRingdownTable        **coincHead
    )
{
  INT4                          removeThisCoinc  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincRingdownTable           *thisCoinc     = NULL;
  CoincRingdownTable           *prevCoinc     = NULL;
  CoincRingdownTable           *otherCoinc    = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALRemoveRepeatedRingdownCoincs", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  thisCoinc = *coincHead;

  while( thisCoinc )
  {
    /* look up the first single ringdown */
    for ( firstEntry = 0; firstEntry < LAL_NUM_IFO; firstEntry++)
    {
      if ( thisCoinc->snglRingdown[firstEntry] )
      {
        LALInfo( status, "Found the first entry in the coinc" );
        break;
      }
    }

    /* get the list of event IDs for this first entry */
    eventIDHead = thisCoinc->snglRingdown[firstEntry]->event_id;

    /* loop over the coincs that firstEntry is a member of and see if
     * thisCoinc is a subset of a higher order coinc */

    removeThisCoinc = 0;

    for( thisID = eventIDHead; thisID; thisID = thisID->next )
    {
      otherCoinc = thisID->coincRingdownTable;

      if( otherCoinc->numIfos >= thisCoinc->numIfos &&
          otherCoinc != thisCoinc )
      {
        /* we have a higher (or equal) coinc, thisCoinc could be a subset 
         * test whether all sngls in thisCoinc are also in otherCoinc */

        for( ifoNumber = firstEntry + 1; ifoNumber < LAL_NUM_IFO; 
            ifoNumber++ )
        {
          if ( thisCoinc->snglRingdown[ifoNumber] && 
              !(thisCoinc->snglRingdown[ifoNumber] == 
                otherCoinc->snglRingdown[ifoNumber]) )
          {
            LALInfo( status, "No Match");
            break;
          }
        }

        if ( ifoNumber == LAL_NUM_IFO )
        {
          LALInfo( status, "Removing lower order coinc");
          removeThisCoinc = 1;
          break;
        }
      }
    }

    if ( removeThisCoinc )
    {
      if ( !prevCoinc )
      {
        *coincHead = thisCoinc->next;
        LALFreeCoincRingdown( status->statusPtr, &thisCoinc );
        thisCoinc = *coincHead;
      }
      else 
      {
        prevCoinc->next = thisCoinc->next;
        LALFreeCoincRingdown( status->statusPtr, &thisCoinc );
        thisCoinc = prevCoinc->next;
      }
    }
    else 
    {
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }
  }

  
  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void
LALFreeCoincRingdown(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr
    )
{
  INITSTATUS( status, "LALFreeCoincRingdown", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  XLALFreeCoincRingdown( coincPtr );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}
    

void
XLALFreeCoincRingdown(
    CoincRingdownTable        **coincPtr
    )
{
  static const char *func = "FreeCoincRingdown";
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  EventIDColumn                *prevID     = NULL;
  EventIDColumn                *thisID     = NULL;
  SnglRingdownTable            *thisSngl   = NULL;

  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( (thisSngl = (*coincPtr)->snglRingdown[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that 
       * points to thisCoinc */
      thisID = thisSngl->event_id;
      prevID = NULL;

      while ( thisID )
      {
        /* test if thisID points to our coinc */ 
        if ( thisID->coincRingdownTable == *coincPtr )
        {
          if ( !prevID )
          {
            thisSngl->event_id = thisID->next;
          }
          else
          {
            prevID->next = thisID->next;
          }
          LALFree(thisID);
          break;
        }
        else
        {
          prevID = thisID;
          thisID = thisID->next;
        }
      }
    }
  }
  LALFree(*coincPtr);
}

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
LALAddSnglRingdownToCoinc(
    LALStatus                  *status,
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    )
/* </lalVerbatim> */
{
  CoincRingdownTable  *coincRingdown = NULL;
  EventIDColumn       *eventId = NULL;
  
  INITSTATUS( status, "LALAddSnglRingdownToCoinc", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglRingdown, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  coincRingdown = *coincPtr;
 
  /* allocate memory for new coinc if it doesn't exist */
  if (! coincRingdown )
  {
    coincRingdown = (CoincRingdownTable *) 
      LALCalloc( 1, sizeof(CoincRingdownTable) );
  }
  
  switch ( (snglRingdown->ifo)[0] ) 
  {
    case 'H':
      if ( !strcmp( snglRingdown->ifo, "H1" ) )
      {
        coincRingdown->snglRingdown[LAL_IFO_H1] = snglRingdown;
      }
      else if (!strcmp( snglRingdown->ifo, "H2" ) )
      {
        coincRingdown->snglRingdown[LAL_IFO_H2] = snglRingdown;
      }
      else
      {
        /* Invalid Hanford Detector */
        ABORT( status, LIGOMETADATAUTILSH_EDET, LIGOMETADATAUTILSH_MSGEDET );
      } 
      break;

    case 'L':
      coincRingdown->snglRingdown[LAL_IFO_L1] = snglRingdown;
      break;

    default:
      /* Invalid Detector Site */
      ABORT( status, LIGOMETADATAUTILSH_EDET, LIGOMETADATAUTILSH_MSGEDET );
  }

  ++(coincRingdown->numIfos);

  /* create an eventId for the single, populate it with the single and coinc */
  if ( ! snglRingdown->event_id )
  {
    eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    snglRingdown->event_id = eventId;
  }
  else
  {
     for( eventId = snglRingdown->event_id; eventId->next; 
         eventId = eventId->next);
     eventId = eventId->next = (EventIDColumn *) 
         LALCalloc( 1, sizeof(EventIDColumn) );
  }
  eventId->snglRingdownTable = snglRingdown;
  eventId->coincRingdownTable = coincRingdown;

  *coincPtr = coincRingdown;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
LALSnglRingdownCoincTest(
    LALStatus                  *status,
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown,
    RingdownAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;
  
  INITSTATUS( status, "LALSnglRingdownCoincTest", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );


  /* Loop over sngl_ringdowns contained in coinc_ringdown */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    thisCoincEntry = coincRingdown->snglRingdown[ifoNumber];

    if ( thisCoincEntry )
    {
      /* snglRingdown entry exists for this IFO, perform coincidence test */
      if ( ifoNumber == XLALIFONumber(snglRingdown->ifo) )
      {
        LALInfo( status, "We already have a coinc from this IFO" );
        accuracyParams->match = 0;
      }

      else
      {
        LALCompareRingdowns ( status->statusPtr, snglRingdown, 
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


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
LALExtractSnglRingdownFromCoinc(
    LALStatus                  *status,
    SnglRingdownTable         **snglPtr,
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )
/* </lalVerbatim> */
{
  SnglRingdownTable  *snglHead = NULL;
  SnglRingdownTable  *thisSngl = NULL;
  SnglRingdownTable  *thisCoincEntry = NULL;
  CoincRingdownTable *thisCoinc = NULL;
  EventIDColumn      *eventId = NULL;
  UINT4               eventNum = 1;
  INT4                j;
  
  INITSTATUS( status, "LALExtractCoincSngls", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( snglPtr, status, 
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( ! (*snglPtr), status, 
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  
  if ( !coincRingdown )
  {
    LALInfo( status, "Empty coincRingdown passed to LALExtractCoincSngls");
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }

  /* loop over the linked list of coinc ringdown */
  for( thisCoinc = coincRingdown; thisCoinc; thisCoinc = thisCoinc->next,
      ++eventNum)
  {
    /* loop over the interferometers */
    for ( j = 0; j < LAL_NUM_IFO; j++)
    {
      thisCoincEntry = thisCoinc->snglRingdown[j];

      if ( thisCoincEntry )
      {
        /* allocate memory for a new sngl ringdown */
        if ( !snglHead )
        {
          thisSngl = snglHead = (SnglRingdownTable *) 
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglRingdownTable *) 
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }

        /* copy thisCoincEntry into our list */
        memcpy( thisSngl, thisCoincEntry, sizeof(SnglRingdownTable) );
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
        eventId->snglRingdownTable = thisSngl;
      }
    }
  }

  *snglPtr = snglHead;
        

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
INT4 XLALCountCoincRingdown( CoincRingdownTable *head )
/* </lalVerbatim> */
{
  INT4 length;
  CoincRingdownTable *event;
  
  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}

#if 0
/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int 
XLALRecreateCoincFromSngls(
    CoincRingdownTable        **coincPtr,
    SnglRingdownTable          *snglRingdown
    )
/* </lalVerbatim> */
{
  static const char *func = "RecreateCoincFromSngls";
  SnglRingdownTable    *thisSngl  = NULL;
  CoincRingdownTable   *thisCoinc = NULL;
  CoincRingdownTable   *prevCoinc = NULL;
  CoincRingdownTable   *coincHead = NULL;
  UINT8                 eventId = 0;
  INT4                  numCoincs = 0;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;
  
    
  if ( !snglRingdown )
  {
    XLALPrintInfo( 
      "XLALRecreateCoincFromSngls: Empty snglRingdown passed as input" );
    return( 0 );
  }

  /* loop over the linked list of sngl ringdowns */
  for( thisSngl = snglringdown; thisSngl; thisSngl = thisSngl->next )
  {
    ifoNumber = XLALIFONumber( thisSngl->ifo );
    thisCoinc = coincHead;
    while ( thisCoinc )
    {
      /* loop over the interferometers to get the event_id*/
      for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
      {
        if ( thisCoinc->snglRingdown[ifoInCoinc] )
        {
          eventId = thisCoinc->snglRingdown[ifoInCoinc]->event_id->id;
          break;
        }
      }
      
      if ( thisSngl->event_id->id == eventId )
      {
        /* thisSngl is part of the coinc, so add it */
        if ( thisCoinc->snglRingdown[ifoNumber] )
        {
          /* already have an event for this ifo */
          XLALPrintError(
              "Already have a single from this ifo with event id %lld",
              eventId);
          /* free memory */
          while ( coincHead )
          {
            thisCoinc = coincHead;
            coincHead = coincHead->next;
            LALFree(thisCoinc);
          }
          XLAL_ERROR(func, XLAL_EDATA);
        }
        else
        {
          thisCoinc->snglRingdown[ifoNumber] = thisSngl;
          thisCoinc->numIfos += 1;
          thisSngl->event_id->coincRingdownTable = thisCoinc;
          break;
        }
      }

      /* proceed to the next coinc */
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }

    if ( thisSngl->event_id->id != eventId )
    {
      /* need to start a new coinc */
      if ( coincHead )
      {
        thisCoinc = prevCoinc->next = 
          LALCalloc( 1, sizeof(CoincRingdownTable) );
      }
      else
      {
        thisCoinc = coincHead = LALCalloc( 1, sizeof(CoincRingdownTable) );
      }
      if ( !thisCoinc )
      {
        /* out of memory: free memory + exit*/
        while ( coincHead )
        {
          thisCoinc = coincHead;
          coincHead = coincHead->next;
          LALFree( thisCoinc );
        }
        XLAL_ERROR(func,XLAL_ENOMEM);
      }

      thisCoinc->snglRingdown[ifoNumber] = thisSngl;
      thisCoinc->numIfos = 1;
      thisSngl->event_id->coincRingdownTable = thisCoinc;
      numCoincs +=1;
    }
  }

  *coincPtr = coincHead;
        
  return( numCoincs );
}

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int 
XLALGenerateCoherentBank(
    SnglRingdownTable         **coherentBank,
    CoincRingdownTable         *coincInput,
    CHAR                       *ifos
    )
/* </lalVerbatim> */
{
  static const char *func = "CreateCoherentBank";
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoMax  = LAL_UNKNOWN_IFO;
  SnglRingdownTable    *bankHead = NULL;
  SnglRingdownTable    *currentTrigger = NULL;
  CoincRingdownTable   *thisCoinc = NULL;
  INT4                  numTmplts = 0;
  
  if ( !coincInput )
  {
    XLALPrintInfo( 
      "XLALGenerateCoherentBank: Empty coincInput passed as input" );
    return( 0 );
  }
  
  for ( thisCoinc = coincInput; thisCoinc; thisCoinc = thisCoinc->next )
  {
    REAL4 max_snr = 0;

    /* loop over the interferometers to get the highest snr*/
    for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
    {
      if (( thisCoinc->snglRingdown[ifoInCoinc] ) && 
        (thisCoinc->snglRingdown[ifoInCoinc]->snr > max_snr) )
      {
        max_snr = thisCoinc->snglRingdown[ifoInCoinc]->snr;
        ifoMax = ifoInCoinc;
      }
    }
    
    for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
    {

      CHAR ifo[LIGOMETA_IFO_MAX];
       
      XLALReturnIFO( ifo, ifoNumber);

      /* decide whether we want a template for this ifo */
      if ( (thisCoinc->snglRingdown[ifoNumber] &&  !ifos) ||
           ( ifos && strstr(ifos,ifo)) )
      {
        numTmplts++;
        
        if( bankHead )
        {
          currentTrigger = currentTrigger->next = 
            LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else 
        {
          bankHead = currentTrigger = 
              LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        if ( !currentTrigger )
        {
          goto error;
        }
        /* copy the info from the loudest trigger */
        memcpy(currentTrigger, thisCoinc->snglRingdown[ifoMax], 
            sizeof(SnglRingdownTable));
        /* terminate the list */
        currentTrigger->next = NULL;
        currentTrigger->event_id = NULL;
        /* set the ifo */
        LALSnprintf( currentTrigger->ifo, LIGOMETA_IFO_MAX, ifo );
        /* set the event id */
        currentTrigger->event_id = LALCalloc( 1, sizeof(EventIDColumn) );
        if ( !(currentTrigger->event_id) )
        {
		      goto error;
        }
        currentTrigger->event_id->id = 
          thisCoinc->snglRingdown[ifoMax]->event_id->id;
        currentTrigger->event_id->snglRingdownTable = currentTrigger;
      }
    }
  }
  
  *coherentBank = bankHead;
  return( numTmplts );

  error:
  while ( bankHead )
  {
    currentTrigger = bankHead;
    bankHead = bankHead->next;
    XLALFreeSnglRingdown( &currentTrigger );
  }
  XLAL_ERROR(func,XLAL_ENOMEM);
  
}
#endif

#if 0  
/* <lalVerbatim file="CoincRingdownUtilsCP"> */
void
XLALRingdownDistanceCut(
    CoincRingdownTable        **coincRingdown,
    RingdownAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  static const char *func = "RingdownDistanceCut";
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincRingdownTable   *thisCoinc = NULL;
  CoincRingdownTable   *prevCoinc = NULL;
  CoincRingdownTable   *coincHead = NULL;

  thisCoinc = *coincRingdown;
  coincHead = NULL;
  
  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 kappaA = 0, kappaB = 0;
    REAL4 epsilonA = 0, epsilonB = 0;
    REAL4 snrA = 0, snrB = 0;
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;
    REAL4 distError = 0;

    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      kappaA = accuracyParams->ifoAccuracy[ifoA].kappa;
      epsilonA = accuracyParams->ifoAccuracy[ifoA].epsilon;
      
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        kappaB = accuracyParams->ifoAccuracy[ifoB].kappa;
        epsilonB = accuracyParams->ifoAccuracy[ifoB].epsilon;

        if( tmpCoinc->snglRingdown[ifoA] && ( kappaA || epsilonA ) 
            && tmpCoinc->snglRingdown[ifoB] && ( kappaB || epsilonB ) )
        {
          /* perform the distance consistency test */
          sigmasqA = tmpCoinc->snglRingdown[ifoA]->sigmasq;
          sigmasqB = tmpCoinc->snglRingdown[ifoB]->sigmasq;
          distA = tmpCoinc->snglRingdown[ifoA]->eff_distance;
          distB = tmpCoinc->snglRingdown[ifoB]->eff_distance;
          snrA = tmpCoinc->snglRingdown[ifoA]->snr;
          snrB = tmpCoinc->snglRingdown[ifoB]->snr;

          if( ( sigmasqA > sigmasqB && 
                fabs(distB - distA)/distA > epsilonB/snrB + kappaB ) ||
              ( sigmasqB > sigmasqA &&
                fabs(distA - distB)/distB > epsilonA/snrA + kappaA ) )
          {
            discardTrigger = 1;
            break;
          }
        }
      }
     
      if ( discardTrigger )
      {
        break;
      }
    }

    if( discardTrigger )
    {
      XLALFreeCoincRingdown( &tmpCoinc );
    }
    else
    {
      if ( ! coincHead )
      {
        coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
    }
  }
  *coincInspiral = coincHead;
}
#endif

#if 0
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
     pointer to the head of a linked list of tables which are in coincs */

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
#endif
