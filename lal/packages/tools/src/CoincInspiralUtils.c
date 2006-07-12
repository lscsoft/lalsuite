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
#include <lal/XLALError.h>

NRCSID( COINCINSPIRALUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{CoincInspiralUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CoincInspiralUtilsCP}
\idx{LALCreateTwoIFOCoincList()}
\idx{LALCreateNIFOCoincList()}
\idx{LALRemoveRepeatedCoincs()}
\idx{LALFreeCoincInspiral()}
\idx{XLALFreeCoincInspiral()}
\idx{LALAddSnglInspiralToCoinc()}
\idx{LALSnglInspiralCoincTest()}
\idx{LALExtractSnglInspiralFromCoinc()}
\idx{XLALRecreateCoincFromSngls()}
\idx{XLALGenerateCoherentBank()}
\idx{XLALInspiralPsi0Psi3CutBCVC()}
\idx{XLALInspiralIotaCutBCVC()}
\idx{LALInspiralDistanceCutCleaning()}
\idx{XLALInspiralDistanceCutBCVC()}
\idx{XLALInspiralDistanceCut()}
\idx{LALCoincCutSnglInspiral()}
\idx{XLALCoincInspiralTimeNS()}
\idx{XLALCoincInspiralStat()}
\idx{XLALClusterCoincInspiralTable()}
\idx{XLALCoincInspiralIfos()}
\idx{XLALCoincInspiralIfosCut()}
\idx{XLALCoincInspiralIdNumber()}
\idx{XLALCoincInspiralSlideCut()}
\idx{XLALInspiralSNRCutBCV2()}
\idx{XLALSNRCutCoincInspiral()}
\idx{XLALSortCoincInspiral()}
\idx{XLALCompareCoincInspiralByTime()}
\idx{XLALCountInspiralTable()}

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

\texttt{LALRemoveRepeatedCoincs()} will remove any lower order coincidences 
if they are contained in a higher order coincidence.  For example, if an H1-L1
double coincident trigger is also part of an H1-H2-L1 triple coincident 
trigger, the double coincident trigger will be removed.  The head of the list 
of coincident triggers is passed and returned as \texttt{coincHead}.

\texttt{XLALFreeCoincInspiral()} \texttt{LALFreeCoincInspiral()} and  free the
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
the triggers in a coincidence.  The distance cut analyzes triggers from two
different instruments.  It determines which instrument was the most sensitive
by comparing the \texttt{sigmasq} values of the two triggers, the instrument 
with the greatest range is designated ifo A, the other ifo B.  It then discards
and triggers for which
%
\begin{equation}
\frac{|distB - distA|}{distA} > \frac{epsilonB}{snrB} + kappaB   
\end{equation}
%

\texttt{LALCoincCutSnglInspiral()} extracts all single inspirals from a
specific ifo which are in coinc inspirals.  The output \texttt{snglPtr} is a
pointer to a linked list of single inspiral tables.  That list contains only
single inspirals from the specified \texttt{ifo} which are found in
coincidence. 

\texttt{XLALCountCoincInspiral()} scans through a linked list of coincidence
inspiral table and counts the number of events. This count is returned 
as \texttt{numTrigs}.



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
  CoincInspiralTable           *lastCoinc     = NULL;
  CoincInspiralTable           *otherCoinc    = NULL;
  CoincInspiralTable           *nIfoCoincHead = NULL;
  CoincInspiralTable           *thisNIfoCoinc = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALCreateNIFOCoincList", COINCINSPIRALUTILSC );
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
        if ( thisCoinc->snglInspiral[firstEntry] )
        {
          LALInfo( status, "Found the first entry in the coinc" );
          break;
        }
      }

      /* get the list of event IDs for this first entry */
      eventIDHead = thisCoinc->snglInspiral[firstEntry]->event_id;

      /* loop over the (N-1) ifo coincs that first entry is a member of 
       * and try to find an N ifo coinc */
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
  if ( lastCoinc )
  {
    lastCoinc->next = nIfoCoincHead;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

void
LALRemoveRepeatedCoincs(
    LALStatus                  *status,
    CoincInspiralTable        **coincHead
    )
{
  INT4                          removeThisCoinc  = 0;
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber          firstEntry = LAL_UNKNOWN_IFO;
  CoincInspiralTable           *thisCoinc     = NULL;
  CoincInspiralTable           *prevCoinc     = NULL;
  CoincInspiralTable           *otherCoinc    = NULL;
  EventIDColumn                *eventIDHead   = NULL;
  EventIDColumn                *thisID        = NULL;


  INITSTATUS( status, "LALRemoveRepeatedCoincs", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* loop over all the coincidences in the list */
  thisCoinc = *coincHead;

  while( thisCoinc )
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

    /* loop over the coincs that firstEntry is a member of and see if
     * thisCoinc is a subset of a higher order coinc */

    removeThisCoinc = 0;

    for( thisID = eventIDHead; thisID; thisID = thisID->next )
    {
      otherCoinc = thisID->coincInspiralTable;

      if( otherCoinc->numIfos >= thisCoinc->numIfos &&
          otherCoinc != thisCoinc )
      {
        /* we have a higher (or equal) coinc, thisCoinc could be a subset 
         * test whether all sngls in thisCoinc are also in otherCoinc */

        for( ifoNumber = firstEntry + 1; ifoNumber < LAL_NUM_IFO; 
            ifoNumber++ )
        {
          if ( thisCoinc->snglInspiral[ifoNumber] && 
              !(thisCoinc->snglInspiral[ifoNumber] == 
                otherCoinc->snglInspiral[ifoNumber]) )
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
        LALFreeCoincInspiral( status->statusPtr, &thisCoinc );
        thisCoinc = *coincHead;
      }
      else 
      {
        prevCoinc->next = thisCoinc->next;
        LALFreeCoincInspiral( status->statusPtr, &thisCoinc );
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
LALFreeCoincInspiral(
    LALStatus                  *status,
    CoincInspiralTable        **coincPtr
    )
{
  INITSTATUS( status, "LALFreeCoincInspiral", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  XLALFreeCoincInspiral( coincPtr );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}
    

int
XLALFreeCoincInspiral(
    CoincInspiralTable        **coincPtr
    )
{
  InterferometerNumber          ifoNumber  = LAL_UNKNOWN_IFO;
  EventIDColumn                *prevID     = NULL;
  EventIDColumn                *thisID     = NULL;
  SnglInspiralTable            *thisSngl   = NULL;
  SimInspiralTable             *thisSim    = NULL;

  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( (thisSngl = (*coincPtr)->snglInspiral[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that 
       * points to thisCoinc */
      thisID = thisSngl->event_id;
      prevID = NULL;

      while ( thisID )
      {
        /* test if thisID points to our coinc */ 
        if ( thisID->coincInspiralTable == *coincPtr )
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
  if ( (thisSim = (*coincPtr)->simInspiral) )
  {
    /* loop over the list of eventID's until we get to the one that 
     * points to thisCoinc */
    thisID = thisSim->event_id;
    prevID = NULL;

    while ( thisID )
    {
      /* test if thisID points to our coinc */ 
      if ( thisID->coincInspiralTable == *coincPtr )
      {
        if ( !prevID )
        {
          thisSim->event_id = thisID->next;
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
  LALFree(*coincPtr);
  return( 0 );
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
  INITSTATUS( status, "LALExtractCoincSngls", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  *snglPtr = XLALExtractSnglInspiralFromCoinc( coincInspiral, gpsStartTime, 
      slideNum );
        

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
SnglInspiralTable *
XLALExtractSnglInspiralFromCoinc(
    CoincInspiralTable         *coincInspiral,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )
/* </lalVerbatim> */
{
  static const char *func = "ExtractSnglInspiralFromCoinc";
  SnglInspiralTable  *snglHead = NULL;
  SnglInspiralTable  *thisSngl = NULL;
  SnglInspiralTable  *thisCoincEntry = NULL;
  CoincInspiralTable *thisCoinc = NULL;
  EventIDColumn      *eventId = NULL;
  UINT4               eventNum = 1;
  INT4                j;
  
  if ( !coincInspiral )
  {
    XLALPrintInfo( 
      "XLALExtractSnglInspiralFromCoinc: Empty coincInspiral passed as input" 
      );
    return( NULL );
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
        if ( thisCoincEntry->event_id->id )
        {
          /* event id number exists, use it */
          eventId->id = thisCoincEntry->event_id->id;
        }
        else if ( gpsStartTime )
        {
          eventId->id = LAL_INT8_C(1000000000) * 
            (INT8) gpsStartTime->gpsSeconds + (INT8) eventNum;
        }
        else 
        {
          XLALPrintError(
              "Event does not have id and no GPS start time given" );
          while ( snglHead )
          {
            thisSngl = snglHead;
            snglHead = snglHead->next;
            XLALFreeSnglInspiral( &thisSngl );
          }
          XLAL_ERROR_NULL(func,XLAL_EIO);
        }

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

  return( snglHead );
        
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int 
XLALRecreateCoincFromSngls(
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable          *snglInspiral
    )
/* </lalVerbatim> */
{
  static const char *func = "RecreateCoincFromSngls";
  SnglInspiralTable    *thisSngl  = NULL;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;
  UINT8                 eventId = 0;
  INT4                  numCoincs = 0;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;
  
    
  if ( !snglInspiral )
  {
    XLALPrintInfo( 
      "XLALRecreateCoincFromSngls: Empty snglInspiral passed as input" );
    return( 0 );
  }

  /* loop over the linked list of sngl inspirals */
  for( thisSngl = snglInspiral; thisSngl; thisSngl = thisSngl->next )
  {
    ifoNumber = XLALIFONumber( thisSngl->ifo );
    thisCoinc = coincHead;
    while ( thisCoinc )
    {
      /* loop over the interferometers to get the event_id*/
      for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
      {
        if ( thisCoinc->snglInspiral[ifoInCoinc] )
        {
          eventId = thisCoinc->snglInspiral[ifoInCoinc]->event_id->id;
          break;
        }
      }
      
      if ( thisSngl->event_id->id == eventId )
      {
        /* thisSngl is part of the coinc, so add it */
        if ( thisCoinc->snglInspiral[ifoNumber] )
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
          thisCoinc->snglInspiral[ifoNumber] = thisSngl;
          thisCoinc->numIfos += 1;
          thisSngl->event_id->coincInspiralTable = thisCoinc;
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
          LALCalloc( 1, sizeof(CoincInspiralTable) );
      }
      else
      {
        thisCoinc = coincHead = LALCalloc( 1, sizeof(CoincInspiralTable) );
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

      thisCoinc->snglInspiral[ifoNumber] = thisSngl;
      thisCoinc->numIfos = 1;
      thisSngl->event_id->coincInspiralTable = thisCoinc;
      numCoincs +=1;
    }
  }
  /* discard any coincs which are from a single ifo */
  thisCoinc = coincHead; 
  coincHead = NULL;
  numCoincs = 0;

  while ( thisCoinc )
  {
    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( tmpCoinc->numIfos > 1 )
    {
      /* ifos match so keep tmpCoinc */
      if ( ! coincHead  )
      {
        coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
      ++numCoincs;
    }
    else
    {
      /* discard tmpCoinc */
      XLALFreeCoincInspiral( &tmpCoinc );
    }
  }

  *coincPtr = coincHead;
        
  return( numCoincs );
}

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int 
XLALGenerateCoherentBank(
    SnglInspiralTable         **coherentBank,
    CoincInspiralTable         *coincInput,
    CHAR                       *ifos
    )
/* </lalVerbatim> */
{
  static const char *func = "CreateCoherentBank";
  InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoMax  = LAL_UNKNOWN_IFO;
  SnglInspiralTable    *bankHead = NULL;
  SnglInspiralTable    *currentTrigger = NULL;
  CoincInspiralTable   *thisCoinc = NULL;
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
      if (( thisCoinc->snglInspiral[ifoInCoinc] ) && 
        (thisCoinc->snglInspiral[ifoInCoinc]->snr > max_snr) )
      {
        max_snr = thisCoinc->snglInspiral[ifoInCoinc]->snr;
        ifoMax = ifoInCoinc;
      }
    }
    
    for (ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
    {

      CHAR ifo[LIGOMETA_IFO_MAX];
       
      XLALReturnIFO( ifo, ifoNumber);

      /* decide whether we want a template for this ifo */
      if ( (thisCoinc->snglInspiral[ifoNumber] &&  !ifos) ||
           ( ifos && strstr(ifos,ifo)) )
      {
        numTmplts++;
        
        if( bankHead )
        {
          currentTrigger = currentTrigger->next = 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        else 
        {
          bankHead = currentTrigger = 
              LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        if ( !currentTrigger )
        {
          goto error;
        }
        /* copy the info from the loudest trigger */
        memcpy(currentTrigger, thisCoinc->snglInspiral[ifoMax], 
            sizeof(SnglInspiralTable));
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
          thisCoinc->snglInspiral[ifoMax]->event_id->id;
        currentTrigger->event_id->snglInspiralTable = currentTrigger;
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
    XLALFreeSnglInspiral( &currentTrigger );
  }
  XLAL_ERROR(func,XLAL_ENOMEM);
  
}
  
/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALInspiralPsi0Psi3CutBCVC(
    CoincInspiralTable        **coincInspiral
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;

  INT4  discardTrigger = 0;
  REAL4 distA = 0, distB = 0;
  REAL4 sigA, sigB;
  REAL4 psi0A, psi0B, psi3A, psi3B, snr, snrA, snrB, x, y, X, Y, theta=0.040;

  thisCoinc = *coincInspiral;
  coincHead = NULL;
   
  /* loop over the coincindent triggers */
  while( thisCoinc )
  {
    discardTrigger=0;

    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;
      

    /* loop over all IFO combinations */
    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        if( tmpCoinc->snglInspiral[ifoA] 
            && tmpCoinc->snglInspiral[ifoB]  )
        {
          /* perform the distance consistency test */
          psi0B = tmpCoinc->snglInspiral[ifoB]->psi0;
          psi0A = tmpCoinc->snglInspiral[ifoA]->psi0;
          psi3B = tmpCoinc->snglInspiral[ifoB]->psi3;
          psi3A = tmpCoinc->snglInspiral[ifoA]->psi3;
          snrA = tmpCoinc->snglInspiral[ifoA]->snr;
          snrB = tmpCoinc->snglInspiral[ifoB]->snr;
          if (snrA<snrB )
          {
            snr = snrB;
          }
          else
          {
            snr = snrA;
          }
        
          x = (psi0A - psi0B) / snr ;
          y = (psi3A - psi3B) / snr ;
          X = x * cos(theta) - y * sin(theta);
          Y = x * sin(theta) + y * cos(theta);

          if ( ((X*X/5000/5000) + (Y*Y/50/50)) >  1 )     
          {
            discardTrigger = 1;
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
        XLALFreeCoincInspiral( &tmpCoinc );
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
    
/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALInspiralIotaCutBCVC(
CoincInspiralTable        **coincInspiral,
InspiralAccuracyList       *accuracyParams
    
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;

  INT4  discardTrigger = 0;
  REAL4 snrA, snrB, sigA, sigB, distA, distB;
  REAL4 iota, iotaCutH1H2, iotaCutH1L1; 

  thisCoinc = *coincInspiral;
  coincHead = NULL;
   

  /* loop over the coincindent triggers */
  while( thisCoinc )
  {
    CoincInspiralTable *tmpCoinc = thisCoinc;

    discardTrigger=0;
    thisCoinc = thisCoinc->next;

    iotaCutH1H2=accuracyParams->iotaCutH1H2;
    iotaCutH1L1=accuracyParams->iotaCutH1L1;      

    /* loop over all IFO combinations */
    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        /*epsilonB = accuracyParams->ifoAccuracy[ifoB].epsilon;*/

        if( tmpCoinc->snglInspiral[ifoA] 
            && tmpCoinc->snglInspiral[ifoB]  )
        {
          /* perform the distance consistency test */
          sigA = tmpCoinc->snglInspiral[ifoA]->sigmasq;
          sigB = tmpCoinc->snglInspiral[ifoB]->sigmasq;
          snrA = tmpCoinc->snglInspiral[ifoA]->snr;
          snrB = tmpCoinc->snglInspiral[ifoB]->snr;
          distA = sigA*sigA/snrA;
          distB = sigB*sigB/snrB;
          iota=2* fabs(distA-distB)/(distA+distB);

	  /* check the iota value */
          if( ( ((ifoA == LAL_IFO_H1)  && (ifoB == LAL_IFO_H2) )
		|| (  (ifoA == LAL_IFO_H2)  && (ifoB == LAL_IFO_H1))  )
		 && iota>iotaCutH1H2 )
	    {
              discardTrigger = 1; 
            }
          if( ( ((ifoA == LAL_IFO_H1)  && (ifoB == LAL_IFO_L1) )
		|| (  (ifoA == LAL_IFO_L1)  && (ifoB == LAL_IFO_H1))  )
		 && iota>iotaCutH1L1 )
	    {
              discardTrigger = 1; 
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
        XLALFreeCoincInspiral( &tmpCoinc );
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
LALInspiralDistanceCutCleaning(
    LALStatus                  *status,
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams,
    REAL4			snrThreshold,
    SummValueTable            **summValueList,
    LALSegList                 *vetoSegsH1, 
    LALSegList                 *vetoSegsH2
    )
/* </lalVerbatim> */
{
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;
  REAL4 dH1, dH2, snrH1, snrH2; 
  REAL4 iotaCut;
  INITSTATUS( status, "LALInspiralDistanceCutCleaning", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );
  
  thisCoinc = *coincInspiral;
  coincHead = NULL;
 
  if (!vetoSegsH1 || !vetoSegsH2 )
  XLALPrintWarning("LALInspiralDistanceCutCleaning: Warning, no veto list provided. you should provide h1 and h2 ones. ");



  /*  compute the iota accuracy  */
  iotaCut  =  1/((2-accuracyParams->ifoAccuracy[LAL_IFO_H1].kappa)
	/(2+accuracyParams->ifoAccuracy[LAL_IFO_H1].kappa));

  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 snrA = 0, snrB = 0;
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;

    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;
    
    if( tmpCoinc->snglInspiral[LAL_IFO_H1] && !tmpCoinc->snglInspiral[LAL_IFO_H2])
    {
      snrH1 = tmpCoinc->snglInspiral[LAL_IFO_H1]->snr;
      LALDistanceScanSummValueTable(status->statusPtr, summValueList,
   tmpCoinc->snglInspiral[LAL_IFO_H1]->end_time, "H1",  &dH1);
      
  LALDistanceScanSummValueTable(status->statusPtr, summValueList,
   tmpCoinc->snglInspiral[LAL_IFO_H1]->end_time, "H2",  &dH2);

      
      /* iota =1 */
      if (dH2/dH1*snrH1 > snrThreshold * iotaCut)
      {
        if (vetoSegsH2)
  {
    if (!XLALSegListSearch( vetoSegsH2, 
    &(tmpCoinc->snglInspiral[LAL_IFO_H1]->end_time)))
        { 
      discardTrigger =1;
    }
        }
  else
  {
  discardTrigger = 1;
  }
      }
    }

    if( tmpCoinc->snglInspiral[LAL_IFO_H2] && !tmpCoinc->snglInspiral[LAL_IFO_H1])
    {
      snrH2 = tmpCoinc->snglInspiral[LAL_IFO_H2]->snr;
      LALDistanceScanSummValueTable(status->statusPtr, summValueList,
  tmpCoinc->snglInspiral[LAL_IFO_H2]->end_time, "H1",  &dH1);
      LALDistanceScanSummValueTable(status->statusPtr, summValueList, 
  tmpCoinc->snglInspiral[LAL_IFO_H2]->end_time, "H2",  &dH2);
      /* iota = 1 */
      if (dH1/dH2*snrH2 > snrThreshold *iotaCut)
      {
        if (vetoSegsH1)
  {
    if (!XLALSegListSearch( vetoSegsH1, 
    &(tmpCoinc->snglInspiral[LAL_IFO_H2]->end_time)))
        { 
      discardTrigger =1;
    }
  } 
  else
  {
    discardTrigger = 1; 
        }
      }
    }

/*    discardTrigger=0;*/
    if( discardTrigger )
    {
      XLALFreeCoincInspiral( &tmpCoinc );
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
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALInspiralDistanceCutBCVC(
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;

  REAL4 kappaA=0, kappaB=0, epsilonA=0, epsilonB=0;
  thisCoinc = *coincInspiral;
  coincHead = NULL;
  
  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 snrA = 0, snrB = 0;
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;

    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      kappaA = accuracyParams->ifoAccuracy[ifoA].kappa;
      epsilonA = accuracyParams->ifoAccuracy[ifoA].epsilon;
      
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        kappaB = accuracyParams->ifoAccuracy[ifoB].kappa;
        epsilonB = accuracyParams->ifoAccuracy[ifoB].epsilon;

        if( tmpCoinc->snglInspiral[ifoA] && ( kappaA || epsilonA ) 
            && tmpCoinc->snglInspiral[ifoB] && ( kappaB || epsilonB ) )
        {
          /* perform the distance consistency test */
          sigmasqA = tmpCoinc->snglInspiral[ifoA]->sigmasq;
          sigmasqB = tmpCoinc->snglInspiral[ifoB]->sigmasq;
          snrA = tmpCoinc->snglInspiral[ifoA]->snr;
          snrB = tmpCoinc->snglInspiral[ifoB]->snr;
          distA = sigmasqA*sigmasqA/snrA;
          distB = sigmasqB*sigmasqB/snrB;
          if( ( fabs(distB - distA)/(distA + distB) > epsilonB/snrB + kappaB ) ||
              ( fabs(distA - distB)/(distB + distA) > epsilonA/snrA + kappaA ) )
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
      XLALFreeCoincInspiral( &tmpCoinc );
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALInspiralDistanceCut(
    CoincInspiralTable        **coincInspiral,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoA = LAL_UNKNOWN_IFO;  
  InterferometerNumber  ifoB = LAL_UNKNOWN_IFO;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;

  thisCoinc = *coincInspiral;
  coincHead = NULL;
  
  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 kappaA = 0, kappaB = 0;
    REAL4 epsilonA = 0, epsilonB = 0;
    REAL4 snrA = 0, snrB = 0;
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;

    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      kappaA = accuracyParams->ifoAccuracy[ifoA].kappa;
      epsilonA = accuracyParams->ifoAccuracy[ifoA].epsilon;
      
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        kappaB = accuracyParams->ifoAccuracy[ifoB].kappa;
        epsilonB = accuracyParams->ifoAccuracy[ifoB].epsilon;

        if( tmpCoinc->snglInspiral[ifoA] && ( kappaA || epsilonA ) 
            && tmpCoinc->snglInspiral[ifoB] && ( kappaB || epsilonB ) )
        {
          /* perform the distance consistency test */
          sigmasqA = tmpCoinc->snglInspiral[ifoA]->sigmasq;
          sigmasqB = tmpCoinc->snglInspiral[ifoB]->sigmasq;
          distA = tmpCoinc->snglInspiral[ifoA]->eff_distance;
          distB = tmpCoinc->snglInspiral[ifoB]->eff_distance;
          snrA = tmpCoinc->snglInspiral[ifoA]->snr;
          snrB = tmpCoinc->snglInspiral[ifoB]->snr;

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
      XLALFreeCoincInspiral( &tmpCoinc );
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


INT8 
XLALCoincInspiralTimeNS (
    CoincInspiralTable         *coincInspiral
    )
{
  static const char *func = "XLALCoincInspiralTimeNS";
  InterferometerNumber  ifoNumber;
  INT8 endTime = 0;
  
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( coincInspiral->snglInspiral[ifoNumber] )
    {
      endTime = XLALGPStoINT8( 
          &(coincInspiral->snglInspiral[ifoNumber]->end_time) );
      return(endTime);
    }
  }
  XLAL_ERROR(func,XLAL_EIO);
}

REAL4 
XLALCoincInspiralStat(
    CoincInspiralTable         *coincInspiral,
    CoincInspiralStatistic      coincStat,
    CoincInspiralBittenLParams *bittenLParams
    )
{
  InterferometerNumber  ifoNumber;
  SnglInspiralTable    *snglInspiral;
  REAL4                 statValues[LAL_NUM_IFO];
  REAL4 statValue = 0;
  INT4  i;
  
  if( coincStat == no_stat )
  {
    return(0);
  }

  /* for bittenL only*/
  if( coincStat == bitten_l)
  {
    for ( i = 0; i < LAL_NUM_IFO ; i++)
    {
      statValues[i] = 1e9; /* sufficiently high values */
    }
  }
 

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( (snglInspiral = coincInspiral->snglInspiral[ifoNumber]) )
    {
      if ( coincStat == snrsq )
      {        
        statValue += snglInspiral->snr * snglInspiral->snr;
      }
      else if ( coincStat == effective_snrsq )
      {
        REAL4 tmp_snr = snglInspiral->snr;
        REAL4 tmp_chisq = snglInspiral->chisq;
        /* XXX Assuming that chisq_dof contains the number of bins, not dof */
        REAL4 tmp_bins = snglInspiral->chisq_dof;
        
        statValue += tmp_snr * tmp_snr / 
          sqrt ( tmp_chisq/(2*tmp_bins-2) * (1+tmp_snr*tmp_snr/250) ) ;
      }

      else if ( coincStat == bitten_l )
      {
        statValues[ifoNumber] = bittenLParams->param_a[ifoNumber] 
                * snglInspiral->snr 
                - bittenLParams->param_b[ifoNumber];
        statValue += snglInspiral->snr * snglInspiral->snr ;          
      }
      else if ( coincStat == s3_snr_chi_stat )
      {
        REAL4 tmp_snr = snglInspiral->snr;
        REAL4 tmp_chisq = snglInspiral->chisq;
        
        statValue += tmp_snr * tmp_snr * tmp_snr * tmp_snr / 
          ( tmp_chisq * ( 250 + tmp_snr * tmp_snr ) );
      }

    }
  }

  /*    for the bitten L case only , we need to compare different 
        values and keep the minimum one */
  if ( coincStat == bitten_l )
  {
    statValue = sqrt(statValue);
    
    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if ( (snglInspiral = coincInspiral->snglInspiral[ifoNumber]) )
      {
        if (statValues[ifoNumber] < statValue)
        {
         statValue = statValues[ifoNumber];
        }
      }
    }
  }

  return( statValue );
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALClusterCoincInspiralTable (
    CoincInspiralTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralBittenLParams *bittenLParams
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALClusterCoincInspiralTable";
  CoincInspiralTable     *thisCoinc = NULL;
  CoincInspiralTable     *prevCoinc = NULL;
  CoincInspiralTable     *nextCoinc = NULL;
  int                     numCoincClust = 0;

  if ( !coincList )
  {
    XLAL_ERROR(func,XLAL_EIO);
  }

  if ( ! *coincList )
  {
    XLALPrintInfo( 
      "XLALClusterCoincInspiralTable: Empty coincList passed as input" );
    return( 0 );
  }
    
  thisCoinc = (*coincList);
  nextCoinc = (*coincList)->next;
  *coincList = NULL;

  while ( nextCoinc )
  {
    INT8 thisTime = XLALCoincInspiralTimeNS( thisCoinc );
    INT8 nextTime = XLALCoincInspiralTimeNS( nextCoinc );

    /* find events within the cluster window */
    if ( (nextTime - thisTime) < dtimeNS )
    {
      REAL4 thisStat = 
        XLALCoincInspiralStat( thisCoinc, coincStat, bittenLParams);
      REAL4 nextStat = 
        XLALCoincInspiralStat( nextCoinc, coincStat, bittenLParams );
      
      if ( nextStat > thisStat )
      {
        /* displace previous event in cluster */
        if( prevCoinc )
        {
          prevCoinc->next = nextCoinc;
        }
        XLALFreeCoincInspiral( &thisCoinc );
        thisCoinc = nextCoinc;
        nextCoinc = thisCoinc->next;
      }
      else
      {
        /* otherwise just dump next event from cluster */
        thisCoinc->next = nextCoinc->next;
        XLALFreeCoincInspiral ( &nextCoinc );
        nextCoinc = thisCoinc->next;
      }
    }
    else 
    {
      /* otherwise we keep this unique event trigger */
      if ( ! (*coincList) )
      {
        *coincList = thisCoinc;
      }
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
      nextCoinc = thisCoinc->next;
      ++numCoincClust;
    }
  }

  /* store the last event */
  if ( ! (*coincList) )
  {
    *coincList = thisCoinc;
  }
  ++numCoincClust;
  
  return(numCoincClust);
}

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALCompareCoincInspiralByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  INT8 ta, tb;

  ta = XLALCoincInspiralTimeNS ( aPtr );
  tb = XLALCoincInspiralTimeNS ( bPtr );

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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALSortCoincInspiral (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                   i;
  INT4                   numEvents = 0;
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
     XLALPrintInfo( 
      "XLALSortCoincInspiral: Empty coincInspiral passed as input" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (CoincInspiralTable **) 
    LALCalloc( numEvents, sizeof(CoincInspiralTable *) );
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

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALCoincInspiralIfos (
    CoincInspiralTable  *coincInspiral,
    char                *ifos    
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  int                   ifosMatch  = 1;
  CHAR                  ifo[LIGOMETA_IFO_MAX];
    
  if ( !coincInspiral )
  {
    return ( 0 );
  }
  
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    XLALReturnIFO( ifo, ifoNumber);

    /* check that the coinc is of the correct type */
    if ( (coincInspiral->snglInspiral[ifoNumber] &&  !strstr(ifos,ifo)) ||
        (!coincInspiral->snglInspiral[ifoNumber] &&  strstr(ifos,ifo)) )
    {
      ifosMatch = 0;
      break;
    }
  }
  return( ifosMatch );
}  
      
         
int
XLALCoincInspiralIfosCut(
    CoincInspiralTable **coincHead,
    char                *ifos    
    )
{
  CoincInspiralTable    *prevCoinc = NULL;
  CoincInspiralTable    *thisCoinc = NULL;
  int                    numCoinc = 0;
  
  thisCoinc = *coincHead; 
  *coincHead = NULL;
  
  while ( thisCoinc )
  {
    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincInspiralIfos( tmpCoinc, ifos ) )
    {
      /* ifos match so keep tmpCoinc */
      if ( ! *coincHead  )
      {
        *coincHead = tmpCoinc;
      }
      else
      {
        prevCoinc->next = tmpCoinc;
      }
      tmpCoinc->next = NULL;
      prevCoinc = tmpCoinc;
      ++numCoinc;
    }
    else
    {
      /* discard tmpCoinc */
      XLALFreeCoincInspiral( &tmpCoinc );
    }
  }
  
  return( numCoinc );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
UINT8
XLALCoincInspiralIdNumber (
    CoincInspiralTable  *coincInspiral
    )
/* </lalVerbatim> */
{
  static const char *func = "CoincInspiralIdNumber";
  SnglInspiralTable    *thisSngl = NULL;
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;

  if ( !coincInspiral )
  {
    XLAL_ERROR(func,XLAL_EIO);
  }
  
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    EventIDColumn *thisID = NULL;
    if ( (thisSngl = coincInspiral->snglInspiral[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that 
       * points to thisCoinc */
      thisID = thisSngl->event_id;

      while ( thisID )
      {
        /* test if thisID points to our coinc */ 
        if ( thisID->coincInspiralTable == coincInspiral )
        {
          return( thisID->id );
          break;
        }
      }
    }
  }
  /* should never get here */
  XLALPrintError( "Unable to find id associated to this event" );
  XLAL_ERROR(func,XLAL_EIO);
}


CoincInspiralTable *
XLALCoincInspiralSlideCut(
    CoincInspiralTable **coincHead,
    int                  slideNum    
    )
{
  CoincInspiralTable    *prevCoinc      = NULL;
  CoincInspiralTable    *thisCoinc      = NULL;
  CoincInspiralTable    *slideHead      = NULL;
  CoincInspiralTable    *thisSlideCoinc = NULL;
  
  UINT8 idNumber = 0;
  
  if( slideNum < 0 )
  {
    slideNum = 5000 - slideNum;
  }
  
  thisCoinc = *coincHead; 
  *coincHead = NULL;
  
  while ( thisCoinc )
  {
    idNumber = XLALCoincInspiralIdNumber( thisCoinc );

    if ( (int) ((idNumber % 1000000000) / 100000) == slideNum )
    {
      /* add thisCoinc to the slideCoinc list */
      if ( slideHead )
      {
        thisSlideCoinc = thisSlideCoinc->next = thisCoinc;
      }
      else
      {
        slideHead = thisSlideCoinc = thisCoinc;
      }

      /* remove from coincHead list */
      if ( prevCoinc )
      {
        prevCoinc->next = thisCoinc->next;
      }

      thisCoinc = thisCoinc->next;
      thisSlideCoinc->next = NULL;
    }
    else
    {
      /* move along the list */
      if( ! *coincHead )
      {
        *coincHead = thisCoinc;
      }
      
      prevCoinc = thisCoinc;
      thisCoinc = thisCoinc->next;
    }
  }
  return( slideHead );
}

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALInspiralSNRCutBCV2(
    CoincInspiralTable        **coincInspiral
    )
/* </lalVerbatim> */
{
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *prevCoinc = NULL;
  CoincInspiralTable   *coincHead = NULL;

  thisCoinc = *coincInspiral;
  coincHead = NULL;

  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    REAL4 snrH1 = 0, snrH2 = 0;

    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if( tmpCoinc->snglInspiral[LAL_IFO_H1] && tmpCoinc->snglInspiral[LAL_IFO_H2] 
      && tmpCoinc->snglInspiral[LAL_IFO_H1]->snr < tmpCoinc->snglInspiral[LAL_IFO_H2]->snr)
    {
       discardTrigger = 1;  
    }

    if( discardTrigger )
    {
      XLALFreeCoincInspiral( &tmpCoinc );
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
INT4 XLALCountCoincInspiral( CoincInspiralTable *head )
/* </lalVerbatim> */
{
  INT4 length;
  CoincInspiralTable *event;
 
  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
CoincInspiralTable *
XLALStatCutCoincInspiral (
    CoincInspiralTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralBittenLParams *bittenLParams,
    REAL4                       statCut
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable    *prevEvent = NULL;


  thisEvent = eventHead;
  eventHead = NULL;
  
  while ( thisEvent )
  {
    CoincInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;
    
    if ( XLALCoincInspiralStat(tmpEvent,coincStat,bittenLParams) >= statCut )
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
      XLALFreeCoincInspiral ( &tmpEvent );
    }
  }
  return( eventHead );
}


