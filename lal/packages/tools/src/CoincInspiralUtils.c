/*
*  Copyright (C) 2007 Stas Babak, Alexander Dietz, Drew Keppel, Gareth Jones, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>
#include <lal/CoincInspiralEllipsoid.h>

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
\idx{XLALCoincInspiralIfosDiscard()}
\idx{XLALCoincInspiralIdNumber()}
\idx{XLALCoincInspiralSlideCut()}
\idx{XLALInspiralSNRCutBCV2()}
\idx{XLALSNRCutCoincInspiral()}
\idx{XLALSortCoincInspiral()}
\idx{XLALCompareCoincInspiralByTime()}
\idx{XLALCompareCoincInspiralByEffectiveSnr()}
\idx{XLALRateCalcCoincInspiral()}
\idx{XLALRateErrorCalcCoincInspiral()}
\idx{XLALRateStatCutCoincInspiral()}
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

\texttt{LALAddSnglInspiralToCoinc()} and \texttt{XLALAddSnglInspiralToCoinc()}
add a pointer to a single inspiral table to a coinc inspiral table.  Upon
entry, if \texttt{coincPtr} points to a \texttt{NULL} coinc inspiral table, the
table is created before a pointer to the single inspiral table is added.
Additionally, an \texttt{eventId} table is created for the single inspiral
table.  This points to both the single and coinc inspirals.  If an
\texttt{eventId} already exists for the single inspiral, another eventId table
is added to the linked list.  The linked list of \texttt{eventId}s associated
to a single inspiral table allow us to easily determine which coincident events
each single is a part of.

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

\texttt{XLALCompleteCoincInspiral()} scans through a linked list of coincidence
inspirals and checks whether the coincs contain a trigger from every ifo in the
ifoList with a non-zero value.  If a trigger does not exist, it is added at the
appropriate time for the appropriate ifo, with zero snr.  The code returns a
linked list of new single inspirals which were created in the process of
completing the coincs.

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
    CHECKSTATUSPTR( status );

    /* set next trigger for comparison */
    currentTrigger[1] = currentTrigger[0]->next;
    LALGPStoINT8( status->statusPtr, &currentTriggerNS[1],
          &(currentTrigger[1]->end_time) );
    CHECKSTATUSPTR( status );

    while ( (currentTriggerNS[1] - currentTriggerNS[0]) < maxTimeDiff )
    {
      /* check that triggers pass coincidence test */
      LALCompareInspirals( status->statusPtr, currentTrigger[0],
          currentTrigger[1], accuracyParams );
      CHECKSTATUSPTR( status );

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
	BEGINFAIL (status) {
	  XLALFreeCoincInspiral( &coincHead );
	} ENDFAIL (status);
        LALAddSnglInspiralToCoinc( status->statusPtr, &thisCoinc,
            currentTrigger[1] );
	BEGINFAIL (status) {
	  XLALFreeCoincInspiral( &coincHead );
	} ENDFAIL (status);

        ++numEvents;

      }

      /* scroll on to the next sngl inspiral */

      if ( (currentTrigger[1] = currentTrigger[1]->next) )
      {
        LALGPStoINT8( status->statusPtr, &currentTriggerNS[1],
            &(currentTrigger[1]->end_time) );
	BEGINFAIL (status) {
	  XLALFreeCoincInspiral( &coincHead );
	} ENDFAIL (status);
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


  void (*func)(CoincInspiralTable *, SnglInspiralTable *, InspiralAccuracyList *) = NULL;


  INITSTATUS( status, "LALCreateNIFOCoincList", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* Choose the appropriate comparison function */
  if ( accuracyParams->test == ellipsoid )
  {
    func = XLALSnglInspiralCoincTestEllipsoid;
  }
  else
  {
    func = XLALSnglInspiralCoincTest;
  }

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
              (*func)( thisCoinc, otherCoinc->snglInspiral[ifoNumber],
                                    accuracyParams );
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
	      BEGINFAIL (status) {
		XLALFreeCoincInspiral( &nIfoCoincHead );
	      } ENDFAIL (status);

              /* add the triggers from the (N-1) coinc to the new N coinc */
              for( ifoNum = 0; ifoNum < LAL_NUM_IFO; ifoNum++ )
              {
                if( thisCoinc->snglInspiral[ifoNum] )
                {
                  LALAddSnglInspiralToCoinc( status->statusPtr, &thisNIfoCoinc,
                      thisCoinc->snglInspiral[ifoNum] );
		  BEGINFAIL (status) {
		    XLALFreeCoincInspiral( &nIfoCoincHead );
		  } ENDFAIL (status);
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
	BEGINFAIL (status) {
	  XLALFreeCoincInspiral( coincHead );
	} ENDFAIL (status);
        thisCoinc = *coincHead;
      }
      else
      {
        prevCoinc->next = thisCoinc->next;
        LALFreeCoincInspiral( status->statusPtr, &thisCoinc );
	BEGINFAIL (status) {
	  XLALFreeCoincInspiral( coincHead );
	} ENDFAIL (status);
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

  INITSTATUS( status, "LALAddSnglInspiralToCoinc", COINCINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglInspiral, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );


  *coincPtr = XLALAddSnglInspiralToCoinc(*coincPtr, snglInspiral);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALAddSnglInspiralToCoinc(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALAddSnglInspiralToCoinc";
  EventIDColumn     *eventId = NULL;

  /* allocate memory for new coinc if it doesn't exist */
  if (! coincInspiral )
  {
    coincInspiral = (CoincInspiralTable *)
      LALCalloc( 1, sizeof(CoincInspiralTable) );
    if ( !coincInspiral )
    {
      LALFree( coincInspiral );
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
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
        XLALPrintError( "Invalid ifo in input snglInspiral" );
        XLAL_ERROR_NULL(func,XLAL_EIO);
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
      XLALPrintError( "Invalid ifo in input snglInspiral" );
      XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  ++(coincInspiral->numIfos);

  /* create an eventId for the single, populate it with the single and coinc */
  if ( ! snglInspiral->event_id )
  {
    eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
    snglInspiral->event_id = eventId;
  }
  else
  {
    for( eventId = snglInspiral->event_id; eventId->next;
        eventId = eventId->next);
    eventId = eventId->next = (EventIDColumn *)
        LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
  }
  eventId->snglInspiralTable = snglInspiral;
  eventId->coincInspiralTable = coincInspiral;

  return coincInspiral;
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

  INITSTATUS( status, "LALSnglInspiralCoincTest", COINCINSPIRALUTILSC );

  XLALSnglInspiralCoincTest( coincInspiral, snglInspiral, accuracyParams );

  RETURN( status );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALSnglInspiralCoincTest(
    CoincInspiralTable         *coincInspiral,
    SnglInspiralTable          *snglInspiral,
    InspiralAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisCoincEntry;
  INT4                  match = 1;
  INT4                  ifoNumber = 0;


  /* Loop over sngl_inspirals contained in coinc_inspiral */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    thisCoincEntry = coincInspiral->snglInspiral[ifoNumber];

    if ( thisCoincEntry )
    {
      /* snglInspiral entry exists for this IFO, perform coincidence test */
      if ( ifoNumber == XLALIFONumber(snglInspiral->ifo) )
      {
        XLALPrintInfo( "We already have a coinc from this IFO" );
        accuracyParams->match = 0;
      }

      else
      {
        XLALCompareInspirals ( snglInspiral,
            thisCoincEntry, accuracyParams );
      }
      /* set match to zero if no match.  Keep same if match */
      match *= accuracyParams->match;
    }
  }
  /* returm errorParams->match to be 1 if we match, zero otherwise */
  accuracyParams->match = match;
  if ( accuracyParams->match == 0 )
    XLALPrintInfo( "Coincidence test failed" );
  if ( accuracyParams->match == 1 )
    XLALPrintInfo( "Coincidence test passed" );

  return;
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
XLALCreateCoincSlideTable(
    CoincInspiralSlideTable   **slideTableHead,
    INT4                        numSlides
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALCreateCoincSlideTable";
  CoincInspiralSlideTable  *thisSlideTable = NULL;
  CoincInspiralSlideTable  *prevSlideTable = NULL;
  INT4                      idx = 0;
  INT4                      slideNum = 0;

  for ( idx = 0; idx < 2*numSlides; idx++ )
  {
    if ( idx < numSlides )
    {
      slideNum = idx + 1;
    }
    else
    {
      slideNum = idx - 2*numSlides;
    }

    if ( *slideTableHead )
    {
      prevSlideTable = thisSlideTable;
      thisSlideTable->next = (CoincInspiralSlideTable*)
          LALCalloc( 1, sizeof(CoincInspiralSlideTable) );
      thisSlideTable = thisSlideTable->next;
    }
    else
    {
      *slideTableHead = thisSlideTable = (CoincInspiralSlideTable*)
          LALCalloc( 1, sizeof(CoincInspiralSlideTable) );
    }

    if ( !thisSlideTable )
    {
      /* out of memory: free memory + exit*/
      while ( *slideTableHead )
      {
        thisSlideTable = *slideTableHead;
        *slideTableHead = (*slideTableHead)->next;
        LALFree( thisSlideTable );
      }

      XLAL_ERROR( func,XLAL_ENOMEM );
    }

    thisSlideTable->coincInspiral = NULL;
    thisSlideTable->slideNum = slideNum;
    thisSlideTable->slideTimeAnalyzed = 0;
    thisSlideTable->currentRate = 0;
    thisSlideTable->next = NULL;
  }

  return( numSlides );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
REAL4
XLALSetupCoincSlideTable(
    CoincInspiralSlideTable    *slideTableHead,
    CoincInspiralTable         *coincSlideHead,
    char                       *timeAnalyzedFileName,
    REAL4                       timeModifier,
    INT4                        numSlides
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALCreateCoincSlideTable";
  CoincInspiralSlideTable  *thisSlideTable = NULL;
  INT4                      idx = 0;
  INT4                      slideNum = 0;
  INT4                      thisSlideNum = 0;
  FILE                     *timeAnalyzedFp = NULL;
  int                       readVal = 0;
  REAL4                     zeroLagTimeAnalyzed = 0;
  REAL4                     thisSlideTimeAnalyzed = 0;

  timeAnalyzedFp = fopen( timeAnalyzedFileName, "r" );
  readVal = fscanf( timeAnalyzedFp, "%i %f\n", &thisSlideNum,
      &zeroLagTimeAnalyzed );
  if ( readVal != 2 )
  {
    /* should never get here */
    XLALPrintError( "Error reading time analyzed file %s",
        timeAnalyzedFileName );
    XLAL_ERROR( func,XLAL_EIO );
  }

  if ( thisSlideNum )
  {
    /* should never get here */
    fclose( timeAnalyzedFp );

    XLALPrintError( "Have no analyzed time associated with zero lag" );
    XLAL_ERROR( func,XLAL_EIO );
  }

  thisSlideTable = slideTableHead;
  while ( thisSlideTable )
  {
    if ( idx < numSlides )
    {
      slideNum = idx + 1;
    }
    else
    {
      slideNum = idx - 2*numSlides;
    }
    idx++;

    readVal = fscanf( timeAnalyzedFp, "%i %f\n", &thisSlideNum,
        &thisSlideTimeAnalyzed );
    if ( readVal != 2 )
    {
      /* should never get here */
      fclose( timeAnalyzedFp );

      XLALPrintError( "Error reading time analyzed file %s",
          timeAnalyzedFileName );
      XLAL_ERROR( func,XLAL_EIO );
    }

    if ( thisSlideNum != slideNum )
    {
      /* should never get here */
      fclose( timeAnalyzedFp );

      XLALPrintError( "Have no analyzed time associated with time slide %d",
          thisSlideNum );
      XLAL_ERROR( func,XLAL_EIO );
    }

    thisSlideTable->coincInspiral = XLALCoincInspiralSlideCut(
            &coincSlideHead, slideNum);
    thisSlideTable->slideNum = slideNum;
    thisSlideTable->slideTimeAnalyzed = thisSlideTimeAnalyzed / timeModifier;
    thisSlideTable->currentRate = 0;
    thisSlideTable = thisSlideTable->next;
  }

  if ( coincSlideHead )
  {
    /* should never get here */
    fclose( timeAnalyzedFp );

    XLALPrintError(
        "Have triggers not associated with a specified time slide" );
    XLAL_ERROR( func,XLAL_EBADLEN );
  }

  readVal = fscanf( timeAnalyzedFp, "%i %f\n", &thisSlideNum,
      &thisSlideTimeAnalyzed );
  if ( readVal != EOF )
  {
    /* should never get here */
    fclose( timeAnalyzedFp );

    XLALPrintError( "Too many lines in time analyzed file %s",
        timeAnalyzedFileName );
    XLAL_ERROR( func,XLAL_EIO );
  }

  fclose( timeAnalyzedFp );

  return( zeroLagTimeAnalyzed * timeModifier );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALRecreateCoincFromSngls(
    CoincInspiralTable        **coincPtr,
    SnglInspiralTable         **snglInspiral
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
  /* InterferometerNumber  ifoInCoinc = LAL_UNKNOWN_IFO; */


  if ( !snglInspiral )
  {
    XLALPrintInfo(
      "XLALRecreateCoincFromSngls: Empty snglInspiral passed as input" );
    return( 0 );
  }

  /* sort single inspiral trigger list according to event_id */
  *snglInspiral = XLALSortSnglInspiral( *snglInspiral,
      LALCompareSnglInspiralByID);

  /* loop over the linked list of sngl inspirals */
  for( thisSngl = *snglInspiral; thisSngl; thisSngl = thisSngl->next )
  {
    ifoNumber = XLALIFONumber( thisSngl->ifo );
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
      }
    }
    else
    {
      /* need to start a new coinc */
      if ( coincHead )
      {
	prevCoinc=thisCoinc;
	thisCoinc->next =
	  LALCalloc( 1, sizeof(CoincInspiralTable) );
	thisCoinc = thisCoinc->next;
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

      eventId=thisSngl->event_id->id;
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
        snprintf( currentTrigger->ifo, LIGOMETA_IFO_MAX, ifo );
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
  REAL4 psi0A, psi0B, psi3A, psi3B, snrC, snrA, snrB, x, y, X, Y, theta=0.040;

  thisCoinc = *coincInspiral;
  coincHead = NULL;

  /* loop over the coincindent triggers */
  while( thisCoinc )
  {
    CoincInspiralTable *tmpCoinc = thisCoinc;
    discardTrigger=0;

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
            snrC = snrB;
          }
          else
          {
            snrC = snrA;
          }

          x = (psi0A - psi0B) / snrC ;
          y = (psi3A - psi3B) / snrC ;
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
    SummValueTable             *summValueList,
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


  /*  compute the iota accuracy  */
  iotaCut  =  1/((2-accuracyParams->iotaCutH1H2)
	/(2+accuracyParams->iotaCutH1H2));

  while( thisCoinc )
  {
    INT4  discardTrigger = 0;
    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if( tmpCoinc->snglInspiral[LAL_IFO_H1] && !tmpCoinc->snglInspiral[LAL_IFO_H2])
    {
      snrH1 = tmpCoinc->snglInspiral[LAL_IFO_H1]->snr;
      dH1 = 0;
      dH2 = 0;
      LALDistanceScanSummValueTable(status->statusPtr, summValueList,
				    tmpCoinc->snglInspiral[LAL_IFO_H1]->end_time, "H1",  &dH1);
      CHECKSTATUSPTR( status );

      LALDistanceScanSummValueTable(status->statusPtr, summValueList,
				    tmpCoinc->snglInspiral[LAL_IFO_H1]->end_time, "H2",  &dH2);
      CHECKSTATUSPTR( status );

      /* iota =1 */
      if (dH2/dH1*snrH1 > snrThreshold * iotaCut)
	{
	  if ( vetoSegsH2->initMagic == SEGMENTSH_INITMAGICVAL )
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
        dH1 = 0;
        dH2 = 0;
	LALDistanceScanSummValueTable(status->statusPtr, summValueList,
				      tmpCoinc->snglInspiral[LAL_IFO_H2]->end_time, "H1",  &dH1);
	CHECKSTATUSPTR( status );
	LALDistanceScanSummValueTable(status->statusPtr, summValueList,
				      tmpCoinc->snglInspiral[LAL_IFO_H2]->end_time, "H2",  &dH2);
	CHECKSTATUSPTR( status );
	/* iota = 1 */
	if (dH1/dH2*snrH2 > snrThreshold *iotaCut)
	  {
	    if ( vetoSegsH1->initMagic == SEGMENTSH_INITMAGICVAL )
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
                2*fabs(distB - distA)/(distA + distB) >
                  epsilonB/snrB + kappaB ) ||
              ( sigmasqB > sigmasqA &&
                2*fabs(distA - distB)/(distA + distB) >
                  epsilonA/snrA + kappaA ) )
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
      CHECKSTATUSPTR( status );
    }
  }

  *eventHead = eventList;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


INT8
XLALCoincInspiralTimeNS (
    const CoincInspiralTable         *coincInspiral
    )
{
  static const char *func = "XLALCoincInspiralTimeNS";
  InterferometerNumber  ifoNumber;
  INT8 endTime = 0;

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( coincInspiral->snglInspiral[ifoNumber] )
    {
      endTime = XLALGPSToINT8NS(
          &(coincInspiral->snglInspiral[ifoNumber]->end_time) );
      return(endTime);
    }
  }
  XLAL_ERROR(func,XLAL_EIO);
}

REAL4
XLALCoincInspiralStat(
    const CoincInspiralTable   *coincInspiral,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )
{
  InterferometerNumber  ifoNumber;
  SnglInspiralTable    *snglInspiral;
  REAL4                 statValues[LAL_NUM_IFO];
  REAL4 statValue = 0;
  INT4  i;
  INT4  ifoCounter = 0;
  /* This replaces the 250 in the effective snr formula. */
  REAL4 eff_snr_denom_fac = bittenLParams->eff_snr_denom_fac;

  if( coincStat == no_stat )
  {
    return(0);
  }

  /* for bittenL only*/
  if( coincStat == bitten_l || coincStat == bitten_lsq)
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
      /* count the number of IFOs for this coincidence */
      ifoCounter++;

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
          sqrt ( tmp_chisq/(2*tmp_bins-2) * (1+tmp_snr*tmp_snr/eff_snr_denom_fac) ) ;
      }
      else if ( coincStat == bitten_l || coincStat == bitten_lsq)
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
      else if ( coincStat == ifar )
      {
        statValue = 1/snglInspiral->alpha;
      }

    }
  }

  /*    for the bitten L case only , we need to compare different
        values and keep the minimum one */
  if ( coincStat == bitten_l || coincStat == bitten_lsq)
  {
    statValue = sqrt(statValue);

    if (coincStat == bitten_l || ifoCounter<3) {
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
  }

  return( statValue );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALClusterCoincInspiralTable (
    CoincInspiralTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
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

int
XLALCompareCoincInspiralByEffectiveSnr (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  REAL4 ta, tb;

  CoincInspiralStatistic coincStat = effective_snrsq;
  CoincInspiralStatParams    bittenLParams;
  memset( &bittenLParams, 0, sizeof(CoincInspiralStatParams   ) );
  /* Default value of denom fac is 250 to preserve old functionality */
  bittenLParams.eff_snr_denom_fac = 250.0;
  ta = XLALCoincInspiralStat(aPtr,coincStat,&bittenLParams);
  tb = XLALCoincInspiralStat(bPtr,coincStat,&bittenLParams);

  if ( ta > tb )
  {
    return -1;
  }
  else if ( ta < tb )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

int  XLALComputeAndStoreEffectiveSNR( CoincInspiralTable *head, CoincInspiralStatistic *stat, CoincInspiralStatParams *par)
  {
  while (head)
    {
    head->stat = XLALCoincInspiralStat(head, *stat, par);
    head = head->next;
    }
  return 0;
  }

int
XLALCompareCoincInspiralByStat (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const CoincInspiralTable *aPtr = *((const CoincInspiralTable * const *)a);
  const CoincInspiralTable *bPtr = *((const CoincInspiralTable * const *)b);
  REAL4 ta, tb;
  ta = aPtr->stat;
  tb = bPtr->stat;

  if ( ta > tb )
  {
    return -1;
  }
  else if ( ta < tb )
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALSortCoincInspiralByStat (
    CoincInspiralTable  *eventHead,
    int(*comparfunc)    (const void *, const void *),
    CoincInspiralStatParams *statParams,
    CoincInspiralStatistic *stat
    )
{
/* </lalVerbatim> */
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

  XLALComputeAndStoreEffectiveSNR( eventHead, stat, statParams);

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

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALCoincInspiralIfosCut(
    CoincInspiralTable **coincHead,
    char                *ifos
    )
/* </lalVerbatim> */
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
int
XLALCoincInspiralIfosDiscard(
    CoincInspiralTable **coincHead,
    char                *ifos
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *prevCoinc = NULL;
  CoincInspiralTable    *thisCoinc = NULL;
  int                    numCoinc = 0;

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )   {
    CoincInspiralTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincInspiralIfos( tmpCoinc, ifos ) )
    {
      /* ifos match so discard tmpCoinc */
      XLALFreeCoincInspiral( &tmpCoinc );
    }
    else
    {
      /* keep tmpCoinc */
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALStatCutCoincInspiral (
    CoincInspiralTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
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

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
int
XLALCalcExpFitNLoudestBackground (
    CoincInspiralTable         *coincSlideHead,
    int                         fitNum,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                      *fitStat,
    REAL4                      *fitA,
    REAL4                      *fitB
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALCalcExpFitAboveNLoudestBackground";

  CoincInspiralTable    *thisSlideEvent = coincSlideHead;
  int idx = 0;
  REAL4 Delta = 0;
  REAL4 X0 = 0;
  REAL4 X1 = 0;
  REAL4 X2 = 0;
  REAL4 Y0 = 0;
  REAL4 Y1 = 0;
  REAL4 thisStat = 0;

  for (idx = 1, thisSlideEvent = coincSlideHead; idx <= fitNum; idx++,
      thisSlideEvent = thisSlideEvent->next )
  {
    if ( ! thisSlideEvent )
    {
      /* should never get here */
      XLALPrintError( "Not enough Background Triggers: have %d, need %d",
          idx - 1, fitNum );
      XLAL_ERROR(func,XLAL_ERANGE);
    }

    thisStat = XLALCoincInspiralStat(thisSlideEvent, coincStat, bittenLParams);
    X0 += idx;
    X1 += thisStat * idx;
    X2 += pow(thisStat, 2.0) * idx;
    Y0 += idx * log(idx);
    Y1 += thisStat * idx * log(idx);
  }

  *fitStat = thisStat;

  Delta = X0 * X2 - X1 * X1;

  *fitA = X2 * Y0 - X1 * Y1;
  *fitA /= Delta;
  *fitA = exp(*fitA);

  *fitB = X0 * Y1 - X1 * Y0;
  *fitB /= Delta;

  return( fitNum );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
REAL4
XLALRateCalcCoincInspiral (
    CoincInspiralTable         *coincZeroHead,
    CoincInspiralTable         *coincSlideHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *thisSlideEvent = coincSlideHead;
  CoincInspiralTable    *thisEvent = coincZeroHead;

  REAL4  thisStat;
  REAL4  thisSlideStat;
  REAL4  thisRate = 0.;
  REAL4  loudestRate = -1;

  while ( thisEvent )
  {
    thisStat = XLALCoincInspiralStat(thisEvent, coincStat, bittenLParams);

    if ( (fitStat > 0) && (thisStat >= fitStat) )
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      thisRate = fitA * exp(fitB * thisStat);

      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate/timeAnalyzed;
      }
    }
    else
    {
      while ( thisSlideEvent )
      {
        thisSlideStat =
            XLALCoincInspiralStat(thisSlideEvent, coincStat, bittenLParams);

        if ( thisSlideStat >= thisStat )
        {
          /* count this slide coinc towards the rate */
          thisRate += 1.;
          thisSlideEvent = thisSlideEvent->next;
        }
        else
        {
          /* no more slide triggers with stat>=thisStat so break */
          break;
        }
      }
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
      }
    }
    }

    if ( loudestRate < 0 ) loudestRate = thisRate/timeAnalyzed;
    thisEvent = thisEvent->next;
  }

  return( loudestRate );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
REAL4
XLALRateErrorCalcCoincInspiral (
    CoincInspiralTable         *coincZeroHead,
    CoincInspiralSlideTable    *slideHeads,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    int                         numSlides,
    REAL4                       timeAnalyzed,
    REAL4                       fitStat,
    REAL4                       fitA,
    REAL4                       fitB
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALRateErrorCalcCoincInspiral";

  CoincInspiralSlideTable    *headSlideHeads = NULL;
  CoincInspiralSlideTable    *thisSlideHead = NULL;
  CoincInspiralSlideTable    *thisHeadSlideHead = NULL;
  CoincInspiralSlideTable    *tmpSlideHead = NULL;
  CoincInspiralTable         *thisEvent = coincZeroHead;

  REAL4  thisStat;
  REAL4  thisSlideStat;
  REAL4  thisRate = 0.;
  REAL4  thisRateNum = 0.;
  REAL4  thisRateDenom = 0.;
  REAL4  thisRateError = 0.;
  REAL4  loudestRate = -1;

  headSlideHeads = slideHeads;
  XLALCreateCoincSlideTable( &thisHeadSlideHead, numSlides );
  thisSlideHead = thisHeadSlideHead;

  while ( thisSlideHead && slideHeads )
  {
    thisSlideHead->coincInspiral = slideHeads->coincInspiral;
    thisSlideHead->slideTimeAnalyzed = slideHeads->slideTimeAnalyzed;
    thisSlideHead->slideNum = slideHeads->slideNum;
    thisSlideHead = thisSlideHead->next;
    slideHeads = slideHeads->next;
  }

  while ( thisEvent )
  {
    thisStat = XLALCoincInspiralStat(thisEvent, coincStat, bittenLParams);

    if ( (fitStat > 0) && (thisStat >= fitStat) )
    {
      /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
       * for all ifos in thisEvent
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      thisRate = fitA * exp(fitB * thisStat);
      thisRateError = pow(thisRate, 0.5);
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
        {
          thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
          thisEvent->snglInspiral[ifoNumber]->alpha1 =
              thisRateError;
        }
      }
    }
    else
    {
      thisRate = 0.;
      thisRateNum = 0.;
      thisRateDenom = 0.;
      thisRateError = 0.;

      thisSlideHead = thisHeadSlideHead;
      while ( thisSlideHead )
      {
        while ( thisSlideHead->coincInspiral )
        {
          thisSlideStat = XLALCoincInspiralStat(thisSlideHead->coincInspiral,
              coincStat, bittenLParams);

          if ( thisSlideStat >= thisStat )
          {
            /* count this slide coinc towards the rate */
            thisSlideHead->currentRate += 1;
            thisSlideHead->coincInspiral = thisSlideHead->coincInspiral->next;
          }
          else
          {
            /* no more slide triggers in this slide number with
             * stat>=thisStat so break */
            break;
          }
        }
        /* add this slide's current rate to thisRate */
        thisRateNum += thisSlideHead->currentRate;
        thisRateDenom += thisSlideHead->slideTimeAnalyzed;

        /* move on to the next slide number */
        thisSlideHead = thisSlideHead->next;
      }

      thisRate = timeAnalyzed * thisRateNum / thisRateDenom;

      thisSlideHead = thisHeadSlideHead;
      while ( thisSlideHead )
      {
        /* calculate error on thisRate */
        thisRateError += pow( timeAnalyzed * thisSlideHead->currentRate / \
            thisSlideHead->slideTimeAnalyzed
            - thisRate, 2.0 ) / \
            (thisRateDenom / thisSlideHead->slideTimeAnalyzed);

        /* move on to the next slide number */
        thisSlideHead = thisSlideHead->next;
      }
      thisRateError = pow(thisRateError, 0.5);

      {
        /* put thisRate in sngl_inspiral->alpha for thisEvent[ifo]
         * for all ifos in thisEvent
         * FIXME in the future */
        InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

        for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
        {
          if ( thisEvent->snglInspiral[ifoNumber] )
          {
            thisEvent->snglInspiral[ifoNumber]->alpha = thisRate;
            thisEvent->snglInspiral[ifoNumber]->alpha1 = thisRateError;
          }
        }
      }
    }

    if ( loudestRate < 0 ) loudestRate = thisRate;
    thisEvent = thisEvent->next;
  }

  if ( thisEvent )
  {
    /* should never get here */
    XLALPrintError( "Have events where FAR not calculated" );
    XLAL_ERROR(func,XLAL_EIO);
  }

  /* free the CoincInspiralSlideTable thisSlideHead */
  thisSlideHead = thisHeadSlideHead;
  while ( thisSlideHead )
  {
    tmpSlideHead = thisSlideHead;
    thisSlideHead = thisSlideHead->next;
    LALFree( tmpSlideHead );
  }

  return( loudestRate );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALRateStatCutCoincInspiral (
    CoincInspiralTable         *eventZeroHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut,
    REAL4                       rateCut
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable    *prevEvent = NULL;
  REAL4                  thisRate;
  REAL4                  thisStat;

  thisEvent = eventZeroHead;
  eventZeroHead = NULL;

  while ( thisEvent )
  {
    CoincInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;
    thisRate = 0;

    {
      /* get thisRate from sngl_inspiral->alpha for thisEvent[ifo]
       * FIXME in the future */
      InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;

      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( thisEvent->snglInspiral[ifoNumber] )
        {
          thisRate = thisEvent->snglInspiral[ifoNumber]->alpha;
        }
      }
    }
    thisStat = XLALCoincInspiralStat(tmpEvent,coincStat,bittenLParams);
    if ( thisRate > rateCut )
    {
      /* discard the remaining templates */
      XLALFreeCoincInspiral ( &tmpEvent );
      while ( thisEvent )
      {
        tmpEvent = thisEvent;
        thisEvent = thisEvent->next;
        XLALFreeCoincInspiral ( &tmpEvent );
      }
    }
    else if ( ( thisRate == rateCut ) &&
        (  thisStat < statCut ) )
    {
      /* discard the remaining templates */
      XLALFreeCoincInspiral ( &tmpEvent );
      while ( thisEvent )
      {
        tmpEvent = thisEvent;
        thisEvent = thisEvent->next;
        XLALFreeCoincInspiral ( &tmpEvent );
      }
    }
    else
    {
      /* keep this template */
      if ( ! eventZeroHead  )
      {
        eventZeroHead = tmpEvent;
      }
      else
      {
        prevEvent->next = tmpEvent;
      }
      tmpEvent->next = NULL;
      prevEvent = tmpEvent;
    }
  }
  return( eventZeroHead );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
SnglInspiralTable *
XLALCompleteCoincInspiral (
    CoincInspiralTable         *eventHead,
    int                         ifoList[LAL_NUM_IFO]
    )
/* </lalVerbatim> */
{
  static const char     *func = "XLALCompleteCoincInspiral";
  CoincInspiralTable    *thisCoinc = NULL;
  SnglInspiralTable     *snglHead  = NULL;
  SnglInspiralTable     *thisSngl   = NULL;
  InterferometerNumber   ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber   ifoNum  = LAL_UNKNOWN_IFO;

  for ( thisCoinc = eventHead; thisCoinc; thisCoinc = thisCoinc->next )
  {
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if ( ifoList[ifoNumber] && !thisCoinc->snglInspiral[ifoNumber] )
      {
        /* we need to add a trigger for this ifo with zero snr,
         * but correct end time */
        if ( !snglHead )
        {
          snglHead = thisSngl = (SnglInspiralTable *)
              LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglInspiralTable *)
              LALCalloc( 1, sizeof(SnglInspiralTable) );
        }
        /* check that the sngl was allocated successfully */
        if ( !thisSngl )
        {
          while ( snglHead )
          {
            thisSngl = snglHead;
            snglHead = snglHead->next;
            LALFree(thisSngl);
          }
          XLAL_ERROR_NULL(func,XLAL_ENOMEM);
        }

        /* populate the ifo field */
        XLALReturnIFO(thisSngl->ifo,ifoNumber);
        XLALPrintInfo( "Appending a zero snr trigger for %s\n", thisSngl->ifo);

        /* obtain the end time */
        ifoNum = 0;
        while (!thisCoinc->snglInspiral[ifoNum]) ifoNum++;
        thisSngl->end_time = thisCoinc->snglInspiral[ifoNum]->end_time;

        /* add sngl to coinc */
        thisCoinc = XLALAddSnglInspiralToCoinc( thisCoinc, thisSngl );
      }
    }
  }
  return( snglHead );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALPlayTestCoincInspiral(
    CoincInspiralTable         *eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *coincEventList = NULL;
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable    *prevEvent = NULL;

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
      CoincInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALCoincInspiralTimeNS( tmpEvent );
      isPlay = XLALINT8NanoSecIsPlayground( &triggerTime );

      if ( ( (*dataType == playground_only)  && isPlay ) ||
          ( (*dataType == exclude_play) && ! isPlay) )
      {
        /* keep this trigger */
        if ( ! coincEventList  )
        {
          coincEventList = tmpEvent;
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
        XLALFreeCoincInspiral ( &tmpEvent );
      }
    }
    eventHead = coincEventList;
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincInspiralTable *
XLALMeanMassCut(
    CoincInspiralTable         *eventHead,
    char                       *massCut,
    REAL4                      massRangeLow,
    REAL4                      massRangeHigh,
    REAL4                      mass2RangeLow,
    REAL4                      mass2RangeHigh
    )
/* </lalVerbatim> */
{
  CoincInspiralTable    *coincEventList = NULL;
  CoincInspiralTable    *thisEvent = NULL;
  CoincInspiralTable    *prevEvent = NULL;

  InterferometerNumber ifoNumber = 0;
  REAL4 meanMass = 0;
  REAL4 meanMass2 = 0;
  INT4 numIfos = 0;
  INT4 numTriggers;
  INT4 massBOOL;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  while ( thisEvent )
  {
    CoincInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;
    numIfos = 0;
    meanMass = 0;
    meanMass2 = 0;

    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if( tmpEvent->snglInspiral[ifoNumber] )
      {
        if ( ! strcmp(massCut,"mchirp") )
        {
          meanMass += tmpEvent->snglInspiral[ifoNumber]->mchirp;
        }
        else if ( ! strcmp(massCut,"mtotal") )
        {
          meanMass += tmpEvent->snglInspiral[ifoNumber]->mass1 +
                      tmpEvent->snglInspiral[ifoNumber]->mass2;
        }
        else if ( ! strcmp(massCut,"mcomp") )
        {
          meanMass += tmpEvent->snglInspiral[ifoNumber]->mass1;
          meanMass2 += tmpEvent->snglInspiral[ifoNumber]->mass2;
        }
        numIfos += 1;
      }
    }

    meanMass = meanMass/numIfos;
    meanMass2 = meanMass2/numIfos;

    if ( ! strcmp(massCut,"mcomp") )
    {
      if ( ( meanMass >= massRangeLow ) && ( meanMass < massRangeHigh ) &&
           ( meanMass2 >= mass2RangeLow ) && ( meanMass2 < mass2RangeHigh ) )
      {
        massBOOL = 1;
      }
      else
      {
        massBOOL = 0;
      }
    }
    else
    {
      if ( ( meanMass >= massRangeLow ) && ( meanMass < massRangeHigh ) )
      {
        massBOOL = 1;
      }
      else
      {
        massBOOL = 0;
      }
    }

    if ( massBOOL )
    {
      /* keep this trigger */
      if ( ! coincEventList  )
      {
        coincEventList = tmpEvent;
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
      XLALFreeCoincInspiral ( &tmpEvent );
    }
  }

  eventHead = coincEventList;
  return(eventHead);
}

/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALPopulateAccuracyParams(
       InspiralAccuracyList  *accuracyParams
)
/* </lalVerbatim> */
{

  const CHAR *func = "XLALPopulateAccuracyParams";

  INT4 ifoNumber, ifoTwo;
  LALDetector aDet, bDet;


  /* check that the accuracyParams structure is allocated */
  if ( accuracyParams == NULL )
  {
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  }

  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, ifoNumber );

    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, ifoTwo );

      /* compute maximum light travel time */
      accuracyParams->lightTravelTime[ ifoNumber][ ifoTwo ] =
	XLALLightTravelTime( &aDet, &bDet );
    }
  }

  /* set the exttrig flag */
  accuracyParams->exttrig=0;
}



/* <lalVerbatim file="CoincInspiralUtilsCP"> */
void
XLALPopulateAccuracyParamsExt(
       InspiralAccuracyList  *accuracyParams,
       const LIGOTimeGPS     *gpstime,
       const REAL8            ra_deg,
       const REAL8            dec_deg
)
/* </lalVerbatim> */
{

  const CHAR *func = "XLALPopulateAccuracyParamsExt";

  INT4 ifoNumber, ifoTwo;
  REAL8 timeDelay;
  REAL8 ra_radians, dec_radians;
  LALDetector aDet, bDet;

  /* check that the accuracyParams structure is allocated */
  if ( accuracyParams == NULL )
  {
    XLAL_ERROR_VOID( func, XLAL_EFAULT );
  }

  /* check the values given */
  if (ra_deg<0 || ra_deg > 360)
  {
    XLALPrintError("Right ascension value outside [0; 360]. Value given: %f\n", ra_deg);
    XLAL_ERROR_VOID( func, XLAL_EDATA );
  }
  if (dec_deg<-90 || dec_deg>90)
  {
    XLALPrintError("Declination value outside [-90; 90]. Value given: %f\n", dec_deg);
    XLAL_ERROR_VOID( func, XLAL_EDATA );
  }

  /* convert position */
  ra_radians  =  ra_deg * LAL_PI_180;
  dec_radians = dec_deg * LAL_PI_180;


  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, ifoNumber );

    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, ifoTwo );

      /* compute signal travel time  */
      timeDelay=-XLALArrivalTimeDiff( aDet.location, bDet.location,
				      ra_radians, dec_radians, gpstime);

      accuracyParams->lightTravelTime[ ifoNumber][ ifoTwo ] =
	(INT8) 1e9*timeDelay;
    }
  }

  /* set the exttrig flag */
  accuracyParams->exttrig=1;

}

