/*
*  Copyright (C) 2007 Lisa M. Goggin
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
int
XLALCoincRingdownIfosDiscard(
    CoincRingdownTable **coincHead,
    char                *ifos
    )
/* </lalVerbatim> */
{
  CoincRingdownTable    *prevCoinc = NULL;
  CoincRingdownTable    *thisCoinc = NULL;
  int                    numCoinc = 0;

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )   {
    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincRingdownIfos( tmpCoinc, ifos ) )
    {
      /* ifos match so discard tmpCoinc */
      XLALFreeCoincRingdown( &tmpCoinc );
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
  /*static const char *func = "FreeCoincRingdown";*/
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
  /*
  CoincRingdownTable  *coincRingdown = NULL;
  EventIDColumn       *eventId = NULL;
  */

  INITSTATUS( status, "LALAddSnglRingdownToCoinc", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( coincPtr, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );
  ASSERT( snglRingdown, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  *coincPtr = XLALAddSnglRingdownToCoinc(*coincPtr, snglRingdown);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
CoincRingdownTable *
XLALAddSnglRingdownToCoinc(
    CoincRingdownTable         *coincRingdown,
    SnglRingdownTable          *snglRingdown
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALAddSnglRingdownToCoinc";
  EventIDColumn     *eventId = NULL;

  /* allocate memory for new coinc if it doesn't exist */
  if (! coincRingdown )
  {
    coincRingdown = (CoincRingdownTable *)
      LALCalloc( 1, sizeof(CoincRingdownTable) );
    if ( !coincRingdown )
    {
      LALFree( coincRingdown );
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
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
        XLALPrintError( "Invalid ifo in input snglInspiral" );
        XLAL_ERROR_NULL(func,XLAL_EIO);
      }
      break;

    case 'L':
      coincRingdown->snglRingdown[LAL_IFO_L1] = snglRingdown;
      break;

    default:
      /* Invalid Detector Site */
      XLALPrintError( "Invalid ifo in input snglInspiral" );
      XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  ++(coincRingdown->numIfos);

  /* create an eventId for the single, populate it with the single and coinc */
  if ( ! snglRingdown->event_id )
  {
    eventId = (EventIDColumn *) LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
    snglRingdown->event_id = eventId;
  }
  else
  {
     for( eventId = snglRingdown->event_id; eventId->next;
         eventId = eventId->next);
     eventId = eventId->next = (EventIDColumn *)
         LALCalloc( 1, sizeof(EventIDColumn) );
    if ( !eventId )
    {
      LALFree(eventId);
      XLAL_ERROR_NULL(func,XLAL_ENOMEM);
    }
  }
  eventId->snglRingdownTable = snglRingdown;
  eventId->coincRingdownTable = coincRingdown;

  return coincRingdown;
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
  INITSTATUS( status, "LALExtractCoincSngls", COINCRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  *snglPtr = XLALExtractSnglRingdownFromCoinc( coincRingdown, gpsStartTime,
      slideNum );


  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="CoincRingdownUtilsCP"> */
SnglRingdownTable *
XLALExtractSnglRingdownFromCoinc(
    CoincRingdownTable         *coincRingdown,
    LIGOTimeGPS                *gpsStartTime,
    INT4                        slideNum
    )
/* </lalVerbatim> */
{
  static const char *func = "ExtractSnglRingdownFromCoinc";
  SnglRingdownTable  *snglHead = NULL;
  SnglRingdownTable  *thisSngl = NULL;
  SnglRingdownTable  *thisCoincEntry = NULL;
  CoincRingdownTable *thisCoinc = NULL;
  EventIDColumn      *eventId = NULL;
  UINT4               eventNum = 1;
  INT4                j;

  if ( !coincRingdown )
  {
    XLALPrintInfo(
        "XLALExtractSnglRingdownFromCoinc: Empty coincRingdown passed as input"
        );
    return( NULL );
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
        memcpy( thisSngl, thisCoincEntry, sizeof(SnglRingdownTable) );                        thisSngl->next = NULL;


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
            XLALFreeSnglRingdown( &thisSngl );
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
        eventId->snglRingdownTable = thisSngl;
      }
    }
  }

  return( snglHead );

}


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int
XLALCoincRingdownIfos (
    CoincRingdownTable  *coincRingdown,
    char                *ifos
    )
/* </lalVerbatim> */
{
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;
  int                   ifosMatch  = 1;
  CHAR                  ifo[LIGOMETA_IFO_MAX];

  if ( !coincRingdown )
  {
    return ( 0 );
  }

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    XLALReturnIFO( ifo, ifoNumber);

    /* check that the coinc is of the correct type */
    if ( (coincRingdown->snglRingdown[ifoNumber] &&  !strstr(ifos,ifo)) ||
        (!coincRingdown->snglRingdown[ifoNumber] &&  strstr(ifos,ifo)) )
    {
      ifosMatch = 0;
      break;
    }
  }
  return( ifosMatch );
}


int
XLALCoincRingdownIfosCut(
    CoincRingdownTable **coincHead,
    char                *ifos
    )
{
  CoincRingdownTable    *prevCoinc = NULL;
  CoincRingdownTable    *thisCoinc = NULL;
  int                    numCoinc = 0;

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )
  {
    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    if ( XLALCoincRingdownIfos( tmpCoinc, ifos ) )
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
      XLALFreeCoincRingdown( &tmpCoinc );
    }
  }

  return( numCoinc );
}


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
UINT8
XLALCoincRingdownIdNumber (
    CoincRingdownTable  *coincRingdown
    )
/* </lalVerbatim> */
{
  static const char *func = "CoincRingdownIdNumber";
  SnglRingdownTable    *thisSngl = NULL;
  InterferometerNumber  ifoNumber  = LAL_UNKNOWN_IFO;

  if ( !coincRingdown )
  {
    XLAL_ERROR(func,XLAL_EIO);
  }

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    EventIDColumn *thisID = NULL;
    if ( (thisSngl = coincRingdown->snglRingdown[ifoNumber]) )
    {
      /* loop over the list of eventID's until we get to the one that
       * points to thisCoinc */
      thisID = thisSngl->event_id;

      while ( thisID )
      {
        /* test if thisID points to our coinc */
        if ( thisID->coincRingdownTable == coincRingdown )
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


CoincRingdownTable *
XLALCoincRingdownSlideCut(
    CoincRingdownTable **coincHead,
    int                  slideNum
    )
{
  CoincRingdownTable    *prevCoinc      = NULL;
  CoincRingdownTable    *thisCoinc      = NULL;
  CoincRingdownTable    *slideHead      = NULL;
  CoincRingdownTable    *thisSlideCoinc = NULL;

  UINT8 idNumber = 0;

  if( slideNum < 0 )
  {
    slideNum = 5000 - slideNum;
  }

  thisCoinc = *coincHead;
  *coincHead = NULL;

  while ( thisCoinc )
  {
    idNumber = XLALCoincRingdownIdNumber( thisCoinc );

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





/* <lalVerbatim file="CoincRingdownUtilsCP"> */
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


/* <lalVerbatim file="CoincInspiralUtilsCP"> */
CoincRingdownTable *
XLALStatCutCoincRingdown (
    CoincRingdownTable         *eventHead,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams,
    REAL4                       statCut
    )
/* </lalVerbatim> */
{
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable    *prevEvent = NULL;


  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    CoincRingdownTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( XLALCoincRingdownStat(tmpEvent,coincStat,bittenLParams) >= statCut )
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
      XLALFreeCoincRingdown ( &tmpEvent );
    }
  }
  return( eventHead );
}



/* <lalVerbatim file="CoincRingdownUtilsCP"> */
SnglRingdownTable *
XLALCompleteCoincRingdown (
    CoincRingdownTable         *eventHead,
    int                         ifoList[LAL_NUM_IFO]
    )
/* </lalVerbatim> */
{
  static const char     *func = "XLALCompleteCoincRingdown";
  CoincRingdownTable    *thisCoinc = NULL;
  SnglRingdownTable     *snglHead  = NULL;
  SnglRingdownTable     *thisSngl   = NULL;
  InterferometerNumber   ifoNumber  = LAL_UNKNOWN_IFO;
  InterferometerNumber   ifoNum  = LAL_UNKNOWN_IFO;

  for ( thisCoinc = eventHead; thisCoinc; thisCoinc = thisCoinc->next )
  {
    for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      if ( ifoList[ifoNumber] && !thisCoinc->snglRingdown[ifoNumber] )
      {
        /* we need to add a trigger for this ifo with zero snr,
         * but correct end time */
        if ( !snglHead )
        {
          snglHead = thisSngl = (SnglRingdownTable *)
              LALCalloc( 1, sizeof(SnglRingdownTable) );
        }
        else
        {
          thisSngl = thisSngl->next = (SnglRingdownTable *)
              LALCalloc( 1, sizeof(SnglRingdownTable) );
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
        while (!thisCoinc->snglRingdown[ifoNum]) ifoNum++;
        thisSngl->start_time = thisCoinc->snglRingdown[ifoNum]->start_time;

        /* add sngl to coinc */
        thisCoinc = XLALAddSnglRingdownToCoinc( thisCoinc, thisSngl );
      }
    }
  }
  return( snglHead );
}



/* <lalVerbatim file="CoincRingdownUtilsCP"> */
CoincRingdownTable *
XLALPlayTestCoincRingdown(
    CoincRingdownTable         *eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  CoincRingdownTable    *coincEventList = NULL;
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable    *prevEvent = NULL;

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
      CoincRingdownTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALCoincRingdownTimeNS( tmpEvent );
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
        XLALFreeCoincRingdown ( &tmpEvent );
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



/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int
XLALRecreateRingdownCoincFromSngls(
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
  for( thisSngl = snglRingdown; thisSngl; thisSngl = thisSngl->next )
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

#if 0

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


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
CoincRingdownTable *
XLALRingdownDistanceCut(
    CoincRingdownTable        **coincRingdown,
    RingdownAccuracyList       *accuracyParams
    )
/* </lalVerbatim> */
{
  /*static const char *func = "RingdownDistanceCut";*/
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
    REAL4 distA = 0, distB = 0;
    REAL8 sigmasqA = 0, sigmasqB = 0;

    CoincRingdownTable *tmpCoinc = thisCoinc;
    thisCoinc = thisCoinc->next;

    for ( ifoA = 0; ifoA < LAL_NUM_IFO; ifoA++ )
    {
      kappaA = accuracyParams->ifoAccuracy[ifoA].kappa;
      for ( ifoB = ifoA + 1; ifoB < LAL_NUM_IFO; ifoB++ )
      {
        kappaB = accuracyParams->ifoAccuracy[ifoB].kappa;

        if( tmpCoinc->snglRingdown[ifoA] && kappaA
            && tmpCoinc->snglRingdown[ifoB] && kappaB  )
        {
          /* perform the distance consistency test */
          sigmasqA = tmpCoinc->snglRingdown[ifoA]->sigma_sq;
          sigmasqB = tmpCoinc->snglRingdown[ifoB]->sigma_sq;
          distA = tmpCoinc->snglRingdown[ifoA]->eff_dist;
          distB = tmpCoinc->snglRingdown[ifoB]->eff_dist;

          if( ( sigmasqA > sigmasqB && distA/distB > kappaA ) ||
              ( sigmasqB > sigmasqA && distB/distA > kappaB ) )
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
  *coincRingdown = coincHead;
  return( coincHead );
}


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

INT8
XLALCoincRingdownTimeNS (
    const CoincRingdownTable         *coincRingdown
    )
{
  static const char *func = "XLALCoincRingdownTimeNS";
  InterferometerNumber  ifoNumber;
  INT8 startTime = 0;

  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( coincRingdown->snglRingdown[ifoNumber] )
    {
      startTime = XLALGPStoINT8(
          &(coincRingdown->snglRingdown[ifoNumber]->start_time) );
      return(startTime);
    }
  }
  XLAL_ERROR(func,XLAL_EIO);
}


REAL4
XLALCoincRingdownStat(
    CoincRingdownTable         *coincRingdown,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )
{
  InterferometerNumber  ifoNumber;
  SnglRingdownTable    *snglRingdown;
  REAL4                 statValues[LAL_NUM_IFO];
  REAL4 statValue = 0;
  INT4  i;
  INT4  ifoCounter = 0;

  if( coincStat == no_stat )
  {
    return(0);
  }

  /* for bittenL only*/
  if( coincStat == bitten_l )
  {
    for ( i = 0; i < LAL_NUM_IFO ; i++)
    {
      statValues[i] = 1e9; /* sufficiently high values */
    }
  }


  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
  {
    if ( (snglRingdown = coincRingdown->snglRingdown[ifoNumber]) )
    {
      /* count the number of IFOs for this coincidence */
      ifoCounter++;

      if ( coincStat == snrsq )
      {
        statValue += snglRingdown->snr * snglRingdown->snr;
      }

      else if ( coincStat == bitten_l )
      {
        statValues[ifoNumber] = bittenLParams->param_a[ifoNumber]
                * snglRingdown->snr
                + bittenLParams->param_b[ifoNumber];
        statValue += snglRingdown->snr ;
      }

    }
  }

  /*    for the bitten L case only , we need to compare different
        values and keep the minimum one */
  if ( coincStat == bitten_l )
  {
    if (coincStat == bitten_l || ifoCounter<3) {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = coincRingdown->snglRingdown[ifoNumber]) )
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


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int
XLALClusterCoincRingdownTable (
    CoincRingdownTable        **coincList,
    INT8                        dtimeNS,
    CoincInspiralStatistic      coincStat,
    CoincInspiralStatParams    *bittenLParams
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALClusterCoincRingdownTable";
  CoincRingdownTable     *thisCoinc = NULL;
  CoincRingdownTable     *prevCoinc = NULL;
  CoincRingdownTable     *nextCoinc = NULL;
  int                     numCoincClust = 0;
  REAL4 thisStat = 0;
  REAL4 nextStat = 0;
  InterferometerNumber  ifoNumber;
  SnglRingdownTable    *snglRingdown;

  CoincInspiralStatistic  tripleStat = snrsq;

  if ( !coincList )
  {
    XLAL_ERROR(func,XLAL_EIO);
  }

  if ( ! *coincList )
  {
    XLALPrintInfo(
      "XLALClusterCoincRingdownTable: Empty coincList passed as input" );
  return( 0 );
  }

  thisCoinc = (*coincList);
  nextCoinc = (*coincList)->next;
  *coincList = NULL;

  while ( nextCoinc )
  {
    INT8 thisTime = XLALCoincRingdownTimeNS( thisCoinc );
    INT8 nextTime = XLALCoincRingdownTimeNS( nextCoinc );

    /* find events within the cluster window */
    if ( (nextTime - thisTime) < dtimeNS )
    {
      if ( thisCoinc->numIfos > 2)
      {
         thisStat = XLALCoincRingdownStat( thisCoinc, tripleStat, bittenLParams );
      }
      else
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, coincStat, bittenLParams );
      }

      if ( nextCoinc->numIfos > 2)
      {
        nextStat = XLALCoincRingdownStat( nextCoinc, tripleStat, bittenLParams );
      }
      else
      {
        nextStat = XLALCoincRingdownStat( nextCoinc, coincStat, bittenLParams );
      }

      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = thisCoinc->snglRingdown[ifoNumber]) )
        {
          thisCoinc->snglRingdown[ifoNumber]->epsilon=thisStat;
        }
        if ( (snglRingdown = nextCoinc->snglRingdown[ifoNumber]) )
        {
          nextCoinc->snglRingdown[ifoNumber]->epsilon=nextStat;
        }
      }

      if ( nextStat > thisStat )
      {
        /* displace previous event in cluster */
        if( prevCoinc )
        {
          prevCoinc->next = nextCoinc;
        }
        XLALFreeCoincRingdown( &thisCoinc );
        thisCoinc = nextCoinc;
        nextCoinc = thisCoinc->next;
      }
      else
      {
        /* otherwise just dump next event from cluster */
        thisCoinc->next = nextCoinc->next;
        XLALFreeCoincRingdown ( &nextCoinc );
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
      if ( thisCoinc->numIfos > 2)
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, tripleStat, bittenLParams );
      }
      else
      {
        thisStat = XLALCoincRingdownStat( thisCoinc, coincStat, bittenLParams );
      }

      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if ( (snglRingdown = thisCoinc->snglRingdown[ifoNumber]) )
        {
          thisCoinc->snglRingdown[ifoNumber]->epsilon=thisStat;
        }
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


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
int
XLALCompareCoincRingdownByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const CoincRingdownTable *aPtr = *((const CoincRingdownTable * const *)a);
  const CoincRingdownTable *bPtr = *((const CoincRingdownTable * const *)b);
  INT8 ta, tb;

  ta = XLALCoincRingdownTimeNS ( aPtr );
  tb = XLALCoincRingdownTimeNS ( bPtr );

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


/* <lalVerbatim file="CoincRingdownUtilsCP"> */
CoincRingdownTable *
XLALSortCoincRingdown (
    CoincRingdownTable  *eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                   i;
  INT4                   numEvents = 0;
  CoincRingdownTable    *thisEvent = NULL;
  CoincRingdownTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
     XLALPrintInfo(
      "XLALSortCoincRingdown: Empty coincRingdown passed as input" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (CoincRingdownTable **)
    LALCalloc( numEvents, sizeof(CoincRingdownTable *) );
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

