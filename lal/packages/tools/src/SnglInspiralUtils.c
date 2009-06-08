/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
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
\idx{LALFreeSnglInspiral()}
\idx{XLALFreeSnglInspiral()}
\idx{LALSortSnglInspiral()}
\idx{XLALSortSnglInspiral()}
\idx{LALCompareSnglInspiralByMass()}
\idx{LALCompareSnglInspiralByPsi()}
\idx{LALCompareSnglInspiralByTime()}
\idx{LALCompareSnglInspiral()}
\idx{LALCompareInspirals()}
\idx{LALClusterSnglInspiralTable()}
\idx{XLALClusterSnglInspiralTable()}
\idx{LALTimeCutSingleInspiral()}
\idx{XLALTimeCutSingleInspiral()}
\idx{LALSNRCutSingleInspiral()}
\idx{XLALSNRCutSingleInspiral()}
\idx{XLALRsqCutSingleInspiral()}
\idx{LALBCVCVetoSingleInspiral()}
\idx{LALalphaFCutSingleInspiral()}
\idx{LALIfoCutSingleInspiral()}
\idx{XLALIfoCutSingleInspiral()}
\idx{LALIfoCountSingleInspiral()}
\idx{LALTimeSlideSingleInspiral()}
\idx{LALPlayTestSingleInspiral()}
\idx{XLALPlayTestSingleInspiral()}
\idx{LALCreateTrigBank()}
\idx{LALIncaCoincidenceTest()}
\idx{LALTamaCoincidenceTest()}
\idx{XLALMaxSnglInspiralOverIntervals(()}
\idx{XLALCountSnglInspiral()}

\subsubsection*{Description}

The function \texttt{LALFreeSnglInspiral()} and \idx{XLALFreeSnglInspiral()}
free the memory associated to a single inspiral table.  The single inspiral
table may point to a linked list of EventIDColumns.  Thus, it is necessary to
free all event ids associated with the single inspiral.

The function \texttt{LALSortSnglInspiral()} and \texttt{XLALSortSnglInspiral()}
sorts a list of single inspiral tables.  The function simply calls qsort with
the appropriate comparison function, \texttt{comparfunc}.  It then ensures that
the head of the sorted list is returned.  There then follow several comparison
functions for single inspiral tables.  \texttt{LALCompareSnglInspiralByMass ()}
first compares the \texttt{mass1} entry of the two inspiral tables, returning 1
if the first mass is larger and -1 if the second is larger.  In the case that
the \texttt{mass1} fields are equal, a similar comparsion is performed on
\texttt{mass2}.  If these also agree, 0 is returned.
\texttt{LALCompareSnglInspiralByPsi()} compares the \texttt{Psi0} and
\texttt{Psi3} fields in two single inspiral tables.  The function is analogous
to the mass comparison described above.  \texttt{LALCompareSnglInspiralByTime}
compares the end times of two single inspiral tables, returnng 1 if the first
time is larger, 0 if equal and -1 if the second time is larger.

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

\texttt{LALTimeCutSingleInspiral()} and
\texttt{XLALTimeCutSingleInspiral()}takes in a linked list of single inspiral
tables and returns only those which occur after the given \texttt{startTime}
and before the \texttt{endTime}.

\texttt{LALSNRCutSingleInspiral()} and \texttt{XLALSNRCutSingleInspiral()}
take in a linked list of single inspiral tables and returns only those
triggers which have snr values above a specific snrCut.

\texttt{XLALRsqCutSingleInspiral()} performs the R-squared veto on a linked
list of single inspiral tables.  Triggers whose snr is less than
\texttt{rsqSnrMax} and whose \texttt{rsqveto\_duration} is greater than
\texttt{rsqVetoThresh} or \texttt{(optional)} whose snr is greater than
\texttt{rsqSnrMax} and whose \texttt{rsqveto\_duration} is greater than
$\mathtt{rsqAboveSnrCoeff} \times \mathtt{snr}^{\mathtt{rsqAboveSnrPow}}$

\texttt{XLALVetoSingleInspiral()} takes in a linked list of single inspiral
tables and a list of segments and returns only those triggers which do not lie
in within the \texttt{vetoSegs}.

\texttt{LALBCVCVetoSingleInspiral()} takes in a linked list of single inspiral
tables and returns only those triggers which have alphaF/SNR values below a
specific threshold and alphaF value between alphaF-hi and alphaF-lo values.  It
is relevant for the BCVC or BCVU search only.

\texttt{LALalphaFCutSingleInspiral()} takes in a linked list of single
inspiral tables and returns only those triggers which have alphaF values below
a specific alphaFcut. It is relevant for the BCV search only.

\texttt{LALIfoCutSingleInspiral()} scans through a linked list of single
inspiral tables and returns those which are from the requested \texttt{ifo}.
On input, \texttt{eventHead} is a pointer to the head of a linked list of
single inspiral tables.  On output, this list contains only single inspirals
from the requested \texttt{ifo}.  \texttt{XLALIfoCutSingleInspiral()} works
similarly, although slightly differently.  This function returns the list of
events from the specified \texttt{ifo}, while on completion,
\texttt{eventHead} contains the list of events from \textit{other} ifos.

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

\texttt{LALPlayTestSingleInspiral()} and \idx{XLALPlayTestSingleInspiral()}
test whether single inspiral events occured in playground or non-playground
times.  It then returns the requested subset of events which occurred in the
times specified by \texttt{dataType} which must be one of
\texttt{playground\_only}, \texttt{exclude\_play} or \texttt{all\_data}.

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
  INITSTATUS( status, "LALFreeSnglInspiral", SNGLINSPIRALUTILSC );
  XLALFreeSnglInspiral( eventHead );
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
  CoincInspiralTable   *thisCoinc;
  InterferometerNumber  ifoNumber;

  while ( (eventId = (*eventHead)->event_id) )
  {
    /* free any associated event_id's */
    (*eventHead)->event_id = (*eventHead)->event_id->next;

    if( (thisCoinc = eventId->coincInspiralTable) )
    {
      /* this Sngl is still part of a coinc, set pointer to NULL */
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( *eventHead == thisCoinc->snglInspiral[ifoNumber] )
        {
          thisCoinc->snglInspiral[ifoNumber] = NULL;
        }
      }
    }
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
  INITSTATUS( status, "LALSortSnglInspiral", SNGLINSPIRALUTILSC );

  *eventHead = XLALSortSnglInspiral ( *eventHead, comparfunc );

  RETURN( status );
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALSortSnglInspiral (
    SnglInspiralTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    XLALPrintInfo(
      "XLALSortSnglInspiral: Empty SnglInspiralTable passed as input\n" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (SnglInspiralTable **)
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );
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
int
LALCompareSnglInspiralByID (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  const SnglInspiralTable *aPtr = *((const SnglInspiralTable * const *)a);
  const SnglInspiralTable *bPtr = *((const SnglInspiralTable * const *)b);

  if ( aPtr->event_id->id > bPtr->event_id->id )
  {
    return 1;
  }
  else if ( aPtr->event_id->id < bPtr->event_id->id )
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

  INITSTATUS( status, "LALCompareInspirals", SNGLINSPIRALUTILSC );

  XLALCompareInspirals( aPtr, bPtr, params );

  RETURN( status );
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALCompareInspirals (
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
  REAL4   dtau0, dtau3;
  InterferometerNumber ifoaNum,  ifobNum;
  SnglInspiralAccuracy aAcc, bAcc;

  /* static const char *func = "XLALCompareInspirals"; */

  params->match = 1;

  /* check that triggers come from different IFOs */
  if( strcmp(aPtr->ifo, bPtr->ifo) )
  {
    XLALPrintInfo( "Triggers from different IFOs\n");
    params->match = 1;
  }
  else
  {
    XLALPrintInfo( "Triggers from same IFO\n");
    params->match = 0;
    return params->match;
  }

  ifoaNum = XLALIFONumber( aPtr->ifo );
  ifobNum = XLALIFONumber( bPtr->ifo );

  ta = XLALGPStoINT8( &(aPtr->end_time) );
  tb = XLALGPStoINT8( &(bPtr->end_time) );

  /* compare on trigger time coincidence */
  aAcc = params->ifoAccuracy[ifoaNum];
  bAcc = params->ifoAccuracy[ifobNum];


  if ( params->exttrig &&
       labs( ta - tb + params->lightTravelTime[ifoaNum][ifobNum]) < (aAcc.dt + bAcc.dt) )
  {
    XLALPrintInfo( "Triggers pass time coincidence test\n" );
    params->match = 1;
  }
  else if (  !params->exttrig &&
      labs( ta - tb ) < (aAcc.dt + bAcc.dt)
      + params->lightTravelTime[ifoaNum][ifobNum])
  {
    XLALPrintInfo( "Triggers pass time coincidence test\n" );
    params->match = 1;
  }
  else
  {
    XLALPrintInfo( "Triggers fail time coincidence test\n" );
    params->match = 0;
    return params->match;
  }

  switch ( params->test )
  {
    case psi0_and_psi3:
      dpsi0 = fabs( aPtr->psi0 - bPtr->psi0 );
      dpsi3 = fabs( aPtr->psi3 - bPtr->psi3 );

      /* compare psi0 and psi3 parameters */
      if ( ( dpsi0 <= (aAcc.dpsi0 + bAcc.dpsi0) )
          && ( dpsi3 <= (aAcc.dpsi3 + bAcc.dpsi3) ))
      {
        XLALPrintInfo( "Triggers are coincident in psi0 and psi3\n" );
        params->match = 1;
      }
      else
      {
        XLALPrintInfo( "Triggers are not coincident in psi0 and psi3\n" );
        params->match = 0;
      }
      break;

    case m1_and_m2:
      dmass1 = fabs( aPtr->mass1 - bPtr->mass1 );
      dmass2 = fabs( aPtr->mass2 - bPtr->mass2 );

      /* compare mass1 and mass2 parameters */
      if ( (dmass1 <= (aAcc.dm + bAcc.dm) )
        && (dmass2 <= (aAcc.dm + bAcc.dm) ))
      {
        XLALPrintInfo( "Triggers are coincident in mass1 and mass2\n" );
        params->match = 1;
      }
      else
      {
        XLALPrintInfo( "Triggers are not coincident in mass1 and mass2\n" );
        params->match = 0;
      }
      break;

    case mchirp_and_eta:
      {
      REAL4 dmchirpTest;
      dmchirp = fabs( aPtr->mchirp - bPtr->mchirp );
      deta = fabs( aPtr->eta - bPtr->eta );

      /* compare mchirp and eta parameters */
      if (aAcc.highMass &&
      ((aPtr->mass1 + aPtr->mass2 > aAcc.highMass) ||
      (bPtr->mass1 + bPtr->mass2 > bAcc.highMass)))
        dmchirpTest = aAcc.dmchirpHi + bAcc.dmchirpHi;
      else
        dmchirpTest = aAcc.dmchirp + bAcc.dmchirp;
      if ( (dmchirp <= dmchirpTest)
            && (deta <= (aAcc.deta + bAcc.deta)) )
      {
        XLALPrintInfo( "Triggers are coincident in mchirp and eta\n" );
        params->match = 1;
      }
      else
      {
        XLALPrintInfo( "Triggers fail mchirp, eta coincidence test\n" );
        params->match = 0;
      }
      }
      break;

    case tau0_and_tau3:
      dtau0 = fabs( aPtr->tau0 - bPtr->tau0 );
      dtau3 = fabs( aPtr->tau3 - bPtr->tau3 );

      /* compare tau0 and tau3 parameters */
      if ( (dtau0 <= (aAcc.dtau0 + bAcc.dtau0) )
        && (dtau3 <= (aAcc.dtau3 + bAcc.dtau3) ))
      {
        XLALPrintInfo( "Triggers are coincident in tau0 and tau3\n" );
        params->match = 1;
      }
      else
      {
        XLALPrintInfo( "Triggers are not coincident in tau0 and tau3\n" );
        params->match = 0;
      }
      break;

    default:
      XLALPrintInfo( "error: unknown test\n" );
      params->match = 0;
      break;
  }

  return params->match;
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALClusterSnglInspiralTable (
    LALStatus                  *status,
    SnglInspiralTable         **inspiralEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALClusterSnglInspiralTable", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( inspiralEvent, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  XLALClusterSnglInspiralTable ( inspiralEvent, dtimeNS, clusterchoice );

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

REAL4
XLALSnglInspiralStat(
    SnglInspiralTable         *snglInspiral,
    SnglInspiralClusterChoice  snglStat
    )
{
  REAL4 statValue = 0;

  if ( snglStat == snr )
  {
    statValue = snglInspiral->snr;
  }
  else if ( snglStat == snrsq_over_chisq )
  {
    statValue = snglInspiral->snr * snglInspiral->snr / snglInspiral->chisq;
  }
  else
  {
    statValue = 0;
  }
  return( statValue );
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
int
XLALClusterSnglInspiralTable (
    SnglInspiralTable         **inspiralList,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  static const char *func = "XLALClusterSnglInspiralTable";
  SnglInspiralTable     *thisEvent=NULL;
  SnglInspiralTable     *prevEvent=NULL;
  SnglInspiralTable     *nextEvent=NULL;
  int                    numSnglClust = 0;

  if ( !inspiralList )
  {
    XLAL_ERROR(func,XLAL_EIO);
  }

  if ( ! *inspiralList )
  {
    XLALPrintInfo(
      "XLALClusterSnglInspiralTable: Empty coincList passed as input\n" );
    return( 0 );
  }



  thisEvent = *inspiralList;
  nextEvent = (*inspiralList)->next;
  *inspiralList = NULL;

  while ( nextEvent )
  {
    INT8 thisTime = XLALGPStoINT8( &(thisEvent->end_time) );
    INT8 nextTime = XLALGPStoINT8( &(nextEvent->end_time) );;

    /* find events within the cluster window */
    if ( (nextTime - thisTime) < dtimeNS )
    {
      REAL4 thisStat = XLALSnglInspiralStat( thisEvent, clusterchoice );
      REAL4 nextStat = XLALSnglInspiralStat( nextEvent, clusterchoice );

      if ( nextStat > thisStat )
      {
        /* displace previous event in cluster */
        if( prevEvent )
        {
          prevEvent->next = nextEvent;
        }
        XLALFreeSnglInspiral( &thisEvent );
        thisEvent = nextEvent;
        nextEvent = thisEvent->next;
      }
      else
      {
        /* otherwise just dump next event from cluster */
        thisEvent->next = nextEvent->next;
        XLALFreeSnglInspiral ( &nextEvent );
        nextEvent = thisEvent->next;
      }
    }
    else
    {
      /* otherwise we keep this unique event trigger */
      if ( ! *inspiralList )
      {
        *inspiralList = thisEvent;
      }
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
      nextEvent = thisEvent->next;
      ++numSnglClust;
    }
  }

    /* store the last event */
  if ( ! (*inspiralList) )
  {
    *inspiralList = thisEvent;
  }
  ++numSnglClust;

  return(numSnglClust);
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
  INITSTATUS( status, "LALTimeCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALTimeCutSingleInspiral( *eventHead, startTime, endTime );

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALTimeCutSingleInspiral(
    SnglInspiralTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  INT8                  startTimeNS = XLALGPStoINT8( startTime );
  INT8                  endTimeNS = XLALGPStoINT8( endTime );


  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = eventHead;

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
  eventHead = inspiralEventList;

  return (eventHead);
}


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void LALSNRCutSingleInspiral (
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       snrCut
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALSNRCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALSNRCutSingleInspiral( *eventHead, snrCut );

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALSNRCutSingleInspiral (
    SnglInspiralTable          *eventHead,
    REAL4                       snrCut
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;


  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( tmpEvent->snr >= snrCut )
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
      XLALFreeSnglInspiral ( &tmpEvent );
    }
  }
  return( eventHead );
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALRsqCutSingleInspiral (
    SnglInspiralTable          *eventHead,
    REAL4                       rsqVetoTimeThresh,
    REAL4                       rsqMaxSnr,
    REAL4                       rsqAboveSnrCoeff,
    REAL4                       rsqAboveSnrPow
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  int                   numTriggers = 0;


  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( (tmpEvent->snr <= rsqMaxSnr)
      && (tmpEvent->rsqveto_duration >= rsqVetoTimeThresh) )
    {
      /* discard this event */
      XLALFreeSnglInspiral ( &tmpEvent );
    }
    else if ( ( (tmpEvent->snr > rsqMaxSnr) && (rsqAboveSnrCoeff > 0)
      && (rsqAboveSnrPow > 0) ) && (tmpEvent->rsqveto_duration >=
      rsqAboveSnrCoeff * pow(tmpEvent->snr,rsqAboveSnrPow) ) )
    {
      /* discard this event */
      XLALFreeSnglInspiral ( &tmpEvent );
    }
    else
    {
      /* keep this event */
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
      numTriggers++;
    }
  }
  return( eventHead );
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALVetoSingleInspiral (
    SnglInspiralTable          *eventHead,
    LALSegList                 *vetoSegs,
    CHAR 			*ifo
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    /*-- Check the time of this event against the veto segment list --*/
    if ( XLALSegListSearch( vetoSegs, &(thisEvent->end_time) )
	&& (strcmp(thisEvent->ifo, ifo)==0) )
    {
      /*-- This event's end_time falls within one of the veto segments --*/
      /* discard the trigger and move to the next one */
      SnglInspiralTable    *tmpEvent = NULL;
      if ( prevEvent ) prevEvent->next = thisEvent->next;
      tmpEvent = thisEvent;
      thisEvent = thisEvent->next;
      XLALFreeSnglInspiral ( &tmpEvent );
    }
    else
    {
      /* This inspiral trigger does not fall within any veto segment */
      /* keep the trigger and increment the count of triggers */
      if ( ! eventHead ) eventHead = thisEvent;
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }
  return( eventHead );
}

void
LALBCVCVetoSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    SnglInspiralBCVCalphafCut   alphafParams
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  REAL4 alphaF;
  INT4 veto;

  INITSTATUS( status, "LALBCVCVetoSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;

    /* calculate the alphaf-value for this trigger */
    thisEvent = thisEvent->next;
    alphaF = tmpEvent->tau5 * pow( tmpEvent->f_final,(2.0/3.0) );
    veto = 0;

    /* check the alphaf-range for each trigger */
    if (strstr(tmpEvent->ifo, "H1") &&
	( (alphaF < alphafParams.h1_lo) || (alphaF > alphafParams.h1_hi ) ) )
    {
      veto = 1;
    }
    else if (strstr(tmpEvent->ifo, "H2") &&
	( (alphaF < alphafParams.h2_lo) || (alphaF > alphafParams.h2_hi ) ) )
    {
      veto = 1;
    }
    else if (strstr(tmpEvent->ifo, "L1") &&
	( (alphaF < alphafParams.l1_lo) || (alphaF > alphafParams.l1_hi ) ) )
    {
      veto = 1;
    }

    if (  (tmpEvent->psi0 < alphafParams.psi0cut))
    {
      veto =0;
    }

    if ( veto == 0 )
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
LALalphaFCutSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    REAL4                       alphaFhi,
    REAL4                       alphaFlo
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;


  INITSTATUS( status, "LALalphaFCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );


  /* Remove all the triggers that are not in alphaFlo <= alphaF <= alphaFhi */

  thisEvent = *eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( ( (tmpEvent->alpha * pow(tmpEvent->f_final,(2.0/3.0))) <= alphaFhi )
      && ( (tmpEvent->alpha * pow(tmpEvent->f_final,(2.0/3.0))) >= alphaFlo ) )
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
  SnglInspiralTable    *ifoHead   = NULL;
  SnglInspiralTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCutSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  ifoHead = XLALIfoCutSingleInspiral( eventHead, ifo );

  /* free events from other ifos */
  while ( *eventHead )
  {
    thisEvent = *eventHead;
    *eventHead = (*eventHead)->next;

    XLALFreeSnglInspiral( &thisEvent );
  }

  *eventHead = ifoHead;
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglInspiralTable *
XLALIfoCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    char                       *ifo
    )
/* </lalVerbatim> */
{
  static const char *func = "IfoCutSingleInspiral";
  SnglInspiralTable    *prevEvent   = NULL;
  SnglInspiralTable    *thisEvent   = NULL;
  SnglInspiralTable    *ifoHead     = NULL;
  SnglInspiralTable    *thisIfoTrig = NULL;

  /* check that eventHead is non-null */
  if ( ! eventHead )
  {
    XLAL_ERROR_NULL(func,XLAL_EIO);
  }

  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent  = *eventHead;
  *eventHead = NULL;

  while ( thisEvent )
  {
    if ( ! strcmp( thisEvent->ifo, ifo ) )
    {
      /* ifos match so keep this event */
      if (  ifoHead  )
      {
        thisIfoTrig = thisIfoTrig->next = thisEvent;
      }
      else
      {
        ifoHead = thisIfoTrig = thisEvent;
      }

      /* remove from eventHead list */
      if ( prevEvent )
      {
        prevEvent->next = thisEvent->next;
      }

      /* move to next event */
      thisEvent = thisEvent->next;
      /* terminate ifo list */
      thisIfoTrig->next = NULL;
    }
    else
    {
      /* move along the list */
      if ( ! *eventHead )
      {
        *eventHead = thisEvent;
      }

      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  return( ifoHead );
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
  LALTimeSlideSegList(
       LALStatus                  *status,
       LALSegList                 *seglist,
       LIGOTimeGPS                *startTime,
       LIGOTimeGPS                *endTime,
       LIGOTimeGPS                *slideTime
				 )
/* </lalVerbatim> */
{
  INT8       startTimeNS = 0;
  INT8       endTimeNS   = 0;
  INT8       lengthTimeNS= 0;
  INT8       slideNS     = 0;
  INT8       segStartNS, segEndNS;
  INT8       newStartNS, newEndNS;
  LALSeg     tmpSeg;
  LALSegList tmplist;
  INT4       m;
  UINT4      n,i;

  INITSTATUS( status, "LALTimeSlideSegList", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  /* time slide segs in a seglist by a time = slideTime, except those from the
   * instrument skipIfo which are left untouched. If you want to slide
   * all triggers, simply set skipIfo = LAL_UNKNOWN_IFO */

  /* check that seglist exist */
  if ( !seglist ) {
    DETATCHSTATUSPTR (status);
    RETURN(status);
  }

  /* Make sure the segment list has been properly initialized */
  if ( seglist->initMagic != SEGMENTSH_INITMAGICVAL ) {
    DETATCHSTATUSPTR (status);
    RETURN(status);
  }

  if ( startTime )
  {
    LALGPStoINT8( status->statusPtr, &startTimeNS, startTime );
  }

  if ( endTime )
  {
    LALGPStoINT8( status->statusPtr, &endTimeNS, endTime );
  }

  /* calculate the slide time in nanoseconds */
  LALGPStoINT8( status->statusPtr, &slideNS, slideTime);

  /* initialize segment-list */
  XLALSegListInit( &tmplist );

  /* make sure slide time is within segment time */
  m=(int)( slideNS/(endTimeNS-startTimeNS) );
  slideNS-=m*(endTimeNS-startTimeNS);

  /*  check the length of the segment list. */
  n=seglist->length;

  /* loop over the entries in seglist */
  for( i=0; i<n; i++ ) {

    /* and seg time in nanoseconds */
    LALGPStoINT8( status->statusPtr, &segStartNS, &(seglist->segs[i].start) );
    LALGPStoINT8( status->statusPtr, &segEndNS,   &(seglist->segs[i].end) );

    /* do the slide */
    segStartNS += slideNS;
    segEndNS   += slideNS;


    /* check times to lie on the ring... */
    if ( startTimeNS && endTimeNS) {
      lengthTimeNS=endTimeNS - startTimeNS;

      /* print warning if the slide-time is too large...*/
      if (slideNS>lengthTimeNS) {
	fprintf(stdout,
		"WARNING in LALTimeSlideSegList: slide-window LARGER than segment-length.\n");
      }

      /* check the different cases where a segment can be
	 slide out of the ring completely or party.
	 In the latter case, the segment has to be split */
      if ( segStartNS<startTimeNS ) {

	if ( segEndNS<startTimeNS ) {
	  segStartNS += lengthTimeNS;
	  segEndNS   += lengthTimeNS;
	} else {
	  /* split up */
	  newStartNS = segStartNS + lengthTimeNS;
	  newEndNS   = endTimeNS;
	  LALINT8toGPS( status->statusPtr, &(tmpSeg.start), &newStartNS );
	  LALINT8toGPS( status->statusPtr, &(tmpSeg.end),   &newEndNS );
	  XLALSegListAppend( &tmplist, &tmpSeg );

	  segStartNS = startTimeNS;
	}

      } else if (segEndNS>endTimeNS) {

	if ( segStartNS>endTimeNS) {
	  segStartNS -= lengthTimeNS;
	  segEndNS   -= lengthTimeNS;
	} else {
	  /* split up */
	  newStartNS = startTimeNS;
	  newEndNS   = segEndNS - lengthTimeNS;
	  LALINT8toGPS( status->statusPtr, &(tmpSeg.start), &newStartNS );
	  LALINT8toGPS( status->statusPtr, &(tmpSeg.end),   &newEndNS );
	  XLALSegListAppend( &tmplist, &tmpSeg );

	  segEndNS = endTimeNS;
	}
      }
    }

    /* restore segment */
    LALINT8toGPS( status->statusPtr, &(tmpSeg.start), &segStartNS );
    LALINT8toGPS( status->statusPtr, &(tmpSeg.end),   &segEndNS );
    XLALSegListAppend( &tmplist, &tmpSeg );
  }


  /* clear the old list */
  XLALSegListClear( seglist );

  /* loop over all segments in this list */
  for (i=0; i<tmplist.length; i++) {
    XLALSegListAppend( seglist, &(tmplist.segs[i]));
  }

  /* clear the temporary list */
  XLALSegListClear( &tmplist );

  /* sort and clean up the new list */
  XLALSegListCoalesce( seglist );

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* ======================================= */
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

    while ( startTimeNS && trigTimeNS < startTimeNS )
    {
      /* if before startTime, then wrap trigger time */
      trigTimeNS += endTimeNS - startTimeNS;
    }
    while ( endTimeNS && trigTimeNS > endTimeNS )
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
SnglInspiralTable *
XLALPlayTestSingleInspiral(
    SnglInspiralTable          *eventHead,
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

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      SnglInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALGPStoINT8( &(tmpEvent->end_time) );
      isPlay = XLALINT8NanoSecIsPlayground( &triggerTime );

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
        XLALFreeSnglInspiral ( &tmpEvent );
      }
    }
    eventHead = inspiralEventList;
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


/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALPlayTestSingleInspiral(
    LALStatus                  *status,
    SnglInspiralTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALPlayTestSingleInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALPlayTestSingleInspiral(*eventHead, dataType);

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
      LALClusterSnglInspiralTable ( status->statusPtr, &timeCoincHead,
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
SnglInspiralTable *
XLALMassCut(
    SnglInspiralTable         *eventHead,
    char                      *massCut,
    REAL4                      massRangeLow,
    REAL4                      massRangeHigh,
    REAL4                      mass2RangeLow,
    REAL4                      mass2RangeHigh
    )
/* </lalVerbatim> */
{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  REAL4 massParam;
  REAL4 mass2Param;
  INT4 numTriggers;
  INT4 massBOOL;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;
    massParam = 0;
    mass2Param = 0;

    if ( ! strcmp(massCut,"mchirp") )
    {
      massParam = tmpEvent->mchirp;
    }
    else if ( ! strcmp(massCut,"mtotal") )
    {
      massParam = tmpEvent->mass1 + tmpEvent->mass2;
    }
    else if ( ! strcmp(massCut,"mcomp") )
    {
      massParam = tmpEvent->mass1;
      mass2Param = tmpEvent->mass2;
    }

    if ( ! strcmp(massCut,"mcomp") )
    {
      if ( ( massParam >= massRangeLow ) && ( massParam < massRangeHigh ) &&
           ( mass2Param >= mass2RangeLow ) && ( mass2Param < mass2RangeHigh ) )
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
      if ( ( massParam >= massRangeLow ) && ( massParam < massRangeHigh ) )
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
      XLALFreeSnglInspiral ( &tmpEvent );
    }
  }

  eventHead = inspiralEventList;
  return(eventHead);
}

