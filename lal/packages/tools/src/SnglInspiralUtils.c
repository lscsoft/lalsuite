/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brown, D. A.
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
\idx{LALCompareSnglInspiralByTime()}

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
  REAL4 sigmaRatio;

  INITSTATUS( status, "LALCompareSnglInspiral", SNGLINSPIRALUTILSC );
  ATTATCHSTATUSPTR( status );

  params->match = 0;

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->end_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->end_time) );

  /* compate on triggger time coincidence */
  if ( labs( ta - tb ) < params->dt )
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
  else
  {
    LALInfo( status, "Triggers fails time coincidence test" );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
