/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglBurstUtils.c
 *
 * Author: Brown, D. A.  Brady, P. R. and Ray Majumder, S. K.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglBurstUtilsCV">
Author: Brown, D. A. Brady, P. R. and Ray Majumder, S. K
$Id$
</lalVerbatim> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>

NRCSID( SNGLBURSTUTILSC, "$Id$" );

#define NANOSEC  LAL_INT8_C(1000000000)

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglBurstUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglBurstUtilsCP}
\idx{LALSortSnglBurst()}
\idx{LALCompareSnglBurstByMass()}
\idx{LALCompareSnglBurstByTime()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.
 
\vfill{\footnotesize\input{SnglBurstUtilsCV}}

</lalLaTeX>
#endif

static int ModifiedforClustering(
				 SnglBurstTable *prevEvent, 
				 SnglBurstTable *thisEvent
				 )
{
  INT8 ta1, ta2, tb1, tb2, ta_p, tb_p;
  REAL4 fa1, fa2, fb1, fb2; 
  LALStatus status; 
  INT8  epsilon = 10;


  memset( &status, 0, sizeof(LALStatus) );

  /*compute the start, stop & peak times of the prevEvent */
  LALGPStoINT8(&status, &ta1, &(prevEvent->start_time));
  ta2 = ta1 + ( NANOSEC * prevEvent->duration );
  LALGPStoINT8(&status, &ta_p, &(prevEvent->peak_time));

  /* compute the start and stop frequencies of the prevEvent */
  fa1 = prevEvent->central_freq - 0.5 * prevEvent->bandwidth;
  fa2 = fa1 + prevEvent->bandwidth;

  /*compute the start, stop & peak times of the thisEvent */
  LALGPStoINT8(&status, &tb1, &(thisEvent->start_time));
  tb2 = tb1 + ( NANOSEC * thisEvent->duration );
  LALGPStoINT8(&status, &tb_p, &(thisEvent->peak_time));

  /* compute the start and stop frequencies of the thisEvent */
  fb1 = thisEvent->central_freq - 0.5 * thisEvent->bandwidth;
  fb2 = fb1 + thisEvent->bandwidth;

  if ( llabs(ta_p - tb_p) < epsilon )
    {
      if((fb1 >= fa1 && fb1 <= fa2) || (fb2 >= fa1 && fb2 <= fa2) 
	 || (fa1 >= fb1 && fa1 <= fb2) || (fa2 >= fb1 && fa2 <= fb2))
        {
          /* set the lower & higher frequencies of the cluster */ 
          fa1 = fb1 < fa1 ? fb1 : fa1 ;
          fa2 = fb2 > fa2 ? fb2 : fa2 ;
          prevEvent->central_freq = 0.5 * (fa1 + fa2);
          prevEvent->bandwidth = (fa2 - fa1);

	  /*set the start & end times of the cluster */
          ta1 = tb1 < ta1 ? tb1 : ta1;
          ta2 = tb2 > ta2 ? tb2 : ta2 ;
          LALINT8toGPS(&status, &(prevEvent->start_time), &ta1);
          prevEvent->duration = (REAL4)(ta2 - ta1)/1.0e9;
	   
	  /* set the amplitude, snr, conf. & peak time */
	  if ( prevEvent->amplitude < thisEvent->amplitude )
            {
	      prevEvent->amplitude = thisEvent->amplitude;
	      prevEvent->snr = thisEvent->snr;
	      prevEvent->peak_time = thisEvent->peak_time;
	      prevEvent->confidence = thisEvent->confidence;
            }
	  return 1;
	}
      else 
	return 0;
    }
  else 
    return 0;
}




/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALSortSnglBurst(
    LALStatus          *status,
    SnglBurstTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglBurstTable    *thisEvent = NULL;
  SnglBurstTable   **eventHandle = NULL;

  INITSTATUS( status, "LALSortSnglBurst", SNGLBURSTUTILSC );

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
  eventHandle = (SnglBurstTable **) 
    LALCalloc( numEvents, sizeof(SnglBurstTable *) );
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



/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
LALCompareSnglBurstByTime(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  SnglBurstTable *aPtr = *((SnglBurstTable **)a);
  SnglBurstTable *bPtr = *((SnglBurstTable **)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

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


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
LALCompareSnglBurstByTimeAndFreq(
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  SnglBurstTable *aPtr = *((SnglBurstTable **)a);
  SnglBurstTable *bPtr = *((SnglBurstTable **)b);
  INT8 ta, tb;
  REAL4 flowa, flowb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

  flowa = aPtr->central_freq - 0.5 * aPtr->bandwidth;
  flowb = bPtr->central_freq - 0.5 * bPtr->bandwidth;

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta == tb && flowa > flowb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else if ( ta == tb && flowa < flowb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSnglBurst(
    LALStatus             *status,
    SnglBurstTable        *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta1, ta2, tb1, tb2;
  REAL4 fa1, fa2, fb1, fb2;
  REAL4 dm1, dm2;
  REAL4 sigmaRatio;

  INITSTATUS (status, "LALCompareSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  params->match=0;

  LALGPStoINT8( status->statusPtr, &ta1, &(aPtr->start_time) );
  LALGPStoINT8( status->statusPtr, &tb1, &(bPtr->start_time) );

  ta2 = ta1 + ( NANOSEC * aPtr->duration );
  tb2 = tb1 + ( NANOSEC * bPtr->duration );
  fa1 = (aPtr->central_freq) - 0.5*(aPtr->bandwidth);
  fa2 = (aPtr->central_freq) + 0.5*(aPtr->bandwidth);
  fb1 = (bPtr->central_freq) - 0.5*(bPtr->bandwidth);
  fb2 = (bPtr->central_freq) + 0.5*(bPtr->bandwidth);

   if(((tb1 >= ta1 && tb1 <= ta2) || (tb2 >= ta1 && tb2 <= ta2)) || ((ta1 >= tb1 && ta1 <= tb2) || (ta2 >= tb1 && ta2 <= tb2)))
     {
       if((fb1 >= fa1 && fb1 <= fa2) || (fb2 >= fa1 && fb2 <= fa2) || (fa1 >= fb1 && fa1 <= fb2) || (fa2 >= fb1 && fa2 <= fb2))
	 {
	   params->match = 1;
	 }
     }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSimBurstAndSnglBurst(
    LALStatus             *status,
    SimBurstTable         *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta, tb1, tb2;
  REAL4 fa, fb1, fb2;

  INITSTATUS(status, "LALCompareSimBurstAndSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR(status);

  LALGPStoINT8(status->statusPtr, &ta, &aPtr->geocent_peak_time);
  LALGPStoINT8(status->statusPtr, &tb1, &bPtr->start_time);

  fa = aPtr->freq;
  tb2 = tb1 + NANOSEC * bPtr->duration;
  fb1 = bPtr->central_freq - 0.5 * bPtr->bandwidth;
  fb2 = bPtr->central_freq + 0.5 * bPtr->bandwidth;

  params->match = (tb1 < ta) && (ta < tb2) && (fb1 < fa) && (fa < fb2);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALClusterSnglBurstTable (
	LALStatus        *status,
	SnglBurstTable   *burstEvent,
	INT4             *nevents
)
/* </lalVerbatim> */
{
	SnglBurstTable *thisEvent = NULL;
	SnglBurstTable *prevEvent = NULL;
	SnglBurstTable *startEvent = NULL;
	INT4 i, j;
	INT4 numModEvent = 1;

	INITSTATUS (status, "LALClusterSnglBurstTable", SNGLBURSTUTILSC);
	ATTATCHSTATUSPTR (status);

	startEvent = burstEvent;

	for(thisEvent = burstEvent->next; thisEvent; thisEvent = prevEvent->next) {
		prevEvent = startEvent;

		for(i = numModEvent; i > 0; i--) {
			if(ModifiedforClustering(prevEvent,thisEvent)) {
				for(j = i; j > 1; j--)
					prevEvent = prevEvent->next;
				prevEvent->next = thisEvent->next;
				LALFree(thisEvent);
				break;
			}
			else
				prevEvent = prevEvent->next;
		}

		if(i == 0)
			numModEvent++;
	}

	/* count the number of events in the modified list */
	for(*nevents = 1; startEvent; startEvent = startEvent->next)
		*nevents++;

	/* normal exit */
	DETATCHSTATUSPTR (status);
	RETURN (status);
}

#undef NANOSEC
