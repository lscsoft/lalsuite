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
  INT8 ta1, ta2, tb1, tb2;
  REAL4 fa1, fa2, fb1, fb2;
  REAL4 dm1, dm2;
  REAL4 sigmaRatio;

  INITSTATUS (status, "LALCompareSimBurstAndSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  params->match=0;

  LALGPStoINT8( status->statusPtr, &ta1, &(aPtr->geocent_peak_time) );
  LALGPStoINT8( status->statusPtr, &tb1, &(bPtr->start_time) );


  tb2 = tb1 + ( NANOSEC * bPtr->duration );
  fa1 = (aPtr->freq);
  fb1 = (bPtr->central_freq) - 0.5*(bPtr->bandwidth);
  fb2 = (bPtr->central_freq) + 0.5*(bPtr->bandwidth);

  if( (tb1 < ta1) && (ta1 < tb2) )
  {
     if( (fb1 < fa1) && (fa1 < fb2) )
       {
      params->match = 1;
       }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
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
  SnglBurstTable     *thisEvent=NULL,*prevEvent=NULL,*startEvent=NULL;
  INT4 i, j;
  INT4 n = 1;
  INT4 tmpnum;
  REAL4 msp_duration;

  INITSTATUS (status, "LALClusterSnglBurstTable", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  thisEvent = burstEvent->next;
  prevEvent = burstEvent;
  startEvent = burstEvent;
  msp_duration=startEvent->duration;
  while (thisEvent != NULL)
  {
    INT8 tb1, tb2;
    REAL4 fb1, fb2;

    LALGPStoINT8(status->statusPtr, &tb1, &(thisEvent->start_time));
    CHECKSTATUSPTR(status);
    tb2 = tb1 + ( NANOSEC * thisEvent->duration );
    fb1 = thisEvent->central_freq - 0.5 * thisEvent->bandwidth;
    fb2 = fb1 + thisEvent->bandwidth;

    prevEvent = startEvent;

    for (i = n; i > 0; i--)
    {   
      INT8 ta1, ta2;
      REAL4 fa1, fa2;

      /*compute the time in nanosec for each event trigger */
      LALGPStoINT8(status->statusPtr, &ta1, &(prevEvent->start_time));
      CHECKSTATUSPTR(status);
      ta2 = ta1 + ( NANOSEC * prevEvent->duration );

      /* compute the start and stop frequencies */
      fa1 = prevEvent->central_freq - 0.5 * prevEvent->bandwidth;
      fa2 = fa1 + prevEvent->bandwidth;

      /* find overlapping events */
      if (((tb1 >= ta1 && tb1 <= ta2) || (tb2 >= ta1 && tb2 <= ta2)) 
          || ((ta1 >= tb1 && ta1 <= tb2) || (ta2 >= tb1 && ta2 <= tb2)))
      {
        if((fb1 >= fa1 && fb1 <= fa2) || (fb2 >= fa1 && fb2 <= fa2) 
            || (fa1 >= fb1 && fa1 <= fb2) || (fa2 >= fb1 && fa2 <= fb2))
        {
          /* set the lower & higher frequencies of the cluster */ 
          fa1 = fb1 < fa1 ? fb1 : fa1 ;
          fa2 = fb2 > fa2 ? fb2 : fa2 ;
          prevEvent->central_freq = 0.5 * (fa1 + fa2);
          prevEvent->bandwidth = (fa2 - fa1);

          /* compare the confidence and use the most confident event
	   * to determine the peak time.   
	   * However if all but the first line is commented then it is 
	   * most likely that amplitude/snr is being considered for the
	   * determination of peak time
	   */

          if ( prevEvent->confidence > thisEvent->confidence )
          {
            prevEvent->confidence = thisEvent->confidence;
	    /* prevEvent->amplitude = thisEvent->amplitude;
            prevEvent->snr = thisEvent->snr;
            prevEvent->peak_time = thisEvent->peak_time;
            msp_duration = thisEvent->duration;*/
          }
	  /*  else if( prevEvent->confidence == thisEvent->confidence )
          {/*if equal confidence use the one with shortest duration 
            if( msp_duration > thisEvent->duration )
            {
              prevEvent->peak_time = thisEvent->peak_time;
              prevEvent->amplitude = thisEvent->amplitude;
              prevEvent->snr = thisEvent->snr;
              msp_duration = thisEvent->duration;
            }
	    }*/

          /*set the start & end times of the cluster */
          ta1 = tb1 < ta1 ? tb1 : ta1;
          ta2 = tb2 > ta2 ? tb2 : ta2 ;
          LALINT8toGPS(status->statusPtr, &(prevEvent->start_time), &ta1);
          CHECKSTATUSPTR(status);
          prevEvent->duration = (REAL4)(ta2 - ta1)/1.0e9;

          /* If amplitude/snr is used to compare instead of confidence. 
	   * However if want to use confidence then comment the lines 
	   * below and uncomment the appropriate line above.
	   */

	  if ( prevEvent->amplitude < thisEvent->amplitude )
            {
            prevEvent->amplitude = thisEvent->amplitude;
            prevEvent->snr = thisEvent->snr;
            prevEvent->peak_time = thisEvent->peak_time;
            }

          /*find the next event to compare */
          for(j=1; j < i; j++)
          {		   
            prevEvent = prevEvent->next;
          }
          prevEvent->next = thisEvent->next;
          LALFree(thisEvent);
          break;
        }
        else 
        {
          /* otherwise keep this as a unique trigger */
          prevEvent = prevEvent->next;
          if ( i==1 )
          {
            n++;
          } 
        }
      }
      else 
      {
        /* otherwise keep this as a unique trigger */  
        prevEvent = prevEvent->next;
        if ( i==1 )
        {
          n++;
        }
      }
    }
    thisEvent = prevEvent->next;    
  }


  tmpnum = 1;
  /*count the number of events in the modified list */

  while (startEvent != NULL)
  {
    startEvent = startEvent->next;
    tmpnum++;
  }

  *nevents = tmpnum;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

#undef NANOSEC
