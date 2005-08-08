/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterOutputVeto.c
 *
 * Author: A.C. Rodriguez
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterOutputVetoCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpFilterOutputVeto.c}}
\label{ss:FindChirpFilterOutputVeto.c}

Memory management functions for creating and destroying input data and 
workspace memory for findchirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpFilterOutputVetoCP}
\idx{LALFindChirpFilterOutputVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpFilterOutputVeto()} implements a signal
based veto.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpFilterOutputVetoCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpFilterOutputVeto.h>

double rint(double x);

NRCSID (FINDCHIRPFILTEROUTPUTVETOC, "$Id$");

/* <lalVerbatim file="FindChirpFilterOutputVetoCP"> */
void 
LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    FindChirpSegment                   *segment,
    REAL4Vector                        *chisqVec,
    REAL8                               deltaT,
    COMPLEX8Vector                     *qVec,
    REAL4                               qNorm,
    FindChirpFilterOutputVetoParams    *params
    )
/* </lalVerbatim> */
{

  UINT4                 x; /* for loops */
  REAL4		        rsqvetoWindow; /* the r^2 veto window */
  UINT4		        eventIndex; /* store the event as an index */
  UINT4		        windowIndex; /* r^2 veto window index */ 
  UINT4                 timeaboversqThresh; /* time spent above the 
                                               r^2 veto threshold */
  REAL4		        rsqvetoThresh; /* the r^2 veto threshold */	
  SnglInspiralTable    *event = NULL;
  event = *eventList;

  INITSTATUS( status, "LALFindChirpFilterFilterOutputVeto", 
      FINDCHIRPFILTEROUTPUTVETOC );
  ATTATCHSTATUSPTR( status );
 
  /*
   *    
   *   check that the arguments are reasonable
   *          
   */
 	 

  /* check that the filterOutputVeto parameter structure exist */	
  ASSERT( params->rsqvetoWindow, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
  ASSERT( params->rsqvetoThresh, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
	
  /* check that the workspace vectors exist */
  ASSERT( chisqVec->data, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL);
	
  /* check that deltaT is reasonable */	
  ASSERT( deltaT > 0, status,
      FINDCHIRPFILTEROUTPUTVETOH_EDTZO, FINDCHIRPFILTEROUTPUTVETOH_MSGEDTZO );
			
  /* make sure that the segment structure exists */
  ASSERT( segment, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL );

  /* make sure that the segment structure contains some input */
  ASSERT( segment->data->epoch.gpsSeconds, status, 
      FINDCHIRPFILTEROUTPUTVETOH_ENULL, FINDCHIRPFILTEROUTPUTVETOH_MSGENULL );
	
  
  /* If an event exists, begin the r^2 veto calculations */  
  while ( event )
  {
    /* convert the time of the event into an index */
    eventIndex = (UINT4) floor( ( (REAL8) (event->end_time.gpsSeconds - 
            segment->data->epoch.gpsSeconds) + 
          (REAL8) (event->end_time.gpsNanoSeconds + 1) * 1.0e-9 ) / deltaT );

    /* convert the window duration into an index */
				    
    windowIndex = (UINT4) floor((REAL8) (((params->rsqvetoWindow) 
                                    / 1000.0) / deltaT));
											
     /* Initialize output to zero */ 
				
     timeaboversqThresh = 0;
	     
     /* Count the number of time samples above the given r^2 veto threshold */
     for( x = eventIndex - windowIndex; x <= eventIndex ; ++x )
	{if(chisqVec->data[x] >= ((REAL4) params->rsqvetoThresh))
           ++timeaboversqThresh;
	}
		
     /* Convert the output to milliseconds and record the computed values in the 
	sngl_inspiral event */
		
     event->rsqveto_duration = (REAL4) timeaboversqThresh * deltaT;

     event = event->next;
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
