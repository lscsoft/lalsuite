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
void LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    COMPLEX8Vector                     *qVec,
    REAL4                               qNorm,
    FindChirpFilterOutputVetoParams    *params
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALFindChirpFilterFilterOutputVeto", 
      FINDCHIRPFILTEROUTPUTVETOC );
  ATTATCHSTATUSPTR( status );

int FindChirpFilterOutputVetoEngine(
    SnglInspiralTable                  *event,
    REAL4TimeSeries                    *input,
    REAL4                               threshold,
    REAL4                               window
    )
{
  /* store the input values in sample points */
	
  UINT4 eventIndex;
  UINT4 windowIndex;

  /* store the output */
  UINT4 duration = 0;
  UINT4 crossings = 0;

  /* convert the time of the event into an index */
	
#if 0

	eventIndex = (UINT4) (event->end_time.gpsSeconds - 
	                     input->segment->data->epoch.gpsSeconds) / input->deltaT;
		
#endif

  eventIndex = 262152;

  /* convert the window duration into an index */
  
	windowIndex = (UINT4) (window / input->deltaT);

  /* compute the duration and the number of threshold crossings */
  
#if 0

  /* count the number of time samples above the given threshold  */ 
		  
			for(x = eventIndex - windowIndex; x <= eventIndex ; ++x)
				{if(input[x] >= threshold)
				 ++duration;
				}
				
	/* count the number of upward-going crossing at threshold */ 
	  
			for(x = eventIndex - windowIndex + 1; x <= eventIndex ; ++x)
				 {if(input[x-1] <= threshold 
					&&  input[x] >= threshold)
					++crossings;
				}
				
			
#endif

	/* Set temporary value for testing, convert duration units to milliseconds */

  duration = 103 * input->deltaT * 1000;
	
	/* Set temporary value for testing */
	
  crossings = 99;

  /* record the computed values in the sngl_inspiral event */
	
#if 0
  event->event_duration = (REAL4) duration ; 
  event->threshold_crossings = (REAL4) crossings; 
#endif

  /* normal exit */
	
  return 0;


  /* JC: DOESN'T DO MUCH... MAKE IT DO SOMETHING */
  eventList = NULL;
  qVec = NULL;
  qNorm = 0;
  params = NULL;
	

}

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
