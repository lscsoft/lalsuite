/* <lalVerbatim file="GPStoFloatCV">
   
Author: Berukoff, S.J.  <steveb@aei-potsdam.mpg.de>
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{GPStoFloat.c}}
\label{ss:GPStoFloat.c}

Converts between \texttt{LIGOTimeGPS} and \texttt{REAL8} formats, and also
to/from \texttt{LALTimeInterval} and \texttt{REAL8} formats.


\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoFloatCP}
\idx{LALGPStoFloat()}
\idx{LALFloatToGPS()}
\idx{LALFloatToInterval()}
\idx{LALIntervalToFloat()}

\subsubsection*{Description}

This modules contains two routines, one of which converts from 
\texttt{LIGOTimeGPS} to \texttt{REAL8}, and the other, from 
\texttt{REAL8} to \texttt{LIGOTimeGPS}.  Accuracy is on par with what one
expects from a typical IEEE-compliant machine epsilon; thus, conversion 
into the \texttt{REAL8} values incurs an error of approximately 1.e-7.

\begin{itemize}
  \item \texttt{LALGPStoFloat()} converts a \texttt{LIGOTimeGPS} to \texttt{REAK8}
  \item \texttt{LALFloatToGPS()} converts a time in \texttt{REAL8} to a \texttt{LIGOTimeGPS}
  \item \texttt{LALFloatToInterval()} converts a time interval in \texttt{REAL8} to a \texttt{LALTimeInterval}
  \item \texttt{LALIntervalToFloat()} converts a time interval in \texttt{LALTimeInterval} to \texttt{REAL8}
\end{itemize}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoFloatCV}}

</lalLaTeX> */

#include <math.h>

#include "Date.h"

NRCSID(GPSTOFLOATC, "$Id$");

/* D'oh. rint(3) is not in the ANSI std. Thanks Jolien. */
#define rint(x) floor(0.5 + (x))

static const INT4 oneBillion = 1000000000;

/* <lalVerbatim file="GPStoFloatCP"> */
void
LALGPStoFloat( LALStatus         *stat,
               REAL8             *p_flt_time, /* output - floating point GPS seconds */
               const LIGOTimeGPS *p_gps_time) /* input - GPS seconds */
{  /* </lalVerbatim> */

  INT4 secs, nanosecs;

  INITSTATUS(stat, "GPStoFloat", GPSTOFLOATC);
  ATTATCHSTATUSPTR(stat);

  ASSERT(p_gps_time != NULL, stat, DATEH_ENULLINPUT,
         DATEH_MSGENULLINPUT);

  ASSERT(p_flt_time != NULL, stat, DATEH_ENULLOUTPUT,
         DATEH_MSGENULLOUTPUT);

  secs        = p_gps_time->gpsSeconds;
  nanosecs    = p_gps_time->gpsNanoSeconds;
  *p_flt_time = (REAL8)secs + (REAL8)nanosecs * 1.e-9;

  DETATCHSTATUSPTR(stat);  
  RETURN(stat);
}  /* END: GPStoFloat() */



/* <lalVerbatim file="GPStoFloatCP"> */
void
LALFloatToGPS( LALStatus   *stat,
               LIGOTimeGPS *p_gps_time,  /* output - GPS time */
               const REAL8 *p_flt_time)  /* input - floating point GPS seconds */
{  /* </lalVerbatim> */
  INT4 secs, ns;
  REAL8 inTime;

  INITSTATUS(stat, "LALFloatToGPS", GPSTOFLOATC);
  ATTATCHSTATUSPTR(stat);
  
  ASSERT(p_flt_time != NULL, stat, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT(p_gps_time != NULL, stat, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /* be generous and allow negative times */
  if ( *p_flt_time < 0)
    inTime = - *p_flt_time;
  else
    inTime = *p_flt_time;

  secs = (INT4) inTime;
  ns  = (INT4) ( ( inTime - secs ) * oneBillion + 0.5);	 /* round properly! */

  /* flip sign if time was negative */
  if ( *p_flt_time < 0)
    {
      secs *= -1;
      ns *= -1;
    }

  /* return values */
  p_gps_time->gpsSeconds     = secs;
  p_gps_time->gpsNanoSeconds = ns;
  
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}  /* END: FloatToGPS() */

  


/* <lalVerbatim file="GPStoFloatCP"> */
void LALFloatToInterval(LALStatus *status,
                        LALTimeInterval *pInterval,  /* output: deltaT in LALTimeInterval format */
                        const REAL8 *pDeltaT) /* input: time interval in floating point */
{  /* </lalVerbatim> */
  INITSTATUS(status, "LALFloatToInterval", GPSTOFLOATC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pInterval, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pDeltaT, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  pInterval->seconds     = (INT4)(*pDeltaT);
  pInterval->nanoSeconds = (INT4)rint((*pDeltaT -
                                       (REAL8)(pInterval->seconds)) *
                                      (REAL8)oneBillion);
  
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: FloatToInterval() */




/* <lalVerbatim file="GPStoFloatCP"> */
void LALIntervalToFloat(LALStatus *status,
                        REAL8 *pDeltaT,    /* output: floating point time interval */
                        const LALTimeInterval *pInterval) /* input: LALTimeInterval */
{ /* </lalVerbatim> */
  INITSTATUS(status, "LALIntervalToFloat", GPSTOFLOATC);
  ATTATCHSTATUSPTR(status);

  ASSERT(pDeltaT, status, DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);
  ASSERT(pInterval, status, DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  *pDeltaT = (REAL8)(pInterval->seconds) +
    ((REAL8)(pInterval->nanoSeconds) / (REAL8)oneBillion);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
