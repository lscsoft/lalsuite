/* <lalVerbatim file="GPStoFloatCV">
   
Author: Berukoff, S.J.  <steveb@aei-potsdam.mpg.de>
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{GPStoFloat.c}}
\label{ss:GPStoFloat.c}

Converts between \texttt{LIGOTimeGPS} and \texttt{REAL8} formats.


\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoFloatCP}
\idx{GPStoFloat()}

\subsubsection*{Description}

This modules contains two routines, one of which converts from 
\texttt{LIGOTimeGPS} to \texttt{REAL8}, and the other, from 
\texttt{REAL8} to \texttt{LIGOTimeGPS}.  Accuracy is on par with what one
expects from a typical IEEE-compliant machine epsilon; thus, conversion 
into the \texttt{REAL8} values incurs an error of approximately 1.e-7.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoFloatCV}}

</lalLaTeX> */

#include <math.h>

#include "Date.h"

NRCSID(GPSTOFLOATC, "$Id$");

/* <lalVerbatim file="GPStoFloatCP"> */
void
LALGPStoFloat( LALStatus   *stat,
               REAL8       *p_flt_time,
               LIGOTimeGPS *p_gps_time)
{  /* </lalVerbatim> */

  INT4 secs, nanosecs;

  INITSTATUS(stat, "GPStoFloat", GPSTOFLOATC);

  ASSERT(p_gps_time != NULL, stat, DATEH_ENULLINPUT,
         DATEH_MSGENULLINPUT);

  ASSERT(p_flt_time != NULL, stat, DATEH_ENULLOUTPUT,
         DATEH_MSGENULLOUTPUT);

  secs        = p_gps_time->gpsSeconds;
  nanosecs    = p_gps_time->gpsNanoSeconds;
  *p_flt_time = (REAL8)secs + (REAL8)nanosecs * 1.e-9;

  RETURN(stat);
}  /* END: GPStoFloat() */



/* <lalVerbatim file="GPStoFloatCP"> */
void
LALFloatToGPS( LALStatus   *stat,
               LIGOTimeGPS *p_gps_time,
               REAL8       *p_flt_time)
{  /* </lalVerbatim> */

  /* Don't blame me for these obtuse variable names */

  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;
  
  INITSTATUS(stat, "GPSFloatConversion", GPSTOFLOATC);
  
  ASSERT(p_flt_time != NULL, stat, DATEH_ENULLINPUT,
         DATEH_MSGENULLINPUT);

  temp0 = floor(*p_flt_time);     /* this is tgps.S */
  temp1 = (*p_flt_time) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2); 
  temp4 = (temp2-temp3) * 0.1;

  p_gps_time->gpsSeconds     = (INT4)temp0;
  p_gps_time->gpsNanoSeconds = (INT4)temp4;
  
  RETURN(stat);
}  /* END: FloatToGPS() */

  
