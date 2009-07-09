/*
*  Copyright (C) 2007 Duncan Brown, David Chin, Jolien Creighton, Kipp Cannon, Peter Shawhan, Reinhard Prix
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

/* <lalVerbatim file="LMST1CV">
Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$
</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LMST1.c}}
\label{ss:LMST1.c}

Routines related to Local and Greenwich Mean Sidereal Time (LMST1 and
GMST1) computations.

\subsection*{Prototypes}
\vspace{0.1in}

\input{LMST1CP}
\idx{LALGPStoGMST1()}
\idx{LALGPStoLMST1()}


\subsubsection*{Description}

The routines in this module compute Mean Sidereal Time in a choice of
units: seconds, hours, degrees, or radians. LMST1 is offset from GMST1
by the longitude of the observing post.

\begin{itemize}
\item  \texttt{LALGPStoGMST1()} computes GMST1 given a GPS time.
\end{itemize}


All the routines will output GMST1 or LMST1 in the units and leap second
accuracy specified by the \texttt{pUnitsAndAcc} argument.  The sidereal
time produced is within $\tilde 1$ sidereal second of values published in
the Almanac.

\subsubsection*{Algorithms}

The basic definitions and formulae are from~\cite{esaa:1992}, Ch. 2,
Sec. 24, and also Sec. B of the Astronomical Almanac. The formula
computes GMST1 for 0h UT1.  To compute GMST1 for other times, a simple
linear interpolation is done.  Since 1984-Jan-01,
GMST has been related to UT1 as follows:
%
\begin{displaymath}
  \mathrm{GMST\>of}\>0^{h}\mathrm{UT1} = 24110^{s}.54841 +
  8640184^{s}.812866\,T_{u} + 0^{s}.093104\,T^{2}_{u} -
  6.2\times10^{-6}\,T^{3}_{u}
\end{displaymath}
%
where $T_{u} = d_{u} / 36525$, $d_{u}$ is the number of days of
Universal Time elapsed since JD 2451545.0 UT1 (January 1, 2000, at 12:00 UT1),
taking on values of $\pm 0.5, \pm 1.5, \pm 2.5,
\pm 3.5, \ldots$



\subsubsection*{Uses}

\subsubsection*{Notes}

From~\cite{esaa:1992}:

\begin{quote}
  Universal Time (UT) is the measure of time used as the basis for all
  civil time-keeping.  It conforms closely to the mean diurnal motion
  of the Sun.  The apparent diurnal motion of the Sun involves both
  the nonuniform diurnal rotation of the Earth and the motion of the
  Earth in its orbit around the Sun.  Although it would be possible to
  define a system of time measurement in terms of the hour angle of
  the Sun, such a system could never be related precisely to sidereal
  time and could not, therefore, be determined by observations of star
  transits.  As a result, Universal Time is directly related to
  sidereal time by means of a numerical formula.  It does not refer to
  the motion of the Earth and is not precisely related to the hour
  angle of the Sun.

  Universal Time at any instant can be derived from observations of
  the diurnal motion of the stars or radio sources.  The uncorrected
  observed rotational timescale, which is dependent on the place of
  observation, is designated UT0.  Correcting this timescale for the
  shift in longitude of the observing station caused by polar motion
  produces the UT1 timescale, which is independent of observing
  location, but is influenced by the slightly variable rotation of the
  Earth.
\end{quote}

A useful resource is \url{http://aa.usno.navy.mil/}.

</lalLaTeX> */


#include <lal/LALRCSID.h>

NRCSID (LMST1C, "$Id$");

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include "date_value.h"
#include <lal/XLALError.h>

#define INFOSTR_LEN 256


/*
 * Compute GMST1 in requested units given GPS time
 */
/* <lalVerbatim file="LMST1CP"> */
void
LALGPStoGMST1( LALStatus         *status,
               REAL8             *p_gmst,   /* output - GMST1 */
               const LIGOTimeGPS *p_gps,    /* input - GPS time */
               const LALMSTUnitsAndAcc *pUnitsAndAcc) /* GMST1 units and accuracy */
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALGPStoGMST1", LMST1C);
  ATTATCHSTATUSPTR(status);

  /*
   * Check pointer to input variables
   */
  ASSERT (p_gps != (LIGOTimeGPS *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);
  ASSERT (pUnitsAndAcc != (LALMSTUnitsAndAcc *)NULL, status,
          DATEH_ENULLINPUT, DATEH_MSGENULLINPUT);

  /*
   * Check pointer to output variable
   */
  ASSERT (p_gmst != (REAL8 *)NULL, status,
          DATEH_ENULLOUTPUT, DATEH_MSGENULLOUTPUT);

  XLALPrintDeprecationWarning("LALGPStoGMST1", "XLALGreenwichMeanSiderealTime");

  /*
   * Compute GMST for GPS on given date in seconds
   */
  *p_gmst = fmod(XLALGreenwichMeanSiderealTime(p_gps) / (2.0 * LAL_PI) * SECS_PER_DAY, SECS_PER_DAY);
  if(*p_gmst < 0.0)
    *p_gmst += SECS_PER_DAY;

  /*
   * Convert GMST to requested units
   */
  switch (pUnitsAndAcc->units)
    {
    case MST_SEC:
      break;
    case MST_HRS:
      *p_gmst /= (REAL8)SECS_PER_HOUR;
      break;
    case MST_DEG:
      *p_gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR;
      break;
    case MST_RAD:
      *p_gmst /= (REAL8)SECS_PER_HOUR / (REAL8)DEGS_PER_HOUR * 180. /
        (REAL8)LAL_PI;
      break;
    default:
      break;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* END: LALGPStoGMST1() */
