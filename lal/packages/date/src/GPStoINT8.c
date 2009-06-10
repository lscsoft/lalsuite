/*
*  Copyright (C) 2007 Duncan Brown, David Chin, Jolien Creighton, Kipp Cannon
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
 * File Name: GPStoINT8.c
 *
 * Author: Brown, D. A., and Creighton, T. D.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="GPStoINT8CV">
Author: Brown D. A., and Creighton, T. D.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{GPStoINT8.c}}
\label{ss:GPStoINT8.c}

Converts between \texttt{LIGOTimeGPS} and \texttt{INT8} formats.

\subsection*{Prototypes}
\vspace{0.1in}
\input{GPStoINT8CP}
\idx{LALGPStoINT8()}
\idx{LALINT8toGPS()}

\subsubsection*{Description}

This modules contains two obsolete LAL routines.  One routine converts from
\texttt{LIGOTimeGPS} to \texttt{INT8} nanoseconds, and the other from
\texttt{INT8} nanoseconds to \texttt{LIGOTimeGPS}.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoINT8CV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/XLALError.h>

NRCSID( GPSTOINT8C, "$Id$" );

/* <lalVerbatim file="GPStoINT8CP"> */
void
LALINT8toGPS (
    LALStatus          *status,
    LIGOTimeGPS        *output,
    const INT8         *input
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALINT8toGPS", GPSTOINT8C );

  XLALPrintDeprecationWarning("LALINT8toGPS", "XLALINT8NSToGPS");

  XLALINT8NSToGPS( output, *input );

  RETURN( status );
}

/*----------------------------------------------------------------------*/
/* <lalVerbatim file="GPStoINT8CP"> */
void
LALGPStoINT8 (
    LALStatus          *status,
    INT8               *output,
    const LIGOTimeGPS  *input
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALGPStoINT8", GPSTOINT8C );

  XLALPrintDeprecationWarning("LALGPStoINT8", "XLALGPSToINT8NS");

  *output = XLALGPSToINT8NS( input );

  RETURN( status );
}

