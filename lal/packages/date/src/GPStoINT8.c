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
\idx{XLALGPStoINT8()}
\idx{LALINT8toGPS()}
\idx{XLALINT8toGPS()}

\subsubsection*{Description}

This modules contains two LAL routines and their XLAL counterparts.  One
pair of routines converts from \texttt{LIGOTimeGPS} to \texttt{INT8}
nanoseconds, and the other, from \texttt{INT8} nanoseconds to
\texttt{LIGOTimeGPS}.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GPStoINT8CV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/Date.h>

NRCSID( GPSTOINT8C, "$Id$" );

/* <lalVerbatim file="GPStoINT8CP"> */
LIGOTimeGPS *
XLALINT8toGPS ( 
    LIGOTimeGPS *output,
    INT8 input
)
/* </lalVerbatim> */
{
  INT8 s = input / LAL_INT8_C(1000000000);
  
  if(output) {
    output->gpsSeconds = (INT4)( s );
    output->gpsNanoSeconds = (INT4)( input - LAL_INT8_C(1000000000)*s );
  }

  return( output );
}

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

  XLALINT8toGPS( output, *input );

  RETURN( status );
}

/*----------------------------------------------------------------------*/
/* <lalVerbatim file="GPStoINT8CP"> */
INT8
XLALGPStoINT8 ( 
    const LIGOTimeGPS *input
)
/* </lalVerbatim> */
{
  return( (INT8) input->gpsNanoSeconds 
    + LAL_INT8_C(1000000000) * (INT8) input->gpsSeconds );
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
  
  *output = XLALGPStoINT8( input );

  RETURN( status );
}

