/******************************** <lalVerbatim file="GetEarthTimesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{GetEarthTimes.c}}
\label{ss:GetEarthTimes.c}

Computes the next sidereal midnight and autumnal equinox.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GetEarthTimesCP}
\index{\texttt{LALGetEarthTimes()}}

\subsubsection*{Description}

This function takes a GPS time from the parameter field
\verb@times->epoch@ and uses it to assign the fields
\verb@times->tAutumn@ and \verb@times->tMidnight@, which are
\verb@REAL8@ representations of the time in seconds from
\verb@times->epoch@ to the next autumnal equinox or sidereal midnight,
respectively.  This routine was written under the \verb@PulsarTimes.h@
header because these quantities are vital for performing pulsar
timing: they characterize the Earth's orbital and rotational phase,
and hence the Dopple modulation on an incoming signal.  See the
\verb@PulsarTimes.h@ header for more information about the
\verb@PulsarTimesParamStruc@ structure.

\subsubsection*{Algorithm}

At present this function is just a stub.  It assumes that the zero of
GPS time was sidereal midnight on the autumnal equinox, which is
clearly not the case.

When assigning the fields of \verb@*times@, it is up to the user to
choose a \verb@times->epoch@ that is close to the actual times that
are being considered.  This is important, since many computations use
a \verb@REAL8@ time variable whose origin is the time
\verb@times->epoch@.  If this is too far from the times of interest,
the \verb@REAL8@ time variables may suffer loss of precision.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GetEarthTimesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/PulsarTimes.h>

NRCSID(GETEARTHTIMESC,"$Id$");

/* <lalVerbatim file="GetEarthTimesCP"> */
void
LALGetEarthTimes( LALStatus *stat, PulsarTimesParamStruc *times )
{ /* </lalVerbatim> */
  REAL8 t; /* The GPS time as a floating-point number, in s. */

  INITSTATUS(stat,"GetEarthTimes",GETEARTHTIMESC);

  /* Make sure the parameters exist. */
  ASSERT(times,stat,PULSARTIMESH_ENUL,PULSARTIMESH_MSGENUL);
  t=(REAL8)times->epoch.gpsSeconds+
    (1.0e-9)*times->epoch.gpsNanoSeconds;

  /* I don't have the actual ephemeris, so for now, assume that GPS
     zero time is midnight on the autumnal equinox. */
  times->tAutumn=fmod(t,LAL_TWOPI*LAL_YRSID_SI);
  times->tMidnight=fmod(t,LAL_TWOPI*LAL_DAYSID_SI);

  RETURN(stat);
}
