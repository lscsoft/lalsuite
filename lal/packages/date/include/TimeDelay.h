/*
<lalVerbatim file="TimeDelayHV">

Author: Chin, David <dwchin@umich.edu> 1-734-730-1274
$Id$
   
</lalVerbatim> 
*/

/* 
<lalLaTeX>

\section{Header \texttt{TimeDelay.h}}
\label{s:TimeDelay.h}

Provides routine to compute time delay between two detectors.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/TimeDelay.h>
\end{verbatim}

This header provides prototypes of routines to compute the difference in
time for a signal to arrive at two detectors.  The routine is a direct
translation of the Maple worksheet by Anderson, \emph{et al.}, available at
\verb+http://dirac.utb.edu/~warren/unprot/beam_patterns.tar.gz+.

</lalLaTeX> 
*/


#ifndef _TIMEDELAY_H
#define _TIMEDELAY_H

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>

#ifdef __cplusplus
extern "C"
{
#endif

NRCSID( TIMEDELAYH, "$Id$" );

/*
<lalLaTeX>

\subsection*{Error conditions}

</lalLaTeX>
*/

/*
<lalErrTable>
*/
#define TIMEDELAYH_ENUL 1

#define TIMEDELAYH_MSGENUL "Unexpected null pointer in arguments"
/*
</lalErrTable>
*/

/*
<lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
TwoDetectors
\end{verbatim}
\index{\texttt{TwoDetectors}}

\noindent This structure stores two \verb@LALDetector@ structures.  The
fields are:
\begin{description}
\item[\texttt{LALPlaceAndGPS *detector1}] The first detector
\item[\texttt{LALPlaceAndGPS *detector2}] The second detector
\item[\texttt{SkyPosition *source}] The source location (equatorial
    co\"{o}dinates in decimal radians
\end{description}

</lalLaTeX>
*/

typedef struct
tagTwoDetsTimeAndASource
{
  LALPlaceAndGPS *det_and_time1; /* the first detector and detection time */
  LALPlaceAndGPS *det_and_time2; /* the second detector and detection time */
  SkyPosition    *source;        /* source Equatorial location (lon=RA, lat=dec)
                                    in decimal radians */
}
TwoDetsTimeAndASource;

/*
<lalLaTeX>
\vfill{\footnotesize\input{TimeDelayHV}}
</lalLaTeX>
*/

/*
 * Function prototypes
 */

/*
<lalLaTeX>
\newpage\input{TimeDelayC}
</lalLaTeX>
*/

void
LALTimeDelay( LALStatus                   *stat,
              REAL8                       *delay,
              const TwoDetsTimeAndASource *two_detectors_time_and_source );

#ifdef __cplusplus
}
#endif

#endif /* !defined _TIMEDELAY_H */
