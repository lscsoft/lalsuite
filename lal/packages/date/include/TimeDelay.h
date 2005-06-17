/*
<lalVerbatim file="TimeDelayHV">

Author: David Chin <dwchin@umich.edu> 1-734-709-9119
$Id$
   
</lalVerbatim> */

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

</lalLaTeX> */


#ifndef _TIMEDELAY_H
#define _TIMEDELAY_H

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>

#ifdef __cplusplus
extern "C"
{
#endif

NRCSID( TIMEDELAYH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error conditions}

</lalLaTeX> */

/* <lalErrTable> */
#define TIMEDELAYH_ENUL 1

#define TIMEDELAYH_MSGENUL "Unexpected null pointer in arguments"
/* </lalErrTable> */

/* <lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
TwoDetsTimeAndASource
\end{verbatim}
\idx[Type]{TwoDetsTimeAndASource}

\noindent This structure stores two pointers to \verb+LALPlaceAndGPS+
structures, and a pointer to a \verb+SkyPosition+ structure.  The
fields are:

\begin{description}
\item{\verb+LALPlaceAndGPS *p_det_and_time1+} The first detector and GPS
\item{\verb+LALPlaceAndGPS *p_det_and_time2+} The second detector and GPS
\item{\verb+SkyPosition *p_source+} The source location (equatorial
    co\"{o}dinates in decimal radians)
\end{description}

</lalLaTeX> */

typedef struct
tagTwoDetsTimeAndASource
{
  LALPlaceAndGPS *p_det_and_time1; /* the first detector and detection time */
  LALPlaceAndGPS *p_det_and_time2; /* the second detector and detection time */
  SkyPosition    *p_source;        /* source Equatorial location
                                    * (lon=RA, lat=dec) in decimal
                                    * radians */ 
}
TwoDetsTimeAndASource;


/* <lalLaTeX>

\begin{verbatim}
DetTimeAndASource
\end{verbatim}
\idx[Type]{DetTimeAndASource}

\noindent This structure stores one pointer to a \verb+LALPlaceAndGPS+
structure, and a pointer to a \verb+SkyPosition+ structure.  The
fields are:

\begin{description}
\item{\verb+LALPlaceAndGPS *p_det_and_time+} The detector and GPS
\item{\verb+SkyPosition *p_source+}  The source location (equatorial
    co\"{o}dinates in decimal radians)
\end{description}

</lalLaTeX> */
  
typedef struct
tagDetTimeAndASource
{
  LALPlaceAndGPS *p_det_and_time; /* detector and detection time */
  SkyPosition    *p_source;       /* source Equatorial location
                                   * (lon=RA, lat=dec) in decimal
                                   * radians */ 
}
DetTimeAndASource;
    

/* <lalLaTeX>
\vfill{\footnotesize\input{TimeDelayHV}}
</lalLaTeX> */

/*
 * Function prototypes
 */

/* <lalLaTeX>
\newpage\input{TimeDelayC}
</lalLaTeX> */

void
LALTimeDelay( LALStatus                   *status,
              REAL8                       *p_delay,
              const TwoDetsTimeAndASource *p_two_detectors_time_and_source );

INT8
XLALLightTravelTime ( const LALDetector *aDet,
                      const LALDetector *bDet 
                     );


/* <lalLaTeX>
\newpage\input{TimeDelayFromEarthCenterC}
</lalLaTeX> */

void
LALTimeDelayFromEarthCenter( LALStatus               *status,
                             REAL8                   *p_delay,
                             const DetTimeAndASource *p_det_time_and_source );

#ifdef __cplusplus
}
#endif

#endif /* !defined _TIMEDELAY_H */
