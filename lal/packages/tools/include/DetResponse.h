/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Teviet Creighton
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

/* <lalVerbatim file="DetResponseHV">

Author: David Chin <dwchin@umich.edu> +1-734-709-9119
$Id$

</lalVerbatim>*/

/*
<lalLaTeX>

\section{Header \texttt{DetResponse.h}}
\label{sec:DetResponse.h}

Provides routines to compute gravitational wave detector response to
polarized planar gravitational wave originating from a given source,
detected at a given time.


\subsection{Synopsis}
\label{ss:Synopsis}

\begin{verbatim}
#include <lal/DetResponse.h>
\end{verbatim}

</lalLaTeX>
*/

#ifndef _DETRESPONSE_H
#define _DETRESPONSE_H

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>

#ifdef __cplusplus
extern "C"
{
#endif

NRCSID( DETRESPONSEH, "$Id$" );

/*
<lalLaTeX>
\subsection*{Error conditions}
</lalLaTeX>
*/

/*
<lalErrTable>
*/
#define DETRESPONSEH_ENULLINPUT  1
#define DETRESPONSEH_ENULLOUTPUT 2
#define DETRESPONSEH_ESRCNOTEQUATORIAL 3

#define DETRESPONSEH_MSGENULLINPUT "Input is NULL"
#define DETRESPONSEH_MSGENULLOUTPUT "Output is NULL"
#define DETRESPONSEH_MSGESRCNOTEQUATORIAL "Source coordinates not in Equatorial system"

/*
</lalErrTable>
*/

/*
<lalLaTeX>
\subsection*{Types and Structures}
</lalLaTeX>
*/

/* <lalLaTeX>

\vfill{\footnotesize\input{DetResponseHV}}

</lalLaTeX> */

/* <lalLaTeX>
\subsubsection*{Structure \texttt{LALSource}}
\idx[Type]{LALSource}

This structure contains gravitational wave source position (in Equatorial
co\"{o}rdinates), and orientation angle.  The orientation is measured
counter-clockwise with respect to the ``line of ascending nodes'',
\textit{i.e.} counter-clockwise with respect to a line perpendicular to the
source's meridian and extending westwards.  For a source in the Northern
celestial hemisphere, and an observer in the Northern hemisphere standing
such that they are facing South, this angle is measured counter-clockwise
from a 3 o'clock position (pointing West) at the source.  The polarization
convention is chosen such that if the source orientation were zero, the
source would be a pure $+$-polarized source.

The fields are:
\begin{description}
\item[\texttt{CHAR *name}]  Name of source
\item[\texttt{SkyPosition equatorialCoords}] Equatorial co\"{o}rdinates of
  source
\item[\texttt{REAL8 orientation}] Orientation angle ($\psi$) of source:
  counter-clockwise angle $x$-axis makes with a line perpendicular to
  meridian of source in Westward direction (\textit{i.e.} North of West),
  in decimal radians.
\end{description}

</lalLaTeX> */

typedef struct
tagLALSource
{
  CHAR         name[LALNameLength];  /* name of source, e.g. catalog number */
  SkyPosition  equatorialCoords;     /* equatorial coordinates of source,
                                        in decimal RADIANS */
  REAL8        orientation;          /* counter-clockwise angle x-axis makes
                                        with a line perpendicular to meridian
                                        of object in Westward direction, in
                                        decimal RADIANS (i.e. North of West) */
}
LALSource;

/* <lalLaTeX>
\subsubsection*{Structure \texttt{LALDetAndSource}}
\idx[Type]{LALDetAndSource}

This structure aggregates a pointer to a \texttt{LALDetector} and a
\texttt{LALSource}.  Its sole function is to allow the user to pass
detector and source parameters to the functions
\texttt{LALComputeDetAMResponse()} and
\texttt{LALComputeDetAMResponseSeries()}.

The fields are:
\begin{description}
\item[\texttt{LALDetector *pDetector}] Pointer to \texttt{LALDetector}
  object containing information about the detector
\item[\texttt{LALSource *pSource}] Pointer to \texttt{LALSource} object
  containing information about the source
\end{description}
</lalLaTeX> */

typedef struct
tagLALDetAndSource
{
  LALDetector  *pDetector;
  LALSource    *pSource;
}
LALDetAndSource;

/* <lalLaTeX>
\subsubsection{Structure \texttt{LALDetAMResponse}}
\idx[Type]{LALDetAMResponse}

This structure encapsulates the detector AM (beam pattern) coefficients for
one source at one instance in time. The fields are:

\begin{description}
\item[\texttt{REAL4 plus}] Detector response to $+$-polarized gravitational
  radiation
\item[\texttt{REAL4 cross}] Detector response to $\times$-polarized
  gravitational radiation
\item[\texttt{REAL4 scalar}] Detector response to scalar gravitational
  radiation (NB: ignored at present -- scalar response computation is not
  yet implemented)
\end{description}

</lalLaTeX> */

typedef struct
tagLALDetAMResponse
{
  REAL4 plus;
  REAL4 cross;
  REAL4 scalar;
}
LALDetAMResponse;

/* <lalLaTeX>
\subsubsection{Structure \texttt{LALDetAMResponseSeries}}
\idx[Type]{LALDetAMResponseSeries}

This structure aggregates together three \texttt{REAL4Vector}s containing
time series of detector AM response.  Since these quantities are
dimensionless, they cannot be accurately stored in a \texttt{TimeSeries}
structure.  However, \texttt{REAL4Vector}s may be conveniently converted to
\texttt{TimeSeries} format.

\begin{description}
\item[\texttt{REAL4TimeSeries *pPlus}] Pointer to a \texttt{REAL4TimeSeries}
  containing detector response to $+$-polarized gravitational radiation
  over a span of time
\item[\texttt{REAL4TimeSeries *pCross}] Pointer to a \texttt{REAL4TimeSeries}
  containing detector response to $\times$-polarized gravitational radiation
  over a span of time
\item[\texttt{REAL4TimeSeries *pScalar}] Pointer to a \texttt{REAL4TimeSeries}
  containing detector response to scalar gravitational radiation
  over a span of time. (NB: This is unused for the moment. Response to
  scalar gravitational radiation is not yet implemented.)
\end{description}
</lalLaTeX> */

typedef struct
tagLALDetAMResponseSeries
{
  REAL4TimeSeries *pPlus;
  REAL4TimeSeries *pCross;
  REAL4TimeSeries *pScalar;
}
LALDetAMResponseSeries;

/* <lalLaTeX>
\subsubsection*{Structure \texttt{LALGPSandAcc}}
\idx[Type]{LALGPSandAcc}

This structure aggregates GPS time and leap second accuracy requirement
for converting GPS time to sidereal time (implicitly used in
\texttt{LALComputeDetAMResponse()}).

\begin{description}
\item[\texttt{LIGOTimeGPS gps}] The GPS time
\item[\texttt{LALLeapSecAccuracy accuracy}] The accuracy parameter
\end{description}

</lalLaTeX> */
typedef struct
tagLALGPSandAcc
{
  LIGOTimeGPS        gps;      /* GPS time */
  LALLeapSecAccuracy accuracy; /* required accuracy in leap second
                                  handling */
}
LALGPSandAcc;


/* <lalLaTeX>
\subsubsection*{Structure \texttt{LALTimeIntervalAndNSample}}
\idx[Type]{LALTimeIntervalAndNSample}

This structure encapsulates time and sampling information for computing a
\texttt{LALDetAMResponseSeries}. Its fields correspond to some fields of the
\texttt{TimeSeries} structures for easy conversion.

\begin{description}
\item[\texttt{LIGOTimeGPS epoch}] The start time $t_0$ of the time series
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$, in seconds
\item[\texttt{UINT4 nSample}] The total number of samples to be computed
\item[\texttt{LALLeapSecAccuracy accuracy}] The required accuracy for handling leap seconds
\end{description}

</lalLaTeX> */

typedef struct
tagLALTimeIntervalAndNSample
{
  LIGOTimeGPS     epoch;
  REAL8           deltaT;    /* sampling interval */
  UINT4           nSample;   /* number of samples */
  LALLeapSecAccuracy accuracy; /* accuracy for handling leap-seconds */
}
LALTimeIntervalAndNSample;



/*
<lalLaTeX>
\vfill{\footnotesize\input{DetResponseHV}}
</lalLaTeX>
*/

/*
 * Function prototypes
 */


/*
<lalLaTeX>
\newpage\input{DetResponseC}
</lalLaTeX>
*/


void
LALComputeDetAMResponse( LALStatus             *status,
                         LALDetAMResponse      *pResponse,
                         const LALDetAndSource *pDetAndSrc,
                         const LALGPSandAcc    *pGPSandAcc);

void XLALComputeDetAMResponse(
	double *fplus,
	double *fcross,
	REAL4 D[3][3],
	const double ra,
	const double dec,
	const double psi,
	const double gmst
);

/*
 * Gives a time series of the detector's response to plus and cross
 * polarization
 */
void
LALComputeDetAMResponseSeries( LALStatus                      *status,
                               LALDetAMResponseSeries         *pResponseSeries,
                               const LALDetAndSource          *pDetAndSource,
                               const LALTimeIntervalAndNSample *pTimeInfo);

int XLALComputeDetAMResponseSeries(
	REAL4TimeSeries **fplus,
	REAL4TimeSeries **fcross,
	REAL4 D[3][3],
	const double ra,
	const double dec,
	const double psi,
	const LIGOTimeGPS *start,
	const double deltaT,
	const int n
);


#ifdef __cplusplus
}
#endif

#endif /* !defined _DETRESPONSE_H */



