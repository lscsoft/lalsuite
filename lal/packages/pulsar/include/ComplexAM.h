/*
*  Copyright (C) 2007 John Whelan
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

/*** <lalVerbatim file="ComplexAMHV">
Author: Whelan, J.T.
$Id$
*** </lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{ComplexAM.h}}
\label{s:ComplexAM.h}
Computes filter components for amplitude demodulation.

\subsection*{Synposis}
\begin{verbatim}
#include <lal/ComplexAM.h>
\end{verbatim}

\noindent  structures and prototypes associated with complex AM coefficients

</lalLaTeX> */


#ifndef _COMPLEXAM_H
#define _COMPLEXAM_H

#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/ComputeFstat.h>
#include "LALBarycenter.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  NRCSID (COMPLEXAMH, "$Id: ComplexAM.h");


/*----- Error-codes -----*/
#define COMPLEXAMC_ENULL 		1
#define COMPLEXAMC_ENONULL 		2
#define COMPLEXAMC_EINPUT   		3
#define COMPLEXAMC_EMEM   		4
#define COMPLEXAMC_EXLAL		5
#define COMPLEXAMC_EIEEE		6

#define COMPLEXAMC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPLEXAMC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPLEXAMC_MSGEINPUT   	"Invalid input"
#define COMPLEXAMC_MSGEMEM   	"Out of memory. Bad."
#define COMPLEXAMC_MSGEXLAL		"XLAL function call failed"
#define COMPLEXAMC_MSGEIEEE		"Floating point failure"

  /* <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
struct AMCoeffs
\end{verbatim}
\index{\texttt{AMCoeffs}}

\noindent This structure contains the AM co\"{e}fficients $a$ and $b$
in the case of a complex detector tensor, and some relevant scalar
products. That is:

\begin{description}
\item[\texttt{COMPLEX8Vector *a}]  The $a$ co\"{e}fficient evaluated at the relevant times
\item[\texttt{COMPLEX8Vector *b}]  The $b$ co\"{e}fficient evaluated at the relevant times
\item[\texttt{REAL4 A}]  The scalar product $(a||a)$
\item[\texttt{REAL4 B}]  The scalar product $(b||b)$
\item[\texttt{REAL4 C}]  The scalar product $(a||b)$
\item[\texttt{REAL4 E}]  The scalar product $(a||ib)$
\item[\texttt{REAL4 D}]  The quantity $AB-C^{2}-E^{2}$
\end{description}

</lalLaTeX> */

typedef struct CmplxAMCoeffsTag
{
  COMPLEX8Vector     *a;          /**< the a coefficient evaluated at the relevant times */
  COMPLEX8Vector     *b;          /**< the b coefficient evaluated at the relevant times  */
  REAL4               A;          /**< the scalar product (a||a) */
  REAL4               B;          /**< the scalar product (b||b) */
  REAL4               C;          /**< the scalar product (a||b) */
  REAL4               E;          /**< the scalar product (a||ib) */
  REAL4               D;          /**< the quantity AB-C^2-E^2    */
} CmplxAMCoeffs;

/** The 'detector tensor' for a GW-detector: symmetric 3x3 matrix, storing only the upper triangle.
 * The coordinate-system is SSB-fixed Cartesian coordinates, in particular EQUATORIAL coords for 
 * Earth-based detectors and ECLIPTIC coords for LISA.
 */
typedef struct 
{
  COMPLEX8 d11;   COMPLEX8 d12;   COMPLEX8 d13;
                  COMPLEX8 d22;   COMPLEX8 d23;
                                  COMPLEX8 d33;
} CmplxDetectorTensor;
  
/* ----- Output types for LALGetDetectorStates() */
/** State-info about position, velocity and LMST of a detector together 
 * with corresponding EarthState.
 */
typedef struct
{
  LIGOTimeGPS tGPS;		/**< GPS timestamps corresponding to this entry */
  REAL8 rDetector[3];		/**< Cartesian coords of detector position in ICRS J2000. Units=sec */
  REAL8 vDetector[3];		/**< Cart. coords. of detector velocity, in dimensionless units (v/c)*/
  CmplxDetectorTensor detT;		/**< Detector-tensor components in SSB-fixed, Cartesian coordinates */
  REAL8 LMST;			/**< local mean sidereal time at the detector-location in radians */
  EarthState earthState;	/**< EarthState information */
} CmplxDetectorState;


/** Timeseries of CmplxDetectorState's, representing the detector-info at different timestamps.
 * In addition to the standard 'vector'-fields we also store the detector-info in here.
 */
typedef struct
{
  UINT4 length;			/**< total number of entries */
  CmplxDetectorState *data;		/**< array of CmplxDetectorState entries */
  LALDetector detector;		/**< detector-info corresponding to this timeseries */
  CoordinateSystem system; 	/**< The coordinate system used for detector's position/velocity and detector-tensor */
} CmplxDetectorStateSeries;

/** Multi-IFO time-series of CmplxDetectorStates */
typedef struct
{
  UINT4 length;			/**< number of detectors */
  CmplxDetectorStateSeries **data;	/**< vector of pointers to CmplxDetectorStateSeries */
  LIGOTimeGPS startTime;	/**< (earliest) startTime of the observation */
  REAL8 Tspan;			/**< total spanned duration of the observation */
} MultiCmplxDetectorStateSeries;

  /* <lalLaTeX>
\newpage\input{ComplexAMHV}
% \newpage\input{ComplexAMC}
</lalLaTeX> */

/*---------- exported prototypes [API] ----------*/
void
LALGetCmplxAMCoeffs(LALStatus *,
		    CmplxAMCoeffs *coeffs,
		    const CmplxDetectorStateSeries *DetectorStates,
		    SkyPosition skypos);

#ifdef __cplusplus
}
#endif

#endif /* _COMPLEXAM_H */
