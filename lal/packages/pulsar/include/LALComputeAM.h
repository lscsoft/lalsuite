/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Steve Berukoff
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

/*** <lalVerbatim file="LALComputeAMHV">
Author: Berukoff, S.J.
$Id$
*** </lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{LALComputeAM.h}}
\label{s:LALComputeAM.h}
Computes filter components for amplitude demodulation.

\subsection*{Synposis}
\begin{verbatim}
#include <lal/LALComputeAM.h>
\end{verbatim}

\noindent  In order to compute the optimal statistic for pulsar searches, one must take account of the various modulations that change the emitted, (fairly) simple sinusoid into a non-trivial function of parameters.  The frequency evolution of the signal (spindown effects, Doppler modulation, etc.) have already been accounted for; this routine filters the amplitude modulation effects.

</lalLaTeX> */


#ifndef _LALCOMPUTEAM_H
#define _LALCOMPUTEAM_H

#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include "LALBarycenter.h"

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID (LALCOMPUTEAMH, "$Id: LALComputeAM.h");

  /* Author-defined error codes */
  /* <lalLaTeX>
     \subsection*{Error conditions}
     \vspace{0.1in}
     \input{LALComputeAMHErrorTable}
     </lalLaTeX> */

  /* <lalErrTable file="LALComputeAMHErrorTable"> */
#define LALCOMPUTEAMH_ENOTS        1
#define LALCOMPUTEAMH_EBCERR       2
#define LALCOMPUTEAMH_EESERR       3
#define LALCOMPUTEAMH_EEPH         4
#define LALCOMPUTEAMH_EDAS         5
#define LALCOMPUTEAMH_EFRD         6

#define LALCOMPUTEAMH_MSGENOTS    "Input LIGOTimeGPS Vector is wrong size or NULL"
#define LALCOMPUTEAMH_MSGEBCERR   "Baryinput pointer is invalid"
#define LALCOMPUTEAMH_MSGEESERR   "EarthState structure invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEEPH     "Ephemeris Table invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEDAS     "Detector and source information invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEFRD     "Detector geometry information invalid, or pointer NULL"
  /* </lalErrTable> */

  /* <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
struct AMCoeffs
\end{verbatim}
\index{\texttt{AMCoeffs}}

\noindent This structure contains the output of the routine: $a(t)$, $b(t)$, and the scalar products therein.  That is:

\begin{description}
\item[\texttt{REAL4Vector *a}]  The function $a(t)$
\item[\texttt{REAL4Vector *b}]  The function $b(t)$
\item[\texttt{REAL4 A}]  The scalar product $(a||a)$
\item[\texttt{REAL4 B}]  The scalar product $(b||b)$
\item[\texttt{REAL4 C}]  The scalar product $(a||b)$
\item[\texttt{REAL4 D}]  The quantity $AB-C^{2}$
\end{description}

\begin{verbatim}
struct CmplxAMCoeffs
\end{verbatim}
\index{\texttt{CmplxAMCoeffs}}

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

\begin{verbatim}
struct AMCoeffsParams
\end{verbatim}
\index{\texttt{AMCoeffsParams}}

\noindent This structure contains the parameters for the routine.  They include:

\begin{description}
\item[\texttt{BarycenterInput *baryinput}]  Parameters from $LALBarycenter()$
\item[\texttt{EarthState *earth}]  The state of the earth at time t
\item[\texttt{EphemerisDate *edat}]  Pointer to the ephemerides
\item[\texttt{LALDetAndSource *das}]  Detector and source information
\item[\texttt{LALFrDetector}]  Detector geometry information
\item[\texttt{REAL4 polAngle}]  Polarization angle
\item[\texttt{LALLeapSecAccuracy leapAcc}] Leap sec accuracy
\end{description}
</lalLaTeX> */

typedef struct AMCoeffsTag
{
  REAL4Vector     *a;          /**< the function a(t)         */
  REAL4Vector     *b;          /**< the function b(t)         */
  REAL4           A;           /**< the scalar product (a||a) */
  REAL4           B;           /**< the scalar product (b||b) */
  REAL4           C;           /**< the scalar product (a||b) */
  REAL4           D;           /**< the quantity AB-C^2       */
} AMCoeffs;

typedef struct CmplxAMCoeffsTag
{
  COMPLEX8Vector     *a;          /**< the a coefficient evaluated at the relevant times */
  COMPLEX8Vector     *b;          /**< the b coefficient evaluated at the relevant times  */
} CmplxAMCoeffs;

typedef struct AMCoeffsParamsTag
{
  BarycenterInput      *baryinput;  /**< data from Barycentring routine */
  EarthState           *earth;      /**< from LALBarycenter()           */
  EphemerisData        *edat;       /**< the ephemerides                */
  LALDetAndSource      *das;        /**< det and source information     */
  LALFrDetector        *det;        /**< detector geometry              */
  REAL4                polAngle;    /**< polarization angle             */
  LALLeapSecAccuracy  leapAcc;     /**< accuracy in def of leap sec */
} AMCoeffsParams;


  /* <lalLaTeX>
\newpage\input{LALComputeAMHV}
\newpage\input{LALComputeAMC}
</lalLaTeX> */

void LALComputeAM (LALStatus          *status,
		   AMCoeffs           *coe,
		   LIGOTimeGPS        *ts,
		   AMCoeffsParams     *params);

#ifdef __cplusplus
}
#endif

#endif /* _LALCOMPUTEAM_H */
