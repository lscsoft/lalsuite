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

/** Struct holding the "antenna-pattern" matrix \f$\mathcal{M}_{\mu\nu} \equiv \left( \mathbf{h}_\mu|\mathbf{h}_\nu\right)\f$,
 * in terms of the multi-detector scalar product. This matrix can be shown to be expressible, in the case of complex AM co\"{e}fficients, as
 * \f{equation}
 * \mathcal{M}_{\mu\nu} = \frac{1}{2}\,\mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c} A_d & C_d & 0 & -E_d \\ C_d & B_d & E_d & 0 \\ 0 & E_d & A_d & C_d \\ -E_d & 0 & C_d & B_d \\ \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the multi-detector noise-floor, and
 * \f{equation}
 * A_d \equiv \sum_{X,\alpha} \mathrm{Re} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\,,\quad
 * B_d \equiv \sum_{X,\alpha} \mathrm{Re} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha \,,\quad
 * C_d \equiv \sum_{X,\alpha} \mathrm{Re} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * E_d \equiv \sum_{X,\alpha} \mathrm{Im} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}^X_\alpha = \sqrt{w^X_\alpha}\,a^X_\alpha\f$, 
 * \f$\widehat{b}^X_\alpha = \sqrt{w^X_\alpha}\,b^X_\alpha\f$, and noise-weights 
 * \f$w^X_\alpha \equiv {S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 * 
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd, \a Ed and the normalization-factor \a Sinv_Tsft separately 
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct {
  REAL8 Ad; 		/**<  \f$A_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\f$ */
  REAL8 Bd; 		/**<  \f$B_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Cd; 		/**<  \f$C_d \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Ed; 		/**<  \f$E_d \equiv \mathrm{Im} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ */
} CmplxAntennaPatternMatrix;

/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) and atenna-pattern matrix M_mu_nu */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  CmplxAMCoeffs **data;	/**< noise-weighted am-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
  CmplxAntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
} MultiCmplxAMCoeffs;

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

void
LALGetMultiCmplxAMCoeffs (LALStatus *, 
		     MultiCmplxAMCoeffs **multiAMcoef,
		     const MultiCmplxDetectorStateSeries *multiDetStates,
		     SkyPosition pos );

int
XLALWeighMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

/* destructors */
void XLALDestroyMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef );

#ifdef __cplusplus
}
#endif

#endif /* _COMPLEXAM_H */
