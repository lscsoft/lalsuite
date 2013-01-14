/*
*  Copyright (C) 2007 John Whelan, Reinhard Prix
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
#ifndef _LALCOMPUTEAM_H
#define _LALCOMPUTEAM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALComputeAM_h Header LALComputeAM.h
 * \ingroup pkg_pulsarCoh
 * \author S.J. Berukoff, Reinhard Prix, John Whelan
 * \date 2007

 * \brief Header-file for computing antenna-pattern components for amplitude demodulation.
 *
 * <tt>\#include <lal/LALComputeAM.h></tt>
 *
 * In order to compute the optimal statistic for pulsar searches, one must take account of the
 * various modulations that change the emitted, (fairly) simple sinusoid into a non-trivial function
 * of parameters.  The frequency evolution of the signal (spindown effects, Doppler modulation, etc.)
 * have already been accounted for; this routine filters the amplitude modulation effects.
 */
/*@{*/

/*---------- exported INCLUDES ----------*/
#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/DetectorStates.h>

/* ---------- exported defines and macros -------------------- */

/** \name Error codes */
/*@{*/
#define LALCOMPUTEAMH_ENOTS        1
#define LALCOMPUTEAMH_EBCERR       2
#define LALCOMPUTEAMH_EESERR       3
#define LALCOMPUTEAMH_EEPH         4
#define LALCOMPUTEAMH_EDAS         5
#define LALCOMPUTEAMH_EFRD         6
#define LALCOMPUTEAMH_ENULL 	   7
#define LALCOMPUTEAMH_EINPUT   	   8
#define LALCOMPUTEAMH_ENONULL 	   9
#define LALCOMPUTEAMH_EMEM   	  10

#define LALCOMPUTEAMH_MSGENOTS    "Input LIGOTimeGPS Vector is wrong size or NULL"
#define LALCOMPUTEAMH_MSGEBCERR   "Baryinput pointer is invalid"
#define LALCOMPUTEAMH_MSGEESERR   "EarthState structure invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEEPH     "Ephemeris Table invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEDAS     "Detector and source information invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGEFRD     "Detector geometry information invalid, or pointer NULL"
#define LALCOMPUTEAMH_MSGENULL 	  "Arguments contained an unexpected null pointer"
#define LALCOMPUTEAMH_MSGEINPUT   "Invalid input"
#define LALCOMPUTEAMH_MSGENONULL  "Output pointer is non-NULL"
#define LALCOMPUTEAMH_MSGEMEM     "Out of memory. Bad."
/*@}*/

/* ---------- exported data types -------------------- */

/** This structure contains the per-SFT (weighted) antenna-pattern functions
 * \f$\widehat{a}_{X\alpha}, \widehat{b}_{X\alpha}\f$,
 * with \f$\alpha\f$ the SFT-index, and \f$X\f$ the IFO index. The per-IFO summed antenna-pattern coefficients
 * for detector X are
 * \f$\widehat{A}_X,\widehat{B}_X,\widehat{C}_X\f$ and their determinant \f$\widehat{D}_X=\widehat{A}_X \widehat{B}_X - \widehat{C}_X^2\f$.
 *
 * \p See Sec.4.1 in CFSv2 notes (https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=T0900149&version=4)
 */
typedef struct tagAMCoeffs
{
  REAL4Vector     *a;          /**< (weighted) per-SFT \f$X\alpha\f$ antenna-pattern function \f$\widehat{a}_{X\alpha}\f$ */
  REAL4Vector     *b;          /**< (weighted) per-SFT \f$X\alpha\f$ antenna-pattern function \f$\widehat{b}_{X\alpha}\f$ */
  REAL4           A;           /**< summed antenna-pattern matrix coefficient: \f$\widehat{A}_X = \sum_{\alpha} \widehat{a}^2_{X\alpha}\f$ */
  REAL4           B;           /**< summed antenna-pattern matrix coefficient: \f$\widehat{B}_X = \sum_{\alpha} \widehat{b}^2_{X\alpha}\f$ */
  REAL4           C;           /**< summed antenna-pattern matrix coefficient: \f$\widehat{C}_X = \sum_{\alpha} \widehat{a}_{X\alpha}\,\widehat{b}_{X\alpha}\f$ */
  REAL4           D;           /**< determinant \f$\widehat{D}_X = \widehat{A}_X \widehat{B}_X - \widehat{C}_X^2\f$  */
} AMCoeffs;

/** This structure contains the parameters for the routine.  They include:
 */
typedef struct tagAMCoeffsParams
{
  BarycenterInput      *baryinput;  /**< data from Barycentring routine */
  EarthState           *earth;      /**< from LALBarycenter()           */
  EphemerisData        *edat;       /**< the ephemerides                */
  LALDetAndSource      *das;        /**< det and source information     */
  LALFrDetector        *det;        /**< detector geometry              */
  REAL4                polAngle;    /**< polarization angle             */
} AMCoeffsParams;


/** Struct holding the "antenna-pattern" matrix \f$\mathcal{M}_{\mu\nu} \equiv \left( \mathbf{h}_\mu|\mathbf{h}_\nu\right)\f$,
 * in terms of the multi-detector scalar product. This matrix can be shown to be expressible as
 * \f{equation}
 * \mathcal{M}_{\mu\nu} = \mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c}
 * \widehat{A} & \widehat{C} & 0 & 0 \\
 * \widehat{C} & \widehat{B} & 0 & 0 \\
 * 0 & 0 & \widehat{A} & \widehat{C} \\
 * 0 & 0 & \widehat{C} & \widehat{B} \\
 * \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the (single-sided!)
 * multi-detector noise-floor, and
 * \f{equation}
 * \widehat{A} \equiv \sum_{X,\alpha} \widehat{a}^2_{X\alpha} \,,\quad
 * \widehat{B} \equiv \sum_{X,\alpha} \widehat{b}^2_{X\alpha} \,,\quad
 * \widehat{C} \equiv \sum_{X,\alpha} \widehat{a}_{X\alpha} \widehat{b}_{X\alpha} \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}_{X\alpha} = \sqrt{w_{X\alpha}}\,a_{X\alpha}\f$,
 * \f$\widehat{b}_{X\alpha} = \sqrt{w_{X\alpha}}\,b_{X\alpha}\f$, and noise-weights
 * \f$w_{X\alpha} \equiv S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 *
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd and the normalization-factor \a Sinv_Tsft separately
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct tagAntennaPatternMatrix {
  REAL8 Ad; 		/**<  \f$\widehat{A} \equiv \sum_{X,\alpha} \widehat{a}^2_{X\alpha} \f$ */
  REAL8 Bd; 		/**<  \f$\widehat{B} \equiv \sum_{X,\alpha} \widehat{b}^2_{X\alpha} \f$ */
  REAL8 Cd; 		/**<  \f$\widehat{C} \equiv \sum_{X,\alpha} \widehat{a}_{X\alpha} \widehat{b}_{X\alpha}\f$ */
  REAL8 Dd; 		/**<  determinant \f$\widehat{D} \equiv \widehat{A} \widehat{B} - \widehat{C}^2 \f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (wrt single-sided PSD!) */
} AntennaPatternMatrix;

/** Multi-IFO container for antenna-pattern coefficients \f$a_{X\alpha}, b_{X\alpha}\f$ and atenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
typedef struct tagMultiAMCoeffs {
  UINT4 length;			/**< number of IFOs */
  AMCoeffs **data;		/**< noise-weighted AM-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
  AntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
} MultiAMCoeffs;


/*---------- exported Global variables ----------*/
/* empty init-structs for the types defined in here */
extern const AntennaPatternMatrix empty_AntennaPatternMatrix;
extern const MultiAMCoeffs empty_MultiAMCoeffs;


/*---------- exported prototypes [API] ----------*/

void LALComputeAM (LALStatus *, AMCoeffs *coe, LIGOTimeGPS *ts, AMCoeffsParams *params);

void LALGetAMCoeffs(LALStatus *, AMCoeffs *coeffs, const DetectorStateSeries *DetectorStates, SkyPosition skypos);

void LALNewGetAMCoeffs(LALStatus *, AMCoeffs *coeffs, const DetectorStateSeries *DetectorStates, SkyPosition skypos);
int XLALComputeAntennaPatternCoeffs ( REAL8 *ai, REAL8 *bi, const SkyPosition *skypos, const LIGOTimeGPS *tGPS, const LALDetector *site, const EphemerisData *edat );
void LALGetMultiAMCoeffs (LALStatus *, MultiAMCoeffs **multiAMcoef, const MultiDetectorStateSeries *multiDetStates, SkyPosition pos );
int XLALWeightMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

AMCoeffs *XLALComputeAMCoeffs ( const DetectorStateSeries *DetectorStates, SkyPosition skypos );
MultiAMCoeffs *XLALComputeMultiAMCoeffs ( const MultiDetectorStateSeries *multiDetStates, const MultiNoiseWeights *multiWeights, SkyPosition skypos );

/* creators and destructors for AM-vectors */
AMCoeffs *XLALCreateAMCoeffs ( UINT4 numSteps );

void XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef );
void XLALDestroyAMCoeffs ( AMCoeffs *amcoef );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _LALCOMPUTEAM_H */
