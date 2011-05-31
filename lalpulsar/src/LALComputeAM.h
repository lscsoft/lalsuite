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

/**
 * \author S.J. Berukoff, Reinhard Prix, John Whelan
 * \date 2007
 * \ingroup pulsarAntenna
 * \file
 * \brief Header-file for computing antenna-pattern components for amplitude demodulation.
 *
 * <tt>\#include <lal/LALComputeAM.h></tt>
 *
 * In order to compute the optimal statistic for pulsar searches, one must take account of the
 * various modulations that change the emitted, (fairly) simple sinusoid into a non-trivial function
 * of parameters.  The frequency evolution of the signal (spindown effects, Doppler modulation, etc.)
 * have already been accounted for; this routine filters the amplitude modulation effects.
 */

#ifndef _LALCOMPUTEAM_H
#define _LALCOMPUTEAM_H

#ifdef __cplusplus
extern "C" {
#endif


/*---------- exported INCLUDES ----------*/
#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/DetectorStates.h>

/* ---------- exported defines and macros -------------------- */
NRCSID (LALCOMPUTEAMH, "$Id: LALComputeAM.h");

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

/** This structure contains the output of the routine: a(t), b(t),
 * and the scalar products therein.  That is:
 */
typedef struct AMCoeffsTag
{
  REAL4Vector     *a;          /**< the function a(t)         */
  REAL4Vector     *b;          /**< the function b(t)         */
  REAL4           A;           /**< the scalar product (a||a) */
  REAL4           B;           /**< the scalar product (b||b) */
  REAL4           C;           /**< the scalar product (a||b) */
  REAL4           D;           /**< the quantity AB-C^2       */
} AMCoeffs;

/** This structure contains the parameters for the routine.  They include:
 */
typedef struct AMCoeffsParamsTag
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
 * \mathcal{M}_{\mu\nu} = \mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c} A_d & C_d & 0 & 0 \\ C_d & B_d & 0 & 0 \\ 0 & 0 & A_d & C_d \\ 0 & 0 & C_d & B_d \\ \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the (single-sided!)
 * multi-detector noise-floor, and
 * \f{equation}
 * A_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{a}^X_\alpha\,,\quad
 * B_d \equiv \sum_{X,\alpha} \widehat{b}^X_\alpha \widehat{b}^X_\alpha \,,\quad
 * C_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{b}^X_\alpha \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}^X_\alpha = \sqrt{w^X_\alpha}\,a^X_\alpha\f$,
 * \f$\widehat{b}^X_\alpha = \sqrt{w^X_\alpha}\,b^X_\alpha\f$, and noise-weights
 * \f$w^X_\alpha \equiv {S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 *
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd and the normalization-factor \a Sinv_Tsft separately
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct {
  REAL8 Ad; 		/**<  \f$A_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{a}^X_\alpha\f$ */
  REAL8 Bd; 		/**<  \f$B_d \equiv \sum_{X,\alpha} \widehat{b}^X_\alpha \widehat{b}^X_\alpha\f$ */
  REAL8 Cd; 		/**<  \f$C_d \equiv \sum_{X,\alpha} \widehat{a}^X_\alpha \widehat{b}^X_\alpha\f$ */
  REAL8 Dd; 		/**<  determinant \f$D_d \equiv A_d B_d - C_d^2 \f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (wrt single-sided PSD!) */
} AntennaPatternMatrix;

/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) and atenna-pattern matrix M_mu_nu */
typedef struct {
  UINT4 length;		/**< number of IFOs */
  AMCoeffs **data;	/**< noise-weighted am-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
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
int XLALWeighMultiAMCoeffs (  MultiAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

AMCoeffs *XLALComputeAMCoeffs ( const DetectorStateSeries *DetectorStates, SkyPosition skypos );
MultiAMCoeffs *XLALComputeMultiAMCoeffs ( const MultiDetectorStateSeries *multiDetStates, const MultiNoiseWeights *multiWeights, SkyPosition skypos );

/* creators and destructors for AM-vectors */
AMCoeffs *XLALCreateAMCoeffs ( UINT4 numSteps );

void XLALDestroyMultiAMCoeffs ( MultiAMCoeffs *multiAMcoef );
void XLALDestroyAMCoeffs ( AMCoeffs *amcoef );


#ifdef __cplusplus
}
#endif

#endif /* _LALCOMPUTEAM_H */
