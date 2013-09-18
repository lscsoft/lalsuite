/*
*  Copyright (C) 2007 John Whelan, Reinhard Prix
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
 * \author Whelan, J.T., Reinhard Prix
 * \date 2007
 * \ingroup pulsarTODO
 * \file
 * \brief structures and prototypes associated with complex AM coefficients
 *
 */

#ifndef _COMPLEXAM_H
#define _COMPLEXAM_H

#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/LALBarycenter.h>
#include <lal/DetectorStates.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- exported defines and macros -------------------- */

/** \name Error codes */
/*@{*/
#define COMPLEXAMC_ENULL 		1
#define COMPLEXAMC_ENONULL 		2
#define COMPLEXAMC_EINPUT   		3
#define COMPLEXAMC_EMEM   		4
#define COMPLEXAMC_EXLAL		5
#define COMPLEXAMC_EIEEE		6
#define COMPLEXAMC_ERAALISA		7

#define COMPLEXAMC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPLEXAMC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPLEXAMC_MSGEINPUT   	"Invalid input"
#define COMPLEXAMC_MSGEMEM   	"Out of memory. Bad."
#define COMPLEXAMC_MSGEXLAL	"XLAL function call failed"
#define COMPLEXAMC_MSGEIEEE	"Floating point failure"
#define COMPLEXAMC_MSGERAALISA	"RAA response only available for LISA"
/*@}*/

/* ---------- exported data types -------------------- */

/**
 * This structure contains the AM coefficients \f$a\f$ and \f$b\f$
 * in the case of a complex detector tensor, and some relevant scalar
 * products. That is:
 */
typedef struct tagCmplxAMCoeffs
{
  COMPLEX8Vector     *a;          /**< the a coefficient evaluated at the relevant times */
  COMPLEX8Vector     *b;          /**< the b coefficient evaluated at the relevant times  */
} CmplxAMCoeffs;

/**
 * Convenience container for precomputed pi f L/c  and skyposition vector
 */
typedef struct tagFreqSkypos_t
{
  REAL4 Freq;		/**< signal frequency */
  REAL8 skyposV[3];	/**< unit vector pointing to skyposition of source */
  SymmTensor3 ePlus;	/**< ePlus polarization tensor (skypos-dependent) */
  SymmTensor3 eCross;	/**< eCross polarization tensor (skypos-dependent) */
} FreqSkypos_t;

/**
 * Struct holding the "antenna-pattern" matrix \f$\mathcal{M}_{\mu\nu} \equiv \left( \mathbf{h}_\mu|\mathbf{h}_\nu\right)\f$,
 * in terms of the multi-detector scalar product. This matrix can be shown to be expressible, in the case of complex AM co\"{e}fficients, as
 * \f{equation}{
 * \mathcal{M}_{\mu\nu} = \mathcal{S}^{-1}\,T_\mathrm{SFT}\,\left( \begin{array}{c c c c}
 * \widehat{A} & \widehat{C} & 0 & -\widehat{E} \\
 * \widehat{C} & \widehat{B} & \widehat{E} & 0 \\
 * 0 & \widehat{E} & \widehat{A} & \widehat{C} \\
 * -\widehat{E} & 0 & \widehat{C} & \widehat{B} \\
 * \end{array}\right)\,,
 * \f}
 * where (here) \f$\mathcal{S} \equiv \frac{1}{N_\mathrm{SFT}}\sum_{X,\alpha} S_{X\alpha}\f$ characterizes the multi-detector noise-floor, and
 * \f{equation}{
 * \widehat{A} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\,,\quad
 * \widehat{B} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha \,,\quad
 * \widehat{C} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * \widehat{E} \equiv \mathrm{Im} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha \,,
 * \f}
 * and the noise-weighted atenna-functions \f$\widehat{a}^X_\alpha = \sqrt{w^X_\alpha}\,a^X_\alpha\f$,
 * \f$\widehat{b}^X_\alpha = \sqrt{w^X_\alpha}\,b^X_\alpha\f$, and noise-weights
 * \f$w^X_\alpha \equiv S^{-1}_{X\alpha}/{\mathcal{S}^{-1}}\f$.
 *
 * \note One reason for storing the un-normalized \a Ad, \a Bd, \a Cd, \a Ed and the normalization-factor \a Sinv_Tsft separately
 * is that the former are of order unity, while \a Sinv_Tsft is very large, and it has numerical advantages for parameter-estimation
 * to use that fact.
 */
typedef struct tagCmplxAntennaPatternMatrix {
  REAL8 Ad; 		/**<  \f$\widehat{A} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{a}^X_\alpha\f$ */
  REAL8 Bd; 		/**<  \f$\widehat{B} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{b}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Cd; 		/**<  \f$\widehat{C} \equiv \mathrm{Re} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Ed; 		/**<  \f$\widehat{E} \equiv \mathrm{Im} \sum_{X,\alpha} \widehat{a}^X_\alpha{}^* \widehat{b}^X_\alpha\f$ */
  REAL8 Dd; 		/**<  determinant \f$\widehat{D} \equiv \widehat{A} \widehat{B} - \widehat{C}^2 -\widehat{E}^2 \f$ */
  REAL8 Sinv_Tsft;	/**< normalization-factor \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (wrt single-sided PSD!) */
} CmplxAntennaPatternMatrix;


/** Multi-IFO container for antenna-pattern coefficients a^X(t), b^X(t) and atenna-pattern matrix M_mu_nu */
typedef struct tagMultiCmplxAMCoeffs {
  UINT4 length;		/**< number of IFOs */
  CmplxAMCoeffs **data;	/**< noise-weighted am-coeffs \f$\widehat{a}_{X\alpha}\f$, and \f$\widehat{b}_{X\alpha}\f$ */
  CmplxAntennaPatternMatrix Mmunu;	/**< antenna-pattern matrix \f$\mathcal{M}_{\mu\nu}\f$ */
} MultiCmplxAMCoeffs;

/*---------- exported prototypes [API] ----------*/
void LALGetCmplxAMCoeffs( LALStatus *, CmplxAMCoeffs *coeffs, const DetectorStateSeries *DetectorStates, const FreqSkypos_t *freq_skypos );
void LALGetMultiCmplxAMCoeffs( LALStatus *, MultiCmplxAMCoeffs **multiAMcoef, const MultiDetectorStateSeries *multiDetStates, PulsarDopplerParams doppler );

int XLALWeightMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

/* destructors */
void XLALDestroyMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef );

#ifdef __cplusplus
}
#endif

#endif /* _COMPLEXAM_H */
