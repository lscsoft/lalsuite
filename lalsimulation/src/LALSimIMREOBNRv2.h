/*
*  Copyright (C) 2010 Craig Robinson 
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
 * \author Craig Robinson
 *
 * \brief File containing most of the structures and prototypes which are
 * used in the generation of the EOBNRv2 waveform. These structures and
 * prototypes are used internally within the waveform generation code, 
 * and shouldn't be needed to generate the waveforms from outside the
 * library. Therefore, people generating EOBNRv2 waveforms should only
 * need to include LALSimIMR.h.
 */

#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>

#ifndef _LALSIMREOBNRv2_H
#define _LALSIMREOBNRv2_H

/**
 * The maximum possible l we have 
 */
#define LALEOB_MAX_MULTIPOLE 8

/**
 *  Structure containing the coefficients for EOBNRv2 A potential function.
 *  The elements in the structure are labelled as follows:
 *  aN, where a is denotes whether the parameter is in the numerator (n)
 *  or denominator (d); and N is the power of r which will multiply this
 *  coefficient. For example, the coefficient of r^5 in the numerator
 *  will be called n5.
 */
typedef struct
tagEOBACoefficients
{
  REAL8 n4;
  REAL8 n5;
  REAL8 d0;
  REAL8 d1;
  REAL8 d2;
  REAL8 d3;
  REAL8 d4;
  REAL8 d5;
}
EOBACoefficients;

/**
 *  Structure containing the coefficients for calculating the factorized
 *  waveform. The coefficients are precomputed in the function
 *  XLALCalcFacWaveformCoefficients()
 */
typedef struct
tagFacWaveformCoeffs
{
  REAL8 delta22vh3;
  REAL8 delta22vh6;
  REAL8 delta22vh8;
  REAL8 delta22vh9;
  REAL8 delta22v5;

  REAL8 rho22v2;
  REAL8 rho22v3;
  REAL8 rho22v4;
  REAL8 rho22v5;
  REAL8 rho22v6;
  REAL8 rho22v6l;
  REAL8 rho22v7;
  REAL8 rho22v8;
  REAL8 rho22v8l;
  REAL8 rho22v10;
  REAL8 rho22v10l;

  REAL8 delta21vh3;
  REAL8 delta21vh6;
  REAL8 delta21vh7;
  REAL8 delta21vh9;
  REAL8 delta21v5;
  REAL8 delta21v7;

  REAL8 rho21v1;
  REAL8 rho21v2;
  REAL8 rho21v3;
  REAL8 rho21v4;
  REAL8 rho21v5;
  REAL8 rho21v6;
  REAL8 rho21v6l;
  REAL8 rho21v7;
  REAL8 rho21v7l;
  REAL8 rho21v8;
  REAL8 rho21v8l;
  REAL8 rho21v10;
  REAL8 rho21v10l;

  REAL8 delta33vh3;
  REAL8 delta33vh6;
  REAL8 delta33vh9;
  REAL8 delta33v5;
  REAL8 delta33v7;

  REAL8 rho33v2;
  REAL8 rho33v3;
  REAL8 rho33v4;
  REAL8 rho33v5;
  REAL8 rho33v6;
  REAL8 rho33v6l;
  REAL8 rho33v7;
  REAL8 rho33v8;
  REAL8 rho33v8l;

  REAL8 delta32vh3;
  REAL8 delta32vh4;
  REAL8 delta32vh6;
  REAL8 delta32vh9;

  REAL8 rho32v;
  REAL8 rho32v2;
  REAL8 rho32v3;
  REAL8 rho32v4;
  REAL8 rho32v5;
  REAL8 rho32v6;
  REAL8 rho32v6l;
  REAL8 rho32v8;
  REAL8 rho32v8l;

  REAL8 delta31vh3;
  REAL8 delta31vh6;
  REAL8 delta31vh7;
  REAL8 delta31vh9;
  REAL8 delta31v5;

  REAL8 rho31v2;
  REAL8 rho31v3;
  REAL8 rho31v4;
  REAL8 rho31v5;
  REAL8 rho31v6;
  REAL8 rho31v6l;
  REAL8 rho31v7;
  REAL8 rho31v8;
  REAL8 rho31v8l;

  REAL8 delta44vh3;
  REAL8 delta44vh6;
  REAL8 delta44v5;

  REAL8 rho44v2;
  REAL8 rho44v3;
  REAL8 rho44v4;
  REAL8 rho44v5;
  REAL8 rho44v6;
  REAL8 rho44v6l;

  REAL8 delta43vh3;
  REAL8 delta43vh4;
  REAL8 delta43vh6;

  REAL8 rho43v;
  REAL8 rho43v2;
  REAL8 rho43v4;
  REAL8 rho43v5;
  REAL8 rho43v6;
  REAL8 rho43v6l;

  REAL8 delta42vh3;
  REAL8 delta42vh6;

  REAL8 rho42v2;
  REAL8 rho42v3;
  REAL8 rho42v4;
  REAL8 rho42v5;
  REAL8 rho42v6;
  REAL8 rho42v6l;

  REAL8 delta41vh3;
  REAL8 delta41vh4;
  REAL8 delta41vh6;

  REAL8 rho41v;
  REAL8 rho41v2;
  REAL8 rho41v4;
  REAL8 rho41v5;
  REAL8 rho41v6;
  REAL8 rho41v6l;

  REAL8 delta55vh3;
  REAL8 delta55v5;
  REAL8 rho55v2;
  REAL8 rho55v3;
  REAL8 rho55v4;
  REAL8 rho55v5;
  REAL8 rho55v6;

  REAL8 delta54vh3;
  REAL8 delta54vh4;
  REAL8 rho54v2;
  REAL8 rho54v3;
  REAL8 rho54v4;

  REAL8 delta53vh3;
  REAL8 rho53v2;
  REAL8 rho53v3;
  REAL8 rho53v4;
  REAL8 rho53v5;

  REAL8 delta52vh3;
  REAL8 delta52vh4;
  REAL8 rho52v2;
  REAL8 rho52v3;
  REAL8 rho52v4;

  REAL8 delta51vh3;
  REAL8 rho51v2;
  REAL8 rho51v3;
  REAL8 rho51v4;
  REAL8 rho51v5;

  REAL8 delta66vh3;
  REAL8 rho66v2;
  REAL8 rho66v3;
  REAL8 rho66v4;

  REAL8 delta65vh3;
  REAL8 rho65v2;
  REAL8 rho65v3;

  REAL8 delta64vh3;
  REAL8 rho64v2;
  REAL8 rho64v3;
  REAL8 rho64v4;

  REAL8 delta63vh3;
  REAL8 rho63v2;
  REAL8 rho63v3;

  REAL8 delta62vh3;
  REAL8 rho62v2;
  REAL8 rho62v3;
  REAL8 rho62v4;

  REAL8 delta61vh3;
  REAL8 rho61v2;
  REAL8 rho61v3;

  REAL8 delta77vh3;
  REAL8 rho77v2;
  REAL8 rho77v3;

  REAL8 rho76v2;

  REAL8 delta75vh3;
  REAL8 rho75v2;
  REAL8 rho75v3;

  REAL8 rho74v2;

  REAL8 delta73vh3;
  REAL8 rho73v2;
  REAL8 rho73v3;

  REAL8 rho72v2;

  REAL8 delta71vh3;
  REAL8 rho71v2;
  REAL8 rho71v3;

  REAL8 rho88v2;
  REAL8 rho87v2;
  REAL8 rho86v2;
  REAL8 rho85v2;
  REAL8 rho84v2;
  REAL8 rho83v2;
  REAL8 rho82v2;
  REAL8 rho81v2;
}
FacWaveformCoeffs;

/**
 * Structure containing all the terms of the Newtonian multipole which
 * are constant over the course of the evolution, and can therefore be
 * pre-computed. They are stored in a two-dimensional array, which is
 * indexed as values[l][m]. Since m has to be <= l, this structure
 * is larger than it needs to be; but it makes the coding a bit neater...
 */
typedef
struct tagNewtonMultipolePrefixes
{
  COMPLEX16 values[LALEOB_MAX_MULTIPOLE+1][LALEOB_MAX_MULTIPOLE+1];
}
NewtonMultipolePrefixes;

/**
 * The coefficients which are used in calculating the non-quasicircular
 * correction to the EOBNRv2 model. The precise definitions of these 
 * coefficients and their use can be found in DCC document T1100433.
 */
typedef struct
tagEOBNonQCCoeffs
{
  REAL8 a1;
  REAL8 a2;
  REAL8 a3;
  REAL8 a4;
  REAL8 b1;
  REAL8 b2;
} EOBNonQCCoeffs;

/**
 * Structure containing all the parameters needed for the EOB waveform.
 * It contains eta, the pre-computed parameters for the A potential function,
 * and the pre-computed parameters for the factorized waveform
 */

typedef
struct tagEOBParams
{
  REAL8 eta;
  REAL8 omega;
  REAL8 m1;
  REAL8 m2;
  EOBACoefficients        *aCoeffs;
  FacWaveformCoeffs       *hCoeffs;
  EOBNonQCCoeffs          *nqcCoeffs;
  NewtonMultipolePrefixes *prefixes;
}
EOBParams;

/**
 * This function calculates the factorized flux in the EOB dynamics for
 * the EOBNR (and potentially subsequent) models. The flux function
 * is found in Phys.Rev.D79:064004,2009.
 */
static REAL8 XLALSimIMREOBFactorizedFlux(
                      REAL8Vector  *values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8  omega,   /**<< Angular frequency omega */
                      EOBParams    *ak,     /**<< Structure containing pre-computed parameters */
                      const INT4   lMax     /**<< Maximum l to include when calculating flux (between 2 and 8) */
                     );

/**
 * Function which computes the various coefficients in the Newtonian
 * multipole. The definition of this can be found in Pan et al, 
 * arXiv:1106.1021v1 [gr-qc]. Note that, although this function gets passed
 * masses measured in Solar masses (which is fine since it is a static function),
 * the units of mass won't matter so long as they are consistent. This is because
 * all uses of the mass within the function are normalized by the total mass.
 */
static int XLALSimIMREOBComputeNewtonMultipolePrefixes(
                NewtonMultipolePrefixes *prefix, /**<< Structure containing the coefficients (populated in function) */
                const REAL8             m1,      /**<< Mass of first component */
                const REAL8             m2       /**<< Nass of second component */
                );

/**
 * This function calculates the Newtonian multipole part of the
 * factorized waveform. This is defined in Pan et al, arXiv:1106.1021v1 [gr-qc].
 */
static int
XLALSimIMREOBCalculateNewtonianMultipole(
                            COMPLEX16 *multipole, /**<< Newtonian multipole (returned) */
                            REAL8 x,              /**<< Dimensionless parameter \f$\equiv v^2\f$ */
                            REAL8 r,              /**<< Orbital separation (units of total mass M */
                            REAL8 phi,            /**<< Orbital phase (in radians) */
                            UINT4  l,             /**<< Mode l */
                            INT4  m,              /**<< Mode m */
                            EOBParams *params     /**<< Pre-computed coefficients, parameters, etc. */
                            );

/**
 * Function which calculates the various coefficients used in the generation
 * of the factorized waveform. These coefficients depend only on the symmetric
 * mass ratio eta. It should be noted that this function calculates the 
 * coefficients used in calculating the flux. For generating the waveforms
 * themselves, the coefficients have additional terms added which are calculated
 * using XLALModifyFacWaveformCoefficients(). THe non-spinning parts of these
 * coefficients can be found in Pan et al, arXiv:1106.1021v1 [gr-qc].
 */
static int XLALSimIMREOBCalcFacWaveformCoefficients(
          FacWaveformCoeffs * const coeffs, /**<< Structure containing coefficients (populated in function) */
          const REAL8               eta     /**<< Symmetric mass ratio */
          );

/**
 * Computes the factorized waveform according to the prescription
 * given in Pan et al, arXiv:1106.1021v1 [gr-qc], for a given
 * mode l,m, for the given values of the dynamics at that point.
 * The function returns XLAL_SUCCESS if everything works out properly,
 * otherwise XLAL_FAILURE will be returned.
 */
static int XLALSimIMREOBGetFactorizedWaveform( 
                                COMPLEX16   * restrict hlm,    /**<< The value of hlm (populated by the function) */
                                REAL8Vector * restrict values, /**<< Vector containing dynamics r, phi, pr, pphi for a given point */
                                const REAL8 v,                 /**<< Velocity (in geometric units) */
                                const INT4  l,                 /**<< Mode l */
                                const INT4  m,                 /**<< Mode m */
                                EOBParams   * restrict params  /**<< Structure containing pre-computed coefficients, etc. */
                                );


/**
 * This function generates the quasinormal mode frequencies for a black
 * hole ringdown. At present, this function works for the 22, 21, 33, 44
 * and 55 modes, and includes 8 overtones. The final frequencies are
 * computed by interpolating the data found on the webpage of 
 * Emanuele Berti, http://www.phy.olemiss.edu/~berti/qnms.html
 */
static INT4 XLALSimIMREOBGenerateQNMFreqV2(
        COMPLEX16Vector          *modefreqs, /**<<The complex frequencies of the overtones (scaled by total mass) */
        const REAL8              mass1,      /**<<The mass of the first component of the system (in Solar masses) */
        const REAL8              mass2,      /**<<The mass of the second component of the system (in Solar masses) */
        UINT4                    l,          /**<<The l value of the mode in question */
        UINT4                    m,          /**<<The m value of the mode in question */
        UINT4                   nmodes       /**<<The number of overtones that should be included (max 8) */
        );

/**
 * The main workhorse function for performing the ringdown attachment for EOB
 * models EOBNRv2 and later. This is the function which gets called by the 
 * code generating the full IMR waveform once generation of the inspiral part
 * has been completed.
 * The ringdown is attached using the hybrid comb matching detailed in 
 * Buonanno et al, arXiv:1106.1021v1 [gr-qc]. Further details of the
 * implementation of the found in the DCC document T1100433.
 */
static int XLALSimIMREOBHybridAttachRingdown(
      REAL8Vector       *signal1,     /**<<Real part of inspiral waveform to which we attach the ringdown */
      REAL8Vector       *signal2,     /**<<Imaginary part of inspiral waveform to which we attach the ringdown */
      const INT4        l,            /**<< Current mode l */
      const INT4        m,            /**<< Current mode m */
      const REAL8       dt,           /**<< Sample time step (in seconds) */
      const REAL8       mass1,        /**<< First component mass (in Solar masses) */
      const REAL8       mass2,        /**<< Second component mass (in Solar masses) */
      REAL8Vector       *timeVec,     /**<< Vector containing the time values */
      REAL8Vector       *matchrange   /**<< Time values chosen as points for performing comb matching */
      );

/**
 * Compute the time offset which should be used in computing the
 * non-quasicircular correction and performing the ringdown attachment.
 * These numbers were tuned to numerical relativity simulations, and
 * are taken from Pan et al, arXiv:1106.1021v1 [gr-qc] 
 */
 static REAL8 XLALSimIMREOBGetNRPeakDeltaT(
                         INT4 l,    /**<< Mode l */
                         INT4 m,    /**<< Mode m */
                         REAL8 eta  /**<< Symmetric mass ratio */
                         );

/**
 * This function computes the coefficients a1, a2, etc. used in the
 * non-quasicircular correction. The details of the calculation of these
 * coefficients are found in the DCC document T1100433. */
static int XLALSimIMREOBCalculateNQCCoefficients(
                 EOBNonQCCoeffs * restrict coeffs,    /**<< NQC coefficients (populated by function) */
                 REAL8Vector    * restrict amplitude, /**<< Amplitude of waveform as function of time */
                 REAL8Vector    * restrict phase,     /**<< Phase of waveform (in radians) as function of time */
                 REAL8Vector    * restrict q1,        /**<< Function of dynamics (see DCC document for details) */
                 REAL8Vector    * restrict q2,        /**<< Function of dynamics (see DCC document for details) */
                 REAL8Vector    * restrict q3,        /**<< Function of dynamics (see DCC document for details) */
                 REAL8Vector    * restrict p1,        /**<< Function of dynamics (see DCC document for details) */
                 REAL8Vector    * restrict p2,        /**<< Function of dynamics (see DCC document for details) */
                 INT4                      l,         /**<< Mode l */
                 INT4                      m,         /**<< Mode m */
                 REAL8                     timePeak,  /**<< Time for which we reach the peak frequency */
                 REAL8                     deltaT,    /**<< Sampling interval */
                 REAL8                     eta        /**<< Symmetric mass ratio */
                 );

/**
 * For the 2,2 mode, there are fits available for the NQC coefficients.
 * This function provides the values of these coefficients, so the 
 * correction can be used in the dynamics prior to finding the more
 * accurate NQC values later on.
 */
static int XLALSimIMREOBGetCalibratedNQCCoeffs( 
                                EOBNonQCCoeffs *coeffs, /**<< Structure for NQC coefficients (populated in function) */
                                INT4            l,      /**<< Mode l */
                                INT4            m,      /**<< Mode m */
                                REAL8           eta     /**<< Symmetric mass ratio */
                                );

/**
 * This function calculates the non-quasicircular correction to apply to 
 * the waveform. The form of this correction can be found in  Pan et al, 
 * arXiv:1106.1021v1 [gr-qc], and also in the DCC document T1100433. Note
 * that when calling this function, the NQC coefficients should already 
 * have been pre-computed.
 */
static int  XLALSimIMREOBNonQCCorrection(
                      COMPLEX16      * restrict nqc,    /**<< The NQC correction (populated in function) */
                      REAL8Vector    * restrict values, /**<< Dynamics r, phi, pr, pphi */
                      const REAL8               omega,  /**<< Angular frequency */
                      EOBNonQCCoeffs * restrict coeffs  /**<< NQC coefficients */
                     );

/**
 * Function to calculate the EOB effective Hamiltonian for the
 * given values of the dynamical variables. The coefficients in the
 * A potential function should already have been computed.
 * Note that the pr used here is the tortoise co-ordinate.
 */
static
REAL8 XLALEffectiveHamiltonian( const REAL8 eta,          /**<< Symmetric mass ratio */
                                const REAL8 r,            /**<< Orbital separation */
                                const REAL8 pr,           /**<< Tortoise co-ordinate */
                                const REAL8 pp,           /**<< Momentum pphi */
                                EOBACoefficients *aCoeffs /**<< Pre-computed coefficients in A function */
                              );

/**
 * This function calculates the EOB A function which using the pre-computed
 * coefficients which should already have been calculated.
 */
static
REAL8 XLALCalculateEOBA( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                         EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                       );

/**
 * Calculated the derivative of the EOB A function with respect to 
 * r, using the pre-computed A coefficients
 */
static
REAL8 XLALCalculateEOBdAdr( const REAL8 r,                     /**<< Orbital separation (in units of total mass M) */
                            EOBACoefficients * restrict coeffs /**<< Pre-computed coefficients for the A function */
                          );
#endif /* _LALSIMEOBNRv2_H */
