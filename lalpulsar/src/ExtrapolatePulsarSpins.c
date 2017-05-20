/*
 * Copyright (C) 2014 Karl Wette
 * Copyright (C) 2005, 2006 Reinhard Prix
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

/* ========== REVIEW information ==================================================
 * This file has been fully reviewed on 13 May, 2009
 * Reviewers: Peter Shawhan, Teviet Creighton, Francesco Salemi, M.A. Papa
 * Minutes of this review:
 * https://www.lsc-group.phys.uwm.edu/twiki/pub/CW/HierarchicalSearchReview/meeting_20090513.txt
 *
 * If you want to modify any existing functions in this file, please submit
 * your patch for review to pulgroup@gravity.phys.uwm.edu
 * ================================================================================
 */

/*********************************************************************************/
/**
 * \file
 * \ingroup ExtrapolatePulsarSpins
 * \author Reinhard Prix
 *
 * \brief Defines functions to extrapolate the pulsar spin-paramters
 * \f$\{f, \dot{f},\ddot{f},...\}\f$ from one SSB epoch to another.
 *
 */

/*---------- INCLUDES ---------- */
#include <math.h>

#include <lal/LALError.h>
#include <lal/Date.h>

#include "ExtrapolatePulsarSpins.h"

/*---------- local DEFINES ----------*/
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

#define TRUE (1==1)
#define FALSE (1==0)

/*---------- main functions ---------- */

/**
 * Initialise a #PulsarSpinRange struct from two #PulsarSpins structs
 */
int XLALInitPulsarSpinRangeFromSpins ( PulsarSpinRange *range,		/**< [out] output spin range */
                                       const LIGOTimeGPS *refTime,	/**< [in] reference time */
                                       const PulsarSpins fkdot1,	/**< [in] input spins */
                                       const PulsarSpins fkdot2		/**< [in] input spins */
  )
{
  XLAL_CHECK( range != NULL, XLAL_EFAULT );
  XLAL_CHECK( refTime != NULL, XLAL_EFAULT );
  XLAL_INIT_MEM( *range );
  range->refTime = *refTime;
  for ( size_t k = 0; k < PULSAR_MAX_SPINS; ++k ) {
    range->fkdot[k] = MYMIN( fkdot1[k], fkdot2[k] );
    range->fkdotBand[k] = fabs( fkdot1[k] - fkdot2[k] );
  }
  return XLAL_SUCCESS;
}

/**
 * General pulsar-spin extrapolation function: given a "spin-range" (ie spins + spin-bands) \c range0
 * at time \f$\tau_0\f$, propagate the whole spin-range to time \f$\tau_1\f$.
 *
 * \note \c *range1 is allowed to point to the same spin-range as \c *range0: the input will be overwritten
 * with the output.
 *
 * \note The output-range is in the 'canonical' order of \f$[ f^{(k)}, f^{(k)} + \Delta f^{(k)}]\f$,
 * where \f$\Delta f^{(k)} \ge 0\f$.
 *
 */
int
XLALExtrapolatePulsarSpinRange ( PulsarSpinRange *range1,		/**< [out] output spin range */
                                 const PulsarSpinRange *range0,		/**< [in] input spin range */
                                 const REAL8 dtau 			/**< [in] time difference \f$\tau_1 - \tau_0\f$ to extrapolate \c range0 to */
                                 )
{
  UINT4 k, l;
  PulsarSpinRange inRange;

  XLAL_CHECK( range1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( range0 != NULL, XLAL_EFAULT );

  /* ----- make a copy of input range, because we allow input == output range, so
   * the input can get overwritten */
  memmove(&inRange, range0, sizeof(inRange));

  for ( l = 0; l < PULSAR_MAX_SPINS; l ++ )
    {
      REAL8 flmin = 0, flmax = 0;
      REAL8 kfact = 1, dtau_powk = 1;	/* values for k=0 */

      for ( k = 0; k < PULSAR_MAX_SPINS - l; k ++ )
	{
	  REAL8 fkltauk0 = inRange.fkdot[k+l] * dtau_powk;
	  REAL8 fkltauk1 = fkltauk0 + inRange.fkdotBand[k+l] * dtau_powk;

	  REAL8 fkltauk_min = MYMIN ( fkltauk0, fkltauk1 );
	  REAL8 fkltauk_max = MYMAX ( fkltauk0, fkltauk1 );

	  flmin += fkltauk_min / kfact;
	  flmax += fkltauk_max / kfact;

	  kfact *= (k+1);
	  dtau_powk *= dtau;

	} /* for k < PULSAR_MAX_SPINS */

      range1->fkdot[l]     = flmin;
      range1->fkdotBand[l] = flmax - flmin;

    } /* for l < PULSAR_MAX_SPINS */

  /* set proper epoch for output */
  range1->refTime = range0->refTime;
  XLALGPSAdd(&range1->refTime, dtau);

  return XLAL_SUCCESS;

} /* XLALExtrapolatePulsarSpinRange() */


/**
 * Extrapolate the Pulsar spin-parameters \f$\{f, \dot{f},\ddot{f},...\}\f$
 * (\c fkdot0) from the initial reference-epoch \f$\tau_0\f$
 * to the new reference-epoch \f$\tau_1\f$.
 *
 * This is equivalent to XLALExtrapolatePulsarSpins(), but uses the fixed-size array-type
 * ::PulsarSpins instead, which is easier to handle and avoids any dynamic-memory hassles.
 *
 * \note This can be called with <tt>fkdot1 == fkdot0</tt>, in which case the input will
 * be correctly replaced by the output.
 */
int
XLALExtrapolatePulsarSpins ( PulsarSpins fkdot1,		/**< [out] output spin-parameter array */
			     const PulsarSpins fkdot0,		/**< [in] input spin-parameter array */
			     REAL8 dtau 			/**< [in] time difference \f$\tau_1 - \tau_0\f$ to extrapolate \c fkdot0 to */
			     )
{
  UINT4 numSpins = sizeof(PulsarSpins) / sizeof(fkdot0[0]); 	/* fixed size array */
  UINT4 k, l;
  REAL8 kfact, dtauk;
  PulsarSpins inSpins;

  /* if dtau is zero, just copy fkdot0 to fkdot1 */
  if ( dtau == 0.0 ) {
    memmove ( fkdot1, fkdot0, sizeof(PulsarSpins) );
    return XLAL_SUCCESS;
  }

  /* keep a local copy of input to allow the input- and output- pointers to be identical */
  memcpy ( inSpins, fkdot0, sizeof(PulsarSpins) );

  for ( l = 0; l < numSpins; l ++ )
    fkdot1[l] = 0;

  kfact = 1;
  dtauk = 1;	/* values of k! and (dTau)^k at k=0 */
  for ( k = 0; k < numSpins; k ++ )
    {
      REAL8 kcoef  = dtauk / kfact;
      for ( l=0; l < numSpins - k ; l ++ )
	fkdot1[l] += inSpins[ k + l ] * kcoef;

      kfact *= (k + 1.0);
      dtauk *= dtau;

    } /* for k < numSpins */

  return XLAL_SUCCESS;

} /* XLALExtrapolatePulsarSpins() */


/**
 * Extrapolate phase \f$\phi_0\f$ from \f$\tau_0\f$ to \f$\tau_1\f$, given the spins \c fkdot1 at \f$\tau_1\f$.
 * Returns \f$\phi_1\f$ in the range \f$[0, 2\pi]\f$.
 */
int
XLALExtrapolatePulsarPhase ( REAL8 *phi1,			/**< [out] output phase at \f$\tau_1\f$ */
                             const PulsarSpins fkdot1,		/**< [in] spin-params at reference \f$\tau_1\f$ */
                             const REAL8 phi0,			/**< [in] input phase at \f$\tau_0\f$ */
                             const REAL8 dtau			/**< [in] time difference \f$\tau_1 - \tau_0\f$ to extrapolate \c phi0 to */
                             )
{
  UINT4 numSpins = PULSAR_MAX_SPINS;
  UINT4 k;
  UINT4 kFact;
  REAL8 dtauk;
  REAL8 frac_cycles;
  REAL8 dummy, phi;

  XLAL_CHECK( fkdot1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( phi1 != NULL, XLAL_EFAULT );

  kFact = 1;
  dtauk = 1.0;
  frac_cycles = 0;

  for ( k=0; k < numSpins; k++ )
    {
      kFact *= (k+1);
      dtauk *= -dtau;
      frac_cycles += modf ( fkdot1[k] * dtauk / kFact, &dummy );
    }

  phi = fmod ( phi0 - LAL_TWOPI * frac_cycles, LAL_TWOPI );
  if ( phi < 0 )
    phi += LAL_TWOPI;

  (*phi1) = phi;

  return XLAL_SUCCESS;

} /* XLALExtrapolatePulsarPhase() */


/**
 * Determines a frequency band which covers the frequency evolution of a band of CW signals between two GPS times.
 * The calculation accounts for the spin evolution of the signals, and the maximum possible Dopper modulation
 * due to detector motion, and (for binary signals) binary orbital motion.
 */
int
XLALCWSignalCoveringBand ( REAL8 *minCoverFreq,                          /**< [out] Minimum frequency of the covering band */
                           REAL8 *maxCoverFreq,                          /**< [out] Maximum frequency of the covering band */
                           const LIGOTimeGPS *time1,                     /**< [in] One end of the GPS time range */
                           const LIGOTimeGPS *time2,                     /**< [in] The other end of the GPS time range */
                           const PulsarSpinRange *spinRange,             /**< [in] Frequency and spindown range of the CW signals */
                           const REAL8 binaryMaxAsini,                   /**< [in] Maximum projected semi-major axis a*sini/c (= 0 for isolated sources) */
                           const REAL8 binaryMinPeriod,                  /**< [in] Minimum orbital period (s); must be 0 for isolated signals */
                           const REAL8 binaryMaxEcc                      /**< [in] Maximal binary eccentricity: must be 0 for isolated signals */
                           )
{
  // Check input
  XLAL_CHECK( minCoverFreq != NULL, XLAL_EFAULT );
  XLAL_CHECK( maxCoverFreq != NULL, XLAL_EFAULT );
  XLAL_CHECK( time1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( time2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( spinRange != NULL, XLAL_EFAULT );
  XLAL_CHECK( binaryMaxAsini >= 0, XLAL_EINVAL );

  XLAL_CHECK( (binaryMaxAsini > 0)  || ((binaryMinPeriod == 0) && (binaryMaxEcc == 0)), XLAL_EINVAL );	// if isolated: required P=0 and e=0
  XLAL_CHECK( (binaryMaxAsini == 0) || ((binaryMinPeriod > 0) && (binaryMaxEcc >= 0) && (binaryMaxEcc<1)), XLAL_EINVAL );	// if binary: P>0, 0<=e<1

  // Extrapolate frequency and spindown range to each end of GPS time range
  PulsarSpinRange spin1, spin2;
  const REAL8 DeltaT1 = XLALGPSDiff(time1, &spinRange->refTime);
  const REAL8 DeltaT2 = XLALGPSDiff(time2, &spinRange->refTime);
  XLAL_CHECK( XLALExtrapolatePulsarSpinRange(&spin1, spinRange, DeltaT1) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALExtrapolatePulsarSpinRange(&spin2, spinRange, DeltaT2) == XLAL_SUCCESS, XLAL_EFUNC );

  // Determine the minimum and maximum frequencies covered
  *minCoverFreq = MYMIN( spin1.fkdot[0], spin2.fkdot[0] );
  *maxCoverFreq = MYMAX( spin1.fkdot[0] + spin1.fkdotBand[0], spin2.fkdot[0] + spin2.fkdotBand[0] );

  // Extra frequency range needed due to detector motion, per unit frequency
  // * Maximum value of the time derivative of the diurnal and (Ptolemaic) orbital phase, plus 5% for luck
  REAL8 extraPerFreq = 1.05 * LAL_TWOPI / LAL_C_SI * ( (LAL_AU_SI/LAL_YRSID_SI) + (LAL_REARTH_SI/LAL_DAYSID_SI) );

  // Extra frequency range needed due to binary orbital motion, per unit frequency
  // Upper bound on maximum value, derived from time derivative of binary-CW phase,
  // see https://bugs.ligo.org/redmine/issues/1567
  if ( binaryMaxAsini > 0 )
    {
      REAL8 maxOmega = LAL_TWOPI / binaryMinPeriod;
      extraPerFreq += maxOmega * binaryMaxAsini / ( 1.0 - binaryMaxEcc );
  }

  // Expand frequency range
  (*minCoverFreq) *= 1.0 - extraPerFreq;
  (*maxCoverFreq) *= 1.0 + extraPerFreq;

  return XLAL_SUCCESS;

} /* XLALCWSignalCoveringBand() */
