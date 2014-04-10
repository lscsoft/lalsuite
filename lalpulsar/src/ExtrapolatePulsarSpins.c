/*
 * Copyright (C) 2014 Karl Wette (XLALified)
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
 * \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$ from one SSB epoch to another.
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
 * General pulsar-spin extrapolation function: given a "spin-range" (ie spins + spin-bands) at time
 * \f$\tau_0\f$, propagate the whole spin-range to time \f$\tau_1\f$.
 *
 * NOTE: *range1 is allowed to point to the same spin-range as *range0: the input will be overwritten
 * with the output.
 *
 * NOTE2: The output-range is in the 'canonical' order of \f$[ f^{(k)}, f^{(k)} + \Delta f^{(k)}]\f$,
 * where \f$\Delta f^{(k)} \ge 0\f$.
 *
 */
int
XLALExtrapolatePulsarSpinRange(
  PulsarSpinRange *range1,			/**< output spin range */
  const PulsarSpinRange *range0,		/**< input spin range */
  const REAL8 dtau 				/**< time difference (range1.refTime - range0.refTime) to extrapolate spin range to */
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
 * Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 * (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epoch0)
 * to the new reference-epoch \f$\tau_1\f$ (\a epoch1).
 *
 * This is equivalent to XLALExtrapolatePulsarSpins(), but uses the fixed-size array-type
 * 'PulsarSpins' = REAL8[PULSAR_MAX_SPINS] instead, which is easier to handle and avoids
 * any dynamic-memory hassles.
 *
 * NOTE: this can be called with fkdotOut == fkdotIn, in which case the input will be correctly
 * replaced by the output.
 */
int
XLALExtrapolatePulsarSpins ( PulsarSpins fkdotOut,		/**< output fkdot array */
			     const PulsarSpins fkdotIn,		/**< inptut fkdot array */
			     REAL8 DeltaTau 			/**< time-difference (tRefOut - tRefIn) to extrapolate fkdot to */
			     )
{
  UINT4 numSpins = sizeof(PulsarSpins) / sizeof(fkdotIn[0]); 	/* fixed size array */
  UINT4 k, l;
  REAL8 kfact, dtauk;
  PulsarSpins inSpins;

  /* if DeltaTau is zero, just copy fkdotIn to fkdotOut */
  if ( DeltaTau == 0.0 ) {
    memmove ( fkdotOut, fkdotIn, sizeof(PulsarSpins) );
    return XLAL_SUCCESS;
  }

  /* keep a local copy of input to allow the input- and output- pointers to be identical */
  memcpy ( inSpins, fkdotIn, sizeof(PulsarSpins) );

  for ( l = 0; l < numSpins; l ++ )
    fkdotOut[l] = 0;

  kfact = 1;
  dtauk = 1;	/* values of k! and (dTau)^k at k=0 */
  for ( k = 0; k < numSpins; k ++ )
    {
      REAL8 kcoef  = dtauk / kfact;
      for ( l=0; l < numSpins - k ; l ++ )
	fkdotOut[l] += inSpins[ k + l ] * kcoef;

      kfact *= (k + 1.0);
      dtauk *= DeltaTau;

    } /* for k < numSpins */

  return XLAL_SUCCESS;

} /* XLALExtrapolatePulsarSpins() */


/**
 * Extrapolate phase phi0 from epoch0 to epoch1, given the spins fkdot1 at epoch1
 * Returns phi1 in the range [0, 2pi]
 */
int
XLALExtrapolatePulsarPhase(
  REAL8 *phi1,			/**< [out] phase at epoch1 */
  PulsarSpins fkdot1,		/**< [in] spin-params at reference epoch1 */
  LIGOTimeGPS epoch1, 		/**< [in] GPS SSB-time of epoch1 */
  REAL8 phi0,			/**< [in] initial phase at epoch 0 */
  LIGOTimeGPS epoch0		/**< [in] GPS SSB-time of reference-epoch */
  )
{
  UINT4 numSpins = PULSAR_MAX_SPINS;
  UINT4 k;
  UINT4 kFact;
  REAL8 dTau, dTauk;
  REAL8 frac_cycles;
  REAL8 dummy, phi;

  XLAL_CHECK( fkdot1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( phi1 != NULL, XLAL_EFAULT );

  kFact = 1;
  dTau = XLALGPSDiff( &epoch0, &epoch1 );
  dTauk = 1.0;
  frac_cycles = 0;

  for ( k=0; k < numSpins; k++ )
    {
      kFact *= (k+1);
      dTauk *= dTau;
      frac_cycles += modf ( fkdot1[k] * dTauk / kFact, &dummy );
    }

  phi = fmod ( phi0 - LAL_TWOPI * frac_cycles, LAL_TWOPI );
  if ( phi < 0 )
    phi += LAL_TWOPI;

  (*phi1) = phi;

  return XLAL_SUCCESS;

} /* LALExtrapolatePulsarPhase() */


void
LALExtrapolatePulsarSpinRange(  LALStatus *status,
				PulsarSpinRange *range1,
				LIGOTimeGPS epoch1,
				const PulsarSpinRange *range0 )
{
  REAL8 dtau;
  INITSTATUS(status);
  ASSERT ( range1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( range0, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  dtau = XLALGPSDiff( &epoch1, &range0->refTime );
  if ( XLALExtrapolatePulsarSpinRange( range1, range0, dtau ) != XLAL_SUCCESS ) {
    ABORT ( status,  EXTRAPOLATEPULSARSPINS_EXLAL,  EXTRAPOLATEPULSARSPINS_MSGEXLAL );
  }
  RETURN(status);
}


void
LALExtrapolatePulsarSpins (LALStatus   *status,
			   PulsarSpins  fkdot1,
			   LIGOTimeGPS  epoch1,
			   const PulsarSpins  fkdot0,
			   LIGOTimeGPS  epoch0
			   )
{
  REAL8 dtau;
  INITSTATUS(status);
  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  dtau = XLALGPSDiff( &epoch1, &epoch0 );
  if ( XLALExtrapolatePulsarSpins ( fkdot1, fkdot0, dtau ) ) {
    ABORT ( status,  EXTRAPOLATEPULSARSPINS_EXLAL,  EXTRAPOLATEPULSARSPINS_MSGEXLAL );
  }
  RETURN(status);
}
