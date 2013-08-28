/*
 *
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

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define TRUE (1==1)
#define FALSE (1==0)

/*---------- main functions ---------- */

/**
 * General pulsar-spin extraploation function: given a "spin-range" (ie spins + spin-bands) at time
 * \f$\tau_0\f$, propagate the whole spin-range to time \f$\tau_1\f$.
 *
 * NOTE: *range1 is allowed to point to the same spin-range as *range0: the input will be overwritten
 * with the output.
 *
 * NOTE2: The output-range is in the 'canonical' order of \f$[ f^{(k)}, f^{(k)} + \Delta f^{(k)}]\f$,
 * where \f$\Delta f^{(k)} \ge 0\f$.
 *
 * NOTE3: This function works correctly for both epoch1 > epoch0 and epoch1 <= epoch0 !
 */
void
LALExtrapolatePulsarSpinRange(  LALStatus *status,
				PulsarSpinRange *range1,
				LIGOTimeGPS epoch1,
				const PulsarSpinRange *range0 )
{
  UINT4 numSpins = PULSAR_MAX_SPINS; 			/* fixed size array */
  UINT4 k, l;
  PulsarSpinRange inRange;
  REAL8 dtau;

  INITSTATUS(status);

  ASSERT ( range1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( range0, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  /* ----- make a copy of input range, because we allow input == output range, so
   * the input can get overwritten */
  for ( k = 0 ; k < numSpins; k ++ )
    {
      inRange.fkdot[k]     = range0->fkdot[k];
      inRange.fkdotBand[k] = range0->fkdotBand[k];
    } /* for k < numSpins */

  /* ----- translate each spin-value \f$\f^{(l)}\f$ from epoch0 to epoch1 */
  dtau = GPS2REAL8(epoch1) - GPS2REAL8(range0->refTime);

  for ( l = 0; l < numSpins; l ++ )
    {
      REAL8 flmin = 0, flmax = 0;
      REAL8 kfact = 1, dtau_powk = 1;	/* values for k=0 */

      for ( k = 0; k < numSpins - l; k ++ )
	{
	  REAL8 fkltauk0 = inRange.fkdot[k+l] * dtau_powk;
	  REAL8 fkltauk1 = fkltauk0 + inRange.fkdotBand[k+l] * dtau_powk;

	  REAL8 fkltauk_min = MYMIN ( fkltauk0, fkltauk1 );
	  REAL8 fkltauk_max = MYMAX ( fkltauk0, fkltauk1 );

	  flmin += fkltauk_min / kfact;
	  flmax += fkltauk_max / kfact;

	  kfact *= (k+1);
	  dtau_powk *= dtau;

	} /* for k < numSpins */

      range1->fkdot[l]     = flmin;
      range1->fkdotBand[l] = flmax - flmin;

    } /* for l < numSpins */

  /* set proper epoch for output */
  range1->refTime = epoch1;

  RETURN( status );

} /* ExtrapolatePulsarSpinRange() */


/**
 * Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 * (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epoch0)
 * to the new reference-epoch \f$\tau_1\f$ (\a epoch1).
 *
 * This is equivalent to LALExtrapolatePulsarSpins(), but uses the fixed-size array-type
 * 'PulsarSpins' = REAL8[PULSAR_MAX_SPINS] instead, which is easier to handle and avoids
 * any dynamic-memory hassles.
 *
 * NOTE: this can be called with fkdot1 == fkdot0, in which case the input will be correctly
 * replaced by the output.
 */
void
LALExtrapolatePulsarSpins (LALStatus   *status,		/**< pointer to LALStatus structure */
			   PulsarSpins  fkdot1,		/**< [out] spin-parameters at epoch1 */
			   LIGOTimeGPS  epoch1, 	/**< [in] GPS SSB-time of new epoch1 */
			   const PulsarSpins  fkdot0,	/**< [in] spin-params at reference epoch0 */
			   LIGOTimeGPS  epoch0		/**< [in] GPS SSB-time of reference-epoch0 */
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

} /* LALExtrapolatePulsarSpins() */

/**
 * Lightweight API to extrapolate PulsarSpins by a time-difference 'DeltaTau'.
 *
 * NOTE: This allows fkdotIn to point to the same memory as fkdotOut !
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
void
LALExtrapolatePulsarPhase (LALStatus *status,		/**< pointer to LALStatus structure */
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

  INITSTATUS(status);

  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( phi1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

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

  RETURN (status);

} /* LALExtrapolatePulsarPhase() */
