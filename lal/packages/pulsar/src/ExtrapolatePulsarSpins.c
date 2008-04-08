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

/*********************************************************************************/
/** \file
 * \ingroup ExtrapolatePulsarSpins
 * \author Reinhard Prix
 * \date $Date$
 *
 * \brief Defines functions to extrapolate the pulsar spin-paramters
 * \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$ from one SSB epoch to another.
 *
 *********************************************************************************/

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

NRCSID( EXTRAPOLATEPULSARSPINSC, "$Id$");

/*---------- main functions ---------- */

/** General pulsar-spin extraploation function: given a "spin-range" (ie spins + spin-bands) at time
 * \f$\tau_0\f$, propagate the whole spin-range to time \f$\tau_1\f$.
 *
 * NOTE: *range1 is allowed to point to the same spin-range as *range0: the input will be overwritten
 *       with the output.
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

  INITSTATUS( status, "LALExtrapolatePulsarSpinRange", EXTRAPOLATEPULSARSPINSC );

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

	} /* for l < numSpins */

      range1->fkdot[l]     = flmin;
      range1->fkdotBand[l] = flmax - flmin;

    } /* for l < numSpins */

  /* set proper epoch for output */
  range1->refTime = epoch1;

  RETURN( status );

} /* ExtrapolatePulsarSpinRange() */


/** Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 *  (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epoch0)
 *  to the new reference-epoch \f$\tau_1\f$ (\a epoch1).
 *
 * This is equivalent to LALExtrapolatePulsarSpins(), but uses the fixed-size array-type
 * 'PulsarSpins' = REAL8[PULSAR_MAX_SPINS] instead, which is easier to handle and avoids
 * any dynamic-memory hassles.
 *
 * NOTE: this can be called with fkdot1 == fkdot0, in which case the input will be correctly
 * replaced by the output.
 */
void
LALExtrapolatePulsarSpins (LALStatus   *status,
			   PulsarSpins  fkdot1,		/**< [out] spin-parameters at epoch1 */
			   LIGOTimeGPS  epoch1, 	/**< [in] GPS SSB-time of new epoch1 */
			   const PulsarSpins  fkdot0,	/**< [in] spin-params at reference epoch0 */
			   LIGOTimeGPS  epoch0		/**< [in] GPS SSB-time of reference-epoch0 */
			   )
{
  UINT4 numSpins = PULSAR_MAX_SPINS; 			/* fixed size array */
  UINT4 k, l;
  REAL8 kfact, dtau, dtauk;
  PulsarSpins inSpins;

  INITSTATUS( status, "LALExtrapolatePulsarSpins", EXTRAPOLATEPULSARSPINSC );

  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  /* keep a local copy of input to allow the input- and output- pointers to be identical */
  memcpy ( inSpins, fkdot0, sizeof(PulsarSpins) );

  /* ----- translate each spin-value \f$\f^{(l)}\f$ from epoch0 to epoch1 */
  dtau = XLALDeltaFloatGPS( &epoch1, &epoch0 );

  for ( l = 0; l < numSpins; l ++ )
    fkdot1[l] = 0;

  kfact = 1;
  dtauk = 1;	/* values of k! and (dTau)^k at k=0 */
  for ( k = 0; k < numSpins; k ++ )
    {
      REAL8 kcoef  = dtauk / kfact;
      for ( l=0; l < numSpins - k ; l ++ )
	fkdot1[l] += inSpins[ k + l ] * kcoef;

      kfact *= (k+1);
      dtauk *= dtau;

    } /* for k < numSpins */

  RETURN(status);

} /* LALExtrapolateSpins2() */


/** Extrapolate phase phi0 from epoch0 to epoch1, given the spins fkdot1 at epoch1
 * Returns phi1 in the range [0, 2pi] */
void
LALExtrapolatePulsarPhase (LALStatus *status,
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

  INITSTATUS( status, "LALExtrapolatePulsarPhase", EXTRAPOLATEPULSARSPINSC );

  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( phi1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  kFact = 1;
  dTau = XLALDeltaFloatGPS( &epoch0, &epoch1 );
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

/* ============================== DEACTIVATED STUFF ============================== */
#if 0
/* ===== DEACTIVATED: use new functions built on PulsarSpins-type instead ===== */

/** Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 *  (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epoch0)
 *  to the new reference-epoch \f$\tau_1\f$ (\a epoch1).
 *
 * This is a simple wrapper to LALExtrapolatePulsarSpinRange() for the simpler case of zero-band ranges.
 * Both fkdot1 and fkdot0 must be allocated and have the same length, but ARE allowed to
 * point to the SAME vector, in which case the input gets overwritten with the output.
 */
void
LALExtrapolatePulsarSpins (LALStatus *status,
			   REAL8Vector *fkdot1,		/**< [out] spin-parameters at epoch1 */
			   LIGOTimeGPS  epoch1, 	/**< [in] GPS SSB-time of new epoch */
			   const REAL8Vector *fkdot0,	/**< [in] spin-params at reference epoch0 */
			   LIGOTimeGPS epoch0		/**< [in] GPS SSB-time of reference-epoch */
			   )
{
  UINT4 numSpins;
  UINT4 k;
  LALPulsarSpinRange *range = NULL;

  INITSTATUS( status, "LALExtrapolatePulsarSpins", EXTRAPOLATEPULSARSPINSC );
  ATTATCHSTATUSPTR ( status );

  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( fkdot0, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  numSpins = fkdot0->length;
  ASSERT ( numSpins == fkdot1->length, status, EXTRAPOLATEPULSARSPINS_ENUMSPINS, EXTRAPOLATEPULSARSPINS_MSGENUMSPINS);
  /* create temporary Range-variable for calling LALExtrapolatePulsarSpinRange() */
  if ( ( range = XLALCreatePulsarSpinRange ( numSpins )) == NULL ) {
    ABORT ( status,  EXTRAPOLATEPULSARSPINS_EMEM,  EXTRAPOLATEPULSARSPINS_MSGEMEM );
  }
  for ( k = 0; k < numSpins; k ++ )
    {
      range->fkdot->data[k] = fkdot0->data[k];
      range->fkdotBand->data[k] = 0;
    } /* for k < numSpins */
  range->epoch = epoch0;

  /* we use the (inRange == outRange) feature and have LALExtrapolatePulsarSpinRange() overwrite the input-range */
  TRY ( LALExtrapolatePulsarSpinRange ( status->statusPtr, range, epoch1, range ), status );

  /* by definition the output-range also has zero-bands => return spins at epoch1 */
  for ( k = 0; k < numSpins; k ++ )
    {
      fkdot1->data[k] = range->fkdot->data[k];
      if ( range->fkdotBand->data[k] != 0 )  /* be paranoid */
	{
	  LALPrintError ( "\nSomething gone wrong in LALExtrapolatePulsarSpinRange(): zero-band range should always remain zero-band!\n\n");
	  ABORT ( status, EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL );
	}
    }  /* for k < numSpins */

  /* we're done, free intermediate range variable */
  XLALDestroyPulsarSpinRange ( range );

  DETATCHSTATUSPTR (status );
  RETURN(status);
} /* LALExtrapolateSpins() */


/** Simple function allocating a PulsarSpinRange with \a numSpins spin-values */
LALPulsarSpinRange *
XLALCreatePulsarSpinRange ( UINT4 numSpins )
{
  LALPulsarSpinRange *range = NULL;

  if ( numSpins == 0 )
    XLAL_ERROR_NULL( "XLALCreatePulsarSpinRange", XLAL_EINVAL );

  if ( ( range = LALCalloc ( 1, sizeof( *range) )) == NULL )
    XLAL_ERROR_NULL( "XLALCreatePulsarSpinRange", XLAL_ENOMEM );

  if ( ( range->fkdot = XLALCreateREAL8Vector ( numSpins ) ) == NULL )
    {
      LALFree( range );
      XLAL_ERROR_NULL( "XLALCreatePulsarSpinRange", XLAL_ENOMEM );
    }

  if ( ( range->fkdotBand = XLALCreateREAL8Vector ( numSpins ) ) == NULL )
    {
      XLALDestroyREAL8Vector ( range->fkdot );
      LALFree( range );
      XLAL_ERROR_NULL( "XLALCreatePulsarSpinRange", XLAL_ENOMEM );
    }

  return ( range );

} /* XLALCreatePulsarSpinRange() */

/** Free memory of a LALPulsarSpinRange */
void
XLALDestroyPulsarSpinRange ( LALPulsarSpinRange *range )
{
  if ( !range )
    return;	/* no error */

  XLALDestroyREAL8Vector ( range->fkdot );
  XLALDestroyREAL8Vector ( range->fkdotBand );
  LALFree( range );

  return;

} /* XLALDestroyPulsarSpinRange() */


#endif
/* ===== DEACTIVATED CODE ===== */

