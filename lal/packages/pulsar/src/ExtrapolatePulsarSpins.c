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
 * \brief Defines function to extrapolate the pulsar spin-paramters 
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

/** General pulsar-spin extraploation function: given a "spin-range" (ie spins + spin-bands) at time
 * \f$\tau_0\f$, propagate the whole spin-range to time \f$\tau_1\f$.
 *
 * NOTE: the output spin-range \a range1 must be allocated, have the same length as the input spin-range 
 * \a range0 and IS allowed to be identical with \a range0, in which case it will be overwritten with the
 * new range !
 *
 * NOTE2: The output-range is in the 'canonical' order of \f$[ f^{(k)}, f^{(k)} + \Delta f^{(k)}]\f$,
 * where \f$\Delta f^{(k)} \ge 0\f$.
 * 
 * NOTE3: This function is supposed to work correctly for both epoch1 > epoch0 and epoch1 <= epoch0 !
 */
void
LALExtrapolatePulsarSpinRange(  LALStatus *status,
			     LALPulsarSpinRange *range1,
			     LIGOTimeGPS epoch1,
			     const LALPulsarSpinRange *range0 )

{
  UINT4 numSpins;
  UINT4 k, l;
  LALPulsarSpinRange *inRange = NULL;
  REAL8 dtau;

  INITSTATUS( status, "LALExtrapolatePulsarSpins", EXTRAPOLATEPULSARSPINSC );  

  ASSERT ( range1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( range0, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( range0->fkdot, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( range0->fkdotBand, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  numSpins = range0->fkdot->length;
  ASSERT ( numSpins > 0, status, EXTRAPOLATEPULSARSPINS_ENUMSPINS, EXTRAPOLATEPULSARSPINS_MSGENUMSPINS);
  ASSERT ( range0->fkdotBand->length == numSpins, status, EXTRAPOLATEPULSARSPINS_ENUMSPINS, EXTRAPOLATEPULSARSPINS_MSGENUMSPINS);
  ASSERT ( range1->fkdot->length == numSpins, status, EXTRAPOLATEPULSARSPINS_ENUMSPINS, EXTRAPOLATEPULSARSPINS_MSGENUMSPINS);
  ASSERT ( range1->fkdotBand->length == numSpins, status, EXTRAPOLATEPULSARSPINS_ENUMSPINS, EXTRAPOLATEPULSARSPINS_MSGENUMSPINS);

  /* ----- make a copy of input range, because we allow input == output range, so 
   * the input can get overwritten */
  if ( ( inRange = XLALCreatePulsarSpinRange ( numSpins ) ) == NULL ) {
    ABORT ( status,  EXTRAPOLATEPULSARSPINS_EMEM,  EXTRAPOLATEPULSARSPINS_MSGEMEM );
  }
  for ( k = 0 ; k < numSpins; k ++ )
    {
      inRange->fkdot->data[k] = range0->fkdot->data[k];
      inRange->fkdotBand->data[k] = range0->fkdotBand->data[k];
    } /* for k < numSpins */
    
  
  /* ----- translate each spin-value \f$\f^{(l)}\f$ from epoch0 to epoch1 */
  dtau = GPS2REAL8(epoch1) - GPS2REAL8(range0->epoch); 

  for ( l = 0; l < numSpins; l ++ )
    {
      REAL8 flmin = 0, flmax = 0;
      REAL8 kfact = 1, dtau_powk = 1;	/* values for k=0 */

      for ( k = 0; k < numSpins - l; k ++ )
	{
	  REAL8 fkltauk0 = inRange->fkdot->data[k+l] * dtau_powk;
	  REAL8 fkltauk1 = fkltauk0 + inRange->fkdotBand->data[k+l] * dtau_powk;

	  REAL8 fkltauk_min = MYMIN ( fkltauk0, fkltauk1 );
	  REAL8 fkltauk_max = MYMAX ( fkltauk0, fkltauk1 );

	  flmin += fkltauk_min / kfact;
	  flmax += fkltauk_max / kfact;

	  kfact *= (k+1);
	  dtau_powk *= dtau;

	} /* for k < numSpins - l */

      range1->fkdot->data[l] = flmin;
      range1->fkdotBand->data[l] = flmax - flmin;

    } /* for l < numSpins */
  
  /* set proper epoch for output */
  range1->epoch = epoch1;

  /* free intermediate copy of input */
  XLALDestroyPulsarSpinRange ( inRange );

  RETURN( status );

} /* ExtrapolatePulsarSpinRange() */

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

/** Extrapolate phase phi0 from epoch0 to epoch1, given the spins fkdot1 at epoch1 
 * Returns phi1 in the range [0, 2pi] */
void
LALExtrapolatePulsarPhase (LALStatus *status,
			   REAL8 *phi1,			/**< [out] phase at epoch1 */
			   const REAL8Vector *fkdot1,	/**< [in] spin-params at reference epoch1 */
			   LIGOTimeGPS  epoch1, 	/**< [in] GPS SSB-time of epoch1 */
			   REAL8 phi0,			/**< [in] initial phase at epoch 0 */
			   LIGOTimeGPS epoch0		/**< [in] GPS SSB-time of reference-epoch */
			   )
{
  UINT4 numSpins;
  UINT4 k;
  UINT4 kFact;
  REAL8 dTau, dTauk;
  REAL8 frac_cycles;
  REAL8 dummy, phi;

  INITSTATUS( status, "LALExtrapolatePulsarPhase", EXTRAPOLATEPULSARSPINSC );

  ASSERT ( fkdot1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT ( phi1, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);

  numSpins = fkdot1->length;

  kFact = 1;
  dTau = XLALDeltaFloatGPS( &epoch0, &epoch1 );
  dTauk = 1.0;
  frac_cycles = 0;
  for ( k=0; k < numSpins; k++ )
    {
      kFact *= (k+1);
      dTauk *= dTau;
      frac_cycles += modf ( fkdot1->data[k] * dTauk / kFact, &dummy );
    }

  phi = fmod ( phi0 - LAL_TWOPI * frac_cycles, LAL_TWOPI );
  if ( phi < 0 )
    phi += LAL_TWOPI;

  (*phi1) = phi;
  
  RETURN (status);

} /* LALExtrapolatePulsarPhase() */
