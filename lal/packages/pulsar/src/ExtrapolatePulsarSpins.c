/*
 * 
 * Copyright (C) 2005 Reinhard Prix
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
#include <lal/LALError.h>
#include <lal/Date.h>

#include "ExtrapolatePulsarSpins.h"

/*---------- local DEFINES ----------*/
#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )


NRCSID( EXTRAPOLATEPULSARSPINSC, "$Id$");

/*---------- main functions ---------- */

/** Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 *  (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epochOld)
 *  to the new reference-epoch \f$\tau_1\f$ (\a epochNew).
 *
 * \note 1) both reference epochs \a epochOld, and \a epochNew refer to arrival
 * time in the solar-system barycenter (SSB).
 *
 * 2) both spin-vectors \a fkdotNew and \a fkdotOld have to be 
 * allocated and have the same length.
 *
 * 3) fkdotNew is NOT allowed to point to the same vector as fkdotOld !
 */
int
XLALExtrapolatePulsarSpins (REAL8Vector *fkdotNew,	/**< [out] spin-parameters at epochNew */
			    LIGOTimeGPS epochNew, 	/**< [in] GPS SSB-time of new reference-epoch */
			    const REAL8Vector *fkdotOld,/**< [in] spin-params at old reference epoch */
			    LIGOTimeGPS epochOld	/**< [in] GPS SSB-time of old reference-epoch */
			    )
{
  UINT4 s, len;
  REAL8 deltaT;	

  /*---------- check input ---------- */
  if ( !fkdotNew || !fkdotOld )
    {
      LALPrintError ("\nNULL Input received!\n\n");
      XLAL_ERROR ( "XLALExtrapolatePulsarSpins", XLAL_EINVAL);
    }
  len = fkdotNew->length;
  if ( len != fkdotOld->length )
    {
      LALPrintError ("\nSpindown-parameters fkdotNew and fkdotOld must have same length!\n\n");
      XLAL_ERROR ( "XLALExtrapolatePulsarSpins", XLAL_EINVAL);
    }
  if ( fkdotNew == fkdotOld )
    {
      LALPrintError ("\nSorry but fkdotNew MUST be different from fkdotOld!\n\n");
      XLAL_ERROR ( "XLALExtrapolatePulsarSpins", XLAL_EINVAL);
    }

  /* ---------- apply the pulsar spindown-model ---------- */
  deltaT = XLALDeltaFloatGPS (&epochNew, &epochOld);	/* tau_1 - tau_0 */

  for ( s = 0; s < len; s++ )
    {
      UINT4 j;
      REAL8 Tas;
      UINT4 jfact;
      fkdotNew->data[s] = 0;

      Tas = 1;
      jfact = 1;
      for ( j = 0; j < len - s; j++ )
	{
	  fkdotNew->data[s] += fkdotOld->data[j+s] * Tas / jfact;
	  Tas *= deltaT;
	  jfact *= (j + 1);
	} /* for j < spdnOrder */

    } /* for s < spdnOrder */

  return XLAL_SUCCESS;

} /* XLALExtrapolatePulsarSpins() */


/** Given a spin-vector, and a vector of spin-bands at a given reference-time,
 * return the frequency-interval spanned at time t0
 */
void
LALExtrapolatePulsarSpinRange (LALStatus *status,
			       LALFrequencyInterval *spinRange,	/**< [out] final spin-range at t0 */
			       const REAL8Vector *fkdotRef,	/**< spin-values at tRef */
			       const REAL8Vector *fkdotBand,	/**< band of spin-values */
			       LIGOTimeGPS tRef,		/**< reference-time tRef */
			       LIGOTimeGPS t0 )			/**< 'final' time t0 */
{
  REAL8Vector *fkmax = NULL;
  REAL8Vector *fkmin = NULL;
  REAL8Vector *fkbuf = NULL;
  UINT4 len, k;
  INT4 ret1=0, ret2=0;
  REAL8 fMax, fMin;

  INITSTATUS( status, "LALExtrapolatePulsarSpinRange", EXTRAPOLATEPULSARSPINSC );

  ASSERT (spinRange, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT (fkdotRef, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT (fkdotBand, status, EXTRAPOLATEPULSARSPINS_ENULL, EXTRAPOLATEPULSARSPINS_MSGENULL);
  ASSERT (fkdotRef->length, status, EXTRAPOLATEPULSARSPINS_EINPUT,  EXTRAPOLATEPULSARSPINS_MSGEINPUT);
  ASSERT (fkdotRef->length == fkdotBand->length, status,  
	  EXTRAPOLATEPULSARSPINS_EINPUT,  EXTRAPOLATEPULSARSPINS_MSGEINPUT);

  len = fkdotRef->length;

  /* allocate temporary storage */
  if ( (fkmax = XLALCreateREAL8Vector (len )) == NULL ) {
    LALPrintError ("XLALCreateREAL8Vector() failed!\n");
    ABORT (status, EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL);
  }
  if ( (fkmin = XLALCreateREAL8Vector (len )) == NULL ) {
    XLALDestroyREAL8Vector (fkmax);
    LALPrintError ("XLALCreateREAL8Vector() failed!\n");
    ABORT (status, EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL);
  }
  if ( (fkbuf = XLALCreateREAL8Vector (len )) == NULL ) {
    XLALDestroyREAL8Vector (fkmax);
    XLALDestroyREAL8Vector (fkmin);
    LALPrintError ("XLALCreateREAL8Vector() failed!\n");
    ABORT (status, EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL);
  }
  
  for ( k = 0; k < len; k++ )
    {
      fkmax->data[k] = MYMAX( fkdotRef->data[k], fkdotRef->data[k] + fkdotBand->data[k] ); 
      fkmin->data[k] = MYMIN( fkdotRef->data[k], fkdotRef->data[k] + fkdotBand->data[k] ); 
    }
  
  ret1 = XLALExtrapolatePulsarSpins (fkbuf, t0, fkmax, tRef );
  fMax = fkbuf->data[0];	/* we're only interested in the final frequency! */

  if ( ret1 == 0 )
    ret2 = XLALExtrapolatePulsarSpins (fkbuf, t0, fkmin, tRef );
  fMin = fkbuf->data[0];

  /* Free memory */
  XLALDestroyREAL8Vector (fkmax);
  XLALDestroyREAL8Vector (fkmin);
  XLALDestroyREAL8Vector (fkbuf);

  /* retarded error-checking ;) */
  if ( ret1 || ret2 ) {
    int code = xlalErrno;
    XLALClearErrno(); 
    LALPrintError ("\nSomething failed in XLALExtrapolatePulsarSpins()! (xlalErrno = %d)\n\n", code);
    ABORT (status,  EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL);
  }

  /* Return result */
  spinRange->FreqMax = fMax;
  spinRange->FreqMin = fMin;

  RETURN (status);
} /* LALExtrapolatePulsarSpinRange() */
				

/** Simple LAL-wrapper to XLALExtrapolatePulsarSpins() for more convenience.. */
void
LALExtrapolatePulsarSpins (LALStatus *status,
			   REAL8Vector *fkdotNew,	/**< [out] spin-parameters at epochNew */
			   LIGOTimeGPS epochNew, 	/**< [in] GPS SSB-time of new reference-epoch */
			   const REAL8Vector *fkdotOld,	/**< [in] spin-params at old reference epoch */
			   LIGOTimeGPS epochOld		/**< [in] GPS SSB-time of old reference-epoch */
			   )
{
  INITSTATUS( status, "LALExtrapolatePulsarSpins", EXTRAPOLATEPULSARSPINSC );

  if (XLALExtrapolatePulsarSpins (fkdotNew, epochNew, fkdotOld, epochOld) != 0 )
    {
      int code = xlalErrno;
      XLALClearErrno(); 
      LALPrintError ("\nERROR: XLALExtrapolateSpins() failed (xlalErrno = %d)!\n\n", code);
      ABORT (status,  EXTRAPOLATEPULSARSPINS_EXLAL, EXTRAPOLATEPULSARSPINS_MSGEXLAL);
    }

  RETURN(status);
} /* LALExtrapolateSpins() */

