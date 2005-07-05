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
/** \author R. Prix
 * \date July 2005 
 * \file 
 * \brief  Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 * from one SSB epoch to another.
 *
 *
 * The extrapolation method used is simply the canonical pulsar spindown-model:
 * \f[ f(\tau_1) = f(\tau_0) + {\stackrel{.}{f}(\tau_0) \over 1!} \,\Delta\tau 
 *     + {\ddot{f}(\tau_0) \over 2!} \,\Delta\tau^2 + ...\f]
 * where \f$\Delta\tau \equiv \tau_1 - \tau_0\f$, and therefore
 * \f[
 * f^{(1)}(\tau_1) = f^{(1)}(\tau_0) + {f^{(2)}(\tau_0) \over 1!}\,\Delta\tau 
 * + {f^{(3)}(\tau_0)\over 2!}\,\Delta\tau^2 +... 
 * \f]
 * \f[
 * f^{(2)}(\tau_1) = f^{(2)}(\tau_0) + {f^{(3)}(\tau_0) \over 1!}\,\Delta\tau 
 *  + {f^{(4)}(\tau_0) \over 2!}\, \Delta\tau^2 + ..
 * \f]
 *
 * so generally:
 * \f[
 * f^{(s)}(\tau_1) = \sum_{j=0} { f^{(j+s)}(\tau_0) \over j! }\, \Delta\tau^j
 * \f]
 * 
 *********************************************************************************/

/*---------- includes ---------- */
#include <lal/LALError.h>
#include <lal/Date.h>

#include <lal/ExtrapolatePulsarSpins.h>

/*---------- main functions ---------- */

/** Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 *  (\a fkdotOld) from the initial reference-epoch \f$\tau_0\f$ (\a epochOld)
 *  to the new reference-epoch \f$\tau_1\f$ (\a epochNew).
 *
 * \note both reference epochs \a epochOld, and \a epochNew refer to arrival
 * time in the solar-system barycenter (SSB).
 *
 * Note also that both spin-vectors \a fkdotNew and \a fkdotOld have to be 
 * allocated and have the same length.
 *
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
