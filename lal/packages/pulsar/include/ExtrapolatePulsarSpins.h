/*
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

/** \defgroup ExtrapolatePulsarSpins
 *  \ingroup pulsarCommon
 *  \author Reinhard Prix
 *  \date $Date$
 *  \brief  Extrapolate the Pulsar spin-paramters \f$\{f, \stackrel{.}{f},\ddot{f},...\}\f$
 *  from one SSB epoch to another.
 * 
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
 */

/** \file 
 * \ingroup ExtrapolatePulsarSpins
 * \author Reinhard Prix
 * \date $Date$
 * \brief Contains prototype for XLALExtrapolatePulsarSpins().
 */

#ifndef _EXTRAPOLATEPULSARSPINS_H  /* Double-include protection. */
#define _EXTRAPOLATEPULSARSPINS_H

#include <lal/AVFactories.h>

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

NRCSID( EXTRAPOLATEPULSARSPINS, "$Id$");

int XLALExtrapolatePulsarSpins (REAL8Vector *fkdotNew, LIGOTimeGPS epochNew,
				const REAL8Vector *fkdotOld, LIGOTimeGPS epochOld);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
