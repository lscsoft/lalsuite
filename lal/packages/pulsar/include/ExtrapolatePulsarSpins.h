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
 *  \brief  Extrapolate the Pulsar spin-paramters 
 *  \f$\{f^{(k)}\}\equiv\{f, \stackrel{.}{f},\ddot{f},...\}\f$, and "spin-ranges" 
 * \f$\{ f^{(k)}, \Delta f^{(k)} \}\f$ from one SSB epoch to another.
 * 
 * The central function of this module is LALExtrapolatePulsarSpinRange(), which extrapolates 
 * a complete "spin range" (defined as LALPulsarSpinRange) from one epoch to another.  
 * A "spin-range" contains an epoch, and \em two vectors, \f$f^{(k)}\f$ and \f$\Delta f^{(k)}\f$
 * (where "canonical" ordering refers to \f$\Delta f^{(k)} >= 0\f$ for all k.
 *
 *
 * The extrapolation is defined by the pulsar spindown-model:
 * \f[ f(\tau_1) = f(\tau_0) + {\stackrel{.}{f}(\tau_0) \over 1!} \,\Delta\tau 
 *     + {\ddot{f}(\tau_0) \over 2!} \,\Delta\tau^2 + ...
 *  = \sum_{k=0}^s {f^{(k)}(\tau_0) \over k!}\, \Delta\tau^k\,,
 * \f]
 * where \f[\Delta\tau \equiv \tau_1 - \tau_0\f\,,\f] 
 * and therefore generally
 *
 * \f[
 * f^{(l)}(\tau_1) = \sum_{k=0}^{s - l} { f^{(k+l)}(\tau_0) \over k! }\, \Delta\tau^k\,.
 * \f]
 *
 * This expression is used to extrapolate a whole "spin-range", namely at each spindown-order \f$(l)\f$ 
 * the extrapolated range is given by 
 * \f[
 * \min\left[ f^{(l)}(\tau_1) \right] = \sum_{k=0}^{s - l} {1\over k!} \min\left[ f^{(k+l)}(\tau_0) \, \Delta\tau^k \right]\,.
 * \f]
 *
 * \f[
 * \max\left[ f^{(l)}(\tau_1) \right] = \sum_{k=0}^{s - l} {1\over k!} \max\left[ f^{(k+l)}(\tau_0) \, \Delta\tau^k \right]\,.
 * \f]
 * 
 * This ensures that the range will be correctly extrapolated even if \f$\tau_1 < \tau_0\f$, i.e. \f$\Delta\tau < 0\f$.
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

NRCSID( EXTRAPOLATEPULSARSPINSH, "$Id$");

/*---------- exported INCLUDES ----------*/

/*---------- exported DEFINES ----------*/

/*----- Error-codes -----*/

#define EXTRAPOLATEPULSARSPINS_ENULL 		1
#define EXTRAPOLATEPULSARSPINS_ENONULL		2
#define EXTRAPOLATEPULSARSPINS_EMEM		3
#define EXTRAPOLATEPULSARSPINS_EINPUT		4
#define EXTRAPOLATEPULSARSPINS_ELIST		5
#define EXTRAPOLATEPULSARSPINS_EXLAL		6
#define EXTRAPOLATEPULSARSPINS_ENUMSPINS	7

#define EXTRAPOLATEPULSARSPINS_MSGENULL 	"Arguments contained an unexpected null pointer"
#define EXTRAPOLATEPULSARSPINS_MSGENONULL	"Output pointer is not NULL"
#define EXTRAPOLATEPULSARSPINS_MSGEMEM		"Out of memory"
#define EXTRAPOLATEPULSARSPINS_MSGEINPUT	"Invald input parameter"
#define EXTRAPOLATEPULSARSPINS_MSGELIST		"Error occurred in list-handling ..."
#define EXTRAPOLATEPULSARSPINS_MSGEXLAL		"(X)LAL sub-routine failed"
#define EXTRAPOLATEPULSARSPINS_MSGENUMSPINS	"Inconsistent number of spins"

/*---------- exported TYPES ----------*/

/** Contains a "spin-range", ie spins \f$f^{(k)}\f$ and corresponding bands \f$\Delta f^{(k)}\f$
 *  at a given (SSB) reference GPS-time \f$\tau\f$. 
 * "Canonical" ordering refers to \f$\Delta f^{(k)} >= 0\f$ for all k.
 */
typedef struct
{
  LIGOTimeGPS epoch;		/**< SSB reference GPS-time at which spin-range is defined */
  REAL8Vector *fkdot;		/**< Vector of spin-values \f$f^{(k)}\f$ */
  REAL8Vector *fkdotBand;	/**< Vector of spin-bands \f$\Delta f^{(k)}\f$, MUST be same length as fkdot */
} LALPulsarSpinRange;
  



/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
LALPulsarSpinRange *XLALCreatePulsarSpinRange ( UINT4 numSpins );
void XLALDestroyPulsarSpinRange ( LALPulsarSpinRange *range );

void LALExtrapolatePulsarSpinRange( LALStatus *, LALPulsarSpinRange *range1, LIGOTimeGPS epoch1, const LALPulsarSpinRange *range0 );

void LALExtrapolatePulsarSpins (LALStatus *, REAL8Vector *fkdot1, LIGOTimeGPS epoch1,
				const REAL8Vector *fkdot0, LIGOTimeGPS epoch0);

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
