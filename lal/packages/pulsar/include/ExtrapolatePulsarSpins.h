/*
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

/** \defgroup ExtrapolatePulsarSpins Extrapolate Pulsar Spins and Phase
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
 *
 * The initial-phase extrapolation in LALExtrapolatePulsarPhase() proceeds in the other direction, extrapolating
 * \f$\phi(\tau_0)\f$ to \f$\phi(\tau_1)\f$, where the spins are given at \f$\tau_1\f$, i.e. \f$f^{(k)}(\tau_1)\f$.
 *
 * By using the above equations, one can arrive at the following expression
 * \f[
 * \phi(\tau_1) = \phi(\tau_0) - \sum_{k=0}^s \frac{f^{(k)}(\tau_1)}{(k+1)!}\,\Delta\tau^{k+1}\,,
 * \f]
 * where <b>NOW</b> we have defined
 * \f[
 * \Delta\tau \equiv \tau_0 - \tau_1\,.
 * \f]
 *
 * This function is used in LALEstimatePulsarAmplitudeParams() to propagate the estimated initial phase
 * from the internal reference time back to the user-input reference time.
 */

/** \file
 * \ingroup ExtrapolatePulsarSpins
 * \author Reinhard Prix
 * \date $Date$
 * \brief Contains prototype for XLALExtrapolatePulsarSpins().
 */

#ifndef _EXTRAPOLATEPULSARSPINS_H  /* Double-include protection. */
#define _EXTRAPOLATEPULSARSPINS_H

#include <lal/PulsarDataTypes.h>
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

/*---------- exported Global variables ----------*/

/*---------- exported prototypes [API] ----------*/
void LALExtrapolatePulsarSpinRange( LALStatus *, PulsarSpinRange *range1, LIGOTimeGPS epoch1,  const PulsarSpinRange *range0 );

void LALExtrapolatePulsarSpins (LALStatus *, PulsarSpins fkdot1, LIGOTimeGPS epoch1, const PulsarSpins fkdot0, LIGOTimeGPS epoch0 );

void LALExtrapolatePulsarPhase (LALStatus *, REAL8 *phi1, PulsarSpins fkdot1, LIGOTimeGPS epoch1, REAL8 phi0, LIGOTimeGPS epoch0 );

int XLALExtrapolatePulsarSpins ( PulsarSpins fkdotOut, const PulsarSpins fkdotIn, REAL8 DeltaTau );

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
