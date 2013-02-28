/*
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _THRESHOLDS_H
#define _THRESHOLDS_H


#include <lal/LALDatatypes.h>


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

  /** \defgroup Thresholds_h Header Thresholds.h
   * \ingroup pkg_burstsearch
   *
   * \brief UNDOCUMENTED
   */
  /*@{*/

REAL8 XLALChisqCdf(
	REAL8 chi2,
	REAL8 dof
);


REAL8 XLALOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
);


REAL8 XLALlnOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
);


REAL8 XLALNoncChisqCdf (
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
);


REAL8 XLALChi2Threshold(
	REAL8 dof,
	REAL8 falseAlarm
);


REAL8 XLALRhoThreshold(
	REAL8 chi2,
	REAL8 dof,
	REAL8 falseDismissal
);

  /*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */


#endif
