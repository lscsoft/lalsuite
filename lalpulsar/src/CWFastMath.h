//
// Copyright (C) 2013 Karl Wette
// Copyright (C) 2005 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

#ifndef _CWFASTMATH_H
#define _CWFASTMATH_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup CWFastMath_h Header CWFastMath.h
/// \ingroup pkg_pulsarCommon
/// \authors Reinhard Prix, Karl Wette
///
/// \brief Various functions for performing fast math in CW functions.
///

// @{

///
/// Calculate sin(x) and cos(x) to roughly 1e-7 precision using a lookup-table and Taylor-expansion.
///
/// \note This function will fail for arguments larger than |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
///
/// Returns XLAL_SUCCESS or XLAL_FAILURE.
///
int
XLALSinCosLUT(
  REAL4 *sinx,
  REAL4 *cosx,
  REAL8 x
  );

///
/// Calculate sin(2*pi*x) and cos(2*pi*x) to roughly 1e-7 precision using a lookup-table and
/// Taylor-expansion.
///
/// \note This function will fail for arguments larger than |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
///
/// Returns XLAL_SUCCESS or XLAL_FAILURE.
///
int
XLALSinCos2PiLUT(
  REAL4 *sin2pix,
  REAL4 *cos2pix,
  REAL8 x
  );

// @}

#ifdef  __cplusplus
}
#endif

#endif // _CWFASTMATH_H
