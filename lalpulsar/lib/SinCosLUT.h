//
// Copyright (C) 2011 Bernd Machenschalk
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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301  USA
//

#ifndef _SINCOSLUT_H
#define _SINCOSLUT_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

///
/// \defgroup SinCosLUT_h Header SinCosLUT.h
/// \ingroup lalpulsar_general
/// \authors Reinhard Prix, Karl Wette
///
/// \brief fast non-vector FPU version of SinCos used in various CW codes
///

/// @{
void XLALSinCosLUTInit (void);

int XLALSinCosLUT ( REAL4 *sinx, REAL4 *cosx, REAL8 x );
int XLALSinCos2PiLUT ( REAL4 *sin2pix, REAL4 *cos2pix, REAL8 x );
int XLALSinCos2PiLUTtrimmed ( REAL4 *s, REAL4 *c, REAL8 x );
/// @}

#ifdef  __cplusplus
}
#endif

#endif // _SINCOSLUT_H
