//
// Copyright (C) 2015 Reinhard Prix
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

// ---------- INCLUDES ----------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <config.h>

#include <lal/LALConstants.h>

#define IN_VECTORMATH
#include <lal/VectorMath.h>

// ---------- local prototypes ----------
static int XLALVectorFuncf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length, float (*f)(float x) );

// ========== function definitions ==========
// ---------- Vanilla FPU versions of vector mathfun() ----------
static int
XLALVectorFuncf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length, float (*f)(float x) )
{
  for ( UINT4 i = 0; i < length; i ++ )
    {
      out[i] = (*f) ( in[i] );
    } // for i < length
  return XLAL_SUCCESS;
} // XLALVectorFuncf_FPU()

int
XLALVectorSinCosf_FPU ( REAL4 *sinx, REAL4 *cosx, const REAL4 *x, UINT4 length )
{
  for ( UINT4 i = 0; i < length; i ++ )
    {
      sinx[i] = sinf ( x[i] );
      cosx[i] = cosf ( x[i] );
    } // for i < length
  return XLAL_SUCCESS;
} // XLALVectorSinCosf_FPU()

int
XLALVectorSinCosf2PI_FPU ( REAL4 *sin2pix, REAL4 *cos2pix, const REAL4 *x, UINT4 length )
{
  for ( UINT4 i = 0; i < length; i ++ )
    {
      sin2pix[i] = sinf ( (REAL4)LAL_TWOPI * x[i] );
      cos2pix[i] = cosf ( (REAL4)LAL_TWOPI * x[i] );
    } // for i < N
  return XLAL_SUCCESS;
} // XLALVectorSinCosf2PI_FPU()

// convenience wrappers for specific functions
int
XLALVectorSinf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_FPU ( out, in, length, sinf );
}
int
XLALVectorCosf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_FPU ( out, in, length, cosf );
}
int
XLALVectorExpf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_FPU ( out, in, length, expf );
}
int
XLALVectorLogf_FPU ( REAL4 *out, const REAL4 *in, UINT4 length )
{
  return XLALVectorFuncf_FPU ( out, in, length, logf );
}
