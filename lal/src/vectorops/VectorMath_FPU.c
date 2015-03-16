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
static int XLALVectorFuncf_FPU ( REAL4Vector *out, const REAL4Vector *in, float (*f)(float x) );

// ========== function definitions ==========
// ---------- Vanilla FPU versions of vector mathfun() ----------
static int
XLALVectorFuncf_FPU ( REAL4Vector *out, const REAL4Vector *in, float (*f)(float x) )
{
  XLAL_CHECK ( (out != NULL) && (in != NULL) && (f != NULL), XLAL_EINVAL );
  XLAL_CHECK ( in->length == out->length, XLAL_EINVAL );

  UINT4 N = in->length;
  for ( UINT4 i = 0; i < N; i ++ )
    {
      out->data[i] = (*f) ( in->data[i] );
    } // for i < N
  return XLAL_SUCCESS;
} // XLALVectorFuncf_FPU()

int
XLALVectorSinCosf_FPU ( REAL4Vector *sinx, REAL4Vector *cosx, const REAL4Vector *x )
{
  XLAL_CHECK ( (sinx != NULL) && (cosx != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sinx->length == cosx->length) && ( x->length == sinx->length), XLAL_EINVAL );

  UINT4 N = x->length;
  for ( UINT4 i = 0; i < N; i ++ )
    {
      sinx->data[i] = sinf ( x->data[i] );
      cosx->data[i] = cosf ( x->data[i] );
    } // for i < N
  return XLAL_SUCCESS;
} // XLALVectorSinCosf_FPU()

int
XLALVectorSinCosf2PI_FPU ( REAL4Vector *sin2pix, REAL4Vector *cos2pix, const REAL4Vector *x )
{
  XLAL_CHECK ( (sin2pix != NULL) && (cos2pix != NULL) && (x != NULL), XLAL_EINVAL );
  XLAL_CHECK ( (sin2pix->length == cos2pix->length) && ( x->length == sin2pix->length), XLAL_EINVAL );

  UINT4 N = x->length;
  for ( UINT4 i = 0; i < N; i ++ )
    {
      sin2pix->data[i] = sinf ( (REAL4)LAL_TWOPI * x->data[i] );
      cos2pix->data[i] = cosf ( (REAL4)LAL_TWOPI * x->data[i] );
    } // for i < N
  return XLAL_SUCCESS;
} // XLALVectorSinCosf2PI_FPU()

// convenience wrappers for specific functions
int
XLALVectorSinf_FPU ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_FPU ( out, in, sinf );
}
int
XLALVectorCosf_FPU ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_FPU ( out, in, cosf );
}
int
XLALVectorExpf_FPU ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_FPU ( out, in, expf );
}
int
XLALVectorLogf_FPU ( REAL4Vector *out, const REAL4Vector *in )
{
  return XLALVectorFuncf_FPU ( out, in, logf );
}
