/*
*  Copyright (C) 2007 Jolien Creighton, Alicia Sintes Olives
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

#include <complex.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

/**
 * \addtogroup VectorPolar_c
 * \author T. D. Creighton, A. M. Sintes
 *
 * \brief Convert complex vector components from rectangular coordinates to polar coordinates.
 *
 * Let \c u be an object of type ::COMPLEX8Vector, and let
 * \c a and \c b be objects of type ::REAL4Vector.
 *
 * The \ref LALCVectorAbs "LALCVectorAbs( &status, &a, &u )" function computes
 * the magnitude of a complex vector \c u:<br>
 * <tt>a.data[i] = sqrt(u.data[i].re^2 + v.data[i].im^2 )</tt>.
 *
 * The \ref LALCVectorAngle "LALCVectorAngle( &status, &a, &u )" function computes
 * the phase angle of a complex vector \c u
 * in the interval \f$[-\pi, \pi]\f$ radians:<br>
 * <tt>a.data[i] = atan2( u.data[i].im, v.data[i].re)</tt>.
 *
 * The \ref LALUnwrapREAL4Angle "LALUnwrapREAL4Angle( &status, &a, &b )" function
 * corrects the radian phase angles of a real vector  \c b
 * by adding multiples of
 * \f$\pm\pi\f$ when the absolute jumps between consecutive
 * angle elements are greater than \f$\pi\f$ radians.
 * This function detects branch cut crossings, but it can be
 * fooled by sparse, rapidly changing phase values.
 *
 * The double-precision functions are similar.
 *
 * ### Algorithm ###
 *
 * The algorithm for LALUnwrapREAL4Angle() and LALUnwrapREAL8Angle()
 * (Inspired from the MATLAP function unwrap):
 * \code
 *
 * a = in->data;
 * b = out->data;
 * n = out->length;
 *
 * cumsum = 0.0;
 * phaseI = *a;
 * *b = phaseI;
 * --n;
 *
 * while (n-- > 0)
 * {
 *   ++a;
 *   ++b;
 *   phaseII = *a;
 *   diffph = phaseII - phaseI;
 *   phaseI = phaseII;
 *
 *   cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
 *
 *   *b= phaseII + cumsum;
 * }
 *
 * \endcode
 *
 * ### Notes ###
 *
 * For the LALUnwrapREAL4Angle() and LALUnwrapREAL8Angle() functions, \c a,
 * and \c b should  not point to the same memory location (<tt>a != b</tt>).
 */
/*@{*/

/** computes the magnitudes of a vector of complex numbers */
int XLALCOMPLEX8VectorAbs( REAL4Vector *out, const COMPLEX8Vector *in )
{
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cabsf( in->data[i] );
	return 0;
}

/** computes the magnitudes of a vector of complex numbers */
int XLALCOMPLEX16VectorAbs( REAL8Vector *out, const COMPLEX16Vector *in )
{
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cabs( in->data[i] );
	return 0;
}


/** computes the arguments of a vector of complex numbers */
int XLALCOMPLEX8VectorArg( REAL4Vector *out, const COMPLEX8Vector *in )
{
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cargf( in->data[i] );
	return 0;
}

/** computes the arguments of a vector of complex numbers */
int XLALCOMPLEX16VectorArg( REAL8Vector *out, const COMPLEX16Vector *in )
{
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = carg( in->data[i] );
	return 0;
}


/**
 * corrects the radian phase angles of a real vector by adding multiples of
 * pi when the absolute jumps between consecutive angle elements are greater
 * pi radians
 */
int XLALREAL4VectorUnwrapAngle( REAL4Vector *out, const REAL4Vector *in )
{
	REAL4 prev;
	REAL4 diff;
	INT4  wrap;
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data || in->length == 0 )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
       	wrap = 0;
	prev = out->data[0] = in->data[0];
	for ( i = 1; i < in->length; ++i ) {
		diff  = in->data[i] - prev;
		prev  = in->data[i];
		wrap += (diff < -LAL_PI) - (diff > LAL_PI);
		out->data[i] = in->data[i] + wrap * LAL_TWOPI;
	}
	return 0;
}

/**
 * corrects the radian phase angles of a real vector by adding multiples of
 * pi when the absolute jumps between consecutive angle elements are greater
 * pi radians
 */
int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in )
{
	REAL8 prev;
	REAL8 diff;
	INT4  wrap;
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( XLAL_EFAULT );
	if ( ! out->data || ! in->data || in->length == 0 )
		XLAL_ERROR( XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( XLAL_EBADLEN );
       	wrap = 0;
	prev = out->data[0] = in->data[0];
	for ( i = 1; i < in->length; ++i ) {
		diff  = in->data[i] - prev;
		prev  = in->data[i];
		wrap += (diff < -LAL_PI) - (diff > LAL_PI);
		out->data[i] = in->data[i] + wrap * LAL_TWOPI;
	}
	return 0;
}


/** UNDOCUMENTED */
void
LALCVectorAbs(
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCOMPLEX8VectorAbs( out, in );

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALZVectorAbs(
    LALStatus             *status,
    REAL8Vector           *out,
    const COMPLEX16Vector *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCOMPLEX16VectorAbs( out, in );

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALCVectorAngle (
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCOMPLEX8VectorArg( out, in );

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALZVectorAngle (
    LALStatus              *status,
    REAL8Vector            *out,
    const COMPLEX16Vector  *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCOMPLEX16VectorArg( out, in );

  RETURN (status);
}

/** UNDOCUMENTED */
void
LALUnwrapREAL4Angle (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPSH_ESAME, VECTOROPSH_MSGESAME);

  XLALREAL4VectorUnwrapAngle( out, in );

  RETURN (status);
}

/** UNDOCUMENTED */
void
LALUnwrapREAL8Angle (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in
    )
{
  INITSTATUS(status);

  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in->data,  status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);

  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPSH_ESAME, VECTOROPSH_MSGESAME);

  XLALREAL8VectorUnwrapAngle( out, in );

  RETURN (status);
}
/*@}*/
