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

/**** <lalVerbatim file="VectorPolarCV">
 * Author: T. D. Creighton, A. M. Sintes
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{VectorPolar.c}}
 *
 * Convert complex vector components from rectangular coordinates to polar
 * coordinates.
 *
 * \subsubsection*{Prototypes}
 * \input{VectorPolarCP}
 * \idx{LALCVectorAbs()}
 * \idx{LALCVectorAngle()}
 * \idx{LALUnwrapREAL4Angle()}
 * \idx{LALZVectorAbs()}
 * \idx{LALZVectorAngle()}
 * \idx{LALUnwrapREAL8Angle()}
 *
 * \subsubsection*{Description}
 *
 * Let \texttt{u} be an object of type \texttt{COMPLEX8Vector}, and let
 * \texttt{a} and \texttt{b} be objects of type \texttt{REAL4Vector}.
 *
 * The \verb:LALCVectorAbs( &status, &a, &u ): function computes
 * the magnitude of a complex vector \texttt{u}.
 * $\mbox{\texttt{a.data[i]}}=\mbox{\texttt{sqrt}} (
 * \mbox{\texttt{u.data[i].re}}^2 + \mbox{\texttt{v.data[i].im}}^2 ) $.
 *
 * The \verb:LALCVectorAngle( &status, &a, &u ): function computes
 * the phase angle of a complex vector \texttt{u}
 * in the interval $[-\pi, \pi]$ radians.\\
 * $\mbox{\texttt{a.data[i]}}=\mbox{\texttt{atan2}} (
 * \mbox{\texttt{u.data[i].im}}, \mbox{\texttt{v.data[i].re}} ) $.
 *
 * The \verb:LALUnwrapREAL4Angle( &status, &a, &b ): function
 * corrects the radian phase angles of a real vector  \texttt{b}
 * by adding multiples of
 * $\pm\pi$ when the absolute jumps between consecutive
 * angle elements are greater than $\pi$ radians.
 * This function detects branch cut crossings, but it can be
 * fooled by sparse, rapidly changing phase values.
 *
 * The double-precision functions are similar.
 *
 * \subsubsection*{Algorithm}
 *
 *
 * The algorithm for LALUnwrapREAL4Angle and LALUnwrapREAL8Angle
 * (Inspired from the MATLAP function unwrap):
 * \begin{verbatim}
 *
 *   a = in->data;
 *   b = out->data;
 *   n = out->length;
 *
 *   cumsum = 0.0;
 *   phaseI = *a;
 *   *b = phaseI;
 *   --n;
 *
 *   while (n-- > 0)
 *   {
 *     ++a;
 *     ++b;
 *     phaseII = *a;
 *     diffph = phaseII - phaseI;
 *     phaseI = phaseII;
 *
 *     cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
 *
 *     *b= phaseII + cumsum;
 *   }
 *
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * For the LALUnwrapREAL4Angle and LALUnwrapREAL8Angle functions, \texttt{a},
 * and \texttt{b} should  not point to the same memory location (\texttt{a !=
 * b}).
 *
 * \vfill{\footnotesize\input{VectorPolarCV}}
 *
 **** </lalLaTeX> */

#include <math.h>
#define LAL_USE_COMPLEX_SHORT_MACROS 1
#include <lal/LALComplex.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

NRCSID (VECTORPOLARC, "$Id$");


/** computes the magnitudes of a vector of complex numbers */
int XLALCOMPLEX8VectorAbs( REAL4Vector *out, const COMPLEX8Vector *in )
{
	static const char *func = "XLALCOMPLEX8VectorAbs";
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cabsf( in->data[i] );
	return 0;
}

/** computes the magnitudes of a vector of complex numbers */
int XLALCOMPLEX16VectorAbs( REAL8Vector *out, const COMPLEX16Vector *in )
{
	static const char *func = "XLALCOMPLEX16VectorAbs";
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cabs( in->data[i] );
	return 0;
}


/** computes the arguments of a vector of complex numbers */
int XLALCOMPLEX8VectorArg( REAL4Vector *out, const COMPLEX8Vector *in )
{
	static const char *func = "XLALCOMPLEX8VectorArg";
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = cargf( in->data[i] );
	return 0;
}

/** computes the arguments of a vector of complex numbers */
int XLALCOMPLEX16VectorArg( REAL8Vector *out, const COMPLEX16Vector *in )
{
	static const char *func = "XLALCOMPLEX16VectorArg";
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
	for ( i = 0; i < in->length; ++i )
		out->data[i] = carg( in->data[i] );
	return 0;
}


/** corrects the radian phase angles of a real vector by adding multiples of
 * pi when the absolute jumps between consecutive angle elements are greater
 * pi radians */
int XLALREAL4VectorUnwrapAngle( REAL4Vector *out, const REAL4Vector *in )
{
	static const char *func = "XLALREAL4VectorUnwrapAngle";
	REAL4 prev;
	REAL4 diff;
	INT4  wrap;
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data || in->length == 0 )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
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

/** corrects the radian phase angles of a real vector by adding multiples of
 * pi when the absolute jumps between consecutive angle elements are greater
 * pi radians */
int XLALREAL8VectorUnwrapAngle( REAL8Vector *out, const REAL8Vector *in )
{
	static const char *func = "XLALREAL8VectorUnwrapAngle";
	REAL8 prev;
	REAL8 diff;
	INT4  wrap;
	UINT4 i;
	if ( ! out || ! in )
		XLAL_ERROR( func, XLAL_EFAULT );
	if ( ! out->data || ! in->data || in->length == 0 )
		XLAL_ERROR( func, XLAL_EINVAL );
	if ( out->length != in->length )
		XLAL_ERROR( func, XLAL_EBADLEN );
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


/* <lalVerbatim file="VectorPolarCP"> */
void
LALCVectorAbs(
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALCVectorAbs", VECTORPOLARC);

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


/* <lalVerbatim file="VectorPolarCP"> */
void
LALZVectorAbs(
    LALStatus             *status,
    REAL8Vector           *out,
    const COMPLEX16Vector *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALZVectorAbs", VECTORPOLARC);

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


/* <lalVerbatim file="VectorPolarCP"> */
void
LALCVectorAngle (
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALCVectorAngle", VECTORPOLARC);

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


/* <lalVerbatim file="VectorPolarCP"> */
void
LALZVectorAngle (
    LALStatus              *status,
    REAL8Vector            *out,
    const COMPLEX16Vector  *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALZVectorAngle", VECTORPOLARC);

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

/* <lalVerbatim file="VectorPolarCP"> */
void
LALUnwrapREAL4Angle (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALUnwrapREAL4Angle", VECTORPOLARC );

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

/* <lalVerbatim file="VectorPolarCP"> */
void
LALUnwrapREAL8Angle (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALUnwrapREAL8Angle", VECTORPOLARC );

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

