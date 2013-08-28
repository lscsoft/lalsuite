/*
*  Copyright (C) 2007 Jolien Creighton
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
#include <lal/VectorOps.h>

/**
 * \addtogroup VectorMultiply_c
 * \author J. D. E. Creighton, T. D. Creighton, A. M. Sintes
 *
 * \brief Multiply two vectors.
 *
 * Let \c u, \c v, and \c w be objects of type
 * ::COMPLEX8Vector, and let \c a, \c b, and \c c be
 * objects of type ::REAL4Vector.
 *
 * The \ref LALCCVectorMultiply "LALCCVectorMultiply( &status, &w, &u, &v )" function computes:<br>
 * <tt>w.data[i]= u.data[i] x v.data[i]</tt>
 *
 * The \ref LALCCVectorMultiplyConjugate "LALCCVectorMultiplyConjugate( &status, &w, &u, &v )" function computes:<br>
 * <tt>w.data[i]=u.data[i] x v.data[i]*</tt>.
 *
 * The \ref LALCCVectorDivide "LALCCVectorDivide( &status, &w, &u, &v )" function computes:<br>
 * <tt>w.data[i]= u.data[i] / v.data[i]</tt>
 *
 * The \ref LALSCVectorMultiply "LALSCVectorMultiply( &status, &w, &a, &v )" function computes:<br>
 * <tt>w.data[i]=a.data[i] x v.data[i]</tt>
 *
 * The \ref LALSSVectorMultiply "LALSSVectorMultiply( &status, &c, &a, &b )" function computes:<br>
 * <tt>c.data[i]=a.data[i] x b.data[i]</tt>
 *
 * The double-precison multiply routines (with \c D or \c Z names) work similarly.
 *
 * \heading{Algorithm}
 *
 * The algorithm for complex division is described in
 * Sec. 5.4 of Ref. [\ref ptvf1992].  The formula used is:
 * \f[
 * \frac{a + ib}{c + id} = \left\{
 * \begin{array}{ll}
 * \frac{[a + b(d/c)] + i[b - a(d/c)]}{c + d(d/c)} & |c| \ge |d| \\
 * \frac{[a(c/d) + b] + i[b(c/d) - a]}{c(c/d) + d} & |c| < |d|.
 * \end{array}
 * \right.
 * \f]
 */

/*@{*/

COMPLEX8Vector * XLALCCVectorDivide(
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ / *b++;

  return out;
}


COMPLEX16Vector * XLALZZVectorDivide(
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ / *b++;

  return out;
}


COMPLEX8Vector * XLALCCVectorMultiply(
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * *b++;

  return out;
}


COMPLEX16Vector * XLALZZVectorMultiply(
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * *b++;

  return out;
}


COMPLEX8Vector * XLALCCVectorMultiplyConjugate(
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * conjf(*b++);

  return out;
}


COMPLEX16Vector * XLALZZVectorMultiplyConjugate(
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * conj(*b++);

  return out;
}


COMPLEX8Vector * XLALSCVectorMultiply(
    COMPLEX8Vector       *out,
    const REAL4Vector    *in1,
    const COMPLEX8Vector *in2
    )
{
  REAL4    *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * *b++;

  return out;
}

COMPLEX16Vector * XLALDZVectorMultiply(
    COMPLEX16Vector       *out,
    const REAL8Vector     *in1,
    const COMPLEX16Vector *in2
    )
{
  REAL8     *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = *a++ * *b++;

  return out;
}


REAL4Vector * XLALSSVectorMultiply(
    REAL4Vector          *out,
    const REAL4Vector    *in1,
    const REAL4Vector    *in2
    )
{
  REAL4 *a;
  REAL4 *b;
  REAL4 *c;
  INT4   n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = (*a++)*(*b++);

  return out;
}


REAL8Vector * XLALDDVectorMultiply(
    REAL8Vector          *out,
    const REAL8Vector    *in1,
    const REAL8Vector    *in2
    )
{
  REAL8 *a;
  REAL8 *b;
  REAL8 *c;
  INT4   n;

  if ( ! out || ! in1 || !in2 || ! out->data || ! in1->data || ! in2->data )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( in1->length != out->length || in2->length != out->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
    *c++ = (*a++)*(*b++);

  return out;
}


/*
 *
 * LAL Routines.
 *
 */



/** UNDOCUMENTED */
void
LALCCVectorDivide (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCCVectorDivide(out, in1, in2);

  RETURN (status);
}



void
LALZZVectorDivide (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALZZVectorDivide(out, in1, in2);

  RETURN (status);
}



/** UNDOCUMENTED */
void
LALCCVectorMultiply (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCCVectorMultiply(out, in1, in2);

  RETURN (status);
}



void
LALZZVectorMultiply (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALZZVectorMultiply(out, in1, in2);

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALCCVectorMultiplyConjugate (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALCCVectorMultiplyConjugate(out, in1, in2);

  RETURN (status);
}



void
LALZZVectorMultiplyConjugate (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALZZVectorMultiplyConjugate(out, in1, in2);

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALSCVectorMultiply (
    LALStatus            *status,
    COMPLEX8Vector       *out,
    const REAL4Vector    *in1,
    const COMPLEX8Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALSCVectorMultiply(out, in1, in2);

  RETURN (status);
}



void
LALDZVectorMultiply (
    LALStatus             *status,
    COMPLEX16Vector       *out,
    const REAL8Vector     *in1,
    const COMPLEX16Vector *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALDZVectorMultiply(out, in1, in2);

  RETURN (status);
}


/** UNDOCUMENTED */
void
LALSSVectorMultiply (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in1,
    const REAL4Vector    *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALSSVectorMultiply(out, in1, in2);

  RETURN (status);
}



void
LALDDVectorMultiply (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in1,
    const REAL8Vector    *in2
    )
{
  INITSTATUS(status);

  ASSERT (out, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (in2, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in1->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);
  ASSERT (in2->data, status, VECTOROPSH_ENULL, VECTOROPSH_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPSH_ESIZE, VECTOROPSH_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPSH_ESZMM, VECTOROPSH_MSGESZMM);

  XLALDDVectorMultiply(out, in1, in2);

  RETURN (status);
}
/*@}*/
