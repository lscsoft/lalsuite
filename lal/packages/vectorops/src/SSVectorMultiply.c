/*----------------------------------------------------------------------- 
 * 
 * File Name: SSVectorMultiply.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * SSVectorMultiply
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Multiply two real vectors.
 * 
 * DIAGNOSTICS 
 * Null pointer, invalid input size, size mismatch.
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _VECTOROPS_H
#include "VectorOps.h"
#ifndef _VECTOROPS_H
#define _VECTOROPS_H
#endif
#endif

NRCSID (SSVECTORMULTIPLYC, "$Id$");

void
SSVectorMultiply (
    Status               *status,
    REAL4Vector          *out,
    const REAL4Vector    *in1,
    const REAL4Vector    *in2
    )
{
  REAL4 *a;
  REAL4 *b;
  REAL4 *c;
  INT4   n;

  INITSTATUS (status, SSVECTORMULTIPLYC);

  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (in2, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in1->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in2->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE);
  ASSERT (in1->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);
  ASSERT (in2->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in1->data;
  b = in2->data;
  c = out->data;
  n = out->length;

  while (n-- > 0)
  {
    *c++ = (*a++)*(*b++);
  }

  RETURN (status);
}

