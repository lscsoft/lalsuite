/*----------------------------------------------------------------------- 
 * 
 * File Name: DDVectorMultiply.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * DDVectorMultiply
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Multiply two double-precision real vectors.
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

NRCSID(DDVECTORMULTIPLYC,"$Id$");

void
DDVectorMultiply (
    Status               *status,
    REAL8Vector          *out,
    const REAL8Vector    *in1,
    const REAL8Vector    *in2
    )
{
  REAL8 *a;
  REAL8 *b;
  REAL8 *c;
  INT4   n;

  INITSTATUS (status, DDVECTORMULTIPLYC);

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

