/*----------------------------------------------------------------------- 
 * 
 * File Name: DZVectorMultiply.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * DZVectorMultiply
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Multiply a double-precision real vector by a double-precision
 * complex vectors.
 * 
 * DIAGNOSTICS 
 * Null pointer, invalid input size, size mismatch.
 *
 * CALLS
 * 
 * NOTES
 * 
 *----------------------------------------------------------------------- */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/VectorOps.h>

NRCSID(DZVECTORMULTIPLYC,"$Id$");

void
LALDZVectorMultiply (
    LALStatus                *status,
    COMPLEX16Vector       *out,
    const REAL8Vector     *in1,
    const COMPLEX16Vector *in2
    )
{
  REAL8     *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "LALDZVectorMultiply", DZVECTORMULTIPLYC);

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
    REAL8 fac = *a;
    REAL8 br  = b->re;
    REAL8 bi  = b->im;

    c->re = fac*br;
    c->im = fac*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}

