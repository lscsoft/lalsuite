/*----------------------------------------------------------------------- 
 * 
 * File Name: ZZVectorDivide.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * ZZVectorDivide
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Divide a double-precision complex vector by another
 * double-precision complex vector.
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
#include "LALStdlib.h"
#include "VectorOps.h"

NRCSID(ZZVECTORDIVIDEC,"$Id$");

void
ZZVectorDivide (
    Status                *status,
    COMPLEX16Vector       *out,
    const COMPLEX16Vector *in1,
    const COMPLEX16Vector *in2
    )
{
  COMPLEX16 *a;
  COMPLEX16 *b;
  COMPLEX16 *c;
  INT4       n;

  INITSTATUS (status, "ZZVectorDivide", ZZVECTORDIVIDEC);

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
    REAL8 ar = a->re;
    REAL8 ai = a->im;
    REAL8 br = b->re;
    REAL8 bi = b->im;

    if (fabs(br) > fabs(bi))
    {
      REAL8 ratio = bi/br;
      REAL8 denom = br + ratio*bi;

      c->re = (ar + ratio*ai)/denom;
      c->im = (ai - ratio*ar)/denom;
    }
    else
    {
      REAL8 ratio = br/bi;
      REAL8 denom = bi + ratio*br;

      c->re = (ar*ratio + ai)/denom;
      c->im = (ai*ratio - ar)/denom;
    }

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}

