/*----------------------------------------------------------------------- 
 * 
 * File Name: CCVectorMultiplyConjugate.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * CCVectorMultiplyConjugate
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Multiply a complex vector by the complex conjugate of another vector.
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

#include <math.h>
#include "LALStdlib.h"
#include "VectorOps.h"

NRCSID (CCVECTORMULTIPLYCONJUGATEC, "$Id$");

void
CCVectorMultiplyConjugate (
    Status               *status,
    COMPLEX8Vector       *out,
    const COMPLEX8Vector *in1,
    const COMPLEX8Vector *in2
    )
{
  COMPLEX8 *a;
  COMPLEX8 *b;
  COMPLEX8 *c;
  INT4      n;

  INITSTATUS (status, "CCVectorMultiplyConjugate", CCVECTORMULTIPLYCONJUGATEC);

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
    REAL4 ar = a->re;
    REAL4 ai = a->im;
    REAL4 br = b->re;
    REAL4 bi = b->im;

    c->re = ar*br + ai*bi;
    c->im = ai*br - ar*bi;

    ++a;
    ++b;
    ++c;
  }

  RETURN (status);
}

