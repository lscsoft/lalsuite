/*----------------------------------------------------------------------- 
 * 
 * File Name: ZVectorAbs.c
 * 
 * Authors: Creighton, T. D., Sintes, A. M. 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * ZVectorAbs
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Computes the magnitudes of a complex vector.
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
#include "LALConstants.h"
#include "VectorOps.h"

NRCSID(ZVECTORABSC,"$Id$");

void
LALZVectorAbs(
    LALStatus               *status,
    REAL8Vector          *out,
    const COMPLEX16Vector *in
    )
{
  COMPLEX16 *a;
  REAL8    *b;
  INT4      n;

  INITSTATUS (status, "LALZVectorAbs", ZVECTORABSC);
  
  /* Make sure the arguments are not NULL: */
  ASSERT (out, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in,  status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);

  /* Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  ASSERT (in->data,  status, VECTOROPS_ENULL, VECTOROPS_MSGENULL);
  
 /* Make sure the size is correct (invalid input size): */
  ASSERT (out->length > 0, status, VECTOROPS_ESIZE, VECTOROPS_MSGESIZE); 
  
  /* Make sure the lengths are the same (size mismatch): */
  ASSERT (in->length == out->length, status,
          VECTOROPS_ESZMM, VECTOROPS_MSGESZMM);

  a = in->data;
  b = out->data;
  n = out->length;

  while (n-- > 0)
  {
    REAL8 ar = a->re;
    REAL8 ai = a->im;

    if (fabs(ar) > fabs(ai))
    {
      REAL8 ratio = ai/ar;
      *b=fabs(ar)*sqrt(1.0 + ratio*ratio);
    }
    else
    {
      if (fabs(ai) < LAL_REAL8_MIN)
      {
         *b=0.0;   /* to avoid NaN */
      }
      else
      {
         REAL8 ratio = ar/ai;
         *b=fabs(ai)*sqrt(1.0 + ratio*ratio);
      }
    }

    ++a;
    ++b;
  }

  RETURN (status);
}

