/*----------------------------------------------------------------------- 
 * 
 * File Name: CVectorAbs.c
 * 
 * Authors: Creighton, T. D., Sintes, A. M. 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * CVectorAbs
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
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

NRCSID(CVECTORABSC,"$Id$");

void
LALCVectorAbs(
    LALStatus               *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{
  COMPLEX8 *a;
  REAL4    *b;
  INT4      n;

  INITSTATUS (status, "LALCVectorAbs", CVECTORABSC);
  
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
    REAL4 ar = a->re;
    REAL4 ai = a->im;

    if (fabs(ar) > fabs(ai))
    {
      REAL4 ratio = ai/ar;
      *b=fabs(ar)*sqrt(1.0 + ratio*ratio);
    }
    else
    {
      if (fabs(ai) < LAL_REAL4_MIN)
      {
         *b=0.0;   /* to avoid NaN */
      }
      else
      {
         REAL4 ratio = ar/ai;
         *b=fabs(ai)*sqrt(1.0 + ratio*ratio);
      }
    }

    ++a;
    ++b;
  }

  RETURN (status);
}

