/*----------------------------------------------------------------------- 
 * 
 * File Name: ZVectorAngle.c
 * 
 * Author:   Sintes, A. M. 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * ZVectorAngle
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Computes the phase angle of a complex vector 
 * in the interval [-pi,pi] radians.
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

NRCSID(ZVECTORANGLEC,"$Id$");

void
LALZVectorAngle (
    LALStatus                 *status,
    REAL8Vector            *out,
    const COMPLEX16Vector  *in
    )
{
  COMPLEX16  *a;
  REAL8      *b;
  INT4        n;

  INITSTATUS (status, "LALZVectorAngle", ZVECTORANGLEC);
  
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
    
    
    if ((fabs(ar)<LAL_REAL8_MIN) && (fabs(ai)<LAL_REAL8_MIN))
    {
      *b=0.0;   /* to avoid NaN when ar=ai=0.0 */
    }
    else
    {
      *b= atan2(ai,ar);
    }


    ++a;
    ++b;
  }

  RETURN (status);
}

