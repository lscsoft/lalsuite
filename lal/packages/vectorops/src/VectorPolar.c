/*----------------------------------------------------------------------- 
 * 
 * File Name: VectorPolar.c
 * 
 * Authors: Creighton, T. D., Sintes, A. M. 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 */ 

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

NRCSID (VECTORPOLARC, "$Id$");

void
LALCVectorAbs(
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{
  COMPLEX8 *a;
  REAL4    *b;
  INT4      n;

  INITSTATUS (status, "LALCVectorAbs", VECTORPOLARC);
  
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


void
LALZVectorAbs(
    LALStatus             *status,
    REAL8Vector           *out,
    const COMPLEX16Vector *in
    )
{
  COMPLEX16 *a;
  REAL8    *b;
  INT4      n;

  INITSTATUS (status, "LALZVectorAbs", VECTORPOLARC);
  
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


void
LALCVectorAngle (
    LALStatus            *status,
    REAL4Vector          *out,
    const COMPLEX8Vector *in
    )
{
  COMPLEX8 *a;
  REAL4    *b;
  INT4      n;

  INITSTATUS (status, "LALCVectorAngle", VECTORPOLARC);
  
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
    
    if ((fabs(ar)<LAL_REAL4_MIN) && (fabs(ai)<LAL_REAL4_MIN))
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


void
LALZVectorAngle (
    LALStatus              *status,
    REAL8Vector            *out,
    const COMPLEX16Vector  *in
    )
{
  COMPLEX16  *a;
  REAL8      *b;
  INT4        n;

  INITSTATUS (status, "LALZVectorAngle", VECTORPOLARC);
  
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


void
LALUnwrapREAL4Angle (
    LALStatus            *status,
    REAL4Vector          *out,
    const REAL4Vector    *in
    )
{
  REAL4    *a;
  REAL4    *b;
  INT4      n;
  REAL4     cumsum;
  REAL4     phaseI;
  REAL4     phaseII;
  REAL4     diffph;
  
  INITSTATUS (status, "LALUnwrapREAL4Angle", VECTORPOLARC );
  
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
          
  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPS_ESAME, VECTOROPS_MSGESAME);


  a = in->data;
  b = out->data;
  n = out->length;
  
  cumsum = 0.0;
  phaseI = *a;
  *b = phaseI;
  --n;

  while (n-- > 0)
  {
    ++a;
    ++b;
    phaseII = *a;
    diffph = phaseII - phaseI;
    phaseI = phaseII;
    
    cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
    
    *b= phaseII + cumsum;
  }

  RETURN (status);
}


void
LALUnwrapREAL8Angle (
    LALStatus            *status,
    REAL8Vector          *out,
    const REAL8Vector    *in
    )
{
  REAL8    *a;
  REAL8    *b;
  INT4      n;
  REAL8     cumsum;
  REAL8     phaseI;
  REAL8     phaseII;
  REAL8     diffph;
  
  INITSTATUS (status, "LALUnwrapREAL8Angle", VECTORPOLARC );
  
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
          
  /* Make sure the arguments are not pointing to the same memory location */
  ASSERT (out != in, status, VECTOROPS_ESAME, VECTOROPS_MSGESAME);


  a = in->data;
  b = out->data;
  n = out->length;
  
  cumsum = 0.0;
  phaseI = *a;
  *b = phaseI;
  --n;

  while (n-- > 0)
  {
    ++a;
    ++b;
    phaseII = *a;
    diffph = phaseII - phaseI;
    phaseI = phaseII;
    
    cumsum += LAL_TWOPI*( (diffph < - LAL_PI) - (diffph > LAL_PI) );
    
    *b= phaseII + cumsum;
  }

  RETURN (status);
}

