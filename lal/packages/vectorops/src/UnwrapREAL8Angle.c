/*----------------------------------------------------------------------- 
 * 
 * File Name: UnwrapREAL8Angle.c
 * 
 * Author:   Sintes, A. M. 
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * UnwrapREAL8Angle
 * 
 * SYNOPSIS 
 *
 *
 * DESCRIPTION 
 * Corrects the radian phase angles by adding multiples of
 * + or - pi when the absolute jumps between consecutive
 * angle elements are greater than pi radians.
 * This function detects branch cut crossings, but it can be 
 * fooled by sparse, rapidly changing phase values.
 * 
 * DIAGNOSTICS 
 * Null pointer, invalid input size, size mismatch, invalid pointer.
 *
 * CALLS
 * 
 * NOTES
 * Inspired from the MATLAP function unwrap
 *----------------------------------------------------------------------- 
 */


#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/VectorOps.h>

NRCSID(UNWRAPREAL8ANGLEC,"$Id$");

void
LALUnwrapREAL8Angle (
    LALStatus               *status,
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
  
  INITSTATUS (status, "LALUnwrapREAL8Angle", UNWRAPREAL8ANGLEC );
  
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

