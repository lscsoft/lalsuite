/*----------------------------------------------------------------------- 
 * 
 * File Name: VectorOpsTest.c
 * 
 * Author: Creighton, J. D. E., Sintes, A. M.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 *   main()
 *
 * SYNOPSIS 
 * 
 * DESCRIPTION 
 *   Test suite for vector operations
 * 
 * DIAGNOSTICS
 * 
 * CALLS
 *   CCVectorDivide
 *   CCVectorMultiply
 *   CCVectorMultiplyConjugate
 *   SCVectorMultiply
 *   SSVectorMultiply
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>

NRCSID (MAIN, "$Id$");

int lalDebugLevel = 2;

int
main ( void )
{
  const int size = 8;
  COMPLEX8Vector *z1 = NULL;
  COMPLEX8Vector *z2 = NULL;
  COMPLEX8Vector *z3 = NULL;
  REAL4Vector    *x1 = NULL;
  REAL4Vector    *x2 = NULL;
  REAL4Vector    *x3 = NULL;
  REAL4Vector    *y_1 = NULL;
  REAL4Vector    *y2 = NULL;
  REAL4Vector    *y3 = NULL;
  static LALStatus   status;
  INT4            i;

  LALCCreateVector(&status, &z1, size);
  LALCCreateVector(&status, &z2, size);
  LALCCreateVector(&status, &z3, size);
  LALSCreateVector(&status, &x1, size);
  LALSCreateVector(&status, &x2, size);
  LALSCreateVector(&status, &x3, size);
  LALSCreateVector(&status, &y_1, size/2);
  y2         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y2->data   = NULL;
  y2->length = size;
  y3         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y3->data   = (REAL4 *)LALMalloc(size*sizeof(REAL4));
  y3->length = 0;

  for (i = 0; i < size; ++i)
  {
    z1->data[i].re = 1 + i;
    z1->data[i].im = 2 + i*i;
    z2->data[i].re = 3 + i + i*i*i;
    z2->data[i].im = 4 + i*i + i*i*i;
    x1->data[i]    = 5 + i + i*i;
    x2->data[i]    = 6 + i + i*i + i*i*i;
  }

  printf("\n");
  LALCCVectorMultiply(&status, z3, z1, z2);
  for (i = 0; i < size; ++i)
    printf("(% 6.0f,% 6.0f) x (% 6.0f,% 6.0f) = (% 6.0f,% 6.0f)\n",
           z1->data[i].re, z1->data[i].im,
           z2->data[i].re, z2->data[i].im,
           z3->data[i].re, z3->data[i].im);

  printf("\n");
  LALCCVectorMultiplyConjugate(&status, z3, z1, z2);
  for (i = 0; i < size; ++i)
    printf("(% 6.0f,% 6.0f) x (% 6.0f,% 6.0f)* = (% 6.0f,% 6.0f)\n",
           z1->data[i].re, z1->data[i].im,
           z2->data[i].re, z2->data[i].im,
           z3->data[i].re, z3->data[i].im);

  printf("\n");
  LALCCVectorDivide(&status, z3, z1, z2);
  for (i = 0; i < size; ++i)
    printf("(% 6.0f,% 6.0f) / (% 6.0f,% 6.0f) = (% 9.6f,% 9.6f)\n",
           z1->data[i].re, z1->data[i].im,
           z2->data[i].re, z2->data[i].im,
           z3->data[i].re, z3->data[i].im);

  printf("\n");
  LALSCVectorMultiply(&status, z3, x1, z1);
  for (i = 0; i < size; ++i)
    printf("% 6.0f x (% 6.0f,% 6.0f) = (% 6.0f,% 6.0f)\n",
           x1->data[i],
           z1->data[i].re, z1->data[i].im,
           z3->data[i].re, z3->data[i].im);

  printf("\n");
  LALSSVectorMultiply(&status, x3, x1, x2);
  for (i = 0; i < size; ++i)
    printf("% 6.0f x % 6.0f = % 6.0f\n",
           x1->data[i], x2->data[i], x3->data[i]);

  printf("\n");
  LALSSVectorMultiply(&status, x3, x1, NULL);
  LALSSVectorMultiply(&status, x3, y2, x2);
  LALSSVectorMultiply(&status, y3, x1, x2);
  LALSSVectorMultiply(&status, x3, x1, y_1);

  LALCDestroyVector(&status, &z1);
  LALCDestroyVector(&status, &z2);
  LALCDestroyVector(&status, &z3);
  LALSDestroyVector(&status, &x1);
  LALSDestroyVector(&status, &x2);
  LALSDestroyVector(&status, &x3);
  LALSDestroyVector(&status, &y_1);
  LALFree(y2);
  LALFree(y3->data);
  LALFree(y3);

  x1 = x2 = x3 = y_1 = y2 = y3 = NULL;
  z1 = z2 = z3 = NULL;
  
  
  LALCCreateVector(&status, &z1, size);
  
  LALSCreateVector(&status, &x1, size);
  LALSCreateVector(&status, &x2, size);
  LALSCreateVector(&status, &x3, size);
  
  
   for (i = 0; i < size; ++i)
  {
    z1->data[i].re = (12.0 + i) *cos(LAL_PI/3.0*i);
    z1->data[i].im = (12.0 + i )*sin(LAL_PI/3.0*i);
  }
  
   printf("\n"); 
   LALCVectorAbs(&status, x1, z1);  
   for (i = 0; i < size; ++i)
    printf(" Abs(% f,%f)  = %f \n",
           z1->data[i].re, z1->data[i].im,
           x1->data[i]);
           
    LALCVectorAngle(&status, x2, z1);
    for (i = 0; i < size; ++i)    
     printf(" Angle(%f,%f)  = %f \n",
           z1->data[i].re, z1->data[i].im,
           x2->data[i]);
 
    LALUnwrapREAL4Angle(&status, x3, x2); 
     for (i = 0; i < size; ++i)    
     printf(" Unwrap Phase Angle ( %f )  = %f \n",
           x2->data[i],
           x3->data[i]);
  
  
  LALSCreateVector(&status, &y_1, size/2);
  
  y2         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y2->data   = NULL;
  y2->length = size;
  
  y3         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y3->data   = (REAL4 *)LALMalloc(size*sizeof(REAL4));
  y3->length = 0;

  printf("\n");
  
    LALCVectorAbs(&status, x1, NULL);
    LALCVectorAbs(&status, NULL, z1);
    LALCVectorAbs(&status, y_1, z1);
    LALCVectorAbs(&status, y2, z1);
    LALCVectorAbs(&status, y3, z1);
    
    
    LALCVectorAngle(&status, x2, NULL);
    LALCVectorAngle(&status, NULL, z1);
    LALCVectorAngle(&status, y_1, z1);
    LALCVectorAngle(&status, y2, z1);
    LALCVectorAngle(&status, y3, z1);
    
    LALUnwrapREAL4Angle(&status, x3, NULL);   
    LALUnwrapREAL4Angle(&status, NULL, x2);   
    LALUnwrapREAL4Angle(&status, y_1, x2);   
    LALUnwrapREAL4Angle(&status, y2, x2);   
    LALUnwrapREAL4Angle(&status, y3, x2);   
    LALUnwrapREAL4Angle(&status, x2, x2);   
 

  LALCDestroyVector(&status, &z1);
 
  LALSDestroyVector(&status, &x1);
  LALSDestroyVector(&status, &x2);
  LALSDestroyVector(&status, &x3);
  
  LALSDestroyVector(&status, &y_1);
  LALFree(y2);
  LALFree(y3->data);
  LALFree(y3);

  return 0;
}
