/*
*  Copyright (C) 2013, 2014 Evan Goetz
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>

#include "../vectormath.h"

#define Relfloaterr(dx,x) (fabsf(x)>0 ? fabsf((dx)/(x)) : fabsf(dx) )
#define Relerr(dx,x) (fabs(x)>0 ? fabs((dx)/(x)) : fabs(dx) )

int main(void)
{

   //struct timespec st,end,st2,end2;

   INT4 ii, length = 100000;
   REAL4VectorAligned *floatvalues0 = NULL, *floatvalues1 = NULL, *floatvalues2 = NULL, *floatvalues3 = NULL;
   alignedREAL8Vector *doublevalues1 = NULL, *doublevalues2 = NULL, *doublevalues3 = NULL;
   alignedREAL4VectorArray *floatvalues = NULL;
   XLAL_CHECK( (floatvalues0 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatvalues1 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatvalues2 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatvalues3 = XLALCreateREAL4VectorAligned(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doublevalues1 = createAlignedREAL8Vector(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doublevalues2 = createAlignedREAL8Vector(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doublevalues3 = createAlignedREAL8Vector(length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatvalues = createAlignedREAL4VectorArray(2, length, 32)) != NULL, XLAL_EFUNC );

   for (ii=0; ii<length; ii++) {
      floatvalues1->data[ii] = (REAL4)(ii-length/2)*2.0e-3;
      doublevalues1->data[ii] = (REAL8)(ii-length/2)*2.0e-3;
      floatvalues2->data[ii] = (REAL4)(ii)*1.0e-3;
      doublevalues2->data[ii] = (REAL8)(ii)*1.0e-3;
      floatvalues3->data[ii] = (REAL4)(ii-length/2)*2.0e-4;
      doublevalues3->data[ii] = (REAL8)(ii-length/2)*2.0e-4;
   }
   memcpy(floatvalues->data[0]->data, floatvalues1->data, sizeof(REAL4)*length);
   memcpy(floatvalues->data[1]->data, floatvalues2->data, sizeof(REAL4)*length);
   memcpy(floatvalues0->data, floatvalues1->data, sizeof(REAL4)*length);

   REAL4VectorAligned *floatresult_vecsum = NULL, *floatresult_vecmult = NULL, *floatresult_addscalar = NULL, *floatresult_scale = NULL;
   alignedREAL8Vector *doubleresult_exp = NULL, *doubleresult_addscalar = NULL, *doubleresult_scale = NULL;
   alignedREAL4VectorArray *arraysumresult = NULL;
   XLAL_CHECK( (doubleresult_exp = createAlignedREAL8Vector(doublevalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_vecsum = XLALCreateREAL4VectorAligned(floatvalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_vecmult = XLALCreateREAL4VectorAligned(floatvalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_addscalar = XLALCreateREAL4VectorAligned(floatvalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_scale = XLALCreateREAL4VectorAligned(floatvalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_addscalar = createAlignedREAL8Vector(doublevalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_scale = createAlignedREAL8Vector(doublevalues1->length, 32)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (arraysumresult = createAlignedREAL4VectorArray(2, length, 32)) != NULL, XLAL_EFUNC );
   memset(arraysumresult->data[0]->data, 0, sizeof(REAL4)*length);
   memset(arraysumresult->data[1]->data, 0, sizeof(REAL4)*length);

   //clock_gettime(CLOCK_REALTIME, &st);

   XLAL_CHECK( sse_exp_REAL8Vector(doubleresult_exp, doublevalues3) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorAddREAL4(floatvalues1->data, floatvalues1->data, floatvalues1->data, floatvalues1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorMultiplyREAL4(floatresult_vecmult->data, floatvalues1->data, floatvalues2->data, floatvalues1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorShiftREAL4(floatresult_addscalar->data, (REAL4)100.0, floatvalues1->data, floatvalues1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( XLALVectorScaleREAL4(floatresult_scale->data, (REAL4)100.0, floatvalues1->data, floatvalues1->length) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseAddScalarToREAL8Vector(doubleresult_addscalar, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseScaleREAL8Vector(doubleresult_scale, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseSSVectorArraySum(arraysumresult, floatvalues, floatvalues, 0, 1, 0, 1) == XLAL_SUCCESS, XLAL_EFUNC );

   //clock_gettime(CLOCK_REALTIME, &end);

   REAL4 maxfloaterr_vecsum = 0.0, maxfloatrelerr_vecsum = 0.0, maxfloaterr_vecmult = 0.0, maxfloatrelerr_vecmult = 0.0, maxfloaterr_addscalar = 0.0, maxfloatrelerr_addscalar = 0.0, maxfloaterr_scale = 0.0, maxfloatrelerr_scale = 0.0, maxfloaterr_seqsum = 0.0, maxfloatrelerr_seqsum = 0.0;
   REAL8 maxdoubleerr_exp = 0.0, maxdoublerelerr_exp = 0.0, maxdoubleerr_addscalar = 0.0, maxdoublerelerr_addscalar = 0.0, maxdoubleerr_scale = 0.0, maxdoublerelerr_scale = 0.0;
   for (ii=0; ii<length; ii++) {
      REAL8 exp_libm = exp(doublevalues3->data[ii]);
      REAL8 doubleerr = fabs(doubleresult_exp->data[ii] - exp_libm);
      REAL8 doublerelerr = Relerr( doubleerr, exp_libm );
      maxdoubleerr_exp = fmax(doubleerr, maxdoubleerr_exp);
      maxdoublerelerr_exp = fmax(doublerelerr, maxdoublerelerr_exp);

      REAL4 sumval = (REAL4)(floatvalues0->data[ii] + floatvalues0->data[ii]);
      REAL4 floaterr = fabsf(floatvalues1->data[ii] - sumval);
      REAL4 floatrelerr = Relfloaterr( floaterr, sumval );
      maxfloaterr_vecsum = fmaxf(floaterr, maxfloaterr_vecsum);
      maxfloatrelerr_vecsum = fmaxf(floatrelerr, maxfloatrelerr_vecsum);

      REAL4 multval = (REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii]);
      floaterr = fabsf(floatresult_vecmult->data[ii] - multval);
      floatrelerr = Relfloaterr(floaterr, multval);
      maxfloaterr_vecmult = fmaxf(floaterr, maxfloaterr_vecmult);
      maxfloatrelerr_vecmult = fmaxf(floatrelerr, maxfloatrelerr_vecmult);

      sumval = (REAL4)(floatvalues1->data[ii]+(REAL4)100.0);
      REAL8 sumvald = (doublevalues1->data[ii]+100.0);
      floaterr = fabsf(floatresult_addscalar->data[ii] - sumval);
      floatrelerr = Relfloaterr(floaterr, sumval);
      doubleerr = fabs(doubleresult_addscalar->data[ii] - sumvald);
      doublerelerr = Relerr( doubleerr, sumvald );
      maxfloaterr_addscalar = fmaxf(floaterr, maxfloaterr_addscalar);
      maxfloatrelerr_addscalar = fmaxf(floatrelerr, maxfloatrelerr_addscalar);
      maxdoubleerr_addscalar = fmax(doubleerr, maxdoubleerr_addscalar);
      maxdoublerelerr_addscalar = fmax(doublerelerr, maxdoublerelerr_addscalar);

      multval = (REAL4)(floatvalues1->data[ii]*(REAL4)100.0);
      REAL8 multvald = (doublevalues1->data[ii]*100.0);
      floaterr = fabsf(floatresult_scale->data[ii] - multval);
      floatrelerr = Relfloaterr(floaterr, multval);
      doubleerr = fabs(doubleresult_scale->data[ii] - multvald);
      doublerelerr = Relerr( doubleerr, multvald );
      maxfloaterr_scale = fmaxf(floaterr, maxfloaterr_scale);
      maxfloatrelerr_scale = fmaxf(floatrelerr, maxfloatrelerr_scale);
      maxdoubleerr_scale = fmax(doubleerr, maxdoubleerr_scale);
      maxdoublerelerr_scale = fmax(doublerelerr, maxdoublerelerr_scale);

      floaterr = fabsf(arraysumresult->data[0]->data[ii] - (REAL4)(floatvalues0->data[ii]+floatvalues2->data[ii]));
      floatrelerr = Relfloaterr(floaterr, (REAL4)(floatvalues0->data[ii]+floatvalues2->data[ii]));
      maxfloaterr_seqsum = fmaxf(floaterr, maxfloaterr_seqsum);
      maxfloatrelerr_seqsum = fmaxf(floatrelerr, maxfloatrelerr_seqsum);
   }

   fprintf(stderr, "Test results SSE:\n");
   fprintf(stderr, "-----------------\n");
   fprintf(stderr, "Add REAL4Vectors: max error = %g, max relative error = %g\n", maxfloaterr_vecsum, maxfloatrelerr_vecsum);
   fprintf(stderr, "Multiply REAL4Vectors: max error = %g, max relative error = %g\n", maxfloaterr_vecmult, maxfloatrelerr_vecmult);
   fprintf(stderr, "Add scalar to REAL4Vector: max error = %g, max relative error = %g\n", maxfloaterr_addscalar, maxfloatrelerr_addscalar);
   fprintf(stderr, "Add scalar to REAL8Vector: max error = %g, max relative error = %g\n", maxdoubleerr_addscalar, maxdoublerelerr_addscalar);
   fprintf(stderr, "Scale REAL4Vector: max error = %g, max relative error = %g\n", maxfloaterr_scale, maxfloatrelerr_scale);
   fprintf(stderr, "Scale REAL8Vector: max error = %g, max relative error = %g\n", maxdoubleerr_scale, maxdoublerelerr_scale);
   fprintf(stderr, "exp(REAL8Vector): max error = %g, max relative error = %g\n", maxdoubleerr_exp, maxdoublerelerr_exp);
   fprintf(stderr, "Sum vectors of vector array into vector array: max error = %g, max relative error = %g\n", maxfloaterr_seqsum, maxfloatrelerr_seqsum);
   //fprintf(stderr, "Time elapsed: %li\n", (end.tv_sec-st.tv_sec)*GIGA+(end.tv_nsec-st.tv_nsec));

#ifdef __AVX__
   //clock_gettime(CLOCK_REALTIME, &st2);

   XLAL_CHECK( avxSSVectorSum(floatresult_vecsum, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxSSVectorMultiply(floatresult_vecmult, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxAddScalarToREAL4Vector(floatresult_addscalar, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxScaleREAL4Vector(floatresult_scale, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxAddScalarToREAL8Vector(doubleresult_addscalar, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxScaleREAL8Vector(doubleresult_scale, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxSSVectorArraySum(arraysumresult, floatvalues, floatvalues, 0, 1, 0, 1) == XLAL_SUCCESS, XLAL_EFUNC );

   //clock_gettime(CLOCK_REALTIME, &end2);

   //REAL4 maxfloaterr_vecsum = 0.0, maxfloatrelerr_vecsum = 0.0, maxfloaterr_vecmult = 0.0, maxfloatrelerr_vecmult = 0.0, maxfloaterr_addscalar = 0.0, maxfloatrelerr_addscalar = 0.0, maxfloaterr_scale = 0.0, maxfloatrelerr_scale = 0.0, maxfloaterr_seqsum = 0.0, maxfloatrelerr_seqsum = 0.0, maxfloaterr_seqsub = 0.0, maxfloatrelerr_seqsub = 0.0;
   //REAL8 maxdoubleerr_addscalar = 0.0, maxdoublerelerr_addscalar = 0.0, maxdoubleerr_scale = 0.0, maxdoublerelerr_scale = 0.0;
   for (ii=0; ii<length; ii++) {
      REAL4 floaterr = fabsf(floatresult_vecsum->data[ii] - (REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii]));
      REAL4 floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecsum->data[ii]/(REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii])));
      if (floaterr>maxfloaterr_vecsum) maxfloaterr_vecsum = floaterr;
      if (floatrelerr>maxfloatrelerr_vecsum) maxfloatrelerr_vecsum = floatrelerr;

      floaterr = fabsf(floatresult_vecmult->data[ii] - (REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecmult->data[ii]/(REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii])));
      if (floaterr>maxfloaterr_vecmult) maxfloaterr_vecmult = floaterr;
      if (floatrelerr>maxfloatrelerr_vecmult) maxfloatrelerr_vecmult = floatrelerr;

      floaterr = fabsf(floatresult_addscalar->data[ii] - (REAL4)(floatvalues1->data[ii]+(REAL4)100.0));
      REAL8 doubleerr = fabs(doubleresult_addscalar->data[ii] - (doublevalues1->data[ii]+100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_addscalar->data[ii]/(REAL4)(floatvalues1->data[ii]+(REAL4)100.0)));
      REAL8 doublerelerr = fabs(1.0 - doubleresult_addscalar->data[ii]/(doublevalues1->data[ii]+100.0));
      if (floaterr>maxfloaterr_addscalar) maxfloaterr_addscalar = floaterr;
      if (floatrelerr>maxfloatrelerr_addscalar) maxfloatrelerr_addscalar = floatrelerr;
      if (doubleerr>maxdoubleerr_addscalar) maxdoubleerr_addscalar = doubleerr;
      if (doublerelerr>maxdoublerelerr_addscalar) maxdoublerelerr_addscalar = doublerelerr;

      floaterr = fabsf(floatresult_scale->data[ii] - (REAL4)(floatvalues1->data[ii]*(REAL4)100.0));
      doubleerr = fabs(doubleresult_scale->data[ii] - (doublevalues1->data[ii]*100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_scale->data[ii]/(REAL4)(floatvalues1->data[ii]*(REAL4)100.0)));
      doublerelerr = fabs(1.0 - doubleresult_scale->data[ii]/(doublevalues1->data[ii]*100.0));
      if (floaterr>maxfloaterr_scale) maxfloaterr_scale = floaterr;
      if (floatrelerr>maxfloatrelerr_scale) maxfloatrelerr_scale = floatrelerr;
      if (doubleerr>maxdoubleerr_scale) maxdoubleerr_scale = doubleerr;
      if (doublerelerr>maxdoublerelerr_scale) maxdoublerelerr_scale = doublerelerr;

      floaterr = fabsf(arraysumresult->data[ii] - (REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - arraysumresult->data[ii]/(REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii])));
      if (floaterr>maxfloaterr_seqsum) maxfloaterr_seqsum = floaterr;
      if (floatrelerr>maxfloatrelerr_seqsum) maxfloatrelerr_seqsum = floatrelerr;
   }

   fprintf(stderr, "Test results AVX:\n");
   fprintf(stderr, "-----------------\n");
   fprintf(stderr, "Add REAL4Vectors: max error = %g, max relative error = %g\n", maxfloaterr_vecsum, maxfloatrelerr_vecsum);
   fprintf(stderr, "Multiply REAL4Vectors: max error = %g, max relative error = %g\n", maxfloaterr_vecmult, maxfloatrelerr_vecmult);
   fprintf(stderr, "Add scalar to REAL4Vector: max error = %g, max relative error = %g\n", maxfloaterr_addscalar, maxfloatrelerr_addscalar);
   fprintf(stderr, "Add scalar to REAL8Vector: max error = %g, max relative error = %g\n", maxdoubleerr_addscalar, maxdoublerelerr_addscalar);
   fprintf(stderr, "Scale REAL4Vector: max error = %g, max relative error = %g\n", maxfloaterr_scale, maxfloatrelerr_scale);
   fprintf(stderr, "Scale REAL8Vector: max error = %g, max relative error = %g\n", maxdoubleerr_scale, maxdoublerelerr_scale);
   fprintf(stderr, "Sum vectors of vector array into vector array: max error = %g, max relative error = %g\n", maxfloaterr_seqsum, maxfloatrelerr_seqsum);
   //fprintf(stderr, "Time elapsed: %li\n", (end2.tv_sec-st2.tv_sec)*GIGA+(end2.tv_nsec-st2.tv_nsec));
#endif

   XLALDestroyREAL4VectorAligned(floatvalues0);
   XLALDestroyREAL4VectorAligned(floatvalues1);
   destroyAlignedREAL8Vector(doublevalues1);
   XLALDestroyREAL4VectorAligned(floatvalues2);
   destroyAlignedREAL8Vector(doublevalues2);
   destroyAlignedREAL8Vector(doubleresult_exp);
   XLALDestroyREAL4VectorAligned(floatresult_vecsum);
   XLALDestroyREAL4VectorAligned(floatresult_vecmult);
   XLALDestroyREAL4VectorAligned(floatresult_addscalar);
   XLALDestroyREAL4VectorAligned(floatresult_scale);
   destroyAlignedREAL8Vector(doubleresult_addscalar);
   destroyAlignedREAL8Vector(doubleresult_scale);
   destroyAlignedREAL4VectorArray(floatvalues);
   destroyAlignedREAL4VectorArray(arraysumresult);

   return 0;

}




