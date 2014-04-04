

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>

#include "../vectormath.h"

int main(void)
{

   INT4 ii, length = 10000000;
   REAL4Vector *floatvalues1 = NULL, *floatvalues2 = NULL;
   REAL8Vector *doublevalues1 = NULL, *doublevalues2 = NULL;
   XLAL_CHECK( (floatvalues1 = XLALCreateREAL4Vector(length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatvalues2 = XLALCreateREAL4Vector(length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doublevalues1 = XLALCreateREAL8Vector(length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doublevalues2 = XLALCreateREAL8Vector(length)) != NULL, XLAL_EFUNC );
   REAL4VectorSequence *floatvalues = NULL;
   XLAL_CHECK( (floatvalues = XLALCreateREAL4VectorSequence(2, length)) != NULL, XLAL_EFUNC );

   for (ii=0; ii<length; ii++) {
      floatvalues1->data[ii] = (REAL4)(ii-length/2)*1.0e-5;
      doublevalues1->data[ii] = (REAL8)(ii-length/2)*1.0e-5;
      floatvalues2->data[ii] = (REAL4)ii*1.0e-5;
      doublevalues2->data[ii] = (REAL8)ii*1.0e-5;
   }
   memcpy(floatvalues->data, floatvalues1->data, sizeof(REAL4)*length);
   memcpy(&(floatvalues->data[length]), floatvalues2->data, sizeof(REAL4)*length);

   REAL4Vector *floatresult_exp = NULL, *floatresult_vecsum = NULL, *floatresult_vecmult = NULL, *floatresult_addscalar = NULL, *floatresult_scale = NULL, *floatresult_sin2pix = NULL, *floatresult_cos2pix = NULL, *sequencesubtractresult = NULL;
   REAL8Vector *doubleresult_exp = NULL, *doubleresult_addscalar = NULL, *doubleresult_scale = NULL, *doubleresult_sin2pix = NULL, *doubleresult_cos2pix = NULL;
   REAL4VectorSequence *sequencesumresult = NULL;
   XLAL_CHECK( (floatresult_exp = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_exp = XLALCreateREAL8Vector(doublevalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_vecsum = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_vecmult = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_addscalar = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_scale = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_addscalar = XLALCreateREAL8Vector(doublevalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_scale = XLALCreateREAL8Vector(doublevalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_sin2pix = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (floatresult_cos2pix = XLALCreateREAL4Vector(floatvalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_sin2pix = XLALCreateREAL8Vector(doublevalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (doubleresult_cos2pix = XLALCreateREAL8Vector(doublevalues1->length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (sequencesubtractresult = XLALCreateREAL4Vector(length)) != NULL, XLAL_EFUNC );
   XLAL_CHECK( (sequencesumresult = XLALCreateREAL4VectorSequence(2, length)) != NULL, XLAL_EFUNC );
   memset(sequencesumresult->data, 0, sizeof(REAL4)*2*length);

   XLAL_CHECK( sse_exp_REAL4Vector(floatresult_exp, floatvalues1) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sse_exp_REAL8Vector(doubleresult_exp, doublevalues1) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseSSVectorSum(floatresult_vecsum, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseSSVectorMultiply(floatresult_vecmult, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseAddScalarToREAL4Vector(floatresult_addscalar, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseScaleREAL4Vector(floatresult_scale, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseAddScalarToREAL8Vector(doubleresult_addscalar, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseScaleREAL8Vector(doubleresult_scale, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sse_sin_cos_2PI_LUT_REAL4Vector(floatresult_sin2pix, floatresult_cos2pix, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sse_sin_cos_2PI_LUT_REAL8Vector(doubleresult_sin2pix, doubleresult_cos2pix, doublevalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseSSVectorSequenceSum(sequencesumresult, floatvalues, floatvalues, 0, 1, 0, 1) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( sseSSVectorSequenceSubtract(sequencesubtractresult, floatvalues, floatvalues2, 0) == XLAL_SUCCESS, XLAL_EFUNC );

   REAL4 maxfloatdiff_exp = 0.0, maxfloatrelerr_exp = 0.0, maxfloatdiff_vecsum = 0.0, maxfloatrelerr_vecsum = 0.0, maxfloatdiff_vecmult = 0.0, maxfloatrelerr_vecmult = 0.0, maxfloatdiff_addscalar = 0.0, maxfloatrelerr_addscalar = 0.0, maxfloatdiff_scale = 0.0, maxfloatrelerr_scale = 0.0, maxfloatdiff_sin2pix = 0.0, maxfloatrelerr_sin2pix = 0.0, maxfloatdiff_cos2pix = 0.0, maxfloatrelerr_cos2pix = 0.0, maxfloatdiff_seqsum = 0.0, maxfloatrelerr_seqsum = 0.0, maxfloatdiff_seqsub = 0.0, maxfloatrelerr_seqsub = 0.0;
   REAL8 maxdoublediff_exp = 0.0, maxdoublerelerr_exp = 0.0, maxdoublediff_addscalar = 0.0, maxdoublerelerr_addscalar = 0.0, maxdoublediff_scale = 0.0, maxdoublerelerr_scale = 0.0, maxdoublediff_sin2pix = 0.0, maxdoublerelerr_sin2pix = 0.0, maxdoublediff_cos2pix = 0.0, maxdoublerelerr_cos2pix = 0.0;
   for (ii=0; ii<length; ii++) {
      REAL4 floatdiff = fabsf(floatresult_exp->data[ii] - expf((REAL4)floatvalues1->data[ii]));
      REAL4 floatrelerr = fabsf((REAL4)(1.0 - floatresult_exp->data[ii]/expf((REAL4)floatvalues1->data[ii])));
      REAL8 doublediff = fabs(doubleresult_exp->data[ii] - exp(doublevalues1->data[ii]));
      REAL8 doublerelerr = fabs(1.0 - doubleresult_exp->data[ii]/exp(doublevalues1->data[ii]));
      if (floatdiff>maxfloatdiff_exp) maxfloatdiff_exp = floatdiff;
      if (floatrelerr>maxfloatrelerr_exp) maxfloatrelerr_exp = floatrelerr;
      if (doublediff>maxdoublediff_exp) maxdoublediff_exp = doublediff;
      if (doublerelerr>maxdoublerelerr_exp) maxdoublerelerr_exp = doublerelerr;

      floatdiff = fabsf(floatresult_vecsum->data[ii] - (REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecsum->data[ii]/(REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_vecsum) maxfloatdiff_vecsum = floatdiff;
      if (floatrelerr>maxfloatrelerr_vecsum) maxfloatrelerr_vecsum = floatrelerr;

      floatdiff = fabsf(floatresult_vecmult->data[ii] - (REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecmult->data[ii]/(REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_vecmult) maxfloatdiff_vecmult = floatdiff;
      if (floatrelerr>maxfloatrelerr_vecmult) maxfloatrelerr_vecmult = floatrelerr;

      floatdiff = fabsf(floatresult_addscalar->data[ii] - (REAL4)(floatvalues1->data[ii]+(REAL4)100.0));
      doublediff = fabs(doubleresult_addscalar->data[ii] - (doublevalues1->data[ii]+100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_addscalar->data[ii]/(REAL4)(floatvalues1->data[ii]+(REAL4)100.0)));
      doublerelerr = fabs(1.0 - doubleresult_addscalar->data[ii]/(doublevalues1->data[ii]+100.0));
      if (floatdiff>maxfloatdiff_addscalar) maxfloatdiff_addscalar = floatdiff;
      if (floatrelerr>maxfloatrelerr_addscalar) maxfloatrelerr_addscalar = floatrelerr;
      if (doublediff>maxdoublediff_addscalar) maxdoublediff_addscalar = doublediff;
      if (doublerelerr>maxdoublerelerr_addscalar) maxdoublerelerr_addscalar = doublerelerr;

      floatdiff = fabsf(floatresult_scale->data[ii] - (REAL4)(floatvalues1->data[ii]*(REAL4)100.0));
      doublediff = fabs(doubleresult_scale->data[ii] - (doublevalues1->data[ii]*100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_scale->data[ii]/(REAL4)(floatvalues1->data[ii]*(REAL4)100.0)));
      doublerelerr = fabs(1.0 - doubleresult_scale->data[ii]/(doublevalues1->data[ii]*100.0));
      if (floatdiff>maxfloatdiff_scale) maxfloatdiff_scale = floatdiff;
      if (floatrelerr>maxfloatrelerr_scale) maxfloatrelerr_scale = floatrelerr;
      if (doublediff>maxdoublediff_scale) maxdoublediff_scale = doublediff;
      if (doublerelerr>maxdoublerelerr_scale) maxdoublerelerr_scale = doublerelerr;

      floatdiff = fabsf(floatresult_sin2pix->data[ii] - sinf((REAL4)LAL_TWOPI * floatvalues1->data[ii]));
      doublediff = fabs(doubleresult_sin2pix->data[ii] - sin(LAL_TWOPI * doublevalues1->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_sin2pix->data[ii]/sinf((REAL4)LAL_TWOPI * floatvalues1->data[ii])));
      doublerelerr = fabs(1.0 - doubleresult_sin2pix->data[ii]/sin(LAL_TWOPI * doublevalues1->data[ii]));
      if (floatdiff>maxfloatdiff_sin2pix) maxfloatdiff_sin2pix = floatdiff;
      if (floatrelerr>maxfloatrelerr_sin2pix) maxfloatrelerr_sin2pix = floatrelerr;
      if (doublediff>maxdoublediff_sin2pix) maxdoublediff_sin2pix = doublediff;
      if (doublerelerr>maxdoublerelerr_sin2pix) maxdoublerelerr_sin2pix = doublerelerr;

      floatdiff = fabsf(floatresult_cos2pix->data[ii] - cosf((REAL4)LAL_TWOPI * floatvalues1->data[ii]));
      doublediff = fabs(doubleresult_cos2pix->data[ii] - cos(LAL_TWOPI * doublevalues1->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_cos2pix->data[ii]/cosf((REAL4)LAL_TWOPI * floatvalues1->data[ii])));
      doublerelerr = fabs(1.0 - doubleresult_cos2pix->data[ii]/cos(LAL_TWOPI * doublevalues1->data[ii]));
      if (floatdiff>maxfloatdiff_cos2pix) maxfloatdiff_cos2pix = floatdiff;
      if (floatrelerr>maxfloatrelerr_cos2pix) maxfloatrelerr_cos2pix = floatrelerr;
      if (doublediff>maxdoublediff_cos2pix) maxdoublediff_cos2pix = doublediff;
      if (doublerelerr>maxdoublerelerr_cos2pix) maxdoublerelerr_cos2pix = doublerelerr;

      floatdiff = fabsf(sequencesumresult->data[ii] - (REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - sequencesumresult->data[ii]/(REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_seqsum) maxfloatdiff_seqsum = floatdiff;
      if (floatrelerr>maxfloatrelerr_seqsum) maxfloatrelerr_seqsum = floatrelerr;

      floatdiff = fabsf(sequencesubtractresult->data[ii] - (REAL4)(floatvalues->data[ii]-floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - sequencesubtractresult->data[ii]/(REAL4)(floatvalues->data[ii]-floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_seqsub) maxfloatdiff_seqsub = floatdiff;
      if (floatrelerr>maxfloatrelerr_seqsub) maxfloatrelerr_seqsub = floatrelerr;
   }

   fprintf(stderr, "Test results SSE:\n");
   fprintf(stderr, "-----------------\n");
   fprintf(stderr, "Add REAL4Vectors: max error = %g, max relative error = %g\n", maxfloatdiff_vecsum, maxfloatrelerr_vecsum);
   fprintf(stderr, "Multiply REAL4Vectors: max error = %g, max relative error = %g\n", maxfloatdiff_vecmult, maxfloatrelerr_vecmult);
   fprintf(stderr, "Add scalar to REAL4Vector: max error = %g, max relative error = %g\n", maxfloatdiff_addscalar, maxfloatrelerr_addscalar);
   fprintf(stderr, "Add scalar to REAL8Vector: max error = %g, max relative error = %g\n", maxdoublediff_addscalar, maxdoublerelerr_addscalar);
   fprintf(stderr, "Scale REAL4Vector: max error = %g, max relative error = %g\n", maxfloatdiff_scale, maxfloatrelerr_scale);
   fprintf(stderr, "Scale REAL8Vector: max error = %g, max relative error = %g\n", maxdoublediff_scale, maxdoublerelerr_scale);
   fprintf(stderr, "exp(REAL4Vector): max error = %g, max relative error = %g\n", maxfloatdiff_exp, maxfloatrelerr_exp);
   fprintf(stderr, "exp(REAL8Vector): max error = %g, max relative error = %g\n", maxdoublediff_exp, maxdoublerelerr_exp);
   fprintf(stderr, "sin(2*pi*REAL4Vector): max error = %g, max relative error = %g\n", maxfloatdiff_sin2pix, maxfloatrelerr_sin2pix);
   fprintf(stderr, "sin(2*pi*REAL8Vector): max error = %g, max relative error = %g\n", maxdoublediff_sin2pix, maxdoublerelerr_sin2pix);
   fprintf(stderr, "cos(2*pi*REAL4Vector): max error = %g, max relative error = %g\n", maxfloatdiff_cos2pix, maxfloatrelerr_cos2pix);
   fprintf(stderr, "cos(2*pi*REAL8Vector): max error = %g, max relative error = %g\n", maxdoublediff_cos2pix, maxdoublerelerr_cos2pix);
   fprintf(stderr, "Sum vectors of vector sequence into vector sequence: max error = %g, max relative error = %g\n", maxfloatdiff_seqsum, maxfloatrelerr_seqsum);
   fprintf(stderr, "Subtract vector from vector sequence: max error = %g, max relative error = %g\n", maxfloatdiff_seqsub, maxfloatrelerr_seqsub);

#ifdef __AVX__
   XLAL_CHECK( avxSSVectorSum(floatresult_vecsum, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxSSVectorMultiply(floatresult_vecmult, floatvalues1, floatvalues2) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxAddScalarToREAL4Vector(floatresult_addscalar, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxScaleREAL4Vector(floatresult_scale, floatvalues1, (REAL4)100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxAddScalarToREAL8Vector(doubleresult_addscalar, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxScaleREAL8Vector(doubleresult_scale, doublevalues1, 100.0) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxSSVectorSequenceSum(sequencesumresult, floatvalues, floatvalues, 0, 1, 0, 1) == XLAL_SUCCESS, XLAL_EFUNC );
   XLAL_CHECK( avxSSVectorSequenceSubtract(sequencesubtractresult, floatvalues, floatvalues2, 0) == XLAL_SUCCESS, XLAL_EFUNC );
   REAL4 maxfloatdiff_vecsum = 0.0, maxfloatrelerr_vecsum = 0.0, maxfloatdiff_vecmult = 0.0, maxfloatrelerr_vecmult = 0.0, maxfloatdiff_addscalar = 0.0, maxfloatrelerr_addscalar = 0.0, maxfloatdiff_scale = 0.0, maxfloatrelerr_scale = 0.0, maxfloatdiff_seqsum = 0.0, maxfloatrelerr_seqsum = 0.0, maxfloatdiff_seqsub = 0.0, maxfloatrelerr_seqsub = 0.0;
   REAL8 maxdoublediff_addscalar = 0.0, maxdoublerelerr_addscalar = 0.0, maxdoublediff_scale = 0.0, maxdoublerelerr_scale = 0.0;
   for (ii=0; ii<length; ii++) {
      REAL4 floatdiff = fabsf(floatresult_vecsum->data[ii] - (REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii]));
      REAL4 floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecsum->data[ii]/(REAL4)(floatvalues1->data[ii] + floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_vecsum) maxfloatdiff_vecsum = floatdiff;
      if (floatrelerr>maxfloatrelerr_vecsum) maxfloatrelerr_vecsum = floatrelerr;

      floatdiff = fabsf(floatresult_vecmult->data[ii] - (REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_vecmult->data[ii]/(REAL4)(floatvalues1->data[ii]*floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_vecmult) maxfloatdiff_vecmult = floatdiff;
      if (floatrelerr>maxfloatrelerr_vecmult) maxfloatrelerr_vecmult = floatrelerr;

      floatdiff = fabsf(floatresult_addscalar->data[ii] - (REAL4)(floatvalues1->data[ii]+(REAL4)100.0));
      REAL8 doublediff = fabs(doubleresult_addscalar->data[ii] - (doublevalues1->data[ii]+100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_addscalar->data[ii]/(REAL4)(floatvalues1->data[ii]+(REAL4)100.0)));
      REAL8 doublerelerr = fabs(1.0 - doubleresult_addscalar->data[ii]/(doublevalues1->data[ii]+100.0));
      if (floatdiff>maxfloatdiff_addscalar) maxfloatdiff_addscalar = floatdiff;
      if (floatrelerr>maxfloatrelerr_addscalar) maxfloatrelerr_addscalar = floatrelerr;
      if (doublediff>maxdoublediff_addscalar) maxdoublediff_addscalar = doublediff;
      if (doublerelerr>maxdoublerelerr_addscalar) maxdoublerelerr_addscalar = doublerelerr;

      floatdiff = fabsf(floatresult_scale->data[ii] - (REAL4)(floatvalues1->data[ii]*(REAL4)100.0));
      doublediff = fabs(doubleresult_scale->data[ii] - (doublevalues1->data[ii]*100.0));
      floatrelerr = fabsf((REAL4)(1.0 - floatresult_scale->data[ii]/(REAL4)(floatvalues1->data[ii]*(REAL4)100.0)));
      doublerelerr = fabs(1.0 - doubleresult_scale->data[ii]/(doublevalues1->data[ii]*100.0));
      if (floatdiff>maxfloatdiff_scale) maxfloatdiff_scale = floatdiff;
      if (floatrelerr>maxfloatrelerr_scale) maxfloatrelerr_scale = floatrelerr;
      if (doublediff>maxdoublediff_scale) maxdoublediff_scale = doublediff;
      if (doublerelerr>maxdoublerelerr_scale) maxdoublerelerr_scale = doublerelerr;

      floatdiff = fabsf(sequencesumresult->data[ii] - (REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - sequencesumresult->data[ii]/(REAL4)(floatvalues1->data[ii]+floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_seqsum) maxfloatdiff_seqsum = floatdiff;
      if (floatrelerr>maxfloatrelerr_seqsum) maxfloatrelerr_seqsum = floatrelerr;

      floatdiff = fabsf(sequencesubtractresult->data[ii] - (REAL4)(floatvalues->data[ii]-floatvalues2->data[ii]));
      floatrelerr = fabsf((REAL4)(1.0 - sequencesubtractresult->data[ii]/(REAL4)(floatvalues->data[ii]-floatvalues2->data[ii])));
      if (floatdiff>maxfloatdiff_seqsub) maxfloatdiff_seqsub = floatdiff;
      if (floatrelerr>maxfloatrelerr_seqsub) maxfloatrelerr_seqsub = floatrelerr;
   }

   fprintf(stderr, "Test results AVX:\n");
   fprintf(stderr, "-----------------\n");
   fprintf(stderr, "Add REAL4Vectors: max error = %g, max relative error = %g\n", maxfloatdiff_vecsum, maxfloatrelerr_vecsum);
   fprintf(stderr, "Multiply REAL4Vectors: max error = %g, max relative error = %g\n", maxfloatdiff_vecmult, maxfloatrelerr_vecmult);
   fprintf(stderr, "Add scalar to REAL4Vector: max error = %g, max relative error = %g\n", maxfloatdiff_addscalar, maxfloatrelerr_addscalar);
   fprintf(stderr, "Add scalar to REAL8Vector: max error = %g, max relative error = %g\n", maxdoublediff_addscalar, maxdoublerelerr_addscalar);
   fprintf(stderr, "Scale REAL4Vector: max error = %g, max relative error = %g\n", maxfloatdiff_scale, maxfloatrelerr_scale);
   fprintf(stderr, "Scale REAL8Vector: max error = %g, max relative error = %g\n", maxdoublediff_scale, maxdoublerelerr_scale);
   fprintf(stderr, "Sum vectors of vector sequence into vector sequence: max error = %g, max relative error = %g\n", maxfloatdiff_seqsum, maxfloatrelerr_seqsum);
   fprintf(stderr, "Subtract vector from vector sequence: max error = %g, max relative error = %g\n", maxfloatdiff_seqsub, maxfloatrelerr_seqsub);
#endif

   XLALDestroyREAL4Vector(floatvalues1);
   XLALDestroyREAL8Vector(doublevalues1);
   XLALDestroyREAL4Vector(floatvalues2);
   XLALDestroyREAL8Vector(doublevalues2);
   XLALDestroyREAL4Vector(floatresult_exp);
   XLALDestroyREAL8Vector(doubleresult_exp);
   XLALDestroyREAL4Vector(floatresult_vecsum);
   XLALDestroyREAL4Vector(floatresult_vecmult);
   XLALDestroyREAL4Vector(floatresult_addscalar);
   XLALDestroyREAL4Vector(floatresult_scale);
   XLALDestroyREAL8Vector(doubleresult_addscalar);
   XLALDestroyREAL8Vector(doubleresult_scale);
   XLALDestroyREAL4Vector(floatresult_sin2pix);
   XLALDestroyREAL4Vector(floatresult_cos2pix);
   XLALDestroyREAL8Vector(doubleresult_sin2pix);
   XLALDestroyREAL8Vector(doubleresult_cos2pix);

   return 0;

}




