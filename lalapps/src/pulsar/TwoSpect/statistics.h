/*
 *  Copyright (C) 2011, 2015 Evan Goetz
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

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include <lal/AVFactories.h>
#include "TwoSpectTypes.h"

INT4 sampleREAL4Vector(REAL4Vector *output, REAL4Vector *input, gsl_rng *rng);
REAL4Vector * sampleAlignedREAL4VectorArray_nozerosaccepted(alignedREAL4VectorArray *input, INT4 numberofvectors, INT4 sampleSize, gsl_rng *rng);

REAL8 calcMeanD(REAL8Vector *vector);
REAL8 calcStddevD(REAL8Vector *vector);
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator);

INT4 ks_test_exp(REAL8 *ksvalue, REAL4Vector *vector);
INT4 kuipers_test_exp(REAL8 *kuipervalue, REAL4Vector *vector);

REAL4 calcMean(REAL4Vector *vector);
REAL4 calcMean_ignoreZeros(REAL4Vector *vector);
INT4 calcHarmonicMean(REAL4 *harmonicMean, REAL4Vector *vector, INT4 numfbins, INT4 numffts);
INT4 calcStddev(REAL4 *sigma, REAL4Vector *vector);
INT4 calcStddev_ignoreZeros(REAL4 *sigma, REAL4Vector *vector);
INT4 calcRms(REAL4 *rms, REAL4Vector *vector);
INT4 calcMedian(REAL4 *median, REAL4Vector *vector);

INT4 sort_float_smallest(REAL4Vector *output, REAL4Vector *input);
void sort_double_ascend(REAL8Vector *vector);
void sort_float_ascend(REAL4Vector *vector);

INT4 max_index(REAL4Vector *vector);
INT4 max_index_double(REAL8Vector *vector);
INT4 max_index_in_range(REAL4Vector *vector, INT4 startlocation, INT4 lastlocation);
INT4 min_max_index_INT4Vector(INT4Vector *inputvector, INT4 *min_index_out, INT4 *max_index_out);
INT4 max_index_from_vector_in_REAL4VectorSequence(REAL4VectorSequence *vectorsequence, INT4 vectornum);

INT4 qsort_REAL4_compar(const void *a, const void *b);
INT4 qsort_REAL8_compar(const void *a, const void *b);

#endif
