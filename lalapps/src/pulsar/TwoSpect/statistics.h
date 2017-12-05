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

INT4 sampleREAL4VectorAligned(REAL4VectorAligned *output, const REAL4VectorAligned *input, const gsl_rng *rng);
REAL4VectorAligned * sampleREAL4VectorAlignedArray_nozerosaccepted(const REAL4VectorAlignedArray *input, const UINT4 numberofvectors, const UINT4 sampleSize, const gsl_rng *rng);

REAL8 calcMeanD(const REAL8Vector *vector);
REAL8 calcStddevD(const REAL8Vector *vector);
REAL8 expRandNum(const REAL8 mu, const gsl_rng *ptrToGenerator);
REAL4VectorAligned * expRandNumVector(const UINT4 length, const REAL8 mu, const gsl_rng *ptrToGenerator);

INT4 ks_test_exp(REAL8 *ksvalue, const REAL4VectorAligned *vector);
INT4 kuipers_test_exp(REAL8 *kuipervalue, const REAL4VectorAligned *vector);

REAL4 calcMean(const REAL4VectorAligned *vector);
REAL4 calcMean_ignoreZeros(const REAL4VectorAligned *vector);
INT4 calcHarmonicMean(REAL4 *harmonicMean, const REAL4VectorAligned *vector, UINT4 numfbins, UINT4 numffts);
INT4 calcStddev(REAL4 *sigma, const REAL4VectorAligned *vector);
INT4 calcStddev_ignoreZeros(REAL4 *sigma, const REAL4VectorAligned *vector);
INT4 calcRms(REAL4 *rms, const REAL4VectorAligned *vector);
INT4 calcMedian(REAL4 *median, const REAL4VectorAligned *vector);

INT4 sort_float_smallest(REAL4VectorAligned *output, const REAL4VectorAligned *input);
void sort_double_ascend(REAL8Vector *vector);
void sort_float_ascend(REAL4VectorAligned *vector);

UINT4 max_index(const REAL4VectorAligned *vector);
UINT4 max_index_double(const REAL8Vector *vector);
UINT4 max_index_in_range(const REAL4VectorAligned *vector, const UINT4 startlocation, const UINT4 lastlocation);
INT4 min_max_index_INT4Vector(const INT4Vector *inputvector, UINT4 *min_index_out, UINT4 *max_index_out);

INT4 qsort_REAL4_compar(const void *a, const void *b);
INT4 qsort_REAL8_compar(const void *a, const void *b);

#endif
