/*
 *  Copyright (C) 2011 Evan Goetz
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

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>

REAL4Vector * sampleREAL4Vector(REAL4Vector *input, INT4 sampleSize, gsl_rng *rng);
REAL4Vector * sampleREAL4VectorSequence(REAL4VectorSequence *input, INT4 numberofvectors, INT4 sampleSize, gsl_rng *rng);
REAL4Vector * sampleREAL4VectorSequence_nozerosaccepted(REAL4VectorSequence *input, INT4 numberofvectors, INT4 sampleSize, gsl_rng *rng);

REAL8 calcMeanD(REAL8Vector *vector);
REAL8 calcStddevD(REAL8Vector *vector);
REAL8 expRandNum(REAL8 mu, gsl_rng *ptrToGenerator);
REAL8 ncx2cdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2cdf_withouttinyprob(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2cdf_withouttinyprob_withmatlabchi2cdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 ncx2pdf(REAL8 x, REAL8 dof, REAL8 delta);
REAL8 binodeviance(REAL8 x, REAL8 np);
REAL8 epsval(REAL8 val);
REAL8 ncx2inv(REAL8 p, REAL8 dof, REAL8 delta);
REAL4 ncx2inv_float(REAL8 p, REAL8 dof, REAL8 delta);
REAL8 norminv(REAL8 p, REAL8 mu, REAL8 sigma);
INT4 ks_test_exp(REAL8 *ksvalue, REAL4Vector *vector);
INT4 kuipers_test_exp(REAL8 *kuipervalue, REAL4Vector *vector);
REAL8 twospect_cdf_chisq_P(REAL8 x, REAL8 nu);
REAL8 matlab_cdf_chisq_P(REAL8 x, REAL8 nu);
REAL8 unitGaussianSNR(REAL8 value, REAL8 dof);

REAL4 ncx2cdf_float(REAL4 x, REAL4 dof, REAL4 delta);
REAL4 ncx2cdf_float_withouttinyprob(REAL4 x, REAL4 dof, REAL4 delta);
REAL4 ncx2cdf_float_withouttinyprob_withmatlabchi2cdf(REAL4 x, REAL4 dof, REAL4 delta);
REAL4 epsval_float(REAL4 val);
REAL4 calcMean(REAL4Vector *vector);
REAL4 calcMean_ignoreZeros(REAL4Vector *vector);
INT4 calcHarmonicMean(REAL4 *harmonicMean, REAL4Vector *vector, INT4 numfbins, INT4 numffts);
INT4 calcStddev(REAL4 *sigma, REAL4Vector *vector);
REAL4 calcStddev_ignoreZeros(REAL4Vector *vector);
INT4 calcRms(REAL4 *rms, REAL4Vector *vector);
INT4 calcMedian(REAL4 *median, REAL4Vector *vector);

INT4 sort_float_smallest(REAL4Vector *output, REAL4Vector *input);
void sort_double_ascend(REAL8Vector *vector);
void sort_float_ascend(REAL4Vector *vector);

void sumseries(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown);
void sumseries_eg(REAL8 *computedprob, REAL8 P, REAL8 C, REAL8 E, INT8 counter, REAL8 x, REAL8 dof, REAL8 halfdelta, REAL8 err, INT4 countdown);
void min_max_index_INT4Vector(INT4Vector *inputvector, INT4 *min_index_out, INT4 *max_index_out);

INT4 max_index(REAL4Vector *vector);
INT4 max_index_double(REAL8Vector *vector);
INT4 max_index_in_range(REAL4Vector *vector, INT4 startlocation, INT4 lastlocation);
INT4 max_index_from_vector_in_REAL4VectorSequence(REAL4VectorSequence *vectorsequence, INT4 vectornum);

INT4 qsort_REAL4_compar(const void *a, const void *b);
INT4 qsort_REAL8_compar(const void *a, const void *b);

#endif
