/*
 *  LALInferenceKDE.c:  Bayesian Followup, kernel density estimator.
 *
 *  Copyright (C) 2013 Ben Farr
 *
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
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceKDE.h>


/**
 * Allocate, fill, and tune a Gaussian kernel density estimate given an array of points.
 *
 * Given a set of points, the distribution from which the points were sampled from is
 * estimated by a kernel density estimate (KDE), using a Gaussian kernal.  This routine
 * initializes a new KDE and sets the bandwidth according to Scott's rule.  From this
 * it is possible to both evaluate the estimate of the probability density function at
 * an arbitrary location, and generate more samples from the estimated distribution. A
 * mask can be provided to select a subset of the provided data, otherwise all data is
 * used.
 * @param[in] pts  A REAL8 array containing the points to estimate the distribution of.
 * @param[in] npts The total number of points contained in \a pts.
 * @param[in] dim  The number of dimensions of the points contained in \a pts.
 * @param[in] mask An optional \a npts long 0/1 array to mask points in \a pts.
 * @return A LALInferenceKDE structure containing the KDE of the distribution.
 */
LALInferenceKDE *LALInferenceNewKDE(REAL8 *pts, UINT4 npts, UINT4 dim, UINT4 *mask) {
    UINT4 i,j;

    /* Determine the total number of samples */
    UINT4 count = 0;
    if (mask != NULL) {
        for (i = 0; i < npts; i++) {
            if (mask[i])
                count++;
        }
    } else {
        count = npts;
    }

    /* Initialize KDE */
    LALInferenceKDE *kde = LALInferenceInitKDE(count, dim);

    /* Fill in samples */
    UINT4 idx = 0;
    for (i = 0; i < npts; i++) {
        if (mask == NULL || mask[i]) {
            for (j = 0; j < dim; j++)
                gsl_matrix_set(kde->data, idx, j, pts[i*dim + j]);
            idx++;
        }
    }

    /* Set bandwidth */
    LALInferenceSetKDEBandwidth(kde);

    return kde;
}


/**
 * Allocate, fill, and tune a Gaussian kernel density estimate given a matrix of points.
 *
 * Estimate the underlying distribution of points in a GSL matrix using a Gaussian
 * KDE.
 * @param[in] data A GSL matrix with rows containing samples from a target distribution.
 * @param[in] mask An optional \a npts long 0/1 array to mask points in \a pts.
 * @return A LALInferenceKDE structure containing the KDE of the distribution.
 * \sa LALInferenceNewKDE()
 */
LALInferenceKDE *LALInferenceNewKDEfromMat(gsl_matrix *data, UINT4 *mask) {
    UINT4 i,j;
    UINT4 npts = data->size1;
    UINT4 dim = data->size2;

    /* Determine the total number of samples */
    UINT4 count = 0;
    if (mask != NULL) {
        for (i = 0; i < npts; i++) {
            if (mask[i])
                count++;
        }
    } else {
        count = npts;
    }

    /* Initialize KDE */
    LALInferenceKDE *kde = LALInferenceInitKDE(count, dim);

    /* Fill in samples */
    UINT4 idx = 0;
    for (i = 0; i < npts; i++) {
        if (mask == NULL || mask[i]) {
            for (j = 0; j < dim; j++)
                gsl_matrix_set(kde->data, idx, j, gsl_matrix_get(data, i, j));
            idx++;
        }
    }

    /* Calculate covariances and set bandwidth */
    LALInferenceSetKDEBandwidth(kde);

    return kde;
}


/**
 * Construct an empty KDE structure.
 *
 * Create an empty LALInferenceKDE structure, allocated to handle the given size and dimension.
 * @param[in] npts The number of samples that will be used to estimate the distribution.
 * @param[in] dim  The number of dimensions that the probability density function will be estimated in.
 * @return An allocated, empty LALInferenceKDE structure.
 * \sa LALInferenceKDE, LALInferenceSetKDEBandwidth()
 */
LALInferenceKDE *LALInferenceInitKDE(UINT4 npts, UINT4 dim) {
    LALInferenceKDE *kde = XLALMalloc(sizeof(LALInferenceKDE));
    kde->dim = dim;
    kde->npts = npts;
    kde->mean = gsl_vector_calloc(dim);
    kde->cholesky_decomp_cov = gsl_matrix_calloc(dim, dim);
    kde->cholesky_decomp_cov_lower = gsl_matrix_calloc(dim, dim);

    kde->data = gsl_matrix_alloc(npts, dim);

    return kde;
}


/**
 * Free an allocated KDE structure.
 *
 * Frees all memory allocated for a given KDE structure.
 * @param[in] KDE The KDE structure to be freed.
 * \sa LALInferenceKDE, LALInferenceInitKDE
 */
void LALInferenceDestroyKDE(LALInferenceKDE *kde) {
    if (kde->mean) gsl_vector_free(kde->mean);
    gsl_matrix_free(kde->cholesky_decomp_cov);
    gsl_matrix_free(kde->cholesky_decomp_cov_lower);
    gsl_matrix_free(kde->data);

    if (kde->cov)
        gsl_matrix_free(kde->cov);

    XLALFree(kde);

    return;
}


/**
 * Compute the Cholesky decomposition of a matrix.
 *
 * A wrapper for gsl_linalg_cholesky_decomp() that avoids halting if the matrix
 * is found not to be positive-definite.  This often happens when decomposing a
 * covariance matrix that is poorly estimated due to low sample size.
 * @param matrix The matrix to decompose (in place).
 * @return Status of call to gsl_linalg_cholesky_decomp().
 */
INT4 LALInferenceCholeskyDecompose(gsl_matrix *mat) {
    INT4 status;

    /* Turn off default GSL error handling (i.e. aborting), and catch
     * errors decomposing due to non-positive definite covariance matrices */
    gsl_error_handler_t *default_gsl_error_handler = gsl_set_error_handler_off();

    status = gsl_linalg_cholesky_decomp(mat);
    if (status) {
        if (status != GSL_EDOM) {
            fprintf(stderr, "ERROR: Unexpected problem Cholesky-decomposing matrix.\n");
            exit(-1);
        }
    }

    /* Return to default GSL error handling */
    gsl_set_error_handler(default_gsl_error_handler);
    return status;
}


/**
 * Calculate the bandwidth and normalization factor for a KDE.
 *
 * Use Scott's rule to determine the bandwidth, and corresponding normalization factor, for a KDE.
 * @param[in] KDE The kernel density estimate to estimate the bandwidth of.
 */
void LALInferenceSetKDEBandwidth(LALInferenceKDE *kde) {
    /* Use Scott's Bandwidth method */
    REAL8 det_cov, cov_factor;
    UINT4 i, j;
    INT4 status;

    /* Calculate average and coveriance */
    kde->mean = LALInferenceComputeMean(kde->data);
    kde->cov = LALInferenceComputeCovariance(kde->data);

    cov_factor = pow((REAL8)kde->npts, -1./(REAL8)(kde->dim + 4));
    gsl_matrix_scale(kde->cov, cov_factor*cov_factor);

    gsl_matrix_memcpy(kde->cholesky_decomp_cov, kde->cov);
    status = LALInferenceCholeskyDecompose(kde->cholesky_decomp_cov);

    /* If cholesky decomposition failed, set the normalization to infinity */
    if (status) {
        fprintf(stderr, "Non-positive-definite matrix encountered when setting KDE bandwidth.");
        fprintf(stderr, " (number of samples = %i).\n", (INT4) kde->data->size1);
        kde->log_norm_factor = INFINITY;
        return;
    }

    /* Zero out upper right triangle of decomposed covariance matrix,
     * which contains the transpose */
    gsl_matrix_memcpy(kde->cholesky_decomp_cov_lower, kde->cholesky_decomp_cov);
    for (i = 0 ; i < kde->dim; i++) {
        for (j = i+1 ; j < kde->dim; j++)
            gsl_matrix_set(kde->cholesky_decomp_cov_lower, i, j, 0.);
    }

    det_cov = LALInferenceMatrixDet(kde->cov);
    kde->log_norm_factor = log(kde->npts * sqrt(pow(2*LAL_PI, kde->dim) * det_cov));
}


/**
 * Evaluate the (log) PDF from a KDE at a single point.
 *
 * Calculate the (log) value of the probability density function estimate from
 * a kernel density estimate at a single point.
 * @param[in] KDE   The kernel density estimate to evaluate.
 * @param[in] point An array containing the point to evaluate the PDF at.
 * @return The value of the estimated probability density function at \a point.
 */
REAL8 LALInferenceKDEEvaluatePoint(LALInferenceKDE *kde, REAL8 *point) {
    UINT4 dim = kde->dim;
    UINT4 npts = kde->npts;
    UINT4 i, j;

    /* If the normalization is inifinite, don't bother calculating anything */
    if (isinf(kde->log_norm_factor))
        return -INFINITY;

    gsl_vector_view x = gsl_vector_view_array(point, dim);
    gsl_vector_view d;

    /* Vectors that will hold the difference and transformed distance */
    gsl_vector *diff = gsl_vector_alloc(dim);
    gsl_vector *tdiff = gsl_vector_alloc(dim);

    /* Loop over points in KDE dataset, using the Cholesky decomposition
     * of the covariance to avoid ever inverting the covariance matrix */
    REAL8 energy;
    REAL8* results = XLALMalloc(npts * sizeof(REAL8));
    for (i=0; i<npts; i++) {
        d = gsl_matrix_row(kde->data, i);
        gsl_vector_memcpy(diff, &d.vector);

        gsl_vector_sub(diff, &x.vector);
        gsl_linalg_cholesky_solve(kde->cholesky_decomp_cov, diff, tdiff);
        gsl_vector_mul(diff, tdiff);

        energy = 0.;
        for (j=0; j<dim; j++)
            energy += gsl_vector_get(diff, j);
        results[i] = -energy/2.;
    }

    /* Normalize the result and return */
    REAL8 result = log_add_exps(results, npts) - kde->log_norm_factor;

    gsl_vector_free(diff);
    gsl_vector_free(tdiff);
    XLALFree(results);

    return result;
}


/**
 * Draw a sample from a kernel density estimate.
 *
 * Draw a sample from a distribution, as estimated by a kernel density estimator.
 * @param[in] KDE The kernel density estimate to draw \a point from.
 * @param[in] rng A GSL random number generator to use for random number generation.
 * @return The sample drawn from the distribution.
 */
REAL8 *LALInferenceDrawKDESample(LALInferenceKDE *kde, gsl_rng *rng) {
    UINT4 dim = kde->dim;
    UINT4 j;
    REAL8 val;

    REAL8 *point = XLALMalloc(dim * sizeof(REAL8));
    gsl_vector_view pt = gsl_vector_view_array(point, dim);

    /* Draw a random sample from KDE dataset */
    UINT4 ind = gsl_rng_uniform_int(rng, kde->npts);
    gsl_vector_view d = gsl_matrix_row(kde->data, ind);
    gsl_vector_memcpy(&pt.vector, &d.vector); 

    /* Draw individual parameters from 1D unit Gaussians */
    gsl_vector *unit_draw = gsl_vector_alloc(dim);
    for (j = 0; j < dim; j++) {
        val = gsl_ran_ugaussian(rng);
        gsl_vector_set(unit_draw, j, val);
    }

    /* Scale and shift the uncorrelated unit-width sample */
    gsl_blas_dgemv(CblasNoTrans, 1.0, kde->cholesky_decomp_cov_lower, unit_draw, 1.0, &pt.vector);

    gsl_vector_free(unit_draw);
    return point;
}


/**
 * Calculate the determinant of a matrix.
 *
 * Calculates the determinant of a GSL matrix.
 * @param[in] mat The matrix to calculate the determinant of.
 * @return The determinant of \a mat.
 */
REAL8 LALInferenceMatrixDet(gsl_matrix *mat) {
    REAL8 det=0.0;
    INT4 sign=0;
    INT4 * signum = &sign;
    UINT4 dim = mat->size1;

    gsl_permutation * p = gsl_permutation_calloc(dim);
    gsl_matrix * tmp_ptr = gsl_matrix_calloc(dim, dim);
    gsl_matrix_memcpy(tmp_ptr, mat);

    gsl_linalg_LU_decomp(tmp_ptr, p, signum);
    det = gsl_linalg_LU_det(tmp_ptr, *signum);

    gsl_permutation_free(p);
    gsl_matrix_free(tmp_ptr);

    return det;
}


/**
 * Calculate the mean of a data set.
 *
 * Calculate the mean vector of a data set contained in the rows of a
 * GSL matrix.
 * @param[in] data GSL matrix data set to compute the mean of.
 * @return A vector containing the mean of \a data.
 */
gsl_vector *LALInferenceComputeMean(gsl_matrix *data) {
    UINT4 i;
    UINT4 npts = data->size1;
    UINT4 dim = data->size2;

    gsl_vector *mean = gsl_vector_calloc(dim);
    gsl_vector_view x;

    for (i = 0; i < npts; i++) {
        x = gsl_matrix_row(data, i);
        gsl_vector_add(mean, &x.vector);
    }
    gsl_vector_scale(mean, 1./(REAL8)npts);

    return mean;
}


/* Compute the (GSL) covariance matrix of a list of samples */
/**
 * Calculate the covariance matrix of a data set.
 *
 * Calculate the covariance matrix of a data set contained in the rows of a
 * GSL matrix.
 * @param[in] data GSL matrix data set to compute the covariance matrix of.
 * @return The covariance matrix of \a data.
 */
gsl_matrix *LALInferenceComputeCovariance(gsl_matrix *data) {
    UINT4 i, j, k;
    REAL8 var;

    UINT4 npts = data->size1;
    UINT4 dim = data->size2;

    gsl_vector *adiff = gsl_vector_alloc(npts);
    gsl_vector *bdiff = gsl_vector_alloc(npts);
    gsl_matrix *cov = gsl_matrix_alloc(dim, dim);

    gsl_vector *mean = LALInferenceComputeMean(data);
    for (i = 0; i < dim; i++) {
        gsl_vector_view a = gsl_matrix_column(data, i);
        gsl_vector_memcpy(adiff, &a.vector);
        gsl_vector_add_constant(adiff, -gsl_vector_get(mean, i));

        for (j = 0; j < dim; j++) {
            gsl_vector_view b = gsl_matrix_column(data, j);
            gsl_vector_memcpy(bdiff, &b.vector);
            gsl_vector_add_constant(bdiff, -gsl_vector_get(mean, j));

            gsl_vector_mul(bdiff, adiff);
            var = 0.;
            for (k = 0; k < npts; k++)
                var += gsl_vector_get(bdiff, k);
            var /= (REAL8)(npts-1);

            gsl_matrix_set(cov, i, j, var);
        }
    }

    gsl_vector_free(adiff);
    gsl_vector_free(bdiff);
    gsl_vector_free(mean);
    return cov;
}


/**
 * Determine the log of the sum of an array of exponentials.
 *
 * Utility for calculating the log of the sum of an array of exponentials.  Useful
 * for avoiding overflows.
 * @param[in] vals Array containing the values to be exponentiated, summed, and logged.
 * @param[in] size Number of elements in \a vals.
 * @return The log of the sum of elements in \a vals.
 */
REAL8 log_add_exps(REAL8 *vals, UINT4 size) {
    UINT4 i;

    REAL8 max_comp = 0.;
    for (i = 0; i < size; i++) {
        if (max_comp < vals[i])
            max_comp = vals[i];
    }

    REAL8 result = 0.;
    for (i = 0; i < size; i++)
        result += exp(vals[i]-max_comp);
    result = max_comp + log(result);

    return result;
}
