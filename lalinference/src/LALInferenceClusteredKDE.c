/*
 *  LALInferenceClusteringKDE.c:  Bayesian Followup, clustering-optimized Gaussian kernel density estimator.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>

#include <lal/LALInference.h>
#include <lal/LALInferenceKDE.h>
#include <lal/LALInferenceClusteredKDE.h>


/**
 * Kmeans cluster data, increasing k until the Bayes Information Criteria stop increasing.
 *
 * Run the kmeans clustering algorithm incrementally, calculating the Bayes Information Criteria (BIC)
 * at each step, stopping when the BIC stops increasing.  The BIC is calculated by estimating the
 * distribution in each cluster with a Gaussian kernel density estimate, then weighting that clusters
 * contribution to the likelihood by the fraction of total samples in that cluster.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceIncrementalKmeans(gsl_matrix *data, UINT4 ntrials, gsl_rng *rng) {
    UINT4 k = 1;
    REAL8 best_bic = -INFINITY;
    REAL8 bic;

    LALInferenceKmeans *kmeans;
    LALInferenceKmeans *best_clustering = NULL;

    /* Start at k=1 and increase k until the BIC stops increasing. */
    while (1) {
        kmeans = LALInferenceKmeansRunBestOf(k, data, ntrials, rng);

        /* If kmeans creation failed, the cluster size was likely too small */
        if (kmeans)
            bic = LALInferenceKmeansBIC(kmeans);
        else
            bic = -INFINITY;

        /* Keep going until the BIC decreases */
        if (bic > best_bic) {
            if (best_clustering)
                LALInferenceKmeansDestroy(best_clustering);
            best_clustering = kmeans;
            best_bic = bic;
        } else {
            LALInferenceKmeansDestroy(kmeans);
            break;
        }
        k++;
    }

    return best_clustering;
}


/**
 * Kmeans cluster data, find k that maximizes BIC assuming BIC(k) is concave-down.
 *
 * Run the kmeans clustering algorithm, searching for the k that maximizes the
 *   Bayes Information Criteria (BIC).  The BIC is calculated by estimating the
 *   distribution in each cluster with a Gaussian kernel density estimate, then
 *   weighing that clusters contribution to the likelihood by the fraction of
 *   total samples in that cluster.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceOptimizedKmeans(gsl_matrix *data, UINT4 ntrials, gsl_rng *rng) {
    UINT4 k, low_k = 1, mid_k = 2, high_k = 4;
    REAL8 bic, low_bic, mid_bic, high_bic;
    LALInferenceKmeans *low_kmeans = NULL;
    LALInferenceKmeans *mid_kmeans = NULL;
    LALInferenceKmeans *high_kmeans = NULL;

    /* Calculate starting clusters and BICs */
    low_kmeans = LALInferenceKmeansRunBestOf(low_k, data, ntrials, rng);
    if (!low_kmeans) {
        fprintf(stderr, "Can't build basic KDE.  Perhaps sample size (%i) is too small.\n",
                (INT4)data->size1);
        if (mid_kmeans) LALInferenceKmeansDestroy(mid_kmeans);
        if (high_kmeans) LALInferenceKmeansDestroy(high_kmeans);
        return NULL;
    }
    low_bic = LALInferenceKmeansBIC(low_kmeans);

    mid_kmeans = LALInferenceKmeansRunBestOf(mid_k, data, ntrials, rng);
    if (!mid_kmeans) {
        printf("Can't cluster, even with k=2.  Falling back to standard KDE.\n");
        return low_kmeans;
    }
    mid_bic = LALInferenceKmeansBIC(mid_kmeans);

    /* If high-k fails, try shrinking before giving up */
    high_kmeans = LALInferenceKmeansRunBestOf(high_k, data, ntrials, rng);
    if (!high_kmeans) {
        high_k--;
        high_kmeans = LALInferenceKmeansRunBestOf(high_k, data, ntrials, rng);
        if (!high_kmeans) {
            if (mid_bic > low_bic) {
                LALInferenceKmeansDestroy(low_kmeans);
                return mid_kmeans;
            } else {
                LALInferenceKmeansDestroy(mid_kmeans);
                return low_kmeans;
            }
        }
    }
    high_bic = LALInferenceKmeansBIC(high_kmeans);

    /* Keep doubling the highest sample until the peak is passed */
    while (high_bic > mid_bic) {
        low_k = mid_k;
        mid_k = high_k;
        high_k *= 2;

        low_bic = mid_bic;
        mid_bic = high_bic;

        LALInferenceKmeansDestroy(low_kmeans);
        low_kmeans = mid_kmeans;
        mid_kmeans = high_kmeans;

        while (1) {
            high_kmeans = LALInferenceKmeansRunBestOf(high_k, data, ntrials, rng);
            if (!high_kmeans) {
                high_k = mid_k + (high_k - mid_k)/2;
                if (high_k < mid_k + 1)
                    high_bic = -INFINITY;
                else
                    continue;
            }
            high_bic = LALInferenceKmeansBIC(high_kmeans);
            break;
        }
    }

    /* Perform bisection search to find maximum */
    LALInferenceKmeans *kmeans = NULL;
    while (high_k - low_k > 2) {
        if (high_k - mid_k > mid_k - low_k) {
            k = mid_k + (high_k - mid_k)/2;
            if (kmeans)
                LALInferenceKmeansDestroy(kmeans);
            kmeans = LALInferenceKmeansRunBestOf(k, data, ntrials, rng);
            bic = LALInferenceKmeansBIC(kmeans);

            if (bic > mid_bic) {
                low_k = mid_k;
                mid_k = k;
                low_bic = mid_bic;
                mid_bic = bic;
                LALInferenceKmeansDestroy(low_kmeans);
                low_kmeans = mid_kmeans;
                mid_kmeans = kmeans;
                kmeans = NULL;
            } else {
                high_k = k;
                high_bic = bic;
                LALInferenceKmeansDestroy(high_kmeans);
                high_kmeans = kmeans;
                kmeans = NULL;
            }
        } else {
            k = low_k + (mid_k - low_k)/2;
            if (kmeans)
                LALInferenceKmeansDestroy(kmeans);
            kmeans = LALInferenceKmeansRunBestOf(k, data, ntrials, rng);
            bic = LALInferenceKmeansBIC(kmeans);

            if (bic > mid_bic) {
                high_k = mid_k;
                mid_k = k;
                high_bic = mid_bic;
                mid_bic = bic;

                LALInferenceKmeansDestroy(high_kmeans);
                high_kmeans = mid_kmeans;
                mid_kmeans = kmeans;
                kmeans = NULL;
            } else {
                low_k = k;
                low_bic = bic;
                LALInferenceKmeansDestroy(low_kmeans);
                low_kmeans = kmeans;
                kmeans = NULL;
            }
        }
    }

    /* Cleanup and return the best clustering */
    if (low_kmeans != mid_kmeans)
        LALInferenceKmeansDestroy(low_kmeans);
    if (high_kmeans != mid_kmeans)
        LALInferenceKmeansDestroy(high_kmeans);
    if (kmeans && kmeans != mid_kmeans)
        LALInferenceKmeansDestroy(kmeans);
    return mid_kmeans;
}


/**
 * Xmeans cluster data, splitting centroids in kmeans according to the Bayes Information Criteria.
 *
 * Run an xmeans-style clustering, increasing a kmeans clustering starting from k=1 by splitting each
 * centroid individually and checking for an increase of the BIC over the parent cluster.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceXmeans(gsl_matrix *data, UINT4 ntrials, gsl_rng *rng) {
    UINT4 kmax = 50;
    UINT4 split_size = 2;
    UINT4 c,k;
    REAL8 starting_bic, ending_bic;
    REAL8 old_bic, new_bic;

    LALInferenceKmeans *sub_kmeans = NULL;
    LALInferenceKmeans *new_kmeans = NULL;

    LALInferenceKmeans *kmeans = LALInferenceCreateKmeans(2, data, rng);

    /* If kmeans wasn't setup, return nothing */
    if (!kmeans)
        return NULL;

    LALInferenceKmeansForgyInitialize(kmeans);
    LALInferenceKmeansRun(kmeans);

    /* Keep increasing k by attempting to split centroids.  Stop if the BIC is maximized */
    while (kmeans->k < kmax) {
        if (kmeans->recursive_centroids != NULL)
            gsl_matrix_free(kmeans->recursive_centroids);

        /* Store the BIC of the initial clustering */
        starting_bic = LALInferenceKmeansBIC(kmeans);

        /* For each existing centroid, attempt to split and see if the BIC increases locally */
        for (c = 0; c < kmeans->k; c++) {
            sub_kmeans = LALInferenceKmeansExtractCluster(kmeans, c);

            /* If the cluster wasn't extracted successfully, reject the split.  */
            if (sub_kmeans) {
                old_bic = LALInferenceKmeansBIC(sub_kmeans);

                new_kmeans = LALInferenceKmeansRunBestOf(split_size, sub_kmeans->data, ntrials, rng);
                new_bic = LALInferenceKmeansBIC(new_kmeans);
            } else {
                old_bic = INFINITY;
                new_bic = -INFINITY;
            }

            if (new_bic > old_bic) {
                /* BIC increased, so keep the new centroids */
                for (k = 0; k < new_kmeans->k; k++) {
                    gsl_vector_view centroid = gsl_matrix_row(new_kmeans->centroids, k);
                    accumulate_vectors(&kmeans->recursive_centroids, &centroid.vector);
                }
            } else {
                /* BIC decreased, keep the parent centroid */
                gsl_vector_view centroid = gsl_matrix_row(kmeans->centroids, c);
                accumulate_vectors(&kmeans->recursive_centroids, &centroid.vector);
            }
            LALInferenceKmeansDestroy(sub_kmeans);
            LALInferenceKmeansDestroy(new_kmeans);
        }

        /* Create a self-consistent kmeans from the resulting centroids */
        gsl_matrix *centroids = kmeans->recursive_centroids;
        k = centroids->size1;

        new_kmeans = LALInferenceCreateKmeans(k, data, rng);
        if (new_kmeans) {
            gsl_matrix_memcpy(new_kmeans->centroids, centroids);
            LALInferenceKmeansRun(new_kmeans);

            /* Store the BIC after all centroids have been attempted to be split */
            ending_bic = LALInferenceKmeansBIC(new_kmeans);
        } else {
            /* If new kmeans couldn't be created, reject it */
            ending_bic = -INFINITY;
        }

        if (ending_bic > starting_bic) {
            /* BIC improved, so update and continue */
            LALInferenceKmeansDestroy(kmeans);
            kmeans = new_kmeans;
        } else {
            /* Clustering was unchanged, or resulted in a lower BIC.  End */
            LALInferenceKmeansDestroy(new_kmeans);
            break;
        }
    }

    return kmeans;
}


/**
 * Run kmeans recursively, splitting each centroid recursively, accepting if the BIC increases.
 *
 * Run kmeans recursively, splitting each centroid in two recursively.  At each step,
 * the BIC is checked before and after splitting.  If the split increases the BIC, each
 * centroid is then attempted to be split.  If the BIC decreases, the split is rejected
 * and the centroid is kept.  This is less precise than other methods, but much faster,
 * as the likelihood is computed at few and fewer points as centroids are split.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceRecursiveKmeans(gsl_matrix *data, UINT4 ntrials, gsl_rng *rng) {
    UINT4 k;

    /* Perform recursive splitting.  Return NULL if kmeans creation fails */
    LALInferenceKmeans *split_kmeans = LALInferenceCreateKmeans(1, data, rng);
    if (!split_kmeans)
        return NULL;

    LALInferenceKmeansForgyInitialize(split_kmeans);
    LALInferenceKmeansRecursiveSplit(split_kmeans, ntrials, rng);

    /* Take the final centroids and make a fully self-consistent kmeans */
    gsl_matrix *centroids = split_kmeans->recursive_centroids;
    LALInferenceKmeansDestroy(split_kmeans);

    k = centroids->size1;
    LALInferenceKmeans *kmeans = LALInferenceCreateKmeans(k, data, rng);

    if (!kmeans)
        return NULL;

    gsl_matrix_memcpy(kmeans->centroids, centroids);
    LALInferenceKmeansRun(kmeans);

    return kmeans;
}


/**
 * Recursively split a k=1 kmeans.
 *
 * Split a centroid in 2 recursively, evaluating the BIC to determine if the
 * split is kept.  This is a cheap method, as the BIC becomes cheaper to
 * compute after each split.
 * @param[in] kmeans A k=1 kmeans to be split recursively.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * \sa LALInferenceRecursiveKmeans()
 */
void LALInferenceKmeansRecursiveSplit(LALInferenceKmeans *kmeans, UINT4 ntrials, gsl_rng *rng) {
    LALInferenceKmeansRun(kmeans);
    UINT4 cluster, i;
    UINT4 split_size = 2;
    REAL8 current_bic, new_bic;

    LALInferenceKmeans *new_kmeans = LALInferenceKmeansRunBestOf(split_size, kmeans->data, ntrials, rng);

    current_bic = LALInferenceKmeansBIC(kmeans);
    new_bic = LALInferenceKmeansBIC(new_kmeans);

    /* If the BIC increased, attempt a further splitting of the centroids */
    if (new_bic > current_bic) {
        for (cluster = 0; cluster < new_kmeans->k; cluster++) {
            LALInferenceKmeans *sub_kmeans = LALInferenceKmeansExtractCluster(new_kmeans, cluster);
            if (sub_kmeans) LALInferenceKmeansRecursiveSplit(sub_kmeans, ntrials, rng);

            UINT4 n_new_centroids = sub_kmeans->recursive_centroids->size1;
            for (i = 0; i < n_new_centroids ; i++) {
                gsl_vector_view centroid = gsl_matrix_row(sub_kmeans->recursive_centroids, i);
                accumulate_vectors(&kmeans->recursive_centroids, &centroid.vector);
            }

            LALInferenceKmeansDestroy(sub_kmeans);
        }
    } else {
        /* The split was rejected, so save the old centroid */
        gsl_vector_view centroid = gsl_matrix_row(kmeans->centroids, 0);
        accumulate_vectors(&kmeans->recursive_centroids, &centroid.vector);
    }

    LALInferenceKmeansDestroy(new_kmeans);
}


/**
 * Run the kmeans algorithm until cluster assignments don't change.
 *
 * Starting with some random initialization, points are assigned to the closest centroids,
 * then the centroids are calculated of the new cluster.  This is repeated until the
 * assignments stop changing.
 * @param kmeans The initialized kmeans to run.
 */
void LALInferenceKmeansRun(LALInferenceKmeans *kmeans) {
    UINT4 i;

    while (kmeans->has_changed) {
        kmeans->has_changed = 0;

        LALInferenceKmeansAssignment(kmeans);
        LALInferenceKmeansUpdate(kmeans);
    }

    for (i = 0; i < kmeans->k; i++)
        kmeans->weights[i] = (REAL8)(kmeans->sizes[i]) / (REAL8)(kmeans->npts);
}


/**
 * Run a kmeans several times and return the best.
 *
 * The kmeans is run \a ntrials times, each time with a new random initialization.  The one that results
 * in the highest Bayes Information Criteria (BIC) is returned.
 * @param[in]  k       The number of clusters to use.
 * @param[in]  samples The (unwhitened) data to cluster.
 * @param[in]  ntrials The number of random initialization to run.
 * @param[in]  rng     The GSL random number generator to use for random initializations.
 * @return The kmeans with the highest BIC of \a ntrials attempts.
 */
LALInferenceKmeans *LALInferenceKmeansRunBestOf(UINT4 k, gsl_matrix *samples, UINT4 ntrials, gsl_rng *rng) {
    UINT4 i;
    REAL8 error = -INFINITY;

    LALInferenceKmeans *kmeans;
    LALInferenceKmeans *best_kmeans = NULL;
    REAL8 bic = 0;
    REAL8 best_bic = -INFINITY;
    for (i = 0; i < ntrials; i++) {
        kmeans = LALInferenceCreateKmeans(k, samples, rng);
        if (!kmeans)
            continue;

        LALInferenceKmeansForgyInitialize(kmeans);
        LALInferenceKmeansRun(kmeans);

        /* Assume BIC hasn't changed if error (summed-dist^2 from assigned centroids) hasn't */
        if (error != kmeans->error) {
            error = kmeans->error;
            bic = LALInferenceKmeansBIC(kmeans);
        }

        if (bic > best_bic) {
            if (best_kmeans)
                LALInferenceKmeansDestroy(best_kmeans);
            best_bic = bic;
            best_kmeans = kmeans;
        } else {
            LALInferenceKmeansDestroy(kmeans);
        }
    }

    return best_kmeans;
}

/**
 * Generate a new kmeans struct from a set of data.
 *
 * Whitens data and fills a newly allocated kmeans structure of the requested size.
 * @param[in] k    The requested number of clusters.
 * @param[in] data The unwhitened data to be clustered.
 * @param[in] rng  A GSL random number generator used for initialization.
 */
LALInferenceKmeans *LALInferenceCreateKmeans(UINT4 k, gsl_matrix *data, gsl_rng *rng) {
    INT4 status;

    /* Return nothing if given empty dataset */
    if (!data)
        return NULL;

    LALInferenceKmeans *kmeans = XLALMalloc(sizeof(LALInferenceKmeans));
    kmeans->k = k;
    kmeans->npts = data->size1;
    kmeans->dim = data->size2;
    kmeans->has_changed = 1;
    kmeans->rng = rng;

    /* Allocate GSL matrices and vectors */
    kmeans->mean = gsl_vector_alloc(kmeans->dim);
    kmeans->data = gsl_matrix_alloc(kmeans->npts, kmeans->dim);
    kmeans->chol_dec_cov = gsl_matrix_alloc(kmeans->dim, kmeans->dim);
    kmeans->unwhitened_chol_dec_cov = gsl_matrix_alloc(kmeans->dim, kmeans->dim);
    kmeans->centroids = gsl_matrix_alloc(kmeans->k, kmeans->dim);
    kmeans->recursive_centroids = NULL;

    /* Allocate everything else */
    kmeans->assignments = XLALMalloc(kmeans->npts * sizeof(UINT4));
    kmeans->sizes = XLALMalloc(kmeans->k * sizeof(UINT4));
    kmeans->mask = XLALMalloc(kmeans->npts * sizeof(UINT4));
    kmeans->weights = XLALMalloc(kmeans->k * sizeof(REAL8));
    kmeans->KDEs = NULL;

    /* Set distance and centroid calculators */
    kmeans->dist = &euclidean_dist_squared;
    kmeans->centroid = &euclidean_centroid;

    /* Store the unwhitened mean before whitening.  Needed for further transformation. */
    LALInferenceComputeMean(kmeans->mean, data);

    /* Whiten the data */
    kmeans->data = LALInferenceWhitenSamples(data);
    if (!kmeans->data) {
        fprintf(stderr, "Unable to whiten data.  No proposal has been built.\n");
        LALInferenceKmeansDestroy(kmeans);
        return NULL;
    }

    LALInferenceComputeCovariance(kmeans->unwhitened_chol_dec_cov, data);
    LALInferenceComputeCovariance(kmeans->chol_dec_cov, kmeans->data);

    /* Cholesky decompose the covariance matrix and decorrelate the data */
    status = LALInferenceCholeskyDecompose(kmeans->unwhitened_chol_dec_cov);
    if (status) {
        fprintf(stderr, "Unwhitened covariance matrix not positive-definite.  No proposal has been built.\n");
        return NULL;
    }

    status = LALInferenceCholeskyDecompose(kmeans->chol_dec_cov);
    if (status) {
        fprintf(stderr, "Whitened covariance matrix not positive-definite.  No proposal has been built.\n");
        return NULL;
    }

   /* Zero out upper right triangle of decomposed covariance matrix, which contains the transpose */
    UINT4 i, j;
    for (i = 0; i < kmeans->dim; i++) {
        for (j = i+1; j < kmeans->dim; j++) {
            gsl_matrix_set(kmeans->unwhitened_chol_dec_cov, i, j, 0.);
            gsl_matrix_set(kmeans->chol_dec_cov, i, j, 0.);
        }
    }

    return kmeans;
}


/**
 * Randomly select points from the data as initial centroids for a kmeans run.
 * @param kmeans The kmeans to initialize the centroids of.
 */
void LALInferenceKmeansForgyInitialize(LALInferenceKmeans *kmeans) {
    UINT4 i, j, u;

    gsl_permutation *p = gsl_permutation_alloc(kmeans->npts);

    gsl_permutation_init(p);
    gsl_ran_shuffle(kmeans->rng, p->data, kmeans->npts, sizeof(size_t));

    for (i = 0; i < kmeans->k; i++) {
        u = gsl_permutation_get(p, i);
        for (j = 0; j < kmeans->dim; j++)
            gsl_matrix_set(kmeans->centroids, i, j, gsl_matrix_get(kmeans->data, u, j));
    }

    gsl_permutation_free(p);
}


/**
 * Use the resulting centroids from a random cluster assignment as the initial centroids of a kmeans run.
 * @param kmeans The kmeans to initialize the centroids of.
 */
void LALInferenceKmeansRandomPartitionInitialize(LALInferenceKmeans *kmeans) {
    UINT4 i;
    UINT4 cluster_id;

    for (i = 0; i < kmeans->npts; i++) {
        cluster_id = gsl_rng_uniform_int(kmeans->rng, kmeans->k);
        kmeans->assignments[i] = cluster_id;
    }

    LALInferenceKmeansUpdate(kmeans);
}


/**
 * The assignment step of the kmeans algorithm.
 *
 * Assign all data to the closest centroid and calculate the error, defined
 * as the cumulative sum of the distance between all points and their closest
 * centroid.
 * @param kmeans The kmeans to peform the assignment step on.
 */
void LALInferenceKmeansAssignment(LALInferenceKmeans *kmeans) {
    UINT4 i, j;

    kmeans->error = 0.;

    for (i = 0; i < kmeans->npts; i++) {
        gsl_vector_view x = gsl_matrix_row(kmeans->data, i);
        gsl_vector_view c;

        UINT4 best_cluster = 0;
        REAL8 best_dist = INFINITY;
        REAL8 dist;

        /* Find the closest centroid */
        for (j = 0; j < kmeans->k; j++) {
            c = gsl_matrix_row(kmeans->centroids, j);
            dist = kmeans->dist(&x.vector, &c.vector);

            if (dist < best_dist) {
                best_cluster = j;
                best_dist = dist;
            }
        }

        /* Check if the point's assignment has changed */
        UINT4 current_cluster = kmeans->assignments[i];
        if (best_cluster != current_cluster) {
            kmeans->has_changed = 1;
            kmeans->assignments[i] = best_cluster;
        }
        kmeans->error += best_dist;
    }

    /* Recalculate cluster sizes */
    for (i = 0; i < kmeans->k; i++)
        kmeans->sizes[i] = 0;

    for (i = 0; i < kmeans->npts; i++)
        kmeans->sizes[kmeans->assignments[i]]++;
}


/**
 * The update step of the kmeans algorithm.
 *
 * Based on the current assignments, calculate the new centroid of each cluster.
 * @param kmeans The kmeans to perform the update step on.
 */
void LALInferenceKmeansUpdate(LALInferenceKmeans *kmeans) {
    UINT4 i;

    for (i = 0; i < kmeans->k; i ++) {
        LALInferenceKmeansConstructMask(kmeans, kmeans->mask, i);

        gsl_vector_view c = gsl_matrix_row(kmeans->centroids, i);
        kmeans->centroid(&c.vector, kmeans->data, kmeans->mask);
    }
}



/**
 * Construct a mask to select only the data assigned to a single cluster.
 *
 * Contruct a mask that selects the data from \a kmeans assigned to the centroid indexed
 * in the centroid list by \a centroid_id.
 * @param[in]  kmeans     The kmeans clustering to mask the data of.
 * @param[out] mask       The mask that will select points assigned to \a cluster_id in \a kmeans.
 * @param[in]  cluster_id The index of the centroid in \a kmeans to choose points that are assigned to.
 */
void LALInferenceKmeansConstructMask(LALInferenceKmeans *kmeans, UINT4 *mask, UINT4 cluster_id) {
    UINT4 i;
    for (i = 0; i < kmeans->npts; i++) {
        if (kmeans->assignments[i] == cluster_id)
            mask[i] = 1;
        else
            mask[i] = 0;
    }
}


/**
 * Extract a single cluster from an existing kmeans as a new 1-means.
 *
 * Given the index of the centroid requested, generate a new 1-means structure
 * containing only the points assigned to that cluster.
 * @param[in] kmeans     The parent kmeans to extract the cluster from.
 * @param[in] cluster_id The index of the centroid in \a kmeans corresponding to the cluster to extract.
 * @result A kmeans with k=1 containing only the points assigned to cluster \a cluster_id from \a kmeans.
 */
LALInferenceKmeans *LALInferenceKmeansExtractCluster(LALInferenceKmeans *kmeans, UINT4 cluster_id) {
    /* Construct a mask to select only the cluster requested */
    LALInferenceKmeansConstructMask(kmeans, kmeans->mask, cluster_id);
    gsl_matrix *masked_data = mask_data(kmeans->data, kmeans->mask);

    /* Create the kmeans if cluster isn't empty (i.e. masked_array != NULL) */
    LALInferenceKmeans *sub_kmeans = NULL;
    if (masked_data) {
        sub_kmeans = LALInferenceCreateKmeans(1, masked_data, kmeans->rng);
        gsl_matrix_free(masked_data);
    }

    /* If cluster was extracted successfully, then run it */
    if (sub_kmeans) {
        /* Initialize and assign all points to the only cluster */
        LALInferenceKmeansForgyInitialize(sub_kmeans);
        LALInferenceKmeansRun(sub_kmeans);
    }

    return sub_kmeans;
}


/**
 * Generate a new matrix by masking an existing one.
 *
 * Allocate and fill a new GSL matrix by masking an
 * existing one.
 * @param[in] data A GSL matrix whose rows will be masked.
 * @param[in] mask A 0/1 array with length equal to the number of rows in \a data.
 * @return A new GSL matrix constaining only the masked data.
 */
gsl_matrix *mask_data(gsl_matrix *data, UINT4 *mask) {
    UINT4 i;
    UINT4 N = data->size1;
    UINT4 dim = data->size2;

    /* Calculate how big the resulting matrix will be */
    UINT4 new_N = 0;
    for (i = 0; i < N; i++)
        new_N += mask[i];
 
    /* Return NULL if masked array is empty */
    gsl_matrix *masked_data = NULL;
    if (new_N > 0)
        masked_data = gsl_matrix_alloc(new_N, dim);
    else
        return NULL;

    gsl_vector_view source, target;

    UINT4 idx = 0;
    for (i = 0; i < N; i++) {
        if (mask[i]) {
            source = gsl_matrix_row(data, i);
            target = gsl_matrix_row(masked_data, idx);
            gsl_vector_memcpy(&target.vector, &source.vector);
            idx++;
        }
    }
    return masked_data;
}


/**
 * Add a vector to the end of a matrix, allocating if necessary.
 *
 * Utility for acculating vectors in a GSL matrix, allocating if its
 * the first vector to be accumulated.
 * @param[in] mat_ptr Pointer to the GSL matrix to append to.
 * @param[in] vec The vector to append to the matrix pointed to by \a mat_ptr.
 * \sa append_vec_to_mat()
 */
void accumulate_vectors(gsl_matrix **mat_ptr, gsl_vector *vec) {
    gsl_matrix *mat = *mat_ptr;

    if (mat == NULL) {
        mat = gsl_matrix_alloc(1, vec->size);
        gsl_vector_view row = gsl_matrix_row(mat, 0);
        gsl_vector_memcpy(&row.vector, vec);
        *mat_ptr = mat;
    } else {
        append_vec_to_mat(mat_ptr, vec);
    }
}

/**
 * Add a vector to the end of a matrix.
 *
 * Utility for acculating vectors in a GSL matrix.
 * @param[in] mat_ptr Pointer to the GSL matrix to append to.
 * @param[in] vec The vector to append to the matrix pointed to by \a mat_ptr.
 */
void append_vec_to_mat(gsl_matrix **mat_ptr, gsl_vector *vec) {
    gsl_matrix *mat = *mat_ptr;
    UINT4 rows = mat->size1;
    UINT4 cols = mat->size2;

    gsl_matrix *new_mat = gsl_matrix_alloc(rows+1, cols);

    /* Copy over existing rows */
    gsl_matrix_view sub_mat_view = gsl_matrix_submatrix(new_mat,  0, 0, rows, cols);
    gsl_matrix_memcpy(&sub_mat_view.matrix, mat);

    /* Copy vector into new row */
    gsl_vector_view row_view = gsl_matrix_row(new_mat, rows);
    gsl_vector_memcpy(&row_view.vector, vec);

    gsl_matrix_free(mat);
    *mat_ptr = new_mat;
}


/**
 * Draw a sample from a kmeans-KDE estimate of a distribution.
 *
 * Draw a random sample from the estimated distribution of a set of points.  A cluster
 * from kmeans is picked at random, from which a sample is drawn from the KDE of the
 * cluster.
 * @param[in] kmeans The kmeans to draw a sample from.
 * @returns An array containing a sample drawn from \a kmeans.
 */
REAL8 *LALInferenceKmeansDraw(LALInferenceKmeans *kmeans) {
    /* Draw a cluster at random, using the weights assigned to each cluster */
    REAL8 randomDraw = gsl_rng_uniform(kmeans->rng);

    UINT4 cluster = 0;
    REAL8 cumulativeWeight = kmeans->weights[cluster];
    while (cumulativeWeight < randomDraw) {
        cluster++;
        cumulativeWeight += kmeans->weights[cluster];
    }

    /* Draw a sample from the chosen cluster's KDE */
    REAL8 *white_point = LALInferenceDrawKDESample(kmeans->KDEs[cluster], kmeans->rng);
    gsl_vector_view white_pt = gsl_vector_view_array(white_point, kmeans->dim);

    /* Create an empty array to contain the unwhitened sample */
    REAL8 *point = XLALMalloc(kmeans->dim * sizeof(REAL8));
    gsl_vector_view pt = gsl_vector_view_array(point, kmeans->dim);
    gsl_vector_memcpy(&pt.vector, kmeans->mean);

    /* Recolor the sample */
    gsl_blas_dgemv(CblasNoTrans, 1.0, kmeans->unwhitened_chol_dec_cov, &white_pt.vector, 1.0, &pt.vector);

    XLALFree(white_point);
    return point;
}


/**
 * Evaluate the estimated (log) PDF from a clustered-KDE at a point.
 *
 * Compute the (log) value of the estimated probability density function from
 * the clustered kernel density estimate of a distribution.
 * @param[in] kmeans The kmeans clustering to estimate the PDF from.
 * @param[in] pt     An array containing the point to evaluate the distribution at.
 * @return The estimated value of the PDF at \a pt.
 */
REAL8 LALInferenceKmeansPDF(LALInferenceKmeans *kmeans, REAL8 *pt) {
    /* The vector is first whitened to be consistent with the data stored
     * in the kernel density estimator. */
    gsl_vector *y = gsl_vector_alloc(kmeans->dim);
    gsl_vector_view pt_view = gsl_vector_view_array(pt, kmeans->dim);

    /* Copy the point to a local vector, since it will be overwritten */
    gsl_vector_memcpy(y, &pt_view.vector);

    /* Subtract the mean from the point */
    gsl_vector_sub(y, kmeans->mean);

    /* The Householder solver destroys the matrix, so don't give it the original */
    gsl_matrix *A = gsl_matrix_alloc(kmeans->dim, kmeans->dim);
    gsl_matrix_memcpy(A, kmeans->unwhitened_chol_dec_cov);
    gsl_linalg_HH_svx(A, y);
    
    REAL8 result = LALInferenceWhitenedKmeansPDF(kmeans, y->data);

    gsl_matrix_free(A);
    gsl_vector_free(y);

    return result;
}


/**
 * Evaluate the estimated (log) PDF from a clustered-KDE at an already whitened point.
 *
 * Calculate the probability at a point from the clustered-KDE estimate, assuming the point has
 * already been whitened using the same process as the stored data in kmeans.  This is particularly
 * useful when evaluating the BIC during clustering.
 * @param[in] kmeans The kmeans clustering to estimate the PDF from.
 * @param[in] pt     An array containing the point to evaluate the distribution at.
 * @return The estimated value of the PDF at \a pt.
 */
REAL8 LALInferenceWhitenedKmeansPDF(LALInferenceKmeans *kmeans, REAL8 *pt) {
    UINT4 j;

    if (kmeans->KDEs == NULL)
        LALInferenceKmeansBuildKDE(kmeans);

    REAL8* cluster_pdfs = XLALMalloc(kmeans->k * sizeof(REAL8));
    for (j = 0; j < kmeans->k; j++)
        cluster_pdfs[j] = log(kmeans->weights[j]) + LALInferenceKDEEvaluatePoint(kmeans->KDEs[j], pt);

    REAL8 result = log_add_exps(cluster_pdfs, kmeans->k);
    XLALFree(cluster_pdfs);
    return result;
}


/**
 * Transform a data set to obtain a 0-mean and identity covariance matrix.
 *
 * Determine and execute the transformation of a data set to remove any global correlations
 * in a data set.
 * @param[in] samples The data set to whiten.
 * @return A whitened data set with samples corresponding to \a samples.
 */
gsl_matrix * LALInferenceWhitenSamples(gsl_matrix *samples) {
    UINT4 i, j;
    UINT4 npts = samples->size1;
    UINT4 dim = samples->size2;
    INT4 status;

    gsl_matrix *whitened_samples = gsl_matrix_alloc(npts, dim);
    gsl_matrix_memcpy(whitened_samples, samples);

    /* Calculate the mean and covariance matrix */
    gsl_vector *mean = gsl_vector_alloc(dim);
    gsl_matrix *cov = gsl_matrix_alloc(dim, dim);
    LALInferenceComputeMean(mean, samples);
    LALInferenceComputeCovariance(cov, samples);

    /* Cholesky decompose the covariance matrix and decorrelate the data */
    gsl_matrix *cholesky_decomp_cov = gsl_matrix_alloc(dim, dim);
    gsl_matrix_memcpy(cholesky_decomp_cov, cov);
    status = LALInferenceCholeskyDecompose(cholesky_decomp_cov);

    /* If data couldn't be whitened, return NULL */
    if (status)
        return NULL;

    /* Zero out upper right triangle of decomposed covariance matrix, which contains the transpose */
    for (i = 0; i < dim; i++) {
        for (j = i+1; j < dim; j++)
            gsl_matrix_set(cholesky_decomp_cov, i, j, 0.);
    }

    /* Decorrelate the data */
    gsl_vector_view y;
    gsl_matrix *A = gsl_matrix_alloc(dim, dim);
    for (i=0; i<npts; i++) {
        y = gsl_matrix_row(whitened_samples, i);
        gsl_vector_sub(&y.vector, mean);

        /* The Householder solver destroys the matrix, so don't give it the original */
        gsl_matrix_memcpy(A, cholesky_decomp_cov);
        gsl_linalg_HH_svx(A, &y.vector);
    }

    gsl_matrix_free(cholesky_decomp_cov);
    gsl_matrix_free(A);
    gsl_vector_free(mean);
    gsl_matrix_free(cov);

    return whitened_samples;
}

/**
 * Calculate the squared Euclidean distance bewteen two points.
 *
 * @param[in] x A GSL vector.
 * @param[in] y Another GSL vector.
 * @return The squared Euclidean distance between \a x and \a y.
 */
REAL8 euclidean_dist_squared(gsl_vector *x, gsl_vector *y) {
    UINT4 i;

    REAL8 dist = 0.;
    for (i = 0; i < x->size; i++) {
        REAL8 diff = gsl_vector_get(x, i) - gsl_vector_get(y, i);
        dist += diff * diff;
    }

    return dist;
}


/**
 * Find the centroid of a masked data set.
 *
 * @param[out] centroid The newly determined centroid.
 * @param[in] data The set of points to determin \a centroid from.
 * @param[in] mask The mask to select which samples in \a data to calculate the centroid of.
 */
void euclidean_centroid(gsl_vector *centroid, gsl_matrix *data, UINT4 *mask) {
    UINT4 npts = data->size1;
    UINT4 dim = data->size2;
    UINT4 i;
    UINT4 count = 0;

    for (i = 0; i < dim; i++)
        gsl_vector_set(centroid, i, 0.);

    for (i = 0; i < npts; i++) {
        if (mask[i]) {
            gsl_vector_view x = gsl_matrix_row(data, i);
            gsl_vector_add(centroid, &x.vector);
            count++;
        }
    }

    gsl_vector_scale(centroid, 1./count);
}


/**
 * Build the kernel density estimate from a kmeans clustering.
 *
 * Given the current clustering, estimate the distribution of each cluster
 * using a kernel density estimator.  This will be used to evaluate the Bayes
 * Information Criteria when deciding the optimal clustering.
 * @param kmeans The kmeans to estimate the cluster distributions of.
 */
void LALInferenceKmeansBuildKDE(LALInferenceKmeans *kmeans) {
    UINT4 i;

    if (kmeans->KDEs == NULL)
        kmeans->KDEs = XLALMalloc(kmeans->k * sizeof(LALInferenceKDE*));

    for (i = 0; i < kmeans->k; i++) {
        LALInferenceKmeansConstructMask(kmeans, kmeans->mask, i);
        kmeans->KDEs[i] = LALInferenceNewKDEfromMat(kmeans->data, kmeans->mask);
    }
}


/**
 * Calculate the maximum likelihood of a given kmeans assuming spherical Gaussian clusters.
 *
 * Determine the maximum likelihood estimate (MLE) of the clustering assuming
 * each cluster is drawn from a spherical Gaussian.  This is not currently
 * used, but is the original MLE used by the xmeans clustering algorithm.
 * @param kmeans The kmeans to determine the BIC of.
 * @return The value of the Bayes Information Criterian of the \a kmeans clustering.
 */
REAL8 LALInferenceKmeansMaxLogL(LALInferenceKmeans *kmeans) {
    UINT4 i, j;

    REAL8 n_c;
    REAL8 N = (REAL8) kmeans->npts;
    REAL8 dim = (REAL8) kmeans->dim;
    REAL8 k = (REAL8) kmeans->k;

    /* Calculate the MLE of the variance */
    gsl_vector_view c;
    REAL8 MLE_variance = 0.;
    for (j = 0; j < kmeans->k; j++) {
        LALInferenceKmeansConstructMask(kmeans, kmeans->mask, j);
        c = gsl_matrix_row(kmeans->centroids, j);

        for (i = 0; i < kmeans->npts; i++) {
            if (kmeans->mask[i]) {
                gsl_vector_view x = gsl_matrix_row(kmeans->data, i);
                MLE_variance += euclidean_dist_squared(&x.vector, &c.vector);
            }
        }
    }

    MLE_variance /= N - k;

    /* Use the MLE of the variance to determine the maximum likelihood */
    REAL8 MLE = 0.;
    for (i = 0; i < kmeans->k; i++) {
        n_c = (REAL8) kmeans->sizes[i];
        MLE += -n_c/2. * (log(2.*LAL_PI) - dim * log(MLE_variance) + 2. * (log(n_c) - log(N))) - (n_c - k)/2.;
    }

    return MLE;
}


/**
 * Calculate the BIC using the KDEs of the cluster distributions.
 *
 * Use the kernal density estimate of each clusters distribution to calculate the Bayes
 * Information Criterian (BIC).  The BIC penalizes the maximum likelihood of a given
 * model based on the number of parameters of the model.  This is used to deside the
 * optimal clustering.
 * @param kmeans The kmeans to determine the BIC of.
 * @return The value of the Bayes Information Criterian of the \a kmeans clustering.
 */
REAL8 LALInferenceKmeansBIC(LALInferenceKmeans *kmeans) {

    /* Return -infinity if kmeans isn't allocated */
    if (!kmeans)
        return -INFINITY;

    UINT4 i;
    REAL8 log_l;
    REAL8 k = (REAL8) kmeans->k;
    REAL8 N = (REAL8) kmeans->npts;
    REAL8 d = (REAL8) kmeans->dim;

    log_l = 0.;
    for (i = 0; i < kmeans->npts; i++) {
        gsl_vector_view pt = gsl_matrix_row(kmeans->data, i);
        log_l += LALInferenceWhitenedKmeansPDF(kmeans, (&pt.vector)->data);
    }

    /* Determine the total number of parameters in clustered-KDE */
    /* Account for centroid locations */
    REAL8 nparams = k * d;

    /* One weight for each cluster, minus one for constraint that all sum to unity */
    nparams += k - 1;

    /* Separate kernel covariances for each cluster */
    nparams += k*(d+1)*d/2.0;

    return log_l - nparams/2.0 * log(N);
}


/**
 * Destroy a kmeans instance.
 *
 * Systematically deallocate all memory associated to \a kmeans.
 * @param kmeans The kmeans instance to deallocate.
 */
void LALInferenceKmeansDestroy(LALInferenceKmeans *kmeans) {
    /* Only destroy when there is something to destroy */
    if (kmeans) {
        /* Free GSL matrices and vectors */
        gsl_vector_free(kmeans->mean);
        gsl_matrix_free(kmeans->chol_dec_cov);
        gsl_matrix_free(kmeans->unwhitened_chol_dec_cov);
        gsl_matrix_free(kmeans->centroids);
        if (kmeans->data) gsl_matrix_free(kmeans->data);
        if (kmeans->recursive_centroids) gsl_matrix_free(kmeans->recursive_centroids);

        /* Free non-GSL arrays */
        XLALFree(kmeans->assignments);
        XLALFree(kmeans->mask);
        XLALFree(kmeans->weights);
        XLALFree(kmeans->sizes);

        /* If KDEs are defined, free all of them */
        if (kmeans->KDEs != NULL) {
            UINT4 k;
            for (k=0; k<kmeans->k; k++)
                if (kmeans->KDEs[k]) LALInferenceDestroyKDE(kmeans->KDEs[k]);
            XLALFree(kmeans->KDEs);
        }

        XLALFree(kmeans);
    }
}
