/*
 *  LALInferenceClusteringKDE.c: Clustering-optimized Gaussian
 *                                 kernel density estimator.
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort_vector_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>

#include <lal/LALInference.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceKDE.h>
#include <lal/LALInferenceClusteredKDE.h>

#ifndef _OPENMP
#define omp ignore
#endif



/**
 * Kmeans cluster data, increasing k until the BIC stops increasing.
 *
 * Run the kmeans clustering algorithm incrementally, calculating the Bayes
 *  Information Criteria (BIC) at each step, stopping when the BIC stops
 *  increasing.  The BIC is calculated by estimating the distribution in each
 *  cluster with a Gaussian kernel density estimate, then weighting that
 *  clusters contribution to the likelihood by the fraction of total samples in
 *  that cluster.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceIncrementalKmeans(gsl_matrix *data,
                                                    INT4 ntrials,
                                                    gsl_rng *rng) {
    INT4 k = 1;
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
 * Kmeans cluster data, maximizing BIC assuming BIC(k) is concave-down.
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
LALInferenceKmeans *LALInferenceOptimizedKmeans(gsl_matrix *data,
                                                INT4 ntrials,
                                                gsl_rng *rng) {
    INT4 k, low_k = 1, mid_k = 2, high_k = 4;
    REAL8 bic, low_bic, mid_bic, high_bic;
    LALInferenceKmeans *low_kmeans = NULL;
    LALInferenceKmeans *mid_kmeans = NULL;
    LALInferenceKmeans *high_kmeans = NULL;

    /* Calculate starting clusters and BICs */
    low_kmeans = LALInferenceKmeansRunBestOf(low_k, data, ntrials, rng);
    if (!low_kmeans) {
        fprintf(stderr, "Can't build basic KDE.");
        fprintf(stderr, " Perhaps sample size (%i) is too small.\n",
                (INT4)data->size1);
        if (mid_kmeans) LALInferenceKmeansDestroy(mid_kmeans);
        if (high_kmeans) LALInferenceKmeansDestroy(high_kmeans);
        return NULL;
    }
    low_bic = LALInferenceKmeansBIC(low_kmeans);

    mid_kmeans = LALInferenceKmeansRunBestOf(mid_k, data, ntrials, rng);

    /* Return vanilla KDE if clustering with k > 1 fails */
    if (!mid_kmeans)
        return low_kmeans;

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
            high_kmeans =
                LALInferenceKmeansRunBestOf(high_k, data, ntrials, rng);
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
 * Cluster data xmeans-style, splitting centroids according to the BIC.
 *
 * Run an xmeans-like clustering, increasing a kmeans clustering starting from
 *  k=1 by splitting each centroid individually and checking for an increase of
 *  the BIC over the parent cluster.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceXmeans(gsl_matrix *data,
                                        INT4 ntrials,
                                        gsl_rng *rng) {
    INT4 kmax = 50;
    INT4 split_size = 2;
    INT4 c,k;
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

    /* Keep increasing k by splitting centroids.  Stop if BIC is maximized */
    while (kmeans->k < kmax) {
        if (kmeans->recursive_centroids != NULL)
            gsl_matrix_free(kmeans->recursive_centroids);

        /* Store the BIC of the initial clustering */
        starting_bic = LALInferenceKmeansBIC(kmeans);

        /* For each existing centroid, split and see if BIC increases locally */
        for (c = 0; c < kmeans->k; c++) {
            sub_kmeans = LALInferenceKmeansExtractCluster(kmeans, c);

            /* Reject split if cluster extraction failed. */
            if (sub_kmeans) {
                old_bic = LALInferenceKmeansBIC(sub_kmeans);

                new_kmeans = LALInferenceKmeansRunBestOf(split_size,
                                                            sub_kmeans->data,
                                                            ntrials, rng);
                new_bic = LALInferenceKmeansBIC(new_kmeans);
            } else {
                old_bic = INFINITY;
                new_bic = -INFINITY;
            }

            if (new_bic > old_bic) {
                /* BIC increased, so keep the new centroids */
                for (k = 0; k < new_kmeans->k; k++) {
                    gsl_vector_view centroid =
                        gsl_matrix_row(new_kmeans->centroids, k);

                    accumulate_vectors(&kmeans->recursive_centroids,
                                        &centroid.vector);
                }
            } else {
                /* BIC decreased, keep the parent centroid */
                gsl_vector_view centroid = gsl_matrix_row(kmeans->centroids, c);

                accumulate_vectors(&kmeans->recursive_centroids,
                                    &centroid.vector);
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

            /* Store BIC after all centroids have been attempted to be split */
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
 * Run kmeans, recursively splitting centroids, accepting if the BIC increases.
 *
 * Run kmeans recursively, splitting each centroid in two recursively.  At each
 *  step, the BIC is checked before and after splitting.  If the split
 *  increases the BIC, each centroid is then attempted to be split.  If the BIC
 *  decreases, the split is rejected and the centroid is kept.  This is less
 *  precise than other methods, but much faster, as the likelihood is computed
 *  at few and fewer points as centroids are split.
 * @param[in] data A GSL matrix containing the data to be clustered.
 * @param[in] ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 * @param[in] rng  A GSL random number generator used by the kmeans algorithm.
 * @result A kmeans structure containing the clustering that maximizes the BIC.
 */
LALInferenceKmeans *LALInferenceRecursiveKmeans(gsl_matrix *data,
                                                INT4 ntrials,
                                                gsl_rng *rng) {
    INT4 k;

    /* Perform recursive splitting.  Return NULL if kmeans creation fails */
    LALInferenceKmeans *split_kmeans = LALInferenceCreateKmeans(1, data, rng);
    if (!split_kmeans)
        return NULL;

    LALInferenceKmeansForgyInitialize(split_kmeans);
    LALInferenceKmeansRecursiveSplit(split_kmeans, ntrials, rng);

    /* Take the final centroids and make a fully self-consistent kmeans */
    gsl_matrix *centroids = split_kmeans->recursive_centroids;

    k = centroids->size1;
    LALInferenceKmeans *kmeans = LALInferenceCreateKmeans(k, data, rng);

    if (!kmeans)
        return NULL;

    gsl_matrix_memcpy(kmeans->centroids, centroids);
    LALInferenceKmeansRun(kmeans);

    LALInferenceKmeansDestroy(split_kmeans);
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
void LALInferenceKmeansRecursiveSplit(LALInferenceKmeans *kmeans,
                                        INT4 ntrials,
                                        gsl_rng *rng) {
    LALInferenceKmeansRun(kmeans);
    INT4 cluster, i;
    INT4 split_size = 2;
    REAL8 current_bic, new_bic;

    LALInferenceKmeans *new_kmeans =
        LALInferenceKmeansRunBestOf(split_size, kmeans->data, ntrials, rng);

    current_bic = LALInferenceKmeansBIC(kmeans);
    new_bic = LALInferenceKmeansBIC(new_kmeans);

    /* If the BIC increased, attempt a further splitting of the centroids */
    if (new_bic > current_bic) {
        for (cluster = 0; cluster < new_kmeans->k; cluster++) {
            LALInferenceKmeans *sub_kmeans =
                LALInferenceKmeansExtractCluster(new_kmeans, cluster);

            if (sub_kmeans)
                LALInferenceKmeansRecursiveSplit(sub_kmeans, ntrials, rng);

            INT4 n_new_centroids = sub_kmeans->recursive_centroids->size1;
            for (i = 0; i < n_new_centroids ; i++) {
                gsl_vector_view centroid =
                    gsl_matrix_row(sub_kmeans->recursive_centroids, i);

                accumulate_vectors(&kmeans->recursive_centroids,
                                    &centroid.vector);
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
 * Starting with some random initialization, points are assigned to the closest
 *  centroids, then the centroids are calculated of the new cluster.  This is
 *  repeated until the assignments stop changing.
 * @param kmeans The initialized kmeans to run.
 */
void LALInferenceKmeansRun(LALInferenceKmeans *kmeans) {
    INT4 i;

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
 * The kmeans is run \a ntrials times, each time with a new random
 *  initialization.  The one that results in the highest Bayes Information
 *  Criteria (BIC) is returned.
 * @param[in]  k       The number of clusters to use.
 * @param[in]  samples The (unwhitened) data to cluster.
 * @param[in]  ntrials The number of random initialization to run.
 * @param[in]  rng     The GSL random number generator to use for random
 *                      initializations.
 * @return The kmeans with the highest BIC of \a ntrials attempts.
 */
LALInferenceKmeans *LALInferenceKmeansRunBestOf(INT4 k,
                                                gsl_matrix *samples,
                                                INT4 ntrials,
                                                gsl_rng *rng) {
    INT4 i;

    LALInferenceKmeans *best_kmeans = NULL;
    REAL8 best_bic = -INFINITY;

    /* Don't bother with multiple trials if k=1 */
    if (k == 1)
        ntrials = 1;

    #pragma omp parallel
    {
        LALInferenceKmeans *kmeans;
        REAL8 bic = 0;
        REAL8 error = -INFINITY;

        #pragma omp for
        for (i = 0; i < ntrials; i++) {
            kmeans = LALInferenceCreateKmeans(k, samples, rng);
            if (!kmeans)
                continue;

            LALInferenceKmeansSeededInitialize(kmeans);
            LALInferenceKmeansRun(kmeans);

            /* Assume BIC hasn't changed if
             * error (summed-dist^2 from assigned centroids) hasn't */
            if (error != kmeans->error) {
                error = kmeans->error;
                bic = LALInferenceKmeansBIC(kmeans);
            }

            #pragma omp critical
            {
                if (bic > best_bic) {
                    if (best_kmeans)
                        LALInferenceKmeansDestroy(best_kmeans);
                    best_bic = bic;
                    best_kmeans = kmeans;
                } else {
                    LALInferenceKmeansDestroy(kmeans);
                }
            }
        }
    }

    return best_kmeans;
}

/**
 * Generate a new kmeans struct from a set of data.
 *
 * Whitens data and fills a newly allocated kmeans structure of the
 *  requested size.
 * @param[in] k    The requested number of clusters.
 * @param[in] data The unwhitened data to be clustered.
 * @param[in] rng  A GSL random number generator used for initialization.
 */
LALInferenceKmeans *LALInferenceCreateKmeans(INT4 k,
                                                gsl_matrix *data,
                                                gsl_rng *rng) {
    INT4 i;

    /* Return nothing if given empty dataset */
    if (!data)
        return NULL;

    LALInferenceKmeans *kmeans = XLALCalloc(1, sizeof(LALInferenceKmeans));
    kmeans->k = k;
    kmeans->npts = data->size1;
    kmeans->dim = data->size2;
    kmeans->has_changed = 1;
    kmeans->rng = rng;

    /* Allocate GSL matrices and vectors */
    kmeans->mean = gsl_vector_alloc(kmeans->dim);
    kmeans->cov  = gsl_matrix_alloc(kmeans->dim, kmeans->dim);
    kmeans->std  = gsl_vector_alloc(kmeans->dim);
    kmeans->data = gsl_matrix_alloc(kmeans->npts, kmeans->dim);
    kmeans->centroids = gsl_matrix_alloc(kmeans->k, kmeans->dim);
    kmeans->recursive_centroids = NULL;

    /* Allocate everything else */
    kmeans->assignments = XLALCalloc(kmeans->npts, sizeof(INT4));
    kmeans->sizes = XLALCalloc(kmeans->k, sizeof(INT4));
    kmeans->mask = XLALCalloc(kmeans->npts, sizeof(INT4));
    kmeans->weights = XLALCalloc(kmeans->k, sizeof(REAL8));
    kmeans->KDEs = NULL;

    /* Set distance and centroid calculators */
    kmeans->dist = &euclidean_dist_squared;
    kmeans->centroid = &euclidean_centroid;

    /* Store the unwhitened mean and covariance matrix for later transformations */
    LALInferenceComputeMean(kmeans->mean, data);
    LALInferenceComputeCovariance(kmeans->cov, data);

    for (i = 0; i < kmeans->dim; i++)
        gsl_vector_set(kmeans->std, i, sqrt(gsl_matrix_get(kmeans->cov, i, i)));

    /* Whiten the data */
    kmeans->data = LALInferenceWhitenSamples(data);
    if (!kmeans->data) {
        fprintf(stderr, "Unable to whiten data. No proposal has been built.\n");
        LALInferenceKmeansDestroy(kmeans);
        return NULL;
    }


    return kmeans;
}


/**
 * A 'k-means++'-like seeded initialization of centroids for a kmeans run.
 * @param kmeans The kmeans to initialize the centroids of.
 */
void LALInferenceKmeansSeededInitialize(LALInferenceKmeans *kmeans) {
    INT4 i, j;
    INT4 u = 0;
    REAL8 dist, norm, randomDraw;
    gsl_vector_view c, x;

    gsl_vector *dists = gsl_vector_alloc(kmeans->npts);
    gsl_permutation *p = gsl_permutation_alloc(kmeans->npts);

    /* Choose first centroid randomly from data */
    i = gsl_rng_uniform_int(kmeans->rng, kmeans->npts);

    for (j = 0; j < kmeans->dim; j++)
        gsl_matrix_set(kmeans->centroids, u, j,
                        gsl_matrix_get(kmeans->data, i, j));

    u++;
    while (u < kmeans->k) {
        /* Calculate distances from centroid */
        norm = 0.0;
        c = gsl_matrix_row(kmeans->centroids, u-1);
        for (i = 0; i < kmeans->npts; i++) {
            x = gsl_matrix_row(kmeans->data, i);
            dist = kmeans->dist(&x.vector, &c.vector);

            gsl_vector_set(dists, i, dist);
            norm += dist;
        }

        gsl_vector_scale(dists, 1.0/norm);
        gsl_sort_vector_index(p, dists);

        randomDraw = gsl_rng_uniform(kmeans->rng);

        i = 0;
        dist = gsl_vector_get(dists, p->data[i]);
        while (dist < randomDraw) {
            i++;
            dist += gsl_vector_get(dists, p->data[i]);
        }

        for (j = 0; j < kmeans->dim; j++)
            gsl_matrix_set(kmeans->centroids, u, j,
                            gsl_matrix_get(kmeans->data, i, j));

        u++;
    }

    gsl_permutation_free(p);
}


/**
 * Randomly select points from the data as initial centroids for a kmeans run.
 * @param kmeans The kmeans to initialize the centroids of.
 */
void LALInferenceKmeansForgyInitialize(LALInferenceKmeans *kmeans) {
    INT4 i, j, u;

    gsl_permutation *p = gsl_permutation_alloc(kmeans->npts);

    gsl_permutation_init(p);
    gsl_ran_shuffle(kmeans->rng, p->data, kmeans->npts, sizeof(size_t));

    for (i = 0; i < kmeans->k; i++) {
        u = gsl_permutation_get(p, i);
        for (j = 0; j < kmeans->dim; j++)
            gsl_matrix_set(kmeans->centroids, i, j,
                            gsl_matrix_get(kmeans->data, u, j));
    }

    gsl_permutation_free(p);
}


/**
 * Use the resulting centroids from a random cluster assignment as the initial
 *  centroids of a kmeans run.
 * @param kmeans The kmeans to initialize the centroids of.
 */
void LALInferenceKmeansRandomPartitionInitialize(LALInferenceKmeans *kmeans) {
    INT4 i;
    INT4 cluster_id;

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
 * @param kmeans The kmeans to perform the assignment step on.
 */
void LALInferenceKmeansAssignment(LALInferenceKmeans *kmeans) {
    INT4 i, j;

    kmeans->error = 0.;

    for (i = 0; i < kmeans->npts; i++) {
        gsl_vector_view x = gsl_matrix_row(kmeans->data, i);
        gsl_vector_view c;

        INT4 best_cluster = 0;
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
        INT4 current_cluster = kmeans->assignments[i];
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
    INT4 i;

    for (i = 0; i < kmeans->k; i ++) {
        LALInferenceKmeansConstructMask(kmeans, kmeans->mask, i);

        gsl_vector_view c = gsl_matrix_row(kmeans->centroids, i);
        kmeans->centroid(&c.vector, kmeans->data, kmeans->mask);
    }
}



/**
 * Construct a mask to select only the data assigned to a single cluster.
 *
 * Contruct a mask that selects the data from \a kmeans assigned to the centroid
 *  indexed in the centroid list by \a centroid_id.
 * @param[in]  kmeans     The kmeans clustering to mask the data of.
 * @param[out] mask       The mask that will select points assigned to
 *                          \a cluster_id in \a kmeans.
 * @param[in]  cluster_id The index of the centroid in \a kmeans to choose points
 *                          that are assigned to.
 */
void LALInferenceKmeansConstructMask(LALInferenceKmeans *kmeans,
                                        INT4 *mask,
                                        INT4 cluster_id) {
    INT4 i;
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
 * @param[in] cluster_id The index of the centroid in \a kmeans corresponding
 *                          to the cluster to extract.
 * @result A kmeans with k=1 containing only the points assigned to cluster
 *          \a cluster_id from \a kmeans.
 */
LALInferenceKmeans *LALInferenceKmeansExtractCluster(LALInferenceKmeans *kmeans,
                                                        INT4 cluster_id) {
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
 * Impose boundaries on individual KDEs.
 *
 * Draw samples from each cluster.  If too many samples lie outside of the
 * prior, impose a cyclic/reflective bound for the offending parameter(s) if
 * requested.  If boundaries aren't incountered, box in the cluster using the
 * samples drawn.  This avoids needing to evaluate KDEs that are too far away.
 * @param[in] kmeans             kmeans to cycle through the clusters of.
 * @param[in] params             Parameters to impose bounds on.
 * @param[in] priorArgs          Variables containing prior boundaries.
 * @param[in] cyclic_reflective  Flag to check for cyclic/reflective bounds.
 *
 */
void LALInferenceKmeansImposeBounds(LALInferenceKmeans *kmeans,
                                    LALInferenceVariables *params,
                                    LALInferenceVariables *priorArgs,
                                    INT4 cyclic_reflective) {
    INT4 i, p, dim;
    INT4 n_below, n_above;
    INT4 ndraws, n_thresh;
    INT4 c;
    REAL8 draw;
    REAL8 min, max, threshold;
    REAL8 drawn_min = INFINITY, drawn_max = -INFINITY;
    REAL8 mean, std;
    LALInferenceVariableItem *param=NULL;

    dim = kmeans->dim;

    ndraws = 1000;      // Number of draws from the KDE to check
    threshold = 0.01;   // Fraction outside boundary to make it cyclic/reflective
    n_thresh = (INT4)(threshold * ndraws);

    REAL8 *draws = XLALCalloc(ndraws * dim, sizeof(REAL8));
    for (c = 0; c < kmeans->k; c++) {
        /* Skip empty clusters */
        if (kmeans->npts == 0)
            continue;

        for (i = 0; i < ndraws; i++) {
            /* Draw a whitened sample from the cluster's KDE */
            REAL8 *point =
                LALInferenceDrawKDESample(kmeans->KDEs[c], kmeans->rng);

            memcpy(&(draws[i*dim]), point, dim*sizeof(REAL8));
            XLALFree(point);
        }

        p = 0;
        for (param = params->head; param; param = param->next) {
            if (LALInferenceCheckVariableNonFixed(params, param->name)) {
                /* Whiten prior boundaries */
                LALInferenceGetMinMaxPrior(priorArgs, param->name, &min, &max);
                mean = gsl_vector_get(kmeans->mean, p);
                std = gsl_vector_get(kmeans->std, p);
                min = (min - mean)/std;
                max = (max - mean)/std;

                n_below = 0;
                n_above = 0;
                for (i = 0; i < ndraws; i++) {
                    draw = draws[i*dim + p];
                    if (draw < min)
                        n_below++;
                    else if (draw > max)
                        n_above++;

                    if (draw < drawn_min)
                        drawn_min = draw;
                    else if (draw > drawn_max)
                        drawn_max = draw;
                }

                if (cyclic_reflective && n_below > n_thresh) {
                    /* Impose cyclic boundaries on both sides */
                    if (param->vary == LALINFERENCE_PARAM_CIRCULAR) {
                        kmeans->KDEs[c]->lower_bound_types[p] = param->vary;
                        kmeans->KDEs[c]->upper_bound_types[p] = param->vary;
                    } else {
                        kmeans->KDEs[c]->lower_bound_types[p] = param->vary;
                    }
                } else {
                    min = drawn_min;
                    kmeans->KDEs[c]->lower_bound_types[p] = LALINFERENCE_PARAM_FIXED;
                }

                if (cyclic_reflective && n_above > n_thresh) {
                    /* Impose cyclic boundaries on both sides */
                    if (param->vary == LALINFERENCE_PARAM_CIRCULAR) {
                        kmeans->KDEs[c]->lower_bound_types[p] = param->vary;
                        kmeans->KDEs[c]->upper_bound_types[p] = param->vary;
                    } else {
                        kmeans->KDEs[c]->upper_bound_types[p] = param->vary;
                    }
                } else {
                    max = drawn_max;
                    kmeans->KDEs[c]->upper_bound_types[p] = LALINFERENCE_PARAM_FIXED;
                }

                kmeans->KDEs[c]->lower_bounds[p] = min;
                kmeans->KDEs[c]->upper_bounds[p] = max;

                p++;
            }
        }
    }
}

/**
 * Generate a new matrix by masking an existing one.
 *
 * Allocate and fill a new GSL matrix by masking an
 * existing one.
 * @param[in] data GSL matrix whose rows will be masked.
 * @param[in] mask 0/1 array with length equal to the number of rows in \a data.
 * @return A new GSL matrix constaining only the masked data.
 */
gsl_matrix *mask_data(gsl_matrix *data, INT4 *mask) {
    INT4 i;
    INT4 N = data->size1;
    INT4 dim = data->size2;

    /* Calculate how big the resulting matrix will be */
    INT4 new_N = 0;
    for (i = 0; i < N; i++)
        new_N += mask[i];

    /* Return NULL if masked array is empty */
    gsl_matrix *masked_data = NULL;
    if (new_N > 0)
        masked_data = gsl_matrix_alloc(new_N, dim);
    else
        return NULL;

    gsl_vector_view source, target;

    INT4 idx = 0;
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
    INT4 rows = mat->size1;
    INT4 cols = mat->size2;

    gsl_matrix *new_mat = gsl_matrix_alloc(rows+1, cols);

    /* Copy over existing rows */
    gsl_matrix_view sub_mat_view =
        gsl_matrix_submatrix(new_mat,  0, 0, rows, cols);
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
 * Draw a random sample from the estimated distribution of a set of points.
 *  A cluster from kmeans is picked at random, from which a sample is drawn from
 *  the KDE of the cluster.
 * @param[in] kmeans The kmeans to draw a sample from.
 * @returns An array containing a sample drawn from \a kmeans.
 */
REAL8 *LALInferenceKmeansDraw(LALInferenceKmeans *kmeans) {
    /* Draw a cluster at random, using the weights assigned to each cluster */
    REAL8 randomDraw = gsl_rng_uniform(kmeans->rng);

    INT4 cluster = 0;
    REAL8 cumulativeWeight = kmeans->weights[cluster];
    while (cumulativeWeight < randomDraw) {
        cluster++;
        cumulativeWeight += kmeans->weights[cluster];
    }

    /* Draw a whitened sample from the chosen cluster's KDE */
    REAL8 *point =
        LALInferenceDrawKDESample(kmeans->KDEs[cluster], kmeans->rng);

    /* Recolor the sample */
    gsl_vector_view pt = gsl_vector_view_array(point, kmeans->dim);
    gsl_vector_mul(&pt.vector, kmeans->std);
    gsl_vector_add(&pt.vector, kmeans->mean);

    return point;
}


/**
 * Evaluate the estimated (log) PDF from a clustered-KDE at a point.
 *
 * Compute the (log) value of the estimated probability density function from
 * the clustered kernel density estimate of a distribution.
 * @param[in] kmeans The kmeans clustering to estimate the PDF from.
 * @param[in] pt     Array containing the point to evaluate the distribution at.
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
    gsl_vector_div(y, kmeans->std);

    REAL8 result = LALInferenceWhitenedKmeansPDF(kmeans, y->data);

    gsl_vector_free(y);

    return result;
}


/**
 * Estimate (log) PDF from a clustered-KDE at an already whitened point.
 *
 * Calculate the probability at a point from the clustered-KDE estimate,
 *  assuming the point has already been whitened using the same process as the
 *  stored data in kmeans.  This is particularly useful when evaluating the BIC
 *  during clustering.
 * @param[in] kmeans The kmeans clustering to estimate the PDF from.
 * @param[in] pt     Array containing the point to evaluate the distribution at.
 * @return The estimated value of the PDF at \a pt.
 */
REAL8 LALInferenceWhitenedKmeansPDF(LALInferenceKmeans *kmeans, REAL8 *pt) {
    INT4 j;

    if (kmeans->KDEs == NULL)
        LALInferenceKmeansBuildKDE(kmeans);

    REAL8* cluster_pdfs = XLALCalloc(kmeans->k, sizeof(REAL8));
    for (j = 0; j < kmeans->k; j++)
        cluster_pdfs[j] = log(kmeans->weights[j]) +
                            LALInferenceKDEEvaluatePoint(kmeans->KDEs[j], pt);

    REAL8 result = log_add_exps(cluster_pdfs, kmeans->k);
    XLALFree(cluster_pdfs);
    return result;
}


/**
 * Transform a data set to obtain a 0-mean and identity covariance matrix.
 *
 * Determine and execute the transformation of a data set to remove any global
 *  correlations in a data set.
 * @param[in] samples The data set to whiten.
 * @return A whitened data set with samples corresponding to \a samples.
 */
gsl_matrix * LALInferenceWhitenSamples(gsl_matrix *samples) {
    INT4 i;
    INT4 npts = samples->size1;
    INT4 dim = samples->size2;

    gsl_matrix *whitened_samples = gsl_matrix_alloc(npts, dim);
    gsl_matrix_memcpy(whitened_samples, samples);

    /* Calculate the mean and covariance matrix */
    gsl_vector *mean = gsl_vector_alloc(dim);
    gsl_matrix *cov = gsl_matrix_alloc(dim, dim);
    gsl_vector *std = gsl_vector_alloc(dim);
    LALInferenceComputeMean(mean, samples);
    LALInferenceComputeCovariance(cov, samples);

    for (i = 0; i < dim; i++)
        gsl_vector_set(std, i, sqrt(gsl_matrix_get(cov, i, i)));

    /* Decorrelate the data */
    gsl_vector_view y;
    for (i=0; i<npts; i++) {
        y = gsl_matrix_row(whitened_samples, i);
        gsl_vector_sub(&y.vector, mean);
        gsl_vector_div(&y.vector, std);
    }

    gsl_vector_free(mean);
    gsl_matrix_free(cov);
    gsl_vector_free(std);

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
    size_t i;

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
 * @param[in] data Set of points to determin \a centroid from.
 * @param[in] mask Mask to select which samples in \a data to calculate the
 *                  centroid of.
 */
void euclidean_centroid(gsl_vector *centroid, gsl_matrix *data, INT4 *mask) {
    INT4 npts = data->size1;
    INT4 dim = data->size2;
    INT4 i;
    INT4 count = 0;

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
    INT4 i;

    if (kmeans->KDEs == NULL)
        kmeans->KDEs = XLALCalloc(kmeans->k, sizeof(LALInferenceKDE*));

    for (i = 0; i < kmeans->k; i++) {
        LALInferenceKmeansConstructMask(kmeans, kmeans->mask, i);
        kmeans->KDEs[i] = LALInferenceNewKDEfromMat(kmeans->data, kmeans->mask);
    }
}


/**
 * Calculate max likelihood of a kmeans assuming spherical Gaussian clusters.
 *
 * Determine the maximum likelihood estimate (MLE) of the clustering assuming
 * each cluster is drawn from a spherical Gaussian.  This is not currently
 * used, but is the original MLE used by the xmeans clustering algorithm.
 * @param kmeans The kmeans to determine the BIC of.
 * @return The maximum likelihood estimate of the \a kmeans clustering.
 */
REAL8 LALInferenceKmeansMaxLogL(LALInferenceKmeans *kmeans) {
    INT4 i, j;

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
        MLE += -n_c/2. * (log(2.*LAL_PI) -
                dim * log(MLE_variance) +
                2. * (log(n_c) - log(N))) - (n_c - k)/2.;
    }

    return MLE;
}


/**
 * Calculate the BIC using the KDEs of the cluster distributions.
 *
 * Use the kernal density estimate of each clusters distribution to calculate
 *  the Bayes Information Criterian (BIC).  The BIC penalizes the maximum
 *  likelihood of a given model based on the number of parameters of the model.
 *  This is used to deside the optimal clustering.
 * @param kmeans The kmeans to determine the BIC of.
 * @return The Bayes Information Criterian of the \a kmeans clustering.
 */
REAL8 LALInferenceKmeansBIC(LALInferenceKmeans *kmeans) {

    /* Return -infinity if kmeans isn't allocated */
    if (!kmeans)
        return -INFINITY;

    INT4 i;
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

    /* One weight for each cluster, minus one for constraint that
     *  all sum to unity */
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
        gsl_matrix_free(kmeans->centroids);
        if (kmeans->data) gsl_matrix_free(kmeans->data);
        if (kmeans->recursive_centroids)
            gsl_matrix_free(kmeans->recursive_centroids);

        /* Free non-GSL arrays */
        XLALFree(kmeans->assignments);
        XLALFree(kmeans->mask);
        XLALFree(kmeans->weights);
        XLALFree(kmeans->sizes);

        /* If KDEs are defined, free all of them */
        if (kmeans->KDEs != NULL) {
            INT4 k;
            for (k=0; k<kmeans->k; k++)
                if (kmeans->KDEs[k]) LALInferenceDestroyKDE(kmeans->KDEs[k]);
            XLALFree(kmeans->KDEs);
        }

        XLALFree(kmeans);
    }
}
