#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>

struct tagkmeans;

/**
 * Structure for performing the kmeans clustering algorithm on a set of samples.
 */
typedef struct
tagkmeans
{
    gsl_matrix *data;                   /**< Data to be clustered (typically whitened) */
    UINT4 dim;                          /**< Dimension of data */
    UINT4 npts;                         /**< Number of points being clustered */
    UINT4 k;                            /**< Number of clusters */
    UINT4 has_changed;                  /**< Flag indicating a change to cluster assignmens */

    REAL8 (*dist) (gsl_vector *x, gsl_vector *y);                           /**< Distance funtion */
    void (*centroid) (gsl_vector *centroid, gsl_matrix *data, UINT4 *mask); /**< Find centroid */

    UINT4 *assignments;                 /**< Cluster assignments */
    UINT4 *sizes;                       /**< Cluster sizes */
    UINT4 *mask;                        /**< Mask used to select data from individual clusters */
    REAL8 *weights;                     /**< Fraction of data points in each cluster */
    gsl_vector *mean;                   /**< Mean of unwhitened data (for tranforming points) */
    gsl_matrix *chol_dec_cov;           /**< Cholesky decomposition of the covariance matrix */
    gsl_matrix *centroids;              /**< Array with rows containing cluster centroids */
    gsl_matrix *recursive_centroids;    /**< Matrix used to accumulate tested centroids */
    gsl_rng *rng;                       /**< Random number generator */

    REAL8 error;                        /**< Error of current clustering */

    LALInferenceKDE **KDEs;             /**< Array of KDEs, one for each cluster */
} LALInferenceKmeans;



/* Kmeans cluster data, increasing k until the Bayes Information Criteria stop increasing. */
LALInferenceKmeans *LALInferenceIncrementalKmeans(gsl_matrix *data);

/* Xmeans cluster data, splitting centroids in kmeans according to the Bayes Information Criteria. */
LALInferenceKmeans *LALInferenceXmeans(gsl_matrix *data);

/* Run kmeans recursively, splitting each centroid recursively, accepting if the BIC increases. */
LALInferenceKmeans *LALInferenceRecursiveKmeans(gsl_matrix *data);

/* Recursively split a k=1 kmeans. */
void LALInferenceKmeansRecursiveSplit(LALInferenceKmeans *kmeans);

/* Run the kmeans algorithm until cluster assignments don't change. */
void LALInferenceKmeansRun(LALInferenceKmeans *kmeans);

/* Run a kmeans several times and return the best. */
LALInferenceKmeans *LALInferenceKmeansRunBestOf(UINT4 k, gsl_matrix *samples, UINT4 iter, gsl_rng *rng);

/* Generate a new kmeans struct from a set of data. */
LALInferenceKmeans * LALInferenceCreateKmeans(UINT4 k, gsl_matrix *data, gsl_rng *rng);

/* Randomly select points from the data as initial centroids for a kmeans run. */
void LALInferenceKmeansForgyInitialize(LALInferenceKmeans *kmeans);

/* Use the resulting centroids from a random cluster assignment as the initial centroids of a kmeans run. */
void LALInferenceKmeansRandomPartitionInitialize(LALInferenceKmeans *kmeans);

/* The assignment step of the kmeans algorithm. */
void LALInferenceKmeansAssignment(LALInferenceKmeans *kmeans);

/* The update step of the kmeans algorithm. */
void LALInferenceKmeansUpdate(LALInferenceKmeans *kmeans);

/* Construct a mask to select only the data assigned to a single cluster. */
void LALInferenceKmeansConstructMask(LALInferenceKmeans *kmeans, UINT4 *mask, UINT4 cluster_id);

/* Extract a single cluster from an existing kmeans as a new 1-means. */
LALInferenceKmeans *LALInferenceKmeansExtractCluster(LALInferenceKmeans *kmeans, UINT4 cluster_id);

/* Generate a new matrix by masking an existing one. */
gsl_matrix *mask_data(gsl_matrix *data, UINT4 *mask);

/* Add a vector to the end of a matrix, allocating if necessary. */
void accumulate_vectors(gsl_matrix **mat_ptr, gsl_vector *vec);

/* Add a vector to the end of a matrix. */
void append_vec_to_mat(gsl_matrix **mat_ptr, gsl_vector *vec);

/* Evaluate the estimated PDF from a clustered-KDE at a point. */
REAL8 LALInferenceKmeansPDF(LALInferenceKmeans *kmeans, REAL8 *pt);

/* Evaluate the estimated PDF from a clustered-KDE at an already whitened point. */
REAL8 LALInferenceInternalKmeansPDF(LALInferenceKmeans *kmeans, REAL8 *pt);

/* Calculate the squared Euclidean distance bewteen two points. */
REAL8 euclidean_dist_squared(gsl_vector *x, gsl_vector *y);

/* Find the centroid of a masked data set. */
void euclidean_centroid(gsl_vector *centroid, gsl_matrix *data, UINT4 *mask);

/* Build the kernel density estimate from a kmeans clustering. */
void LALInferenceKmeansBuildKDE(LALInferenceKmeans *kmeans);

/* Calculate the maximum likelihood of a given kmeans assuming spherical Gaussian clusters. */
REAL8 LALInferenceKmeansMaxLogL(LALInferenceKmeans *kmeans);

/* Calculate the BIC using the KDEs of the cluster distributions. */
REAL8 LALInferenceKmeansBIC(LALInferenceKmeans *kmeans);

/* Destroy a kmeans instance. */
void LALInferenceKmeansDestroy(LALInferenceKmeans *kmeans);
