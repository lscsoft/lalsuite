/*
 *  Copyright (C) 2017 Michael Puerrer, Ben Lackey
 *  Reduced Order Model for SEOBNRv4T spin-aligned tidal effective-one-body model
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdbool.h>
#include <alloca.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>


#ifdef LAL_HDF5_ENABLED
#include <lal/H5FileIO.h>
static const char SurDataHDF5[] = "SEOBNRv4T_surrogate_v1.0.0.hdf5";
static const INT4 SurDataHDF5_VERSION_MAJOR = 1;
static const INT4 SurDataHDF5_VERSION_MINOR = 0;
static const INT4 SurDataHDF5_VERSION_MICRO = 0;
#endif

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMREOBNRv2.h"
#include "LALSimUniversalRelations.h"
#include "LALSimInspiralPNCoefficients.c"
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif



#ifdef LAL_PTHREAD_LOCK
static pthread_once_t Surrogate_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

struct tagSurrogatedata_submodel
{
  gsl_matrix *hyp_amp;           // GP hyperparameters log amplitude
  gsl_matrix *hyp_phi;           // GP hyperparameters for dephasing
  gsl_matrix *kinv_dot_y_amp;    // kinv_dot_y for log amplitude
  gsl_matrix *kinv_dot_y_phi;    // kinv_dot_y for dephasing
  gsl_matrix *x_train;           // Training points

  // location of surrogate spline nodes
  gsl_vector *mf_amp;            // location of spline nodes for log amplitude
  gsl_vector *mf_phi;            // location of spline nodes for dephasing

  // location of spline nodes for TaylorF2 at low frequency
  gsl_vector *TF2_mf_amp_cubic;  // cubic spline nodes for TaylorF2 amplitude
  gsl_vector *TF2_mf_phi_cubic;  // cubic spline nodes for TaylorF2 phasing
  gsl_vector *TF2_mf_amp_linear; // linear spline nodes for TaylorF2 amplitude
  gsl_vector *TF2_mf_phi_linear; // linear spline nodes for TaylorF2 phasing

  // 5D parameter space bounds of surrogate
  gsl_vector *q_bounds;          // [q_min, q_max]
  gsl_vector *chi1_bounds;       // [chi1_min, chi1_max]
  gsl_vector *chi2_bounds;       // [chi2_min, chi2_max]
  gsl_vector *lambda1_bounds;    // [lambda1_min, lambda1_max]
  gsl_vector *lambda2_bounds;    // [lambda2_min, lambda2_max]
};
typedef struct tagSurrogatedata_submodel Surrogatedata_submodel;

struct tagSurrogatedata
{
  UINT4 setup;
  Surrogatedata_submodel* sub1;
};
typedef struct tagSurrogatedata Surrogatedata;

static Surrogatedata __lalsim_SurrogateDS_data;

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

/**************** Internal functions **********************/

UNUSED static void save_gsl_frequency_series(const char filename[], gsl_vector *x, gsl_vector *y);

static gsl_vector *gsl_vector_prepend_value(gsl_vector *v, double value);

double kernel(
  gsl_vector *x1,          // array with shape ndim
  gsl_vector *x2,          // array with shape ndim
  gsl_vector *hyperparams, // Hyperparameters
  gsl_vector *work         // work array
);

double gp_predict(
  gsl_vector *xst,          // Point x_* where you want to evaluate the function.
  gsl_vector *hyperparams,  // Hyperparameters for the GPR kernel.
  gsl_matrix *x_train,      // Training set points.
  gsl_vector *Kinv_dot_y,   // The interpolating weights at each training set point.
  gsl_vector *work          // work array for kernel
);


UNUSED static void Surrogate_Init_LALDATA(void);
UNUSED static int Surrogate_Init(const char dir[]);
UNUSED static bool Surrogate_IsSetup(void);

UNUSED static int Surrogatedata_Init(Surrogatedata *romdata, const char dir[]);
UNUSED static void Surrogatedata_Cleanup(Surrogatedata *romdata);

static double xi_of_lambda(double lambda);

static int GPR_evaluation_5D(
  double q,                      // Input: q-value (q >= 1)
  double chi1,                   // Input: chi1-value
  double chi2,                   // Input: chi2-value
  double lambda1,                // Input: lambda1-value
  double lambda2,                // Input: lambda2-value
  gsl_matrix *hyp_amp,           // Input: GPR hyperparameters for log amplitude
  gsl_matrix *hyp_phi,           // Input: GPR hyperparameters for dephasing
  gsl_matrix *kinv_dot_y_amp,    // Input: kinv_dot_y for log amplitude
  gsl_matrix *kinv_dot_y_phi,    // Input: kinv_dot_y for dephasing
  gsl_matrix *x_train,           // Input: GPR training points
  gsl_vector *amp_at_EI_nodes,   // Output: log amplitude at EI nodes (preallocated)
  gsl_vector *phi_at_EI_nodes,   // Output: dephasing at EI nodes     (preallocated)
  gsl_vector *work               // Input: work array
);

UNUSED static int Surrogatedata_Init_submodel(
  UNUSED Surrogatedata_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[]
);

UNUSED static void Surrogatedata_Cleanup_submodel(Surrogatedata_submodel *submodel);

UNUSED static int CheckParameterSpaceBounds(
  Surrogatedata_submodel *sur,
  double q,      // mass-ratio q >= 1
  double chi1,
  double chi2,
  double lambda1,
  double lambda2
);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
UNUSED static int SurrogateCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double chi1,
  double chi2,
  double lambda1,
  double lambda2,
  const REAL8Sequence *freqs, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  SEOBNRv4TSurrogate_spline_order spline_order
);

static size_t NextPow2(const size_t n);

static int TaylorF2Phasing(
  double Mtot,         // Total mass in solar masses
  double q,            // Mass-ration m1/m2 >= 1
  double chi1,         // Dimensionless aligned spin of body 1
  double chi2,         // Dimensionless aligned spin of body 2
  double lambda1,      // Tidal deformability of body 1
  double lambda2,      // Tidal deformability of body 2
  gsl_vector *Mfs,     // Input geometric frequencies
  gsl_vector **PNphase // Output: TaylorF2 phase at frequencies Mfs
);

static int TaylorF2Amplitude1PN(
  double eta,        // Symmetric mass-ratio
  gsl_vector *Mfs,   // Input geometric frequencies
  gsl_vector **PNamp // Output: TaylorF2 amplitude at frequencies Mfs
);

/********************* Definitions begin here ********************/

UNUSED static void save_gsl_frequency_series(const char filename[], gsl_vector *x, gsl_vector *y) {
  FILE *fp = fopen(filename, "w");
  fprintf(stderr, "save_gsl_frequency_series: %s: %zu %zu\n", filename, x->size, y->size);
  assert(x->size == y->size);
  for (size_t i=0; i<x->size; i++) {
    fprintf(fp, "%g\t%g\n", gsl_vector_get(x, i), gsl_vector_get(y, i));
  }
  fclose(fp);
}

static gsl_vector *gsl_vector_prepend_value(gsl_vector *v, double value) {
// Helper function to prepend a value to a gsl_vector
// Returns the augmented gsl_vector
// Deallocates the input gsl_vector
  int n = v->size;
  gsl_vector *vout = gsl_vector_alloc(n+1);

  gsl_vector_set(vout, 0, value);
  for (int i=1; i<=n; i++)
    gsl_vector_set(vout, i, gsl_vector_get(v, i-1));
  gsl_vector_free(v);

  return vout;
}

double kernel(
  gsl_vector *x1,          // parameter space point 1
  gsl_vector *x2,          // parameter space point 2
  gsl_vector *hyperparams, // GPR Hyperparameters
  gsl_vector *work         // work array
)
{
  // Matern covariance function for n-dimensional data.
  //
  // Parameters
  // ----------
  // x1 : array with shape ndim
  // x2 : array with shape ndim
  // hyperparams : array with shape ndim+2 [sigma_f, ls0, ls1, ..., sigma_n]
  //     sigma_f : Approximately the range (ymax-ymin) of values that the data takes.
  //         sigma_f^2 called the signal variance.
  //     sigma_n : Noise term. The uncertainty in the y values of the data.
  //     lsi : Length scales for the variation in dimension i.
  //
  // Returns
  // -------
  // covariance : double

  double sigma_f = gsl_vector_get(hyperparams, 0);
  double sigma_n = gsl_vector_get(hyperparams, hyperparams->size-1);
  gsl_vector ls = gsl_vector_subvector(hyperparams, 1, hyperparams->size-2).vector;

  XLAL_CHECK((x1->size == x2->size) && (x1->size == ls.size), XLAL_EDIMS,
  "kernel(): dimensions of vectors x1, x2 and ls: %zu, %zu, %zu have to be consistent.\n",
  x1->size, x2->size, ls.size);

  // Noise nugget for diagonal elements
  double nugget = 0.0;
  if (gsl_vector_equal(x1, x2))
    nugget = sigma_n*sigma_n;

  gsl_vector_memcpy(work, x1);
  gsl_vector_sub(work, x2);
  gsl_vector_div(work, &ls); // (x1 - x2) / ls
  double r = gsl_blas_dnrm2(work);

  // nu = 5/2 Matern covariance
  double matern = (1.0 + sqrt(5.0)*r + 5.0*r*r/3.0) * exp(-sqrt(5.0)*r);

  // Full covariance
  // Include the nugget to agree with scikit-learn when the points x1, x2 are exactly the same.
  return sigma_f*sigma_f * matern + nugget;
}

double gp_predict(
  gsl_vector *xst,          // Point x_* where you want to evaluate the function.
  gsl_vector *hyperparams,  // Hyperparameters for the GPR kernel.
  gsl_matrix *x_train,      // Training set points.
  gsl_vector *Kinv_dot_y,   // The interpolating weights at each training set point.
  gsl_vector *work
)
{
  // Interpolate the function at the point xst using Gaussian process regression.
  //
  // Parameters
  // ----------
  // xst : array of shape ndim.
  //     Point x_* where you want to evaluate the function.
  // hyperparams : array with shape ndim+2 [sigma_f, ls0, ls1, ..., sigma_n].
  //     Hyperparameters for the GPR kernel.
  // x_train : array of shape (n_train, ndim).
  //     Training set points.
  // Kinv_dot_y : array of shape n_train.
  //     The interpolating weights at each training set point.
  //
  // Returns
  // -------
  // yst : double
  //     Interpolated value at the point xst.

  // Evaluate vector K_*
  int n = x_train->size1;
  gsl_vector *Kst = gsl_vector_alloc(n);
  for (int i=0; i < n; i++) {
    gsl_vector x = gsl_matrix_const_row(x_train, i).vector;
    double ker = kernel(xst, &x, hyperparams, work);
    gsl_vector_set(Kst, i, ker);
  }

  // Evaluate y_*
  double res = 0;
  gsl_blas_ddot(Kst, Kinv_dot_y, &res);
  gsl_vector_free(Kst);
  return res;
}

/** Setup Surrogate model using data files installed in dir */
static int Surrogate_Init(const char dir[]) {
  if(__lalsim_SurrogateDS_data.setup) {
    XLALPrintError("Error: Surrogate data was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }
  Surrogatedata_Init(&__lalsim_SurrogateDS_data, dir);

  if(__lalsim_SurrogateDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the Surrogate model has been initialised */
static bool Surrogate_IsSetup(void) {
  if(__lalsim_SurrogateDS_data.setup)
    return true;
  else
    return false;
}

/** Coordinate transformation */
static double xi_of_lambda(double lambda) {
  const double a = 100.0;
  return log10(lambda / a + 1.0);
}

// Compute amplitude and phase at empirical interpolant nodes from GPRs.
// This entails interpolation over the 5D parameter space (q, chi1, chi2, lambda1, lambda2).
static int GPR_evaluation_5D(
  double q,                      // Input: q-value (q >= 1)
  double chi1,                   // Input: chi1-value
  double chi2,                   // Input: chi2-value
  double lambda1,                // Input: lambda1-value
  double lambda2,                // Input: lambda2-value
  gsl_matrix *hyp_amp,           // Input: GPR hyperparameters for log amplitude
  gsl_matrix *hyp_phi,           // Input: GPR hyperparameters for dephasing
  gsl_matrix *kinv_dot_y_amp,    // Input: kinv_dot_y for log amplitude
  gsl_matrix *kinv_dot_y_phi,    // Input: kinv_dot_y for dephasing
  gsl_matrix *x_train,           // Input: GPR training points
  gsl_vector *amp_at_nodes,      // Output: log amplitude at frequency nodes (preallocated)
  gsl_vector *phi_at_nodes,      // Output: dephasing at frequency nodes (preallocated)
  gsl_vector *work               // Input: workspace
)
{
  // Assemble evaluation point
  gsl_vector *xst = gsl_vector_alloc(5);
  double q_inv = 1.0/q;
  gsl_vector_set(xst, 0, q_inv);
  gsl_vector_set(xst, 1, chi1);
  gsl_vector_set(xst, 2, chi2);
  gsl_vector_set(xst, 3, xi_of_lambda(lambda1));
  gsl_vector_set(xst, 4, xi_of_lambda(lambda2));

  // Evaluate GPR for amplitude spline nodes
  for (size_t i=0; i<amp_at_nodes->size; i++) {
    gsl_vector hyp_amp_i = gsl_matrix_const_row(hyp_amp, i).vector;
    gsl_vector kinv_dot_y_amp_i = gsl_matrix_const_row(kinv_dot_y_amp, i).vector;
    double pred = gp_predict(xst, &hyp_amp_i, x_train, &kinv_dot_y_amp_i, work);
    gsl_vector_set(amp_at_nodes, i, pred);
  }

  // Evaluate GPR for phase spline nodes
  for (size_t i=0; i<phi_at_nodes->size; i++) {
    gsl_vector hyp_phi_i = gsl_matrix_const_row(hyp_phi, i).vector;
    gsl_vector kinv_dot_y_phi_i = gsl_matrix_const_row(kinv_dot_y_phi, i).vector;
    double pred = gp_predict(xst, &hyp_phi_i, x_train, &kinv_dot_y_phi_i, work);
    gsl_vector_set(phi_at_nodes, i, pred);
  }

  gsl_vector_free(xst);

  return XLAL_SUCCESS;
}

/* Set up a new submodel, using data contained in dir */
UNUSED static int Surrogatedata_Init_submodel(
  Surrogatedata_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[]
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(Surrogatedata_submodel));
  else
    Surrogatedata_Cleanup_submodel(*submodel);

#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(SurDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, SurDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  //LALH5File *file = XLALH5FileOpen("/Users/mpuer/Documents/gpsurrogate/src/TEOB-LAL-implementation/TEOBv4_surrogate.hdf5", "r");
  LALH5File *root = XLALH5GroupOpen(file, "/");

  //////////////////////////////////////////////////////////////////////////////
  // load everything we need
  // GP hyperparameters
  ReadHDF5RealMatrixDataset(root, "hyp_amp", & (*submodel)->hyp_amp);
  ReadHDF5RealMatrixDataset(root, "hyp_phi", & (*submodel)->hyp_phi);

  // K^{-1} . y
  ReadHDF5RealMatrixDataset(root, "kinv_dot_y_amp", & (*submodel)->kinv_dot_y_amp);
  ReadHDF5RealMatrixDataset(root, "kinv_dot_y_phi", & (*submodel)->kinv_dot_y_phi);

  // Training points
  ReadHDF5RealMatrixDataset(root, "x_train", & (*submodel)->x_train);

  // Frequency grids for surrogate corrections
  ReadHDF5RealVectorDataset(root, "spline_nodes_amp", & (*submodel)->mf_amp);
  ReadHDF5RealVectorDataset(root, "spline_nodes_phase", & (*submodel)->mf_phi);

  // Frequency grids for TaylorF2
  ReadHDF5RealVectorDataset(root, "TF2_Mf_amp_cubic", & (*submodel)->TF2_mf_amp_cubic);
  ReadHDF5RealVectorDataset(root, "TF2_Mf_phi_cubic", & (*submodel)->TF2_mf_phi_cubic);
  ReadHDF5RealVectorDataset(root, "TF2_Mf_amp_linear", & (*submodel)->TF2_mf_amp_linear);
  ReadHDF5RealVectorDataset(root, "TF2_Mf_phi_linear", & (*submodel)->TF2_mf_phi_linear);

  // Physical domain covered by surrogate
  ReadHDF5RealVectorDataset(root, "q_bounds", & (*submodel)->q_bounds);
  ReadHDF5RealVectorDataset(root, "chi1_bounds", & (*submodel)->chi1_bounds);
  ReadHDF5RealVectorDataset(root, "chi2_bounds", & (*submodel)->chi2_bounds);
  ReadHDF5RealVectorDataset(root, "lambda1_bounds", & (*submodel)->lambda1_bounds);
  ReadHDF5RealVectorDataset(root, "lambda2_bounds", & (*submodel)->lambda2_bounds);

  // Prepend the point [mf_amp[0], 0] to the phase nodes
  double mf_min = gsl_vector_get( (*submodel)->mf_amp, 0); // Follow definition of mf_a in GPSplineSurrogate constructor
  gsl_vector *phi_nodes = gsl_vector_prepend_value((*submodel)->mf_phi, mf_min);
  (*submodel)->mf_phi = phi_nodes;

  XLALFree(path);
  XLALH5FileClose(file);
  ret = XLAL_SUCCESS;
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  return ret;
}

/* Deallocate contents of the given Surrogatedata_submodel structure */
static void Surrogatedata_Cleanup_submodel(Surrogatedata_submodel *submodel) {
  if(submodel->hyp_amp) gsl_matrix_free(submodel->hyp_amp);
  if(submodel->hyp_phi) gsl_matrix_free(submodel->hyp_phi);
  if(submodel->kinv_dot_y_amp) gsl_matrix_free(submodel->kinv_dot_y_amp);
  if(submodel->kinv_dot_y_phi) gsl_matrix_free(submodel->kinv_dot_y_phi);
  if(submodel->x_train) gsl_matrix_free(submodel->x_train);
  if(submodel->mf_amp) gsl_vector_free(submodel->mf_amp);
  if(submodel->mf_phi) gsl_vector_free(submodel->mf_phi);
  if(submodel->TF2_mf_amp_cubic) gsl_vector_free(submodel->TF2_mf_amp_cubic);
  if(submodel->TF2_mf_phi_cubic) gsl_vector_free(submodel->TF2_mf_phi_cubic);
  if(submodel->TF2_mf_amp_linear) gsl_vector_free(submodel->TF2_mf_amp_linear);
  if(submodel->TF2_mf_phi_linear) gsl_vector_free(submodel->TF2_mf_phi_linear);

  if(submodel->q_bounds) gsl_vector_free(submodel->q_bounds);
  if(submodel->chi1_bounds) gsl_vector_free(submodel->chi1_bounds);
  if(submodel->chi2_bounds) gsl_vector_free(submodel->chi2_bounds);
  if(submodel->lambda1_bounds) gsl_vector_free(submodel->lambda1_bounds);
  if(submodel->lambda2_bounds) gsl_vector_free(submodel->lambda2_bounds);
}

/* Set up a new ROM model, using data contained in dir */
int Surrogatedata_Init(
  UNUSED Surrogatedata *romdata,
  UNUSED const char dir[])
{
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the Surrogate model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

#ifdef LAL_HDF5_ENABLED
  // First, check we got the correct version number
  size_t size = strlen(dir) + strlen(SurDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, SurDataHDF5);
  LALH5File *file = XLALH5FileOpen(path, "r");

  XLALPrintInfo("Surrogate metadata\n============\n");
  PrintInfoStringAttribute(file, "Email");
  PrintInfoStringAttribute(file, "Description");
  ret = ROM_check_version_number(file, SurDataHDF5_VERSION_MAJOR,
                                       SurDataHDF5_VERSION_MINOR,
                                       SurDataHDF5_VERSION_MICRO);

  XLALFree(path);
  XLALH5FileClose(file);

  ret = Surrogatedata_Init_submodel(&(romdata)->sub1, dir, "sub1");
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 1 loaded successfully.\n", __func__);

  if(XLAL_SUCCESS==ret)
    romdata->setup=1;
  else
    Surrogatedata_Cleanup(romdata);
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  return (ret);
}

/* Deallocate contents of the given Surrogatedata structure */
static void Surrogatedata_Cleanup(Surrogatedata *romdata) {
  Surrogatedata_Cleanup_submodel((romdata)->sub1);
  XLALFree((romdata)->sub1);
  (romdata)->sub1 = NULL;
  romdata->setup=0;
}

/* Return the closest higher power of 2  */
// Note: NextPow(2^k) = 2^k for integer values k.
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}


static int TaylorF2Phasing(
  double Mtot,         // Total mass in solar masses
  double q,            // Mass-ration m1/m2 >= 1
  double chi1,         // Dimensionless aligned spin of body 1
  double chi2,         // Dimensionless aligned spin of body 2
  double lambda1,      // Tidal deformability of body 1
  double lambda2,      // Tidal deformability of body 2
  gsl_vector *Mfs,     // Input geometric frequencies
  gsl_vector **PNphase // Output: TaylorF2 phase at frequencies Mfs
) {
  XLAL_CHECK(PNphase != NULL, XLAL_EFAULT);
  XLAL_CHECK(*PNphase == NULL, XLAL_EFAULT);
  *PNphase = gsl_vector_alloc(Mfs->size);

  PNPhasingSeries *pn = NULL;
  LALDict *extraParams = XLALCreateDict();
  XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
  XLALSimInspiralWaveformParamsInsertPNTidalOrder(extraParams, LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN);
  XLALSimInspiralWaveformParamsInsertTidalLambda1(extraParams, lambda1);
  XLALSimInspiralWaveformParamsInsertTidalLambda2(extraParams, lambda2);

  // Use the universal relation fit to determine the quadrupole-monopole terms.
  double dquadmon1 = XLALSimUniversalRelationQuadMonVSlambda2Tidal(lambda1) - 1.0;
  double dquadmon2 = XLALSimUniversalRelationQuadMonVSlambda2Tidal(lambda2) - 1.0;

  if ((dquadmon1 > 0) || (dquadmon2 > 0)) {
    XLALPrintInfo("Using quadrupole-monopole terms from fit: %g, %g\n", dquadmon1, dquadmon2);
    XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams, dquadmon1);
    XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams, dquadmon2);
  }

  double m1OverM = q / (1.0+q);
  double m2OverM = 1.0 / (1.0+q);
  double m1 = Mtot * m1OverM * LAL_MSUN_SI;
  double m2 = Mtot * m2OverM * LAL_MSUN_SI;
  // This function is a thin wrapper around XLALSimInspiralPNPhasing_F2
  XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, extraParams);

  // Compute and subtract pn_ss3 term (See LALSimInspiralPNCoefficients.c: XLALSimInspiralPNPhasing_F2()).
  // Rationale: SEOBNRv4T does not contain this ss3 term, therefore 
  // we remove it from the TF2 phasing that is used as a base waveform for the surrogate.
  double m1M = m1OverM;
  double m2M = m2OverM;
  double chi1L = chi1;
  double chi2L = chi2;
  double chi1sq = chi1*chi1;
  double chi2sq = chi2*chi2;
  double qm_def1 = 1 + dquadmon1;
  double qm_def2 = 1 + dquadmon2;
  double eta = q / ((1.0 + q)*(1.0 + q));
  double pn_ss3 = XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(eta)*chi1L*chi2L
    + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m1M))*chi1sq
    + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m2M))*chi2sq;
  pn->v[6] -= pn_ss3 * pn->v[0];

  for (size_t i=0; i < Mfs->size; i++) {
      const double Mf = gsl_vector_get(Mfs, i);
      const double v = cbrt(LAL_PI * Mf);
      const double logv = log(v);
      const double v2 = v * v;
      const double v3 = v * v2;
      const double v4 = v * v3;
      const double v5 = v * v4;
      const double v6 = v * v5;
      const double v7 = v * v6;
      const double v8 = v * v7;
      const double v9 = v * v8;
      const double v10 = v * v9;
      const double v12 = v2 * v10;
      double phasing = 0.0;

      phasing += pn->v[7] * v7;
      phasing += (pn->v[6] + pn->vlogv[6] * logv) * v6;
      phasing += (pn->v[5] + pn->vlogv[5] * logv) * v5;
      phasing += pn->v[4] * v4;
      phasing += pn->v[3] * v3;
      phasing += pn->v[2] * v2;
      phasing += pn->v[1] * v;
      phasing += pn->v[0];

      /* Tidal terms in phasing */
      phasing += pn->v[12] * v12;
      phasing += pn->v[10] * v10;

      phasing /= v5;

      gsl_vector_set(*PNphase, i, -phasing);
  }

  XLALDestroyDict(extraParams);
  XLALFree(pn);

  return XLAL_SUCCESS;
}

// 1PN point-particle amplitude.
// Expression from Eq. (6.10) of arXiv:0810.5336.
static int TaylorF2Amplitude1PN(
  double eta,        // Symmetric mass-ratio
  gsl_vector *Mfs,   // Input geometric frequencies
  gsl_vector **PNamp // Output: TaylorF2 amplitude at frequencies Mfs
) {
  XLAL_CHECK(PNamp != NULL, XLAL_EFAULT);
  XLAL_CHECK(*PNamp == NULL, XLAL_EFAULT);
  *PNamp = gsl_vector_alloc(Mfs->size);

  for (size_t i=0; i < Mfs->size; i++) {
    const double Mf = gsl_vector_get(Mfs, i);
    const double v = cbrt(LAL_PI * Mf);
    const double x = v*v;

    double a00 = sqrt( (5.0*LAL_PI/24.0)*eta );
    double a10 = -323.0/224.0 + 451.0*eta/168.0;
    double amp =  -a00 * pow(x, -7.0/4.0) * (1.0 + a10*x);
    gsl_vector_set(*PNamp, i, amp);
  }

  return XLAL_SUCCESS;
}

static int CheckParameterSpaceBounds(
  Surrogatedata_submodel *sur,
  double q,      // mass-ratio q >= 1
  double chi1,
  double chi2,
  double lambda1,
  double lambda2
) {
  // convert to q >= 1
  double q_max = 1.0 / gsl_vector_get(sur->q_bounds, 0);
  double q_min = 1.0 / gsl_vector_get(sur->q_bounds, 1);

  double chi1_min = gsl_vector_get(sur->chi1_bounds, 0);
  double chi1_max = gsl_vector_get(sur->chi1_bounds, 1);
  double chi2_min = gsl_vector_get(sur->chi2_bounds, 0);
  double chi2_max = gsl_vector_get(sur->chi2_bounds, 1);

  double lambda1_min = gsl_vector_get(sur->lambda1_bounds, 0);
  double lambda1_max = gsl_vector_get(sur->lambda1_bounds, 1);
  double lambda2_min = gsl_vector_get(sur->lambda2_bounds, 0);
  double lambda2_max = gsl_vector_get(sur->lambda2_bounds, 1);

  if (q < q_min  || q > q_max) {
    XLALPrintError("XLAL Error - %s: mass-ratio q (%f) out of bounds: [%f, %f]!\n", __func__, q, q_min, q_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (chi1 < chi1_min || chi1 > chi1_max) {
    XLALPrintError("XLAL Error - %s: aligned-spin chi1 (%f) out of bounds: [%f, %f]!\n", __func__, chi1, chi1_min, chi1_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (chi2 < chi2_min || chi2 > chi2_max) {
    XLALPrintError("XLAL Error - %s: aligned-spin chi2 (%f) out of bounds: [%f, %f]!\n", __func__, chi2, chi2_min, chi2_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (lambda1 < lambda1_min || lambda1 > lambda1_max) {
    XLALPrintError("XLAL Error - %s: tidal deformability lambda1 (%f) out of bounds: [%f, %f]!\n", __func__, lambda1, lambda1_min, lambda1_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (lambda2 < lambda2_min || lambda2 > lambda2_max) {
    XLALPrintError("XLAL Error - %s: tidal deformability lambda2 (%f) out of bounds: [%f, %f]!\n", __func__, lambda2, lambda2_min, lambda2_max);
    XLAL_ERROR( XLAL_EDOM );
  }

  return XLAL_SUCCESS;
}

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi1, chi2).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SurrogateCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef, // orbital reference phase
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double eta,
  double chi1,
  double chi2,
  double lambda1,
  double lambda2,
  const REAL8Sequence *freqs_in, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  SEOBNRv4TSurrogate_spline_order spline_order
  )
{
  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);

  Surrogatedata *romdata=&__lalsim_SurrogateDS_data;
  if (!Surrogate_IsSetup()) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up Surrogate data - check your $LAL_DATA_PATH\n");
  }

  if(*hptilde || *hctilde) {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",
                   (*hptilde), (*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  double Mtot = Mtot_sec / LAL_MTSUN_SI;
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);

  Surrogatedata_submodel *sur;
  sur = romdata->sub1;

  retcode = CheckParameterSpaceBounds(sur, q, chi1, chi2, lambda1, lambda2);
  if(retcode != 0) XLAL_ERROR(retcode);

  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];

  if(fRef==0.0)
    fRef=fLow;

  int N_amp = sur->mf_amp->size;
  int N_phi = sur->mf_phi->size;

  // Point to linear or cubic spline nodes and set order for gsl calls lateron
  gsl_vector *TF2_mf_amp = NULL;
  gsl_vector *TF2_mf_phi = NULL;
  const gsl_interp_type *spline_type = NULL;
  switch (spline_order) {
    case SEOBNRv4TSurrogate_CUBIC:
      TF2_mf_amp = sur->TF2_mf_amp_cubic;
      TF2_mf_phi = sur->TF2_mf_phi_cubic;
      spline_type = gsl_interp_cspline;
      break;
    case SEOBNRv4TSurrogate_LINEAR:
      TF2_mf_amp = sur->TF2_mf_amp_linear;
      TF2_mf_phi = sur->TF2_mf_phi_linear;
      spline_type = gsl_interp_linear;
      break;
    default:
      XLAL_ERROR( XLAL_EINVAL, "Unknown spline order!\n. Must be LINEAR or CUBIC." );
  }
  int N_amp_TF2 = TF2_mf_amp->size;
  int N_phi_TF2 = TF2_mf_phi->size;

  // Allowed geometric frequency ranges for
  // surrogate and supported TaylorF2 low frequency extension
  double Mf_TF2_min = fmax(gsl_vector_get(TF2_mf_amp, 0),
                         gsl_vector_get(TF2_mf_phi, 0));
  double Mf_TF2_max = fmin(gsl_vector_get(TF2_mf_amp, N_amp_TF2-1),
                         gsl_vector_get(TF2_mf_phi, N_phi_TF2-1));
  double Mf_sur_min = fmax(gsl_vector_get(sur->mf_amp, 0),
                           gsl_vector_get(sur->mf_phi, 0));
  double Mf_sur_max = fmin(gsl_vector_get(sur->mf_amp, N_amp-1),
                           gsl_vector_get(sur->mf_phi, N_phi-1));

  // sanity checks: sparse grids for TaylorF2 need to contain the
  // frequency range where the surrogate corrections are applied
  XLAL_CHECK(Mf_TF2_min < Mf_sur_min, XLAL_EFAULT);
  XLAL_CHECK(Mf_TF2_max >= Mf_sur_max, XLAL_EFAULT);

  double fLow_geom = fLow * Mtot_sec;
  double fHigh_geom = fHigh * Mtot_sec;
  double fRef_geom = fRef * Mtot_sec;
  double deltaF_geom = deltaF * Mtot_sec;

  // Enforce allowed geometric frequency range
  if (fLow_geom < Mf_TF2_min)
    XLAL_ERROR(XLAL_EDOM, "Starting frequency Mflow=%g is smaller than lowest frequency in surrogate Mf=%g.\n", fLow_geom, Mf_TF2_min);
  if (fHigh_geom == 0 || fHigh_geom > Mf_sur_max)
    fHigh_geom = Mf_sur_max;
  else if (fHigh_geom < Mf_sur_min)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than surrogate starting frequency %g!\n", fHigh_geom, Mf_TF2_min);
  if (fHigh_geom <= fLow_geom)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than (or equal to) starting frequency %g!\n", fHigh_geom, fLow_geom);
  if (fRef_geom > Mf_sur_max) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is greater than maximal frequency in surrogate Mf=%g. Starting at maximal frequency in ROM.\n", fRef_geom, Mf_sur_max);
    fRef_geom = Mf_sur_max; // If fref > fhigh we reset fref to default value of cutoff frequency.
  }
  if (fRef_geom < Mf_TF2_min) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in surrogate Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_TF2_min);
    fRef_geom = Mf_TF2_min;
  }

  // Evaluate GPR for log amplitude and dephasing
  gsl_vector *sur_amp_at_nodes = gsl_vector_alloc(N_amp);
  gsl_vector *sur_phi_at_nodes_tmp = gsl_vector_alloc(N_phi - 1); // Will prepend a point below

  // Allocate workspace for the kernel function
  gsl_vector *work = gsl_vector_alloc(5);

  retcode |= GPR_evaluation_5D(
    q, chi1, chi2, lambda1, lambda2,
    sur->hyp_amp,
    sur->hyp_phi,
    sur->kinv_dot_y_amp,
    sur->kinv_dot_y_phi,
    sur->x_train,
    sur_amp_at_nodes,
    sur_phi_at_nodes_tmp,
    work
  );

  gsl_vector_free(work);

  if(retcode != 0) {
    gsl_vector_free(sur_amp_at_nodes);
    gsl_vector_free(sur_phi_at_nodes_tmp);
    XLAL_ERROR(retcode);
  }

  // Prepend the point [mf_min, 0] to the phase nodes
  // This has already been done in the setup for mf_phi
  gsl_vector *sur_phi_at_nodes = gsl_vector_prepend_value(sur_phi_at_nodes_tmp, 0.0);

  // Spline the surrogate amplitude and phase corrections in frequency
  gsl_interp_accel *acc_amp_sur = gsl_interp_accel_alloc();
  gsl_spline *spline_amp_sur = gsl_spline_alloc(gsl_interp_cspline, N_amp);
  gsl_spline_init(spline_amp_sur, gsl_vector_const_ptr(sur->mf_amp, 0),
                  gsl_vector_const_ptr(sur_amp_at_nodes, 0), N_amp);

  gsl_interp_accel *acc_phi_sur = gsl_interp_accel_alloc();
  gsl_spline *spline_phi_sur = gsl_spline_alloc(gsl_interp_cspline, N_phi);
  gsl_spline_init(spline_phi_sur, gsl_vector_const_ptr(sur->mf_phi, 0),
                  gsl_vector_const_ptr(sur_phi_at_nodes, 0), N_phi);

  gsl_vector_free(sur_amp_at_nodes);
  gsl_vector_free(sur_phi_at_nodes);


  // Evaluate TaylorF2 amplitude and phase on sparse grids
  gsl_vector *TF2_phi_at_nodes = NULL;
  retcode |= TaylorF2Phasing(Mtot, q, chi1, chi2, lambda1, lambda2, TF2_mf_phi, &TF2_phi_at_nodes);

  gsl_vector *TF2_amp_at_nodes = NULL;
  retcode |= TaylorF2Amplitude1PN(eta, TF2_mf_amp, &TF2_amp_at_nodes);

  if(retcode != 0) {
    gsl_vector_free(TF2_amp_at_nodes);
    gsl_vector_free(TF2_phi_at_nodes);
    XLAL_ERROR(retcode);
  }

  // Spline TaylorF2 amplitude and phase
  // We add the surrogate corrections to the TaylorF2 spline data
  // in the frequency range where the surrogate has support.

  // amplitude
  gsl_interp_accel *acc_amp_TF2 = gsl_interp_accel_alloc();
  gsl_spline *spline_amp_TF2 = gsl_spline_alloc(spline_type, N_amp_TF2);
  gsl_vector *spline_amp_values = gsl_vector_alloc(N_amp_TF2);
  for (int i=0; i<N_amp_TF2; i++) {
    double Mf = gsl_vector_get(TF2_mf_amp, i);
    double amp_i = gsl_vector_get(TF2_amp_at_nodes, i);
    if ((Mf >= Mf_sur_min) & (Mf <= Mf_sur_max))
      amp_i *= exp(gsl_spline_eval(spline_amp_sur, Mf, acc_amp_sur));
    gsl_vector_set(spline_amp_values, i, amp_i);
  }
  gsl_spline_init(spline_amp_TF2, gsl_vector_const_ptr(TF2_mf_amp, 0),
                  gsl_vector_const_ptr(spline_amp_values, 0), N_amp_TF2);

  // phase
  gsl_interp_accel *acc_phi_TF2 = gsl_interp_accel_alloc();
  gsl_spline *spline_phi_TF2 = gsl_spline_alloc(spline_type, N_phi_TF2);
  gsl_vector *spline_phi_values = gsl_vector_alloc(N_phi_TF2);
  for (int i=0; i<N_phi_TF2; i++) {
    double Mf = gsl_vector_get(TF2_mf_phi, i);
    double phi_i = gsl_vector_get(TF2_phi_at_nodes, i);
    if ((Mf >= Mf_sur_min) & (Mf <= Mf_sur_max))
      phi_i += gsl_spline_eval(spline_phi_sur, Mf, acc_phi_sur);
    gsl_vector_set(spline_phi_values, i, phi_i);
  }
  gsl_spline_init(spline_phi_TF2, gsl_vector_const_ptr(TF2_mf_phi, 0),
                  gsl_vector_const_ptr(spline_phi_values, 0), N_phi_TF2);

  gsl_vector_free(TF2_amp_at_nodes);
  gsl_vector_free(TF2_phi_at_nodes);
  gsl_vector_free(spline_amp_values);
  gsl_vector_free(spline_phi_values);

  gsl_spline_free(spline_amp_sur);
  gsl_spline_free(spline_phi_sur);
  gsl_interp_accel_free(acc_amp_sur);
  gsl_interp_accel_free(acc_phi_sur);


  // Now setup LAL datastructures for waveform polarizations
  size_t npts = 0;
  LIGOTimeGPS tC = {0, 0};
  UINT4 offset = 0; // Index shift between freqs and the frequency series
  REAL8Sequence *freqs = NULL;
  if (deltaF > 0)  { // freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0
    /* Set up output array with size closest power of 2 */
    npts = NextPow2(fHigh_geom / deltaF_geom) + 1;
    if (fHigh_geom < fHigh * Mtot_sec) /* Resize waveform if user wants f_max larger than cutoff frequency */
      npts = NextPow2(fHigh * Mtot_sec / deltaF_geom) + 1;

    XLALGPSAdd(&tC, -1. / deltaF);  /* coalesce at t=0 */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, npts);

    // Recreate freqs using only the lower and upper bounds
    // Use fLow, fHigh and deltaF rather than geometric frequencies for numerical accuracy
    double fHigh_temp = fHigh_geom / Mtot_sec;
    UINT4 iStart = (UINT4) ceil(fLow / deltaF);
    UINT4 iStop = (UINT4) ceil(fHigh_temp / deltaF);
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs) {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF_geom;

    offset = iStart;
  } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    npts = freqs_in->length;
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    offset = 0;

    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    if (!freqs) {
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    }
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i] * Mtot_sec;
  }

  if (!(*hptilde) || !(*hctilde))	{
      XLALDestroyREAL8Sequence(freqs);
      gsl_interp_accel_free(acc_phi_TF2);
      gsl_spline_free(spline_phi_TF2);
      gsl_interp_accel_free(acc_amp_TF2);
      gsl_spline_free(spline_amp_TF2);
      XLAL_ERROR(XLAL_EFUNC);
  }
  memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));

  XLALUnitMultiply(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  XLALUnitMultiply(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  COMPLEX16 *pdata=(*hptilde)->data->data;
  COMPLEX16 *cdata=(*hctilde)->data->data;

  REAL8 cosi = cos(inclination);
  REAL8 pcoef = 0.5*(1.0 + cosi*cosi);
  REAL8 ccoef = cosi;

  double amp0 = Mtot * Mtot_sec * LAL_MRSUN_SI / (distance); // amplitude prefactor

  // Evaluate reference phase for setting phiRef correctly
  double phase_change = 0.0;
  phase_change = gsl_spline_eval(spline_phi_TF2, fRef_geom, acc_phi_TF2) - 2*phiRef;

  // The maximum frequency for which we generate waveform data is set to the
  // maximum frequency covered by the surrogate.
  Mf_sur_max = gsl_vector_get(sur->mf_phi, N_phi-1);
  double Mf_final = Mf_sur_max;

  // Assemble waveform from aplitude and phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i];
    if (f > Mf_final) continue; // We're beyond the highest allowed frequency;
    // since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
    int j = i + offset; // shift index for frequency series if needed
    double A = gsl_spline_eval(spline_amp_TF2, f, acc_amp_TF2);
    double phase = gsl_spline_eval(spline_phi_TF2, f, acc_phi_TF2) - phase_change;
    COMPLEX16 htilde = amp0*A * (cos(phase) + I*sin(phase));//cexp(I*phase);
    pdata[j] =      pcoef * htilde;
    cdata[j] = -I * ccoef * htilde;
  }

  /* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */
  UINT4 L = freqs->length;
  // prevent gsl interpolation errors
  if (Mf_final > freqs->data[L-1])
    Mf_final = freqs->data[L-1];
  if (Mf_final < freqs->data[0]) {
    XLALDestroyREAL8Sequence(freqs);
    gsl_interp_accel_free(acc_phi_TF2);
    gsl_spline_free(spline_phi_TF2);
    gsl_interp_accel_free(acc_amp_TF2);
    gsl_spline_free(spline_amp_TF2);
    XLAL_ERROR(XLAL_EDOM, "f_ringdown < f_min");
  }

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // We compute the dimensionless time correction t/M since we use geometric units.
  // We use the Schwarzschild ISCO as a rough proxy for the amplitude peak of the BNS waveform.
  // We used XLALSimNSNSMergerFreq() earlier and it turned out not to be reliable for extreme input parameters. 
  REAL8 Mf_ISCO_Schwrazschild = 1.0 / (pow(6.,3./2.)*LAL_PI);
  REAL8 t_corr = gsl_spline_eval_deriv(spline_phi_TF2, Mf_ISCO_Schwrazschild, acc_phi_TF2) / (2*LAL_PI);

  // Now correct phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i] - fRef_geom;
    int j = i + offset; // shift index for frequency series if needed
    double phase_factor = -2*LAL_PI * f * t_corr;
    COMPLEX16 t_factor = (cos(phase_factor) + I*sin(phase_factor));//cexp(I*phase_factor);
    pdata[j] *= t_factor;
    cdata[j] *= t_factor;
  }

  XLALDestroyREAL8Sequence(freqs);
  gsl_interp_accel_free(acc_phi_TF2);
  gsl_spline_free(spline_phi_TF2);
  gsl_interp_accel_free(acc_amp_TF2);
  gsl_spline_free(spline_amp_TF2);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRv4TSurrogate_c
 *
 * \author Michael Puerrer, Ben Lackey
 *
 * \brief C code for SEOBNRv4T surrogate model
 * See arXiv:1812.08643
 *
 * This is a frequency domain model that approximates the time domain TEOBv4 model.
 *
 * The binary data HDF5 file (SEOBNRv4T_surrogate_v1.0.0.hdf5)
 * will be available at on LIGO clusters in /home/cbc/ and can be downloaded from
 * https://git.ligo.org/lscsoft/lalsuite-extra/blob/master/data/lalsimulation/SEOBNRv4T_surrogate_v1.0.0.hdf5.
 * Make sure the files are in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the surrogate has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv4T model.
 *
 * @note Parameter ranges:
 *   1 <= q <= 3
 *   -0.5 <= chi1 <= 0.5
 *   -0.5 <= chi2 <= 0.5
 *   0 <= lambda1 <= 5000
 *   0 <= lambda2 <= 5000
 *
 *  Aligned component spins chi1, chi2.
 *  Tidal deformabilities of neutron stars lambda1, lambda2.
 *  Mass-ratio q = m1/m2
 *
 * @{
 */


/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv4T_surrogate model.
 *
 * XLALSimIMRSEOBNRv4TSurrogate() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv4TSurrogateFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv4TSurrogateFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv4TSurrogateFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin 2 */
  REAL8 lambda1,                                /**< Dimensionless tidal deformability 1 */
  REAL8 lambda2,                                /**< Dimensionless tidal deformability 2 */
  SEOBNRv4TSurrogate_spline_order spline_order) /**< Spline order in frequency */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    double lambda1temp = lambda1;
    m1SI = m2SI;
    chi1 = chi2;
    lambda1 = lambda2;
    m2SI = m1temp;
    chi2 = chi1temp;
    lambda2 = lambda1temp;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);  /* Symmetric mass-ratio */
  double Mtot_sec = Mtot * LAL_MTSUN_SI;     /* Total mass in seconds */

  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&Surrogate_is_initialized, Surrogate_Init_LALDATA);
#else
  Surrogate_Init_LALDATA();
#endif

  if(!Surrogate_IsSetup()) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up Surrogate data - check your $LAL_DATA_PATH\n");
  }

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SurrogateCore(hptilde, hctilde, phiRef, fRef, distance,
                                inclination, Mtot_sec, eta, chi1, chi2, 
                                lambda1, lambda2, freqs, 0, spline_order);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv4T_surrogate model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv4TSurrogate(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin 2 */
  REAL8 lambda1,                                /**< Dimensionless tidal deformability 1 */
  REAL8 lambda2,                                /**< Dimensionless tidal deformability 2 */
  SEOBNRv4TSurrogate_spline_order spline_order) /**< Spline order in frequency */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    double lambda1temp = lambda1;
    m1SI = m2SI;
    chi1 = chi2;
    lambda1 = lambda2;
    m2SI = m1temp;
    chi2 = chi1temp;
    lambda2 = lambda1temp;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  double Mtot_sec = Mtot * LAL_MTSUN_SI;       /* Total mass in seconds */

  if(fRef==0.0)
    fRef=fLow;

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&Surrogate_is_initialized, Surrogate_Init_LALDATA);
#else
  Surrogate_Init_LALDATA();
#endif

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SurrogateCore(hptilde, hctilde, phiRef, fRef, distance,
                                inclination, Mtot_sec, eta, chi1, chi2, 
                                lambda1, lambda2, freqs, deltaF, spline_order);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */


/** Setup Surrogate model using data files installed in $LAL_DATA_PATH
 */
UNUSED static void Surrogate_Init_LALDATA(void)
{
  if (Surrogate_IsSetup()) return;

  // Expect ROM datafile in a directory listed in LAL_DATA_PATH,
#ifdef LAL_HDF5_ENABLED
#define datafile SurDataHDF5
  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = Surrogate_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find Surrogate data files in $LAL_DATA_PATH\n");
#else
  XLAL_ERROR_VOID(XLAL_EFAILED, "Surrogate requires HDF5 support which is not enabled\n");
#endif
}
