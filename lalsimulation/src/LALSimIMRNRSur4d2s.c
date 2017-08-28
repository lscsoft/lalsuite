/*
 *  Copyright (C) 2017 Jonathan Blackman
 *  NRSur4d2s_FDROM model, Reduced Order Model for NRSur4d2s.
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_complex_math.h>
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
#include <lal/SphericalHarmonics.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


/******* Surrogate model parameter space ********/
static const double Q_MIN = 1.;
static const double Q_MAX = 2.;
static const double CHI1_MAG_MAX = 0.8;
static const double CHI2_MIN = -0.8;
static const double CHI2_MAX = 0.8;

static const char NRSUR4D2S_DATAFILE[] = "NRSur4d2s_FDROM.hdf5";

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t NRSur4d2s_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

struct tagNRSurrogateData_submodel // A submodel is one (l, m) mode.
{
  gsl_vector* cvec_re;      // Flattened real spline coefficients
  gsl_vector* cvec_im;      // Flattened imag spline coefficients
  gsl_matrix *EI_re;        // Empirical interpolation matrix (real part)
  gsl_matrix *EI_im;        // Empirical interpolation matrix (imag part)
};
typedef struct tagNRSurrogateData_submodel NRSurrogateData_submodel;

typedef struct tagSplineData5d
{
  gsl_bspline_workspace *bw_x0; // Spline workspace for the first parameter
  gsl_bspline_workspace *bw_x1; // etc
  gsl_bspline_workspace *bw_x2;
  gsl_bspline_workspace *bw_x3;
  gsl_bspline_workspace *bw_x4;
} SplineData5d;

struct tagNRSurrogateData
{
  UINT4 setup;
  int n_nodes;  // Number of empirical nodes
  int n_freqs;  // Number of frequency samples in empirical interpolant output
  int *nc;    // Number of spline coefficients in each parameter dimension
  double *qvec;       //q points in tensor grid
  double *chi1vec;     //chi1 points in tensor grid
  double *chi1thetavec;
  double *chi1phivec;
  double *chi2vec;
  NRSurrogateData_submodel* mode_2m2;
  NRSurrogateData_submodel* mode_2m1;
  NRSurrogateData_submodel* mode_2p0;
  NRSurrogateData_submodel* mode_2p1;
  NRSurrogateData_submodel* mode_2p2;
  NRSurrogateData_submodel* mode_3m3;
  NRSurrogateData_submodel* mode_3m2;
  NRSurrogateData_submodel* mode_3m1;
  NRSurrogateData_submodel* mode_3p0;
  NRSurrogateData_submodel* mode_3p1;
  NRSurrogateData_submodel* mode_3p2;
  NRSurrogateData_submodel* mode_3p3;

  gsl_vector* fEI;                      // Frequency samples of empirical interpolant output
  SplineData5d* splinedata;               // Data associated with the tensor product grid in parameter space
};
typedef struct tagNRSurrogateData NRSurrogateData;

static NRSurrogateData __lalsim_NRSurrogate_data;

/**************** Internal functions **********************/

static void err_handler(const char *reason, const char *file, int line, int gsl_errno);
static void NRSurrogate_Init_LALDATA(void);
static int NRSurrogate_Init(const char dir[]);
static bool NRSurrogate_IsSetup(void);

static int NRSurrogateData_Init(NRSurrogateData *data, const char dir[]);
static void NRSurrogateData_Cleanup(NRSurrogateData *data);

// Helper function to interpolate all empirical nodes in 5d for one waveform mode
static int TP_Spline_interpolation_5d(
  REAL8 x0,                 // Coordinates onto which we interpolate
  REAL8 x1,
  REAL8 x2,
  REAL8 x3,
  REAL8 x4,
  int n_nodes,              // Input: number of empirical nodes
  int *nc,                // Input: number of spline coefs per dimension
  gsl_vector *cvec_re,      // Input: flattened real spline coefficients
  gsl_vector *cvec_im,      // Input: flattened imag spline coefficients
  SplineData5d *splinedata, // Input: Data associated with the tensor grid in parameter space
  gsl_vector *nodes_re,     // Output: interpolated real empirical nodes
  gsl_vector *nodes_im      // Output: interpolated imag empirical nodes
);

static int NRSurrogateData_Init_submodel(
  NRSurrogateData_submodel **submodel,
  LALH5File *file,
  const int i_mode,
  int n_nodes,
  int n_freqs,
  int *nc
);

static void NRSurrogateData_Cleanup_submodel(NRSurrogateData_submodel *submodel);

/**
 * Evaluates one mode of the surrogate model and adds it to the waveform
 * Real and imaginary parts used due to lack of gsl_vector_complex_conjugate() function
 * Note that here the output waveform is a single complex function containing positive
 * and negative frequency content, from which we may obtain hplus and hcross.
 * The output is evaluated at the empirical interpolation output samples.
 * Output frequencies and waveform are dimensionless.
*/
static void AddModeContribution(
  double q,             // Parameters at which we evaluate the surrogate mode
  double chi1mag,
  double chi1theta,
  double chi1phi,
  double chi2z,
  double spheretheta,   // The polar angle of the direction of GW propagation from the source
  double spherephi,     // The azimuthal angle of the direction of GW propagation from the source
  int swsh_L,           // The mode L
  int swsh_m,           // The mode m
  int n_nodes,
  int *nc,
  gsl_vector *h_re,     // Gets modified - the real part of the waveform
  gsl_vector *h_im,     // Gets modified - the imag part of the waveform
  NRSurrogateData_submodel *modeData, // The surrogate data for this mode
  SplineData5d *splinedata, // Data associated with the tensor product grid in param space
  gsl_vector *nodes_re, // Intermediate result - the real part of the empirical nodes for this mode
  gsl_vector *nodes_im, // Intermediate result - the imag part of the empirical nodes for this mode
  gsl_vector *mode_re,  // Intermediate result - the real part of the evaluated mode
  gsl_vector *mode_im  // Intermediate result - the imag part of the evaluated mode
);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int NRSurrogateCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi1mag,
  double chi1theta,
  double chi1phi,
  double chi2z,
  const REAL8Sequence *freqs /* Frequency points at which to evaluate the waveform (Hz) */
);

static void SplineData5d_Destroy(SplineData5d *splinedata);
static void SplineData5d_Init(
  SplineData5d **splinedata,
  int *nc,   // Number of spline coefs per dimension
  double *x1,  // Parameter grid values
  double *x2,
  double *x3,
  double *x4,
  double *x5
);

static int load_data_sub(
  const int i_mode,
  LALH5File *file,
  gsl_vector *cvec_re,
  gsl_vector *cvec_im,
  gsl_matrix *EI_re,
  gsl_matrix *EI_im
);

// Function which sums over all contributing basis splines
static REAL8 Interpolate_Coefficent_Tensor_5d(
  gsl_vector *v,
  gsl_vector *B0,               // 4 non-zero basis splines
  gsl_vector *B1,
  gsl_vector *B2,
  gsl_vector *B3,
  gsl_vector *B4,
  size_t is0,                   // Index of first of 4 coefficients to evaluate
  size_t is1,
  size_t is2,
  size_t is3,
  size_t is4,
  int *nc                     // Number of spline coefficients in each dimension
);

/********************* Definitions begin here ********************/

static void err_handler(const char *reason, const char *file, int line, int gsl_errno) {
  XLALPrintError("gsl: %s:%d: %s - %d\n", file, line, reason, gsl_errno);
}

// Function which sums over all contributing basis splines
static REAL8 Interpolate_Coefficent_Tensor_5d(
  gsl_vector *v,
  gsl_vector *B0,               // 4 non-zero basis splines
  gsl_vector *B1,
  gsl_vector *B2,
  gsl_vector *B3,
  gsl_vector *B4,
  size_t is0,                   // Index of first of 4 coefficients to evaluate
  size_t is1,
  size_t is2,
  size_t is3,
  size_t is4,
  int *nc                     // Number of coefficients in each dimension
) {
  // Compute coefficient at desired parameters from
  // C = c_{i0, i1, i2, i3, i4} * B1[i0] * B2[i1] * ...
  // while summing over indices i0,...,i4 where the B-splines are nonzero.
  double sum = 0;

  for (int i0=0; i0<4; i0++) {
    int ii0 = is0 + i0;
    for (int i1=0; i1<4; i1++) {
      int ii1 = is1 + i1;
      for (int i2=0; i2<4; i2++) {
        int ii2 = is2 + i2;
        for (int i3=0; i3<4; i3++) {
          int ii3 = is3 + i3;
          for (int i4=0; i4<4; i4++) {
            int ii4 = is4 + i4;
            double c = gsl_vector_get(v, (((ii0*nc[1] + ii1)*nc[2] + ii2)*nc[3] + ii3)*nc[4] + ii4);
            sum += c * gsl_vector_get(B0, i0) * gsl_vector_get(B1, i1) * gsl_vector_get(B2, i2) *
                            gsl_vector_get(B3, i3) * gsl_vector_get(B4, i4);
          }
        }
      }
    }
  }

  return sum;
}

/** Setup NRSurrogate model using data files installed in dir
 */
static int NRSurrogate_Init(const char dir[]) {
  if(__lalsim_NRSurrogate_data.setup) {
    XLALPrintError("Error: NRSurrogate was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  NRSurrogateData_Init(&__lalsim_NRSurrogate_data, dir);

  if(__lalsim_NRSurrogate_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the NRSurrogate model has been initialised */
static bool NRSurrogate_IsSetup(void) {
  if(__lalsim_NRSurrogate_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients for submodel 1
static int load_data_sub(
  const int i_mode,
  LALH5File *file,
  gsl_vector *cvec_re,
  gsl_vector *cvec_im,
  gsl_matrix *EI_re,
  gsl_matrix *EI_im
) {
  // Load H5 data sets for spline coefficients and empirical interpolation matrix

  int ret = XLAL_SUCCESS;
  size_t name_length = 22;
  if (i_mode > 9) name_length = 23;
  char *dataset_name = malloc(name_length);

  snprintf(dataset_name, name_length, "NRSurrogate_cre_%d", i_mode);
  ret |= ReadHDF5RealVectorDataset(file, dataset_name, &cvec_re);

  snprintf(dataset_name, name_length, "NRSurrogate_cim_%d", i_mode);
  ret |= ReadHDF5RealVectorDataset(file, dataset_name, &cvec_im);

  snprintf(dataset_name, name_length, "NRSurrogate_EIr_%d.dat", i_mode);
  ret |= ReadHDF5RealMatrixDataset(file, dataset_name, &EI_re);

  snprintf(dataset_name, name_length, "NRSurrogate_EIi_%d.dat", i_mode);
  ret |= ReadHDF5RealMatrixDataset(file, dataset_name, &EI_im);

  free(dataset_name);
  return(ret);
}

// Setup B-spline basis functions for given points
static void SplineData5d_Init(
  SplineData5d **splinedata,
  int *nc,  // Number of coefficients per dimension
  double *x0, // Parameter grid values
  double *x1,
  double *x2,
  double *x3,
  double *x4
) {
  if(!splinedata) exit(1);
  if(*splinedata) SplineData5d_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData5d));

  // Set up B-spline basis for desired knots
  const size_t nbreak_0 = nc[0]-2;  // nbreak = ncoefs - 2 for cubic splines
  const size_t nbreak_1 = nc[1]-2;
  const size_t nbreak_2 = nc[2]-2;
  const size_t nbreak_3 = nc[3]-2;
  const size_t nbreak_4 = nc[4]-2;

  // Allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bw_x0 = gsl_bspline_alloc(4, nbreak_0);
  gsl_bspline_workspace *bw_x1 = gsl_bspline_alloc(4, nbreak_1);
  gsl_bspline_workspace *bw_x2 = gsl_bspline_alloc(4, nbreak_2);
  gsl_bspline_workspace *bw_x3 = gsl_bspline_alloc(4, nbreak_3);
  gsl_bspline_workspace *bw_x4 = gsl_bspline_alloc(4, nbreak_4);

  // Set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_0 = gsl_vector_alloc(nbreak_0);
  gsl_vector *breakpts_1 = gsl_vector_alloc(nbreak_1);
  gsl_vector *breakpts_2 = gsl_vector_alloc(nbreak_2);
  gsl_vector *breakpts_3 = gsl_vector_alloc(nbreak_3);
  gsl_vector *breakpts_4 = gsl_vector_alloc(nbreak_4);
  for (UINT4 i=0; i<nbreak_0; i++)
    gsl_vector_set(breakpts_0, i, x0[i]);
  for (UINT4 i=0; i<nbreak_1; i++)
    gsl_vector_set(breakpts_1, i, x1[i]);
  for (UINT4 i=0; i<nbreak_2; i++)
    gsl_vector_set(breakpts_2, i, x2[i]);
  for (UINT4 i=0; i<nbreak_3; i++)
    gsl_vector_set(breakpts_3, i, x3[i]);
  for (UINT4 i=0; i<nbreak_4; i++)
    gsl_vector_set(breakpts_4, i, x4[i]);

  // Calculate knots
  gsl_bspline_knots(breakpts_0, bw_x0);
  gsl_bspline_knots(breakpts_1, bw_x1);
  gsl_bspline_knots(breakpts_2, bw_x2);
  gsl_bspline_knots(breakpts_3, bw_x3);
  gsl_bspline_knots(breakpts_4, bw_x4);

  gsl_vector_free(breakpts_0);
  gsl_vector_free(breakpts_1);
  gsl_vector_free(breakpts_2);
  gsl_vector_free(breakpts_3);
  gsl_vector_free(breakpts_4);

  (*splinedata)->bw_x0=bw_x0;
  (*splinedata)->bw_x1=bw_x1;
  (*splinedata)->bw_x2=bw_x2;
  (*splinedata)->bw_x3=bw_x3;
  (*splinedata)->bw_x4=bw_x4;
}

static void SplineData5d_Destroy(SplineData5d *splinedata)
{
  if(!splinedata) return;
  if(splinedata->bw_x0) gsl_bspline_free(splinedata->bw_x0);
  if(splinedata->bw_x1) gsl_bspline_free(splinedata->bw_x1);
  if(splinedata->bw_x2) gsl_bspline_free(splinedata->bw_x2);
  if(splinedata->bw_x3) gsl_bspline_free(splinedata->bw_x3);
  if(splinedata->bw_x4) gsl_bspline_free(splinedata->bw_x4);

  XLALFree(splinedata);
}

// Interpolate empirical nodes over the parameter space.
// The multi-dimensional interpolation is carried out via a tensor product decomposition.
static int TP_Spline_interpolation_5d(
  REAL8 x0,                 // Coordinates onto which we interpolate
  REAL8 x1,
  REAL8 x2,
  REAL8 x3,
  REAL8 x4,
  int n_nodes,              // Input: number of empirical nodes
  int *nc,                // Input: number of coefficients in each dimension
  gsl_vector *cvec_re,      // Input: flattened real spline coefficients
  gsl_vector *cvec_im,      // Input: flattened imag spline coefficients
  SplineData5d *splinedata, // Input: Data associated with the tensor grid in parameter space
  gsl_vector *nodes_re,     // Output: interpolated real empirical nodes
  gsl_vector *nodes_im      // Output: interpolated imag empirical nodes
) {

  // Store nonzero cubic (order k=4) B-spline basis functions in all directions.
  gsl_vector *B0 = gsl_vector_calloc(4);
  gsl_vector *B1 = gsl_vector_calloc(4);
  gsl_vector *B2 = gsl_vector_calloc(4);
  gsl_vector *B3 = gsl_vector_calloc(4);
  gsl_vector *B4 = gsl_vector_calloc(4);

  size_t is0, is1, is2, is3, is4; // First non-zero basis spline
  size_t ie0, ie1, ie2, ie3, ie4; // Last non-zero basis spline

  // Evaluate all potentially nonzero cubic B-spline basis functions and store them in the B vectors.
  // Since the B-splines are of compact support we only need to store a small
  // number of basis functions to avoid computing terms that would be zero anyway.
  // https://www.gnu.org/software/gsl/manual/html_node/Overview-of-B_002dsplines.html#Overview-of-B_002dsplines
  gsl_bspline_eval_nonzero(x0,  B0, &is0, &ie0, splinedata->bw_x0);
  gsl_bspline_eval_nonzero(x1,  B1, &is1, &ie1, splinedata->bw_x1);
  gsl_bspline_eval_nonzero(x2,  B2, &is2, &ie2, splinedata->bw_x2);
  gsl_bspline_eval_nonzero(x3,  B3, &is3, &ie3, splinedata->bw_x3);
  gsl_bspline_eval_nonzero(x4,  B4, &is4, &ie4, splinedata->bw_x4);

  int N = nc[0]*nc[1]*nc[2]*nc[3]*nc[4];  // Size of the data matrix for one empirical node

  // Evaluate the TP spline for all empirical nodes
  for (int k=0; k<n_nodes; k++) {
    // Pick out the coefficient matrix corresponding to the k-th empirical node.
    gsl_vector v_re = gsl_vector_subvector(cvec_re, k*N, N).vector;
    gsl_vector v_im = gsl_vector_subvector(cvec_im, k*N, N).vector;

    REAL8 node_re = Interpolate_Coefficent_Tensor_5d(&v_re,
            B0, B1, B2, B3, B4,
            is0, is1, is2, is3, is4, nc);

    REAL8 node_im = Interpolate_Coefficent_Tensor_5d(&v_im,
            B0, B1, B2, B3, B4,
            is0, is1, is2, is3, is4, nc);
    gsl_vector_set(nodes_re, k, node_re);
    gsl_vector_set(nodes_im, k, node_im);
  }

  gsl_vector_free(B0);
  gsl_vector_free(B1);
  gsl_vector_free(B2);
  gsl_vector_free(B3);
  gsl_vector_free(B4);

  return(0);
}

/* Set up a new ROM submodel, using data contained in the NRSur4d2s.h5 h5 file */
static int NRSurrogateData_Init_submodel(
  NRSurrogateData_submodel **submodel,
  LALH5File *file,
  const int i_mode,
  int n_nodes,
  int n_freqs,
  int *nc
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(NRSurrogateData_submodel));
  else
    NRSurrogateData_Cleanup_submodel(*submodel);

  int N = nc[0]*nc[1]*nc[2]*nc[3]*nc[4]; // Total number of spline coefficients for one empirical node

  // Initalize actual ROM data
  (*submodel)->cvec_re = gsl_vector_calloc(N*n_nodes);
  if ((*submodel)->cvec_re == NULL)
    XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_calloc(%d) failed", N*n_nodes);
  (*submodel)->cvec_im = gsl_vector_calloc(N*n_nodes);
  if ((*submodel)->cvec_im == NULL)
    XLAL_ERROR(XLAL_ENOMEM, "gsl_vector_calloc(%d) failed", N*n_nodes);
  (*submodel)->EI_re = gsl_matrix_calloc(n_nodes, n_freqs);
  if ((*submodel)->EI_re == NULL)
    XLAL_ERROR(XLAL_ENOMEM, "gsl_matrix_calloc(%d, %d) failed", n_nodes, n_freqs);
  (*submodel)->EI_im = gsl_matrix_calloc(n_nodes, n_freqs);
  if ((*submodel)->EI_im == NULL)
    XLAL_ERROR(XLAL_ENOMEM, "gsl_matrix_calloc(%d, %d) failed", n_nodes, n_freqs);

  // Load ROM data for this submodel
  printf("Loading mode %d...\n", i_mode);
  ret=load_data_sub(i_mode, file, (*submodel)->cvec_re, (*submodel)->cvec_im, (*submodel)->EI_re, (*submodel)->EI_im);

  return ret;
}

/* Deallocate contents of the given NRSurrogateData_submodel structure */
static void NRSurrogateData_Cleanup_submodel(NRSurrogateData_submodel *submodel) {
  if(submodel->cvec_re) gsl_vector_free(submodel->cvec_re);
  if(submodel->cvec_im) gsl_vector_free(submodel->cvec_im);
  if(submodel->EI_re) gsl_matrix_free(submodel->EI_re);
  if(submodel->EI_im) gsl_matrix_free(submodel->EI_im);
}

/* Set up a new NRSurrogate model, using data contained in dir */
int NRSurrogateData_Init(NRSurrogateData *data, const char dir[]) {
  int ret = XLAL_FAILURE;

  if(data->setup) {
    XLALPrintError("WARNING: You tried to setup the NRSurrogate model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  printf("Loading NRSur4d2s_FDROM data...\n");

  size_t size = strlen(dir) + strlen(NRSUR4D2S_DATAFILE) + 2;
  char *filename = XLALMalloc(size);
  snprintf(filename, size, "%s/%s", dir, NRSUR4D2S_DATAFILE);

  LALH5File *file = XLALH5FileOpen(filename, "r");
  LALH5Dataset *dset;

  // Parse surrogate model parameters
  int n_nodes, n_freqs, n_q, n_chi1, n_chi1theta, n_chi1phi, n_chi2;
  INT8 tmp_param;
  dset = XLALH5DatasetRead(file, "n_nodes");
  ret = XLALH5DatasetQueryData(&tmp_param, dset);
  n_nodes = (int) tmp_param;
  data->n_nodes = n_nodes;
  dset = XLALH5DatasetRead(file, "n_freqs");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_freqs = (int) tmp_param;
  data->n_freqs = n_freqs;
  dset = XLALH5DatasetRead(file, "n_q");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_q = (int) tmp_param;
  dset = XLALH5DatasetRead(file, "n_chi1");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_chi1 = (int) tmp_param;
  dset = XLALH5DatasetRead(file, "n_chi1theta");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_chi1theta = (int) tmp_param;
  dset = XLALH5DatasetRead(file, "n_chi1phi");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_chi1phi = (int) tmp_param;
  dset = XLALH5DatasetRead(file, "n_chi2");
  ret |= XLALH5DatasetQueryData(&tmp_param, dset);
  n_chi2 = (int) tmp_param;

  data->nc = (int *) malloc(5*sizeof(double));
  data->nc[0] = n_q+2;
  data->nc[1] = n_chi1+2;
  data->nc[2] = n_chi1theta+2;
  data->nc[3] = n_chi1phi+2;
  data->nc[4] = n_chi2+2;

  // Setup parameter vectors
  int i;
  double delta;

  // q
  data->qvec = (double *)malloc(n_q*sizeof(double));
  delta = (Q_MAX - Q_MIN)/(n_q - 1);
  for(i=0; i<n_q; i++) data->qvec[i] = Q_MIN + i*delta;

  // chi1
  data->chi1vec = (double *)malloc(n_chi1*sizeof(double));
  delta = CHI1_MAG_MAX/(n_chi1 - 1);
  for(i=0; i<n_chi1; i++) data->chi1vec[i] = i*delta;

  // chi1theta
  data->chi1thetavec = (double *)malloc(n_chi1theta*sizeof(double));
  delta = LAL_PI/(n_chi1theta - 1);
  for(i=0; i<n_chi1theta; i++) data->chi1thetavec[i] = i*delta;

  // chi1phi
  data->chi1phivec = (double *)malloc(n_chi1phi*sizeof(double));
  delta = 2.0*LAL_PI/(n_chi1phi - 1);
  for(i=0; i<n_chi1phi; i++) data->chi1phivec[i] = i*delta;

  // chi2
  data->chi2vec = (double *)malloc(n_chi2*sizeof(double));
  delta = (CHI2_MAX - CHI2_MIN)/(n_chi2 - 1);
  for(i=0; i<n_chi2; i++) data->chi2vec[i] = CHI2_MIN + i*delta;

  // Setup modes
  gsl_set_error_handler(&err_handler);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_2m2, file, 0, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 0 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_2m1, file, 1, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 1 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_2p0, file, 2, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 2 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_2p1, file, 3, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 3 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_2p2, file, 4, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 4 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3m3, file, 5, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 5 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3m2, file, 6, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 6 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3m1, file, 7, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 7 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3p0, file, 8, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 8 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3p1, file, 9, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 9 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3p2, file, 10, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 10 loaded sucessfully.\n", __func__);

  ret |= NRSurrogateData_Init_submodel(&(data)->mode_3p3, file, 11, n_nodes, n_freqs, data->nc);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : Mode 11 loaded sucessfully.\n", __func__);


  if (ret==XLAL_SUCCESS) printf("Successfully loaded all modes!\n");

  /* Load frequency samples */
  data->fEI = gsl_vector_calloc(n_freqs);
  ret |= ReadHDF5RealVectorDataset(file, "freqs", &(data->fEI));

  XLALH5FileClose(file);

  /* setup splinedata */
  SplineData5d *splinedata=NULL;
  SplineData5d_Init(&splinedata, data->nc, data->qvec, data->chi1vec, data->chi1thetavec, data->chi1phivec, data->chi2vec);
  data->splinedata = splinedata;

  if(XLAL_SUCCESS==ret)
    data->setup=1;
  else
    NRSurrogateData_Cleanup(data);

  printf("Done setting up\n");
  return (ret);
}

/* Deallocate contents of the given NRSurrogate structure */
static void NRSurrogateData_Cleanup(NRSurrogateData *data) {
  NRSurrogateData_Cleanup_submodel((data)->mode_2m2);
  NRSurrogateData_Cleanup_submodel((data)->mode_2m1);
  NRSurrogateData_Cleanup_submodel((data)->mode_2p0);
  NRSurrogateData_Cleanup_submodel((data)->mode_2p1);
  NRSurrogateData_Cleanup_submodel((data)->mode_2p2);
  NRSurrogateData_Cleanup_submodel((data)->mode_3m3);
  NRSurrogateData_Cleanup_submodel((data)->mode_3m2);
  NRSurrogateData_Cleanup_submodel((data)->mode_3m1);
  NRSurrogateData_Cleanup_submodel((data)->mode_3p0);
  NRSurrogateData_Cleanup_submodel((data)->mode_3p1);
  NRSurrogateData_Cleanup_submodel((data)->mode_3p2);
  NRSurrogateData_Cleanup_submodel((data)->mode_3p3);
  XLALFree((data)->mode_2m2);
  XLALFree((data)->mode_2m1);
  XLALFree((data)->mode_2p0);
  XLALFree((data)->mode_2p1);
  XLALFree((data)->mode_2p2);
  XLALFree((data)->mode_3m3);
  XLALFree((data)->mode_3m2);
  XLALFree((data)->mode_3m1);
  XLALFree((data)->mode_3p0);
  XLALFree((data)->mode_3p1);
  XLALFree((data)->mode_3p2);
  XLALFree((data)->mode_3p3);
  (data)->mode_2m2 = NULL;
  (data)->mode_2m1 = NULL;
  (data)->mode_2p0 = NULL;
  (data)->mode_2p1 = NULL;
  (data)->mode_2p2 = NULL;
  (data)->mode_3p3 = NULL;
  (data)->mode_3p2 = NULL;
  (data)->mode_3p1 = NULL;
  (data)->mode_3p0 = NULL;
  (data)->mode_3m1 = NULL;
  (data)->mode_3m2 = NULL;
  (data)->mode_3m3 = NULL;
  if(data->nc) free(data->nc);
  if(data->fEI) gsl_vector_free(data->fEI);
  if(data->qvec) free(data->qvec);
  if(data->chi1vec) free(data->chi1vec);
  if(data->chi1thetavec) free(data->chi1thetavec);
  if(data->chi1phivec) free(data->chi1thetavec);
  if(data->chi2vec) free(data->chi2vec);
  if(data->splinedata) SplineData5d_Destroy(data->splinedata);
  data->setup=0;
}

//Add one spherical harmonic mode to the waveform h, weighted by the SWSH coefficient
static void AddModeContribution(
  double q,             // Parameters at which we evaluate the surrogate mode
  double chi1mag,
  double chi1theta,
  double chi1phi,
  double chi2z,
  double spheretheta,   // The polar angle of the direction of GW propagation from the source
  double spherephi,     // The azimuthal angle of the direction of GW propagation from the source
  int swsh_L,           // The mode L
  int swsh_m,           // The mode m
  int n_nodes,          // Input: the number of empirical nodes
  int *nc,              // Input: the number of spline coefs per dimension
  gsl_vector *h_re,     // Gets modified - the real part of the waveform
  gsl_vector *h_im,     // Gets modified - the imag part of the waveform
  NRSurrogateData_submodel *modeData, // The surrogate data for this mode
  SplineData5d *splinedata, // Data associated with the tensor product grid in param space
  gsl_vector *nodes_re, // Intermediate result - the real part of the empirical nodes for this mode
  gsl_vector *nodes_im, // Intermediate result - the imag part of the empirical nodes for this mode
  gsl_vector *mode_re,  // Intermediate result - the real part of the evaluated mode
  gsl_vector *mode_im  // Intermediate result - the imag part of the evaluated mode
){
  //Interpolate nodes
  TP_Spline_interpolation_5d(q, chi1mag, chi1theta, chi1phi, chi2z, n_nodes, nc, modeData->cvec_re, modeData->cvec_im,
                             splinedata, nodes_re, nodes_im);

  // Evaluate empirical interpolant
  gsl_blas_dgemv(CblasTrans, 1.0, modeData->EI_re, nodes_re, 0., mode_re); // mode_re = EI_re*nodes_re
  gsl_blas_dgemv(CblasTrans, -1., modeData->EI_im, nodes_im, 1., mode_re); // mode_re -= EI_im*nodes_im
  gsl_blas_dgemv(CblasTrans, 1.0, modeData->EI_im, nodes_re, 0., mode_im); // mode_im = EI_im*nodes_re
  gsl_blas_dgemv(CblasTrans, 1.0, modeData->EI_re, nodes_im, 1., mode_im); // mode_im += EI_re*nodes_im

  // Add (spin-weighted spherical harmonic coefficient)*mode to h
  COMPLEX16 swsh_coef = XLALSpinWeightedSphericalHarmonic(spheretheta, spherephi, -2, swsh_L, swsh_m);
  double c_re = creal(swsh_coef);
  double c_im = cimag(swsh_coef);
  gsl_blas_daxpy(c_re, mode_re, h_re);      // h_re += c_re*mode_re
  gsl_blas_daxpy(-1*c_im, mode_im, h_re);   // h_re += -c_im*mode_im
  gsl_blas_daxpy(c_re, mode_im, h_im);      // h_im += c_re*mode_im
  gsl_blas_daxpy(c_im, mode_re, h_im);      // h_im += c_im*mode_re
}

static int NRSurrogateCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi1mag,
  double chi1theta,
  double chi1phi,
  double chi2z,
  const REAL8Sequence *freqs /* Dimensionless frequency points at which to evaluate the waveform */
)
{
  // Initialize waveform at empirical interpolant output frequencies (contains negative and positive frequency content)
  int n_freqs = (&__lalsim_NRSurrogate_data)->n_freqs;
  gsl_vector *h_full_re = gsl_vector_calloc(n_freqs);
  gsl_vector *h_full_im = gsl_vector_calloc(n_freqs);

  // Allocate memory for intermediate results
  int n_nodes = (&__lalsim_NRSurrogate_data)->n_nodes;
  gsl_vector *nodes_re = gsl_vector_calloc(n_nodes);
  gsl_vector *nodes_im = gsl_vector_calloc(n_nodes);
  gsl_vector *mode_re = gsl_vector_calloc(n_freqs);
  gsl_vector *mode_im = gsl_vector_calloc(n_freqs);

  // Add up all modes
  int *nc = (&__lalsim_NRSurrogate_data)->nc;
  double phi_sphere = 0.5*LAL_PI - phiRef;
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 2, -2, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_2m2, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 2, -1, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_2m1, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 2, 0, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_2p0, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 2, 1, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_2p1, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 2, 2, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_2p2, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, -3, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3m3, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, -2, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3m2, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, -1, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3m1, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, 0, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3p0, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, 1, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3p1, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, 2, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3p2, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);
  AddModeContribution(q, chi1mag, chi1theta, chi1phi, chi2z, inclination, phi_sphere, 3, 3, n_nodes, nc, h_full_re, h_full_im,
                      (&__lalsim_NRSurrogate_data)->mode_3p3, (&__lalsim_NRSurrogate_data)->splinedata,
                      nodes_re, nodes_im, mode_re, mode_im);

  gsl_vector_free(mode_re);
  gsl_vector_free(mode_im);
  gsl_vector_free(nodes_re);
  gsl_vector_free(nodes_im);

  // We need h(f) for f > 0 and h(-f) for f < 0
  int n = n_freqs/2;
  int i0 = n_freqs - n_freqs/2;
  gsl_vector h_negf_re = gsl_vector_subvector(h_full_re, 0, n).vector;
  gsl_vector h_negf_im = gsl_vector_subvector(h_full_im, 0, n).vector;
  gsl_vector h_posf_re = gsl_vector_subvector(h_full_re, i0, n).vector;
  gsl_vector h_posf_im = gsl_vector_subvector(h_full_im, i0, n).vector;
  gsl_vector_reverse(&h_negf_re);
  gsl_vector_reverse(&h_negf_im);

  // TODO: There's probably a more memory efficient way to do do the previous step and this
  // hp = 0.5*(h(f) + h(-f)*)
  // hc = i*0.5*(h(f) - h(-f)*)
  gsl_vector *hp_re = gsl_vector_alloc(n);
  gsl_vector *hp_im = gsl_vector_alloc(n);
  gsl_vector *hc_re = gsl_vector_alloc(n);
  gsl_vector *hc_im = gsl_vector_alloc(n);

  gsl_vector_memcpy(hp_re, &h_posf_re); //hp_re = h_posf_re
  gsl_blas_daxpy(1.0, &h_negf_re, hp_re); //hp_re = h_posf_re + h_negf_re
  gsl_blas_dscal(0.5, hp_re); //hp_re = 0.5*(h_posf_re + h_negf_re)

  gsl_vector_memcpy(hp_im, &h_posf_im);
  gsl_blas_daxpy(-1.0, &h_negf_im, hp_im);
  gsl_blas_dscal(0.5, hp_im); //hp_im = 0.5*(h_posf_im - h_negf_im)

  gsl_vector_memcpy(hc_im, &h_posf_re);
  gsl_blas_daxpy(-1.0, &h_negf_re, hc_im);
  gsl_blas_dscal(0.5, hc_im); //hc_im = 0.5*(h_posf_re - h_negf_re)

  gsl_vector_memcpy(hc_re, &h_posf_im);
  gsl_blas_daxpy(1.0, &h_negf_im, hc_re);
  gsl_blas_dscal(-0.5, hc_re); //hp_re = -0.5*(h_posf_im + h_negf_im)

  // Setup interpolation of hp, hc onto hptilde, hctilde at freqs
  gsl_vector fEI = gsl_vector_subvector((&__lalsim_NRSurrogate_data)->fEI, i0, n).vector;
  const double *model_freqs = gsl_vector_const_ptr(&fEI, 0);
  gsl_interp_accel *acc = gsl_interp_accel_alloc(); //All interpolations use the same grid -> use a single accelerator
  gsl_spline *spl_hp_re = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline *spl_hp_im = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline *spl_hc_re = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline *spl_hc_im = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(spl_hp_re, model_freqs, gsl_vector_const_ptr(hp_re, 0), n);
  gsl_spline_init(spl_hp_im, model_freqs, gsl_vector_const_ptr(hp_im, 0), n);
  gsl_spline_init(spl_hc_re, model_freqs, gsl_vector_const_ptr(hc_re, 0), n);
  gsl_spline_init(spl_hc_im, model_freqs, gsl_vector_const_ptr(hc_im, 0), n);

  // Interpolate and dimensionalize the waveform
  double a0 = Mtot_sec * Mtot_sec * LAL_C_SI / distance;
  LIGOTimeGPS tC = {0, 0};
  size_t npts = freqs->length;
  *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, 0, 1.0, &lalStrainUnit, npts);
  *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, 0, 1.0, &lalStrainUnit, npts);
  memset((*hptilde)->data->data, 0, npts * sizeof(COMPLEX16));
  memset((*hctilde)->data->data, 0, npts * sizeof(COMPLEX16));
  double f_min = model_freqs[0];
  double f_max = model_freqs[n-1];

  for (UINT4 i=0; i < npts; i++) {
    double f = freqs->data[i];
    if (f > f_max || f < f_min) continue;
    (*hptilde)->data->data[i] = a0*(gsl_spline_eval(spl_hp_re, f, acc) + I*gsl_spline_eval(spl_hp_im, f, acc));
    (*hctilde)->data->data[i] = a0*(gsl_spline_eval(spl_hc_re, f, acc) + I*gsl_spline_eval(spl_hc_im, f, acc));
  }


  // Keep these around for the vector views until we don't need h_negf_re etc.
  gsl_vector_free(h_full_re);
  gsl_vector_free(h_full_im);

  gsl_vector_free(hp_re);
  gsl_vector_free(hp_im);
  gsl_vector_free(hc_re);
  gsl_vector_free(hc_im);
  gsl_spline_free(spl_hp_re);
  gsl_spline_free(spl_hp_im);
  gsl_spline_free(spl_hc_re);
  gsl_spline_free(spl_hc_im);
  gsl_interp_accel_free(acc);

  return XLAL_SUCCESS;
}

int XLALSimNRSur4d2sFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz), need to be strictly monotonically increasing */
  REAL8 phiRef,                                 /**< Orbital phase at reference frequency*/
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 S1x,                                    /**< x-component of the dimensionless spin of companion 1 */
  REAL8 S1y,                                    /**< y-component of the dimensionless spin of companion 1 */
  REAL8 S1z,                                    /**< z-component of the dimensionless spin of companion 1 */
  REAL8 S2x,                                    /**< x-component of the dimensionless spin of companion 2 */
  REAL8 S2y,                                    /**< y-component of the dimensionless spin of companion 2 */
  REAL8 S2z                                     /**< z-component of the dimensionless spin of companion 2 */
) {
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double tmp = m1SI;
    m1SI = m2SI;
    m2SI = tmp;
    tmp = S1x;
    S1x = S2x;
    S2x = tmp;
    tmp = S1y;
    S1y = S2y;
    S2y = tmp;
    tmp = S1z;
    S1z = S2z;
    S2z = tmp;
    phiRef += LAL_PI;
  }

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */
  double q = mass1/mass2; /* mass ratio (>= 1)*/
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Get model parameter spins
  double chi1mag = sqrt(S1x*S1x + S1y*S1y + S1z*S1z);
  double chi1theta = 0;
  double chi1phi = 0;
  if (chi1mag > 0.) chi1theta = acos(S1z/chi1mag);
  if (fabs(S1y) + fabs(S1x) > 0.) chi1phi = atan2(S1y*cos(phiRef) - S1x*sin(phiRef), S1x*cos(phiRef) + S1y*sin(phiRef));
  if (chi1phi < 0.) chi1phi += 2*LAL_PI;
  double chi2z = S2z;

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur4d2s_is_initialized, NRSurrogate_Init_LALDATA);
#else
  NRSurrogate_Init_LALDATA();
#endif

  if(!NRSurrogate_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up NRSur4d2s data - check your $LAL_DATA_PATH\n");

  int retcode = NRSurrogateCore(hptilde, hctilde,
            phiRef, distance, inclination, Mtot_sec, q, chi1mag, chi1theta, chi1phi, chi2z, freqs);

  return(retcode);
}



/**
 * Compute waveform in LAL format for the NRSur4d2s_FDROM NR surrogate model.
 *
 *  Described in:
 *  https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.104023
 *  https://arxiv.org/abs/1701.00550
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency fHigh to the next power of 2 in
 * the size of the frequency series.
 */
int XLALSimNRSur4d2s(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at 4500M before peak amplitude.
                                                Note that usually this is TWICE the orb phase.*/
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 S1x,    /**< Spins are in the source frame, which is different from the alignment frame.
                See http://software.ligo.org/docs/lalsuite/lalsimulation/group__lalsimulation__inspiral.html
                In the alignment frame, the larger BH is on the positive x-axis. x-component of chi_1 */
  REAL8 S1y,                                    /**< y-component of the dimensionless spin of companion 1 */
  REAL8 S1z,                                    /**< z-component of the dimensionless spin of companion 1 */
  REAL8 S2x,                                    /**< x-component of the dimensionless spin of companion 2 */
  REAL8 S2y,                                    /**< y-component of the dimensionless spin of companion 2 */
  REAL8 S2z                                     /**< z-component of the dimensionless spin of companion 2 */
) {
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2.
    double tmp = m1SI;
    m1SI = m2SI;
    m2SI = tmp;
    tmp = S1x;
    S1x = S2x;
    S2x = tmp;
    tmp = S1y;
    S1y = S2y;
    S2y = tmp;
    tmp = S1z;
    S1z = S2z;
    S2z = tmp;
    phiRef += LAL_PI;
  }

  /* Get mass ratio and dimensionful-constant */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */
  double q = mass1/mass2; /* mass ratio (>= 1)*/

  // Get model parameter spins.
  double chi1mag = sqrt(S1x*S1x + S1y*S1y + S1z*S1z);
  double chi1theta = 0;
  double chi1phi = 0;
  if (chi1mag > 0.) chi1theta = acos(S1z/chi1mag);
  if (fabs(S1y) + fabs(S1x) > 0.) chi1phi = atan2(S1y*cos(phiRef) - S1x*sin(phiRef), S1x*cos(phiRef) + S1y*sin(phiRef));
  if (chi1phi < 0.) chi1phi += 2*LAL_PI;
  double chi2z = S2z;
  if (fabs(S2y) + fabs(S2x) > 0.) XLAL_ERROR(XLAL_FAILURE, "NRsurrogate does not support x or y spin components on the smaller BH\n");

  // Load data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&NRSur4d2s_is_initialized, NRSurrogate_Init_LALDATA);
#else
  NRSurrogate_Init_LALDATA();
#endif

  if(!NRSurrogate_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up NRSurrogate data - check your $LAL_DATA_PATH\n");

  // Compute dimensionless frequency sequence

  if (fHigh == 0.) {
    // lal convention: secret flag to use max model frequency
    fHigh = gsl_vector_get((&__lalsim_NRSurrogate_data)->fEI, (&__lalsim_NRSurrogate_data)->n_freqs - 1)/Mtot_sec;
    }

  UINT4 nf = (UINT4) ceil(fHigh / deltaF);
  UINT4 i0 = (UINT4) ceil(fLow / deltaF);
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(nf);
  double deltaF_dimless = deltaF*Mtot_sec;
  UINT4 i;
  for (i=0; i<i0; i++) {
    // Ensure the waveform is zero below fLow
    freqs->data[i] = 0.0;
  }
  for (i=i0; i<nf; i++) {
    freqs->data[i] = i*deltaF_dimless;
  }
  int retcode = NRSurrogateCore(hptilde,hctilde,
            phiRef, distance, inclination, Mtot_sec, q, chi1mag, chi1theta, chi1phi, chi2z, freqs);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}


static void NRSurrogate_Init_LALDATA(void)
{
  if (NRSurrogate_IsSetup()) return;

  char *path = XLALFileResolvePathLong(NRSUR4D2S_DATAFILE, PKG_DATA_DIR);

  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", NRSUR4D2S_DATAFILE);
  char *dir = dirname(path);
  int ret = NRSurrogate_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find NRSurrogateData data files in $LAL_DATA_PATH\n");
}
