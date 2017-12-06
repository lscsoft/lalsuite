/*
 *  Copyright (C) 2014, 2015, 2016 Michael Puerrer, John Veitch
 *  Reduced Order Model for SEOBNR
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
static const char ROMDataHDF5[] = "SEOBNRv4ROM_v2.0.hdf5";
static const INT4 ROMDataHDF5_VERSION_MAJOR = 2;
static const INT4 ROMDataHDF5_VERSION_MINOR = 0;
static const INT4 ROMDataHDF5_VERSION_MICRO = 0;
#endif

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif



#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv4ROM_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

typedef struct tagSEOBNRROMdataDS_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} SEOBNRROMdataDS_coeff;

struct tagSEOBNRROMdataDS_submodel
{
  gsl_vector* cvec_amp;      // Flattened amplitude projection coefficients
  gsl_vector* cvec_phi;      // Flattened phase projection coefficients
  gsl_matrix *Bamp;          // Reduced SVD basis for amplitude
  gsl_matrix *Bphi;          // Reduced SVD basis for phase
  int nk_amp;                // Number frequency points for amplitude
  int nk_phi;                // Number of frequency points for phase
  gsl_vector *gA;            // Sparse frequency points for amplitude
  gsl_vector *gPhi;          // Sparse frequency points for phase
  gsl_vector *etavec;        // B-spline knots in eta
  gsl_vector *chi1vec;       // B-spline knots in chi1
  gsl_vector *chi2vec;       // B-spline knots in chi2
  int ncx, ncy, ncz;         // Number of points in eta, chi1, chi2
  double eta_bounds[2];      // [eta_min, eta_max]
  double chi1_bounds[2];     // [chi1_min, chi1_max]
  double chi2_bounds[2];     // [chi2_min, chi2_max]
};
typedef struct tagSEOBNRROMdataDS_submodel SEOBNRROMdataDS_submodel;

struct tagSEOBNRROMdataDS
{
  UINT4 setup;
  SEOBNRROMdataDS_submodel* sub1;
  SEOBNRROMdataDS_submodel* sub2;
  SEOBNRROMdataDS_submodel* sub3;
};
typedef struct tagSEOBNRROMdataDS SEOBNRROMdataDS;

static SEOBNRROMdataDS __lalsim_SEOBNRv4ROMDS_data;

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  gsl_bspline_workspace *bwz;
} SplineData;

/**************** Internal functions **********************/

UNUSED static void SEOBNRv4ROM_Init_LALDATA(void);
UNUSED static int SEOBNRv4ROM_Init(const char dir[]);
UNUSED static bool SEOBNRv4ROM_IsSetup(void);

UNUSED static int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[]);
UNUSED static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata);

static int TP_Spline_interpolation_3d(
  REAL8 eta,                // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
//  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  int nk_amp,               // number of SVD-modes == number of basis functions for amplitude
  int nk_phi,               // number of SVD-modes == number of basis functions for phase
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *etavec,     // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi         // Output: interpolated projection coefficients for phase
//  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
);

UNUSED static int SEOBNRROMdataDS_Init_submodel(
  UNUSED SEOBNRROMdataDS_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[]
);

UNUSED static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
UNUSED static int SEOBNRv4ROMCore(
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
  const REAL8Sequence *freqs, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  int nk_max // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
);

UNUSED static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi);
UNUSED static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
UNUSED static void SplineData_Destroy(SplineData *splinedata);
UNUSED static void SplineData_Init(
  SplineData **splinedata,
  int ncx,                // Number of points in eta  + 2
  int ncy,                // Number of points in chi1 + 2
  int ncz,                // Number of points in chi2 + 2
  const double *etavec,   // B-spline knots in eta
  const double *chi1vec,  // B-spline knots in chi1
  const double *chi2vec   // B-spline knots in chi2
);

UNUSED static int SEOBNRv4ROMTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi1,                                   // Aligned spin of companion 1
  REAL8 chi2,                                   // Aligned spin of companion 2
  REAL8 *Mf_ROM_min,                            // Lowest geometric frequency for ROM
  REAL8 *Mf_ROM_max                             // Highest geometric frequency for ROM
);

UNUSED static REAL8 Interpolate_Coefficent_Matrix(
  gsl_vector *v,
  REAL8 eta,
  REAL8 chi,
  int ncx,
  int ncy,
  gsl_bspline_workspace *bwx,
  gsl_bspline_workspace *bwy
);

UNUSED static void GlueAmplitude(
  // INPUTS
  SEOBNRROMdataDS_submodel *submodel_lo,
  SEOBNRROMdataDS_submodel *submodel_hi,
  gsl_vector* amp_f_lo,
  gsl_vector* amp_f_hi,
  double amp_pre_lo,
  double amp_pre_hi,
  const double Mfm,
  // OUTPUTS
  gsl_interp_accel **acc_amp,
  gsl_spline **spline_amp
);

UNUSED static void GluePhasing(
  // INPUTS
  SEOBNRROMdataDS_submodel *submodel_lo,
  SEOBNRROMdataDS_submodel *submodel_hi,
  gsl_vector* phi_f_lo,
  gsl_vector* phi_f_hi,
  const double Mfm,
  // OUTPUTS
  gsl_interp_accel **acc_phi_out,
  gsl_spline **spline_phi_out
);


/********************* Definitions begin here ********************/

/** Setup SEOBNRv4ROM model using data files installed in dir
 */
static int SEOBNRv4ROM_Init(const char dir[]) {
  if(__lalsim_SEOBNRv4ROMDS_data.setup) {
    XLALPrintError("Error: SEOBNRv4ROM data was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  SEOBNRROMdataDS_Init(&__lalsim_SEOBNRv4ROMDS_data, dir);

  if(__lalsim_SEOBNRv4ROMDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv4ROM model has been initialised */
static bool SEOBNRv4ROM_IsSetup(void) {
  if(__lalsim_SEOBNRv4ROMDS_data.setup)
    return true;
  else
    return false;
}

// Setup B-spline basis functions for given points
static void SplineData_Init(
  SplineData **splinedata,
  int ncx,                // Number of points in eta  + 2
  int ncy,                // Number of points in chi1 + 2
  int ncz,                // Number of points in chi2 + 2
  const double *etavec,   // B-spline knots in eta
  const double *chi1vec,  // B-spline knots in chi1
  const double *chi2vec   // B-spline knots in chi2
)
{
  if(!splinedata) exit(1);
  if(*splinedata) SplineData_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData));

  // Set up B-spline basis for desired knots
  const size_t nbreak_x = ncx-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_y = ncy-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_z = ncz-2;  // must have nbreak = n-2 for cubic splines

  // Allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bwx = gsl_bspline_alloc(4, nbreak_x);
  gsl_bspline_workspace *bwy = gsl_bspline_alloc(4, nbreak_y);
  gsl_bspline_workspace *bwz = gsl_bspline_alloc(4, nbreak_z);

  // Set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_x = gsl_vector_alloc(nbreak_x);
  gsl_vector *breakpts_y = gsl_vector_alloc(nbreak_y);
  gsl_vector *breakpts_z = gsl_vector_alloc(nbreak_z);
  for (UINT4 i=0; i<nbreak_x; i++)
    gsl_vector_set(breakpts_x, i, etavec[i]);
  for (UINT4 j=0; j<nbreak_y; j++)
    gsl_vector_set(breakpts_y, j, chi1vec[j]);
  for (UINT4 k=0; k<nbreak_z; k++)
    gsl_vector_set(breakpts_z, k, chi2vec[k]);

  gsl_bspline_knots(breakpts_x, bwx);
  gsl_bspline_knots(breakpts_y, bwy);
  gsl_bspline_knots(breakpts_z, bwz);

  gsl_vector_free(breakpts_x);
  gsl_vector_free(breakpts_y);
  gsl_vector_free(breakpts_z);

  (*splinedata)->bwx=bwx;
  (*splinedata)->bwy=bwy;
  (*splinedata)->bwz=bwz;
}

static void SplineData_Destroy(SplineData *splinedata)
{
  if(!splinedata) return;
  if(splinedata->bwx) gsl_bspline_free(splinedata->bwx);
  if(splinedata->bwy) gsl_bspline_free(splinedata->bwy);
  if(splinedata->bwz) gsl_bspline_free(splinedata->bwz);
  XLALFree(splinedata);
}

// Interpolate projection coefficients for amplitude and phase over the parameter space (q, chi).
// The multi-dimensional interpolation is carried out via a tensor product decomposition.
static int TP_Spline_interpolation_3d(
  REAL8 eta,                // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  int nk_amp,               // number of SVD-modes == number of basis functions for amplitude
  int nk_phi,               // number of SVD-modes == number of basis functions for phase
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *etavec,     // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi        // Output: interpolated projection coefficients for phase
  ) {
  if (nk_max != -1) {
    if (nk_max > nk_amp || nk_max > nk_phi)
      XLAL_ERROR(XLAL_EDOM, "Truncation parameter nk_max %d must be smaller or equal to nk_amp %d and nk_phi %d", nk_max, nk_amp, nk_phi);
    else { // truncate SVD modes
      nk_amp = nk_max;
      nk_phi = nk_max;
    }
  }

  SplineData *splinedata=NULL;
  SplineData_Init(&splinedata, ncx, ncy, ncz, etavec, chi1vec, chi2vec);

  gsl_bspline_workspace *bwx=splinedata->bwx;
  gsl_bspline_workspace *bwy=splinedata->bwy;
  gsl_bspline_workspace *bwz=splinedata->bwz;

  int N = ncx*ncy*ncz;  // Size of the data matrix for one SVD-mode
  // Evaluate the TP spline for all SVD modes - amplitude
  for (int k=0; k<nk_amp; k++) { // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Tensor(&v, eta, chi1, chi2, ncy, ncz, bwx, bwy, bwz);
    gsl_vector_set(c_amp, k, csum);
  }

  // Evaluate the TP spline for all SVD modes - phase
  for (int k=0; k<nk_phi; k++) {  // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Tensor(&v, eta, chi1, chi2, ncy, ncz, bwx, bwy, bwz);
    gsl_vector_set(c_phi, k, csum);
  }

  SplineData_Destroy(splinedata);

  return(0);
}

/* Set up a new ROM submodel, using data contained in dir */
UNUSED static int SEOBNRROMdataDS_Init_submodel(
  SEOBNRROMdataDS_submodel **submodel,
  UNUSED const char dir[],
  UNUSED const char grp_name[]
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(SEOBNRROMdataDS_submodel));
  else
    SEOBNRROMdataDS_Cleanup_submodel(*submodel);

#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  LALH5File *sub = XLALH5GroupOpen(file, grp_name);

  // Read ROM coefficients
  ReadHDF5RealVectorDataset(sub, "Amp_ciall", & (*submodel)->cvec_amp);
  ReadHDF5RealVectorDataset(sub, "Phase_ciall", & (*submodel)->cvec_phi);

  // Read ROM basis functions
  ReadHDF5RealMatrixDataset(sub, "Bamp", & (*submodel)->Bamp);
  ReadHDF5RealMatrixDataset(sub, "Bphase", & (*submodel)->Bphi);

  // Read sparse frequency points
  ReadHDF5RealVectorDataset(sub, "Mf_grid_Amp", & (*submodel)->gA);
  ReadHDF5RealVectorDataset(sub, "Mf_grid_Phi", & (*submodel)->gPhi);

  // Read parameter space nodes
  ReadHDF5RealVectorDataset(sub, "etavec", & (*submodel)->etavec);
  ReadHDF5RealVectorDataset(sub, "chi1vec", & (*submodel)->chi1vec);
  ReadHDF5RealVectorDataset(sub, "chi2vec", & (*submodel)->chi2vec);

  // Initialize other members
  (*submodel)->nk_amp = (*submodel)->gA->size;
  (*submodel)->nk_phi = (*submodel)->gPhi->size;
  (*submodel)->ncx = (*submodel)->etavec->size + 2;
  (*submodel)->ncy = (*submodel)->chi1vec->size + 2;
  (*submodel)->ncz = (*submodel)->chi2vec->size + 2;

  // Domain of definition of submodel
  (*submodel)->eta_bounds[0] = gsl_vector_get((*submodel)->etavec, 0);
  (*submodel)->eta_bounds[1] = gsl_vector_get((*submodel)->etavec, (*submodel)->etavec->size - 1);
  (*submodel)->chi1_bounds[0] = gsl_vector_get((*submodel)->chi1vec, 0);
  (*submodel)->chi1_bounds[1] = gsl_vector_get((*submodel)->chi1vec, (*submodel)->chi1vec->size - 1);
  (*submodel)->chi2_bounds[0] = gsl_vector_get((*submodel)->chi2vec, 0);
  (*submodel)->chi2_bounds[1] = gsl_vector_get((*submodel)->chi2vec, (*submodel)->chi2vec->size - 1);

  XLALFree(path);
  XLALH5FileClose(file);
  ret = XLAL_SUCCESS;
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  return ret;
}

/* Deallocate contents of the given SEOBNRROMdataDS_submodel structure */
static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel) {
  if(submodel->cvec_amp) gsl_vector_free(submodel->cvec_amp);
  if(submodel->cvec_phi) gsl_vector_free(submodel->cvec_phi);
  if(submodel->Bamp) gsl_matrix_free(submodel->Bamp);
  if(submodel->Bphi) gsl_matrix_free(submodel->Bphi);
  if(submodel->gA)   gsl_vector_free(submodel->gA);
  if(submodel->gPhi) gsl_vector_free(submodel->gPhi);
  if(submodel->etavec)  gsl_vector_free(submodel->etavec);
  if(submodel->chi1vec) gsl_vector_free(submodel->chi1vec);
  if(submodel->chi2vec) gsl_vector_free(submodel->chi2vec);
}

/* Set up a new ROM model, using data contained in dir */
int SEOBNRROMdataDS_Init(
  UNUSED SEOBNRROMdataDS *romdata,
  UNUSED const char dir[])
{
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the SEOBNRv4ROM model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  gsl_set_error_handler(&err_handler);

#ifdef LAL_HDF5_ENABLED
  // First, check we got the correct version number
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);
  LALH5File *file = XLALH5FileOpen(path, "r");

  XLALPrintInfo("ROM metadata\n============\n");
  PrintInfoStringAttribute(file, "Email");
  PrintInfoStringAttribute(file, "Description");
  ret = ROM_check_version_number(file, ROMDataHDF5_VERSION_MAJOR,
                                 ROMDataHDF5_VERSION_MINOR,
                                 ROMDataHDF5_VERSION_MICRO);

  XLALFree(path);
  XLALH5FileClose(file);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->sub1, dir, "sub1");
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 1 loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->sub2, dir, "sub2");
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 2 loaded sucessfully.\n", __func__);

  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->sub3, dir, "sub3");
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 3 loaded sucessfully.\n", __func__);

  if(XLAL_SUCCESS==ret)
    romdata->setup=1;
  else
    SEOBNRROMdataDS_Cleanup(romdata);
#else
  XLAL_ERROR(XLAL_EFAILED, "HDF5 support not enabled");
#endif

  return (ret);
}

/* Deallocate contents of the given SEOBNRROMdataDS structure */
static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata) {
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->sub1);
  XLALFree((romdata)->sub1);
  (romdata)->sub1 = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->sub2);
  XLALFree((romdata)->sub2);
  (romdata)->sub2 = NULL;
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->sub3);
  XLALFree((romdata)->sub3);
  (romdata)->sub3 = NULL;
  romdata->setup=0;
}

/* Structure for internal use */
static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi) {
  if(!romdatacoeff) exit(1);
  /* Create storage for structures */
  if(!*romdatacoeff)
    *romdatacoeff=XLALCalloc(1,sizeof(SEOBNRROMdataDS_coeff));
  else
    SEOBNRROMdataDS_coeff_Cleanup(*romdatacoeff);

  (*romdatacoeff)->c_amp = gsl_vector_alloc(nk_amp);
  (*romdatacoeff)->c_phi = gsl_vector_alloc(nk_phi);
}

/* Deallocate contents of the given SEOBNRROMdataDS_coeff structure */
static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff) {
  if(romdatacoeff->c_amp) gsl_vector_free(romdatacoeff->c_amp);
  if(romdatacoeff->c_phi) gsl_vector_free(romdatacoeff->c_phi);
  XLALFree(romdatacoeff);
}

/* Return the closest higher power of 2  */
// Note: NextPow(2^k) = 2^k for integer values k.
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

static void GlueAmplitude(
  // INPUTS
  SEOBNRROMdataDS_submodel *submodel_lo,
  SEOBNRROMdataDS_submodel *submodel_hi,
  gsl_vector* amp_f_lo,
  gsl_vector* amp_f_hi,
  // amp_pre_* can be set to 1 if not using amplitude prefactors
  double amp_pre_lo,
  double amp_pre_hi,
  const double Mfm,
  // OUTPUTS
  gsl_interp_accel **acc_amp,
  gsl_spline **spline_amp
) {
  // First need to find overlaping frequency interval
  int jA_lo;
  // Find index so that Mf < Mfm
  for (jA_lo=0; jA_lo < submodel_lo->nk_amp; jA_lo++) {
    if (gsl_vector_get(submodel_lo->gA, jA_lo) > Mfm) {
      jA_lo--;
      break;
    }
  }

  int jA_hi;
  // Find index so that Mf > Mfm
  for (jA_hi=0; jA_hi < submodel_hi->nk_amp; jA_hi++)
    if (gsl_vector_get(submodel_hi->gA, jA_hi) > Mfm)
      break;

  int nA = 1 + jA_lo + (submodel_hi->nk_amp - jA_hi); // length of the union of frequency points of the low and high frequency models glued at MfM

  gsl_vector *gAU = gsl_vector_alloc(nA); // glued frequency grid
  gsl_vector *amp_f = gsl_vector_alloc(nA); // amplitude on glued frequency grid
  // Note: We don't interpolate the amplitude, but this may already be smooth enough for practical purposes.
  // To improve this we would evaluate both amplitue splines times the prefactor at the matching frequency and correct with the ratio, so we are C^0.
  for (int i=0; i<=jA_lo; i++) {
    gsl_vector_set(gAU, i, gsl_vector_get(submodel_lo->gA, i));
    double A = amp_pre_lo * gsl_vector_get(amp_f_lo, i);
    gsl_vector_set(amp_f, i, A);
  }

  for (int i=jA_lo+1; i<nA; i++) {
    int k = jA_hi - (jA_lo+1) + i;
    gsl_vector_set(gAU, i, gsl_vector_get(submodel_hi->gA, k));
    double A = amp_pre_hi * gsl_vector_get(amp_f_hi, k);
    gsl_vector_set(amp_f, i, A);
  }

  // Setup 1d splines in frequency from glued amplitude grids & data
  *acc_amp = gsl_interp_accel_alloc();
  *spline_amp = gsl_spline_alloc(gsl_interp_cspline, nA);
  //gsl_spline_init(spline_amp, gAU->data, gsl_vector_const_ptr(amp_f,0), nA);
  gsl_spline_init(*spline_amp, gsl_vector_const_ptr(gAU,0),
                  gsl_vector_const_ptr(amp_f,0), nA);

  gsl_vector_free(gAU);
  gsl_vector_free(amp_f);
  gsl_vector_free(amp_f_lo);
  gsl_vector_free(amp_f_hi);
}

// Glue phasing in frequency to C^1 smoothness
static void GluePhasing(
  // INPUTS
  SEOBNRROMdataDS_submodel *submodel_lo,
  SEOBNRROMdataDS_submodel *submodel_hi,
  gsl_vector* phi_f_lo,
  gsl_vector* phi_f_hi,
  const double Mfm,
  // OUTPUTS
  gsl_interp_accel **acc_phi,
  gsl_spline **spline_phi
) {
  // First need to find overlaping frequency interval
  int jP_lo;
  // Find index so that Mf < Mfm
  for (jP_lo=0; jP_lo < submodel_lo->nk_phi; jP_lo++) {
    if (gsl_vector_get(submodel_lo->gPhi, jP_lo) > Mfm) {
      jP_lo--;
      break;
    }
  }

  int jP_hi;
  // Find index so that Mf > Mfm
  for (jP_hi=0; jP_hi < submodel_hi->nk_phi; jP_hi++)
    if (gsl_vector_get(submodel_hi->gPhi, jP_hi) > Mfm)
      break;

  int nP = 1 + jP_lo + (submodel_hi->nk_phi - jP_hi); // length of the union of frequency points of the low and high frequency models glued at MfM
  gsl_vector *gPU = gsl_vector_alloc(nP); // glued frequency grid
  gsl_vector *phi_f = gsl_vector_alloc(nP); // phase on glued frequency grid
  // We need to do a bit more work to glue the phase with C^1 smoothness
  for (int i=0; i<=jP_lo; i++) {
    gsl_vector_set(gPU, i, gsl_vector_get(submodel_lo->gPhi, i));
    double P = gsl_vector_get(phi_f_lo, i);
    gsl_vector_set(phi_f, i, P);
  }

  for (int i=jP_lo+1; i<nP; i++) {
    int k = jP_hi - (jP_lo+1) + i;
    gsl_vector_set(gPU, i, gsl_vector_get(submodel_hi->gPhi, k));
    double P = gsl_vector_get(phi_f_hi, k);
    gsl_vector_set(phi_f, i, P);
  }

  // Set up phase data across the gluing frequency Mfm
  // We need to set up a spline for the low frequency model and evaluate at the designated points for the high frequency model so that we work with the *same* frequeny interval!

  // We could optimize this further by not constructing the whole spline for
  // submodel_lo, but this may be insignificant since the number of points is small anyway.
  gsl_interp_accel *acc_phi_lo = gsl_interp_accel_alloc();
  gsl_spline *spline_phi_lo = gsl_spline_alloc(gsl_interp_cspline, submodel_lo->nk_phi);
  gsl_spline_init(spline_phi_lo, gsl_vector_const_ptr(submodel_lo->gPhi,0),
                  gsl_vector_const_ptr(phi_f_lo,0), submodel_lo->nk_phi);

  const int nn = 15;
  gsl_vector_const_view gP_hi_data = gsl_vector_const_subvector(submodel_hi->gPhi, jP_hi - nn, 2*nn+1);
  gsl_vector_const_view P_hi_data = gsl_vector_const_subvector(phi_f_hi, jP_hi - nn, 2*nn+1);
  gsl_vector *P_lo_data = gsl_vector_alloc(2*nn+1);
  for (int i=0; i<2*nn+1; i++) {
    double P = gsl_spline_eval(spline_phi_lo, gsl_vector_get(&gP_hi_data.vector, i), acc_phi_lo);
    gsl_vector_set(P_lo_data, i, P);
  }

  // Fit phase data to cubic polynomial in frequency
  gsl_vector *cP_lo = Fit_cubic(&gP_hi_data.vector, P_lo_data);
  gsl_vector *cP_hi = Fit_cubic(&gP_hi_data.vector, &P_hi_data.vector);

  double P_lo_derivs[2];
  double P_hi_derivs[2];
  gsl_poly_eval_derivs(cP_lo->data, 4, Mfm, P_lo_derivs, 2);
  gsl_poly_eval_derivs(cP_hi->data, 4, Mfm, P_hi_derivs, 2);

  double delta_omega = P_hi_derivs[1] - P_lo_derivs[1];
  double delta_phi   = P_hi_derivs[0] - P_lo_derivs[0] - delta_omega * Mfm;

  for (int i=jP_lo+1; i<nP; i++) {
    int k = jP_hi - (jP_lo+1) + i;
    double f = gsl_vector_get(submodel_hi->gPhi, k);
    gsl_vector_set(gPU, i, f);
    double P = gsl_vector_get(phi_f_hi, k) - delta_omega * f - delta_phi; // Now correct phase of high frequency submodel
    gsl_vector_set(phi_f, i, P);
  }

  // free some vectors
  gsl_vector_free(P_lo_data);
  gsl_vector_free(cP_lo);
  gsl_vector_free(cP_hi);
  gsl_vector_free(phi_f_lo);
  gsl_vector_free(phi_f_hi);

  // Setup 1d splines in frequency from glued phase grids & data
  *acc_phi = gsl_interp_accel_alloc();
  *spline_phi = gsl_spline_alloc(gsl_interp_cspline, nP);
  //gsl_spline_init(spline_phi, gPU->data, gsl_vector_const_ptr(phi_f,0), nP);
  gsl_spline_init(*spline_phi, gsl_vector_const_ptr(gPU,0),
                  gsl_vector_const_ptr(phi_f,0), nP);

  /**** Finished gluing ****/

  gsl_vector_free(phi_f);
  gsl_vector_free(gPU);
  gsl_spline_free(spline_phi_lo);
  gsl_interp_accel_free(acc_phi_lo);
}



/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi1, chi2).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv4ROMCore(
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
  const REAL8Sequence *freqs_in, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  int nk_max // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
  )
{

  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);

  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4ROMDS_data;
  if (!SEOBNRv4ROM_IsSetup()) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }

  if(*hptilde || *hctilde) {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",
                   (*hptilde), (*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (eta > 0.25)     nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01)     nudge(&eta, 0.01, 1e-6);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 1.0 || chi2 > 1.0) {
    XLALPrintError("XLAL Error - %s: chi1 or chi2 smaller than -1.0 or larger than 1.0!\n"
                   "SEOBNRv4ROM is only available for spins in the range -1 <= a/M <= 1.0.\n",
                   __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta<0.01 || eta > 0.25) {
    XLALPrintError("XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\n"
                   "SEOBNRv4ROM is only available for eta in the range 0.01 <= eta <= 0.25.\n",
                   __func__, eta);
    XLAL_ERROR( XLAL_EDOM );
  }

  /* We always need to glue two submodels together for this ROM */
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM
  submodel_lo = romdata->sub1;

  /* Select high frequency ROM submodel */
  if (chi1 < romdata->sub3->chi1_bounds[0] || eta > romdata->sub3->eta_bounds[1]) // only check the two conditions that apply for this ROM; could be more general, but slower
    submodel_hi = romdata->sub2;
  else
    submodel_hi = romdata->sub3;


  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];

  if(fRef==0.0)
    fRef=fLow;

  /* Convert to geometric units for frequency */
  // lowest allowed geometric frequency for ROM
  double Mf_ROM_min = fmax(gsl_vector_get(submodel_lo->gA, 0),
                           gsl_vector_get(submodel_lo->gPhi,0));
  // highest allowed geometric frequency for ROM
  double Mf_ROM_max = fmin(gsl_vector_get(submodel_hi->gA, submodel_hi->nk_amp-1),
                           gsl_vector_get(submodel_hi->gPhi, submodel_hi->nk_phi-1));
  double fLow_geom = fLow * Mtot_sec;
  double fHigh_geom = fHigh * Mtot_sec;
  double fRef_geom = fRef * Mtot_sec;
  double deltaF_geom = deltaF * Mtot_sec;

  // Enforce allowed geometric frequency range
  if (fLow_geom < Mf_ROM_min)
    XLAL_ERROR(XLAL_EDOM, "Starting frequency Mflow=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
  if (fHigh_geom == 0 || fHigh_geom > Mf_ROM_max)
    fHigh_geom = Mf_ROM_max;
  else if (fHigh_geom < Mf_ROM_min)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than ROM starting frequency %g!\n", fHigh_geom, Mf_ROM_min);
  if (fHigh_geom <= fLow_geom)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than (or equal to) starting frequency %g!\n", fHigh_geom, fLow_geom);
  if (fRef_geom > Mf_ROM_max) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is greater than maximal frequency in ROM Mf=%g. Starting at maximal frequency in ROM.\n", fRef_geom, Mf_ROM_max);
    fRef_geom = Mf_ROM_max; // If fref > fhigh we reset fref to default value of cutoff frequency.
  }
  if (fRef_geom < Mf_ROM_min) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for waveform coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff_lo=NULL;
  SEOBNRROMdataDS_coeff *romdata_coeff_hi=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_lo, submodel_lo->nk_amp, submodel_lo->nk_phi);
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_hi, submodel_hi->nk_amp, submodel_hi->nk_phi);
  REAL8 amp_pre_lo = 1.0; // unused here
  REAL8 amp_pre_hi = 1.0;

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_lo->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_lo->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_lo->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_lo->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_lo->ncx,             // Number of points in eta  + 2
    submodel_lo->ncy,             // Number of points in chi1 + 2
    submodel_lo->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel_lo->etavec, 0),          // B-spline knots in eta
    gsl_vector_const_ptr(submodel_lo->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel_lo->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff_lo->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_lo->c_phi       // Output: interpolated projection coefficients for phase
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
    XLAL_ERROR(retcode);
  }

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_hi->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_hi->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_hi->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_hi->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_hi->ncx,             // Number of points in eta  + 2
    submodel_hi->ncy,             // Number of points in chi1 + 2
    submodel_hi->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel_hi->etavec, 0),         // B-spline knots in eta
    gsl_vector_const_ptr(submodel_hi->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel_hi->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff_hi->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_hi->c_phi       // Output: interpolated projection coefficients for phase
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);
    XLAL_ERROR(retcode);
  }


  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // amp_pts = B_A^T . c_A
  // phi_pts = B_phi^T . c_phi
  gsl_vector* amp_f_lo = gsl_vector_alloc(submodel_lo->nk_amp);
  gsl_vector* phi_f_lo = gsl_vector_alloc(submodel_lo->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_lo->Bamp, romdata_coeff_lo->c_amp, 0.0, amp_f_lo);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_lo->Bphi, romdata_coeff_lo->c_phi, 0.0, phi_f_lo);

  gsl_vector* amp_f_hi = gsl_vector_alloc(submodel_hi->nk_amp);
  gsl_vector* phi_f_hi = gsl_vector_alloc(submodel_hi->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_hi->Bamp, romdata_coeff_hi->c_amp, 0.0, amp_f_hi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_hi->Bphi, romdata_coeff_hi->c_phi, 0.0, phi_f_hi);

  const double Mfm = 0.01; // Gluing frequency: the low and high frequency ROMs overlap here; this is used both for amplitude and phase.

  // Glue amplitude
  gsl_interp_accel *acc_amp;
  gsl_spline *spline_amp;
  GlueAmplitude(submodel_lo, submodel_hi, amp_f_lo, amp_f_hi, amp_pre_lo, amp_pre_hi, Mfm,
    &acc_amp, &spline_amp
  );

  // Glue phasing in frequency to C^1 smoothness
  gsl_interp_accel *acc_phi;
  gsl_spline *spline_phi;
  GluePhasing(submodel_lo, submodel_hi, phi_f_lo, phi_f_hi, Mfm,
    &acc_phi, &spline_phi
  );

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
      gsl_spline_free(spline_amp);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_amp);
      gsl_interp_accel_free(acc_phi);
      SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
      SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);
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

  REAL8 s = 0.5; // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv4
  double Mtot = Mtot_sec / LAL_MTSUN_SI;
  double amp0 = Mtot * Mtot_sec * LAL_MRSUN_SI / (distance); // Correct overall amplitude to undo mass-dependent scaling used in ROM

  // Evaluate reference phase for setting phiRef correctly
  double phase_change = gsl_spline_eval(spline_phi, fRef_geom, acc_phi) - 2*phiRef;

  // Assemble waveform from aplitude and phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i];
    if (f > Mf_ROM_max) continue; // We're beyond the highest allowed frequency; since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
    int j = i + offset; // shift index for frequency series if needed
    double A = gsl_spline_eval(spline_amp, f, acc_amp);
    double phase = gsl_spline_eval(spline_phi, f, acc_phi) - phase_change;
    COMPLEX16 htilde = s*amp0*A * (cos(phase) + I*sin(phase));//cexp(I*phase);
    pdata[j] =      pcoef * htilde;
    cdata[j] = -I * ccoef * htilde;
  }

  /* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */

  // Get SEOBNRv4 ringdown frequency for 22 mode
  double Mf_final = SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(Mtot_sec, eta, chi1,
                                                        chi2, SEOBNRv4);

  UINT4 L = freqs->length;
  // prevent gsl interpolation errors
  if (Mf_final > freqs->data[L-1])
    Mf_final = freqs->data[L-1];
  if (Mf_final < freqs->data[0]) {
    XLALDestroyREAL8Sequence(freqs);
    gsl_spline_free(spline_amp);
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(acc_amp);
    gsl_interp_accel_free(acc_phi);
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);
    XLAL_ERROR(XLAL_EDOM, "f_ringdown < f_min");
  }

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // We compute the dimensionless time correction t/M since we use geometric units.
  REAL8 t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI);

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

  gsl_spline_free(spline_amp);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_amp);
  gsl_interp_accel_free(acc_phi);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRv4ROM_c
 *
 * \author Michael Puerrer
 *
 * \brief C code for SEOBNRv4 reduced order model
 * See CQG 31 195010, 2014, arXiv:1402.4146 for the basic approach.
 * Further details in PRD 93, 064041, 2016, arXiv:1512.02248.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv4 model.
 *
 * The binary data HDF5 file (SEOBNRv4ROM_DS_HI_v1.0.hdf5)
 * will be available at on LIGO clusters in /home/cbc/.
 * Make sure the files are in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv4 model.
 *
 * @note Parameter ranges:
 *   * 0.01 <= eta <= 0.25
 *   * -1 <= chi_i <= 0.99
 *   * Mtot >= 3 Msun
 *
 *  Aligned component spins chi1, chi2.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * This ROM consists of three submodels and glues together one low-mass and 2 high-mass models
 * These submodels and theor boundaries are not explicit in the source, just in the HDF5 data file.
 *
 * @{
 */


/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv4_ROM model.
 *
 * XLALSimIMRSEOBNRv4ROM() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv4ROMFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv4ROMFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv4ROMFrequencySequence(
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
  INT4 nk_max)                                  /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
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
  (void) pthread_once(&SEOBNRv4ROM_is_initialized, SEOBNRv4ROM_Init_LALDATA);
#else
  SEOBNRv4ROM_Init_LALDATA();
#endif

  if(!SEOBNRv4ROM_IsSetup()) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv4ROMCore(hptilde, hctilde, phiRef, fRef, distance,
                                inclination, Mtot_sec, eta, chi1, chi2, freqs,
                                0, nk_max);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv4_ROM model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv4ROM(
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
  INT4 nk_max)                                  /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
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
  (void) pthread_once(&SEOBNRv4ROM_is_initialized, SEOBNRv4ROM_Init_LALDATA);
#else
  SEOBNRv4ROM_Init_LALDATA();
#endif

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv4ROMCore(hptilde, hctilde, phiRef, fRef, distance,
                                inclination, Mtot_sec, eta, chi1, chi2, freqs,
                                deltaF, nk_max);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */

// Auxiliary function to perform setup of phase spline for t(f) and f(t) functions
static int SEOBNRv4ROMTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi1,                                   // Aligned spin of companion 1
  REAL8 chi2,                                   // Aligned spin of companion 2
  REAL8 *Mf_ROM_min,                            // Lowest geometric frequency for ROM
  REAL8 *Mf_ROM_max                             // Highest geometric frequency for ROM
)
{
  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1 + mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  *Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  // 'Nudge' parameter values to allowed boundary values if close by
  if (eta > 0.25)     nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01)     nudge(&eta, 0.01, 1e-6);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 1.0 || chi2 > 1.0) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1.0 or larger than 1.0!\nSEOBNRv4ROM is only available for spins in the range -1 <= a/M <= 1.0.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta < 0.01 || eta > 0.25) {
    XLALPrintError( "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv4ROM is only available for spins in the range 0.01 <= eta <= 0.25.\n", __func__,eta);
    XLAL_ERROR( XLAL_EDOM );
  }

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv4ROM_is_initialized, SEOBNRv4ROM_Init_LALDATA);
#else
  SEOBNRv4ROM_Init_LALDATA();
#endif

  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv4ROMDS_data;
  if (!SEOBNRv4ROM_IsSetup()) {
    XLAL_ERROR(XLAL_EFAILED,
               "Error setting up SEOBNRv4ROM data - check your $LAL_DATA_PATH\n");
  }

  /* We always need to glue two submodels together for this ROM */
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM
  submodel_lo = romdata->sub1;

  /* Select high frequency ROM submodel */
  if (chi1 < romdata->sub3->chi1_bounds[0] || eta > romdata->sub3->eta_bounds[1]) // only check the two conditions that apply for this ROM; could be more general, but slower
    submodel_hi = romdata->sub2;
  else
    submodel_hi = romdata->sub3;

  /* Internal storage for waveform coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff_lo=NULL;
  SEOBNRROMdataDS_coeff *romdata_coeff_hi=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_lo, submodel_lo->nk_amp,
                             submodel_lo->nk_phi);
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_hi, submodel_hi->nk_amp,
                             submodel_hi->nk_phi);

  // lowest allowed geometric frequency for ROM
  *Mf_ROM_min = fmax(gsl_vector_get(submodel_lo->gA, 0),
                     gsl_vector_get(submodel_lo->gPhi,0));
  // highest allowed geometric frequency for ROM
  *Mf_ROM_max = fmin(gsl_vector_get(submodel_hi->gA, submodel_hi->nk_amp-1),
                     gsl_vector_get(submodel_hi->gPhi, submodel_hi->nk_phi-1));


  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  int nk_max = -1; // adjust truncation parameter if speed is an issue
  int retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_lo->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_lo->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_lo->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_lo->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_lo->ncx,             // Number of points in eta  + 2
    submodel_lo->ncy,             // Number of points in chi1 + 2
    submodel_lo->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel_lo->etavec, 0),         // B-spline knots in eta
    gsl_vector_const_ptr(submodel_lo->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel_lo->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff_lo->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_lo->c_phi       // Output: interpolated projection coefficients for phase
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
    XLAL_ERROR(retcode);
  }

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_hi->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_hi->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_hi->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_hi->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_hi->ncx,             // Number of points in eta  + 2
    submodel_hi->ncy,             // Number of points in chi1 + 2
    submodel_hi->ncz,             // Number of points in chi2 + 2
    gsl_vector_const_ptr(submodel_hi->etavec, 0),         // B-spline knots in eta
    gsl_vector_const_ptr(submodel_hi->chi1vec, 0),        // B-spline knots in chi1
    gsl_vector_const_ptr(submodel_hi->chi2vec, 0),        // B-spline knots in chi2
    romdata_coeff_hi->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_hi->c_phi       // Output: interpolated projection coefficients for phase
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);
    XLAL_ERROR(retcode);
  }

  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // phi_pts = B_phi^T . c_phi
  gsl_vector* phi_f_lo = gsl_vector_alloc(submodel_lo->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_lo->Bphi, romdata_coeff_lo->c_phi,
                 0.0, phi_f_lo);

  gsl_vector* phi_f_hi = gsl_vector_alloc(submodel_hi->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_hi->Bphi, romdata_coeff_hi->c_phi,
                 0.0, phi_f_hi);

  const double Mfm = 0.01; // Gluing frequency: the low and high frequency ROMs overlap here; this is used both for amplitude and phase.

  // Glue phasing in frequency to C^1 smoothness
  GluePhasing(submodel_lo, submodel_hi, phi_f_lo, phi_f_hi, Mfm,
    acc_phi, spline_phi
  );

  // Get SEOBNRv4 ringdown frequency for 22 mode
  *Mf_final = SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(*Mtot_sec, eta, chi1, chi2,
                                                  SEOBNRv4);

  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);

  return(XLAL_SUCCESS);
}

/**
 * Compute the 'time' elapsed in the ROM waveform from a given starting frequency until the ringdown.
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv4.
 *
 * The allowed frequency range for the starting frequency in geometric frequency is [0.00053, 0.135].
 * The SEOBNRv4 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv4ROMFrequencyOfTime() for the inverse function.
 */
int XLALSimIMRSEOBNRv4ROMTimeOfFrequency(
  REAL8 *t,         /**< Output: time (s) elapsed from starting frequency to ringdown */
  REAL8 frequency,  /**< Starting frequency (Hz) */
  REAL8 m1SI,       /**< Mass of companion 1 (kg) */
  REAL8 m2SI,       /**< Mass of companion 2 (kg) */
  REAL8 chi1,       /**< Dimensionless aligned component spin 1 */
  REAL8 chi2        /**< Dimensionless aligned component spin 2 */
)
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
  }

  // Set up phase spline
  gsl_spline *spline_phi;
  gsl_interp_accel *acc_phi;
  double Mf_final, Mtot_sec;
  double Mf_ROM_min, Mf_ROM_max;
  int ret = SEOBNRv4ROMTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final,
                                          &Mtot_sec, m1SI, m2SI, chi1, chi2,
                                          &Mf_ROM_min, &Mf_ROM_max);
  if(ret != 0)
    XLAL_ERROR(ret);

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  //XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  double Mf = frequency * Mtot_sec;
  if (Mf < Mf_ROM_min || Mf > Mf_ROM_max) {
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(acc_phi);
    XLAL_ERROR(XLAL_EDOM, "Frequency %g is outside allowed frequency range.\n", frequency);
   }

  // Compute time relative to origin at merger
  double time_M = gsl_spline_eval_deriv(spline_phi, frequency * Mtot_sec, acc_phi) / (2*LAL_PI) - t_corr;
  *t = time_M * Mtot_sec;

  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_phi);

  return(XLAL_SUCCESS);
}

/**
 * Compute the starting frequency so that the given amount of 'time' elapses in the ROM waveform
 * from the starting frequency until the ringdown.
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv4.
 *
 * If the frequency that corresponds to the specified elapsed time is lower than the
 * geometric frequency Mf=0.00053 (ROM starting frequency) or above half of the SEOBNRv4
 * ringdown frequency an error is thrown.
 * The SEOBNRv4 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv4ROMTimeOfFrequency() for the inverse function.
 */
int XLALSimIMRSEOBNRv4ROMFrequencyOfTime(
  REAL8 *frequency,   /**< Output: Frequency (Hz) */
  REAL8 t,            /**< Time (s) at frequency */
  REAL8 m1SI,         /**< Mass of companion 1 (kg) */
  REAL8 m2SI,         /**< Mass of companion 2 (kg) */
  REAL8 chi1,         /**< Dimensionless aligned component spin 1 */
  REAL8 chi2          /**< Dimensionless aligned component spin 2 */
)
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    // Swap m1 and m2
    double m1temp = m1SI;
    double chi1temp = chi1;
    m1SI = m2SI;
    chi1 = chi2;
    m2SI = m1temp;
    chi2 = chi1temp;
  }

  // Set up phase spline
  gsl_spline *spline_phi;
  gsl_interp_accel *acc_phi;
  double Mf_final, Mtot_sec;
  double Mf_ROM_min, Mf_ROM_max;
  int ret = SEOBNRv4ROMTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final,
                                          &Mtot_sec, m1SI, m2SI, chi1, chi2,
                                          &Mf_ROM_min, &Mf_ROM_max);
  if(ret != 0)
    XLAL_ERROR(ret);

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  //XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  // Assume for now that we only care about f(t) *before* merger so that f(t) - f_ringdown >= 0.
  // Assume that we only need to cover the frequency range [f_min, f_ringdown/2].
  int N = 20;
  double log_f_pts[N];
  double log_t_pts[N];
  double log_f_min   = log(Mf_ROM_min * 1.000001); // raise minimum frequency slightly, so exp(log()) doesn't go below it
  double log_f_rng_2 = log(Mf_final/2.0);
  double dlog_f = (log_f_rng_2 - log_f_min) / (N-1);

  // Set up data in log-log space
  for (int i=0; i<N; i++) {
    log_f_pts[i] = log_f_rng_2 - i*dlog_f; // gsl likes the x-values to be monotonically increasing
    // Compute time relative to origin at merger
    double time_M = gsl_spline_eval_deriv(spline_phi, exp(log_f_pts[i]), acc_phi) / (2*LAL_PI) - t_corr;
    log_t_pts[i] = log(time_M * Mtot_sec);
  }

  // Check whether time is in bounds
  double t_rng_2 = exp(log_t_pts[0]);   // time of f_ringdown/2
  double t_min   = exp(log_t_pts[N-1]); // time of f_min
  if (t < t_rng_2 || t > t_min) {
    gsl_spline_free(spline_phi);
    gsl_interp_accel_free(acc_phi);
    XLAL_ERROR(XLAL_EDOM, "The frequency of time %g is outside allowed frequency range.\n", t);
  }

  // create new spline for data
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, N);
  gsl_spline_init(spline, log_t_pts, log_f_pts, N);

  *frequency = exp(gsl_spline_eval(spline, log(t), acc)) / Mtot_sec;

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_phi);

  return(XLAL_SUCCESS);
}


/** Setup SEOBNRv4ROM model using data files installed in $LAL_DATA_PATH
 */
UNUSED static void SEOBNRv4ROM_Init_LALDATA(void)
{
  if (SEOBNRv4ROM_IsSetup()) return;

  // Expect ROM datafile in a directory listed in LAL_DATA_PATH,
#ifdef LAL_HDF5_ENABLED
#define datafile ROMDataHDF5
  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv4ROM_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv4ROM data files in $LAL_DATA_PATH\n");
#else
  XLAL_ERROR_VOID(XLAL_EFAILED, "SEOBNRv4ROM requires HDF5 support which is not enabled\n");
#endif
}
