/*
 *  Copyright (C) 2014 Michael Puerrer, John Veitch
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

/**
 * \author Michael Puerrer, John Veitch
 *
 * \file
 *
 * \brief C code for SEOBNRv1 reduced order model (double spin version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv1 model.
 * Note that SEOBNRv2 supersedes SEOBNRv1.
 *
 * The binary data files are available at https://dcc.ligo.org/T1400701-v1.
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv1 model.
 *
 * @note Parameter ranges:
 *   * q <= 10
 *   * -1 <= chi_i <= 0.6
 *   * Mtot >= 12Msun
 *
 *  Aligned component spins chi1, chi2.
 *  Asymmetric mass-ratio q = max(m1/m2, m2/m1).
 *  Total mass Mtot.
 *
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

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


/********* Input data for spline basis points **************/
#define nk_amp 95  // number of SVD-modes == number of basis functions for amplitude
#define nk_phi 123 // number of SVD-modes == number of basis functions for phase

// Frequency points for amplitude and phase
static const double gA[] = {0.000631238, 0.000669113, 0.00070926, 0.000751815, 0.000796924, \
    0.000844739, 0.000895424, 0.000949149, 0.0010061, 0.00106646, \
    0.00113045, 0.00119828, 0.00127018, 0.00134639, 0.00142717, \
    0.0015128, 0.00160357, 0.00169978, 0.00180177, 0.00190987, \
    0.00202447, 0.00214594, 0.00227469, 0.00241117, 0.00255584, \
    0.00270919, 0.00287175, 0.00304405, 0.00322669, 0.00342029, \
    0.00362551, 0.00384304, 0.00407363, 0.00431804, 0.00457713, \
    0.00485175, 0.00514286, 0.00545143, 0.00577852, 0.00612523, \
    0.00649274, 0.00688231, 0.00729524, 0.00773296, 0.00819694, \
    0.00868875, 0.00921008, 0.00976268, 0.0103484, 0.0109693, 0.0116275, \
    0.0123252, 0.0130647, 0.0138486, 0.0146795, 0.0155602, 0.0164938, \
    0.0174835, 0.0185325, 0.0196444, 0.0208231, 0.0220725, 0.0233968, \
    0.0248006, 0.0262887, 0.027866, 0.029538, 0.0313102, 0.0331889, \
    0.0351802, 0.037291, 0.0395285, 0.0419002, 0.0444142, 0.047079, \
    0.0499038, 0.052898, 0.0560719, 0.0594362, 0.0630024, 0.0667825, \
    0.0707895, 0.0750368, 0.079539, 0.0843114, 0.0893701, 0.0947323, \
    0.100416, 0.106441, 0.112828, 0.119597, 0.126773, 0.13438, 0.142442, \
    0.15};

static const double gPhi[] = {0.00062519, 0.000638554, 0.000652301, 0.000666444, 0.000680997, \
    0.000695976, 0.000711395, 0.000727272, 0.000743622, 0.000760465, \
    0.000777818, 0.000795701, 0.000814135, 0.00083314, 0.000852738, \
    0.000872954, 0.000893812, 0.000915337, 0.000937555, 0.000960495, \
    0.000984187, 0.00100866, 0.00103395, 0.00106009, 0.00108711, \
    0.00111506, 0.00114396, 0.00117387, 0.00120483, 0.00123688, \
    0.00127008, 0.00130446, 0.00134009, 0.00137703, 0.00141533, \
    0.00145506, 0.00149628, 0.00153906, 0.00158349, 0.00162963, \
    0.00167757, 0.0017274, 0.00177922, 0.00183312, 0.00188921, \
    0.00194759, 0.0020084, 0.00207175, 0.00213777, 0.00220662, \
    0.00227844, 0.00235339, 0.00243165, 0.0025134, 0.00259883, \
    0.00268816, 0.0027816, 0.0028794, 0.00298181, 0.0030891, 0.00320158, \
    0.00331954, 0.00344334, 0.00357333, 0.00370991, 0.00385349, \
    0.00400452, 0.0041635, 0.00433095, 0.00450744, 0.00469358, \
    0.00489005, 0.00509755, 0.00531687, 0.00554887, 0.00579446, \
    0.00605465, 0.00633054, 0.00662331, 0.00693427, 0.00726485, \
    0.0076166, 0.00799125, 0.00839066, 0.00881692, 0.00927229, \
    0.00975928, 0.0102807, 0.0108395, 0.0114393, 0.0120836, 0.0127768, \
    0.0135236, 0.0143291, 0.0151992, 0.0161404, 0.0171602, 0.0182667, \
    0.0194694, 0.0207788, 0.0222069, 0.0237674, 0.0254758, 0.0273498, \
    0.0294099, 0.0316794, 0.0341854, 0.0369591, 0.0400368, 0.043461, \
    0.0472811, 0.0515553, 0.0563523, 0.0617535, 0.0678557, 0.0747749, \
    0.0826504, 0.0916509, 0.101981, 0.113893, 0.127695, 0.143771, 0.15};

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv1ROMDoubleSpin_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

typedef struct tagSEOBNRROMdataDS_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} SEOBNRROMdataDS_coeff;

struct tagSEOBNRROMdataDS
{
  UINT4 setup;
  gsl_vector* cvec_amp;
  gsl_vector* cvec_phi;
  gsl_matrix *Bamp;
  gsl_matrix *Bphi;
  gsl_vector* cvec_amp_pre;
};
typedef struct tagSEOBNRROMdataDS SEOBNRROMdataDS;

static SEOBNRROMdataDS __lalsim_SEOBNRv1ROMDS_data;

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  gsl_bspline_workspace *bwz;
  int ncx, ncy, ncz;
} SplineData;

/**************** Internal functions **********************/

static void SEOBNRv1ROMDoubleSpin_Init_LALDATA(void);
static int SEOBNRv1ROMDoubleSpin_Init(const char dir[]);
static bool SEOBNRv1ROMDoubleSpin_IsSetup(void);

static int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[]);
static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata);

static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff);
static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
static void SplineData_Destroy(SplineData *splinedata);
static void SplineData_Init(SplineData **splinedata);

static int read_vector(const char dir[], const char fname[], gsl_vector *v);
static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);

static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);

static int TP_Spline_interpolation_3d(
  REAL8 q,                  // Input: q-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
int SEOBNRv1ROMDoubleSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi1,
  double chi2,
  const REAL8Sequence *freqs_in, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
);


/********************* Definitions begin here ********************/


/** Setup SEOBNRv1ROMDoubleSpin model using data files installed in dir
 */
int SEOBNRv1ROMDoubleSpin_Init(const char dir[]) {
  if(__lalsim_SEOBNRv1ROMDS_data.setup) {
    XLALPrintError("Error: SEOBNRROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  SEOBNRROMdataDS_Init(&__lalsim_SEOBNRv1ROMDS_data, dir);

  if(__lalsim_SEOBNRv1ROMDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv1ROMDoubleSpin model has been initialised */
bool SEOBNRv1ROMDoubleSpin_IsSetup(void) {
  if(__lalsim_SEOBNRv1ROMDS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients
static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients as computed in Mathematica
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv1ROM_DS_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv1ROM_DS_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv1ROM_DS_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv1ROM_DS_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv1ROM_DS_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
}

static void SplineData_Init( SplineData **splinedata )
{
  if(!splinedata) exit(1);
  if(*splinedata) SplineData_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData));

  int ncx = 41+2;     // points in q
  int ncy = 21+2;     // points in chi1
  int ncz = 21+2;     // points in chi2
  (*splinedata)->ncx = ncx;
  (*splinedata)->ncy = ncy;
  (*splinedata)->ncz = ncz;

  // Set up B-spline basis for desired knots
  double qvec[] = {1., 1.125, 1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2., 2.25, 2.5, \
      2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 5.25, 5.5, 5.75, \
      6., 6.25, 6.5, 6.75, 7., 7.25, 7.5, 7.75, 8., 8.25, 8.5, 8.75, 9., \
      9.25, 9.5, 9.75, 10.};
  double chi1vec[] = {-1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.6, -0.5, -0.4, -0.3, \
      -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6};
  double chi2vec[] = {-1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.6, -0.5, -0.4, -0.3, \
      -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6};

  const size_t nbreak_x = ncx-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_y = ncy-2;  // must have nbreak = n-2 for cubic splines
  const size_t nbreak_z = ncz-2;  // must have nbreak = n-2 for cubic splines

  // allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bwx = gsl_bspline_alloc(4, nbreak_x);
  gsl_bspline_workspace *bwy = gsl_bspline_alloc(4, nbreak_y);
  gsl_bspline_workspace *bwz = gsl_bspline_alloc(4, nbreak_z);

  // set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_x = gsl_vector_alloc(nbreak_x);
  gsl_vector *breakpts_y = gsl_vector_alloc(nbreak_y);
  gsl_vector *breakpts_z = gsl_vector_alloc(nbreak_z);
  for (UINT4 i=0; i<nbreak_x; i++)
    gsl_vector_set(breakpts_x, i, qvec[i]);
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
  REAL8 q,                  // Input: q-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
) {

  SplineData *splinedata=NULL;
  SplineData_Init(&splinedata);
  gsl_bspline_workspace *bwx=splinedata->bwx;
  gsl_bspline_workspace *bwy=splinedata->bwy;
  gsl_bspline_workspace *bwz=splinedata->bwz;

  int ncx = splinedata->ncx; // points in q
  int ncy = splinedata->ncy; // points in chi1
  int ncz = splinedata->ncz; // points in chi2
  int N = ncx*ncy*ncz;  // size of the data matrix for one SVD-mode

  // Evaluate the TP spline for all SVD modes - amplitude
  for (int k=0; k<nk_amp; k++) { // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Tensor(&v, q, chi1, chi2, ncy, ncz, bwx, bwy, bwz);
    gsl_vector_set(c_amp, k, csum);
  }

  // Evaluate the TP spline for all SVD modes - phase
  for (int k=0; k<nk_phi; k++) {  // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Tensor(&v, q, chi1, chi2, ncy, ncz, bwx, bwy, bwz);
    gsl_vector_set(c_phi, k, csum);
  }

  // Evaluate the TP spline for the amplitude prefactor
  *amp_pre = Interpolate_Coefficent_Tensor(cvec_amp_pre, q, chi1, chi2, ncy, ncz, bwx, bwy, bwz);

  SplineData_Destroy(splinedata);

  return(0);
}

/* Set up a new ROM model, using data contained in dir */
int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[]) {
  // set up ROM
  int ncx = 41+2;       // points in q
  int ncy = 21+2;       // points in chi1
  int ncz = 21+2;       // points in chi2
  int N = ncx*ncy*ncz;  // size of the data matrix for one SVD-mode

  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup)
  {
    XLALPrintError("WARNING: You tried to setup the SEOBNRv1ROMDoubleSpin model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  gsl_set_error_handler(&err_handler);
  (romdata)->cvec_amp = gsl_vector_alloc(N*nk_amp);
  (romdata)->cvec_phi = gsl_vector_alloc(N*nk_phi);
  (romdata)->Bamp = gsl_matrix_alloc(nk_amp, nk_amp);
  (romdata)->Bphi = gsl_matrix_alloc(nk_phi, nk_phi);
  (romdata)->cvec_amp_pre = gsl_vector_alloc(N);
  ret=load_data(dir, (romdata)->cvec_amp, (romdata)->cvec_phi, (romdata)->Bamp, (romdata)->Bphi, (romdata)->cvec_amp_pre);

  if(XLAL_SUCCESS==ret) romdata->setup=1;
  else SEOBNRROMdataDS_Cleanup(romdata);

  return (ret);
}

/* Deallocate contents of the given SSEOBNRROMdata structure */
void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata) {
  if(romdata->cvec_amp) gsl_vector_free(romdata->cvec_amp);
  if(romdata->cvec_phi) gsl_vector_free(romdata->cvec_phi);
  if(romdata->Bamp) gsl_matrix_free(romdata->Bamp);
  if(romdata->Bphi) gsl_matrix_free(romdata->Bphi);
  if(romdata->cvec_amp_pre) gsl_vector_free(romdata->cvec_amp_pre);
  romdata->setup=0;
}

/* Structure for internal use */
static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff) {

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
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
int SEOBNRv1ROMDoubleSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi1,
  double chi2,
  const REAL8Sequence *freqs_in, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  )
{
  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv1ROMDS_data;
  if(*hptilde || *hctilde) {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (q < 1.0)     nudge(&q, 1.0, 1e-6);
  if (q > 10.0)    nudge(&q, 10.0, 1e-6);
  if (chi1 < -1.0) nudge(&chi1, -1.0, 1e-6);
  if (chi1 > 0.6)  nudge(&chi1, 0.6, 1e-6);
  if (chi2 < -1.0) nudge(&chi2, -1.0, 1e-6);
  if (chi2 > 0.6)  nudge(&chi2, 0.6, 1e-6);

  /* If either spin > 0.6, model not available, exit */
  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 0.6 || chi2 > 0.6 ) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1 or larger than 0.6!\nSEOBNRv1ROMDoubleSpin is only available for spins in the range -1 <= a/M <= 0.6.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (q > 10) {
    XLALPrintError( "XLAL Error - %s: q=%lf larger than 10!\nSEOBNRv1ROMDoubleSpin is only available for q in the range 1 <= q <= 10.\n", __func__, q);
    XLAL_ERROR( XLAL_EDOM );
  }

  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];

  if(fRef==0.0)
    fRef=fLow;

  /* Convert to geometric units for frequency */
  double Mf_ROM_min = fmax(gA[0], gPhi[0]);               // lowest allowed geometric frequency for ROM
  double Mf_ROM_max = fmin(gA[nk_amp-1], gPhi[nk_phi-1]); // highest allowed geometric frequency for ROM
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
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than starting frequency %g!\n", fHigh_geom, fLow_geom);
  if (fRef_geom > Mf_ROM_max)
    fRef_geom = Mf_ROM_max; // If fref > fhigh we reset fref to default value of cutoff frequency.
  if (fRef_geom < Mf_ROM_min) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for w.f. coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (q,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    q,                         // Input: q-value for which projection coefficients should be evaluated
    chi1,                      // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                      // Input: chi2-value for which projection coefficients should be evaluated
    romdata->cvec_amp,         // Input: data for spline coefficients for amplitude
    romdata->cvec_phi,         // Input: data for spline coefficients for phase
    romdata->cvec_amp_pre,     // Input: data for spline coefficients for amplitude prefactor
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode);
  }

  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // amp_pts = B_A^T . c_A
  // phi_pts = B_phi^T . c_phi
  gsl_vector* amp_f = gsl_vector_alloc(nk_amp);
  gsl_vector* phi_f = gsl_vector_alloc(nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, romdata->Bamp, romdata_coeff->c_amp, 0.0, amp_f);
  gsl_blas_dgemv(CblasTrans, 1.0, romdata->Bphi, romdata_coeff->c_phi, 0.0, phi_f);

  // Setup 1d splines in frequency
  gsl_interp_accel *acc_amp = gsl_interp_accel_alloc();
  gsl_spline *spline_amp = gsl_spline_alloc(gsl_interp_cspline, nk_amp);
  gsl_spline_init(spline_amp, gA, gsl_vector_const_ptr(amp_f,0), nk_amp);

  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_spline *spline_phi = gsl_spline_alloc(gsl_interp_cspline, nk_phi);
  gsl_spline_init(spline_phi, gPhi, gsl_vector_const_ptr(phi_f,0), nk_phi);


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
    UINT4 iStart = (UINT4) ceil(fLow_geom / deltaF_geom);
    UINT4 iStop = (UINT4) ceil(fHigh_geom / deltaF_geom);
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF_geom;

    offset = iStart;
  } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    npts = freqs_in->length;
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &tC, fLow, 0, &lalStrainUnit, npts);
    offset = 0;

    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i] * Mtot_sec;
  }


  if (!(*hptilde) || !(*hctilde))
  {
      XLALDestroyREAL8Sequence(freqs);

      gsl_spline_free(spline_amp);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_amp);
      gsl_interp_accel_free(acc_phi);
      gsl_vector_free(amp_f);
      gsl_vector_free(phi_f);
      SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
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

  REAL8 s = 0.5; // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv1
  double Mtot = Mtot_sec / LAL_MTSUN_SI;
  double amp0 = Mtot * amp_pre * Mtot_sec * LAL_MRSUN_SI / (distance); // Correct overall amplitude to undo mass-dependent scaling used in ROM

  // Evaluate reference phase for setting phiRef correctly
  double phase_change = gsl_spline_eval(spline_phi, fRef_geom, acc_phi) - 2*phiRef;

  // Assemble waveform from aplitude and phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i];
    if (f > Mf_ROM_max) continue; // We're beyond the highest allowed frequency; since freqs may not be ordered, we'll just skip the current frequency and leave zero in the buffer
    int j = i + offset; // shift index for frequency series if needed
    double A = gsl_spline_eval(spline_amp, f, acc_amp);
    double phase = gsl_spline_eval(spline_phi, f, acc_phi) - phase_change;
    COMPLEX16 htilde = s*amp0*A * cexp(I*phase);
    pdata[j] =      pcoef * htilde;
    cdata[j] = -I * ccoef * htilde;
  }
  
  /* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */

  // Get SEOBNRv1 ringdown frequency for 22 mode
  // XLALSimInspiralGetFinalFreq wants masses in SI units, so unfortunately we need to convert back
  double Mtot_SI = Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * 1.0/(1.0+q);
  double m2_SI = Mtot_SI * q/(1.0+q);
  double Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi1, 0,0,chi2, SEOBNRv1) * Mtot_sec;

  UINT4 L = freqs->length;
  // prevent gsl interpolation errors
  if (Mf_final > freqs->data[L-1])
    Mf_final = freqs->data[L-1];
  if (Mf_final < freqs->data[0])
  {
      XLALDestroyREAL8Sequence(freqs);

      gsl_spline_free(spline_amp);
      gsl_spline_free(spline_phi);
      gsl_interp_accel_free(acc_amp);
      gsl_interp_accel_free(acc_phi);
      gsl_vector_free(amp_f);
      gsl_vector_free(phi_f);
      SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
      XLAL_ERROR(XLAL_EDOM, "f_ringdown < f_min");
  }

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // We compute the dimensionless time correction t/M since we use geometric units.
  REAL8 t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI);
  XLAL_PRINT_INFO("t_corr [s] = %g\n", t_corr * Mtot_sec);

  // Now correct phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i];
    int j = i + offset; // shift index for frequency series if needed
    pdata[j] *= cexp(-2*LAL_PI * I * f * t_corr);
    cdata[j] *= cexp(-2*LAL_PI * I * f * t_corr);
  }

  XLALDestroyREAL8Sequence(freqs);

  gsl_spline_free(spline_amp);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_amp);
  gsl_interp_accel_free(acc_phi);
  gsl_vector_free(amp_f);
  gsl_vector_free(phi_f);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}

/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv1_ROM_DoubleSpin model.
 *
 * XLALSimIMRSEOBNRv1ROMDoubleSpin() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv1ROMDoubleSpinFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv1ROMDoubleSpinFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv1ROMDoubleSpinFrequencySequence(
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
  REAL8 chi2)                                   /**< Dimensionless aligned component spin 2 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI) {
    /* Swap m1 and m2 */
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
  double q = mass2 / mass1;
  if(q<1.0) q=1./q;
  /* Total mass in seconds */
  double Mtot_sec = Mtot * LAL_MTSUN_SI;

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv1ROMDoubleSpin_is_initialized, SEOBNRv1ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv1ROMDoubleSpin_Init_LALDATA();
#endif

  if(!SEOBNRv1ROMDoubleSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv1ROMDoubleSpin - check your $LAL_DATA_PATH\n");

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv1ROMDoubleSpinCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, q, chi1, chi2, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv1_ROM_DoubleSpin model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv1ROMDoubleSpin(
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
  REAL8 chi2)                                   /**< Dimensionless aligned component spin 2 */
{
  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1SI < m2SI)
  {
    /* Swap m1 and m2 */
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
  double q = mass2 / mass1;
  if(q<1.0) q=1./q;
  /* Total mass in seconds */
  double Mtot_sec = Mtot * LAL_MTSUN_SI;

  if(fRef==0.0)
    fRef=fLow;

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv1ROMDoubleSpin_is_initialized, SEOBNRv1ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv1ROMDoubleSpin_Init_LALDATA();
#endif

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv1ROMDoubleSpinCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, q, chi1, chi2, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** Setup SEOBNRv1ROMDoubleSpin model using data files installed in $LAL_DATA_PATH
 */
void SEOBNRv1ROMDoubleSpin_Init_LALDATA(void)
{
  if (SEOBNRv1ROMDoubleSpin_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "SEOBNRv1ROM_DS_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv1ROMDoubleSpin_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv1ROMDoubleSpin data files in $LAL_DATA_PATH\n");
}
