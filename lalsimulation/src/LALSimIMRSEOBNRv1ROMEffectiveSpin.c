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
#define nk_amp 78  // number of SVD-modes == number of basis functions for amplitude
#define nk_phi 200 // number of SVD-modes == number of basis functions for phase

static const double gA[] = {0.0001, 0.00011, 0.000121, 0.0001331, 0.00014641, 0.000161051, \
  0.000177156, 0.000194872, 0.000214359, 0.000235795, 0.000259374, \
  0.000285312, 0.000313843, 0.000345227, 0.00037975, 0.000417725, \
  0.000459497, 0.000505447, 0.000555992, 0.000611591, 0.00067275, \
  0.000740025, 0.000814027, 0.00089543, 0.000984973, 0.00108347, \
  0.00119182, 0.001311, 0.0014421, 0.00158631, 0.00174494, 0.00191943, \
  0.00211138, 0.00232252, 0.00255477, 0.00281024, 0.00309127, \
  0.00340039, 0.00374043, 0.00411448, 0.00452593, 0.00497852, \
  0.00547637, 0.00602401, 0.00662641, 0.00728905, 0.00801795, \
  0.00881975, 0.00970172, 0.0106719, 0.0117391, 0.012913, 0.0142043, \
  0.0156247, 0.0171872, 0.0189059, 0.0207965, 0.0228762, 0.0251638, \
  0.0276801, 0.0304482, 0.033493, 0.0368423, 0.0405265, 0.0445792, \
  0.0490371, 0.0539408, 0.0593349, 0.0652683, 0.0717952, 0.0789747, \
  0.0868722, 0.0955594, 0.105115, 0.115627, 0.12719, 0.139908, 0.14};

static const double gPhi[] = {0.0001, 0.000101411, 0.000102849, 0.000104314, 0.000105806, \
  0.000107328, 0.000108878, 0.000110459, 0.00011207, 0.000113712, \
  0.000115387, 0.000117095, 0.000118836, 0.000120613, 0.000122424, \
  0.000124272, 0.000126157, 0.000128081, 0.000130044, 0.000132047, \
  0.000134091, 0.000136177, 0.000138307, 0.000140481, 0.000142701, \
  0.000144968, 0.000147283, 0.000149648, 0.000152063, 0.000154531, \
  0.000157052, 0.000159627, 0.00016226, 0.00016495, 0.0001677, \
  0.000170512, 0.000173386, 0.000176325, 0.000179331, 0.000182405, \
  0.00018555, 0.000188768, 0.000192059, 0.000195428, 0.000198876, \
  0.000202405, 0.000206017, 0.000209716, 0.000213504, 0.000217383, \
  0.000221357, 0.000225428, 0.000229598, 0.000233872, 0.000238253, \
  0.000242743, 0.000247346, 0.000252066, 0.000256907, 0.000261871, \
  0.000266965, 0.00027219, 0.000277553, 0.000283057, 0.000288707, \
  0.000294507, 0.000300464, 0.000306582, 0.000312866, 0.000319323, \
  0.000325958, 0.000332778, 0.000339788, 0.000346996, 0.000354409, \
  0.000362034, 0.000369878, 0.000377949, 0.000386256, 0.000394808, \
  0.000403612, 0.00041268, 0.00042202, 0.000431643, 0.00044156, \
  0.000451782, 0.000462321, 0.000473188, 0.000484398, 0.000495963, \
  0.000507897, 0.000520216, 0.000532935, 0.00054607, 0.000559639, \
  0.000573659, 0.000588149, 0.00060313, 0.000618621, 0.000634645, \
  0.000651225, 0.000668385, 0.00068615, 0.000704548, 0.000723607, \
  0.000743356, 0.000763826, 0.000785052, 0.000807068, 0.000829911, \
  0.000853621, 0.000878237, 0.000903805, 0.000930369, 0.00095798, \
  0.000986689, 0.00101655, 0.00104762, 0.00107997, 0.00111365, \
  0.00114874, 0.00118532, 0.00122345, 0.00126323, 0.00130475, \
  0.00134809, 0.00139336, 0.00144067, 0.00149013, 0.00154188, \
  0.00159603, 0.00165273, 0.00171213, 0.0017744, 0.0018397, 0.00190823, \
  0.00198018, 0.00205578, 0.00213524, 0.00221883, 0.00230681, \
  0.00239947, 0.00249712, 0.00260011, 0.0027088, 0.00282359, \
  0.00294492, 0.00307324, 0.00320907, 0.00335297, 0.00350553, \
  0.00366741, 0.00383935, 0.00402211, 0.00421656, 0.00442364, \
  0.0046444, 0.00487997, 0.0051316, 0.00540067, 0.00568873, 0.00599744, \
  0.0063287, 0.00668457, 0.00706737, 0.00747967, 0.00792436, \
  0.00840463, 0.00892411, 0.00948683, 0.0100974, 0.0107608, 0.011483, \
  0.0122706, 0.013131, 0.0140727, 0.0151056, 0.0162408, 0.0174911, \
  0.0188713, 0.0203987, 0.0220931, 0.0239776, 0.0260796, 0.0284307, \
  0.0310686, 0.0340377, 0.0373912, 0.0411922, 0.0455169, 0.0504574, \
  0.0561255, 0.0626582, 0.0702238, 0.0790313, 0.0893416, 0.101483, \
  0.115873, 0.133047, 0.14};

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv1ROMEffectiveSpin_is_initialized = PTHREAD_ONCE_INIT;
#endif

/*************** type definitions ******************/

typedef struct tagSEOBNRROMdata_coeff
{
  gsl_vector* c_amp;
  gsl_vector* c_phi;
} SEOBNRROMdata_coeff;

struct tagSEOBNRROMdata
{
  UINT4 setup;
  gsl_vector* cvec_amp;
  gsl_vector* cvec_phi;
  gsl_matrix *Bamp;
  gsl_matrix *Bphi;
  gsl_vector* cvec_amp_pre;
};
typedef struct tagSEOBNRROMdata SEOBNRROMdata;

static SEOBNRROMdata __lalsim_SEOBNRv1ROMSS_data;

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  int ncx, ncy;
} SplineData;

/**************** Internal functions **********************/

static void SEOBNRv1ROMEffectiveSpin_Init_LALDATA(void);
static int SEOBNRv1ROMEffectiveSpin_Init(const char dir[]);
static bool SEOBNRv1ROMEffectiveSpin_IsSetup(void);

static int SEOBNRROMdata_Init(SEOBNRROMdata *romdata, const char dir[]);
static void SEOBNRROMdata_Cleanup(SEOBNRROMdata *romdata);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv1ROMEffectiveSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi,
  const REAL8Sequence *freqs, /* Frequency points at which to evaluate the waveform (Hz) */
  double deltaF
  /* If deltaF >0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
);

static void SEOBNRROMdata_coeff_Init(SEOBNRROMdata_coeff **romdatacoeff);
static void SEOBNRROMdata_coeff_Cleanup(SEOBNRROMdata_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
static void SplineData_Destroy(SplineData *splinedata);
static void SplineData_Init(SplineData **splinedata);

static int read_vector(const char dir[], const char fname[], gsl_vector *v);
static int read_matrix(const char dir[], const char fname[], gsl_matrix *m);

static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);

static int TP_Spline_interpolation_2d(
  REAL8 q,                  // Input: q-value for which projection coefficients should be evaluated
  REAL8 chi,                // Input: chi-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
);


/********************* Definitions begin here ********************/

/** Setup SEOBNRv1ROMEffectiveSpin model using data files installed in dir
 */
static int SEOBNRv1ROMEffectiveSpin_Init(const char dir[]) {
  if(__lalsim_SEOBNRv1ROMSS_data.setup) {
    XLALPrintError("Error: SEOBNRROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  SEOBNRROMdata_Init(&__lalsim_SEOBNRv1ROMSS_data, dir);

  if(__lalsim_SEOBNRv1ROMSS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv1ROMEffectiveSpin model has been initialised */
static bool SEOBNRv1ROMEffectiveSpin_IsSetup(void) {
  if(__lalsim_SEOBNRv1ROMSS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients
static int load_data(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients as computed in Mathematica
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv1ROM_SS_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv1ROM_SS_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv1ROM_SS_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv1ROM_SS_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv1ROM_SS_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
}

static void SplineData_Init( SplineData **splinedata )
{
  if(!splinedata) exit(1);
  if(*splinedata) SplineData_Destroy(*splinedata);

  (*splinedata)=XLALCalloc(1,sizeof(SplineData));

  const int ncx = 159;    // points in q
  const int ncy = 49;     // points in chi
  (*splinedata)->ncx = ncx;
  (*splinedata)->ncy = ncy;

  // Set up B-spline basis for desired knots
  double qvec[] = {1., 1.125, 1.25, 1.375, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, \
    3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, \
    9., 9.5, 10., 10.5, 11., 11.5, 12., 12.5, 13., 13.5, 14., 14.5, 15., \
    15.5, 16., 16.5, 17., 17.5, 18., 18.5, 19., 19.5, 20., 20.5, 21., \
    21.5, 22., 22.5, 23., 23.5, 24., 24.5, 25., 25.5, 26., 26.5, 27., \
    27.5, 28., 28.5, 29., 29.5, 30., 30.5, 31., 31.5, 32., 32.5, 33., \
    33.5, 34., 34.5, 35., 35.5, 36., 36.5, 37., 37.5, 38., 38.5, 39., \
    39.5, 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51., \
    52., 53., 54., 55., 56., 57., 58., 59., 60., 61., 62., 63., 64., 65., \
    66., 67., 68., 69., 70., 71., 72., 73., 74., 75., 76., 77., 78., 79., \
    80., 81., 82., 83., 84., 85., 86., 87., 88., 89., 90., 91., 92., 93., \
    94., 95., 95.5, 96., 96.5, 97., 97.5, 98., 98.5, 98.75, 99., 99.25, \
    99.5, 99.75, 100.};

  double chivec[] = {-1., -0.975, -0.95, -0.925, -0.9, -0.875, -0.85, -0.825, -0.8, \
    -0.775, -0.75, -0.725, -0.7, -0.675, -0.65, -0.625, -0.6, -0.55, \
    -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., \
    0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, \
    0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6};

  const size_t nbreak_x = ncx-2;  // must have nbreak = n -2 for cubic splines
  const size_t nbreak_y = ncy-2;  // must have nbreak = n -2 for cubic splines

  // allocate a cubic bspline workspace (k = 4)
  gsl_bspline_workspace *bwx = gsl_bspline_alloc(4, nbreak_x);
  gsl_bspline_workspace *bwy = gsl_bspline_alloc(4, nbreak_y);

  // set breakpoints (and thus knots by hand)
  gsl_vector *breakpts_x = gsl_vector_alloc(nbreak_x);
  gsl_vector *breakpts_y = gsl_vector_alloc(nbreak_y);
  for (UINT4 i=0; i<nbreak_x; i++)
    gsl_vector_set(breakpts_x, i, qvec[i]);
  for (UINT4 j=0; j<nbreak_y; j++)
    gsl_vector_set(breakpts_y, j, chivec[j]);

  gsl_bspline_knots(breakpts_x, bwx);
  gsl_bspline_knots(breakpts_y, bwy);

  gsl_vector_free(breakpts_x);
  gsl_vector_free(breakpts_y);

  (*splinedata)->bwx=bwx;
  (*splinedata)->bwy=bwy;
}

static void SplineData_Destroy(SplineData *splinedata)
{
  if(!splinedata) return;
  if(splinedata->bwx) gsl_bspline_free(splinedata->bwx);
  if(splinedata->bwy) gsl_bspline_free(splinedata->bwy);
  XLALFree(splinedata);
}

// Interpolate projection coefficients for amplitude and phase over the parameter space (q, chi).
// The multi-dimensional interpolation is carried out via a tensor product decomposition.
static int TP_Spline_interpolation_2d(
  REAL8 q,                  // Input: q-value for which projection coefficients should be evaluated
  REAL8 chi,                // Input: chi-value for which projection coefficients should be evaluated
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

  int ncx = splinedata->ncx; // points in q
  int ncy = splinedata->ncy; // points in chi
  int N = ncx*ncy;  // size of the data matrix for one SVD-mode

  // Evaluate the TP spline for all SVD modes - amplitude
  for (int k=0; k<nk_amp; k++) { // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_amp, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Matrix(&v, q, chi, ncx, ncy, bwx, bwy);
    gsl_vector_set(c_amp, k, csum);
  }

  // Evaluate the TP spline for all SVD modes - phase
  for (int k=0; k<nk_phi; k++) {  // For each SVD mode
    gsl_vector v = gsl_vector_subvector(cvec_phi, k*N, N).vector; // Pick out the coefficient matrix corresponding to the k-th SVD mode.
    REAL8 csum = Interpolate_Coefficent_Matrix(&v, q, chi, ncx, ncy, bwx, bwy);
    gsl_vector_set(c_phi, k, csum);
  }

  // Evaluate the TP spline for the amplitude prefactor
  *amp_pre = Interpolate_Coefficent_Matrix(cvec_amp_pre, q, chi, ncx, ncy, bwx, bwy);

  SplineData_Destroy(splinedata);

  return(0);
}

/* Set up a new ROM model, using data contained in dir */
static int SEOBNRROMdata_Init(SEOBNRROMdata *romdata, const char dir[]) {
  // set up ROM
  int ncx = 159;    // points in q
  int ncy = 49;     // points in chi
  int N = ncx*ncy;  // size of the data matrix for one SVD-mode

  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup)
  {
    XLALPrintError("WARNING: You tried to setup the SEOBNRv1ROMEffectiveSpin model that was already initialised. Ignoring\n");
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
  else SEOBNRROMdata_Cleanup(romdata);

  return (ret);
}


/* Deallocate contents of the given SEOBNRROMdata structure */
static void SEOBNRROMdata_Cleanup(SEOBNRROMdata *romdata) {
  if(romdata->cvec_amp) gsl_vector_free(romdata->cvec_amp);
  if(romdata->cvec_phi) gsl_vector_free(romdata->cvec_phi);
  if(romdata->Bamp) gsl_matrix_free(romdata->Bamp);
  if(romdata->Bphi) gsl_matrix_free(romdata->Bphi);
  if(romdata->cvec_amp_pre) gsl_vector_free(romdata->cvec_amp_pre);
  romdata->setup=0;
}

/* Structure for internal use */
static void SEOBNRROMdata_coeff_Init(SEOBNRROMdata_coeff **romdatacoeff) {

  if(!romdatacoeff) exit(1);
  /* Create storage for structures */
  if(!*romdatacoeff)
    *romdatacoeff=XLALCalloc(1,sizeof(SEOBNRROMdata_coeff));
  else
    SEOBNRROMdata_coeff_Cleanup(*romdatacoeff);

  (*romdatacoeff)->c_amp = gsl_vector_alloc(nk_amp);
  (*romdatacoeff)->c_phi = gsl_vector_alloc(nk_phi);
}

/* Deallocate contents of the given SEOBNRROMdata_coeff structure */
static void SEOBNRROMdata_coeff_Cleanup(SEOBNRROMdata_coeff *romdatacoeff) {
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
static int SEOBNRv1ROMEffectiveSpinCore(
  COMPLEX16FrequencySeries **hptilde,
  COMPLEX16FrequencySeries **hctilde,
  double phiRef,
  double fRef,
  double distance,
  double inclination,
  double Mtot_sec,
  double q,
  double chi,
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
  SEOBNRROMdata *romdata=&__lalsim_SEOBNRv1ROMSS_data;
  if(*hptilde || *hctilde) {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (q < 1.0)    nudge(&q, 1.0, 1e-6);
  if (q > 100.0)  nudge(&q, 100.0, 1e-6);
  if (chi < -1.0) nudge(&chi, -1.0, 1e-6);
  if (chi > 0.6)  nudge(&chi, 0.6, 1e-6);

  /* If either spin > 0.6, model not available, exit */
  if ( chi < -1.0 || chi > 0.6 ) {
    XLALPrintError( "XLAL Error - %s: chi smaller than -1 or larger than 0.6!\nSEOBNRv1ROMEffectiveSpin is only available for spins in the range -1 <= a/M <= 0.6.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (q > 100) {
    XLALPrintError( "XLAL Error - %s: q=%lf larger than 100!\nSEOBNRv1ROMEffectiveSpin is only available for q in the range 1 <= q <= 100.\n", __func__,q);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (q >= 20 && q <= 40 && chi < -0.75 && chi > -0.9) {
    XLALPrintWarning( "XLAL Warning - %s: q in [20,40] and chi in [-0.8]. The SEOBNRv1 model is not trustworthy in this region!\nSee Fig 15 in CQG 31 195010, 2014 for details.", __func__);
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
  if (fHigh_geom == 0)
    fHigh_geom = Mf_ROM_max;
  else if (fHigh_geom > Mf_ROM_max) {
	  XLALPrintWarning("Maximal frequency Mf_high=%g is greater than highest ROM frequency Mf_ROM_Max=%g. Using Mf_high=Mf_ROM_Max.", fHigh_geom, Mf_ROM_max);
	  fHigh_geom = Mf_ROM_max;
  }
  else if (fHigh_geom < Mf_ROM_min)
    XLAL_ERROR(XLAL_EDOM, "End frequency %g is smaller than starting frequency %g!\n", fHigh_geom, fLow_geom);
  if (fRef_geom > Mf_ROM_max) {
	  XLALPrintWarning("Reference frequency Mf_ref=%g is greater than maximal frequency in ROM Mf=%g. Starting at maximal frequency in ROM.\n", fRef_geom, Mf_ROM_max);
    fRef_geom = Mf_ROM_max; // If fref > fhigh we reset fref to default value of cutoff frequency.
  }
  if (fRef_geom < Mf_ROM_min) {
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fLow_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for w.f. coefficiencts */
  SEOBNRROMdata_coeff *romdata_coeff=NULL;
  SEOBNRROMdata_coeff_Init(&romdata_coeff);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (q,chi) */
  retcode=TP_Spline_interpolation_2d(
    q,                         // Input: q-value for which projection coefficients should be evaluated
    chi,                       // Input: chi-value for which projection coefficients should be evaluated
    romdata->cvec_amp,         // Input: data for spline coefficients for amplitude
    romdata->cvec_phi,         // Input: data for spline coefficients for phase
    romdata->cvec_amp_pre,     // Input: data for spline coefficients for amplitude prefactor
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode, "Parameter-space interpolation failed.");
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
      SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
      XLAL_ERROR(XLAL_EFUNC, "Waveform allocation failed.");
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

  REAL8 s = 1.0/sqrt(2.0); // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv1
  double Mtot = Mtot_sec / LAL_MTSUN_SI;
  double amp0 = Mtot * amp_pre * Mtot_sec * LAL_MRSUN_SI / (distance); // Correct overall amplitude to undo mass-dependent scaling used in single-spin ROM

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
  double Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi, 0,0,chi, SEOBNRv1) * Mtot_sec;

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
      SEOBNRROMdata_coeff_Cleanup(romdata_coeff);
      XLAL_ERROR(XLAL_EDOM, "f_ringdown < f_min");
  }

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // We compute the dimensionless time correction t/M since we use geometric units.
  REAL8 t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI);

  // Now correct phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double f = freqs->data[i] - fRef_geom;
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
  SEOBNRROMdata_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRROM_c
 *
 * @brief Functions for producing SEOBNRv1 and v2 waveforms
 * using reduced order models.
 *
 * @review SEOBNRv1/2_ROM_(Effective/Double)Spin reviewed by Frank Ohme, Sarah Caudill, Michael Puerrer, Ian Harry, John Veitch, Gareth Thomas. Review concluded with git hash 9dc5e84583bfe2707ac20638e7b89bf988d4d482 (July 2015).
 *
 * @{
 *
 * @name SEOBNRv1 Reduced Order Model (Effective Spin)
 *
 * @author Michael Puerrer, John Veitch
 *
 * @brief C code for SEOBNRv1 reduced order model (equal spin version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv1 model with equal spins.
 * Note that SEOBNRv2 supersedes SEOBNRv1.
 *
 * The binary data files are available at https://dcc.ligo.org/T1400701-v1.
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv1 model.
 *
 * @note Due to non-smoothness in SEOBNRv1 at chi1=chi2 ~ -0.8 and 20 <= q <= 40 the
 * ROM deviates from the SEOBNRv1 behavior there. See arXiv:1402.4146, Fig 7,
 * Fig 11, and Fig 13 for details.
 *
 * @note Parameter ranges:
 *   * 1 <= q <= 100
 *   * -1 <= chi <= 0.6
 *   * Mtot >= 1.4Msun
 *
 *  Equal spin chi = chi1 = chi2.
 *  Asymmetric mass-ratio q = max(m1/m2, m2/m1).
 *  Total mass Mtot.
 *
 * @{
 */

/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv1_ROM_EffectiveSpin model.
 *
 * XLALSimIMRSEOBNRv1ROMEffectiveSpin() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv1ROMEffectiveSpinFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv1ROMEffectiveSpinFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv1ROMEffectiveSpinFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz), need to be strictly monotonically increasing */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi)                                    /**< Effective aligned spin */
{

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double q = mass2 / mass1;
  if(q<1.0) q=1./q;
  /* Total mass in seconds */
  double Mtot_sec = Mtot * LAL_MTSUN_SI;

  // Load ROM data if not already loaded
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv1ROMEffectiveSpin_is_initialized, SEOBNRv1ROMEffectiveSpin_Init_LALDATA);
#else
  SEOBNRv1ROMEffectiveSpin_Init_LALDATA();
#endif

  if(!SEOBNRv1ROMEffectiveSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv1ROMEffectiveSpin model - check your $LAL_DATA_PATH\n");

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv1ROMEffectiveSpinCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, q, chi, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv1_ROM_EffectiveSpin model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency fHigh to the next power of 2 in
 * the size of the frequency series.
 */
int XLALSimIMRSEOBNRv1ROMEffectiveSpin(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi)                                    /**< Effective aligned spin */
{

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

  // Load ROM data if not already loaded
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv1ROMEffectiveSpin_is_initialized, SEOBNRv1ROMEffectiveSpin_Init_LALDATA);
#else
  SEOBNRv1ROMEffectiveSpin_Init_LALDATA();
#endif

  if(!SEOBNRv1ROMEffectiveSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv1ROMEffectiveSpin data - check your $LAL_DATA_PATH\n");

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv1ROMEffectiveSpinCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, q, chi, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */
/** @} */

/** Setup SEOBNRv1ROMEffectiveSpin model using data files installed in $LAL_DATA_PATH
 */
static void SEOBNRv1ROMEffectiveSpin_Init_LALDATA(void)
{
  if (SEOBNRv1ROMEffectiveSpin_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "SEOBNRv1ROM_SS_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv1ROMEffectiveSpin_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv1ROMEffectiveSpin data files in $LAL_DATA_PATH\n");
}
