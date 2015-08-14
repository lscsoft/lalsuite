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


/******************************************************************
 * The double-spin SEOBNRv2 ROM consists of a number of submodels *
 * for subdomains of the whole parameter space:                   *
 * These submodels may change in the future.                      *
 *                                                                *
 * "Core" submodel 1                                              *
 * B-spline points: 54x24x24                                      *
 * Frequency points: {133, 139}                                   *
 *                                                                *
 * "Near q=1" submodel 2                                          *
 * B-spline points: 11x60x60                                      *
 * Frequency points: {200, 187}                                   *
 *                                                                *
 * "High q, high chi1" submodel 3                                 *
 * B-spline points: 24x30x24                                      *
 * frequency points: {133, 139}                                   *
 *                                                                *
 *****************************************************************/


/********* Input data for spline basis points **************/

// The frequency points are not the same for all submodels

#define nk_amp_sub1 133  // number of SVD-modes == number of basis functions for amplitude
#define nk_phi_sub1 139  // number of SVD-modes == number of basis functions for phase

// Frequency points for amplitude and phase
static const double gA_sub1[] = {0.0005, 0.000525, 0.00055125, 0.000578812, 0.000607753, 0.000638141, 0.000670048,
  0.00070355, 0.000738728, 0.000775664, 0.000814447, 0.00085517, 0.000897928, 0.000942825,
  0.000989966, 0.00103946, 0.00109144, 0.00114601, 0.00120331, 0.00126348, 0.00132665,
  0.00139298, 0.00146263, 0.00153576, 0.00161255, 0.00169318, 0.00177784, 0.00186673,
  0.00196006, 0.00205807, 0.00216097, 0.00226902, 0.00238247, 0.00250159, 0.00262667,
  0.00275801, 0.00289591, 0.0030407, 0.00319274, 0.00335238, 0.00351999, 0.00369599, 0.00388079,
  0.00407483, 0.00427858, 0.0044925, 0.00471713, 0.00495299, 0.00520063, 0.00546067, 0.0057337,
  0.00602038, 0.0063214, 0.00663747, 0.00696935, 0.00731782, 0.00768371, 0.00806789, 0.00847129,
  0.00889485, 0.00933959, 0.00980657, 0.0102969, 0.0108117, 0.0113523, 0.01192, 0.0125159,
  0.0131417, 0.0137988, 0.0144888, 0.0152132, 0.0159739, 0.0167726, 0.0176112, 0.0184918,
  0.0194163, 0.0203872, 0.0214065, 0.0224768, 0.0236007, 0.0247807, 0.0260198, 0.0273207,
  0.0286868, 0.0301211, 0.0316272, 0.0332085, 0.034869, 0.0366124, 0.038443, 0.0403652,
  0.0423834, 0.0445026, 0.0467277, 0.0490641, 0.0515173, 0.0540932, 0.0567979, 0.0596378,
  0.0626196, 0.0657506, 0.0690382, 0.0724901, 0.0761146, 0.0799203, 0.0839163, 0.0881121,
  0.0925177, 0.0971436, 0.102001, 0.107101, 0.112456, 0.118079, 0.123983, 0.130182, 0.136691,
  0.143525, 0.150702, 0.158237, 0.166149, 0.174456, 0.183179, 0.192338, 0.201955, 0.212052,
  0.222655, 0.233788, 0.245477, 0.257751, 0.270639, 0.28417, 0.298379, 0.3};

static const double gPhi_sub1[] = {0.0005, 0.000509921, 0.000520106, 0.000530563, 0.000541301, 0.000552329, 0.000563659,
  0.000575299, 0.000587261, 0.000599555, 0.000612194, 0.00062519, 0.000638554, 0.000652301,
  0.000666444, 0.000680997, 0.000695976, 0.000711395, 0.000727272, 0.000743622, 0.000760465,
  0.000777818, 0.000795701, 0.000814135, 0.00083314, 0.000852738, 0.000872954, 0.000893812,
  0.000915337, 0.000937555, 0.000960495, 0.000984187, 0.00100866, 0.00103395, 0.00106009,
  0.00108711, 0.00111506, 0.00114396, 0.00117387, 0.00120483, 0.00123688, 0.00127008,
  0.00130446, 0.00134009, 0.00137703, 0.00141533, 0.00145506, 0.00149628, 0.00153906,
  0.00158349, 0.00162963, 0.00167757, 0.0017274, 0.00177922, 0.00183312, 0.00188921, 0.00194759,
  0.0020084, 0.00207175, 0.00213777, 0.00220662, 0.00227844, 0.00235339, 0.00243165, 0.0025134,
  0.00259883, 0.00268816, 0.0027816, 0.0028794, 0.00298181, 0.0030891, 0.00320158, 0.00331954,
  0.00344334, 0.00357333, 0.00370991, 0.00385349, 0.00400452, 0.0041635, 0.00433095, 0.00450744,
  0.00469358, 0.00489005, 0.00509755, 0.00531687, 0.00554887, 0.00579446, 0.00605465,
  0.00633054, 0.00662331, 0.00693427, 0.00726485, 0.0076166, 0.00799125, 0.00839066, 0.00881692,
  0.00927229, 0.00975928, 0.0102807, 0.0108395, 0.0114393, 0.0120836, 0.0127768, 0.0135236,
  0.0143291, 0.0151992, 0.0161404, 0.0171602, 0.0182667, 0.0194694, 0.0207788, 0.0222069,
  0.0237674, 0.0254758, 0.0273498, 0.0294099, 0.0316794, 0.0341854, 0.0369591, 0.0400368,
  0.043461, 0.0472811, 0.0515553, 0.0563523, 0.0617535, 0.0678557, 0.0747749, 0.0826504,
  0.0916509, 0.101981, 0.113893, 0.127695, 0.143771, 0.1626, 0.184787, 0.2111, 0.242523,
  0.280334, 0.3};

#define nk_amp_sub3 nk_amp_sub1
#define nk_phi_sub3 nk_phi_sub1
#define gA_sub3 gA_sub1
#define gPhi_sub3 gPhi_sub1

#define nk_amp_sub2 200
#define nk_phi_sub2 187

static const double gA_sub2[] = {0.0005, 0.00051635, 0.000533235, 0.000550671, 0.000568678, \
0.000587274, 0.000606478, 0.00062631, 0.00064679, 0.00066794, \
0.000689782, 0.000712338, 0.000735631, 0.000759686, 0.000784528, \
0.000810182, 0.000836675, 0.000864034, 0.000892288, 0.000921466, \
0.000951598, 0.000982715, 0.00101485, 0.00104804, 0.00108231, \
0.0011177, 0.00115425, 0.00119199, 0.00123097, 0.00127122, \
0.00131279, 0.00135572, 0.00140005, 0.00144583, 0.00149311, \
0.00154194, 0.00159236, 0.00164443, 0.0016982, 0.00175373, \
0.00181108, 0.0018703, 0.00193146, 0.00199462, 0.00205984, 0.0021272, \
0.00219676, 0.00226859, 0.00234277, 0.00241938, 0.0024985, 0.0025802, \
0.00266457, 0.0027517, 0.00284168, 0.00293461, 0.00303057, \
0.00312967, 0.00323201, 0.00333769, 0.00344684, 0.00355955, \
0.00367594, 0.00379615, 0.00392028, 0.00404848, 0.00418086, \
0.00431757, 0.00445876, 0.00460456, 0.00475513, 0.00491062, \
0.0050712, 0.00523703, 0.00540828, 0.00558513, 0.00576776, \
0.00595637, 0.00615114, 0.00635229, 0.00656, 0.00677452, 0.00699604, \
0.00722481, 0.00746107, 0.00770504, 0.007957, 0.00821719, 0.00848589, \
0.00876338, 0.00904994, 0.00934588, 0.00965149, 0.00996709, 0.010293, \
0.0106296, 0.0109772, 0.0113361, 0.0117068, 0.0120896, 0.012485, \
0.0128932, 0.0133148, 0.0137502, 0.0141999, 0.0146642, 0.0151437, \
0.0156389, 0.0161503, 0.0166784, 0.0172238, 0.017787, 0.0183687, \
0.0189693, 0.0195896, 0.0202302, 0.0208917, 0.0215749, 0.0222804, \
0.023009, 0.0237614, 0.0245384, 0.0253408, 0.0261694, 0.0270251, \
0.0279089, 0.0288215, 0.0297639, 0.0307372, 0.0317423, 0.0327803, \
0.0338522, 0.0349592, 0.0361024, 0.0372829, 0.0385021, 0.0397611, \
0.0410613, 0.042404, 0.0437906, 0.0452225, 0.0467013, 0.0482284, \
0.0498055, 0.0514341, 0.053116, 0.0548529, 0.0566466, 0.058499, \
0.0604119, 0.0623874, 0.0644274, 0.0665342, 0.0687099, 0.0709567, \
0.073277, 0.0756731, 0.0781476, 0.0807031, 0.083342, 0.0860673, \
0.0888817, 0.0917882, 0.0947896, 0.0978893, 0.10109, 0.104396, \
0.10781, 0.111335, 0.114976, 0.118735, 0.122618, 0.126628, 0.130768, \
0.135044, 0.13946, 0.144021, 0.14873, 0.153594, 0.158616, 0.163803, \
0.169159, 0.174691, 0.180403, 0.186302, 0.192395, 0.198686, 0.205183, \
0.211892, 0.218821, 0.225977, 0.233366, 0.240997, 0.248878, 0.257016, \
0.265421, 0.2741, 0.283063, 0.292319, 0.3};

static const double gPhi_sub2[] = {0.000527978, 0.000535297, 0.000542751, 0.000550343, 0.000558078, \
0.000565958, 0.000573987, 0.000582167, 0.000590504, 0.000599, \
0.00060766, 0.000616487, 0.000625485, 0.000634659, 0.000644013, \
0.000653551, 0.000663278, 0.000673198, 0.000683317, 0.000693639, \
0.000704169, 0.000714913, 0.000725876, 0.000737064, 0.000748482, \
0.000760137, 0.000772035, 0.000784181, 0.000796583, 0.000809247, \
0.00082218, 0.00083539, 0.000848883, 0.000862667, 0.000876751, \
0.000891143, 0.00090585, 0.000920881, 0.000936246, 0.000951954, \
0.000968015, 0.000984437, 0.00100123, 0.00101841, 0.00103598, \
0.00105396, 0.00107236, 0.00109118, 0.00111045, 0.00113017, \
0.00115036, 0.00117103, 0.0011922, 0.00121388, 0.00123608, \
0.00125883, 0.00128214, 0.00130603, 0.00133052, 0.00135561, \
0.00138134, 0.00140773, 0.00143478, 0.00146254, 0.00149101, \
0.00152022, 0.0015502, 0.00158097, 0.00161255, 0.00164498, \
0.00167829, 0.00171249, 0.00174763, 0.00178373, 0.00182083, \
0.00185896, 0.00189816, 0.00193847, 0.00197992, 0.00202256, \
0.00206643, 0.00211157, 0.00215802, 0.00220585, 0.0022551, \
0.00230581, 0.00235806, 0.00241188, 0.00246736, 0.00252454, \
0.00258349, 0.00264428, 0.002707, 0.0027717, 0.00283847, 0.00290739, \
0.00297856, 0.00305206, 0.00312798, 0.00320644, 0.00328752, \
0.00337136, 0.00345806, 0.00354774, 0.00364054, 0.00373658, \
0.00383602, 0.00393901, 0.00404569, 0.00415625, 0.00427086, \
0.00438969, 0.00451296, 0.00464086, 0.00477362, 0.00491147, \
0.00505465, 0.00520342, 0.00535805, 0.00551885, 0.00568611, \
0.00586017, 0.00604136, 0.00623006, 0.00642666, 0.00663158, \
0.00684526, 0.00706816, 0.00730079, 0.00754369, 0.00779742, \
0.0080626, 0.00833986, 0.00862992, 0.0089335, 0.0092514, 0.00958447, \
0.00993363, 0.0102998, 0.0106842, 0.0110877, 0.0115118, 0.0119576, \
0.0124265, 0.0129201, 0.0134401, 0.0139881, 0.0145661, 0.0151762, \
0.0158206, 0.0165017, 0.0172222, 0.0179849, 0.0187931, 0.01965, \
0.0205593, 0.0215253, 0.0225522, 0.0236449, 0.0248088, 0.0260497, \
0.0273741, 0.0287889, 0.0303021, 0.0319223, 0.033659, 0.0355228, \
0.0375255, 0.0396801, 0.0420012, 0.044505, 0.0472099, 0.0501361, \
0.0533066, 0.0567473, 0.0604871, 0.0645592, 0.0690008, 0.0738545, \
0.0791686, 0.0849986, 0.091408, 0.0984697, 0.106268, 0.1149, 0.12448, \
0.13514};

/******* B-spline knots over the parameter space *******/
static const double etavec_sub1[] = {0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, \
    0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, \
    0.032, 0.034, 0.036, 0.038, 0.04, 0.042, 0.044, 0.048, 0.05, 0.056, \
    0.06, 0.065, 0.073, 0.082, 0.089, 0.1, 0.11, 0.12, 0.133, 0.15, \
    0.162, 0.18, 0.195, 0.21, 0.22, 0.2225, 0.225, 0.2275, 0.23, 0.2325, \
    0.235, 0.2375, 0.24, 0.25};
static const double chi1vec_sub1[] = {-1., -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0., 0.15, 0.3, 0.45, 0.6, \
    0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99};
static const double chi2vec_sub1[] = {-1., -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0., 0.15, 0.3, 0.45, 0.6, \
    0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99};

static const int ncx_sub1 = 54+2;       // points in eta  + 2
static const int ncy_sub1 = 24+2;       // points in chi1 + 2
static const int ncz_sub1 = 24+2;       // points in chi2 + 2


static const double etavec_sub2[] = {0.242, 0.244, 0.245, 0.246, 0.247, 0.248, 0.249, 0.2495, 0.2498, \
0.2499, 0.25};
static const double chi1vec_sub2[] = {-1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, \
-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., \
0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
0.65, 0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.88, 0.89, 0.9, 0.91, \
0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.97, 0.975, 0.98, 0.983, 0.985, \
0.986, 0.987, 0.988, 0.989, 0.99};
static const double chi2vec_sub2[] = {-1., -0.95, -0.9, -0.85, -0.8, -0.75, -0.7, -0.65, -0.6, -0.55, \
-0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0., \
0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, \
0.65, 0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.88, 0.89, 0.9, 0.91, \
0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.97, 0.975, 0.98, 0.983, 0.985, \
0.986, 0.987, 0.988, 0.989, 0.99};

static const int ncx_sub2 = 11+2;       // points in eta  + 2
static const int ncy_sub2 = 60+2;       // points in chi1 + 2
static const int ncz_sub2 = 60+2;       // points in chi2 + 2


static const double etavec_sub3[] = {0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.016, 0.017, 0.018, 0.019, \
0.02, 0.022, 0.024, 0.026, 0.028, 0.03, 0.036, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};
static const double chi1vec_sub3[] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.825, 0.85, 0.875, 0.88, \
0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.97, 0.975, \
0.98, 0.983, 0.985, 0.986, 0.987, 0.988, 0.989, 0.99};
static const double chi2vec_sub3[] = {-1., -0.9, -0.75, -0.6, -0.45, -0.3, -0.15, 0., 0.15, 0.3, 0.45, \
0.6, 0.7, 0.75, 0.8, 0.85, 0.875, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.99};

static const int ncx_sub3 = 24+2;       // points in eta  + 2
static const int ncy_sub3 = 30+2;       // points in chi1 + 2
static const int ncz_sub3 = 24+2;       // points in chi2 + 2

#ifdef LAL_PTHREAD_LOCK
static pthread_once_t SEOBNRv2ROMDoubleSpin_is_initialized = PTHREAD_ONCE_INIT;
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
  gsl_vector* cvec_amp_pre;  // AMplitude prefactor coefficient
  int nk_amp;                // Number frequency points for amplitude
  int nk_phi;                // Number of frequency points for phase
  const double *gA;          // Sparse frequency points for amplitude
  const double *gPhi;        // Sparse frequency points for phase
  const double *etavec;      // B-spline knots in eta
  const double *chi1vec;     // B-spline knots in chi1
  const double *chi2vec;     // B-spline knots in chi2
  int ncx, ncy, ncz;         // Number of points in eta, chi1, chi2
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

static SEOBNRROMdataDS __lalsim_SEOBNRv2ROMDS_data;

typedef int (*load_dataPtr)(const char*, gsl_vector *, gsl_vector *, gsl_matrix *, gsl_matrix *, gsl_vector *);

typedef struct tagSplineData
{
  gsl_bspline_workspace *bwx;
  gsl_bspline_workspace *bwy;
  gsl_bspline_workspace *bwz;
} SplineData;

/**************** Internal functions **********************/

static void SEOBNRv2ROMDoubleSpin_Init_LALDATA(void);
static int SEOBNRv2ROMDoubleSpin_Init(const char dir[]);
static bool SEOBNRv2ROMDoubleSpin_IsSetup(void);

static int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[]);
static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata);

static int TP_Spline_interpolation_3d(
  REAL8 eta,                // Input: eta-value for which projection coefficients should be evaluated
  REAL8 chi1,               // Input: chi1-value for which projection coefficients should be evaluated
  REAL8 chi2,               // Input: chi2-value for which projection coefficients should be evaluated
  gsl_vector *cvec_amp,     // Input: data for spline coefficients for amplitude
  gsl_vector *cvec_phi,     // Input: data for spline coefficients for phase
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  int nk_amp,               // number of SVD-modes == number of basis functions for amplitude
  int nk_phi,               // number of SVD-modes == number of basis functions for phase
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *etavec,     // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
);

static int SEOBNRROMdataDS_Init_submodel(
  SEOBNRROMdataDS_submodel **submodel,
  const int nk_amp,
  const int nk_phi,
  const double *gA,
  const double *gPhi,
  const double *etavec,
  const double *chi1vec,
  const double *chi2vec,
  const int ncx,
  const int ncy,
  const int ncz,
  const char dir[],
  load_dataPtr load_data
);

static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel);

/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv2ROMDoubleSpinCore(
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
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
);

static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi);
static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff);

static size_t NextPow2(const size_t n);
static void SplineData_Destroy(SplineData *splinedata);
static void SplineData_Init(
  SplineData **splinedata,
  int ncx,                // Number of points in eta  + 2
  int ncy,                // Number of points in chi1 + 2
  int ncz,                // Number of points in chi2 + 2
  const double *etavec,   // B-spline knots in eta
  const double *chi1vec,  // B-spline knots in chi1
  const double *chi2vec   // B-spline knots in chi2
);

static int load_data_sub1(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);
static int load_data_sub2(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);
static int load_data_sub3(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre);

static int SEOBNRv2ROMDoubleSpinTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi1,                                   // Aligned spin of companion 1
  REAL8 chi2                                    // Aligned spin of companion 2
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

/********************* Definitions begin here ********************/

/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in dir
 */
static int SEOBNRv2ROMDoubleSpin_Init(const char dir[]) {
  if(__lalsim_SEOBNRv2ROMDS_data.setup) {
    XLALPrintError("Error: DSEOBNRROMdata was already set up!");
    XLAL_ERROR(XLAL_EFAILED);
  }

  SEOBNRROMdataDS_Init(&__lalsim_SEOBNRv2ROMDS_data, dir);

  if(__lalsim_SEOBNRv2ROMDS_data.setup) {
    return(XLAL_SUCCESS);
  }
  else {
    return(XLAL_EFAILED);
  }
}

/** Helper function to check if the SEOBNRv2ROMDoubleSpin model has been initialised */
static bool SEOBNRv2ROMDoubleSpin_IsSetup(void) {
  if(__lalsim_SEOBNRv2ROMDS_data.setup)
    return true;
  else
    return false;
}

// Read binary ROM data for basis functions and coefficients for submodel 1
static int load_data_sub1(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "Core" submodel 1
  // B-spline points: 54x24x24
  // Frequency points: {133, 139}
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub1_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub1_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub1_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub1_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub1_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
}

// Read binary ROM data for basis functions and coefficients for submodel 2
static int load_data_sub2(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "Near q=1" submodel 2
  // B-spline points: 11x60x60
  // Frequency points: {200, 187}
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub2_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub2_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub2_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub2_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub2_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
}

// Read binary ROM data for basis functions and coefficients for submodel 3
static int load_data_sub3(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "High q, high chi1" submodel 3
  // B-spline points: 24x30x24
  // frequency points: {133, 139}
  int ret = XLAL_SUCCESS;
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub3_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub3_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub3_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_sub3_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_sub3_AmpPrefac_ci.dat", cvec_amp_pre);
  return(ret);
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
  gsl_vector *cvec_amp_pre, // Input: data for spline coefficients for amplitude prefactor
  int nk_amp,               // number of SVD-modes == number of basis functions for amplitude
  int nk_phi,               // number of SVD-modes == number of basis functions for phase
  int ncx,                  // Number of points in eta  + 2
  int ncy,                  // Number of points in chi1 + 2
  int ncz,                  // Number of points in chi2 + 2
  const double *etavec,     // B-spline knots in eta
  const double *chi1vec,    // B-spline knots in chi1
  const double *chi2vec,    // B-spline knots in chi2
  gsl_vector *c_amp,        // Output: interpolated projection coefficients for amplitude
  gsl_vector *c_phi,        // Output: interpolated projection coefficients for phase
  REAL8 *amp_pre            // Output: interpolated amplitude prefactor
) {
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

  // Evaluate the TP spline for the amplitude prefactor
  *amp_pre = Interpolate_Coefficent_Tensor(cvec_amp_pre, eta, chi1, chi2, ncy, ncz, bwx, bwy, bwz);

  SplineData_Destroy(splinedata);

  return(0);
}

/* Set up a new ROM submodel, using data contained in dir */
static int SEOBNRROMdataDS_Init_submodel(
  SEOBNRROMdataDS_submodel **submodel,
  const int nk_amp,
  const int nk_phi,
  const double *gA,
  const double *gPhi,
  const double *etavec,
  const double *chi1vec,
  const double *chi2vec,
  const int ncx,
  const int ncy,
  const int ncz,
  const char dir[],
  load_dataPtr load_data
) {
  int ret = XLAL_FAILURE;

  if(!submodel) exit(1);
  /* Create storage for submodel structures */
  if (!*submodel)
    *submodel = XLALCalloc(1,sizeof(SEOBNRROMdataDS_submodel));
  else
    SEOBNRROMdataDS_Cleanup_submodel(*submodel);

  int N = ncx*ncy*ncz; // Total number of points over parameter space = size of the data matrix for one SVD-mode

  // Initalize actual ROM data
  (*submodel)->cvec_amp = gsl_vector_alloc(N*nk_amp);
  (*submodel)->cvec_phi = gsl_vector_alloc(N*nk_phi);
  (*submodel)->Bamp = gsl_matrix_alloc(nk_amp, nk_amp);
  (*submodel)->Bphi = gsl_matrix_alloc(nk_phi, nk_phi);
  (*submodel)->cvec_amp_pre = gsl_vector_alloc(N);

  // Load ROM data for this submodel
  ret=load_data(dir, (*submodel)->cvec_amp, (*submodel)->cvec_phi, (*submodel)->Bamp, (*submodel)->Bphi, (*submodel)->cvec_amp_pre);

  // Initialize other members
  (*submodel)->nk_amp = nk_amp;
  (*submodel)->nk_phi = nk_phi;
  (*submodel)->gA = gA;
  (*submodel)->gPhi = gPhi;
  (*submodel)->etavec = etavec;
  (*submodel)->chi1vec = chi1vec;
  (*submodel)->chi2vec = chi2vec;
  (*submodel)->ncx = ncx;
  (*submodel)->ncy = ncy;
  (*submodel)->ncz = ncz;

  return ret;
}

/* Deallocate contents of the given SEOBNRROMdataDS_submodel structure */
static void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel) {
  if(submodel->cvec_amp) gsl_vector_free(submodel->cvec_amp);
  if(submodel->cvec_phi) gsl_vector_free(submodel->cvec_phi);
  if(submodel->Bamp) gsl_matrix_free(submodel->Bamp);
  if(submodel->Bphi) gsl_matrix_free(submodel->Bphi);
  if(submodel->cvec_amp_pre) gsl_vector_free(submodel->cvec_amp_pre);
}

/* Set up a new ROM model, using data contained in dir */
int SEOBNRROMdataDS_Init(SEOBNRROMdataDS *romdata, const char dir[]) {
  int ret = XLAL_FAILURE;

  /* Create storage for structures */
  if(romdata->setup) {
    XLALPrintError("WARNING: You tried to setup the SEOBNRv2ROMDoubleSpin model that was already initialised. Ignoring\n");
    return (XLAL_FAILURE);
  }

  gsl_set_error_handler(&err_handler);

  load_dataPtr load_data = &load_data_sub1;
  ret = SEOBNRROMdataDS_Init_submodel(&(romdata)->sub1, nk_amp_sub1, nk_phi_sub1,
          gA_sub1, gPhi_sub1, etavec_sub1, chi1vec_sub1, chi2vec_sub1, ncx_sub1, ncy_sub1, ncz_sub1, dir, load_data);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 1 loaded sucessfully.\n", __func__);

  load_data = &load_data_sub2;
  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->sub2, nk_amp_sub2, nk_phi_sub2,
          gA_sub2, gPhi_sub2, etavec_sub2, chi1vec_sub2, chi2vec_sub2, ncx_sub2, ncy_sub2, ncz_sub2, dir, load_data);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 2 loaded sucessfully.\n", __func__);

  load_data = &load_data_sub3;
  ret |= SEOBNRROMdataDS_Init_submodel(&(romdata)->sub3, nk_amp_sub3, nk_phi_sub3,
          gA_sub3, gPhi_sub3, etavec_sub3, chi1vec_sub3, chi2vec_sub3, ncx_sub3, ncy_sub3, ncz_sub3, dir, load_data);
  if (ret==XLAL_SUCCESS) XLALPrintInfo("%s : submodel 3 loaded sucessfully.\n", __func__);

  if(XLAL_SUCCESS==ret)
    romdata->setup=1;
  else
    SEOBNRROMdataDS_Cleanup(romdata);

  return (ret);
}

/* Deallocate contents of the given SEOBNRROMdataDS structure */
static void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata) {
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->sub1);
  XLALFree((romdata)->sub1);
  (romdata)->sub1 = NULL;
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


/**
 * Core function for computing the ROM waveform.
 * Interpolate projection coefficient data and evaluate coefficients at desired (q, chi1, chi2).
 * Construct 1D splines for amplitude and phase.
 * Compute strain waveform from amplitude and phase.
*/
static int SEOBNRv2ROMDoubleSpinCore(
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
  double deltaF
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  )
{

  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv2ROMDS_data;
  if(*hptilde || *hctilde)
  {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (eta > 0.25)  nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01)  nudge(&eta, 0.01, 1e-6);
  if (chi1 < -1.0) nudge(&chi1, -1.0, 1e-6);
  if (chi1 > 0.99) nudge(&chi1, 0.99, 1e-6);
  if (chi2 < -1.0) nudge(&chi2, -1.0, 1e-6);
  if (chi1 > 0.99) nudge(&chi2, 0.99, 1e-6);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 0.99 || chi2 > 0.99 ) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1 or larger than 0.99!\nSEOBNRv2ROMDoubleSpin is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta<0.01 || eta > 0.25) {
    XLALPrintError( "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv2ROMDoubleSpin is only available for eta in the range 0.01 <= eta <= 0.25.\n", __func__,eta);
    XLAL_ERROR( XLAL_EDOM );
  }

  /* Select ROM submodel */
  SEOBNRROMdataDS_submodel *submodel;
  if (eta >= 0.242)
    submodel = romdata->sub2;
  else if (eta < 0.1 && chi1 > 0.5)
    submodel = romdata->sub3;
  else
    submodel = romdata->sub1;

  /* Find frequency bounds */
  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];

  if(fRef==0.0)
    fRef=fLow;

  /* Convert to geometric units for frequency */
  double Mf_ROM_min = fmax(submodel->gA[0], submodel->gPhi[0]);                                   // lowest allowed geometric frequency for ROM
  double Mf_ROM_max = fmin(submodel->gA[submodel->nk_amp-1], submodel->gPhi[submodel->nk_phi-1]); // highest allowed geometric frequency for ROM
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
    XLALPrintWarning("Reference frequency Mf_ref=%g is smaller than lowest frequency in ROM Mf=%g. Starting at lowest frequency in ROM.\n", fRef_geom, Mf_ROM_min);
    fRef_geom = Mf_ROM_min;
  }

  /* Internal storage for waveform coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff, submodel->nk_amp, submodel->nk_phi);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (q,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    eta,                       // Input: eta-value for which projection coefficients should be evaluated
    chi1,                      // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                      // Input: chi2-value for which projection coefficients should be evaluated
    submodel->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel->cvec_phi,        // Input: data for spline coefficients for phase
    submodel->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel->nk_phi,          // number of SVD-modes == number of basis functions for phase
    submodel->ncx,             // Number of points in eta  + 2
    submodel->ncy,             // Number of points in chi1 + 2
    submodel->ncz,             // Number of points in chi2 + 2
    submodel->etavec,          // B-spline knots in eta
    submodel->chi1vec,         // B-spline knots in chi1
    submodel->chi2vec,         // B-spline knots in chi2
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode, "Parameter-space interpolation failed.");
  }

  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // amp_pts = B_A^T . c_A
  // phi_pts = B_phi^T . c_phi
  gsl_vector* amp_f = gsl_vector_alloc(submodel->nk_amp);
  gsl_vector* phi_f = gsl_vector_alloc(submodel->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Bamp, romdata_coeff->c_amp, 0.0, amp_f);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Bphi, romdata_coeff->c_phi, 0.0, phi_f);

  // Setup 1d splines in frequency
  gsl_interp_accel *acc_amp = gsl_interp_accel_alloc();
  gsl_spline *spline_amp = gsl_spline_alloc(gsl_interp_cspline, submodel->nk_amp);
  gsl_spline_init(spline_amp, submodel->gA, gsl_vector_const_ptr(amp_f,0), submodel->nk_amp);

  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_spline *spline_phi = gsl_spline_alloc(gsl_interp_cspline, submodel->nk_phi);
  gsl_spline_init(spline_phi, submodel->gPhi, gsl_vector_const_ptr(phi_f,0), submodel->nk_phi);


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

  REAL8 s = 0.5; // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv2
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

  // Get SEOBNRv2 ringdown frequency for 22 mode
  // XLALSimInspiralGetFinalFreq wants masses in SI units, so unfortunately we need to convert back
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  double Mtot_SI = Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * 1.0/(1.0+q);
  double m2_SI = Mtot_SI * q/(1.0+q);
  double Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi1, 0,0,chi2, SEOBNRv2) * Mtot_sec;

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
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRROM_c
 *
 * @{
 *
 * @name SEOBNRv2 Reduced Order Model (Double Spin)
 *
 * @author Michael Puerrer, John Veitch
 *
 * @brief C code for SEOBNRv2 reduced order model (double spin version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for details.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv2 model.
 *
 * The binary data files are available at https://dcc.ligo.org/T1400701-v1.
 * Put the untared data into a location in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv2 model.
 *
 * @note Parameter ranges:
 *   * 0.01 <= eta <= 0.25
 *   * -1 <= chi_i <= 0.99
 *   * Mtot >= 12Msun
 *
 *  Aligned component spins chi1, chi2.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */

/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv2_ROM_DoubleSpin model.
 *
 * XLALSimIMRSEOBNRv2ROMDoubleSpin() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz), need to be strictly monotonically increasing */
  REAL8 phiRef,                                 /**< Orbital phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1SI,                                   /**< Mass of companion 1 (kg) */
  REAL8 m2SI,                                   /**< Mass of companion 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin 1 */
  REAL8 chi2                                    /**< Dimensionless aligned component spin 2 */
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

  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1+mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  double Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv2ROMDoubleSpin_is_initialized, SEOBNRv2ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv2ROMDoubleSpin_Init_LALDATA();
#endif

  if(!SEOBNRv2ROMDoubleSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv2ROMDoubleSpin data - check your $LAL_DATA_PATH\n");

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv2ROMDoubleSpinCore(hptilde,hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi1, chi2, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv2_ROM_DoubleSpin model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency fHigh to the next power of 2 in
 * the size of the frequency series.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpin(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Orbital phase at reference frequency*/
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
  (void) pthread_once(&SEOBNRv2ROMDoubleSpin_is_initialized, SEOBNRv2ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv2ROMDoubleSpin_Init_LALDATA();
#endif

  if(!SEOBNRv2ROMDoubleSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv2ROMDoubleSpin data - check your $LAL_DATA_PATH\n");

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv2ROMDoubleSpinCore(hptilde,hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi1, chi2, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/**
 * Compute the 'time' elapsed in the ROM waveform from a given starting frequency until the ringdown.
 * UNREVIEWED!
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * The allowed frequency range for the starting frequency in geometric frequency is [0.00053, 0.135].
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencyOfTime() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinTimeOfFrequency(
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
  int ret = SEOBNRv2ROMDoubleSpinTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final, &Mtot_sec, m1SI, m2SI, chi1, chi2);
  if(ret != 0)
    XLAL_ERROR(ret);

  // ROM frequency bounds in Mf
  double Mf_ROM_min = 0.00053;
  double Mf_ROM_max = 0.135;

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  double Mf = frequency * Mtot_sec;
  if (Mf < Mf_ROM_min || Mf > Mf_ROM_max)
  {
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
 * UNREVIEWED!
 *
 * The notion of elapsed 'time' (in seconds) is defined here as the difference of the
 * frequency derivative of the frequency domain phase between the ringdown frequency
 * and the starting frequency ('frequency' argument). This notion of time is similar to the
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * If the frequency that corresponds to the specified elapsed time is lower than the
 * geometric frequency Mf=0.00053 (ROM starting frequency) or above half of the SEOBNRv2
 * ringdown frequency an error is thrown.
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMDoubleSpinTimeOfFrequency() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinFrequencyOfTime(
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
  int ret = SEOBNRv2ROMDoubleSpinTimeFrequencySetup(&spline_phi, &acc_phi, &Mf_final, &Mtot_sec, m1SI, m2SI, chi1, chi2);
  if(ret != 0)
    XLAL_ERROR(ret);

  // ROM frequency bounds in Mf
  double Mf_ROM_min = 0.00053;

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

  // Assume for now that we only care about f(t) *before* merger so that f(t) - f_ringdown >= 0.
  // Assume that we only need to cover the frequency range [f_min, f_ringdown/2].
  int N = 20;
  double log_f_pts[N];
  double log_t_pts[N];
  double log_f_min   = log(Mf_ROM_min);
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
  if (t < t_rng_2 || t > t_min)
  {
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

/** @} */
/** @} */

// Auxiliary function to perform setup of phase spline for t(f) and f(t) functions
static int SEOBNRv2ROMDoubleSpinTimeFrequencySetup(
  gsl_spline **spline_phi,                      // phase spline
  gsl_interp_accel **acc_phi,                   // phase spline accelerator
  REAL8 *Mf_final,                              // ringdown frequency in Mf
  REAL8 *Mtot_sec,                              // total mass in seconds
  REAL8 m1SI,                                   // Mass of companion 1 (kg)
  REAL8 m2SI,                                   // Mass of companion 2 (kg)
  REAL8 chi1,                                   // Aligned spin of companion 1
  REAL8 chi2                                    // Aligned spin of companion 2
)
{
  /* Get masses in terms of solar mass */
  double mass1 = m1SI / LAL_MSUN_SI;
  double mass2 = m2SI / LAL_MSUN_SI;
  double Mtot = mass1 + mass2;
  double eta = mass1 * mass2 / (Mtot*Mtot);    /* Symmetric mass-ratio */
  *Mtot_sec = Mtot * LAL_MTSUN_SI; /* Total mass in seconds */

  // 'Nudge' parameter values to allowed boundary values if close by
  nudge(&eta, 0.25, 1e-6);
  nudge(&eta, 0.01, 1e-6);
  nudge(&chi1, -1.0, 1e-6);
  nudge(&chi1, 0.99, 1e-6);
  nudge(&chi2, -1.0, 1e-6);
  nudge(&chi2, 0.99, 1e-6);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 0.99 || chi2 > 0.99 ) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1 or larger than 0.99!\nSEOBNRv2ROMDoubleSpin is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta<0.01 || eta > 0.25) {
    XLALPrintError( "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv2ROMDoubleSpin is only available for spins in the range 0.01 <= eta <= 0.25.\n", __func__,eta);
    XLAL_ERROR( XLAL_EDOM );
  }

  // Load ROM data if not loaded already
#ifdef LAL_PTHREAD_LOCK
  (void) pthread_once(&SEOBNRv2ROMDoubleSpin_is_initialized, SEOBNRv2ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv2ROMDoubleSpin_Init_LALDATA();
#endif

  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv2ROMDS_data;

  /* Select ROM submodel */
  SEOBNRROMdataDS_submodel *submodel;
  if (eta >= 0.242)
    submodel = romdata->sub2;
  else if (eta < 0.1 && chi1 > 0.5)
    submodel = romdata->sub3;
  else
    submodel = romdata->sub1;

  /* Internal storage for w.f. coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff, submodel->nk_amp, submodel->nk_phi);
  REAL8 amp_pre;

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  int retcode=TP_Spline_interpolation_3d(
    eta,                       // Input: eta-value for which projection coefficients should be evaluated
    chi1,                      // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                      // Input: chi2-value for which projection coefficients should be evaluated
    submodel->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel->cvec_phi,        // Input: data for spline coefficients for phase
    submodel->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel->nk_phi,          // number of SVD-modes == number of basis functions for phase
    submodel->ncx,             // Number of points in eta  + 2
    submodel->ncy,             // Number of points in chi1 + 2
    submodel->ncz,             // Number of points in chi2 + 2
    submodel->etavec,          // B-spline knots in eta
    submodel->chi1vec,         // B-spline knots in chi1
    submodel->chi2vec,         // B-spline knots in chi2
    romdata_coeff->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);
    XLAL_ERROR(retcode);
  }

  // Compute function values of phase on sparse frequency points by evaluating matrix vector products
  // phi_pts = B_phi^T . c_phi
  gsl_vector* phi_f = gsl_vector_alloc(submodel->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel->Bphi, romdata_coeff->c_phi, 0.0, phi_f);

  // Setup 1d phase spline in frequency
  *acc_phi = gsl_interp_accel_alloc();
  *spline_phi = gsl_spline_alloc(gsl_interp_cspline, submodel->nk_phi);
  gsl_spline_init(*spline_phi, submodel->gPhi, gsl_vector_const_ptr(phi_f,0), submodel->nk_phi);

  // Get SEOBNRv2 ringdown frequency for 22 mode
  double q = (1.0 + sqrt(1.0 - 4.0*eta) - 2.0*eta) / (2.0*eta);
  double Mtot_SI = *Mtot_sec / LAL_MTSUN_SI * LAL_MSUN_SI;
  double m1_SI = Mtot_SI * 1.0/(1.0+q);
  double m2_SI = Mtot_SI * q/(1.0+q);
  *Mf_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0,0,chi1, 0,0,chi2, SEOBNRv2) * (*Mtot_sec);

  gsl_vector_free(phi_f);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff);

  return(XLAL_SUCCESS);
}



/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in $LAL_DATA_PATH
 */
static void SEOBNRv2ROMDoubleSpin_Init_LALDATA(void)
{
  if (SEOBNRv2ROMDoubleSpin_IsSetup()) return;

  // If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  char datafile[] = "SEOBNRv2ROM_DS_sub1_Phase_ciall.dat";

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv2ROMDoubleSpin_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv2ROMDoubleSpin data files in $LAL_DATA_PATH\n");
}
