/*
 *  Copyright (C) 2014, 2015 Michael Puerrer, John Veitch
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
static const char ROMDataHDF5[] = "SEOBNRv2ROM_DS_HI_v1.0.hdf5";
static const INT4 ROMDataHDF5_VERSION_MAJOR = 1;
static const INT4 ROMDataHDF5_VERSION_MINOR = 0;
static const INT4 ROMDataHDF5_VERSION_MICRO = 0;
#endif

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

// Helper functions to read gsl_vector and gsl_matrix data with error checking
UNUSED static int read_vector_test(const char dir[], const char fname[], gsl_vector *v) {
  size_t size = strlen(dir) + strlen(fname) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, fname);

  FILE *f = fopen(path, "rb");
  if (!f)
    XLAL_ERROR(XLAL_EIO, "Could not find ROM data file at path `%s'", path);
  int ret = gsl_vector_fread(f, v);
  if (ret != 0)
     XLAL_ERROR(XLAL_EIO, "Error reading data from `%s'", path);
  fclose(f);

  XLAL_PRINT_INFO("Sucessfully read data file `%s'", path);
  XLALFree(path);
  return(XLAL_SUCCESS);
}

#include "LALSimIMRSEOBNRROMUtilities.c"

#include <lal/LALConfig.h>
#ifdef LAL_PTHREAD_LOCK
#include <pthread.h>
#endif


// This ROM consists of different submodels and also glues together a low mass and 2 high mass models
// Filenames of the submodels; need to load up all 3 of them

/******************************************************************
 * The double-spin SEOBNRv2 ROM consists of a number of submodels *
 * for subdomains of the whole parameter space:                   *
 *                                                                *
 * "Core Low frequency" submodel 1                                *
 * B-spline points: 67x12x12                                      *
 * Frequency points: {183, 242}                                   *
 *                                                                *
 * "High frequency, low chi1" submodel 2                          *
 * B-spline points: 63x37x51                                      *
 * Frequency points: {133, 133}                                   *
 *                                                                *
 * "High frequency, high chi1" submodel 3                         *
 * B-spline points: 97x119x51                                     *
 * frequency points: {133, 133}                                   *
 *                                                                *
 *****************************************************************/


/********* Input data for spline basis points **************/

// The frequency points are different for all submodels.

#define nk_amp_sub1 200  // number of SVD-modes == number of basis functions for amplitude
#define nk_phi_sub1 250  // number of SVD-modes == number of basis functions for phase

// Frequency points for amplitude and phase

// sub1: core low mass ROM
static const double gA_sub1[] = {0.0000985, 0.000102213, 0.000106067, 0.000110066, 0.000114215, \
0.000118521, 0.000122989, 0.000127626, 0.000132437, 0.00013743, \
0.000142611, 0.000147988, 0.000153567, 0.000159357, 0.000165364, \
0.000171598, 0.000178068, 0.000184781, 0.000191747, 0.000198976, \
0.000206477, 0.000214262, 0.000222339, 0.000230721, 0.00023942, \
0.000248446, 0.000257812, 0.000267532, 0.000277618, 0.000288084, \
0.000298945, 0.000310215, 0.00032191, 0.000334046, 0.000346639, \
0.000359708, 0.000373269, 0.000387341, 0.000401944, 0.000417097, \
0.000432822, 0.000449139, 0.000466071, 0.000483642, 0.000501876, \
0.000520796, 0.00054043, 0.000560805, 0.000581947, 0.000603886, \
0.000626653, 0.000650278, 0.000674793, 0.000700233, 0.000726632, \
0.000754026, 0.000782452, 0.000811951, 0.000842561, 0.000874326, \
0.000907288, 0.000941493, 0.000976987, 0.00101382, 0.00105204, \
0.0010917, 0.00113286, 0.00117557, 0.00121989, 0.00126588, 0.0013136, \
0.00136312, 0.00141451, 0.00146784, 0.00152318, 0.0015806, \
0.00164019, 0.00170203, 0.00176619, 0.00183278, 0.00190187, \
0.00197357, 0.00204798, 0.00212519, 0.00220531, 0.00228845, \
0.00237472, 0.00246425, 0.00255715, 0.00265355, 0.00275359, \
0.0028574, 0.00296513, 0.00307691, 0.00319291, 0.00331329, 0.0034382, \
0.00356782, 0.00370232, 0.0038419, 0.00398674, 0.00413704, \
0.00429301, 0.00445485, 0.0046228, 0.00479708, 0.00497793, 0.0051656, \
0.00536034, 0.00556243, 0.00577213, 0.00598974, 0.00621555, \
0.00644988, 0.00669304, 0.00694537, 0.00720721, 0.00747892, \
0.00776087, 0.00805346, 0.00835707, 0.00867214, 0.00899908, \
0.00933834, 0.0096904, 0.0100557, 0.0104348, 0.0108282, 0.0112364, \
0.0116601, 0.0120996, 0.0125558, 0.0130291, 0.0135203, 0.0140301, \
0.014559, 0.0151079, 0.0156774, 0.0162685, 0.0168818, 0.0175182, \
0.0181787, 0.018864, 0.0195752, 0.0203132, 0.021079, 0.0218737, \
0.0226983, 0.023554, 0.024442, 0.0253635, 0.0263197, 0.0273119, \
0.0283416, 0.0294101, 0.0305188, 0.0316694, 0.0328633, 0.0341023, \
0.0353879, 0.036722, 0.0381065, 0.0395431, 0.0410339, 0.0425808, \
0.0441861, 0.0458519, 0.0475806, 0.0493744, 0.0512358, 0.0531674, \
0.0551718, 0.0572517, 0.0594101, 0.0616499, 0.0639741, 0.0663859, \
0.0688887, 0.0714858, 0.0741808, 0.0769774, 0.0798794, 0.0828909, \
0.0860159, 0.0892587, 0.0926237, 0.0961157, 0.0997392, 0.103499, \
0.107401, 0.11145, 0.115652, 0.120012, 0.124537, 0.129232, 0.134104, \
0.139159, 0.144406, 0.14985, 0.15};

static const double gPhi_sub1[] = {0.0000985, 0.0000996191, 0.000100755, 0.000101908, 0.000103079, \
0.000104268, 0.000105476, 0.000106702, 0.000107947, 0.000109211, \
0.000110495, 0.000111799, 0.000113124, 0.00011447, 0.000115838, \
0.000117227, 0.000118638, 0.000120072, 0.000121529, 0.00012301, \
0.000124515, 0.000126045, 0.000127599, 0.00012918, 0.000130786, \
0.000132419, 0.000134079, 0.000135768, 0.000137484, 0.00013923, \
0.000141005, 0.00014281, 0.000144647, 0.000146515, 0.000148415, \
0.000150348, 0.000152314, 0.000154315, 0.000156352, 0.000158424, \
0.000160532, 0.000162679, 0.000164863, 0.000167087, 0.000169351, \
0.000171656, 0.000174003, 0.000176393, 0.000178826, 0.000181305, \
0.000183829, 0.0001864, 0.00018902, 0.000191688, 0.000194407, \
0.000197177, 0.000200001, 0.000202878, 0.00020581, 0.0002088, \
0.000211847, 0.000214954, 0.000218121, 0.000221351, 0.000224645, \
0.000228005, 0.000231431, 0.000234927, 0.000238492, 0.000242131, \
0.000245843, 0.000249632, 0.000253499, 0.000257445, 0.000261474, \
0.000265587, 0.000269787, 0.000274075, 0.000278455, 0.000282928, \
0.000287497, 0.000292165, 0.000296934, 0.000301807, 0.000306787, \
0.000311877, 0.00031708, 0.000322399, 0.000327838, 0.000333398, \
0.000339086, 0.000344902, 0.000350852, 0.00035694, 0.000363169, \
0.000369543, 0.000376066, 0.000382744, 0.00038958, 0.000396579, \
0.000403747, 0.000411088, 0.000418607, 0.000426311, 0.000434203, \
0.000442292, 0.000450582, 0.000459079, 0.000467791, 0.000476725, \
0.000485886, 0.000495283, 0.000504923, 0.000514813, 0.000524964, \
0.000535381, 0.000546076, 0.000557056, 0.000568331, 0.000579912, \
0.000591808, 0.000604031, 0.000616592, 0.000629502, 0.000642774, \
0.00065642, 0.000670454, 0.000684889, 0.000699741, 0.000715023, \
0.000730752, 0.000746943, 0.000763615, 0.000780785, 0.000798472, \
0.000816695, 0.000835474, 0.000854831, 0.000874789, 0.00089537, \
0.0009166, 0.000938503, 0.000961107, 0.000984439, 0.00100853, \
0.00103341, 0.00105911, 0.00108567, 0.00111312, 0.0011415, \
0.00117085, 0.0012012, 0.00123261, 0.00126513, 0.00129879, \
0.00133365, 0.00136976, 0.00140718, 0.00144597, 0.00148619, \
0.00152792, 0.00157121, 0.00161614, 0.0016628, 0.00171126, \
0.00176161, 0.00181395, 0.00186837, 0.00192498, 0.00198389, \
0.00204521, 0.00210908, 0.00217561, 0.00224496, 0.00231727, \
0.00239271, 0.00247143, 0.00255363, 0.0026395, 0.00272923, \
0.00282306, 0.00292121, 0.00302394, 0.00313151, 0.00324421, \
0.00336236, 0.00348627, 0.00361632, 0.00375287, 0.00389633, \
0.00404716, 0.00420582, 0.00437283, 0.00454874, 0.00473414, \
0.00492969, 0.00513608, 0.00535408, 0.00558449, 0.00582823, \
0.00608624, 0.0063596, 0.00664946, 0.00695705, 0.00728377, 0.0076311, \
0.00800069, 0.00839433, 0.00881401, 0.0092619, 0.00974038, 0.0102521, \
0.0108, 0.0113872, 0.0120175, 0.0126946, 0.0134231, 0.0142079, \
0.0150544, 0.0159688, 0.0169581, 0.0180299, 0.0191929, 0.020457, \
0.0218334, 0.0233345, 0.0249749, 0.0267707, 0.0287408, 0.0309065, \
0.0332925, 0.0359272, 0.0388435, 0.0420796, 0.0456801, 0.0496972, \
0.054192, 0.0592368, 0.0649174, 0.0713356, 0.0786136, 0.086898, \
0.0963666, 0.107235, 0.115, 0.119768, 0.13, 0.134291, 0.145, 0.15};

#define nk_amp_sub2 113
#define nk_phi_sub2 113

// sub2: high frequency, low chi1
static const double g_sub2[] = {0.00739041, 0.00753436, 0.00768207, 0.00783364, 0.00798922, \
0.00814894, 0.00831292, 0.00848132, 0.00865428, 0.00883196, \
0.00901452, 0.00920213, 0.00939497, 0.00959321, 0.00979705, \
0.0100067, 0.0102223, 0.0104442, 0.0106725, 0.0109074, 0.0111493, \
0.0113984, 0.0116549, 0.0119192, 0.0121915, 0.012472, 0.0127613, \
0.0130595, 0.013367, 0.0136843, 0.0140116, 0.0143494, 0.0146981, \
0.0150581, 0.0154299, 0.0158141, 0.016211, 0.0166213, 0.0170455, \
0.0174842, 0.017938, 0.0184075, 0.0188935, 0.0193968, 0.0199179, \
0.0204578, 0.0210174, 0.0215974, 0.0221988, 0.0228227, 0.0234701, \
0.0241421, 0.0248398, 0.0255646, 0.0263177, 0.0271005, 0.0279145, \
0.0287613, 0.0296425, 0.0305599, 0.0315153, 0.0325108, 0.0335484, \
0.0346304, 0.0357592, 0.0369373, 0.0381675, 0.0394525, 0.0407956, \
0.0422, 0.0436692, 0.045207, 0.0468174, 0.0485048, 0.0502737, \
0.0521292, 0.0540765, 0.0561214, 0.0582701, 0.0605292, 0.0629058, \
0.0654076, 0.0680429, 0.0708208, 0.0737509, 0.0768437, 0.0801107, \
0.0835641, 0.0872175, 0.0910854, 0.0951836, 0.0995295, 0.104142, \
0.109042, 0.114251, 0.119795, 0.1257, 0.131997, 0.138718, 0.145899, \
0.15358, 0.161804, 0.170621, 0.180084, 0.190254, 0.201196, 0.212986, \
0.225705, 0.239447, 0.254316, 0.270429, 0.287917, 0.3};

#define gA_sub2 g_sub2
#define gPhi_sub2 g_sub2

// sub3: high frequency, high chi1
#define nk_amp_sub3 nk_amp_sub2
#define nk_phi_sub3 nk_phi_sub2
#define gA_sub3 g_sub2
#define gPhi_sub3 g_sub2


/******* B-spline knots over the parameter space *******/
// These knots are different for all submodels.

// sub1: core low frequency ROM
static const double etavec_sub1[] = {0.01, 0.011, 0.012, 0.013, 0.015, 0.017, 0.018, 0.02, 0.021, 0.022, \
0.023, 0.024, 0.025, 0.027, 0.03, 0.035, 0.037, 0.04, 0.042, 0.045, \
0.048, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, \
0.095, 0.1, 0.105, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.18, 0.2, \
0.22, 0.23, 0.235, 0.24, 0.241, 0.242, 0.243, 0.244, 0.245, 0.246, \
0.247, 0.248, 0.2485, 0.2488, 0.249, 0.2491, 0.2492, 0.2493, 0.2494, \
0.2495, 0.2496, 0.2497, 0.2498, 0.2499, 0.24995, 0.25};
static const double chi1vec_sub1[] = {-0.9999, -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 0.9, \
0.98999};
static const double chi2vec_sub1[] = {-0.9999, -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 0.9, \
0.98999};

static const int ncx_sub1 = 67+2;       // points in eta  + 2
static const int ncy_sub1 = 12+2;       // points in chi1 + 2
static const int ncz_sub1 = 12+2;       // points in chi2 + 2

// sub2: high frequency, low chi1
static const double etavec_sub2[] = {0.01, 0.0105, 0.011, 0.0115, 0.0125, 0.015, 0.0175, 0.02, 0.025, \
0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, \
0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, \
0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, \
0.185, 0.19, 0.195, 0.2, 0.205, 0.21, 0.215, 0.22, 0.225, 0.23, \
0.235, 0.24, 0.2425, 0.2435, 0.245, 0.24625, 0.2475, 0.2485, 0.249, \
0.2495, 0.2498, 0.2499, 0.25};
static const double chi1vec_sub2[] = {-0.9999, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, \
-0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, \
-0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0., 0.04, 0.08, 0.12, 0.16, \
0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44};
static const double chi2vec_sub2[] = {-0.9999, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, \
-0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, \
-0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0., 0.04, 0.08, 0.12, 0.16, \
0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, \
0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 0.98999};

static const int ncx_sub2 = 63+2;       // points in eta  + 2
static const int ncy_sub2 = 37+2;       // points in chi1 + 2
static const int ncz_sub2 = 51+2;       // points in chi2 + 2

// sub3: high frequency, high chi1
static const double etavec_sub3[] = {0.01, 0.0105, 0.011, 0.0115, 0.0125, 0.015, 0.0175, 0.02, 0.0225, \
0.025, 0.0275, 0.03, 0.0325, 0.035, 0.0375, 0.04, 0.0425, 0.045, \
0.0475, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.0625, 0.065, 0.0675, \
0.07, 0.0725, 0.075, 0.0775, 0.08, 0.0825, 0.085, 0.0875, 0.09, \
0.0925, 0.095, 0.0975, 0.1, 0.1025, 0.105, 0.1075, 0.11, 0.1125, \
0.115, 0.1175, 0.12, 0.1225, 0.125, 0.1275, 0.13, 0.1325, 0.135, \
0.1375, 0.14, 0.1425, 0.145, 0.1475, 0.15, 0.1525, 0.155, 0.1575, \
0.16, 0.1625, 0.165, 0.1675, 0.17, 0.1725, 0.175, 0.1775, 0.18, \
0.1825, 0.185, 0.1875, 0.19, 0.1925, 0.195, 0.1975, 0.2, 0.2025, \
0.205, 0.2075, 0.21, 0.2125, 0.215, 0.2175, 0.22, 0.2225, 0.225, \
0.2275, 0.23, 0.2325, 0.235, 0.2375, 0.24, 0.2425, 0.245, 0.24625, \
0.2475, 0.2485, 0.249, 0.2495, 0.2498, 0.2499, 0.25};
static const double chi1vec_sub3[] = {0.4, 0.405, 0.41, 0.415, 0.42, 0.425, 0.43, 0.435, 0.44, 0.445, \
0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49, 0.495, 0.5, \
0.505, 0.51, 0.515, 0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, \
0.555, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.6, \
0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.635, 0.64, 0.645, 0.65, \
0.655, 0.66, 0.665, 0.67, 0.675, 0.68, 0.685, 0.69, 0.695, 0.7, \
0.705, 0.71, 0.715, 0.72, 0.725, 0.73, 0.735, 0.74, 0.745, 0.75, \
0.755, 0.76, 0.765, 0.77, 0.775, 0.78, 0.785, 0.79, 0.795, 0.8, \
0.805, 0.81, 0.815, 0.82, 0.825, 0.83, 0.835, 0.84, 0.845, 0.85, \
0.855, 0.86, 0.865, 0.87, 0.875, 0.88, 0.885, 0.89, 0.895, 0.9, \
0.905, 0.91, 0.915, 0.92, 0.925, 0.93, 0.935, 0.94, 0.945, 0.95, \
0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.98999};
static const double chi2vec_sub3[] = {-0.9999, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, \
-0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, \
-0.24, -0.2, -0.16, -0.12, -0.08, -0.04, 0., 0.04, 0.08, 0.12, 0.16, \
0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, 0.6, 0.64, \
0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 0.98999};

static const int ncx_sub3 = 106+2;      // points in eta  + 2
static const int ncy_sub3 = 119+2;      // points in chi1 + 2
static const int ncz_sub3 = 51+2;       // points in chi2 + 2

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
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
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
  double deltaF,
  /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  int nk_max // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
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

static void GluePhasing(
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

/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in dir
 */
int SEOBNRv2ROMDoubleSpin_Init(const char dir[]) {
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
bool SEOBNRv2ROMDoubleSpin_IsSetup(void) {
  if(__lalsim_SEOBNRv2ROMDS_data.setup)
    return true;
  else
    return false;
}


// Read binary ROM data for basis functions and coefficients for submodel 1
static int load_data_sub1(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "Core Low frequency" submodel 1
  // B-spline points: 67x12x12
  // Frequency points: {183, 242}

  // Read data into preallocated gsl_vectors and gsl_matrices
  int ret = XLAL_SUCCESS;
#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  LALH5File *sub1 = XLALH5GroupOpen(file, "sub1");

  // Read ROM coefficients
  ReadHDF5RealVectorDataset(sub1, "Amp_ciall", &cvec_amp);
  ReadHDF5RealVectorDataset(sub1, "Phase_ciall", &cvec_phi);
  ReadHDF5RealVectorDataset(sub1, "AmpPrefac_ci", &cvec_amp_pre);
  // Read ROM basis functions
  ReadHDF5RealMatrixDataset(sub1, "Bamp", &Bamp);
  ReadHDF5RealMatrixDataset(sub1, "Bphase", &Bphi);

  // Check frequency points
  CheckVectorFromHDF5(sub1, "Mf_grid_Amp", gA_sub1, nk_amp_sub1);
  CheckVectorFromHDF5(sub1, "Mf_grid_Phi", gPhi_sub1, nk_phi_sub1);
  // Check parameter space nodes
  CheckVectorFromHDF5(sub1, "etavec", etavec_sub1, ncx_sub1-2);
  CheckVectorFromHDF5(sub1, "chi1vec", chi1vec_sub1, ncy_sub1-2);
  CheckVectorFromHDF5(sub1, "chi2vec", chi2vec_sub1, ncz_sub1-2);

  XLALFree(path);
  XLALH5FileClose(file);
#else
  // fall back to reading gsl binary data
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub1_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub1_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub1_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub1_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub1_AmpPrefac_ci.dat", cvec_amp_pre);
#endif
  return(ret);
}

// Read binary ROM data for basis functions and coefficients for submodel 2
static int load_data_sub2(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "High frequency, low chi1" submodel 2
  // B-spline points: 63x37x51
  // Frequency points: {113, 113}

  // Read data into preallocated gsl_vectors and gsl_matrices
  int ret = XLAL_SUCCESS;
#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  LALH5File *sub2 = XLALH5GroupOpen(file, "sub2");

  // Read ROM coefficients
  ReadHDF5RealVectorDataset(sub2, "Amp_ciall", &cvec_amp);
  ReadHDF5RealVectorDataset(sub2, "Phase_ciall", &cvec_phi);
  ReadHDF5RealVectorDataset(sub2, "AmpPrefac_ci", &cvec_amp_pre);
  // Read ROM basis functions
  ReadHDF5RealMatrixDataset(sub2, "Bamp", &Bamp);
  ReadHDF5RealMatrixDataset(sub2, "Bphase", &Bphi);

  // Check frequency points
  CheckVectorFromHDF5(sub2, "Mf_grid_Amp", gA_sub2, nk_amp_sub2);
  CheckVectorFromHDF5(sub2, "Mf_grid_Phi", gPhi_sub2, nk_phi_sub2);
  // Check parameter space nodes
  CheckVectorFromHDF5(sub2, "etavec", etavec_sub2, ncx_sub2-2);
  CheckVectorFromHDF5(sub2, "chi1vec", chi1vec_sub2, ncy_sub2-2);
  CheckVectorFromHDF5(sub2, "chi2vec", chi2vec_sub2, ncz_sub2-2);

  XLALH5FileClose(file);
  XLALFree(path);
#else
  // fall back to reading gsl binary data
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub2_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub2_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub2_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub2_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub2_AmpPrefac_ci.dat", cvec_amp_pre);
#endif
  return(ret);
}

// Read binary ROM data for basis functions and coefficients for submodel 3
static int load_data_sub3(const char dir[], gsl_vector *cvec_amp, gsl_vector *cvec_phi, gsl_matrix *Bamp, gsl_matrix *Bphi, gsl_vector *cvec_amp_pre) {
  // Load binary data for amplitude and phase spline coefficients and reduced bases as computed in Mathematica
  // "High frequency, high chi1" submodel 3
  // B-spline points: 97x119x51
  // frequency points: {113, 113}

  // Read data into preallocated gsl_vectors and gsl_matrices
  int ret = XLAL_SUCCESS;
#ifdef LAL_HDF5_ENABLED
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);

  LALH5File *file = XLALH5FileOpen(path, "r");
  LALH5File *sub3 = XLALH5GroupOpen(file, "sub3");

  // Read ROM coefficients
  ReadHDF5RealVectorDataset(sub3, "Amp_ciall", &cvec_amp);
  ReadHDF5RealVectorDataset(sub3, "Phase_ciall", &cvec_phi);
  ReadHDF5RealVectorDataset(sub3, "AmpPrefac_ci", &cvec_amp_pre);
  // Read ROM basis functions
  ReadHDF5RealMatrixDataset(sub3, "Bamp", &Bamp);
  ReadHDF5RealMatrixDataset(sub3, "Bphase", &Bphi);

  // Check frequency points
  CheckVectorFromHDF5(sub3, "Mf_grid_Amp", gA_sub3, nk_amp_sub3);
  CheckVectorFromHDF5(sub3, "Mf_grid_Phi", gPhi_sub3, nk_phi_sub3);
  // Check parameter space nodes
  CheckVectorFromHDF5(sub3, "etavec", etavec_sub3, ncx_sub3-2);
  CheckVectorFromHDF5(sub3, "chi1vec", chi1vec_sub3, ncy_sub3-2);
  CheckVectorFromHDF5(sub3, "chi2vec", chi2vec_sub3, ncz_sub3-2);

  XLALH5FileClose(file);
  XLALFree(path);
#else
  // fall back to reading gsl binary data
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub3_Amp_ciall.dat", cvec_amp);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub3_Phase_ciall.dat", cvec_phi);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub3_Bamp_bin.dat", Bamp);
  ret |= read_matrix(dir, "SEOBNRv2ROM_DS_HI_sub3_Bphase_bin.dat", Bphi);
  ret |= read_vector(dir, "SEOBNRv2ROM_DS_HI_sub3_AmpPrefac_ci.dat", cvec_amp_pre);
#endif
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

  (*splinedata)=malloc(sizeof(SplineData));

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
  free(splinedata);
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
  int nk_max,               // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
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
    *submodel = malloc(sizeof(SEOBNRROMdataDS_submodel));
  else
    SEOBNRROMdataDS_Cleanup_submodel(*submodel);

  int N = ncx*ncy*ncz; // Total number of points over parameter space = size of the data matrix for one SVD-mode

  // Initalize actual ROM data
  // Note: At the moment the sizes of the vectors & matrices are hardwired; they could also be read from HDF5.
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
void SEOBNRROMdataDS_Cleanup_submodel(SEOBNRROMdataDS_submodel *submodel) {
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

#ifdef LAL_HDF5_ENABLED
  // First, check we got the correct version number
  size_t size = strlen(dir) + strlen(ROMDataHDF5) + 2;
  char *path = XLALMalloc(size);
  snprintf(path, size, "%s/%s", dir, ROMDataHDF5);
  LALH5File *file = XLALH5FileOpen(path, "r");

  ret = ROM_check_version_number(file, ROMDataHDF5_VERSION_MAJOR, ROMDataHDF5_VERSION_MINOR, ROMDataHDF5_VERSION_MICRO);
  PrintInfoStringAttribute(file, "Email");
  PrintInfoStringAttribute(file, "Description");

  XLALFree(path);
  XLALH5FileClose(file);
#endif

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
void SEOBNRROMdataDS_Cleanup(SEOBNRROMdataDS *romdata) {
  SEOBNRROMdataDS_Cleanup_submodel((romdata)->sub1);
  free((romdata)->sub1);
  (romdata)->sub1 = NULL;
  romdata->setup=0;
}

/* Structure for internal use */
static void SEOBNRROMdataDS_coeff_Init(SEOBNRROMdataDS_coeff **romdatacoeff, int nk_amp, int nk_phi) {
  if(!romdatacoeff) exit(1);
  /* Create storage for structures */
  if(!*romdatacoeff)
    *romdatacoeff=malloc(sizeof(SEOBNRROMdataDS_coeff));
  else
    SEOBNRROMdataDS_coeff_Cleanup(*romdatacoeff);

  (*romdatacoeff)->c_amp = gsl_vector_alloc(nk_amp);
  (*romdatacoeff)->c_phi = gsl_vector_alloc(nk_phi);
}

/* Deallocate contents of the given SEOBNRROMdataDS_coeff structure */
static void SEOBNRROMdataDS_coeff_Cleanup(SEOBNRROMdataDS_coeff *romdatacoeff) {
  if(romdatacoeff->c_amp) gsl_vector_free(romdatacoeff->c_amp);
  if(romdatacoeff->c_phi) gsl_vector_free(romdatacoeff->c_phi);
  free(romdatacoeff);
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
  double amp_pre_lo,
  double amp_pre_hi,
  const double Mfm,
  // OUTPUTS
  gsl_interp_accel **acc_amp,
  gsl_spline **spline_amp
);

static void GlueAmplitude(
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
) {
  // First need to find overlaping frequency interval
  int jA_lo;
  // Find index so that Mf < Mfm
  for (jA_lo=0; jA_lo < submodel_lo->nk_amp; jA_lo++) {
    if (submodel_lo->gA[jA_lo] > Mfm) {
      jA_lo--;
      break;
    }
  }

  int jA_hi;
  // Find index so that Mf > Mfm
  for (jA_hi=0; jA_hi < submodel_hi->nk_amp; jA_hi++)
    if (submodel_hi->gA[jA_hi] > Mfm)
      break;

  int nA = 1 + jA_lo + (submodel_hi->nk_amp - jA_hi); // length of the union of frequency points of the low and high frequency models glued at MfM

  gsl_vector *gAU = gsl_vector_alloc(nA); // glued frequency grid
  gsl_vector *amp_f = gsl_vector_alloc(nA); // amplitude on glued frequency grid
  // Note: We don't interpolate the amplitude, but this may already be smooth enough for practical purposes.
  // To improve this we would evaluate both amplitue splines times the prefactor at the matching frequency and correct with the ratio, so we are C^0.
  for (int i=0; i<=jA_lo; i++) {
    gsl_vector_set(gAU, i, submodel_lo->gA[i]);
    double A = amp_pre_lo * gsl_vector_get(amp_f_lo, i);
    gsl_vector_set(amp_f, i, A);
  }

  for (int i=jA_lo+1; i<nA; i++) {
    int k = jA_hi - (jA_lo+1) + i;
    gsl_vector_set(gAU, i, submodel_hi->gA[k]);
    double A = amp_pre_hi * gsl_vector_get(amp_f_hi, k);
    gsl_vector_set(amp_f, i, A);
  }

  // Setup 1d splines in frequency from glued amplitude grids & data
  *acc_amp = gsl_interp_accel_alloc();
  *spline_amp = gsl_spline_alloc(gsl_interp_cspline, nA);
  //gsl_spline_init(spline_amp, gAU->data, gsl_vector_const_ptr(amp_f,0), nA);
  gsl_spline_init(*spline_amp, gsl_vector_const_ptr(gAU,0), gsl_vector_const_ptr(amp_f,0), nA);

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
    if (submodel_lo->gPhi[jP_lo] > Mfm) {
      jP_lo--;
      break;
    }
  }

  int jP_hi;
  // Find index so that Mf > Mfm
  for (jP_hi=0; jP_hi < submodel_hi->nk_phi; jP_hi++)
    if (submodel_hi->gPhi[jP_hi] > Mfm)
      break;

  int nP = 1 + jP_lo + (submodel_hi->nk_phi - jP_hi); // length of the union of frequency points of the low and high frequency models glued at MfM
  gsl_vector *gPU = gsl_vector_alloc(nP); // glued frequency grid
  gsl_vector *phi_f = gsl_vector_alloc(nP); // phase on glued frequency grid
  // We need to do a bit more work to glue the phase with C^1 smoothness
  for (int i=0; i<=jP_lo; i++) {
    gsl_vector_set(gPU, i, submodel_lo->gPhi[i]);
    double P = gsl_vector_get(phi_f_lo, i);
    gsl_vector_set(phi_f, i, P);
  }

  for (int i=jP_lo+1; i<nP; i++) {
    int k = jP_hi - (jP_lo+1) + i;
    gsl_vector_set(gPU, i, submodel_hi->gPhi[k]);
    double P = gsl_vector_get(phi_f_hi, k);
    gsl_vector_set(phi_f, i, P);
  }

  // Set up phase data across the gluing frequency Mfm
  // We need to set up a spline for the low frequency model and evaluate at the designated points for the high frequency model so that we work with the *same* frequeny interval!

  // We could optimize this further by not constructing the whole spline for
  // submodel_lo, but this may be insignificant since the number of points is small anyway.
  gsl_interp_accel *acc_phi_lo = gsl_interp_accel_alloc();
  gsl_spline *spline_phi_lo = gsl_spline_alloc(gsl_interp_cspline, submodel_lo->nk_phi);
  gsl_spline_init(spline_phi_lo, submodel_lo->gPhi, gsl_vector_const_ptr(phi_f_lo,0), submodel_lo->nk_phi);

  const int nn = 15;
  gsl_vector_const_view gP_hi_data = gsl_vector_const_view_array(submodel_hi->gPhi + jP_hi - nn, 2*nn+1);
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
    double f = submodel_hi->gPhi[k];
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
  gsl_spline_init(*spline_phi, gsl_vector_const_ptr(gPU,0), gsl_vector_const_ptr(phi_f,0), nP);

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
  SEOBNRROMdataDS *romdata=&__lalsim_SEOBNRv2ROMDS_data;
  if(*hptilde || *hctilde)
  {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }
  int retcode=0;

  // 'Nudge' parameter values to allowed boundary values if close by
  if (eta > 0.25)     nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01)     nudge(&eta, 0.01, 1e-6);
  if (chi1 < -0.9999) nudge(&chi1, -0.9999, 1e-4);
  if (chi1 > 0.98999) nudge(&chi1, 0.98999, 1e-4);
  if (chi2 < -0.9999) nudge(&chi2, -0.9999, 1e-4);
  if (chi2 > 0.98999) nudge(&chi2, 0.98999, 1e-4);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 0.99 || chi2 > 0.99) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1.0 or larger than 0.99!\nSEOBNRv2ROMDoubleSpinHI is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta<0.01 || eta > 0.25) {
    XLALPrintError( "XLAL Error - %s: eta (%f) smaller than 0.01 or unphysical!\nSEOBNRv2ROMDoubleSpin is only available for spins in the range 0.01 <= eta <= 0.25.\n", __func__,eta);
    XLAL_ERROR( XLAL_EDOM );
  }

  /* We always need to glue two submodels together for this ROM */
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM
  submodel_lo = romdata->sub1;
  /* Select high frequency ROM submodel */
  if (chi1 < 0.41)
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
  double Mf_ROM_min = fmax(submodel_lo->gA[0], submodel_lo->gPhi[0]);                                   // lowest allowed geometric frequency for ROM
  double Mf_ROM_max = fmin(submodel_hi->gA[submodel_hi->nk_amp-1], submodel_hi->gPhi[submodel_hi->nk_phi-1]); // highest allowed geometric frequency for ROM
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

  /* Internal storage for waveform coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff_lo=NULL;
  SEOBNRROMdataDS_coeff *romdata_coeff_hi=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_lo, submodel_lo->nk_amp, submodel_lo->nk_phi);
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_hi, submodel_hi->nk_amp, submodel_hi->nk_phi);
  REAL8 amp_pre_lo, amp_pre_hi;

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_lo->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_lo->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_lo->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel_lo->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_lo->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_lo->ncx,             // Number of points in eta  + 2
    submodel_lo->ncy,             // Number of points in chi1 + 2
    submodel_lo->ncz,             // Number of points in chi2 + 2
    submodel_lo->etavec,          // B-spline knots in eta
    submodel_lo->chi1vec,         // B-spline knots in chi1
    submodel_lo->chi2vec,         // B-spline knots in chi2
    romdata_coeff_lo->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_lo->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre_lo                   // Output: interpolated amplitude prefactor
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
    submodel_hi->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel_hi->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_hi->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_hi->ncx,             // Number of points in eta  + 2
    submodel_hi->ncy,             // Number of points in chi1 + 2
    submodel_hi->ncz,             // Number of points in chi2 + 2
    submodel_hi->etavec,          // B-spline knots in eta
    submodel_hi->chi1vec,         // B-spline knots in chi1
    submodel_hi->chi2vec,         // B-spline knots in chi2
    romdata_coeff_hi->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_hi->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre_hi                   // Output: interpolated amplitude prefactor
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
    UINT4 iStart = (UINT4) ceil(fLow_geom / deltaF_geom);
    UINT4 iStop = (UINT4) ceil(fHigh_geom / deltaF_geom);
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

  XLALUnitDivide(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);
  XLALUnitDivide(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);

  COMPLEX16 *pdata=(*hptilde)->data->data;
  COMPLEX16 *cdata=(*hctilde)->data->data;

  REAL8 cosi = cos(inclination);
  REAL8 pcoef = 0.5*(1.0 + cosi*cosi);
  REAL8 ccoef = cosi;

  REAL8 s = 0.5; // Scale polarization amplitude so that strain agrees with FFT of SEOBNRv2
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
    COMPLEX16 htilde = s*amp0*A * cexp(I*phase);
    pdata[j] =      pcoef * htilde;
    cdata[j] = -I * ccoef * htilde;
  }

  /* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */
  /* JV: disable this so as not to clutter logs */
  // if (deltaF > 0)
  //   XLAL_PRINT_WARNING("Warning: Depending on specified frequency sequence correction to time of coalescence may not be accurate.\n");

  // Get SEOBNRv2 ringdown frequency for 22 mode
  double Mf_final = SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(Mtot_sec, eta, chi1, chi2, SEOBNRv2);

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
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_lo);
  SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);

  return(XLAL_SUCCESS);
}

/**
 * @addtogroup LALSimIMRSEOBNRv2ROMDoubleSpin_c
 *
 * \author Michael Puerrer
 *
 * \brief C code for SEOBNRv2 reduced order model
 * (double spin high resolution low mass version).
 * See CQG 31 195010, 2014, arXiv:1402.4146 for the basic approach.
 * Further details in PRD 93, 064041, 2016, arXiv:1512.02248.
 *
 * This is a frequency domain model that approximates the time domain SEOBNRv2 model.
 *
 * The binary data HDF5 file (SEOBNRv2ROM_DS_HI_vXYZ.hdf5) and the gsl-binary data files (SEOBNRv2ROM_DS_HI_vXYZ.tar)
 * will be available at on LIGO clusters in /home/cbc/.
 * Make sure the files are in your LAL_DATA_PATH.
 *
 * @note Note that due to its construction the iFFT of the ROM has a small (~ 20 M) offset
 * in the peak time that scales with total mass as compared to the time-domain SEOBNRv2 model.
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
 * @{
 */


/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv2_ROM_DoubleSpin_HI model.
 *
 * XLALSimIMRSEOBNRv2ROMDoubleSpinHI() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence(
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
  UINT4 nk_max)                                 /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
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
  (void) pthread_once(&SEOBNRv2ROMDoubleSpin_is_initialized, SEOBNRv2ROMDoubleSpin_Init_LALDATA);
#else
  SEOBNRv2ROMDoubleSpin_Init_LALDATA();
#endif

  if(!SEOBNRv2ROMDoubleSpin_IsSetup()) XLAL_ERROR(XLAL_EFAILED,"Error setting up SEOBNRv2ROMDoubleSpinHI data - check your $LAL_DATA_PATH\n");

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv2ROMDoubleSpinCore(hptilde,hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi1, chi2, freqs, 0, nk_max);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv2_ROM_DoubleSpin_HI model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinHI(
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
  UINT4 nk_max)                                 /**< Truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1 */
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

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv2ROMDoubleSpinCore(hptilde,hctilde,
            phiRef, fRef, distance, inclination, Mtot_sec, eta, chi1, chi2, freqs, deltaF, nk_max);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

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
  if (eta > 0.25)     nudge(&eta, 0.25, 1e-6);
  if (eta < 0.01)     nudge(&eta, 0.01, 1e-6);
  if (chi1 < -0.9999) nudge(&chi1, -0.9999, 1e-4);
  if (chi1 > 0.98999) nudge(&chi1, 0.98999, 1e-4);
  if (chi2 < -0.9999) nudge(&chi2, -0.9999, 1e-4);
  if (chi2 > 0.98999) nudge(&chi2, 0.98999, 1e-4);

  if ( chi1 < -1.0 || chi2 < -1.0 || chi1 > 0.99 || chi2 > 0.99) {
    XLALPrintError( "XLAL Error - %s: chi1 or chi2 smaller than -1.0 or larger than 0.99!\nSEOBNRv2ROMDoubleSpinHI is only available for spins in the range -1 <= a/M <= 0.99.\n", __func__);
    XLAL_ERROR( XLAL_EDOM );
  }

  if (eta < 0.01 || eta > 0.25) {
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

  /* We always need to glue two submodels together for this ROM */
  SEOBNRROMdataDS_submodel *submodel_hi; // high frequency ROM
  SEOBNRROMdataDS_submodel *submodel_lo; // low frequency ROM
  submodel_lo = romdata->sub1;
  /* Select high frequency ROM submodel */
  if (chi1 < 0.41)
    submodel_hi = romdata->sub2;
  else
    submodel_hi = romdata->sub3;

  /* Internal storage for waveform coefficiencts */
  SEOBNRROMdataDS_coeff *romdata_coeff_lo=NULL;
  SEOBNRROMdataDS_coeff *romdata_coeff_hi=NULL;
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_lo, submodel_lo->nk_amp, submodel_lo->nk_phi);
  SEOBNRROMdataDS_coeff_Init(&romdata_coeff_hi, submodel_hi->nk_amp, submodel_hi->nk_phi);
  REAL8 amp_pre_lo, amp_pre_hi;

  /* Interpolate projection coefficients and evaluate them at (eta,chi1,chi2) */
  int nk_max = -1; // adjust truncation parameter if speed is an issue
  int retcode=TP_Spline_interpolation_3d(
    eta,                          // Input: eta-value for which projection coefficients should be evaluated
    chi1,                         // Input: chi1-value for which projection coefficients should be evaluated
    chi2,                         // Input: chi2-value for which projection coefficients should be evaluated
    submodel_lo->cvec_amp,        // Input: data for spline coefficients for amplitude
    submodel_lo->cvec_phi,        // Input: data for spline coefficients for phase
    submodel_lo->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel_lo->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_lo->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_lo->ncx,             // Number of points in eta  + 2
    submodel_lo->ncy,             // Number of points in chi1 + 2
    submodel_lo->ncz,             // Number of points in chi2 + 2
    submodel_lo->etavec,          // B-spline knots in eta
    submodel_lo->chi1vec,         // B-spline knots in chi1
    submodel_lo->chi2vec,         // B-spline knots in chi2
    romdata_coeff_lo->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_lo->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre_lo                   // Output: interpolated amplitude prefactor
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
    submodel_hi->cvec_amp_pre,    // Input: data for spline coefficients for amplitude prefactor
    submodel_hi->nk_amp,          // number of SVD-modes == number of basis functions for amplitude
    submodel_hi->nk_phi,          // number of SVD-modes == number of basis functions for phase
    nk_max,                       // truncate interpolants at SVD mode nk_max; don't truncate if nk_max == -1
    submodel_hi->ncx,             // Number of points in eta  + 2
    submodel_hi->ncy,             // Number of points in chi1 + 2
    submodel_hi->ncz,             // Number of points in chi2 + 2
    submodel_hi->etavec,          // B-spline knots in eta
    submodel_hi->chi1vec,         // B-spline knots in chi1
    submodel_hi->chi2vec,         // B-spline knots in chi2
    romdata_coeff_hi->c_amp,      // Output: interpolated projection coefficients for amplitude
    romdata_coeff_hi->c_phi,      // Output: interpolated projection coefficients for phase
    &amp_pre_hi                   // Output: interpolated amplitude prefactor
  );

  if(retcode!=0) {
    SEOBNRROMdataDS_coeff_Cleanup(romdata_coeff_hi);
    XLAL_ERROR(retcode);
  }

  // Compute function values of amplitude an phase on sparse frequency points by evaluating matrix vector products
  // phi_pts = B_phi^T . c_phi
  gsl_vector* phi_f_lo = gsl_vector_alloc(submodel_lo->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_lo->Bphi, romdata_coeff_lo->c_phi, 0.0, phi_f_lo);

  gsl_vector* phi_f_hi = gsl_vector_alloc(submodel_hi->nk_phi);
  gsl_blas_dgemv(CblasTrans, 1.0, submodel_hi->Bphi, romdata_coeff_hi->c_phi, 0.0, phi_f_hi);

  const double Mfm = 0.01; // Gluing frequency: the low and high frequency ROMs overlap here; this is used both for amplitude and phase.

  // Glue phasing in frequency to C^1 smoothness
  GluePhasing(submodel_lo, submodel_hi, phi_f_lo, phi_f_hi, Mfm,
    acc_phi, spline_phi
  );

  // Get SEOBNRv2 ringdown frequency for 22 mode
  *Mf_final = SEOBNRROM_Ringdown_Mf_From_Mtot_Eta(*Mtot_sec, eta, chi1, chi2, SEOBNRv2);

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
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * The allowed frequency range for the starting frequency in geometric frequency is [0.00053, 0.135].
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencyOfTime() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency(
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
  double Mf_ROM_min = 0.0000985;
  double Mf_ROM_max = 0.3;

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  double t_corr = gsl_spline_eval_deriv(spline_phi, Mf_final, acc_phi) / (2*LAL_PI); // t_corr / M
  XLAL_PRINT_INFO("t_corr[s] = %g\n", t_corr * Mtot_sec);

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
 * chirp time, but it includes both the inspiral and the merger ringdown part of SEOBNRv2.
 *
 * If the frequency that corresponds to the specified elapsed time is lower than the
 * geometric frequency Mf=0.00053 (ROM starting frequency) or above half of the SEOBNRv2
 * ringdown frequency an error is thrown.
 * The SEOBNRv2 ringdown frequency can be obtained by calling XLALSimInspiralGetFinalFreq().
 *
 * See XLALSimIMRSEOBNRv2ROMDoubleSpinHITimeOfFrequency() for the inverse function.
 */
int XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencyOfTime(
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


/** Setup SEOBNRv2ROMDoubleSpin model using data files installed in $LAL_DATA_PATH
 */
void SEOBNRv2ROMDoubleSpin_Init_LALDATA(void)
{
  if (SEOBNRv2ROMDoubleSpin_IsSetup()) return;

  // For gsl binary data: If we find one ROM datafile in a directory listed in LAL_DATA_PATH,
  // then we expect the remaining datafiles to also be there.
  // For HDF5 there is only one file.
#ifdef LAL_HDF5_ENABLED
#define datafile ROMDataHDF5
#else
  const char datafile[] = "SEOBNRv2ROM_DS_HI_sub1_Phase_ciall.dat";
#endif

  char *path = XLALFileResolvePathLong(datafile, PKG_DATA_DIR);
  if (path==NULL)
    XLAL_ERROR_VOID(XLAL_EIO, "Unable to resolve data file %s in $LAL_DATA_PATH\n", datafile);
  char *dir = dirname(path);
  int ret = SEOBNRv2ROMDoubleSpin_Init(dir);
  XLALFree(path);

  if(ret!=XLAL_SUCCESS)
    XLAL_ERROR_VOID(XLAL_FAILURE, "Unable to find SEOBNRv2ROMDoubleSpin_HI data files in $LAL_DATA_PATH\n");
}
