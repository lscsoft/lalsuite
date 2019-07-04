 /*
 *  Copyright (C) 2019 Edward Jacob Fauchon-Jones, Jonathan E. Thompson, Sebastian Khan
 *  Test code for LALSimIMRPhenomNSBH
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include <lal/Units.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/LALConstants.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#include <lal/LALSimNoise.h>
#include <lal/ComplexFFT.h>

#include <lal/ComplexFFT.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define MYUNUSED(expr) do { (void)(expr); } while (0)

bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon);
bool approximatelyEqual(REAL8 a, REAL8 b, REAL8 epsilon) {
  if (a == 0)
    return fabs(b) < epsilon;
  else if (b == 0)
    return fabs(a) < epsilon;
  else
    return !gsl_fcmp(a, b, epsilon);
}

void print_difference(const char *name, REAL8 u, REAL8 u_expected);
void print_difference(const char *name, REAL8 u, REAL8 u_expected) {
  printf("%s: %-20.17g\t%-20.17g\t%-20.17g\n", name, u, u_expected, u - u_expected);
}

#define COUNT_FGWINKERR 9
#define TOLERANCE_FGWINKERR 1e-8
double input_fGWinKerr[COUNT_FGWINKERR][3] = {
  {8.7173522796064891, 1.0, -0.900}, {13.0760284194097345, 1.0, -0.900},
  {17.4347045592129781, 1.0, -0.900}, {6.0000000000000000, 1.0, 0.000},
  {9.0000000000000000, 1.0, 0.000}, {12.0000000000000000, 1.0, 0.000},
  {2.3208830417618871, 1.0, 0.900}, {3.4813245626428309, 1.0, 0.900},
  {4.6417660835237742, 1.0, 0.900}};
double expected_fGWinKerr[COUNT_FGWINKERR] = {
  0.01281537532531727, 0.00686250007959751, 0.00442721737694587,
  0.02165824447871323, 0.0117892550438441, 0.007657345769747112,
  0.07176032468044616, 0.04304065885628561, 0.0292012072733511};
static void Test_fGWinKerr(void);
static void Test_fGWinKerr(void) {
  printf("\n## Test_fGWinKerr\n\n");
  for (int i=0; i<COUNT_FGWINKERR; i++) {
    double output = XLALSimNSBH_fGWinKerr(
      input_fGWinKerr[i][0],
      input_fGWinKerr[i][1],
      input_fGWinKerr[i][2]);
    print_difference("XLALSimNSBH_fGWinKerr", output, expected_fGWinKerr[i]);
    assert(approximatelyEqual(output, expected_fGWinKerr[i], TOLERANCE_FGWINKERR));
  }
}

#define COUNT_RKERRISCO 7
#define TOLERANCE_RKERRISCO 1e-8
double input_rKerrISCO[COUNT_RKERRISCO][1] = {
  {-0.90}, {-0.60}, {-0.30}, {-0.00}, {0.30}, {0.60}, {0.90}};
double expected_rKerrISCO[COUNT_RKERRISCO] = {
  8.717352279606489, 7.850686185306578, 6.949272527004718, 6, 4.97861683057595,
  3.829069418813151, 2.320883041761887};
static void Test_rKerrISCO(void);
static void Test_rKerrISCO(void) {
  printf("\n## Test_rKerrISCO\n\n");
  for (int i=0; i<COUNT_RKERRISCO; i++) {
    double output = XLALSimNSBH_rKerrISCO(
      input_rKerrISCO[i][0]);
    print_difference("XLALSimNSBH_rKerrISCO", output, expected_rKerrISCO[i]);
    assert(approximatelyEqual(output, expected_rKerrISCO[i], TOLERANCE_RKERRISCO));
  }
}

#define COUNT_XI_TIDE 36
#define TOLERANCE_XI_TIDE 1e-8
double input_xi_tide[COUNT_XI_TIDE][3] = {
  {2.00, -0.90, 0.28}, {2.00, -0.90, 0.34}, {2.00, -0.90, 0.4}, {2.00, 0.00, 0.28},
  {2.00, 0.00, 0.34}, {2.00, 0.00, 0.4}, {2.00, 0.90, 0.28}, {2.00, 0.90, 0.34}, {2.00, 0.90, 0.4},
  {3.00, -0.90, 0.42}, {3.00, -0.90, 0.51}, {3.00, -0.90, 0.6}, {3.00, 0.00, 0.42},
  {3.00, 0.00, 0.51}, {3.00, 0.00, 0.6}, {3.00, 0.90, 0.42}, {3.00, 0.90, 0.51}, {3.00, 0.90, 0.6},
  {4.00, -0.90, 0.56}, {4.00, -0.90, 0.68}, {4.00, -0.90, 0.8}, {4.00, 0.00, 0.56},
  {4.00, 0.00, 0.68}, {4.00, 0.00, 0.8}, {4.00, 0.90, 0.56}, {4.00, 0.90, 0.68}, {4.00, 0.90, 0.8},
  {5.00, -0.90, 0.7}, {5.00, -0.90, 0.85}, {5.00, -0.90, 1}, {5.00, 0.00, 0.7}, {5.00, 0.00, 0.85},
  {5.00, 0.00, 1}, {5.00, 0.90, 0.7}, {5.00, 0.90, 0.85}, {5.00, 0.90, 1}};
double expected_xi_tide[COUNT_XI_TIDE] = {
  2.06871189988092, 2.165647459230484, 2.279729557247444, 1.957702074572125,
  2.005770585734644, 2.063058314326985, 1.874012683952472, 1.884478196440186,
  1.894984431945135, 2.535178557740413, 2.72161470061767, 2.939474255173244,
  2.324039062239011, 2.419350068761182, 2.5369293431522, 2.162563923478413,
  2.178495207637674, 2.19511940553916, 2.992537647891789, 3.286528350419773,
  3.623185704395681, 2.661277126455292, 2.820959995701167, 3.018859346483027,
  2.397505338613715, 2.419802208282962, 2.444234246009384, 3.455945198308276,
  3.868295212767634, 4.330079575392291, 2.991668576796354, 3.230879137923294,
  3.522510111613574, 2.600443557587901, 2.630444772711848, 2.665283137995229};
static void Test_xi_tide(void);
static void Test_xi_tide(void) {
  printf("\n## Test_xi_tide\n\n");
  for (int i=0; i<COUNT_XI_TIDE; i++) {
    double output = XLALSimNSBH_xi_tide(
      input_xi_tide[i][0],
      input_xi_tide[i][1],
      input_xi_tide[i][2]);
    print_difference("XLALSimNSBH_xi_tide", output, expected_xi_tide[i]);
    assert(approximatelyEqual(output, expected_xi_tide[i], TOLERANCE_XI_TIDE));
  }
}

#define COUNT_COMPACTNESS_FROM_LAMBDA 21
#define TOLERANCE_COMPACTNESS_FROM_LAMBDA 1e-8
double input_compactness_from_lambda[COUNT_COMPACTNESS_FROM_LAMBDA][1] = {
  {0.0}, {200.0}, {400.0}, {600.0}, {800.0}, {1000.0}, {1200.0}, {1400.0}, {1600.0}, {1800.0},
  {2000.0}, {2200.0}, {2400.0}, {2600.0}, {2800.0}, {3000.0}, {3200.0}, {3400.0}, {3600.0}, {3800.0},
  {4000.0}};
double expected_compactness_from_lambda[COUNT_COMPACTNESS_FROM_LAMBDA] = {
  0.5, 0.1917006111637932, 0.1726108500082392, 0.1617580970945444,
  0.1541985276023099, 0.1484152311071196, 0.1437420509536809,
  0.1398274999660475, 0.1364636439460052, 0.1335173810600413,
  0.1308984340168845, 0.1285427920411775, 0.1264034435624421,
  0.1244448606165661, 0.1226395498186449, 0.1209658060963831,
  0.1194061990393252, 0.1179465229202867, 0.1165750499338684,
  0.1152819874103816, 0.114059075676274};
static void Test_compactness_from_lambda(void);
static void Test_compactness_from_lambda(void) {
  printf("\n## Test_compactness_from_lambda\n\n");
  for (int i=0; i<COUNT_COMPACTNESS_FROM_LAMBDA; i++) {
    double output = XLALSimNSBH_compactness_from_lambda(
      input_compactness_from_lambda[i][0]);
    print_difference("XLALSimNSBH_compactness_from_lambda", output, expected_compactness_from_lambda[i]);
    assert(approximatelyEqual(output, expected_compactness_from_lambda[i], TOLERANCE_COMPACTNESS_FROM_LAMBDA));
  }
}

#define COUNT_TORUS_MASS_FIT 36
#define TOLERANCE_TORUS_MASS_FIT 1e-8
double input_torus_mass_fit[COUNT_TORUS_MASS_FIT][3] = {
  {2.00, -0.90, 0.14}, {2.00, -0.90, 0.17}, {2.00, -0.90, 0.20}, {2.00, 0.00, 0.14},
  {2.00, 0.00, 0.17}, {2.00, 0.00, 0.20}, {2.00, 0.90, 0.14}, {2.00, 0.90, 0.17}, {2.00, 0.90, 0.20},
  {3.00, -0.90, 0.14}, {3.00, -0.90, 0.17}, {3.00, -0.90, 0.20}, {3.00, 0.00, 0.14},
  {3.00, 0.00, 0.17}, {3.00, 0.00, 0.20}, {3.00, 0.90, 0.14}, {3.00, 0.90, 0.17}, {3.00, 0.90, 0.20},
  {4.00, -0.90, 0.14}, {4.00, -0.90, 0.17}, {4.00, -0.90, 0.20}, {4.00, 0.00, 0.14},
  {4.00, 0.00, 0.17}, {4.00, 0.00, 0.20}, {4.00, 0.90, 0.14}, {4.00, 0.90, 0.17}, {4.00, 0.90, 0.20},
  {5.00, -0.90, 0.14}, {5.00, -0.90, 0.17}, {5.00, -0.90, 0.20}, {5.00, 0.00, 0.14},
  {5.00, 0.00, 0.17}, {5.00, 0.00, 0.20}, {5.00, 0.90, 0.14}, {5.00, 0.90, 0.17}, {5.00, 0.90, 0.20}};
double expected_torus_mass_fit[COUNT_TORUS_MASS_FIT] = {
  0.02349705295506288, 0, 0, 0.1299454661328111, 0.0430073416291199, 0,
  0.2882657031643916, 0.2332155204085186, 0.1778008350569429, 0, 0, 0,
  0.06437920494437788, 0, 0, 0.2941998033123807, 0.2231866136920418,
  0.1517306063389851, 0, 0, 0, 0, 0, 0, 0.2887085776862365, 0.2028602793140871,
  0.1165992019782404, 0, 0, 0, 0, 0, 0, 0.2763968308942356, 0.1765433406768961,
  0.07648328516666991};
static void Test_torus_mass_fit(void);
static void Test_torus_mass_fit(void) {
  printf("\n## Test_torus_mass_fit\n\n");
  for (int i=0; i<COUNT_TORUS_MASS_FIT; i++) {
    double output = XLALSimNSBH_torus_mass_fit(
      input_torus_mass_fit[i][0],
      input_torus_mass_fit[i][1],
      input_torus_mass_fit[i][2]);
    print_difference("XLALSimNSBH_torus_mass_fit", output, expected_torus_mass_fit[i]);
    assert(approximatelyEqual(output, expected_torus_mass_fit[i], TOLERANCE_TORUS_MASS_FIT));
  }
}

int main(int argc, char *argv[]) {
  MYUNUSED(argc);
  MYUNUSED(argv);

  Test_fGWinKerr();
  Test_rKerrISCO();

  Test_xi_tide();
  Test_compactness_from_lambda();
  Test_torus_mass_fit();

  return 0;
}
