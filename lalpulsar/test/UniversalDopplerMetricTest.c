/*
 * Copyright (C) 2009 Reinhard Prix
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

#include <math.h>
#include <sys/times.h>

#include <lal/LALMalloc.h>
#include <lal/LALInitBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/LALConstants.h>
#include <lal/StringVector.h>
#include <lal/LogPrintf.h>

#include <lal/PtoleMetric.h>
#include <lal/UniversalDopplerMetric.h>

/** \author Reinhard Prix
 * \file
 * \ingroup UniversalDopplerMetric_h
 * \brief Tests for exported functions in UniversalDopplerMetric
 *
 */

// ---------- defines --------------------

// ---------- Macros --------------------

// copy 3 components of Euklidean vector
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)
// compute norm of Euklidean 3-vector
#define NORM(v) ( sqrt ( (v)[0]*(v)[0] + (v)[1]*(v)[1] + (v)[2]*(v)[2] )  )
// compute scalar division of Euklidean 3-vector by div
#define DIV_VECT(dst,src,div) do { (dst)[0] = (src)[0]/(div); (dst)[1] = (src)[1]/(div); (dst)[2] = (src)[2]/(div); } while(0)
#define SUB_VECT(dst,src) do { (dst)[0] -= (src)[0]; (dst)[1] -= (src)[1]; (dst)[2] -= (src)[2]; } while(0)
#define MULT_VECT(v,lam) do{ (v)[0] *= (lam); (v)[1] *= (lam); (v)[2] *= (lam); } while(0)

// ---------- global variables --------------------
extern int lalDebugLevel;

static LALStatus empty_LALStatus;
static PtoleMetricIn empty_PtoleMetricIn;

// ---------- local prototypes
static int test_XLALComputeOrbitalDerivatives ( void );
static int test_XLALDopplerFstatMetric ( void );

static REAL8 compare_metrics ( const gsl_matrix *m1, const gsl_matrix *m2 );
static gsl_matrix *convert_old_metric_2_new ( const REAL8Vector *m0, REAL8 Freq );

DopplerMetric *XLALOldDopplerFstatMetric ( const DopplerMetricParams *metricParams, const EphemerisData *edat );

// ---------- function definitions --------------------
/** MAIN function: calls a number of unit-tests
 */
int main( void )
{
  INT4 ret;
  lalDebugLevel = 3;

  ret = test_XLALDopplerFstatMetric();
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "test_XLALDopplerFstatMetric() failed with status = %d.\n", ret );

  ret = test_XLALComputeOrbitalDerivatives();
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "test_XLALComputeOrbitalDerivatives() failed with status = %d.\n", ret );

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  /* all tests passed */
  return XLAL_SUCCESS;

} /* main() */


/**
 * Unit test for 'do-it-all' metric function XLALDopplerFstatMetric()
 *
 * Initially modelled afer testMetricCodes.py script:
 * Check metric codes 'getMetric' 'FstatMetric' and 'FstatMetric_v2' by
 * comparing them against each other.
 * Given that they represent 3 very different implementations of
 * metric calculations, this provides a very powerful consistency test
 *
 */
static int
test_XLALDopplerFstatMetric ( void )
{
  int ret;
  REAL8 tolPh = 0.01;	// 1% tolerance on phase metrics [taken from testMetricCodes.py]

  // ----- load ephemeris
  char earthEphem[] = TEST_DATA_DIR "earth00-19-DE200.dat.gz";
  char sunEphem[]   = TEST_DATA_DIR "sun00-19-DE200.dat.gz";
  EphemerisData *edat = XLALInitBarycenter ( earthEphem, sunEphem );
  XLAL_CHECK ( edat != NULL, XLAL_EFUNC, "XLALInitBarycenter('%s','%s') failed with xlalErrno = %d\n", earthEphem, sunEphem, xlalErrno );

  // ----- set test-parameters ----------
  LIGOTimeGPS startTimeGPS = { 792576013, 0 };
  REAL8 duration = 60000;

  REAL8 Alpha = 1.0;
  REAL8 Delta = 0.5;
  REAL8 Freq  = 100;
  REAL8 f1dot = 0;// -1e-8;

  LALStringVector *detNames = XLALCreateStringVector ( "H1", "L1", "V1",  NULL );
  LALStringVector *sqrtSX   = XLALCreateStringVector ( "1.0", "0.5", "1.5", NULL );
  MultiDetectorInfo detInfo;
  XLAL_CHECK ( XLALParseMultiDetectorInfo ( &detInfo, detNames, sqrtSX ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyStringVector ( detNames );
  XLALDestroyStringVector ( sqrtSX );

  // prepare metric parameters for modern XLALDopplerFstatMetric() and mid-old XLALOldDopplerFstatMetric()
  DopplerCoordinateSystem coordSys = { 4, { DOPPLERCOORD_FREQ, DOPPLERCOORD_ALPHA, DOPPLERCOORD_DELTA, DOPPLERCOORD_F1DOT } };
  PulsarAmplitudeParams Amp = { 0.03, -0.3, 0.5, 0.0 };	// h0, cosi, psi, phi0
  PulsarDopplerParams dop = empty_PulsarDopplerParams;
  dop.refTime  = startTimeGPS;
  dop.Alpha    = Alpha;
  dop.Delta    = Delta;
  dop.fkdot[0] = Freq;
  dop.fkdot[1] = f1dot;

  UINT4 Nseg = 1;
  REAL8 Tseg = duration / Nseg;
  LALSegList segList;
  ret = XLALSegListInitSimpleSegments ( &segList, startTimeGPS, Nseg, Tseg );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALSegListInitSimpleSegments() failed with xlalErrno = %d\n", xlalErrno );

  DopplerMetricParams pars2 = empty_DopplerMetricParams;

  pars2.coordSys      		= coordSys;
  pars2.detMotionType 		= DETMOTION_SPIN_ORBIT;
  pars2.segmentList   		= segList;
  pars2.detInfo 		= detInfo;
  pars2.signalParams.Amp     	= Amp;
  pars2.signalParams.Doppler 	= dop;
  pars2.projectCoord  		= - 1;	// -1==no projection
  pars2.metricType    		= METRIC_TYPE_PHASE;
  pars2.approxPhase   		= 0;

  // ----- prepare call-parameters of ancient LALPulsarMetric()
  LALStatus status = empty_LALStatus;
  PtoleMetricIn pars0 = empty_PtoleMetricIn;

  pars0.spindown = XLALCreateREAL4Vector ( 1 );
  XLAL_CHECK ( pars0.spindown != NULL, XLAL_EFUNC, "XLALCreateREAL4Vector(1) failed.\n");

  pars0.position.longitude 	= Alpha;
  pars0.position.latitude  	= Delta;
  pars0.position.system    	= COORDINATESYSTEM_EQUATORIAL;
  pars0.spindown->data[0]	= f1dot / Freq; /* 'old-style' "f1" = f1dot / Freq !!*/
  pars0.epoch    		= startTimeGPS;
  pars0.duration 		= duration;
  pars0.maxFreq  		= Freq;
  pars0.site     		= &(detInfo.sites[0]);	// use first detector for phase-metric
  pars0.ephemeris 		= edat;
  pars0.metricType 		= LAL_PMETRIC_COH_EPHEM;


  // ========== BEGINNING OF TEST CALLS ==========
  REAL8Vector *metric0 = 0;
  gsl_matrix *g0_ij;
  DopplerMetric *metric1, *metric2;
  REAL8 diff_2_0, diff_2_1, diff_1_0;


  XLALPrintWarning("\n---------- ROUND 1: ephemeris-based (single-IFO) phase-metrics ----------\n");

  pars2.detInfo.length = 1;	// truncate to first detector

  // 0) compute metric using ancient LALPulsarMetric() function (used in lalapps_getMetric)
  LALPulsarMetric ( &status, &metric0, &pars0 );
  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALPulsarMetric() failed with status=%d: '%s'\n", status.statusCode, status.statusDescription );
  XLAL_CHECK ( (g0_ij = convert_old_metric_2_new ( metric0, Freq )) != NULL, XLAL_EFUNC );
  XLALDestroyREAL8Vector ( metric0 ); metric0 = NULL;
  // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
  XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
  // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
  XLAL_CHECK ( (metric2 = XLALDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

  // compare all 3 metrics against each other:
  XLAL_CHECK ( (diff_2_0 = compare_metrics ( metric2->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g2,g0)= %g exceeds tolerance of %g\n", diff_2_0, tolPh );
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( metric2->g_ij, metric1->g_ij )) < tolPh, XLAL_ETOL, "Error(g2,g1)= %g exceeds tolerance of %g\n", diff_2_1, tolPh );
  XLAL_CHECK ( (diff_1_0 = compare_metrics ( metric1->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g1,g0)= %g exceeds tolerance of %g\n", diff_1_0, tolPh );

  gsl_matrix_free ( g0_ij );
  XLALDestroyDopplerMetric ( metric1 );
  XLALDestroyDopplerMetric ( metric2 );

  XLALPrintWarning ("diff_2_0 = %e\n", diff_2_0 );
  XLALPrintWarning ("diff_2_1 = %e\n", diff_2_1 );
  XLALPrintWarning ("diff_1_0 = %e\n", diff_1_0 );


  XLALPrintWarning("\n---------- ROUND 2: Ptolemaic-based (single-IFO) phase-metrics ----------\n");

  pars0.metricType    = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
  pars2.detMotionType = DETMOTION_SPIN_PTOLEORBIT;

  // 0) compute metric using ancient LALPulsarMetric() function (used in lalapps_getMetric)
  LALPulsarMetric ( &status, &metric0, &pars0 );
  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALPulsarMetric() failed with status=%d: '%s'\n", status.statusCode, status.statusDescription );
  XLAL_CHECK ( (g0_ij = convert_old_metric_2_new ( metric0, Freq )) != NULL, XLAL_EFUNC );
  XLALDestroyREAL8Vector ( metric0 ); metric0 = NULL;
  // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
  XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
  // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
  XLAL_CHECK ( (metric2 = XLALDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

  // compare all 3 metrics against each other:
  XLAL_CHECK ( (diff_2_0 = compare_metrics ( metric2->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g2,g0)= %e exceeds tolerance of %e\n", diff_2_0, tolPh );
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( metric2->g_ij, metric1->g_ij )) < tolPh, XLAL_ETOL, "Error(g2,g1)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
  XLAL_CHECK ( (diff_1_0 = compare_metrics ( metric1->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g1,g0)= %e exceeds tolerance of %e\n", diff_1_0, tolPh );

  gsl_matrix_free ( g0_ij );
  XLALDestroyDopplerMetric ( metric1 );
  XLALDestroyDopplerMetric ( metric2 );

  XLALPrintWarning ("diff_2_0 = %e\n", diff_2_0 );
  XLALPrintWarning ("diff_2_1 = %e\n", diff_2_1 );
  XLALPrintWarning ("diff_1_0 = %e\n", diff_1_0 );


  XLALPrintWarning("\n---------- ROUND 3: (ephemeris-base) *full* (multi-IF) F-stat metrics ----------\n");

  pars2.detMotionType = DETMOTION_SPIN_ORBIT;
  pars2.metricType    = METRIC_TYPE_FSTAT;
  pars2.detInfo       = detInfo;	// 3 IFOs

  // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
  XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
  // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
  XLAL_CHECK ( (metric2 = XLALDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

  // compare both metrics against each other:
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( metric2->gF_ij,   metric1->gF_ij ))   < tolPh, XLAL_ETOL, "Error(gF2,gF1)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
  XLALPrintWarning ("gF:   diff_2_1 = %e\n", diff_2_1 );
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( metric2->gFav_ij, metric1->gFav_ij )) < tolPh, XLAL_ETOL, "Error(gFav2,gFav1)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
  XLALPrintWarning ("gFav: diff_2_1 = %e\n", diff_2_1 );

  XLALDestroyDopplerMetric ( metric1 );
  XLALDestroyDopplerMetric ( metric2 );


  XLALPrintWarning("\n---------- ROUND 4: compare analytic {f,f1dot,f2dot,f3dot} phase-metric vs  XLALDopplerFstatMetric() ----------\n");
  pars2.detInfo.length  = 1;	// truncate to 1st detector
  pars2.detMotionType   = DETMOTION_SPIN_ORBIT;
  pars2.metricType      = METRIC_TYPE_PHASE;
  pars2.approxPhase     = 1;	// use same phase-approximation as in analytic solution to improve comparison

  DopplerCoordinateSystem coordSys2 = { 4, { DOPPLERCOORD_FREQ, DOPPLERCOORD_F1DOT, DOPPLERCOORD_F2DOT, DOPPLERCOORD_F3DOT } };
  pars2.coordSys = coordSys2;
  gsl_matrix* gN_ij;

  // a) compute metric at refTime = startTime
  pars2.signalParams.Doppler.refTime = startTimeGPS;
  XLAL_CHECK ( (metric2 = XLALDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gN_ij = XLALNaturalizeMetric ( metric2->g_ij, &pars2 )) != NULL, XLAL_EFUNC );

  REAL8 gStart_ij[] = {   1.0/3,      2.0/3,    6.0/5,    32.0/15,      \
                          2.0/3,    64.0/45,    8.0/3,  512.0/105,      \
                          6.0/5,      8.0/3,   36.0/7,     48.0/5,      \
                          32.0/15,  512.0/105, 48.0/5, 4096.0/225 };
  const gsl_matrix_view gStart = gsl_matrix_view_array ( gStart_ij, 4, 4 );

  // compare natural-units metric against analytic solution
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( gN_ij,   &(gStart.matrix) )) < tolPh, XLAL_ETOL,
               "RefTime=StartTime: Error(g_ij,g_analytic)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
  XLALPrintWarning ("Analytic (refTime=startTime): diff_2_1 = %e\n", diff_2_1 );

  XLALDestroyDopplerMetric ( metric2 );
  gsl_matrix_free ( gN_ij );


  // b) compute metric at refTime = midTime
  pars2.signalParams.Doppler.refTime = startTimeGPS;
  pars2.signalParams.Doppler.refTime.gpsSeconds += duration / 2;

  XLAL_CHECK ( (metric2 = XLALDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gN_ij = XLALNaturalizeMetric ( metric2->g_ij, &pars2 )) != NULL, XLAL_EFUNC );

  REAL8 gMid_ij[] = { 1.0/3,    0,        1.0/5,         0,       \
                      0,        4.0/45,       0,   8.0/105,       \
                      1.0/5,    0,        1.0/7,         0,       \
                      0,        8.0/105,      0,  16.0/225  };
  const gsl_matrix_view gMid = gsl_matrix_view_array ( gMid_ij, 4, 4 );

  // compare natural-units metric against analytic solution
  XLAL_CHECK ( (diff_2_1 = compare_metrics ( gN_ij,   &(gMid.matrix) )) < tolPh, XLAL_ETOL,
               "RefTime=MidTime: Error(g_ij,g_analytic)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
  XLALPrintWarning ("Analytic (refTime=midTime):   diff_2_1 = %e\n\n", diff_2_1 );

  XLALDestroyDopplerMetric ( metric2 );
  gsl_matrix_free ( gN_ij );

  // ----- clean up memory
  XLALSegListClear ( &segList );
  XLALDestroyREAL4Vector ( pars0.spindown );
  XLALDestroyEphemerisData ( edat );

  return XLAL_SUCCESS;

} /* test_XLALDopplerFstatMetric() */

// convert 'old-style' REAL8Vector metric into a new 'gsl_matrix' type
// also converts 'old-style' spindown f1 = f1dot / Freq, back into 'f1dot' coordinates
static gsl_matrix *
convert_old_metric_2_new ( const REAL8Vector *m0, REAL8 Freq )
{
  XLAL_CHECK_NULL ( m0, XLAL_EINVAL, "Invalid NULL input 'm0'\n");
  XLAL_CHECK_NULL ( Freq > 0, XLAL_EINVAL, "Invalid input Freq = %g must be > 0\n", Freq );
  XLAL_CHECK_NULL ( m0->length > 0, XLAL_EINVAL, "Invalid zero-length metric input 'm0'\n");

  UINT4 len = m0->length;
  REAL4 dim0 = -0.5 + sqrt ( 2 * len + 0.25 );

  UINT4 dim = round ( dim0 );
  XLAL_CHECK_NULL ( fabs (dim0 - 1.0 * dim) < 1e-4, XLAL_EDOM, "Vector length '%d' invalid: must be of the form (n^2 + n)/2, where n is a positive integer\n", len );

  gsl_matrix *gij = gsl_matrix_calloc ( dim, dim );
  XLAL_CHECK_NULL ( gij != NULL, XLAL_ENOMEM, "gsl_matrix_calloc(%d,%d) failed\n", dim, dim );

  for ( UINT4 i = 0; i < dim; i ++ )
    {
      for ( UINT4 j = i; j < dim; j ++ )
        {
          REAL8 el_ij = m0->data[ PMETRIC_INDEX ( i, j ) ];
          if ( i == 3 ) el_ij /= Freq;
          if ( j == 3 ) el_ij /= Freq;
          gsl_matrix_set ( gij, i, j, el_ij );
          gsl_matrix_set ( gij, j, i, el_ij );

        } // for j = i ... dim
    } // for i < dim

  return gij;

} // convert_old_metric_2_new()


// flexible comparison function for 2 metrics
// returns maximal relative deviation, measured in terms of diagnal entries,
// ie  dg_ij = (g1_ij - g2_ij) / sqrt( g_ii * g_jj )
// this should always be well-defined, as we deal with positive-definite square matrices
static REAL8
compare_metrics ( const gsl_matrix *g1_ij, const gsl_matrix *g2_ij )
{
  XLAL_CHECK_REAL8 ( g1_ij, XLAL_EINVAL, "Invalid NULL input metric 'g1_ij'\n");
  XLAL_CHECK_REAL8 ( g2_ij, XLAL_EINVAL, "Invalid NULL input metric 'g2_ij'\n");

  UINT4 dim = g1_ij->size1;
  XLAL_CHECK_REAL8 ( g1_ij->size1 == g1_ij->size2, XLAL_EDOM, "Input matrix 'g1_ij' is not square, got %d x %d\n", g1_ij->size1, g1_ij->size2 );
  XLAL_CHECK_REAL8 ( g2_ij->size1 == g2_ij->size2, XLAL_EDOM, "Input matrix 'g2_ij' is not square, got %d x %d\n", g2_ij->size1, g2_ij->size2 );
  XLAL_CHECK_REAL8 ( g1_ij->size1 == g2_ij->size1, XLAL_EDOM, "Input metrics have different sizes: g1_ij = %d-dim, g2_ij = %d-dim\n", g1_ij->size1, g2_ij->size1 );

  //XLALfprintfGSLmatrix ( stderr, "%.15e", g1_ij );
  //XLALfprintfGSLmatrix ( stderr, "%.15e", g2_ij );

  REAL8 errmax = 0;
  for ( UINT4 i = 0; i < dim; i ++ )
    {
      for ( UINT4 j = 0; j < dim; j ++ )
        {
          REAL8 norm = sqrt ( gsl_matrix_get ( g1_ij, i, i ) * gsl_matrix_get ( g1_ij, j, j ) );
          REAL8 e1 = gsl_matrix_get ( g1_ij, i, j ) / norm;
          REAL8 e2 = gsl_matrix_get ( g2_ij, i, j ) / norm;
          REAL8 base;
          if ( e2 == 0 )
            base = 1;
          else
            base = e2;
          REAL8 reldiff = fabs ( e1 - e2 ) / base;

          errmax = fmax ( errmax, reldiff );
        } // for j < dim
    } // for i < dim

  return errmax;

} /* compare_metrics() */


/** Unit test function for XLALComputeOrbitalDerivatives()
 */
static int
test_XLALComputeOrbitalDerivatives ( void )
{
  // ----- load an example ephemeris, describing a pure cicular 2D
  // orbit w period of one year
  CHAR earthEphem[] = TEST_DATA_DIR "circularEphem.dat";
  CHAR sunEphem[]   = TEST_DATA_DIR "sun00-19-DE405.dat";

  EphemerisData *edat = XLALInitBarycenter ( earthEphem, sunEphem );
  XLAL_CHECK ( edat != NULL, XLAL_EFUNC, "XLALInitBarycenter('%s','%s') failed with xlalErrno = %d\n", earthEphem, sunEphem, xlalErrno );

  /* Unit test for XLALComputeOrbitalDerivatives() */
  // ----------------------------------------
  // the tests consists in using a purely 2D circular orbit to
  // compute the derivatives, so we can analytically check them!

  // one line from the middle of cicularEphem.dat:
  // 700244800    498.41220200278  24.31154155846971  0    -4.84039532365306e-06  9.923320107134072e-05  0    -1.975719726623028e-11 -9.63716218195963e-13 0

  // we can use this to check derivatives '0-'2'
  vect3Dlist_t *rOrb_n = NULL;
  LIGOTimeGPS t0 = { 700244800, 0 };
  vect3D_t r_0 = {498.41220200278, 24.31154155846971,  0 };
  vect3D_t r_1 = {-4.84039532365306e-06,   9.923320107134072e-05,  0 };
  vect3D_t r_2 = { -1.975719726623028e-11, -9.63716218195963e-13,  0 };
  UINT4 maxorder = 4;

  rOrb_n = XLALComputeOrbitalDerivatives ( maxorder, &t0, edat );
  XLAL_CHECK ( rOrb_n != NULL, XLAL_EFUNC, "XLALComputeOrbitalDerivatives() failed with xlalErrno = %d!.\n", xlalErrno );

  // ----- first check differences to known first derivatives 0 - 2:
  vect3D_t diff;
  REAL8 reldiff;
  REAL8 reltol = 1e-6;
  UINT4 i;

  /* order = 0 */
  for (i=0; i<3; i++) diff[i] = r_0[i] - rOrb_n->data[0][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  XLAL_CHECK ( reldiff <= reltol, XLAL_ETOL, "Relative error %g on 0th order r_0 exceeds tolerance of %g.\n", reldiff, reltol );

  /* order = 1 */
  for (i=0; i<3; i++) diff[i] = r_1[i] - rOrb_n->data[1][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  XLAL_CHECK ( reldiff <= reltol, XLAL_ETOL, "Relative error %g on 1st order r_1 exceeds tolerance of %g.\n", reldiff, reltol );

  /* order = 1 */
  for (i=0; i<3; i++) diff[i] = r_2[i] - rOrb_n->data[2][i];
  reldiff = sqrt ( NORM(diff) / NORM(r_0) );
  XLAL_CHECK ( reldiff <= reltol, XLAL_ETOL, "Relative error %g on 2n order r_2 exceeds tolerance of %g.\n", reldiff, reltol );

  // ----- second check derivatives against known analytic results
  REAL8 RorbC = LAL_AU_SI / LAL_C_SI;
  REAL8 Om = LAL_TWOPI / LAL_YRSID_SI;
  vect3D_t nR, nV;
  REAL8 norm;
  vect3D_t rTest_n[maxorder+1];

  norm = NORM( r_0 );
  DIV_VECT (nR, r_0, norm );
  norm = NORM( r_1 );
  DIV_VECT (nV, r_1, norm );

  /* r_0 */
  COPY_VECT(rTest_n[0], nR );
  MULT_VECT(rTest_n[0], RorbC );
  /* r_1 */
  COPY_VECT(rTest_n[1], nV );
  MULT_VECT(rTest_n[1], RorbC*Om );
  /* r_2 */
  COPY_VECT(rTest_n[2], nR );
  MULT_VECT(rTest_n[2], -RorbC*Om*Om );
  /* r_3 */
  COPY_VECT(rTest_n[3], nV );
  MULT_VECT(rTest_n[3], -RorbC*Om*Om*Om );
  /* r_4 */
  COPY_VECT(rTest_n[4], nR );
  MULT_VECT(rTest_n[4], RorbC*Om*Om*Om*Om );


  UINT4 n;
  reltol = 1e-2;	/* 1% error should be enough to enforce order 0-2 are ~1e-5 - 1e-8, order 3-4 are ~ 5e-3*/
  for ( n=0; n <= maxorder; n ++ )
    {
      for (i=0; i<3; i++) {
        diff[i] = rTest_n[n][i] - rOrb_n->data[n][i];
      }
      reldiff = sqrt ( NORM(diff) / NORM(rTest_n[n]) );
      XLALPrintInfo ("order %d: relative difference = %g\n", n, reldiff );
      XLAL_CHECK ( reldiff <= reltol, XLAL_ETOL,
                   "Relative error %g on r_%d exceeds tolerance of %g.\n"
                   "rd[%d] = {%g, %g, %g},  rTest[%d] = {%g, %g, %g}\n",
                   reldiff, n, reltol,
                   n, rOrb_n->data[n][0], rOrb_n->data[n][1], rOrb_n->data[n][2],
                   n, rTest_n[n][0], rTest_n[n][1], rTest_n[n][2] );
    } /* for n <= maxorder */

  /* free memory */
  XLALDestroyVect3Dlist ( rOrb_n );
  XLALDestroyEphemerisData ( edat );

  return XLAL_SUCCESS;

} /* test_XLALComputeOrbitalDerivatives() */
