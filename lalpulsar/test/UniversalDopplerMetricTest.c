/*
 * Copyright (C) 2015 Karl Wette
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
#include <lal/MetricUtils.h>

#include <lal/GSLHelpers.h>

/**
 * \author Reinhard Prix
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

// ---------- local types --------------------

typedef enum {
  OLDMETRIC_TYPE_PHASE = 0,        /**< compute phase metric only */
  OLDMETRIC_TYPE_FSTAT = 1,        /**< compute full F-metric only */
  OLDMETRIC_TYPE_ALL   = 2,        /**< compute both F-metric and phase-metric */
  OLDMETRIC_TYPE_LAST
} OldMetricType_t;

typedef struct tagOldDopplerMetric
{
  DopplerMetricParams meta;             /**< "meta-info" describing/specifying the type of Doppler metric */

  gsl_matrix *g_ij;                     /**< symmetric matrix holding the phase-metric, averaged over segments */
  gsl_matrix *g_ij_seg;                 /**< the phase-metric for each segment, concatenated by column: [g_ij_1, g_ij_2, ...] */

  gsl_matrix *gF_ij;                    /**< full F-statistic metric gF_ij, including antenna-pattern effects (see \cite Prix07) */
  gsl_matrix *gFav_ij;                  /**< 'average' Fstat-metric */
  gsl_matrix *m1_ij, *m2_ij, *m3_ij;    /**< Fstat-metric sub components */

  gsl_matrix *Fisher_ab;                /**< Full 4+n dimensional Fisher matrix, ie amplitude + Doppler space */

  double maxrelerr_gPh;                 /**< estimate for largest relative error in phase-metric component integrations */
  double maxrelerr_gF;                  /**< estimate for largest relative error in Fmetric component integrations */

  REAL8 rho2;                           /**< signal SNR rho^2 = A^mu M_mu_nu A^nu */
} OldDopplerMetric;

// ---------- local prototypes
static int test_XLALComputeOrbitalDerivatives ( void );
static int test_XLALComputeDopplerMetrics ( void );

static gsl_matrix *convert_old_metric_2_new ( const REAL8Vector *m0, REAL8 Freq );

OldDopplerMetric *XLALOldDopplerFstatMetric ( const OldMetricType_t metricType, const DopplerMetricParams *metricParams, const EphemerisData *edat );
void XLALDestroyOldDopplerMetric ( OldDopplerMetric *metric );
int XLALAddOldDopplerMetric ( OldDopplerMetric **metric1, const OldDopplerMetric *metric2 );
int XLALScaleOldDopplerMetric ( OldDopplerMetric *m, REAL8 scale );

// ---------- function definitions --------------------
/**
 * MAIN function: calls a number of unit-tests
 */
int main( void )
{
  INT4 ret;

  ret = test_XLALComputeDopplerMetrics();
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "test_XLALComputeDopplerMetrics() failed with status = %d.\n", ret );

  ret = test_XLALComputeOrbitalDerivatives();
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "test_XLALComputeOrbitalDerivatives() failed with status = %d.\n", ret );

  /* check for memory leaks */
  LALCheckMemoryLeaks();

  /* all tests passed */
  return XLAL_SUCCESS;

} /* main() */


/**
 * Unit test for metric functions XLALComputeDopplerPhaseMetric()
 * and XLALComputeDopplerFstatMetric()
 *
 * Initially modelled afer testMetricCodes.py script:
 * Check metric codes 'getMetric' 'FstatMetric' and 'FstatMetric_v2' by
 * comparing them against each other.
 * Given that they represent 3 very different implementations of
 * metric calculations, this provides a very powerful consistency test
 *
 */
static int
test_XLALComputeDopplerMetrics ( void )
{
  int ret;
  const REAL8 tolPh = 0.01;	// 1% tolerance on phase metrics [taken from testMetricCodes.py]

  // ----- load ephemeris
  const char earthEphem[] = TEST_DATA_DIR "earth00-19-DE200.dat.gz";
  const char sunEphem[]   = TEST_DATA_DIR "sun00-19-DE200.dat.gz";
  EphemerisData *edat = XLALInitBarycenter ( earthEphem, sunEphem );
  XLAL_CHECK ( edat != NULL, XLAL_EFUNC, "XLALInitBarycenter('%s','%s') failed with xlalErrno = %d\n", earthEphem, sunEphem, xlalErrno );

  // ----- set test-parameters ----------
  const LIGOTimeGPS startTimeGPS = { 792576013, 0 };
  const REAL8 Tseg = 60000;

  const REAL8 Alpha = 1.0;
  const REAL8 Delta = 0.5;
  const REAL8 Freq  = 100;
  const REAL8 f1dot = 0;// -1e-8;

  LALStringVector *detNames = XLALCreateStringVector ( "H1", "L1", "V1",  NULL );
  LALStringVector *sqrtSX   = XLALCreateStringVector ( "1.0", "0.5", "1.5", NULL );

  MultiLALDetector multiIFO;
  XLAL_CHECK ( XLALParseMultiLALDetector ( &multiIFO, detNames ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyStringVector ( detNames );

  MultiNoiseFloor multiNoiseFloor;
  XLAL_CHECK ( XLALParseMultiNoiseFloor ( &multiNoiseFloor, sqrtSX, multiIFO.length ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyStringVector ( sqrtSX );

  // prepare metric parameters for modern XLALComputeDopplerFstatMetric() and mid-old XLALOldDopplerFstatMetric()
  const DopplerCoordinateSystem coordSys = { 4, { DOPPLERCOORD_FREQ, DOPPLERCOORD_ALPHA, DOPPLERCOORD_DELTA, DOPPLERCOORD_F1DOT } };
  const PulsarAmplitudeParams Amp = { 0.03, -0.3, 0.5, 0.0 };	// h0, cosi, psi, phi0
  const PulsarDopplerParams dop = {
    .refTime  = startTimeGPS,
    .Alpha    = Alpha,
    .Delta    = Delta,
    .fkdot    = { Freq, f1dot },
  };

  LALSegList XLAL_INIT_DECL(segList);
  ret = XLALSegListInitSimpleSegments ( &segList, startTimeGPS, 1, Tseg );
  XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALSegListInitSimpleSegments() failed with xlalErrno = %d\n", xlalErrno );

  const DopplerMetricParams master_pars2 = {
    .coordSys			= coordSys,
    .detMotionType		= DETMOTION_SPIN | DETMOTION_ORBIT,
    .segmentList		= segList,
    .multiIFO			= multiIFO,
    .multiNoiseFloor		= multiNoiseFloor,
    .signalParams		= { .Amp = Amp, .Doppler = dop },
    .projectCoord		= - 1,	// -1==no projection
    .approxPhase		= 0,
  };

  // ----- prepare call-parameters of ancient LALPulsarMetric()
  LALStatus XLAL_INIT_DECL(status);
  REAL4 master_pars0_spindown_data[1] = { f1dot / Freq }; /* 'old-style' "f1" = f1dot / Freq !!*/
  REAL4Vector master_pars0_spindown = { .length = 1, .data = master_pars0_spindown_data };
  const PtoleMetricIn master_pars0 = {
    .spindown		= &master_pars0_spindown,
    .position		= { .longitude = Alpha, .latitude = Delta, .system = COORDINATESYSTEM_EQUATORIAL },
    .epoch		= startTimeGPS,
    .duration		= Tseg,
    .maxFreq		= Freq,
    .site		= &(multiIFO.sites[0]),	// use first detector for phase-metric
    .ephemeris		= edat,
    .metricType		= LAL_PMETRIC_COH_EPHEM,
  };


  // ========== BEGINNING OF TEST CALLS ==========


  XLALPrintWarning("\n---------- ROUND 1: ephemeris-based, single-IFO phase metrics ----------\n");
  {
    REAL8Vector *metric0 = 0;
    gsl_matrix *g0_ij;
    OldDopplerMetric *metric1;
    DopplerPhaseMetric *metric2P;
    REAL8 diff_2_0, diff_2_1, diff_1_0;

    PtoleMetricIn pars0 = master_pars0;
    DopplerMetricParams pars2 = master_pars2;

    pars2.multiIFO.length = 1;	// truncate to first detector
    pars2.multiNoiseFloor.length = 1;	// truncate to first detector

    // 0) compute metric using ancient LALPulsarMetric() function (used in lalapps_getMetric)
    LALPulsarMetric ( &status, &metric0, &pars0 );
    XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALPulsarMetric() failed with status=%d: '%s'\n", status.statusCode, status.statusDescription );
    XLAL_CHECK ( (g0_ij = convert_old_metric_2_new ( metric0, Freq )) != NULL, XLAL_EFUNC );
    XLALDestroyREAL8Vector ( metric0 ); metric0 = NULL;
    // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
    XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( OLDMETRIC_TYPE_PHASE, &pars2, edat )) != NULL, XLAL_EFUNC );
    // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
    XLAL_CHECK ( (metric2P = XLALComputeDopplerPhaseMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

    // compare all 3 metrics against each other:
    XLAL_CHECK ( (diff_2_0 = XLALCompareMetrics ( metric2P->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g2,g0)= %g exceeds tolerance of %g\n", diff_2_0, tolPh );
    XLALPrintWarning ("diff_2_0 = %e\n", diff_2_0 );
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( metric2P->g_ij, metric1->g_ij )) < tolPh, XLAL_ETOL, "Error(g2,g1)= %g exceeds tolerance of %g\n", diff_2_1, tolPh );
    XLALPrintWarning ("diff_2_1 = %e\n", diff_2_1 );
    XLAL_CHECK ( (diff_1_0 = XLALCompareMetrics ( metric1->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g1,g0)= %g exceeds tolerance of %g\n", diff_1_0, tolPh );
    XLALPrintWarning ("diff_1_0 = %e\n", diff_1_0 );

    gsl_matrix_free ( g0_ij );
    XLALDestroyOldDopplerMetric ( metric1 );
    XLALDestroyDopplerPhaseMetric ( metric2P );
  }


  XLALPrintWarning("\n---------- ROUND 2: Ptolemaic-based, single-IFO phase metrics ----------\n");
  {
    REAL8Vector *metric0 = 0;
    gsl_matrix *g0_ij;
    OldDopplerMetric *metric1;
    DopplerPhaseMetric *metric2P;
    REAL8 diff_2_0, diff_2_1, diff_1_0;

    PtoleMetricIn pars0 = master_pars0;
    DopplerMetricParams pars2 = master_pars2;

    pars2.multiIFO.length = 1;	// truncate to first detector
    pars2.multiNoiseFloor.length = 1;	// truncate to first detector

    pars0.metricType    = LAL_PMETRIC_COH_PTOLE_ANALYTIC;
    pars2.detMotionType = DETMOTION_SPIN | DETMOTION_PTOLEORBIT;

    // 0) compute metric using ancient LALPulsarMetric() function (used in lalapps_getMetric)
    LALPulsarMetric ( &status, &metric0, &pars0 );
    XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALPulsarMetric() failed with status=%d: '%s'\n", status.statusCode, status.statusDescription );
    XLAL_CHECK ( (g0_ij = convert_old_metric_2_new ( metric0, Freq )) != NULL, XLAL_EFUNC );
    XLALDestroyREAL8Vector ( metric0 ); metric0 = NULL;
    // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
    XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( OLDMETRIC_TYPE_PHASE, &pars2, edat )) != NULL, XLAL_EFUNC );
    // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
    XLAL_CHECK ( (metric2P = XLALComputeDopplerPhaseMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

    // compare all 3 metrics against each other:
    XLAL_CHECK ( (diff_2_0 = XLALCompareMetrics ( metric2P->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g2,g0)= %g exceeds tolerance of %g\n", diff_2_0, tolPh );
    XLALPrintWarning ("diff_2_0 = %e\n", diff_2_0 );
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( metric2P->g_ij, metric1->g_ij )) < tolPh, XLAL_ETOL, "Error(g2,g1)= %g exceeds tolerance of %g\n", diff_2_1, tolPh );
    XLALPrintWarning ("diff_2_1 = %e\n", diff_2_1 );
    XLAL_CHECK ( (diff_1_0 = XLALCompareMetrics ( metric1->g_ij, g0_ij )) < tolPh, XLAL_ETOL, "Error(g1,g0)= %g exceeds tolerance of %g\n", diff_1_0, tolPh );
    XLALPrintWarning ("diff_1_0 = %e\n", diff_1_0 );

    gsl_matrix_free ( g0_ij );
    XLALDestroyOldDopplerMetric ( metric1 );
    XLALDestroyDopplerPhaseMetric ( metric2P );
  }


  XLALPrintWarning("\n---------- ROUND 3: ephemeris-based, multi-IFO F-stat metrics ----------\n");
  {
    OldDopplerMetric *metric1;
    DopplerFstatMetric *metric2F;
    REAL8 diff_2_1;

    DopplerMetricParams pars2 = master_pars2;

    pars2.detMotionType = DETMOTION_SPIN | DETMOTION_ORBIT;
    pars2.multiIFO      = multiIFO;	// 3 IFOs
    pars2.multiNoiseFloor = multiNoiseFloor;// 3 IFOs

    // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
    XLAL_CHECK ( (metric1 = XLALOldDopplerFstatMetric ( OLDMETRIC_TYPE_FSTAT, &pars2, edat )) != NULL, XLAL_EFUNC );
    // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
    XLAL_CHECK ( (metric2F = XLALComputeDopplerFstatMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

    // compare both metrics against each other:
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( metric2F->gF_ij,   metric1->gF_ij ))   < tolPh, XLAL_ETOL, "Error(gF2,gF1)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
    XLALPrintWarning ("gF:   diff_2_1 = %e\n", diff_2_1 );
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( metric2F->gFav_ij, metric1->gFav_ij )) < tolPh, XLAL_ETOL, "Error(gFav2,gFav1)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
    XLALPrintWarning ("gFav: diff_2_1 = %e\n", diff_2_1 );

    XLALDestroyOldDopplerMetric ( metric1 );
    XLALDestroyDopplerFstatMetric ( metric2F );
  }


  XLALPrintWarning("\n---------- ROUND 4: compare analytic {f,f1dot,f2dot,f3dot} phase metric vs XLALComputeDopplerPhaseMetric() ----------\n");
  {
    DopplerPhaseMetric *metric2P;
    REAL8 diff_2_1;

    DopplerMetricParams pars2 = master_pars2;

    pars2.multiIFO.length  = 1;	// truncate to 1st detector
    pars2.multiNoiseFloor.length  = 1;	// truncate to 1st detector
    pars2.detMotionType   = DETMOTION_SPIN | DETMOTION_ORBIT;
    pars2.approxPhase     = 1;	// use same phase-approximation as in analytic solution to improve comparison

    DopplerCoordinateSystem coordSys2 = { 4, { DOPPLERCOORD_FREQ, DOPPLERCOORD_F1DOT, DOPPLERCOORD_F2DOT, DOPPLERCOORD_F3DOT } };
    pars2.coordSys = coordSys2;
    gsl_matrix* gN_ij;

    // a) compute metric at refTime = startTime
    pars2.signalParams.Doppler.refTime = startTimeGPS;
    XLAL_CHECK ( (metric2P = XLALComputeDopplerPhaseMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
    gN_ij = NULL;
    XLAL_CHECK ( XLALNaturalizeMetric ( &gN_ij, NULL, metric2P->g_ij, &pars2 ) == XLAL_SUCCESS, XLAL_EFUNC );

    REAL8 gStart_ij[] = {   1.0/3,      2.0/3,    6.0/5,    32.0/15,      \
                            2.0/3,    64.0/45,    8.0/3,  512.0/105,      \
                            6.0/5,      8.0/3,   36.0/7,     48.0/5,      \
                            32.0/15,  512.0/105, 48.0/5, 4096.0/225 };
    const gsl_matrix_view gStart = gsl_matrix_view_array ( gStart_ij, 4, 4 );

    // compare natural-units metric against analytic solution
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( gN_ij,   &(gStart.matrix) )) < tolPh, XLAL_ETOL,
                 "RefTime=StartTime: Error(g_ij,g_analytic)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
    XLALPrintWarning ("Analytic (refTime=startTime): diff_2_1 = %e\n", diff_2_1 );

    XLALDestroyDopplerPhaseMetric ( metric2P );
    gsl_matrix_free ( gN_ij );

    // b) compute metric at refTime = midTime
    pars2.signalParams.Doppler.refTime = startTimeGPS;
    pars2.signalParams.Doppler.refTime.gpsSeconds += Tseg / 2;

    XLAL_CHECK ( (metric2P = XLALComputeDopplerPhaseMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );
    gN_ij = NULL;
    XLAL_CHECK ( XLALNaturalizeMetric ( &gN_ij, NULL, metric2P->g_ij, &pars2 ) == XLAL_SUCCESS, XLAL_EFUNC );

    REAL8 gMid_ij[] = { 1.0/3,    0,        1.0/5,         0,       \
                        0,        4.0/45,       0,   8.0/105,       \
                        1.0/5,    0,        1.0/7,         0,       \
                        0,        8.0/105,      0,  16.0/225  };
    const gsl_matrix_view gMid = gsl_matrix_view_array ( gMid_ij, 4, 4 );

    // compare natural-units metric against analytic solution
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( gN_ij,   &(gMid.matrix) )) < tolPh, XLAL_ETOL,
                 "RefTime=MidTime: Error(g_ij,g_analytic)= %e exceeds tolerance of %e\n", diff_2_1, tolPh );
    XLALPrintWarning ("Analytic (refTime=midTime):   diff_2_1 = %e\n\n", diff_2_1 );

    XLALDestroyDopplerPhaseMetric ( metric2P );
    gsl_matrix_free ( gN_ij );
  }


  XLALPrintWarning("\n---------- ROUND 5: ephemeris-based, single-IFO, segment-averaged phase metrics ----------\n");
  {
    OldDopplerMetric *metric1;
    DopplerPhaseMetric *metric2P;
    REAL8 diff_2_1;

    DopplerMetricParams pars2 = master_pars2;

    pars2.detMotionType = DETMOTION_SPIN | DETMOTION_ORBIT;
    pars2.multiIFO.length = 1;	// truncate to first detector
    pars2.multiNoiseFloor.length = 1;	// truncate to first detector
    pars2.approxPhase = 1;

    const UINT4 Nseg = 10;
    LALSegList XLAL_INIT_DECL(NsegList);
    ret = XLALSegListInitSimpleSegments ( &NsegList, startTimeGPS, Nseg, Tseg );
    XLAL_CHECK ( ret == XLAL_SUCCESS, XLAL_EFUNC, "XLALSegListInitSimpleSegments() failed with xlalErrno = %d\n", xlalErrno );
    pars2.segmentList = NsegList;

    LALSegList XLAL_INIT_DECL(segList_k);
    LALSeg segment_k;
    XLALSegListInit( &segList_k );	// prepare single-segment list containing segment k
    segList_k.arraySize = 1;
    segList_k.length = 1;
    segList_k.segs = &segment_k;

    // 1) compute metric using old FstatMetric code, now wrapped into XLALOldDopplerFstatMetric()
    metric1 = NULL;
    for (UINT4 k = 0; k < Nseg; ++k) {
      // setup 1-segment segment-list pointing k-th segment
      DopplerMetricParams pars2_k = pars2;
      pars2_k.segmentList = segList_k;
      pars2_k.segmentList.segs[0] = pars2.segmentList.segs[k];
      // XLALOldDopplerFstatMetric() does not agree numerically with UniversalDopplerMetric when using refTime != startTime
      pars2_k.signalParams.Doppler.refTime = pars2_k.segmentList.segs[0].start;

      OldDopplerMetric *metric1_k;   // per-segment coherent metric
      XLAL_CHECK ( (metric1_k = XLALOldDopplerFstatMetric ( OLDMETRIC_TYPE_PHASE, &pars2_k, edat )) != NULL, XLAL_EFUNC );

      // manually correct reference time of metric1_k->g_ij; see Prix, "Frequency metric for CW searches" (2014-08-17), p. 4
      const double dt = XLALGPSDiff( &(pars2_k.signalParams.Doppler.refTime), &(pars2.signalParams.Doppler.refTime) );
      const double gFF = gsl_matrix_get( metric1_k->g_ij, 0, 0 );
      const double gFA = gsl_matrix_get( metric1_k->g_ij, 0, 1 );
      const double gFD = gsl_matrix_get( metric1_k->g_ij, 0, 2 );
      const double gFf = gsl_matrix_get( metric1_k->g_ij, 0, 3 );
      const double gAf = gsl_matrix_get( metric1_k->g_ij, 1, 3 );
      const double gDf = gsl_matrix_get( metric1_k->g_ij, 2, 3 );
      const double gff = gsl_matrix_get( metric1_k->g_ij, 3, 3 );
      gsl_matrix_set( metric1_k->g_ij, 0, 3, gFf + gFF*dt ); gsl_matrix_set( metric1_k->g_ij, 3, 0, gsl_matrix_get( metric1_k->g_ij, 0, 3 ) );
      gsl_matrix_set( metric1_k->g_ij, 1, 3, gAf + gFA*dt ); gsl_matrix_set( metric1_k->g_ij, 3, 1, gsl_matrix_get( metric1_k->g_ij, 1, 3 ) );
      gsl_matrix_set( metric1_k->g_ij, 2, 3, gDf + gFD*dt ); gsl_matrix_set( metric1_k->g_ij, 3, 2, gsl_matrix_get( metric1_k->g_ij, 2, 3 ) );
      gsl_matrix_set( metric1_k->g_ij, 3, 3, gff + 2*gFf*dt + gFF*dt*dt );

      XLAL_CHECK ( XLALAddOldDopplerMetric ( &metric1, metric1_k ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyOldDopplerMetric ( metric1_k );
    }
    XLAL_CHECK ( XLALScaleOldDopplerMetric ( metric1, 1.0 / Nseg ) == XLAL_SUCCESS, XLAL_EFUNC );

    // 2) compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
    XLAL_CHECK ( (metric2P = XLALComputeDopplerPhaseMetric ( &pars2, edat )) != NULL, XLAL_EFUNC );

    GPMAT( metric1->g_ij, "%0.4e" );
    GPMAT( metric2P->g_ij, "%0.4e" );

    // compare both metrics against each other:
    XLAL_CHECK ( (diff_2_1 = XLALCompareMetrics ( metric2P->g_ij, metric1->g_ij )) < tolPh, XLAL_ETOL, "Error(g2,g1)= %g exceeds tolerance of %g\n", diff_2_1, tolPh );
    XLALPrintWarning ("diff_2_1 = %e\n", diff_2_1 );

    XLALDestroyOldDopplerMetric ( metric1 );
    XLALDestroyDopplerPhaseMetric ( metric2P );

    XLALSegListClear ( &NsegList );
  }


  XLALPrintWarning("\n---------- ROUND 6: directed binary orbital metric ----------\n");
  {
    REAL8 Period = 68023.70496;
    REAL8 Omega = LAL_TWOPI / Period;
    REAL8 asini = 1.44;
    REAL8 tAsc = 897753994;
    REAL8 argp = 0;
    LIGOTimeGPS tP; XLALGPSSetREAL8 ( &tP, tAsc + argp / Omega );

    const PulsarDopplerParams dopScoX1 = {
      .refTime  = startTimeGPS,
      .Alpha    = Alpha,
      .Delta    = Delta,
      .fkdot    = { Freq },
      .asini    = asini,
      .period   = Period,
      .tp       = tP
    };
    REAL8 TspanScoX1 = 20 * 19 * 3600;	// 20xPorb for long-segment regime
    LALSegList XLAL_INIT_DECL(segListScoX1);
    XLAL_CHECK ( XLALSegListInitSimpleSegments ( &segListScoX1, startTimeGPS, 1, TspanScoX1 ) == XLAL_SUCCESS, XLAL_EFUNC );
    REAL8 tMid = XLALGPSGetREAL8(&startTimeGPS) + 0.5 * TspanScoX1;
    REAL8 DeltaMidAsc = tMid - tAsc;
    const DopplerCoordinateSystem coordSysScoX1 = { 6, { DOPPLERCOORD_FREQ, DOPPLERCOORD_ASINI, DOPPLERCOORD_TASC, DOPPLERCOORD_PORB, DOPPLERCOORD_KAPPA, DOPPLERCOORD_ETA } };
    DopplerMetricParams pars_ScoX1 = {
      .coordSys			= coordSysScoX1,
      .detMotionType		= DETMOTION_SPIN | DETMOTION_ORBIT,
      .segmentList		= segListScoX1,
      .multiIFO			= multiIFO,
      .multiNoiseFloor		= multiNoiseFloor,
      .signalParams		= { .Amp = Amp, .Doppler = dopScoX1 },
      .projectCoord		= - 1,	// -1==no projection
      .approxPhase		= 1,
    };
    pars_ScoX1.multiIFO.length = 1;	// truncate to first detector
    pars_ScoX1.multiNoiseFloor.length = 1;	// truncate to first detector

    // compute metric using modern UniversalDopplerMetric module: (used in lalapps_FstatMetric_v2)
    DopplerPhaseMetric *metric_ScoX1;
    XLAL_CHECK ( (metric_ScoX1 = XLALComputeDopplerPhaseMetric ( &pars_ScoX1, edat )) != NULL, XLAL_EFUNC );

    // compute analytic metric computed from Eq.(47) in Leaci,Prix PRD91, 102003 (2015):
    gsl_matrix *g0_ij;
    XLAL_CHECK ( (g0_ij = gsl_matrix_calloc ( 6, 6 )) != NULL, XLAL_ENOMEM, "Failed to gsl_calloc a 6x6 matrix\n");
    gsl_matrix_set ( g0_ij, 0, 0, pow ( LAL_PI * TspanScoX1, 2 ) / 3.0 );
    gsl_matrix_set ( g0_ij, 1, 1, 2.0 * pow ( LAL_PI * Freq, 2 ) );
    gsl_matrix_set ( g0_ij, 2, 2, 2.0 * pow ( LAL_PI * Freq * asini * Omega, 2 ) );
    gsl_matrix_set ( g0_ij, 3, 3, 0.5 * pow ( Omega, 4 ) * pow ( Freq * asini, 2 ) * ( pow ( TspanScoX1, 2 ) / 12.0 + pow ( DeltaMidAsc, 2 ) ) );
    REAL8 gPAsc = LAL_PI * pow ( Freq * asini, 2 ) * pow ( Omega, 3 ) * DeltaMidAsc;
    gsl_matrix_set ( g0_ij, 2, 3, gPAsc );
    gsl_matrix_set ( g0_ij, 3, 2, gPAsc );
    gsl_matrix_set ( g0_ij, 4, 4, 0.5 * pow ( LAL_PI * Freq * asini, 2 ) );
    gsl_matrix_set ( g0_ij, 5, 5, 0.5 * pow ( LAL_PI * Freq * asini, 2 ) );

    GPMAT ( metric_ScoX1->g_ij, "%0.4e" );
    GPMAT ( g0_ij, "%0.4e" );

    // compare metrics against each other
    REAL8 diff, tolScoX1 = 0.05;
    XLAL_CHECK ( (diff = XLALCompareMetrics ( metric_ScoX1->g_ij, g0_ij )) < tolScoX1, XLAL_ETOL, "Error(gNum,gAn)= %g exceeds tolerance of %g\n", diff, tolScoX1 );
    XLALPrintWarning ("diff_Num_An = %e\n", diff );

    gsl_matrix_free ( g0_ij );
    XLALDestroyDopplerPhaseMetric ( metric_ScoX1 );
    XLALSegListClear ( &segListScoX1 );
  }


  // ----- clean up memory
  XLALSegListClear ( &segList );
  XLALDestroyEphemerisData ( edat );

  return XLAL_SUCCESS;

} /* test_XLALComputeDopplerMetrics() */

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


/**
 * Unit test function for XLALComputeOrbitalDerivatives()
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
  XLALPrintInfo ("\nComparing Earth orbital derivatives:\n");
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
