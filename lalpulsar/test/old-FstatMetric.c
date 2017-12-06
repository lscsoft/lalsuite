/*
 * Copyright (C) 2006, 2013 Reinhard Prix
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
 * \file
 * \ingroup UniversalDopplerMetric_h
 * \author Reinhard Prix
 * \brief
 * The only purpose of this file is to serve as a backwards-comparison
 * check for XLALDopplerFstatMetric(). This used to be a standalone-code
 * 'lalapps_FstatMetric', and was XLALified and moved into the test-directory,
 * main() was wrapped into the forwards-compatible function XLALOldDopplerFstatMetric()
 * and called in UniversalDopplerMetricTest for comparison.
 */

/* ---------- includes ---------- */
#include <math.h>


#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


#include <lal/UserInput.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/LALInitBarycenter.h>
#include <lal/AVFactories.h>
#include <lal/SkyCoordinates.h>
#include <lal/ComputeFstat.h>
#include <lal/PulsarTimes.h>
#include <lal/SFTutils.h>

#include <lal/FlatPulsarMetric.h>
#include <lal/ComputeFstat.h>
#include <lal/UniversalDopplerMetric.h>

/* ----- compile switches ----- */
/* uncomment the following to turn off range-checking in GSL vector-functions */
/* #define GSL_RANGE_CHECK_OFF 1*/

/*---------- local defines ---------- */

#define NUM_SPINS 	2
#define METRIC_DIM 	2 + NUM_SPINS

/* ----- Macros ----- */
/** Simple Euklidean scalar product for two 3-dim vectors in cartesian coords */
#define SCALAR(u,v) ((u)[0]*(v)[0] + (u)[1]*(v)[1] + (u)[2]*(v)[2])

/** copy 3 components of Euklidean vector */
#define COPY_VECT(dst,src) do { (dst)[0] = (src)[0]; (dst)[1] = (src)[1]; (dst)[2] = (src)[2]; } while(0)

#define SQ(x) ((x) * (x))

/* ---------- local types ---------- */
typedef enum {
  PHASE_NONE = -1,
  PHASE_FULL = 0,
  PHASE_ORBITAL,
  PHASE_SPIN,
  PHASE_PTOLE,
  PHASE_LAST
} PhaseType_t;

/** a 'point' in the "Doppler parameter space" {alpha, delta, fkdot } */
typedef struct
{
  SkyPosition skypos;
  REAL8Vector *fkdot;
} DopplerPoint;

typedef struct
{
  REAL8Vector *dAlpha;
  REAL8Vector *dDelta;
  REAL8Vector *dFreq;
  REAL8Vector *df1dot;
} PhaseDerivs;

typedef struct
{
  UINT4 length;		/**< number of IFOs */
  PhaseDerivs **data;	/**< phase-derivs array */
} MultiPhaseDerivs;

typedef struct {
  REAL8 pos[3];
  REAL8 vel[3];
} PosVel_t;

typedef struct
{
  const EphemerisData *edat;		/**< ephemeris data (from LALInitBarycenter()) */
  LIGOTimeGPS startTime;		/**< start time of observation */
  LIGOTimeGPS refTime;			/**< reference time for spin-parameters  */
  DopplerPoint dopplerPoint;		/**< sky-position and spins */
  MultiDetectorStateSeries *multiDetStates;/**< pos, vel and LMSTs for detector at times t_i */
  MultiAMCoeffs *multiAMcoe;         	/**< Amplitude Modulation coefficients a,b(t)*/
  MultiPhaseDerivs *multidPhi;		/**< Phase-derivatives d_i phi(t) */
  REAL8Vector *GLweights;		/**< Gauss-Legendre Integration-weights */
  MultiNoiseWeights *multiNoiseWeights;	/**< noise-weights for the different IFOs */
  REAL8 Al1, Al2, Al3;			/**< amplitude factors alpha1, alpha2, alpha3 */
  REAL8 Ad, Bd, Cd, Dd;
  PhaseType_t phaseType;
} ConfigVariables;

/* ---------- global variables ----------*/

/* ---------- local prototypes ---------- */
int InitCode ( ConfigVariables *cfg, const DopplerMetricParams *metricParams, const EphemerisData *edat );
MultiPhaseDerivs * getMultiPhaseDerivs ( const MultiDetectorStateSeries *detStates, const DopplerPoint *dopplerPoint, PhaseType_t type );

int computeFstatMetric ( gsl_matrix *gF_ij, gsl_matrix *gFav_ij,
			 gsl_matrix *m1_ij, gsl_matrix *m2_ij, gsl_matrix *m3_ij,
			 ConfigVariables *cfg );

int computePhaseMetric ( gsl_matrix *g_ij, const PhaseDerivs *dphi, const REAL8Vector *GLweights );

int project_metric( gsl_matrix *ret_ij, gsl_matrix *g_ij, const UINT4 coordinate );
int outer_product ( gsl_matrix *ret_ij, const gsl_vector *u_i, const gsl_vector *v_j );
int symmetrize ( gsl_matrix *mat );
REAL8 quad_form ( const gsl_matrix *mat, const gsl_vector *vec );


void getPtolePosVel( PosVel_t *posvel, REAL8 tGPS, REAL8 tAutumn );

void XLALDestroyMultiPhaseDerivs ( MultiPhaseDerivs *mdPhi );

void gauleg(double x1, double x2, double x[], double w[], int n);

DopplerMetric *XLALOldDopplerFstatMetric ( const DopplerMetricParams *metricParams, const EphemerisData *edat );

/*============================================================
 * FUNCTION definitions
 *============================================================*/

/**
 * The only purpose of this function is to serve as a backwards-comparison
 * check for XLALDopplerFstatMetric(). This is why it has been moved
 * into the test-directory
 *
 * This is basically a wrapper of the 'main()' function from the old
 * standalone 'lalapps_FstatMetric' code, providing an API compatible
 * with XLALDopplerFstatMetric().
 *
 */
DopplerMetric *
XLALOldDopplerFstatMetric ( const DopplerMetricParams *metricParams,  	/**< input parameters determining the metric calculation */
                            const EphemerisData *edat			/**< ephemeris data */
                            )
{
  XLAL_CHECK_NULL ( metricParams != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( edat != NULL, XLAL_EINVAL );

  XLAL_CHECK_NULL ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL );
  UINT4 Nseg = metricParams->segmentList.length;
  XLAL_CHECK_NULL ( Nseg == 1, XLAL_EINVAL, "Segment list must only contain Nseg=1 segments, got Nseg=%d", Nseg );

  MetricType_t metricType = metricParams->metricType;
  XLAL_CHECK_NULL ( metricType < METRIC_TYPE_LAST, XLAL_EDOM );

  LIGOTimeGPS *startTime = &(metricParams->segmentList.segs[0].start);
  LIGOTimeGPS *endTime   = &(metricParams->segmentList.segs[0].end);
  REAL8 duration = XLALGPSDiff( endTime, startTime );

  const DopplerCoordinateSystem *coordSys = &(metricParams->coordSys);
  XLAL_CHECK_NULL ( coordSys->dim == METRIC_DIM, XLAL_EINVAL );
  XLAL_CHECK_NULL ( coordSys->coordIDs[0] == DOPPLERCOORD_FREQ, XLAL_EDOM );
  XLAL_CHECK_NULL ( coordSys->coordIDs[1] == DOPPLERCOORD_ALPHA, XLAL_EDOM );
  XLAL_CHECK_NULL ( coordSys->coordIDs[2] == DOPPLERCOORD_DELTA, XLAL_EDOM );
  XLAL_CHECK_NULL ( coordSys->coordIDs[3] == DOPPLERCOORD_F1DOT, XLAL_EDOM );

  DopplerMetric *metric;
  XLAL_CHECK_NULL ( (metric = XLALCalloc ( 1, sizeof(*metric) )) != NULL, XLAL_ENOMEM );

  ConfigVariables XLAL_INIT_DECL(config);

  /* basic setup and initializations */
  metric->gF_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  metric->gFav_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  metric->m1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  metric->m2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  metric->m3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  metric->g_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  XLAL_CHECK_NULL ( InitCode( &config, metricParams, edat ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ----- compute phase derivatives ----- */
  config.multidPhi = getMultiPhaseDerivs ( config.multiDetStates, &(config.dopplerPoint), config.phaseType );
  XLAL_CHECK_NULL ( config.multidPhi != NULL, XLAL_EFUNC, "getMultiPhaseDerivs() failed.\n" );

  if ( (metricType == METRIC_TYPE_FSTAT) || (metricType == METRIC_TYPE_ALL) )
    {
      XLAL_CHECK_NULL ( computeFstatMetric ( metric->gF_ij, metric->gFav_ij, metric->m1_ij, metric->m2_ij, metric->m3_ij, &config ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  if ( (metricType == METRIC_TYPE_PHASE) || (metricType == METRIC_TYPE_ALL) )
    {
      if ( metricParams->detMotionType == (DETMOTION_SPINXY | DETMOTION_ORBIT) )
        {
          XLAL_CHECK_NULL ( XLALFlatMetricCW ( metric->g_ij, config.refTime, config.startTime, duration, edat ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      else
        {
          XLAL_CHECK_NULL ( computePhaseMetric ( metric->g_ij, config.multidPhi->data[0], config.GLweights) == XLAL_SUCCESS, XLAL_EFUNC );
        }
    } // endif metricType==PHASE || ALL

  // ----- Free internal memory
  XLALDestroyMultiPhaseDerivs ( config.multidPhi );
  XLALDestroyMultiAMCoeffs ( config.multiAMcoe );
  XLALDestroyMultiDetectorStateSeries ( config.multiDetStates );
  XLALDestroyREAL8Vector ( config.GLweights );
  for ( UINT4 X=0; X < config.multiNoiseWeights->length; X++ ) {
    XLALDestroyREAL8Vector ( config.multiNoiseWeights->data[X] );
  }
  XLALFree ( config.multiNoiseWeights->data );
  XLALFree ( config.multiNoiseWeights );
  XLALDestroyREAL8Vector ( config.dopplerPoint.fkdot );

  return metric;

} /* XLALOldDopplerFstatMetric() */


/* calculate Fstat-metric components m1_ij, m2_ij, m3_ij, and metrics gF_ij, gFav_ij,
 * which must be allocated 4x4 matrices.
 *
 * Return 0 = OK, -1 on error.
 */
int
computeFstatMetric ( gsl_matrix *gF_ij, gsl_matrix *gFav_ij,
		     gsl_matrix *m1_ij, gsl_matrix *m2_ij, gsl_matrix *m3_ij,
		     ConfigVariables *cfg )
{
  UINT4 i;
  REAL8 AMA;

  gsl_matrix *P1_ij, *P2_ij, *P3_ij;
  gsl_vector *a2_dPhi_i, *b2_dPhi_i, *ab_dPhi_i;

  gsl_matrix *Q1_ij, *Q2_ij, *Q3_ij;

  gsl_vector *dPhi_i;
  gsl_matrix *mat1, *mat2;
  gsl_vector *vec;

  gsl_matrix *a2_a2, *a2_b2, *a2_ab;
  gsl_matrix *b2_b2, *b2_ab;
  gsl_matrix *ab_ab;

  /* ----- check input ----- */
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL );
  UINT4 numDet = cfg->multiDetStates->length;

  XLAL_CHECK ( (gF_ij != NULL) && ( gFav_ij != NULL ) && (m1_ij != NULL) && (m2_ij != NULL) && (m3_ij != NULL), XLAL_EINVAL );

  XLAL_CHECK ( (gF_ij->size1 == gF_ij->size2) && (gF_ij->size1 == METRIC_DIM), XLAL_EINVAL );
  XLAL_CHECK ( (gFav_ij->size1 == gFav_ij->size2) && (gFav_ij->size1 == METRIC_DIM), XLAL_EINVAL );
  XLAL_CHECK ( (m1_ij->size1 == m1_ij->size2) && (m1_ij->size1 == METRIC_DIM), XLAL_EINVAL );
  XLAL_CHECK ( (m2_ij->size1 == m2_ij->size2) && (m2_ij->size1 == METRIC_DIM), XLAL_EINVAL );
  XLAL_CHECK ( (m3_ij->size1 == m3_ij->size2) && (m3_ij->size1 == METRIC_DIM), XLAL_EINVAL );

  XLAL_CHECK ( cfg->multiAMcoe->length == numDet, XLAL_EINVAL );
  XLAL_CHECK ( cfg->multidPhi->length == numDet, XLAL_EINVAL );

  /* ----- allocate matrices/vectors ----- */
  dPhi_i = gsl_vector_calloc ( METRIC_DIM );

  mat1 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  mat2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  vec = gsl_vector_calloc (  METRIC_DIM );

  P1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  P2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  P3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  a2_dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  b2_dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  ab_dPhi_i = gsl_vector_calloc ( METRIC_DIM );

  Q1_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  Q2_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  Q3_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  a2_a2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  a2_b2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  a2_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  b2_b2 = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  b2_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  ab_ab = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  /* ----- calculate averages ----- */
  REAL8 A = 0, B = 0, C = 0;
  for ( UINT4 X=0; X < numDet; X ++ )
    {
      UINT4 numSteps = cfg->multiDetStates->data[X]->length;
      AMCoeffs *amcoe = cfg->multiAMcoe->data[X];
      PhaseDerivs *dPhi = cfg->multidPhi->data[X];
      gsl_matrix *P1_Xij, *P2_Xij, *P3_Xij;
      gsl_vector *a2_dPhi_Xi, *b2_dPhi_Xi, *ab_dPhi_Xi;

      P1_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
      P2_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
      P3_Xij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

      a2_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );
      b2_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );
      ab_dPhi_Xi = gsl_vector_calloc ( METRIC_DIM );

      for ( i = 0; i < numSteps; i ++ )
	{
	  REAL8 a = amcoe->a->data[i];
	  REAL8 b = amcoe->b->data[i];
	  REAL8 a2 = SQ(a);
	  REAL8 b2 = SQ(b);
	  REAL8 ab = a * b;

          A += a2;
          B += b2;
          C += ab;

	  gsl_vector_set ( dPhi_i, 0, dPhi->dFreq->data[i] );
	  gsl_vector_set ( dPhi_i, 1, dPhi->dAlpha->data[i]);
	  gsl_vector_set ( dPhi_i, 2, dPhi->dDelta->data[i]);
	  gsl_vector_set ( dPhi_i, 3, dPhi->df1dot->data[i]);

	  outer_product(mat1, dPhi_i, dPhi_i );

	  /* ----- P1_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, a2 );
	  gsl_matrix_add ( P1_Xij, mat2 );

	  /* ----- P2_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, b2 );
	  gsl_matrix_add ( P2_Xij, mat2 );

	  /* ----- P3_ij ----- */
	  gsl_matrix_memcpy (mat2, mat1);
	  gsl_matrix_scale ( mat2, ab );
	  gsl_matrix_add ( P3_Xij, mat2 );

	  /* ----- a2_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, a2 );
	  gsl_vector_add ( a2_dPhi_Xi, vec );

	  /* ----- b2_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, b2 );
	  gsl_vector_add ( b2_dPhi_Xi, vec );

	  /* ----- ab_dPhi_i ----- */
	  gsl_vector_memcpy ( vec, dPhi_i );
	  gsl_vector_scale ( vec, ab );
	  gsl_vector_add ( ab_dPhi_Xi, vec );

	} /* for i < numSteps */

      gsl_matrix_add ( P1_ij, P1_Xij );
      gsl_matrix_add ( P2_ij, P2_Xij );
      gsl_matrix_add ( P3_ij, P3_Xij );

      gsl_vector_add ( a2_dPhi_i, a2_dPhi_Xi );

      gsl_vector_add ( b2_dPhi_i, b2_dPhi_Xi );

      gsl_vector_add ( ab_dPhi_i, ab_dPhi_Xi );

      gsl_matrix_free ( P1_Xij );
      gsl_matrix_free ( P2_Xij );
      gsl_matrix_free ( P3_Xij );

      gsl_vector_free ( a2_dPhi_Xi );
      gsl_vector_free ( b2_dPhi_Xi );
      gsl_vector_free ( ab_dPhi_Xi );

    } /* for X < numDet */

  /* ---------- composite quantities ---------- */
  REAL8 D = A * B - C*C;
  cfg->Ad = A;
  cfg->Bd = B;
  cfg->Cd = C;
  cfg->Dd = D;

  outer_product (a2_a2, a2_dPhi_i, a2_dPhi_i );
  outer_product (a2_b2, a2_dPhi_i, b2_dPhi_i );
  outer_product (a2_ab, a2_dPhi_i, ab_dPhi_i );

  outer_product (b2_b2, b2_dPhi_i, b2_dPhi_i );
  outer_product (b2_ab, b2_dPhi_i, ab_dPhi_i );

  outer_product (ab_ab, ab_dPhi_i, ab_dPhi_i );

  /* ----- Q1_ij ----- */
  gsl_matrix_memcpy ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q1_ij, mat1 );	/*  = (A/D)<ab_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_a2 );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q1_ij, mat1 );	/*  + (B/D)<a2_dPhi_i><a2_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_ab );
  gsl_matrix_scale ( mat1, - 2.0 * C/D );
  gsl_matrix_add ( Q1_ij, mat1 );	/*  -2(C/D)<a2_dPhi_i><ab_dPhi_j> */

  symmetrize ( Q1_ij );		/* (i,j) */

  /* ----- Q2_ij ----- */
  gsl_matrix_memcpy ( mat1, b2_b2 );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q2_ij, mat1 );	/*  = (A/D)<b2_dPhi_i><b2_dPhi_j> */

  gsl_matrix_memcpy ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q2_ij, mat1 );	/*  + (B/D)<ab_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, b2_ab );
  gsl_matrix_scale ( mat1, - 2.0 * C/D );
  gsl_matrix_add ( Q2_ij, mat1 );	/*  -2(C/D)<b2_dPhi_i><ab_dPhi_j> */

  symmetrize ( Q2_ij );		/* (i,j) */

  /* ----- Q3_ij ----- */
  gsl_matrix_memcpy ( mat1, b2_ab );
  gsl_matrix_scale ( mat1, A/D );
  gsl_matrix_memcpy ( Q3_ij, mat1 );	/*  = (A/D)<b2_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_ab );
  gsl_matrix_scale ( mat1, B/D );
  gsl_matrix_add ( Q3_ij, mat1 );	/*  + (B/D)<a2_dPhi_i><ab_dPhi_j> */

  gsl_matrix_memcpy ( mat1, a2_b2 );
  gsl_matrix_add ( mat1, ab_ab );
  gsl_matrix_scale ( mat1, - 1.0 * C/D );
  gsl_matrix_add ( Q3_ij, mat1 );	/*  -(C/D)(<a2_dPhi_i><b2_dPhi_j> + <ab_dPhi_i><ab_dPhi_j>)*/

  symmetrize ( Q3_ij );		/* (i,j) */

  /* ===== final matrics: m1_ij, m2_ij, m3_ij ===== */
  gsl_matrix_memcpy ( m1_ij, P1_ij );
  gsl_matrix_sub ( m1_ij, Q1_ij );

  gsl_matrix_memcpy ( m2_ij, P2_ij );
  gsl_matrix_sub ( m2_ij, Q2_ij );

  gsl_matrix_memcpy ( m3_ij, P3_ij );
  gsl_matrix_sub ( m3_ij, Q3_ij );


  /* ===== full F-metric gF_ij */
  AMA = A  * cfg->Al1 + B  * cfg->Al2 + 2.0 * C  * cfg->Al3;

  gsl_matrix_memcpy (gF_ij, m1_ij );
  gsl_matrix_scale ( gF_ij, cfg->Al1 );	/* alpha1 m^1_ij */

  gsl_matrix_memcpy (mat1, m2_ij );
  gsl_matrix_scale ( mat1, cfg->Al2 );
  gsl_matrix_add ( gF_ij, mat1 );	/* + alpha2 m^2_ij */

  gsl_matrix_memcpy (mat1, m3_ij );
  gsl_matrix_scale ( mat1, 2.0 * cfg->Al3 );
  gsl_matrix_add ( gF_ij, mat1 );	/* + 2 * alpha3 m^3_ij */

  gsl_matrix_scale ( gF_ij, 1.0 / AMA );

  /* ===== averaged F-metric gFav_ij */
  gsl_matrix_memcpy (gFav_ij, m1_ij );
  gsl_matrix_scale ( gFav_ij, B );	/* B m^1_ij */

  gsl_matrix_memcpy (mat1, m2_ij );
  gsl_matrix_scale ( mat1, A );
  gsl_matrix_add ( gFav_ij, mat1 );	/* + A m^2_ij */

  gsl_matrix_memcpy (mat1, m3_ij );
  gsl_matrix_scale ( mat1, - 2.0 * C );
  gsl_matrix_add ( gFav_ij, mat1 );	/* - 2C m^3_ij */

  gsl_matrix_scale ( gFav_ij, 1.0 / ( 2.0 * D ) ); /* 1/ (2D) */


  /* ----- free memory ----- */
  gsl_vector_free ( dPhi_i );

  gsl_matrix_free ( mat1 );
  gsl_matrix_free ( mat2 );
  gsl_vector_free ( vec );

  gsl_matrix_free ( P1_ij );
  gsl_matrix_free ( P2_ij );
  gsl_matrix_free ( P3_ij );

  gsl_vector_free ( a2_dPhi_i );
  gsl_vector_free ( b2_dPhi_i );
  gsl_vector_free ( ab_dPhi_i );

  gsl_matrix_free ( Q1_ij );
  gsl_matrix_free ( Q2_ij );
  gsl_matrix_free ( Q3_ij );

  gsl_matrix_free ( a2_a2 );
  gsl_matrix_free ( a2_b2 );
  gsl_matrix_free ( a2_ab );

  gsl_matrix_free ( b2_b2 );
  gsl_matrix_free ( b2_ab );

  gsl_matrix_free ( ab_ab );


  return 0;

} /* computeFstatMetric() */

/* calculate pure Phase-metric gij, which must be allocated 4x4 matrix.
 *
 * Return 0 = OK, -1 on error.
 */
int
computePhaseMetric ( gsl_matrix *g_ij, const PhaseDerivs *dPhi, const REAL8Vector *GLweights )
{
  UINT4 i;
  UINT4 numSteps = GLweights->length;

  gsl_vector *dPhi_i;
  gsl_matrix *dPhi_i_dPhi_j;

  gsl_matrix *aPhi_ij;
  gsl_vector *aPhi_i;
  gsl_matrix *aPhi_i_aPhi_j;

  /* ----- check input ----- */
  if ( !g_ij  )
    return -1;

  if ( (g_ij->size1 != METRIC_DIM) || (g_ij->size2 != METRIC_DIM) )
    return -1;

  if ( dPhi->dAlpha->length != numSteps )
    return -1;

  /* ----- allocate matrices/vectors ----- */
  dPhi_i = gsl_vector_calloc ( METRIC_DIM );
  dPhi_i_dPhi_j = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  aPhi_i = gsl_vector_calloc ( METRIC_DIM );
  aPhi_ij = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );
  aPhi_i_aPhi_j = gsl_matrix_calloc ( METRIC_DIM, METRIC_DIM );

  /* ----- calculate averages using Gauss-Legendre integration ----- */
  for ( i = 0; i < numSteps; i ++ )
    {
      REAL8 wi = GLweights->data[i];
      gsl_vector_set ( dPhi_i, 0, dPhi->dFreq->data[i] );
      gsl_vector_set ( dPhi_i, 1, dPhi->dAlpha->data[i]);
      gsl_vector_set ( dPhi_i, 2, dPhi->dDelta->data[i]);
      gsl_vector_set ( dPhi_i, 3, dPhi->df1dot->data[i]);

      outer_product(dPhi_i_dPhi_j, dPhi_i, dPhi_i );

      /* Gauss-Legendre integration: weighted sum */
      gsl_vector_scale ( dPhi_i, wi );
      gsl_matrix_scale ( dPhi_i_dPhi_j, wi );

      /* ----- av_dPhi_i ----- */
      gsl_vector_add ( aPhi_i, dPhi_i );

      /* ----- av_dPhi_i_dPhi_j ----- */
      gsl_matrix_add ( aPhi_ij, dPhi_i_dPhi_j );

    } /* for i < numSteps */

  /* ---------- composite quantities ---------- */
  outer_product (aPhi_i_aPhi_j, aPhi_i, aPhi_i );

  /* ===== final answer: m1_ij, m2_ij, m3_ij ===== */
  gsl_matrix_memcpy ( g_ij, aPhi_ij );
  gsl_matrix_sub ( g_ij, aPhi_i_aPhi_j );

  /* ----- free memory ----- */
  gsl_vector_free ( dPhi_i );
  gsl_matrix_free ( dPhi_i_dPhi_j );
  gsl_matrix_free ( aPhi_ij );
  gsl_vector_free ( aPhi_i );
  gsl_matrix_free ( aPhi_i_aPhi_j );

  return 0;

} /* computePhaseMetric() */

/**
 * basic initializations: set-up 'ConfigVariables'
 * Taken from FstatMetric where it parsed user-input into ConfigVariables,
 * now basically just translates from modern-API 'metricParams' into old-API 'ConfigVariables'
 *
 */
int
InitCode ( ConfigVariables *cfg,
           const DopplerMetricParams *metricParams,
           const EphemerisData *edat
           )
{
  XLAL_CHECK ( cfg != NULL, XLAL_EINVAL );
  XLAL_CHECK ( metricParams != NULL, XLAL_EINVAL );
  XLAL_CHECK ( edat != NULL, XLAL_EINVAL );

  XLAL_CHECK ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL );
  XLAL_CHECK ( metricParams->segmentList.length == 1, XLAL_EDOM );	// only 1 segment allowed

  LIGOTimeGPSVector *GLtimestamps = NULL;
  REAL8Vector *detWeights;

  cfg->edat = edat;

  cfg->startTime = metricParams->segmentList.segs[0].start;
  cfg->refTime = metricParams->signalParams.Doppler.refTime;

  REAL8 startTimeREAL8 = XLALGPSGetREAL8 ( &(metricParams->segmentList.segs[0].start) );
  REAL8 endTimeREAL8   = XLALGPSGetREAL8 ( &(metricParams->segmentList.segs[0].end) );
  REAL8 duration = XLALGPSDiff ( &(metricParams->segmentList.segs[0].end), &(metricParams->segmentList.segs[0].start) );

  /* NOTE: internally we always deal with a fixed number of spins */
  XLAL_CHECK ( (cfg->dopplerPoint.fkdot = XLALCreateREAL8Vector ( NUM_SPINS )) != NULL, XLAL_EFUNC );

  cfg->dopplerPoint.fkdot->data[0] = metricParams->signalParams.Doppler.fkdot[0];
  cfg->dopplerPoint.fkdot->data[1] = metricParams->signalParams.Doppler.fkdot[1];

  cfg->dopplerPoint.skypos.system    = COORDINATESYSTEM_EQUATORIAL;
  cfg->dopplerPoint.skypos.longitude = metricParams->signalParams.Doppler.Alpha;
  cfg->dopplerPoint.skypos.latitude  = metricParams->signalParams.Doppler.Delta;

  /* ----- construct Gauss-Legendre timestamps and corresponding weights
   * for computing the integrals by Gauss-Legendre quadrature
   */
  UINT4 numSteps = 2000;	// just kept the default value from lalapps_FstatMetric, which was used in testMetricCodes.py

  REAL8Vector *ti;	/* temporary */
  REAL8 Tinv = 1.0 / duration;

  XLAL_CHECK ( (cfg->GLweights = XLALCreateREAL8Vector ( numSteps )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (ti = XLALCreateREAL8Vector ( numSteps )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (GLtimestamps = XLALCreateTimestampVector ( numSteps )) != NULL, XLAL_EFUNC );

  /* compute Gauss-Legendre roots, timestamps and associated weights */
  gauleg ( startTimeREAL8, endTimeREAL8, ti->data, cfg->GLweights->data, numSteps );

  /* convert REAL8-times into LIGOTimeGPS-times */
  for ( UINT4 i=0; i < numSteps; i ++ )
    {
      XLALGPSSetREAL8 ( &(GLtimestamps->data[i]), ti->data[i] );
      cfg->GLweights->data[i] *= Tinv;
    }

  XLALDestroyREAL8Vector ( ti );

  /* ----- initialize IFOs and (Multi-)DetectorStateSeries  ----- */
  UINT4 numDet = metricParams->multiIFO.length;

  XLAL_CHECK ( (cfg->multiDetStates = XLALCalloc ( 1, sizeof( *(cfg->multiDetStates) ))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( (cfg->multiDetStates->data = XLALCalloc (numDet, sizeof( *(cfg->multiDetStates->data)))) != NULL, XLAL_ENOMEM );

  cfg->multiDetStates->length = numDet;
  cfg->multiDetStates->Tspan = duration;

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      const LALDetector *ifo = &(metricParams->multiIFO.sites[X]);
      /* obtain detector positions and velocities, together with LMSTs */
      cfg->multiDetStates->data[X] = XLALGetDetectorStates ( GLtimestamps, ifo, edat, 0 );
      XLAL_CHECK ( cfg->multiDetStates->data[X] != NULL, XLAL_EFUNC );
    } /* for X < numDet */

  /* ----- get relative detector-weights from user ---------- */
  XLAL_CHECK ( (detWeights = XLALCreateREAL8Vector ( numDet )) != NULL, XLAL_EFUNC );

  for ( UINT4 X=0; X < numDet ; X ++ )
    {
      detWeights->data[X] = metricParams->multiNoiseFloor.sqrtSn[X];
    }

  /* ---------- combine relative detector-weights with GL-weights ----------*/
  MultiNoiseWeights *tmp;
  XLAL_CHECK ( (tmp = XLALCalloc(1, sizeof(*tmp))) != NULL, XLAL_ENOMEM );
  tmp->length = numDet;
  XLAL_CHECK ( (tmp->data = XLALCalloc( numDet, sizeof(*(tmp->data)))) != NULL, XLAL_ENOMEM );

  for ( UINT4 X = 0; X < numDet; X ++ )
    {
      /* create k^th weights vector */
      XLAL_CHECK ( (tmp->data[X] = XLALCreateREAL8Vector ( numSteps )) != NULL, XLAL_EFUNC );

      /* set all weights to detectorWeight * GLweight */
      for ( UINT4 i = 0; i < numSteps; i ++ )
        {
          tmp->data[X]->data[i] = detWeights->data[X] * cfg->GLweights->data[i];
        } // for i < numSteps

    } /* for X < numDet */

  cfg->multiNoiseWeights = tmp;

  cfg->multiAMcoe = XLALComputeMultiAMCoeffs ( cfg->multiDetStates, cfg->multiNoiseWeights, cfg->dopplerPoint.skypos );
  XLAL_CHECK ( cfg->multiAMcoe != NULL, XLAL_EFUNC );

  /* ----- compute amplitude-factors alpha1, alpha2, alpha3 ----- */
  REAL8 cosi = metricParams->signalParams.Amp.cosi;
  REAL8 psi  = metricParams->signalParams.Amp.psi;

  REAL8 Aplus = 0.5 * ( 1.0 + SQ(cosi) );
  REAL8 Across = cosi;
  REAL8 cos2psi = cos(2.0 * psi );
  REAL8 sin2psi = sin(2.0 * psi );
  cfg->Al1 = SQ(Aplus) * SQ( cos2psi ) + SQ(Across) * SQ(sin2psi);
  cfg->Al2 = SQ(Aplus) * SQ( sin2psi ) + SQ(Across) * SQ(cos2psi);
  cfg->Al3 = ( SQ(Aplus) - SQ(Across) ) * sin2psi * cos2psi ;

  // ----- translate 'new-API' DetectorMotionType into 'old-API' PhaseType_t
  switch ( (int)metricParams->detMotionType )
    {
    case DETMOTION_SPIN | DETMOTION_ORBIT:
      cfg->phaseType = PHASE_FULL;
      break;

    case DETMOTION_ORBIT:
      cfg->phaseType = PHASE_ORBITAL;
      break;

    case DETMOTION_SPIN:
      cfg->phaseType = PHASE_SPIN;
      break;

    case DETMOTION_SPIN | DETMOTION_PTOLEORBIT:
      cfg->phaseType = PHASE_PTOLE;
      break;

    default:
      XLAL_ERROR ( XLAL_EDOM, "Can't deal with detMotionType = %d, not analog exists for XLALOldDopplerFstatMetric()\n", metricParams->detMotionType );
      break;
    } // switch(detMotionType)

  /* free temporary memory */
  XLALDestroyREAL8Vector ( detWeights );
  XLALDestroyTimestampVector ( GLtimestamps);

  return XLAL_SUCCESS;

} /* InitFStat() */


/**
 * calculate the phase-derivatives \f$\partial_i \phi \f$ for the
 * time-series detStates and the given doppler-point.
 * Has the option of using only the orbital part of the phase (PHASE_ORBITAL)
 * or the full-phase (PHASE_FULL).
 *
 * returned PhaseDerivs is allocated in here.
 */
MultiPhaseDerivs *
getMultiPhaseDerivs ( const MultiDetectorStateSeries *multiDetStates,
                      const DopplerPoint *dopplerPoint,
                      PhaseType_t phaseType
                      )
{
  XLAL_CHECK_NULL ( multiDetStates != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( dopplerPoint != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( dopplerPoint->skypos.system == COORDINATESYSTEM_EQUATORIAL, XLAL_EDOM );
  XLAL_CHECK_NULL ( dopplerPoint->fkdot->length == 2, XLAL_EDOM );

  UINT4 numDet = multiDetStates->length;
  REAL8 Alpha = dopplerPoint->skypos.longitude;
  REAL8 Delta = dopplerPoint->skypos.latitude;

  REAL8 vn[3];		/* unit-vector pointing to source in Cart. equatorial coord. */
  vn[0] = cos(Delta) * cos(Alpha);
  vn[1] = cos(Delta) * sin(Alpha);
  vn[2] = sin(Delta);

  LIGOTimeGPS refTimeGPS = multiDetStates->data[0]->data[0].tGPS;	/* use 1st detectors startTime as refTime */
  REAL8 refTime = XLALGPSGetREAL8 ( &refTimeGPS );

  /* get tAutumn */
  PulsarTimesParamStruc XLAL_INIT_DECL(times);
  times.epoch = refTimeGPS;
  XLAL_CHECK_NULL ( XLALGetEarthTimes ( &(times.epoch), &(times.tMidnight), &(times.tAutumn) ) == XLAL_SUCCESS, XLAL_EFUNC );

  MultiPhaseDerivs *mdPhi = NULL;
  XLAL_CHECK_NULL( (mdPhi = XLALCalloc ( 1, sizeof( *mdPhi ))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( (mdPhi->data = XLALCalloc ( numDet, sizeof( *(mdPhi->data) ))) != NULL, XLAL_ENOMEM );

  mdPhi->length = numDet;

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      UINT4 numStepsX = multiDetStates->data[X]->length;
      DetectorStateSeries *detStatesX = multiDetStates->data[X];
      PhaseDerivs *dPhi;

      XLAL_CHECK_NULL( (dPhi = XLALCalloc ( 1, sizeof ( *dPhi ) )) != NULL, XLAL_ENOMEM );

      dPhi->dFreq  = XLALCreateREAL8Vector ( numStepsX );
      dPhi->dAlpha = XLALCreateREAL8Vector ( numStepsX );
      dPhi->dDelta = XLALCreateREAL8Vector ( numStepsX );
      dPhi->df1dot = XLALCreateREAL8Vector ( numStepsX );

      mdPhi->data[X] = dPhi;

      for (UINT4 i=0; i < numStepsX; i++ )
	{
	  REAL8 ti, dT, taui;
	  REAL8 rX[3];	/* radius-vector to use: full, orbital or spin-only */
	  REAL8 rDet[3];	/* vector from earth center to detector */
	  REAL8 fi;		/* instantaneous intrinsic frequency in SSB */
	  PosVel_t posvel;

	  ti = XLALGPSGetREAL8 ( &(detStatesX->data[i].tGPS) ) - refTime;

	  /* compute detector-vector relative to earth's center */
	  rDet[0] = detStatesX->data[i].rDetector[0] - detStatesX->data[i].earthState.posNow[0];
	  rDet[1] = detStatesX->data[i].rDetector[1] - detStatesX->data[i].earthState.posNow[1];
	  rDet[2] = detStatesX->data[i].rDetector[2] - detStatesX->data[i].earthState.posNow[2];

	  switch ( phaseType )
	    {
	    case PHASE_FULL:
	      COPY_VECT( rX, detStatesX->data[i].rDetector );
	      break;
	    case PHASE_ORBITAL:
	      COPY_VECT (rX, detStatesX->data[i].earthState.posNow );
	      break;
	    case PHASE_SPIN: /* just given for completeness, probably not very useful */
	      COPY_VECT ( rX, rDet );
	      break;
	    case PHASE_PTOLE: /* use Ptolemaic orbital approximation */
	      getPtolePosVel( &posvel, ti, times.tAutumn );
	      COPY_VECT ( rX, posvel.pos );
	      /* add on the detector-motion due to the Earth's spin */

	      rX[0] += rDet[0];
	      rX[1] += rDet[1];
	      rX[2] += rDet[2];

	      break;
	    default:
	      XLAL_ERROR_NULL ( XLAL_EINVAL, "Unknown phase-type specified '%d'\n", phaseType);
	      break;
	    } /* switch(phaseType) */

	  /* correct for time-delay from SSB to detector */
	  dT = SCALAR(vn, rX );
	  taui = ( ti + dT );

	  fi = dopplerPoint->fkdot->data[0] + taui * dopplerPoint->fkdot->data[1];

	  /* phase-derivatives */
	  dPhi->dFreq->data[i] 	= LAL_TWOPI * taui;

	  dPhi->dAlpha->data[i]	= LAL_TWOPI * fi * cos(Delta)*( - rX[0] * sin(Alpha) + rX[1] * cos(Alpha) );

	  dPhi->dDelta->data[i] = LAL_TWOPI * fi * ( - rX[0] * cos(Alpha) * sin(Delta)
                                                     - rX[1] * sin(Alpha) * sin(Delta)
                                                     + rX[2] * cos(Delta) );

	  dPhi->df1dot->data[i] = LAL_TWOPI * 0.5 * SQ( taui );

	} /* for i < numStepsX */

    } /* for X < numDet */

  return mdPhi;

} /* getMultiPhaseDerivs() */

/**
 * Calculate the projected metric onto the subspace of 'c' given by
 * ret_ij = g_ij - ( g_ic * g_jc / g_cc ) , where c is the value of the projected coordinate
 * The output-matrix ret must be allocated
 *
 * return 0 = OK, -1 on error.
 */
int
project_metric( gsl_matrix *ret_ij, gsl_matrix * g_ij, const UINT4 c )
{
  UINT4 i,j;

  if ( !ret_ij )
    return -1;
  if ( !g_ij )
    return -1;
  if ( (ret_ij->size1 != ret_ij->size2) )
    return -1;
  if ( (g_ij->size1 != g_ij->size2) )
    return -1;

  for ( i=0; i < ret_ij->size1; i ++) {
    for ( j=0; j < ret_ij->size2; j ++ ) {
      if ( i==c || j==c ) {
	gsl_matrix_set ( ret_ij, i, j, 0.0 );
      }
      else {
	gsl_matrix_set ( ret_ij, i, j, ( gsl_matrix_get(g_ij, i, j) - (gsl_matrix_get(g_ij, i, c) * gsl_matrix_get(g_ij, j, c) / gsl_matrix_get(g_ij, c, c)) ));
      }
    }
  }
  return 0;
}


/**
 * Calculate the outer product ret_ij of vectors u_i and v_j, given by
 * ret_ij = u_i v_j
 * The output-matrix ret must be allocated and have dimensions len(u) x len(v)
 *
 * return 0 = OK, -1 on error.
 */
int
outer_product (gsl_matrix *ret_ij, const gsl_vector *u_i, const gsl_vector *v_j )
{
  UINT4 i, j;

  if ( !ret_ij || !u_i || !v_j )
    return -1;

  if ( (ret_ij->size1 != u_i->size) || ( ret_ij->size2 != v_j->size) )
    return -1;


  for ( i=0; i < ret_ij->size1; i ++)
    for ( j=0; j < ret_ij->size2; j ++ )
      gsl_matrix_set ( ret_ij, i, j, gsl_vector_get(u_i, i) * gsl_vector_get(v_j, j) );

  return 0;

} /* outer_product() */


/* symmetrize the input-matrix 'mat' (which must be quadratic)
 */
int
symmetrize ( gsl_matrix *mat )
{
  gsl_matrix *tmp;

  if ( !mat )
    return -1;
  if ( mat->size1 != mat->size2 )
    return -1;

  tmp = gsl_matrix_calloc ( mat->size1, mat->size2 );

  gsl_matrix_transpose_memcpy ( tmp, mat );

  gsl_matrix_add (mat, tmp );

  gsl_matrix_scale ( mat, 0.5 );

  gsl_matrix_free ( tmp );

  return 0;

} /* symmetrize() */


/* compute the quadratic form = vec.mat.vec
 */
REAL8
quad_form ( const gsl_matrix *mat, const gsl_vector *vec )
{
  UINT4 i, j;
  REAL8 ret = 0;

  if ( !mat || !vec )
    return 0;
  if ( (mat->size1 != mat->size2) || ( mat->size1 != vec->size ) )
    return 0;

  for (i=0; i < mat->size1; i ++ )
    for (j=0; j < mat->size2; j ++ )
      ret += gsl_vector_get(vec, i) * gsl_matrix_get (mat, i, j) * gsl_vector_get(vec,j);

  return ret;

} /* quad_form() */


/**
 * Get Ptolemaic position and velocity at time tGPS
 * cut-down version of LALDTBaryPtolemaic()
 */

void
getPtolePosVel( PosVel_t *posvel, REAL8 tGPS, REAL8 tAutumnGPS )
{
  REAL8 rev;   /* Earth revolution angle, in radians. */

  /* Some local constants. */
  REAL8 tRev = LAL_AU_SI / LAL_C_SI;
  REAL8 vRev = LAL_TWOPI * tRev / LAL_YRSID_SI;
  REAL8 cosi = cos(LAL_IEARTH);
  REAL8 sini = sin(LAL_IEARTH);

  if ( !posvel )
    return;

  rev = LAL_TWOPI * ( tGPS - tAutumnGPS ) /LAL_YRSID_SI;

  /* Get detector position. */
  posvel->pos[0] = tRev * cos(rev);
  posvel->pos[1] = tRev * sin(rev) * cosi;
  posvel->pos[2]=  tRev * sin(rev) * sini;

  /* Get detector velocity. */
  posvel->vel[0] = -vRev * sin(rev);
  posvel->vel[1] =  vRev * cos(rev) * cosi;
  posvel->vel[2] =  vRev * cos(rev) * sini;

  return;

} /* getPtolePosVel() */

void
XLALDestroyMultiPhaseDerivs ( MultiPhaseDerivs *mdPhi )
{
  UINT4 numDet, X;

  if ( !mdPhi )
    return;
  numDet = mdPhi->length;

  if ( numDet && !mdPhi->data )
    XLAL_ERROR_VOID ( XLAL_EINVAL );

  for ( X=0; X < numDet; X ++ )
    {
      PhaseDerivs *dPhiX = mdPhi->data[X];
      if ( !dPhiX )
	continue;
      XLALDestroyREAL8Vector ( dPhiX->dAlpha);
      XLALDestroyREAL8Vector ( dPhiX->dDelta );
      XLALDestroyREAL8Vector ( dPhiX->dFreq );
      XLALDestroyREAL8Vector ( dPhiX->df1dot );

      LALFree ( dPhiX );

    } /* for X < numDet */

  LALFree ( mdPhi->data );
  LALFree ( mdPhi );

  return;

} /* XLALDestroyMultiPhaseDerivs() */

/* ================================================================================*/
/* taken from Numerical recipes:
 * function to compute roots and weights for order-N Gauss-Legendre integration
 * modified to use standard C-conventions for arrays (indexed 0... N-1)
 */
#define EPS 3.0e-11

/* Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
 * arrays x[0..n-1] and w[0..n-1] of length n, containing the abscissas and weights of the Gauss-
 * Legendre n-point quadrature formula.
 */
void
gauleg(double x1, double x2, double x[], double w[], int n)
{
     int m,j,i;
     /* High precision is a good idea for this routine. */
     double z1,z,xm,xl,pp,p3,p2,p1;

     /* The roots are symmetric in the interval, so
      * we only have to find half of them. */
     m=(n+1)/2;

     xm=0.5*(x2+x1);
     xl=0.5*(x2-x1);
     /* Loop over the desired roots. */
     for (i=1;i<=m;i++) {
       z=cos( LAL_PI *(i-0.25)/(n+0.5));
       /* Starting with the above approximation to the ith root,
	* we enter the main loop of refinement by Newton's method.
	*/
       do {
	 p1=1.0;
	 p2=0.0;
	 /* Loop up the recurrence relation to get the
	  * Legendre polynomial evaluated at z.
	  */
	 for (j=1;j<=n;j++) {

	   p3=p2;
	   p2=p1;
	   p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	 }
	 /* p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
	  * by a standard relation involving also p2, the polynomial of one lower order.
	  */
	 pp=n*(z*p1-p2)/(z*z-1.0);
	 z1=z;
	 /* Newton's method. */
	 z=z1-p1/pp;
       } while (fabs(z-z1) > EPS);
       /* Scale the root to the desired interval,  and put in its symmetric counterpart. */
       x[i - 1]=xm-xl*z;	/*RP: changed to C-convention */

       x[n+1-i - 1]=xm+xl*z; 	/*RP: changed to C-convention */

       /* Compute the weight and its symmetric counterpart. */
       w[i-1]=2.0*xl/((1.0-z*z)*pp*pp); /*RP: changed to C-convention */

       w[n+1-i - 1]=w[i - 1];		/*RP: changed to C-convention */
     }

} /* gauleg() */

