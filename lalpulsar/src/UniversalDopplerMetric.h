/*
 * Copyright (C) 2012 Karl Wette
 * Copyright (C) 2008, 2009 Reinhard Prix
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

#ifndef _UNIVERSALDOPPLERMETRIC_H  /* Double-include protection. */
#define _UNIVERSALDOPPLERMETRIC_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup UniversalDopplerMetric_h Header UniversalDopplerMetric.h
 * \ingroup pkg_pulsarMetric
 * \author Reinhard Prix, Karl Wette
 *
 * Function to compute the full F-statistic metric, including
 * antenna-pattern functions from multi-detector, derived in \ref Prix07.
 *
 */
/*@{*/

/*---------- INCLUDES ----------*/
#include <math.h>

/* gsl includes */
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>

#include <lal/LALDatatypes.h>
#include <lal/LALBarycenter.h>
#include <lal/XLALGSL.h>
#include <lal/DetectorStates.h>
#include <lal/SFTutils.h>
#include <lal/Segments.h>

/*---------- exported types ----------*/
  /** 2D vector */
typedef REAL8 vect2D_t[2];
  /** 3D vector */
typedef REAL8 vect3D_t[3];
  /** 3x3 matrix, useful for spatial 3D vector operations */
typedef REAL8 mat33_t[3][3];

/** variable-length list of 2D-vectors */
typedef struct tagvect2Dlist_t {
  UINT4 length;			/**< number of elements */
  vect2D_t *data;		/**< array of 2D vectors */
} vect2Dlist_t;

/** variable-length list of 3D vectors */
typedef struct tagvect3Dlist_t {
  UINT4 length;			/**< number of elements */
  vect3D_t *data;		/**< array of 3D vectors */
} vect3Dlist_t;


/** Small Container to hold two 3D vectors: position and velocity */
typedef struct tagPosVel3D_t {
  vect3D_t pos;
  vect3D_t vel;
} PosVel3D_t;


/** Different types of detector-motion to use in order to compute the Doppler-metric,
 * the most 'realistic' obviously being DETMOTION_EPHEMORBIT_SPIN, which includes
 * both orbital and spin motion, and uses the Earth-ephemeris.
 */
typedef enum {
  DETMOTION_SPIN_ORBIT,		/**< full ephemeris-based detector motion (orbit+spin) */
  DETMOTION_ORBIT,		/**< ephemeris-based, purely orbital detector-motion, no Earth spin */
  DETMOTION_SPIN,		/**< purely Earth-spin detector motion (no orbit) */

  DETMOTION_SPIN_PTOLEORBIT,	/**< ptole-orbital motion (on a circle) + Earth spin */
  DETMOTION_PTOLEORBIT,		/**< pure "Ptolemaic" orbital motion, no Earth spin */

  DETMOTION_ORBIT_SPINZ,	/**< orbital motion plus *only* z-component of Earth spin-motion wrt to ecliptic plane */
  DETMOTION_ORBIT_SPINXY,	/**< orbital motion plus *only* x+y component of Earth spin-motion in the ecliptic */

  DETMOTION_LAST
} DetectorMotionType;


typedef enum {
  METRIC_TYPE_PHASE = 0,	/**< compute phase metric only */
  METRIC_TYPE_FSTAT = 1,	/**< compute full F-metric only */
  METRIC_TYPE_ALL   = 2,	/**< compute both F-metric and phase-metric */
  METRIC_TYPE_LAST
} MetricType_t;


/** enum listing symbolic 'names' for all Doppler Coordinates
 * supported by the metric codes in FstatMetric
 */
typedef enum {
  DOPPLERCOORD_NONE = -1,	/**< No Doppler component */

  DOPPLERCOORD_FREQ,		/**< Frequency [Units: Hz]. */
  DOPPLERCOORD_F1DOT,		/**< First spindown [Units: Hz/s]. */
  DOPPLERCOORD_F2DOT,		/**< Second spindown [Units: Hz/s^2]. */
  DOPPLERCOORD_F3DOT,		/**< Third spindown [Units: Hz/s^3]. */

  DOPPLERCOORD_GC_NU0,		/**< Global correlation frequency [Units: Hz]. Activates 'reduced' detector position. */
  DOPPLERCOORD_GC_NU1,		/**< Global correlation first spindown [Units: Hz/s]. Activates 'reduced' detector position. */
  DOPPLERCOORD_GC_NU2,		/**< Global correlation second spindown [Units: Hz/s^2]. Activates 'reduced' detector position. */
  DOPPLERCOORD_GC_NU3,		/**< Global correlation third spindown [Units: Hz/s^3]. Activates 'reduced' detector position. */

  DOPPLERCOORD_ALPHA,		/**< Right ascension [Units: radians]. Uses 'reduced' detector position. */
  DOPPLERCOORD_DELTA,		/**< Declination [Units: radians]. Uses 'reduced' detector position. */

  DOPPLERCOORD_N2X_EQU,		/**< X component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */
  DOPPLERCOORD_N2Y_EQU,		/**< Y component of contrained sky position in equatorial coordinates [Units: none]. Uses 'reduced' detector position. */

  DOPPLERCOORD_N2X_ECL,		/**< X component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */
  DOPPLERCOORD_N2Y_ECL,		/**< Y component of contrained sky position in ecliptic coordinates [Units: none]. Uses 'reduced' detector position. */

  DOPPLERCOORD_N3X_EQU,		/**< X component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
  DOPPLERCOORD_N3Y_EQU,		/**< Y component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
  DOPPLERCOORD_N3Z_EQU,		/**< Z component of unconstrained super-sky position in equatorial coordinates [Units: none]. */

  DOPPLERCOORD_N3X_ECL,		/**< X component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
  DOPPLERCOORD_N3Y_ECL,		/**< Y component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
  DOPPLERCOORD_N3Z_ECL,		/**< Z component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */

  DOPPLERCOORD_N3SX_EQU,	/**< X spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */
  DOPPLERCOORD_N3SY_EQU,	/**< Y spin-component of unconstrained super-sky position in equatorial coordinates [Units: none]. */

  DOPPLERCOORD_N3OX_ECL,	/**< X orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */
  DOPPLERCOORD_N3OY_ECL,	/**< Y orbit-component of unconstrained super-sky position in ecliptic coordinates [Units: none]. */

  DOPPLERCOORD_LAST
} DopplerCoordinateID;

#define DOPPLERMETRIC_MAX_DIM 60	/**< should be large enough for a long time ... */
/** type describing a Doppler coordinate system:
 * lists the number of dimensions and the symbolic names of the coordinates.
 */
typedef struct tagDopplerCoordinateSystem
{
  UINT4 dim;						/**< number of dimensions covered */
  DopplerCoordinateID coordIDs[DOPPLERMETRIC_MAX_DIM];	/**< coordinate 'names' */
} DopplerCoordinateSystem;

/** meta-info specifying a Doppler-metric
 */
typedef struct tagDopplerMetricParams
{
  DopplerCoordinateSystem coordSys;		/**< number of dimensions and coordinate-IDs of Doppler-metric */
  DetectorMotionType detMotionType;		/**< the type of detector-motion assumed: full spin+orbit, pure orbital, Ptole, ... */

  LALSegList segmentList;			/**< segment list: Nseg segments of the form (startGPS endGPS numSFTs) */
  MultiDetectorInfo detInfo;			/**< detectors (and their noise-weights) to compute metric for */

  PulsarParams signalParams;			/**< parameter-space point to compute metric for (doppler + amplitudes) */
  INT4 projectCoord;				/**< project metric onto subspace orthogonal to this axis (-1 = none, 0 = 1st coordinate, etc) */

  MetricType_t metricType;			/**< switch controlling which types of metric to compute: 0 = PhaseMetric g_ij, 1 = Fmetrics gF.., 2=BOTH */
  BOOLEAN approxPhase;				/**< use an approximate phase-model, neglecting Roemer delay in spindown coordinates */
} DopplerMetricParams;



/** Struct to hold the 'atoms', ie weighted phase-derivative averages like \f$\langle a^2 \partial_i \phi \partial_j \phi\rangle>\f$
 *  from which the F-metric is computed, but also the full Fisher-matrix. The noise-weighted average is defined as
 * \f$\langle Q\rangle \equiv \frac{1}{T} \, \sum_X w^X\, \int_0^T Q\, dt \f$, where \f$w^X\f$ is the noise-weight for detector X,
 * and \f$T\f$ is the observation time, see \ref Prix07 for details.
 */
typedef struct tagFmetricAtoms_t
{
  REAL8 a_a;		/**< \f$ \langle a^2 \rangle = A \f$ */
  REAL8 a_b;		/**< \f$ \langle a b \rangle = C \f$ */
  REAL8 b_b;		/**< \f$ \langle b^2 \rangle = B \f$ */

  gsl_vector *a_a_i;	/**< \f$ \langle a^2 \, \partial_i\phi \rangle \f$ */
  gsl_vector *a_b_i;	/**< \f$ \langle a\, b \, \partial_i\phi \rangle \f$ */
  gsl_vector *b_b_i;	/**< \f$ \langle b^2 \, \partial_i\phi \rangle \f$ */

  gsl_matrix *a_a_i_j;	/**< \f$ \langle a^2 \, \partial_i\phi \, \partial_j\phi \rangle \f$ */
  gsl_matrix *a_b_i_j;	/**< \f$ \langle a\,b \, \partial_i\phi \, \partial_j\phi \rangle \f$ */
  gsl_matrix *b_b_i_j;	/**< \f$ \langle b^2 \, \partial_i\phi \, \partial_j\phi \rangle \f$ */

  double maxrelerr;	/**< estimate for largest relative error in metric component integrations */
} FmetricAtoms_t;


/** struct to hold a DopplerMetric, including meta-info on the number of
 * dimensions, the coordinate-system and type of metric.
 */
typedef struct tagDopplerMetric
{
  DopplerMetricParams meta;		/**< "meta-info" describing/specifying the type of Doppler metric */

  gsl_matrix *g_ij;			/**< symmetric matrix holding the usual Phase-metric */
  double maxrelerr_gPh;			/**< estimate for largest relative error in phase-metric component integrations */

  gsl_matrix *gF_ij;			/**< full F-statistic metric gF_ij, including antenna-pattern effects (see \ref Prix07) */
  gsl_matrix *gFav_ij;			/**< 'average' Fstat-metric */
  gsl_matrix *m1_ij, *m2_ij, *m3_ij;	/**< Fstat-metric sub components */

  gsl_matrix *Fisher_ab;		/**< Full 4+n dimensional Fisher matrix, ie amplitude + Doppler space */

  double maxrelerr_gF;			/**< estimate for largest relative error in Fmetric component integrations */

  REAL8 rho2;				/**< signal SNR rho^2 = A^mu M_mu_nu A^nu */
} DopplerMetric;


/*---------- Global variables ----------*/
extern const PosVel3D_t empty_PosVel3D_t;
extern const DopplerMetricParams empty_DopplerMetricParams;
extern const MultiDetectorInfo empty_MultiDetectorInfo;
extern const DopplerCoordinateSystem empty_DopplerCoordinateSystem;
#define empty_vect3D_t {0,0,0}

/*---------- exported prototypes [API] ----------*/
gsl_matrix *
XLALDopplerPhaseMetric ( const DopplerMetricParams *metricParams,
			 const EphemerisData *edat,
                         double *relerr_max
			 );

DopplerMetric*
XLALDopplerFstatMetric ( const DopplerMetricParams *metricParams,
			 const EphemerisData *edat
			 );

DopplerMetric*
XLALDopplerFstatMetricCoh ( const DopplerMetricParams *metricParams,
                            const EphemerisData *edat
                            );

FmetricAtoms_t*
XLALComputeAtomsForFmetric ( const DopplerMetricParams *metricParams,
			     const EphemerisData *edat
			     );


int
XLALDetectorPosVel (PosVel3D_t *spin_posvel,
		    PosVel3D_t *orbit_posvel,
		    const LIGOTimeGPS *tGPS,
		    const LALDetector *site,
		    const EphemerisData *edat,
		    DetectorMotionType special
		    );

int XLALPtolemaicPosVel ( PosVel3D_t *posvel, const LIGOTimeGPS *tGPS );
gsl_matrix *XLALProjectMetric ( const gsl_matrix * g_ij, const UINT4 c );

void XLALequatorialVect2ecliptic ( vect3D_t out, const vect3D_t in );
void XLALeclipticVect2equatorial ( vect3D_t out, const vect3D_t in );
void XLALmatrix33_in_vect3 ( vect3D_t out, mat33_t mat, const vect3D_t in );

vect3Dlist_t *
XLALComputeOrbitalDerivatives ( UINT4 maxorder, const LIGOTimeGPS *tGPS, const EphemerisData *edat );

FmetricAtoms_t* XLALCreateFmetricAtoms ( UINT4 dim );
void XLALDestroyFmetricAtoms ( FmetricAtoms_t *atoms );

void XLALDestroyDopplerMetric ( DopplerMetric *metric );

int XLALParseDetectorMotionString ( const CHAR *detMotionString );
int XLALParseDopplerCoordinateString ( const CHAR *coordName );
int XLALDopplerCoordinateNames2System ( DopplerCoordinateSystem *coordSys, const LALStringVector *coordNames );

const CHAR *XLALDetectorMotionName ( DetectorMotionType detType );
const CHAR *XLALDopplerCoordinateName ( DopplerCoordinateID coordID );
const CHAR *XLALDopplerCoordinateHelp ( DopplerCoordinateID coordID );
CHAR *XLALDopplerCoordinateHelpAll ( void );

gsl_matrix* XLALNaturalizeMetric( const gsl_matrix* g_ij, const DopplerMetricParams *metricParams );

gsl_matrix *XLALDiagNormalizeMetric ( const gsl_matrix * g_ij );
int XLALAddDopplerMetric ( DopplerMetric **metric1, const DopplerMetric *metric2 );
int XLALScaleDopplerMetric ( DopplerMetric *m, REAL8 scale );
// destructor for vect3Dlist_t type
void XLALDestroyVect3Dlist ( vect3Dlist_t *list );

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
