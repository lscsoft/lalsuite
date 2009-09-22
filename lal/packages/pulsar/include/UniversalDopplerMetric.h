/*
 * Copyright (C) 2008 Reinhard Prix
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
 *
 * \author{Reinhard Prix}
 *
 * Function to compute the full F-statistic metric, including
 * antenna-pattern functions from multi-detector, derived in \ref Prix07.
 *
 * Revision: $Id$
 *
 */

#ifndef _UNIVERSALDOPPLERMETRIC_H  /* Double-include protection. */
#define _UNIVERSALDOPPLERMETRIC_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

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

/*---------- DEFINES ----------*/
NRCSID( UNIVERSALDOPPLERMETRICH, "$Id$");

/*----- Error-codes -----*/


/*---------- exported types ----------*/

/** Small Container to hold two 3D vectors: position and velocity */
typedef struct {
  REAL8 pos[3];
  REAL8 vel[3];
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

/** Array of symbolic 'names' for various detector-motions
 */
#ifdef IN_UNIVERSALDOPPLERMETRICC
const CHAR *DetectorMotionNames[] = {
  "spin+orbit",
  "orbit",
  "spin",

  "spin+ptoleorbit",
  "ptoleorbit",

  "orbit+spin_Z",
  "orbit+spin_XY",

  "NONE"
};
#endif

/** enum listing symbolic 'names' for all Doppler Coordinates
 * supported by the metric codes in FstatMetric
 */
typedef enum {
  DOPPLERCOORD_NONE = -1,		/**< used to denote 'empty', i.e. no Doppler component */
  DOPPLERCOORD_FREQ_SI = 0,		/**< frequency in Hz */
  DOPPLERCOORD_F1DOT_SI,		/**< f1dot = dFreq/dt in Hz/s */
  DOPPLERCOORD_F2DOT_SI,		/**< f2dot = d2Freq/dt2 in Hz/s^2 */
  DOPPLERCOORD_F3DOT_SI,		/**< f3dot = d3Freq/dt3 in Hz/s^3 */

  DOPPLERCOORD_ALPHA_RAD,		/**< right-ascencion (longitude) in radians, using coord-system of ephemeris-file */
  DOPPLERCOORD_DELTA_RAD,		/**< declination (latitude) in radians,  using coord-system of ephemeris-file */

  DOPPLERCOORD_FREQ_NAT,		/**< frequency in "natural units": om0 = 2pi f * Tspan */
  DOPPLERCOORD_F1DOT_NAT,		/**< f1dot in "natural units":     om1 = 2pi f1dot/2! * Tspan^2 */
  DOPPLERCOORD_F2DOT_NAT,		/**< f2dot in "natural units":     om2 = 2pi f2dot/3! * Tspan^3 */
  DOPPLERCOORD_F3DOT_NAT,		/**< f3dot in "natural units":     om3 = 2pi f3dot/4! * Tspan^4 */

  DOPPLERCOORD_LAST
} DopplerCoordinateID;

#ifdef IN_UNIVERSALDOPPLERMETRICC
/** Array of Doppler coordinate names, which *MUST*
 * correspond to the order in which they are listed in
 * DopplerCoordinateID
 */
const CHAR *DopplerCoordinateNames[] = {
  "Freq",
  "f1dot",
  "f2dot",
  "f3dot",

  "Alpha",
  "Delta",

  "Freq_Nat",
  "f1dot_Nat",
  "f2dot_Nat",
  "f3dot_Nat",

  "NONE"
};

/** Array of help-strings explaining the meaning/conventions of the Doppler coordinate names,
 * NOTE: this *MUST* correspond to the order in which they are listed in DopplerCoordinateID.
 *
 * NOTE2: It's important to also specify the "coordinate-set" this coordinate is meant to belong to,
 * in the sense of which other coordinates need to be held constant in partial derivatives wrt to this coordinate!
 *
 */
const CHAR *DopplerCoordinateNamesHelp[] = {
  "Signal frequency in SSB [Units:Hz]. Coordinate-set: {fkdot, sky}.",
  "First frequency-derivative dFreq/dtau in SSB [Units:Hz/s]. Coordinate-set: {fkdot, sky}.",
  "Second frequency-derivative d2Freq/dtau^2 in SSB [Units:Hz/s^2]. Coordinate-set: {fkdot, sky}.",
  "Third frequency-derivative d3Freq/dtau^3 in SSB [Units:Hz/s^3]. Coordinate-set: {fkdot, sky}.",

  "Sky-position: Right-ascension (longitude) wrt ephemeris coord-system [Units:rad]. Coordinate-set: {fkdot, Alpha, Delta}.",
  "Sky-position: Declination (latitude) wrt ephemeris coord-system [Units:rad]. Coordinate-set: {fkdot, Alpha, Delta}.",

  "Same as Freq, but in 'natural units': Freq_Nat = 2 pi Freq Tspan [Units:1]",
  "Same as f1dot, but in 'natural units': f1dot_Nat = 2 pi f1dot/2! Tspan^2 [Units:1]",
  "Same as f2dot, but in 'natural units': f2dot_Nat = 2 pi f2dot/3! Tspan^3 [Units:1]",
  "Same as f3dot, but in 'natural units': f3dot_Nat = 2 pi f3dot/4! Tspan^4 [Units:1]",

  "NONE"
};
#endif


#define DOPPLERMETRIC_MAX_DIM 60	/**< should be large enough for a long time ... */
/** type describing a Doppler coordinate system:
 * lists the number of dimensions and the symbolic names of the coordinates.
 */
typedef struct
{
  UINT4 dim;						/**< number of dimensions covered */
  DopplerCoordinateID coordIDs[DOPPLERMETRIC_MAX_DIM];	/**< coordinate 'names' */
} DopplerCoordinateSystem;

#define DOPPLERMETRIC_MAX_DETECTORS 60	/**< should be way large enough forever */
/** type describing a set of detectors and their relative noise-weights
 * This is only used for full multi-IFO Fstatistic-metrics
 */
typedef struct
{
  UINT4 length;						/**< number N of detectors */
  LALDetector sites[DOPPLERMETRIC_MAX_DETECTORS]; 	/**< array of N detectors */
  REAL8 detWeights[DOPPLERMETRIC_MAX_DETECTORS];	/**< array of N detector noise-weights: must satisfy \f$\sum_{i=1}^N w_i = 1\f$ */
} MultiDetectorInfo;

/**< meta-info specifying a Doppler-metric
 */
typedef struct
{
  DopplerCoordinateSystem coordSys;		/**< number of dimensions and coordinate-IDs of Doppler-metric */
  DetectorMotionType detMotionType;		/**< the type of detector-motion assumed: full spin+orbit, pure orbital, Ptole, ... */

  LIGOTimeGPS startTime;			/**< startTime of the observation */
  REAL8 Tspan;					/**< total spanned duration of the observation */
  MultiDetectorInfo detInfo;			/**< detectors (and their noise-weights) to compute metric for */

  PulsarParams signalParams;			/**< parameter-space point to compute metric for (doppler + amplitudes) */
} DopplerMetricParams;



/** Struct to hold the 'atoms', ie weighted phase-derivative averages like \f$\langle a^2 \partial_i \phi \partial_j \phi\rangle>\f$
 *  from which the F-metric is computed, but also the full Fisher-matrix. The noise-weighted average is defined as
 * \f$\langle Q\rangle \equiv \frac{1}{T} \, \sum_X w^X\, \int_0^T Q\, dt \f$, where \f$w^X\f$ is the noise-weight for detector X,
 * and \f$T\f$ is the observation time, see \ref Prix07 for details.
 */
typedef struct
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
typedef struct
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
extern DopplerMetricParams empty_DopplerMetricParams;
extern MultiDetectorInfo empty_MultiDetectorInfo;

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


FmetricAtoms_t*
XLALComputeAtomsForFmetric ( const DopplerMetricParams *metricParams,
			     const EphemerisData *edat
			     );


int
XLALDetectorPosVel ( PosVel3D_t *pos_vel3D,
		    const LIGOTimeGPS *tGPS,
		    const LALDetector *site,
		    const EphemerisData *edat,
		    DetectorMotionType special
		    );

int
XLALAmplitudeParams2Vect ( PulsarAmplitudeVect *Amu, const PulsarAmplitudeParams *Amp );


FmetricAtoms_t* XLALCreateFmetricAtoms ( UINT4 dim );
void XLALDestroyFmetricAtoms ( FmetricAtoms_t *atoms );

void XLALDestroyDopplerMetric ( DopplerMetric *metric );

DetectorMotionType XLALParseDetectorMotionString ( const CHAR *detMotionString );
DopplerCoordinateID XLALParseDopplerCoordinateString ( const CHAR *coordName );
int XLALDopplerCoordinateNames2System ( DopplerCoordinateSystem *coordSys, const LALStringVector *coordNames );

const CHAR *XLALDetectorMotionName ( DetectorMotionType detType );
const CHAR *XLALDopplerCoordinateName ( DopplerCoordinateID coordID );
const CHAR *XLALDopplerCoordinateHelp ( DopplerCoordinateID coordID );
CHAR *XLALDopplerCoordinateHelpAll ( void );
int XLALParseMultiDetectorInfo ( MultiDetectorInfo *detInfo, const LALStringVector *detNames, const LALStringVector *detWeights );


#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
