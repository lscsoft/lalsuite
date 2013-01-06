/*
*  Copyright (C) 2012 Miroslav Shaltev, R Prix
*  Copyright (C) 2007 Curt Cutler, Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#ifndef _LALBARYCENTER_H    /* Protect against double-inclusion */
#define _LALBARYCENTER_H

#include <stdio.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup LALBarycenter_h Header LALBarycenter.h
 * \ingroup pkg_pulsarCommon
 * \author Curt Cutler
 * \date 2001
 *
 * \brief Provides routines for transforming from arrival time at detector (GPS) to pulse emission time (TDB); ie
 * for ``barycentering'' the measured astronomical time series.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/LALInitBarycenter.h>
 * #include <lal/LALBarycenter.h>
 * \endcode
 *
 */
/*@{*/

/** \name Error codes */
/*@{*/
#define LALBARYCENTERH_ENULL  2	/**< Null input to Barycenter routine. */
#define LALBARYCENTERH_EOUTOFRANGEE  4	/**< tgps not in range of earth.dat file */
#define LALBARYCENTERH_EOUTOFRANGES  8	/**< tgps not in range of sun.dat file */
#define LALBARYCENTERH_EBADSOURCEPOS 16	/**< source position not in standard range */
#define LALBARYCENTERH_EXLAL 	     32	/**< XLAL function failed. */
/*@}*/

/** \cond DONT_DOXYGEN */
#define LALBARYCENTERH_MSGENULL  "Null input to Barycenter routine."
#define LALBARYCENTERH_MSGEOUTOFRANGEE  "tgps not in range of earth.dat file"
#define LALBARYCENTERH_MSGEOUTOFRANGES  "tgps not in range of sun.dat file"
#define LALBARYCENTERH_MSGEBADSOURCEPOS "source position not in standard range"
#define LALBARYCENTERH_MSGEXLAL 	"XLAL function failed."
/** \endcond */

/** \brief Enumerated type denoting the time system type to be produced in
 * the solar system barycentring routines.
 *
 * The type denotes the time system in which solar system barycentred times
 * should be given. \c TIMECORRECTION_TDB and \c TIMECORRECTION_TEMPO will mean
 * times are in the Barycentric Dynamical Time (TDB) system, where the
 * conversion has been performed using a time correction ephemeris look-up
 * table as used by TEMPO2 (\c TIMECORRECTION_TEMPO is so-called because the
 * pulsar timing software TEMPO uses the TDB time system by default); \c
 * TIMECORRECTION_ORIGINAL will mean times are in the TDB system, but with the
 * conversion performed using the original \c XLALBarycenterEarth function; and,
 * \c TIMECORRECTION_TCB and \c TIMECORRECTION_TEMPO2 will mean times are in the
 * Coordinate Barycentric Time (TCB) system, where the conversion has been
 * performed using a time correction ephemeris look-up table as used by TEMPO2
 * (\c TIMECORRECTION_TEMPO2 is so-called because the pulsar timing software
 * TEMPO2 uses the TCB time system by default). */
typedef enum{
  TIMECORRECTION_NONE = 0,
  TIMECORRECTION_TDB,
  TIMECORRECTION_TCB,
  TIMECORRECTION_TEMPO,
  TIMECORRECTION_TEMPO2,
  TIMECORRECTION_ORIGINAL,
  TIMECORRECTION_LAST
} TimeCorrectionType;

/** \brief Enumerated type denoting the JPL solar system ephemeris to be used
 * in calculating barycentre time corrections.
 */
typedef enum {
  EPHEM_NONE = 0,
  EPHEM_DE200,
  EPHEM_DE405,
  EPHEM_DE414,
  EPHEM_DE421,
  EPHEM_LAST
} EphemerisType;

/** \name Constants from Irwin and Fukushima, A&A, 348, 1999 (taken from TEMPO2)
 * used for ephemeris conversions. */
/*@{*/
#define IFTE_JD0  2443144.5003725 /**< Epoch of TCB, TCG and TT in Julian Days */
#define IFTE_MJD0 43144.0003725 /**< Epoch of TCB, TCG and TT in Modified Julian Days */
#define IFTE_TEPH0 -65.564518e-6 /**< Equation 17 of Irwin and Fukushima. */
#define IFTE_LC 1.48082686742e-8 /**< Equation 20 of Irwin and Fukushima. */
#define IFTE_KM1 1.55051979176e-8 /**< Value of K-1, defined using the IAU definition of L_B = 1.55051976772e-8 and K=1/(1-L_B) (see TEMPO2). */
#define IFTE_K (((long double)1.0) + ((long double)IFTE_KM1)) /**< Factor relating ephemeris units for time and distance to corresponding SI units, from Eq. 2 of Irwin and Fukushima. */
/*@}*/

#define JPL_AU_DE405 149597870.6910000 	/**< Definition of 1 AU from the JPL DE405 ephemeris in km */
#define JPL_AU_DE200 149597870.6600000 	/**< Definition of 1 AU from the JPL DE200 ephemeris in km */
#define CURT_AU 149597870.6600 		/**< 1 AU from create_solar_system_barycenter.c as used in Curt's original routines */

/** \brief This structure contains
 * two pointers to the ephemeris data files containing arrays
 * of center-of-mass positions for the Earth and Sun, respectively.
 *
 * The tables are derived from the JPL ephemeris.
 *
 * Files tabulate positions for one calendar year
 * (actually, a little more than one year, to deal
 * with overlaps).  The first line of each table summarizes
 * what is in it. Subsequent lines give the time (GPS) and the
 * Earth's position \f$(x,y,z)\f$,
 * velocity \f$(v_x, v_y, v_z)\f$, and acceleration \f$(a_x, a_y, a_z)\f$
 * at that instant.  All in units of seconds; e.g. positions have
 * units of seconds, and accelerations have units 1/sec.
 *
 */
typedef struct tagEphemerisFilenames
{
  CHAR *earthEphemeris;         /**< File containing Earth's position.  */
  CHAR *sunEphemeris;           /**< File containing Sun's position. */
}
EphemerisFilenames;

/** Structure holding a REAL8 time, and a position, velocity and acceleration vector. */
typedef struct tagPosVelAcc
{
  REAL8 gps;            /**< REAL8 timestamp */
  REAL8 pos[3];         /**< position-vector */
  REAL8 vel[3];         /**< velocity-vector */
  REAL8 acc[3];         /**< acceleration-vector */
}
PosVelAcc;

/** This structure contains all information about the
 * center-of-mass positions of the Earth and Sun, listed at regular
 * time intervals.
 */
typedef struct tagEphemerisData
{
  EphemerisFilenames ephiles; /**< Names of the two files containing positions of
                               * Earth and Sun, respectively at evenly spaced times. */
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(1D_ARRAY(PosVelAcc, ephemE, INT4, nentriesE));
  SWIGLAL(1D_ARRAY(PosVelAcc, ephemS, INT4, nentriesS));
#endif /* SWIG */
  INT4  nentriesE;      /**< The number of entries in Earth ephemeris table. */
  INT4  nentriesS;      /**< The number of entries in Sun ephemeris table. */

  REAL8 dtEtable;       /**< The spacing in sec between consecutive instants in Earth ephemeris table.*/
  REAL8 dtStable;       /**< The spacing in sec between consecutive instants in Sun ephemeris table.*/

  PosVelAcc *ephemE;    /**< Array containing pos,vel,acc of earth, as extracted from earth
                         * ephem file. Units are sec, 1, 1/sec respectively */
  PosVelAcc *ephemS;    /**< Array with pos, vel and acc for the sun (see ephemE) */

  EphemerisType etype;  /**< The ephemeris type e.g. DE405 */
}
EphemerisData;


/** This structure will contain a vector of time corrections
 * used during conversion from TT to TDB/TCB/Teph */
typedef struct tagTimeCorrectionData{
  CHAR *timeEphemeris;   /**< File containing the time ephemeris */

  UINT4  nentriesT;      /**< The number of entries in Time ephemeris table. */
  REAL8 dtTtable;        /**< The spacing in sec between consecutive instants in Time ephemeris table.*/
  REAL8 *timeCorrs;      /**< Array of time delays for converting TT to TDB/TCB from the Time table (seconds).*/
  REAL8 timeCorrStart;   /**< The initial GPS time of the time delay table. */
} TimeCorrectionData;


/** Basic output structure of LALBarycenterEarth.c.
 */
typedef struct tagEarthState
{
  REAL8  einstein;      /**<  the einstein delay equiv TDB - TDT or TCB - TDT */
  REAL8 deinstein;      /**< d(einstein)/d(tgps) */

  REAL8 posNow[3];      /**< Cartesian coords of Earth's center at tgps,
                         * extrapolated from JPL DE405 ephemeris; units= sec */
  REAL8 velNow[3];      /**< dimensionless velocity of Earth's center at tgps,
                         * extrapolated from JPL DE405 ephemeris */

  REAL8 gmstRad;        /**< Greenwich Mean Sidereal Time (GMST) in radians, at tgps */
  REAL8 gastRad;        /**< Greenwich Apparent Sidereal Time, in radians, at tgps;
                         * Is basically the angle thru which Earth has spun at
                         * given time; gast is like gmst, but has
                         * additional correction for short-term nutation */

  REAL8 tzeA;           /**< variable describing effect of lunisolar precession, at tgps */
  REAL8 zA;             /**< variable describing effect of lunisolar precession, at tgps */
  REAL8 thetaA;         /**< variable describing effect of lunisolar precession, at tgps */
  REAL8 delpsi;         /**< variable describing effect of Earth nutation, at tgps*/
  REAL8 deleps;         /**< variable describing effect of Earth nutation, at tgps*/

  REAL8 se[3];          /**< vector that points from Sun to Earth at instant tgps,
                         * in DE405 coords; units = sec */
  REAL8 dse[3];         /**< d(se[3])/d(tgps); Dimensionless */
  REAL8 rse;            /**< length of vector se[3]; units = sec */
  REAL8 drse;           /**< d(rse)/d(tgps); dimensionless */

  TimeCorrectionType ttype; /**< Time correction type */
}
EarthState;

/** Basic input structure to LALBarycenter.c.
 */
typedef struct tagBarycenterInput
{
  LIGOTimeGPS  tgps;    /**< input GPS arrival time. I use tgps (lower case)
                         * to remind that here the LAL structure is a
                         * field in the larger structure BarycenterInput.
                         * I use tGPS as an input structure (by itself) to
                         * LALBarycenterEarth */

  LALDetector site;     /**< detector site info.  <b>NOTE:</b>
                         * the <tt>site.location</tt> field must be modified
                         * to give the detector location in units of
                         * <em>seconds</em> (i.e. divide the values normally
                         * stored in <tt>site.location</tt> by <tt>LAL_C_SI</tt> */

  REAL8 alpha;          /**<  Source right ascension in ICRS J2000 coords (radians). */
  REAL8 delta;          /**< Source declination in ICRS J2000 coords (radians) */
  REAL8 dInv;           /**< 1/(distance to source), in 1/sec.
                         * This is needed to correct Roemer delay for very
                         * nearby sources; correction is about 10 microsec for
                         * source at 100 pc */
}
BarycenterInput;

/*Curt: probably best to take 1.0 OUT of tDot--ie., output tDot-1.
But most users would immediately add back the one anyway.
*/

/*Curt: rem te is ``time pulse would arrive at a GPS clock
way out in empty space, if you renormalized  and zero-ed the latter
to give, on average, the same arrival time as the GPS clock on Earth'' */

/**  Basic output structure produced by LALBarycenter.c.
 */
typedef struct tagEmissionTime
{
  REAL8 deltaT;         /**< \f$t_e\f$(TDB) - \f$t_a\f$(GPS)
                         * + (light-travel-time from source to SSB) */

  LIGOTimeGPS te;       /**< pulse emission time (TDB); also sometimes called
                         * ``arrival time (TDB) of same wavefront at SSB'' */
  REAL8 tDot;           /**< d(emission time in TDB)/d(arrival time in GPS)  */

  REAL8 rDetector[3];   /**< Cartesian coords (0=x,1=y,2=z) of detector position
                         * at $t_a$ (GPS), in ICRS J2000 coords. Units = sec. */

  REAL8 vDetector[3];   /* Cartesian coords (0=x,1=y,2=z) of detector velocity
                         * at \f$t_a\f$ (GPS), in ICRS J2000 coords. Dimensionless. */

  REAL8 roemer;         /**<  the Roemer delay */
  REAL8 droemer;        /**<  d(Roemer)/d(tgps) */

  REAL8 shapiro;        /**<  the Shapiro delay */
  REAL8 dshapiro;       /**<  d(Shapiro)/d(tgps) */

  REAL8 erot;           /**< Earth rotation delay */
  REAL8 derot;          /**< d(erot)/d(tgps) */
}
EmissionTime;


/// internal (opaque) buffer type for optimized Barycentering function
typedef struct tagBarycenterBuffer BarycenterBuffer;

/* Function prototypes. */
int XLALBarycenterEarth ( EarthState *earth, const LIGOTimeGPS *tGPS, const EphemerisData *edat);
int XLALBarycenter ( EmissionTime *emit, const BarycenterInput *baryinput, const EarthState *earth);
int XLALBarycenterOpt ( EmissionTime *emit, const BarycenterInput *baryinput, const EarthState *earth, BarycenterBuffer **buffer);

/* Function that uses time delay look-up tables to calculate time delays */
int XLALBarycenterEarthNew ( EarthState *earth,
                             const LIGOTimeGPS *tGPS,
                             const EphemerisData *edat,
                             const TimeCorrectionData *tdat,
                             TimeCorrectionType ttype );

/* Function to calculate positions */
void precessionMatrix( REAL8 prn[3][3],
                       REAL8 mjd,
                       REAL8 dpsi,
                       REAL8 deps );
void observatoryEarth( REAL8 obsearth[3],
                       const LALDetector det,
                       const LIGOTimeGPS *tgps,
                       REAL8 gmst,
                       REAL8 dpsi,
                       REAL8 deps );

// deprecated LAL interface

void LALBarycenterEarth ( LALStatus *status, EarthState *earth, const LIGOTimeGPS *tGPS, const EphemerisData *edat);
void LALBarycenter ( LALStatus *status, EmissionTime *emit, const BarycenterInput *baryinput, const EarthState *earth);


/*@}*/

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */
