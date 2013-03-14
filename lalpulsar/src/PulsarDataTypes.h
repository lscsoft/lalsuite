/*
 * Copyright (C) 2004, 2005 Reinhard Prix
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

#ifndef _PULSARDATATYPES_H  /* Double-include protection. */
#define _PULSARDATATYPES_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \author Reinhard Prix
 * \date 2005
 * \defgroup PulsarDataTypes_h Header PulsarDataTypes.h
 * \ingroup pkg_pulsarCommon
 * \brief Some common useful data-types for pulsar-searches.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/PulsarDataTypes.h>
 * \endcode
 *
 */
/*@{*/

#include <gsl/gsl_matrix.h>

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#include <lal/SFTutils.h>

/** maximal number of spin-parameters (Freq + spindowns) we can handle */
#define PULSAR_MAX_SPINS	7

/** maximal number of detectors we can handle (for static arrays of detector quantities) */
#ifndef PULSAR_MAX_DETECTORS	// allow this value to be overridden for e.g. E@H apps
#define PULSAR_MAX_DETECTORS	10
#endif

/** Type defining the orbital parameters of a binary pulsar */
typedef struct tagBinaryOrbitParams {
  LIGOTimeGPS tp;         /**< time of observed periapsis passage (in SSB) */
  REAL8 argp;            /**< argument of periapsis (radians) */
  REAL8 asini;            /**< projected, normalized orbital semi-major axis (s) */
  REAL8 ecc;              /**< orbital eccentricity */
  REAL8 period;           /**< orbital period (sec) */
} BinaryOrbitParams;

/** Type containing the JKS 'amplitude parameters' {h0, cosi, phi0, psi} */
typedef struct tagPulsarAmplitudeParams {
  REAL8 h0;	/**< overall signal amplitude */
  REAL8 cosi;	/**< cos(iota) of inclination angle iota of spin-axis wrt line of sight */
  REAL8 psi;	/**< polarization angle psi */
  REAL8 phi0;	/**< initial signal-phase (at some reference time) */
} PulsarAmplitudeParams;

/** Struct for 'canonical' coordinates in amplitude-params space A^mu = {A1, A2, A3, A4} */
typedef REAL8 PulsarAmplitudeVect[4];

/** Typedef for fixed-size array holding GW frequency and derivatives fk = d^k Freq/dt^k|(tau_ref) */
typedef REAL8 PulsarSpins[PULSAR_MAX_SPINS];

/** Contains a "spin-range", ie spins \f$f^{(k)}\f$ and corresponding bands \f$\Delta f^{(k)}\f$
 *  at a given (SSB) reference GPS-time \f$\tau\f$.
 * "Canonical" ordering refers to \f$\Delta f^{(k)} >= 0\f$ for all k.
 */
typedef struct tagPulsarSpinRange
{
  LIGOTimeGPS refTime;		/**< SSB reference GPS-time at which spin-range is defined */
  PulsarSpins fkdot;		/**< Vector of spin-values \f$f^{(k)}\f$ */
  PulsarSpins fkdotBand;	/**< Vector of spin-bands \f$\Delta f^{(k)}\f$, MUST be same length as fkdot */
} PulsarSpinRange;

/** Type containing the 'Doppler-parameters' affecting the time-evolution of the phase */
typedef struct tagPulsarDopplerParams {
  LIGOTimeGPS refTime;	/**< reference time of pulsar parameters (in SSB!) */
  REAL8 Alpha;		/**< skyposition: RA (longitude) in equatorial coords and radians */
  REAL8 Delta;		/**< skyposition: DEC (latitude) in equatorial coords and radians */
  PulsarSpins fkdot;	/**< intrinsic spins: [Freq, f1dot, f2dot, ... ] where fkdot = d^kFreq/dt^k */
  BinaryOrbitParams *orbit;	/**< binary orbital parameters (or NULL if isolated pulsar) */
  REAL8 dFreq;		/**< if used for a frequency band: frequency increment between bins */
  UINT4 numFreqBins;	/**< if used for a frequency band: number of frequency bins */
} PulsarDopplerParams;

// ---------- transient-CW related types ----------
/** Struct to define parameters of a 'transient window' to be applied to obtain transient signals */
typedef enum {
  TRANSIENT_NONE = 0,		/**< Note: in this case the window-parameters will be ignored, and treated as rect={data},
                                 * i.e. a simple rectangular window covering all the data => this should always reproduce the
                                 * standard F-statistic computation. */
  TRANSIENT_RECTANGULAR = 1,	/**< standard rectangular window covering [t0, t0+tau] */
  TRANSIENT_EXPONENTIAL,	/**< exponentially decaying window e^{-t0/tau} starting at t0.
                                 * Note: we'll truncate this at some small (eg 3x) e-folding TRANSIENT_EXP_EFOLDING */
  TRANSIENT_LAST
} transientWindowType_t;

/** Struct defining one transient window instance */
typedef struct tagtransientWindow_t
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  UINT4 t0;			/**< GPS start-time 't0' */
  UINT4 tau;			/**< transient timescale tau in seconds */
} transientWindow_t;


// ---------- 'integrated' types describing a complete CW signal ----------

/** Type defining the parameters of a pulsar-source of CW Gravitational waves */
typedef struct tagPulsarParams {
  PulsarAmplitudeParams Amp;		/**< 'Amplitude-parameters': h0, cosi, phi0, psi */
  PulsarDopplerParams   Doppler;	/**< 'Phase-evolution parameters': {skypos, fkdot, orbital params } */
  transientWindow_t     Transient;	/**< Transient window-parameters (start-time, duration, window-type) */
} PulsarParams;

/** Type containing a "candidate": parameter-space point with estimated errors and Fstat-value/significance */
typedef struct tagPulsarCandidate {
  PulsarAmplitudeParams Amp, dAmp;	/**< amplitude-parameters and error-estimates */
  PulsarDopplerParams Doppler, dDoppler;/**< Doppler-parameters and error-bars */
  REAL8 significance;			/**< a (user-chosen) measure of 'significance': Fstat, Hough-count,... */
  gsl_matrix *AmpFisherMatrix;		/**< Fisher-matrix of amplitude-subspace: has more info than dAmp! */
} PulsarCandidate;


/** DEPRECATED Type defining the parameters of a pulsar-source of Gravitational waves.
 * \note this type is obsolete and should no longer be used,
 * however, it's too entrenched in LALGeneratePulsarSignal() et al, and codes using it,
 * so we can't easily get rid of it and keep it for now....
 * \deprecated Use #PulsarParams instead.
 */
typedef struct tagPulsarSourceParams {
   LIGOTimeGPS refTime;	/**< reference time of pulsar parameters (in SSB!) */
   SkyPosition position; /**< source location (in radians) */
   REAL4 psi;            /**< polarization angle (radians) at tRef */
   REAL4 aPlus; 	/**< plus-polarization amplitude at tRef */
   REAL4 aCross;  	/**< cross-polarization amplitude at tRef */
   REAL8 phi0;           /**< initial phase (radians) at tRef */
   REAL8 f0;             /**< WAVE-frequency(!) at tRef (in Hz) */
   REAL8Vector *spindown;/**< wave-frequency spindowns at tRef (NOT f0-normalized!) */
} PulsarSourceParams;

/*---------- Global variables ----------*/
/** \name empty struct initializers */
/*@{*/
extern const BinaryOrbitParams empty_BinaryOrbitParams;
extern const PulsarAmplitudeParams empty_PulsarAmplitudeParams;
extern const PulsarSpinRange empty_PulsarSpinRange;
extern const PulsarDopplerParams empty_PulsarDopplerParams;
extern const PulsarParams empty_PulsarParams;
extern const PulsarCandidate empty_PulsarCandidate;

/* initializers for frequently used non-pulsar types */
extern const LIGOTimeGPS empty_LIGOTimeGPS;
extern const BarycenterInput empty_BarycenterInput;
extern const EphemerisData empty_EphemerisData;
extern const SkyPosition empty_SkyPosition;
/*@}*/

/*@}*/

#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
