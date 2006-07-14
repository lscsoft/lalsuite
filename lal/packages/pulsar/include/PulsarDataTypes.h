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

/**
 * \author Reinhard Prix
 * \date 2005
 * \file 
 * \ingroup pulsarCommon
 * \brief Some common useful data-types for pulsar-searches.
 *
 * $Id$
 *
 */

#ifndef _PULSARDATATYPES_H  /* Double-include protection. */
#define _PULSARDATATYPES_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

#include "SFTutils.h"

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PULSARDATATYPESH, "$Id$");

/** Type defining the orbital parameters of a binary pulsar */
typedef struct {
  LIGOTimeGPS orbitEpoch; /**< time of periapsis passage (in SSB) */
  REAL8 omega;            /**< argument of periapsis (radians) */
  REAL8 rPeriNorm;        /**< projected, normalized periapsis (s) */
  REAL8 oneMinusEcc;      /**< 1 - orbital eccentricity */
  REAL8 angularSpeed;     /**< angular speed at periapsis (Hz) */
} BinaryOrbitParams;
  
/** Type containing the JKS 'amplitude parameters' {h0, cosi, phi0, psi} */
typedef struct {
  REAL8 h0;	/**< overall signal amplitude */
  REAL8 cosi;	/**< cos(iota) of inclination angle iota of spin-axis wrt line of sight */
  REAL8 psi;	/**< polarization angle psi */
  REAL8 phi0;	/**< initial signal-phase (at some reference time) */
} PulsarAmplitudeParams;

/** Type containing the 'Doppler-parameters' affecting the time-evolution of the phase */
typedef struct {
  LIGOTimeGPS refTime;	/**< reference time of pulsar parameters (in SSB!) */
  REAL8 Alpha;		/**< skyposition: RA (longitude) in equatorial coords and radians */ 
  REAL8 Delta;		/**< skyposition: DEC (latitude) in equatorial coords and radians */ 
  REAL8 Freq;		/**< intrinsic signal-frequency at refTime */
  REAL8 f1dot;		/**< dFreq/dt in pulsar frame */
  REAL8 f2dot;		/**< d^2Freq/dt^2 in pulsar frame */
  REAL8 f3dot;		/**< d^3Freq/dt^3 in pulsar frame */
  BinaryOrbitParams *orbit;	/**< binary orbital parameters (or NULL if isolated pulsar) */
} PulsarDopplerParams;

/** Type defining the parameters of a pulsar-source of Gravitational waves */
typedef struct {
  PulsarAmplitudeParams Amp;	/**< 'Amplitude-parameters': h0, cosi, phi0, psi */
  PulsarDopplerParams Doppler;	/**< 'Doppler-parameters': {skypos, fkdot, orbital params } */
} PulsarParams;

/** Type containing a "candidate": parameter-space point with estimated errors and Fstat-value/significance */
typedef struct {
  PulsarAmplitudeParams Amp, dAmp;	/**< amplitude-parameters and error-estimates */
  PulsarDopplerParams Doppler, dDoppler;/**< Doppler-parameters and error-bars */
  REAL8 significance;			/**< a (user-chosen) measure of 'significance': Fstat, Hough-count,... */
} PulsarCandidate;


/** [OBSOLETE] Type defining the parameters of a pulsar-source of Gravitational waves.
 * NOTE: this type is obsolete and should no longer be used [==> use 'PulsarParams' instead]
 * however, it's too entrenched in the the GeneratePulsarSignal() functions and codes using it,
 * so we can't easily get rid of it any more, so we keep it for now....
 */
typedef struct {
   LIGOTimeGPS tRef;	/**< reference time of pulsar parameters (in SSB!) */
   SkyPosition position;	/**< source location (in radians) */
   REAL4 psi;            /**< polarization angle (radians) at tRef */
   REAL4 aPlus; 		/**< plus-polarization amplitude at tRef */
   REAL4 aCross;  	/**< cross-polarization amplitude at tRef */
   REAL8 phi0;           /**< initial phase (radians) at tRef */
   REAL8 f0;             /**< WAVE-frequency(!) at tRef (in Hz) */
   REAL8Vector *spindown;/**< wave-frequency spindowns at tRef (NOT f0-normalized!) */
} PulsarSourceParams;

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
