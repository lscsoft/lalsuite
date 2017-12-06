/*
 * Copyright (C) 2011 J. Creighton
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

#ifndef _LALSIMBLACKHOLERINGDOWN_TIGER_H
#define _LALSIMBLACKHOLERINGDOWN_TIGER_H

#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

REAL8 XLALSimSphericalHarmonicPlus(UINT4 l, INT4 m, REAL8 iota);
REAL8 XLALSimSphericalHarmonicCross(UINT4 l, INT4 m, REAL8 iota);

REAL8 XLALSimRingdownQNMAmplitudes(INT4 l, INT4 m, REAL8 eta, REAL8 chiEff);

COMPLEX16 XLALSimRingdownFitOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a);

REAL8 XLALQNMFreqOfOmega(COMPLEX16 omega, REAL8 mtot);
REAL8 XLALQNMTauOfOmega(COMPLEX16 omega, REAL8 mtot);

/*REAL8 XLALSimSphericalHarmonicPlusFD(UINT4 l, INT4 m, REAL8 iota);
REAL8 XLALSimSphericalHarmonicCrossFD(UINT4 l, INT4 m, REAL8 iota);

REAL8 XLALSimRingdownQNMAmplitudesFD(INT4 l, INT4 m, REAL8 eta, REAL8 chiEff);

COMPLEX16 XLALSimRingdownFitOmegaFD(UINT4 l, INT4 m, UINT4 n, REAL8 a);

REAL8 XLALQNMFreqOfOmegaFD(COMPLEX16 omega, REAL8 mtot);
REAL8 XLALQNMTauOfOmegaFD(COMPLEX16 omega, REAL8 mtot); */

int XLALSimBlackHoleRingdownTigerAllModes(
    REAL8TimeSeries **hplus,      /**< plus-polarization waveform [returned] */
    REAL8TimeSeries **hcross,     /**< cross-polarization waveform [returned] */
    const LIGOTimeGPS *t0,                /**< start time of ringdown */
    REAL8 phi0,                   /**< initial phase of ringdown (rad) */
    REAL8 deltaT,                 /**< sampling interval (s) */
    REAL8 mass,                   /**< black hole mass (kg) */
    REAL8 a,      /**< black hole dimensionless spin parameter */
    REAL8 eta,         /**< symmetric mass ratio of progenitor */
    REAL8 spin1[3],    /**< initial spin for 1st component */
    REAL8 spin2[3],    /**< initial spin for 2nd component */
    REAL8 chiEff,      /**< effective spin parameter for initial spins */
    REAL8 distance,               /**< distance to source (m) */
    REAL8 inclination,            /**< inclination of source's spin axis (rad) */
    LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */
                                         );

int XLALSimBlackHoleRingdownTiger(
                                  REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
                                  REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
                                  SphHarmTimeSeries *hlms, /**< Head of linked list of waveform modes */
                                  const LIGOTimeGPS *t0,		/**< start time of ringdown */
                                  REAL8 phi0,			/**< initial phase of ringdown (rad) */
                                  REAL8 deltaT,			/**< sampling interval (s) */
                                  REAL8 mass,			/**< black hole mass (kg) */
                                  REAL8 a,	/**< black hole dimensionless spin parameter */
                                  //	double fractional_mass_loss,	/**< fraction of mass radiated in this mode */
                                  REAL8 eta,          /**< symmetric mass ratio of progenitor */
                                  REAL8 spin1[3],     /**< initial spin for 1st component */
                                  REAL8 spin2[3],     /**< initial spin for 2nd component */
                                  REAL8 chiEff,
                                  REAL8 distance,		/**< distance to source (m) */
                                  REAL8 inclination,		/**< inclination of source's spin axis (rad) */
                                  LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */ );

int XLALSimBlackHoleRingdownModeTiger(
                                      COMPLEX16TimeSeries **hlmmode, /**< complex waveform for lm mode */
                                      const LIGOTimeGPS *t0,         /**< start time of ringdown */
                                      REAL8 phi0,         /**< initial phase of ringdown (rad) */
                                      REAL8 deltaT,       /**< sampling interval (s) */
                                      REAL8 mass,         /**< black hole mass (kg) */
                                      REAL8 a,            /**< black hole dimensionless spin parameter */
                                      //	double fractional_mass_loss,	/**< fraction of mass radiated in this mode */
                                      REAL8 distance,     /**< distance to source (m) */
                                      REAL8 inclination,  /**< inclination of source's spin axis (rad) */
                                      REAL8 eta,          /**< symmetric mass ratio of progenitor */
                                      REAL8 spin1[3],     /**< initial spin for 1st component */
                                      REAL8 spin2[3],     /**< initial spin for 2nd component */
                                      REAL8 chiEff,
                                      UINT4 l,            /**< polar mode number */
                                      INT4 m,             /**< azimuthal mode number */
                                      REAL8 dfreq,        /**< relative shift in the real frequency parameter */
                                      REAL8 dtau          /**< relative shift in the damping time parameter */
                                      );
int XLALSimBlackHoleRingdownTigerFD(
                COMPLEX16FrequencySeries **hptilde,    /**< FD plus polarization */
                COMPLEX16FrequencySeries **hctilde,    /**< FD cross polarization */
                REAL8 phi0,                            /**<initial phase of ringdown*/
                REAL8 deltaF,                         /**<sampling interval (Hz)*/
                REAL8 fEnd,
                REAL8 fStart,
                REAL8 mass,                           /**<blackhole mass*/
                REAL8 a,                              /**<blackhole dimensionless spin parameter*/
                REAL8 eta,                            /**<symmetric mass ratio*/
                REAL8 chiEff,
                REAL8 distance,                       /**distance(m)*/
                REAL8 inclination,                    /**< inclination of source's spin axis (rad)*/
                LALSimInspiralTestGRParam *nonGRparams  /**< testing GR parameters */
                );

int XLALSimBlackHoleRingdownModeTigerFD(
                COMPLEX16FrequencySeries **hptilde_lm,    /**< FD plus polarization*/
                COMPLEX16FrequencySeries **hctilde_lm,    /**< FD cross polarization*/
                const REAL8 deltaF,                    /**< Frequency resolution */
                const REAL8 fStart,                    /**< Start GW frequency (Hz)*/
                const REAL8 fEnd,                      /**< Highest GW frequency (Hz)*/
                REAL8 phi0,                       /**< initial phase of ringdown */
                REAL8 mass,                       /**< black hole mass (kg) */
                REAL8 a,  /**< black hole dimensionless spin parameter */
                REAL8 distance,           /**< distance to source (m) */
                REAL8 inclination,                /**< inclination of source's spin axis (rad)*/
                REAL8 eta,         /**< symmetric mass ratio of progenitor */
                UINT4 l,
                INT4 m,
                REAL8 chiEff,      /**< effective spin parameter for initial spins */
                REAL8 dfreq,         /**< relative shift in thereal frequency parameter */
                REAL8 dtau          /**< relative shift in the damping time parameter */
                );
#endif
