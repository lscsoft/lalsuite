/*
 *  Copyright (C) 2012 Chris Pankow
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
 *
 * \author Chris Pankow
 *
 * \file
 *
 * \brief Functions to take an arbitrary waveform time series and impose the
 * effects of causing the viewing angle to precess about a cone of L around J.
 * The cone currently has a constant opening angle.
 *
 * */

/* Considerations:
 * Don't have full time series for alpha, beta, and gamma unless we have a
 * specific need for it. It's wasteful of memory.
 */

#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALDatatypes.h>
#include <lal/SphericalHarmonics.h>

/**
 * Takes in the h_lm spherical harmonic decomposed modes and rotates the modes
 * by Euler angles alpha, beta, and gamma using the Wigner D matricies.
 * 
 * e.g.
 *
 * \f$\tilde{h}_{l,m}(t) = D^l_{m,m'} h_{l,m'}(t)\f$
 */
int XLALSimInspiralPrecessionRotateModes(
				SphHarmTimeSeries* h_lm, /**< spherical harmonic decomposed modes, modified in place */
				REAL8TimeSeries* alpha, /**< alpha Euler angle time series */
				REAL8TimeSeries* beta, /**< beta Euler angle time series */
				REAL8TimeSeries* gam /**< gamma Euler angle time series */
);

/**
 * Takes in the l=2, abs(m)=2 decomposed modes as a strain time series and
 * imposes the effect of a constant cone of precession. This is accomplished
 * by taking the axis of the binary rotational plane and rotating the Y_lm
 * such that it appears to "precess" around a fixed J direction.
 * Note that h_2_2 and h_22 are modified in place.
 *
 * Future revisions will change the first two pointers to this:
 * COMPLEX16TimeSeries** h_lm
 *
 * and add 
 * unsigned int l
 *
 * Thus the input h_lm will be considered a list of pointers to the h_lm for a
 * given l and the appropriate action will be taken for *all* of the submodes.
 */
int XLALSimInspiralConstantPrecessionConeWaveformModes(
				COMPLEX16TimeSeries** h_2_2, /**< (2,-2) mode, modified in place */
				COMPLEX16TimeSeries** h_22, /**< (2,2) mode, modified in place */
				double precess_freq, /**< Precession frequency in Hz */
				double a, /**< Opening angle of precession cone in rads  */
				double phi_precess, /**< initial phase in cone of L around J */
				double alpha_0, /**< azimuth btwn center of cone and line of sight */
				double beta_0 /**< zenith btwn center of cone and line of sight */
);

/**
 * Takes in the l=2, abs(m)=2 decomposed modes as a strain time series and
 * imposes the effect of a constant cone of precession. The result is returned
 * in the physical waveforms hp, hx, after they have been resummed from the 
 * modified h_22 waveforms.
 *
 * NOTE: the modes h_2_2 and h_22 will be modified in place
 */
int XLALSimInspiralConstantPrecessionConeWaveform(
				REAL8TimeSeries* hp, /**< Output precessing plus polarization */
				REAL8TimeSeries* hx, /**< Output precessing cross polarization*/
				COMPLEX16TimeSeries* h_2_2, /**< Input non-precessing (2,-2) mode - modified in place */
				COMPLEX16TimeSeries* h_22, /**< Input non-precessing (2,2) mode - modified in place */
				double precess_freq, /**< Precession frequency in Hz */
				double a, /**< Opening angle of precession cone in rads  */
				double phi_precess, /**< initial phase in cone of L around J */
				double alpha_0, /**< azimuth btwn center of cone and line of sight */
				double beta_0 /**< zenith btwn center of cone and line of sight */
);

