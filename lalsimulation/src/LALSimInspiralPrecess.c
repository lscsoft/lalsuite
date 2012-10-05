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

/* Considerations...
 * 1. Instead of h_lm, have the user give the mode decomposition for a specific
 * l, and provide that l to the function. Then everything proceeds as before
 * with the assumption that the double pointer provided has 2l+1 values to 
 * operate on.
 * 
 * 2. Don't have full time series for alpha, beta, and gamma unless we have a
 * specific need for it. It's wasteful of memory.
 */

#include "LALSimInspiralPrecess.h"

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
				COMPLEX16TimeSeries* h_2_2, /**< (2,-2) mode, modified in place */
				COMPLEX16TimeSeries* h_22, /**< (2,2) mode, modified in place */
				double precess_freq, /**< Precession frequency in Hz */
				double a, /**< Opening angle of precession cone in rads  */
				double phi_precess, /**< initial phase in cone of L around J */
				double alpha_0, /**< azimuth btwn center of cone and line of sight */
				double beta_0 /**< zenith btwn center of cone and line of sight */
) {

		// Error checking
		// Since we need at least three points to do any of the numerical work
		// we intend to, if the waveform is two points or smaller, we're in
		// trouble. I don't expect this to be a problem.
		if( h_2_2->data->length <= 2 ){
			XLALPrintError( "XLAL Error - %s: Waveform length is too small to evolve accurately.", __func__);
			XLAL_ERROR( XLAL_EBADLEN );
		}
        if( h_2_2->data->length != h_22->data->length ){
            XLALPrintError( "XLAL Error - %s: Input (2,2) and (2,-2) modes have different length.", __func__);
            XLAL_ERROR( XLAL_EBADLEN );
        }

		unsigned int i;
		double omg_p = 2*LAL_PI*precess_freq;
		double t=0;
		int l = 2, mp;
		// h_lm samples
		complex double *x_lm = XLALCalloc( (2*l+1), sizeof(complex double) );

		// time evolved Euler angles
		REAL8TimeSeries* alpha = XLALCreateREAL8TimeSeries(
			"euler angle alpha",
			&(h_22->epoch),
			h_22->f0,
			h_22->deltaT,
			&(h_22->sampleUnits),
			h_22->data->length
		);
		REAL8TimeSeries* beta = XLALCreateREAL8TimeSeries(
			"euler angle beta",
			&(h_22->epoch),
			h_22->f0,
			h_22->deltaT,
			&(h_22->sampleUnits),
			h_22->data->length
		);
		REAL8TimeSeries* gam = XLALCreateREAL8TimeSeries(
			"euler angle gamma",
			&(h_22->epoch),
			h_22->f0,
			h_22->deltaT,
			&(h_22->sampleUnits),
			h_22->data->length
		);

		// Minimal rotation constraint
		// \gamma(t) = \int \cos(\beta(t)) \alpha'(t) dt
		// Uses the second order finite difference to estimate dalpha/dt
		// Then the trapezoid rule for the integration
		for(i=0; i<alpha->data->length; i++){
			t = h_22->deltaT*i;
			alpha->data->data[i] = a*sin(omg_p * t + phi_precess) + alpha_0;
			beta->data->data[i] = a*cos(omg_p * t + phi_precess) + beta_0;
		}

		// NOTE: The step size cancels out between the difference and sum and
		// thus does not appear below.
		// two point forward difference
		double dalpha_0 = alpha->data->data[1] - alpha->data->data[0];
		// three point central difference
		double dalpha_1 = 0.5*(alpha->data->data[2] - alpha->data->data[0]);
		gam->data->data[0] = 0.;
		gam->data->data[1] =
				cos(beta->data->data[0])*dalpha_0 +
				cos(beta->data->data[1])*dalpha_1;
		for(i=2; i<gam->data->length-1; i++){
			// three point central difference
			dalpha_0 = dalpha_1;
			dalpha_1 = 0.5*(alpha->data->data[i+1] - alpha->data->data[i-1]);
			// Two point numerical integration over the interval
			gam->data->data[i] = 
				gam->data->data[i-1] +
				cos(beta->data->data[i-1])*dalpha_0 +
				cos(beta->data->data[i])*dalpha_1;
		}
		// Use two point backward difference for last point
        dalpha_0 = dalpha_1;
        dalpha_1 = alpha->data->data[i] - alpha->data->data[i-1];
		gam->data->data[i] = gam->data->data[i-1] +
                cos(beta->data->data[i-1])*dalpha_0 +
                cos(beta->data->data[i])*dalpha_1;

		// Rotate waveform
		// TODO: Make a loop on m inside the mp loop and sum across all
		// submodes. Note that this requires a change to the function signature
		// to input all modes for a given l.
		for(i=0; i<h_22->data->length; i++){
			x_lm[4] = h_22->data->data[i];
			x_lm[0] = h_2_2->data->data[i];

			h_22->data->data[i] = 0;
			h_2_2->data->data[i] = 0;

			for(mp=-l; mp<=l; mp++){
				h_2_2->data->data[i] += 
				x_lm[mp+2] * XLALWignerDMatrix( 2, mp, -2, alpha->data->data[i], beta->data->data[i], gam->data->data[i] );
				h_22->data->data[i] += 
				x_lm[mp+2] * XLALWignerDMatrix( 2, mp, 2, alpha->data->data[i], beta->data->data[i], gam->data->data[i] );
			}
		}	

		XLALDestroyREAL8TimeSeries( alpha );
		XLALDestroyREAL8TimeSeries( beta );
		XLALDestroyREAL8TimeSeries( gam );

		return XLAL_SUCCESS;
}

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
) {
		int ret = XLALSimInspiralConstantPrecessionConeWaveformModes( 
						h_2_2, h_22, 
						precess_freq, a, phi_precess,
						alpha_0, beta_0 );
        if( ret != XLAL_SUCCESS )
            XLAL_ERROR( XLAL_EFUNC );

		if( !hp ){
			XLALDestroyREAL8TimeSeries( hp );
			hp = XLALCreateREAL8TimeSeries(
				"h_+ precessed waveform",
				&(h_22->epoch),
				h_22->f0,
				h_22->deltaT,
				&(h_22->sampleUnits),
				h_22->data->length
			);
		}
		if( !hx ){
			XLALDestroyREAL8TimeSeries( hx );
			hx = XLALCreateREAL8TimeSeries(
				"h_x precessed waveform",
				&(h_22->epoch),
				h_22->f0,
				h_22->deltaT,
				&(h_22->sampleUnits),
				h_22->data->length
			);
		}

		unsigned int i;
		complex double x_t = 0I;
		// FIXME: Should these be fixed?
		double view_th = 0.0, view_ph = 0.0;
		// Reconstitute the waveform from the h_lm
		for(i=0; i<h_22->data->length; i++){
			x_t = h_22->data->data[i] * XLALSpinWeightedSphericalHarmonic( view_th, view_ph, -2, 2, 2 );
			x_t += h_2_2->data->data[i] * XLALSpinWeightedSphericalHarmonic( view_th, view_ph, -2, 2, -2 );
			hp->data->data[i] = crealf( x_t );
			hx->data->data[i] = cimagf( x_t );
		}

		// User should do this.
		//XLALDestroyCOMPLEX16TimeSeries( h_22 );
		//XLALDestroyCOMPLEX16TimeSeries( h_2_2 );

        return XLAL_SUCCESS;
}
