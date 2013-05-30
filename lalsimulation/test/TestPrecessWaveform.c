/*
 *  Copyright (C) 2012 Chris Pankow, Evan Ochsner
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
 * \brief Testing constant precession code on TaylorT1 waveform.
 *
 * */


#include <lal/LALSimInspiralPrecess.h>

int main(void){

		FILE* h_ref = fopen("h_ref.txt", "w");
		FILE* h_rot = fopen("h_rot.txt", "w");
		REAL8TimeSeries *hp = NULL, *hx = NULL;

		int ret;
		unsigned int i;
		
		// Waveform parameters
		REAL8 m1 = 2.0*LAL_MSUN_SI, m2 = 5.0*LAL_MSUN_SI;
		REAL8 f_min = 40.0, f_ref = 0., dist = 1e6*LAL_PC_SI;
		REAL8 lambda1 = 0.0, lambda2 = 0.0;
		REAL8 phi = 0.0, dt = 1/16384.0;
		REAL8 inclination = LAL_PI_4, psi = 0.;
		LALSimInspiralWaveformFlags *waveFlags = NULL;
		LALSimInspiralTestGRParam *nonGRparams = NULL;
		int Lmax = 5, amplitudeOrder = -1, phaseOrder = -1;
		Approximant approximant = TaylorT1;

		// Parameters define a constant precession cone
		REAL8 precess_freq = 10.; // Freq. of L's motion about cone (Hz)
		REAL8 cone_opening = LAL_PI_4; // Opening angle of precession cone
		REAL8 cone_azimuth = 0.; // Initial azimuthal angle of L on its cone
		REAL8 J_azimuth = 0.;//azimuth btwn center of cone (J) and line of sight
		REAL8 J_zenith = inclination; // zenith angle btwn J and line of sight

		// Generate all available waveform modes
		SphHarmTimeSeries *ts = XLALSimInspiralChooseTDModes(
			phi, dt,
			m1, m2,
			f_min, f_ref, 
			dist,
			lambda1, lambda2,
			waveFlags,
			nonGRparams,
			amplitudeOrder, phaseOrder,
			Lmax,
			approximant
		);

		// Generate the unrotated polarizations from the modes
		ret = XLALSimInspiralPolarizationsFromSphHarmTimeSeries(&hp, &hx, ts,
				inclination, psi);
		if( ret != XLAL_SUCCESS ) XLAL_ERROR( XLAL_EFUNC );

		// Write out unrotated polarizations
		REAL8 t0 = XLALGPSGetREAL8(&(hp->epoch));
		for(i=0; i<hp->data->length; i++)
			fprintf( h_ref, "%.16g %.16g %.16g\n", t0 + i * hp->deltaT,
					hp->data->data[i], hx->data->data[i] );

		// Transform waveform so L moves on a constant precession cone
		ret = XLALSimInspiralConstantPrecessionConeWaveform(
				&hp, &hx,
				ts,
				precess_freq,
				cone_opening, cone_azimuth,
				J_azimuth, J_zenith );
        if( ret != XLAL_SUCCESS ) XLAL_ERROR( XLAL_EFUNC );

		XLALDestroySphHarmTimeSeries( ts );

		// Write out rotated waveform
		for(i=0; i<hp->data->length; i++)
			fprintf( h_rot, "%.16g %.16g %.16g\n", t0 + i * hp->deltaT,
					hp->data->data[i], hx->data->data[i] );

		// We're done.
		XLALDestroyREAL8TimeSeries( hp );
		XLALDestroyREAL8TimeSeries( hx );

		fclose(h_ref);
		fclose(h_rot);

		return 0;
}
