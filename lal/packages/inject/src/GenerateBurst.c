/*
 * Copyright (C) 2007 Jolien Creighton, Patrick Brady, Saikat Ray-Majumder,
 * Xavier Siemens, Teviet Creighton, Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <string.h>
#include <gsl/gsl_rng.h>
#include <lal/Date.h>
#include <lal/GenerateBurst.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/LALSimBurst.h>
#include <lal/LALSimulation.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/TimeSeries.h>


NRCSID(GENERATEBURSTC, "$Id$");


/*
 * ============================================================================
 *
 *                              sim_burst Nexus
 *
 * ============================================================================
 */


/*
 * Generate the + and x time series for a single sim_burst table row.
 * Note:  only the row pointed to by sim_burst is processed, the linked
 * list is not iterated over.
 */


static int XLALGenerateSimBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	const SimBurst *sim_burst,
	double delta_t
)
{
	static const char func[] = "XLALGenerateSimBurst";

	if(!strcmp(sim_burst->waveform, "BTLWNB")) {
		/* E_{GW}/r^{2} is in M_{sun} / pc^{2}, so we multiply by
		 * (M_{sun} c^2) to convert to energy/pc^{2}, and divide by
		 * (distance/pc)^{2} to convert to energy/distance^{2},
		 * which is then multiplied by (4 G / c^3) to convert to a
		 * value of \int \dot{h}^{2} \diff t.  From the values of
		 * the LAL constants, the total factor multiplying
		 * egw_over_rsquared is 1.8597e-21. */

		double int_hdot_squared_dt = sim_burst->egw_over_rsquared * LAL_MSUN_SI * 4 * LAL_G_SI / LAL_C_SI / LAL_PC_SI / LAL_PC_SI;

		/* the waveform number is interpreted as the seed for GSL's
		 * Mersenne twister random number generator */
		gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
		if(!rng)
			XLAL_ERROR(func, XLAL_ENOMEM);
		gsl_rng_set(rng, sim_burst->waveform_number);

		if(XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(hplus, hcross, sim_burst->duration, sim_burst->frequency, sim_burst->bandwidth, int_hdot_squared_dt, delta_t, rng)) {
			gsl_rng_free(rng);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		gsl_rng_free(rng);
	} else if(!strcmp(sim_burst->waveform, "StringCusp")) {
		if(XLALGenerateStringCusp(hplus, hcross, sim_burst->amplitude, sim_burst->frequency, delta_t))
			XLAL_ERROR(func, XLAL_EFUNC);
	} else if(!strcmp(sim_burst->waveform, "SineGaussian")) {
		if(XLALSimBurstSineGaussian(hplus, hcross, sim_burst->q, sim_burst->frequency, sim_burst->hrss, sim_burst->pol_ellipse_e, sim_burst->pol_ellipse_angle, delta_t))
			XLAL_ERROR(func, XLAL_EFUNC);
	} else
		/* unrecognized waveform */
		XLAL_ERROR(func, XLAL_EINVAL);

	/* done */

	return 0;
}


/*
 * Convenience wrapper to iterate over the entries in a sim_burst linked
 * list and inject them into a time series.  Passing NULL for the response
 * disables it (input time series is strain).
 */


int XLALBurstInjectSignals(
	REAL8TimeSeries *series,
	const SimBurst *sim_burst,
	const COMPLEX16FrequencySeries *response
)
{
	static const char func[] = "XLALBurstInjectSignals";
	/* to be deduced from the time series' channel name */
	LALDetector detector;
	/* + and x time series for injection waveform */
	REAL8TimeSeries *injection_hplus, *injection_hcross;
	/* injection time series as added to detector's */
	REAL8TimeSeries *injection_h;
	int i;

	/* turn the first two characters of the channel name into a
	 * detector */

	for(i = 0; i < LAL_NUM_DETECTORS; i++) {
		if(!strncmp(series->name, lalCachedDetectors[i].frDetector.prefix, 2)) {
			detector = lalCachedDetectors[i];
			break;
		}
	}
	if(i >= LAL_NUM_DETECTORS) {
		XLALPrintError("can't identify detector from channel '%s'", series->name);
		XLAL_ERROR(func, XLAL_EDATA);
	}

	/* iterate over injections */

	for(; sim_burst; sim_burst = sim_burst->next) {
		/* construct the h+ and hx time series for the injection
		 * waveform */

		if(XLALGenerateSimBurst(&injection_hplus, &injection_hcross, sim_burst, series->deltaT))
			XLAL_ERROR(func, XLAL_EFUNC);

		/* project the wave strain onto the detector's response
		 * tensor to produce the injection strain as seen in the
		 * detector. */

		injection_h = XLALSimDetectorStrainREAL8TimeSeries(injection_hplus, injection_hcross, sim_burst->ra, sim_burst->dec, sim_burst->psi, &detector, &sim_burst->time_geocent_gps);
		XLALDestroyREAL8TimeSeries(injection_hplus);
		XLALDestroyREAL8TimeSeries(injection_hcross);
		if(!injection_h)
			XLAL_ERROR(func, XLAL_EFUNC);

		/* add the injection strain time series to the detector
		 * data */

		if(XLALAddInjectionREAL8TimeSeries(series, injection_h, response)) {
			XLALDestroyREAL8TimeSeries(injection_h);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALDestroyREAL8TimeSeries(injection_h);
	}

	return 0;
}
