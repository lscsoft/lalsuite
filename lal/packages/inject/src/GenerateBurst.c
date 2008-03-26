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


int XLALGenerateSimBurst(
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

		XLALPrintInfo("%s(): BTLWNB: f = %.16g Hz, df = %.16g Hz, dt = %.16g s, hdot^2 = %.16g\n", func, sim_burst->frequency, sim_burst->bandwidth, sim_burst->duration, int_hdot_squared_dt);
		if(XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(hplus, hcross, sim_burst->duration, sim_burst->frequency, sim_burst->bandwidth, int_hdot_squared_dt, delta_t, rng)) {
			gsl_rng_free(rng);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		gsl_rng_free(rng);
	} else if(!strcmp(sim_burst->waveform, "StringCusp")) {
		XLALPrintInfo("%s(): string cusp: A = %.16g, fhigh = %.16g Hz\n", func, sim_burst->amplitude, sim_burst->frequency);
		if(XLALGenerateStringCusp(hplus, hcross, sim_burst->amplitude, sim_burst->frequency, delta_t))
			XLAL_ERROR(func, XLAL_EFUNC);
	} else if(!strcmp(sim_burst->waveform, "SineGaussian")) {
		XLALPrintInfo("%s(): sine-Gaussian: f = %.16g Hz, Q = %.16g, hrss = %.16g\n", func, sim_burst->frequency, sim_burst->q, sim_burst->hrss);
		if(XLALSimBurstSineGaussian(hplus, hcross, sim_burst->q, sim_burst->frequency, sim_burst->hrss, sim_burst->pol_ellipse_e, sim_burst->pol_ellipse_angle, delta_t))
			XLAL_ERROR(func, XLAL_EFUNC);
	} else {
		/* unrecognized waveform */
		XLALPrintError("%s(): error: unrecognized waveform\n");
		XLAL_ERROR(func, XLAL_EINVAL);
	}

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
	REAL8TimeSeries *hplus, *hcross;
	/* injection time series as added to detector's */
	REAL8TimeSeries *h;
	/* skip injections whose geocentre times are more than this many
	 * seconds outside of the target time series */
	const double injection_window = 600.0;
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
		XLALPrintError("%s(): error: can't identify instrument from channel '%s'\n", func, series->name);
		XLAL_ERROR(func, XLAL_EDATA);
	}
	XLALPrintInfo("%s(): channel name is '%s', instrument appears to be '%s'\n", func, series->name, detector.frDetector.prefix);

	/* iterate over injections */

	for(; sim_burst; sim_burst = sim_burst->next) {
		/* skip injections whose "times" are too far outside of the
		 * target time series */

		if(XLALGPSDiff(&series->epoch, &sim_burst->time_geocent_gps) > injection_window || XLALGPSDiff(&sim_burst->time_geocent_gps, &series->epoch) > (series->data->length * series->deltaT + injection_window))
			continue;

		/* construct the h+ and hx time series for the injection
		 * waveform.  in the time series produced by this function,
		 * t = 0 */

		if(XLALGenerateSimBurst(&hplus, &hcross, sim_burst, series->deltaT))
			XLAL_ERROR(func, XLAL_EFUNC);

		/* add the time of the injection at the geocentre to the
		 * start times of the h+ and hx time series.  after this,
		 * their epochs mark the start of those time series at the
		 * geocentre. */

		XLALGPSAddGPS(&hcross->epoch, &sim_burst->time_geocent_gps);
		XLALGPSAddGPS(&hplus->epoch, &sim_burst->time_geocent_gps);

		/* project the wave onto the detector to produce the strain
		 * in the detector. */

		h = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, sim_burst->ra, sim_burst->dec, sim_burst->psi, &detector);
		XLALDestroyREAL8TimeSeries(hplus);
		XLALDestroyREAL8TimeSeries(hcross);
		if(!h)
			XLAL_ERROR(func, XLAL_EFUNC);

		/* add the injection strain time series to the detector
		 * data */

		if(XLALSimAddInjectionREAL8TimeSeries(series, h, response)) {
			XLALDestroyREAL8TimeSeries(h);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALDestroyREAL8TimeSeries(h);
	}

	/* done */

	return 0;
}
