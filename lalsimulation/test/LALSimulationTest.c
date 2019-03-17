#include <math.h>
#include <stdio.h>

#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/TimeDelay.h>
#include <lal/LALSimulation.h>
#include <gsl/gsl_sf_trig.h>

static LIGOTimeGPS gps_zero = LIGOTIMEGPSZERO;

static REAL8TimeSeries *new_series(double deltaT, unsigned length, double fill)
{
	REAL8TimeSeries *new = XLALCreateREAL8TimeSeries("blah", &gps_zero, 0.0, deltaT, &lalDimensionlessUnit, length);
	unsigned i;
	for(i = 0 ; i < new->data->length; i++)
		new->data->data[i] = fill;
	return new;
}

static REAL8TimeSeries *copy_series(const REAL8TimeSeries *src)
{
	return XLALCutREAL8TimeSeries(src, 0, src->data->length);
}

static void add_circular_polarized_sine(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross,  LIGOTimeGPS epoch, double ampl, double freq)
{
	double t0 = XLALGPSDiff(&hplus->epoch, &epoch);
	unsigned i;

	for(i = 0; i < hplus->data->length; i++) {
		double t = t0 + i * hplus->deltaT;
		hplus->data->data[i] = ampl * cos(LAL_TWOPI * freq * t);
		hcross->data->data[i] = ampl * sin(LAL_TWOPI * freq * t);
	}
}

static double modelsignal(LIGOTimeGPS t, LIGOTimeGPS epoch, double ampl, double freq, double armlen, double xcos, double ycos, double fxplus, double fxcross, double fyplus, double fycross)
{
	double t0 = XLALGPSDiff(&t, &epoch);

	double Txplus = armlen * (1. + xcos) / LAL_C_SI;
	double Txminus = armlen * (1. - xcos) / LAL_C_SI;
	double Typlus = armlen * (1. + ycos) / LAL_C_SI;
	double Tyminus = armlen * (1. - ycos) / LAL_C_SI;

	double t1x = t0 - Txplus / 2.0;
	double t2x = t0 + Txminus / 2.0;
	double signalxplus = ampl * fxplus * (gsl_sf_sinc(freq * Txminus) * cos(LAL_TWOPI * freq * t1x) + gsl_sf_sinc(freq * Txplus) * cos(LAL_TWOPI * freq * t2x)) / 2.0;
	double signalxcross = ampl * fxcross * (gsl_sf_sinc(freq * Txminus) * sin(LAL_TWOPI * freq * t1x) + gsl_sf_sinc(freq * Txplus) * sin(LAL_TWOPI * freq * t2x)) / 2.0;

	double t1y = t0 - Typlus / 2.0;
	double t2y = t0 + Tyminus / 2.0;
	double signalyplus = ampl * fyplus * (gsl_sf_sinc(freq * Tyminus) * cos(LAL_TWOPI * freq * t1y) + gsl_sf_sinc(freq * Typlus) * cos(LAL_TWOPI * freq * t2y)) / 2.0;
	double signalycross = ampl * fycross * (gsl_sf_sinc(freq * Tyminus) * sin(LAL_TWOPI * freq * t1y) + gsl_sf_sinc(freq * Typlus) * sin(LAL_TWOPI * freq * t2y)) / 2.0;

	return signalxplus + signalxcross + signalyplus + signalycross;
}

/* add the signal the detector observes when hplus = ampl * cos(2 * pi * f * (t - epoch)) and hcross = ampl * sin(2 * pi * f * (t - epoch)) */
static void compute_answer(REAL8TimeSeries *mdl, LIGOTimeGPS epoch, double ampl, double freq, REAL8 right_ascension, REAL8 declination, REAL8 psi, const LALDetector *detector)
{
	double armlen = XLAL_REAL8_FAIL_NAN;
	double xcos = XLAL_REAL8_FAIL_NAN;
	double ycos = XLAL_REAL8_FAIL_NAN;
	double fxplus = XLAL_REAL8_FAIL_NAN;
	double fxcross = XLAL_REAL8_FAIL_NAN;
	double fyplus = XLAL_REAL8_FAIL_NAN;
	double fycross = XLAL_REAL8_FAIL_NAN;
	unsigned i;

	for(i = 0; i < mdl->data->length; i++) {
		LIGOTimeGPS t = mdl->epoch;
		XLALGPSAdd(&t, i * mdl->deltaT);

		XLALComputeDetAMResponseParts(&armlen, &xcos, &ycos, &fxplus, &fyplus, &fxcross, &fycross, detector, right_ascension, declination, psi, XLALGreenwichMeanSiderealTime(&t));

		XLALGPSAdd(&t, -XLALTimeDelayFromEarthCenter(detector->location, right_ascension, declination, &t));

		mdl->data->data[i] = modelsignal(t, epoch, ampl, freq, armlen, xcos, ycos, fxplus, fxcross, fyplus, fycross);
	}
}

static REAL8TimeSeries *error(const REAL8TimeSeries *s1, const REAL8TimeSeries *s0)
{
	REAL8TimeSeries *result = copy_series(s1);
	unsigned i;

	for(i = 0; i < s1->data->length; i++)
		result->data->data[i] -= s0->data->data[i];

	return result;
}


static double RMS(const REAL8TimeSeries *s)
{
	double rms = 0.;
	unsigned i;

	for(i = 0; i < s->data->length; i++)
		rms += s->data->data[i] * s->data->data[i];

	return sqrt(rms / s->data->length);
}


static void minmax(const REAL8TimeSeries *s, double *min, double *max)
{
	unsigned i;

	*min = *max = s->data->data[0];
	for(i = 1; i < s->data->length; i++) {
		if(s->data->data[i] < *min)
			*min = s->data->data[i];
		if(s->data->data[i] > *max)
			*max = s->data->data[i];
	}
}


static void check_result(const REAL8TimeSeries *model, const REAL8TimeSeries *result, double rms_bound, double residual_min, double residual_max)
{
	REAL8TimeSeries *err = error(model, result);
	double rms, min, max;

	rms = RMS(err);
	minmax(err, &min, &max);

	fprintf(stderr, "error vector: RMS=%g, min=%g, max=%g\n", rms, min, max);
	XLALDestroyREAL8TimeSeries(err);
	if(rms > rms_bound || min < residual_min || max > residual_max) {
		fprintf(stderr, "error vector larger than allowed\n");
		exit(1);
	}
}


int main(void)
{
	REAL8TimeSeries *hplus, *hcross, *dst, *short_dst, *mdl;
	LALDetector detector = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	double f, ampl, dt;
	unsigned length_origin, start_mdl, length_mdl;
	REAL8 right_ascension = 0.0, declination = 0.0, psi = 0.0;

	ampl = 1.0;
	f = 100.0;
	dt = 1.0 / (f * 4.0);
	length_origin = 1024 * 3;
	start_mdl = 1024;
	length_mdl = 1024;

	hplus = new_series(dt, length_origin, 0.0);
	hcross = copy_series(hplus);

	add_circular_polarized_sine(hplus, hcross, hplus->epoch, ampl, f);
	dst = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, right_ascension, declination, psi, &detector);
	short_dst = XLALCutREAL8TimeSeries(dst, start_mdl, length_mdl);

	fprintf(stderr, "injecting unit amplitude %g Hz circular polarized monochromatic GWs sampled at %g Hz into LHO data\n", f, 1/ dt);

	mdl = copy_series(short_dst);
	compute_answer(mdl, hplus->epoch,  ampl, f, right_ascension, declination, psi, &detector);

	check_result(mdl, short_dst, 0.0012, -0.0016, 0.0016);

	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(dst);
	XLALDestroyREAL8TimeSeries(short_dst);
	XLALDestroyREAL8TimeSeries(mdl);

	detector = lalCachedDetectors[LAL_ET1_DETECTOR];
	f = 10000.0;
	dt = 1.0 / (f * 4.0);
	length_origin = 1024 * 3;
	start_mdl = 1024;
	length_mdl = 1024;

	hplus = new_series(dt, length_origin, 0.0);
	hcross = copy_series(hplus);

	add_circular_polarized_sine(hplus, hcross, hplus->epoch, ampl, f);
	dst = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, right_ascension, declination, psi, &detector);
	short_dst = XLALCutREAL8TimeSeries(dst, start_mdl, length_mdl);

	fprintf(stderr, "injecting unit amplitude %g Hz circular polarized monochromatic GWs sampled at %g Hz into ET data\n", f, 1/ dt);

	mdl = copy_series(short_dst);
	compute_answer(mdl, hplus->epoch,  ampl, f, right_ascension, declination, psi, &detector);

	check_result(mdl, short_dst, 0.0004, -0.0006, 0.0006);

	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(dst);
	XLALDestroyREAL8TimeSeries(short_dst);
	XLALDestroyREAL8TimeSeries(mdl);

	exit(0);
}
