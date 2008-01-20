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


#include <math.h>
#include <gsl/gsl_rng.h>
#include <lal/LALSimBurst.h>
#include <lal/LALSimulation.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/GenerateBurst.h>


/* FIXME:  which of these are still needed? */
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/Inject.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>


NRCSID(GENERATEBURSTC, "$Id$");


/*
 * ============================================================================
 *
 *                              sim_burst Nexus
 *
 * ============================================================================
 */


/*
 * Convenience wrapper to iterate over the entries in a sim_burst linked
 * list and inject them into a strain time series.
 */


int XLALBurstInjectSignals(LALDetector *detector, REAL8TimeSeries *h, const SimBurstTable *sim_burst)
{
	static const char func[] = "XLALBurstInjectSignals";
	/* + and x time series for injection waveform */
	REAL8TimeSeries *injection_hplus, *injection_hcross;
	/* injection time series as added to detector's */
	REAL8TimeSeries *injection_h;

	for(; sim_burst; sim_burst = sim_burst->next) {
		/* construct the h+ and hx time series for the injection
		 * waveform */

		if(!strcmp(sim_burst->waveform, "BTLWNB")) {
			/* hrss --> int \dot{h}^2 dt, freq --> f_{0},
			 * dtplus --> duration, dtminus --> bandwidth,
			 * zm_number --> seed */
			gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
			if(!rng)
				XLAL_ERROR(func, XLAL_ENOMEM);
			gsl_rng_set(rng, sim_burst->zm_number);
			if(XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(&injection_hplus, &injection_hcross, sim_burst->dtplus, sim_burst->freq, sim_burst->dtminus, sim_burst->hrss, h->deltaT, rng)) {
				gsl_rng_free(rng);
				XLAL_ERROR(func, XLAL_EFUNC);
			}
			gsl_rng_free(rng);
		} else if(!strcmp(sim_burst->waveform, "StringCusp")) {
			/* hpeak --> amplitude, freq --> f_{high} */
			if(XLALGenerateStringCusp(&injection_hplus, &injection_hcross, sim_burst->hpeak, sim_burst->freq, h->deltaT))
				XLAL_ERROR(func, XLAL_EFUNC);
		} else
			/* unrecognized waveform */
			XLAL_ERROR(func, XLAL_EINVAL);

		/* project the wave strain onto the detector's response
		 * tensor to produce the injection strain as seen in the
		 * detector.  longitude --> right ascension, latitude -->
		 * declination, polarization --> psi, geocent_peak_time -->
		 * "time" of injection at geocentre */

		injection_h = XLALSimDetectorStrainREAL8TimeSeries(injection_hplus, injection_hcross, sim_burst->longitude, sim_burst->latitude, sim_burst->polarization, detector, &sim_burst->geocent_peak_time);
		XLALDestroyREAL8TimeSeries(injection_hplus);
		XLALDestroyREAL8TimeSeries(injection_hcross);
		if(!injection_h)
			XLAL_ERROR(func, XLAL_EFUNC);

		/* add the injection strain time series to the detector
		 * data */

		if(XLALAddInjectionREAL8TimeSeries(h, injection_h, NULL)) {
			XLALDestroyREAL8TimeSeries(injection_h);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALDestroyREAL8TimeSeries(injection_h);
	}

	return 0;
}


/*
 * ============================================================================
 *
 *                   Legacy Code --- Please Update to XLAL!
 *
 * ============================================================================
 */


static int XLALGenerateBurst(CoherentGW *output, SimBurstTable *simBurst, BurstParamStruc *params)
{
	static const char func[] = "XLALGenerateBurst";
	UINT4 n, i;		/* number of and index over samples */
	REAL8 t, dt, duration;	/* time, interval */
	REAL8 t0, tau, gtime;	/* central time, decay time, gaussian time */
	REAL8 f0;		/* initial frequency */
	REAL8 twopif0;		/* 2*pi*f0 */
	REAL4 hpeak;		/* peak strain for burst */
	REAL4 *fData;		/* pointer to frequency data */
	REAL8 *phiData;		/* pointer to phase data */
	REAL4 *aData;		/* pointer to frequency data */
	LIGOTimeGPS startTime;	/* start time of injection */

	/* Set up some constants to avoid repeated dereferencing.  notice
	 * the factor of 2 in the definition of n confusingly makes
	 * injections twice as long as the variable duration */
	duration = simBurst->dtplus + simBurst->dtminus;
	dt = params->deltaT;
	n = 2.0 * duration / dt;
	if(!n)
		XLAL_ERROR(func, XLAL_EINVAL);

	/* start time of data is peak time duration */
	startTime = simBurst->geocent_peak_time;
	XLALGPSAdd(&startTime, -duration);

	/* Generic burst parameters */
	hpeak = simBurst->hpeak;
	tau = simBurst->tau;
	f0 = simBurst->freq;
	twopif0 = LAL_TWOPI * f0;

	/* Allocate output structures. */
	output->a = LALCalloc(1, sizeof(*output->a));
	output->f = XLALCreateREAL4TimeSeries("Burst frequency", &startTime, 0.0, params->deltaT, &lalHertzUnit, n);
	output->phi = XLALCreateREAL8TimeSeries("Burst phase", &startTime, 0.0, params->deltaT, &lalDimensionlessUnit, n);
	if(!output->a || !output->f || !output->phi) {
		LALFree(output->a);
		XLALDestroyREAL4TimeSeries(output->f);
		XLALDestroyREAL8TimeSeries(output->phi);
		output->a = NULL;
		output->f = NULL;
		output->phi = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	output->a->data = XLALCreateREAL4VectorSequence(n, 2);
	if(!output->a->data) {
		LALFree(output->a);
		XLALDestroyREAL4TimeSeries(output->f);
		XLALDestroyREAL8TimeSeries(output->phi);
		output->a = NULL;
		output->f = NULL;
		output->phi = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* Set output structure metadata fields. */
	output->position.longitude = simBurst->longitude;
	output->position.latitude = simBurst->latitude;
	output->position.system = params->system;
	output->psi = simBurst->polarization;
	output->a->epoch = startTime;
	output->a->deltaT = params->deltaT;
	output->a->sampleUnits = lalStrainUnit;
	LALSnprintf(output->a->name, LALNameLength, "Burst amplitudes");

	/* Fill frequency and phase arrays. */
	fData = output->f->data->data;
	phiData = output->phi->data->data;
	aData = output->a->data->data;

	/* this depends on the waveform type */
	if(!(strcmp(simBurst->waveform, "SineGaussian"))) {
		/* find the peak time as a REAL8 relative to start of segment */
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* construct the signal */
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*(fData++) = f0;
			*(phiData++) = twopif0 * (t - t0);
			*(aData++) = hpeak * exp(-gtime * gtime);
			*(aData++) = 0.0;
		}
	} else if(!(strcmp(simBurst->waveform, "Gaussian"))) {
		/* find the peak time as a REAL8 relative to start of segment */
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* construct the signal */
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*(fData++) = 0.0;
			*(phiData++) = 0.0;
			*(aData++) = hpeak * exp(-gtime * gtime);
			*(aData++) = 0.0;
		}
	} else if(!(strcmp(simBurst->waveform, "Ringdown"))) {
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*fData++ = f0;
			*phiData++ = twopif0 * (t - t0);
			if(gtime > 0)
				*aData++ = hpeak * exp(-gtime);
			else
				*aData++ = 0;
			*aData++ = 0;
		}
	} else if(!(strcmp(simBurst->waveform, "Ringup"))) {
		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);
		for(i = 0; i < n; i++) {
			t = i * dt;
			gtime = (t - t0) / tau;
			*fData++ = f0;
			*phiData++ = twopif0 * (t - t0);
			if(gtime < 0)
				*aData++ = hpeak * exp(gtime);
			else
				*aData++ = 0;
			*aData++ = 0;
		}
	} else if(!(strcmp(simBurst->waveform, "StringCusp"))) {
		/* I use the simburst table as follows: The duration is still
		   dtplus+dtminus; hpeak is the amplitude of the cusp, not the
		   value of the strain at the peak, t0 is the central time of the
		   cusp, which introduces a phase into the waveform. The low
		   frequency cutoff will be fixed at 1Hz there's nothing special
		   about 1Hz except that it low compared to the ferquecny at which
		   we should be high-passing the data; the high frequency cutoff
		   is given by f0 */
		REAL4Vector *vector;
		COMPLEX8Vector *vtilde;
		RealFFTPlan *rplan;
		REAL4 dfreq = 1 / (2 * duration);	/* the factor of two here is becaus the length of injections 
							   is actually twice the value of the variable duration */
		REAL4 flow = 1;

		t0 = XLALGPSDiff(&(simBurst->geocent_peak_time), &startTime);

		/* create vector to store h(t), create vector that will
		 * hold frequency domain template, create fft plan */
		vector = XLALCreateREAL4Sequence(n);
		vtilde = XLALCreateCOMPLEX8Sequence(n / 2 + 1);
		rplan = XLALCreateReverseREAL4FFTPlan(n, 0);
		if(!vector || !vtilde || !rplan) {
			XLALDestroyREAL4Sequence(vector);
			XLALDestroyCOMPLEX8Sequence(vtilde);
			XLALDestroyREAL4FFTPlan(rplan);
			XLAL_ERROR(func, XLAL_EFUNC);
		}

		/* Set the FD template */
		for(i = 0; i < vtilde->length - 1; i++) {
			REAL4 freq = i * dfreq;
			vtilde->data[i].re = hpeak * pow((sqrt(1 + pow(flow, 2) * pow(freq, -2))), -8) * pow(freq, -4.0 / 3.0);

			if(freq >= f0)
				vtilde->data[i].re *= exp(1 - freq / f0);

			vtilde->data[i].im = vtilde->data[i].re * sin(-LAL_TWOPI * freq * duration);
			vtilde->data[i].re = vtilde->data[i].re * cos(-LAL_TWOPI * freq * duration);
		}

		/* set dc to zero */
		vtilde->data[0].re = 0;
		vtilde->data[0].im = 0;
		/* set nyquist to zero */
		vtilde->data[vtilde->length - 1].re = 0;
		vtilde->data[vtilde->length - 1].im = 0;

		/* Reverse FFT */
		if(XLALREAL4ReverseFFT(vector, vtilde, rplan))
			XLAL_ERROR(func, XLAL_EFUNC);

		/* multiply times dfreq to make sure units are correct */
		for(i = 0; i < vector->length; i++)
			vector->data[i] *= dfreq;

		/* make sure injection starts precisely at 0 */
		for(i = 0; i < vector->length; i++)
			vector->data[i] -= vector->data[0];

		for(i = 0; i < n; i++) {
			*fData++ = 0.0;
			*phiData++ = 0.0;
			*aData++ = vector->data[i];
			*aData++ = 0;
		}

		/* free the data */
		XLALDestroyREAL4Sequence(vector);
		XLALDestroyCOMPLEX8Sequence(vtilde);
		XLALDestroyREAL4FFTPlan(rplan);
	} else if(!(strcmp(simBurst->waveform, "warren"))) {
		/* set everything to 0 */
		for(i = 0; i < n; i++) {
			*(fData++) = 0.0;
			*(phiData++) = 0.0;
			*(aData++) = 0.0;
			*(aData++) = 0.0;
		}
	} else {
		/* unknown waveform */
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	return 0;
}


void LALBurstInjectSignals(LALStatus * stat, REAL4TimeSeries * series, SimBurstTable * injections, COMPLEX8FrequencySeries * resp, INT4 calType)
{
	UINT4 k;
	INT4 injStartTime;
	INT4 injStopTime;
	DetectorResponse detector;
	COMPLEX8Vector *unity = NULL;
	CoherentGW waveform;
	BurstParamStruc burstParam;
	REAL4TimeSeries signal;
	SimBurstTable *simBurst = NULL;
	LALDetector *tmpDetector = NULL /*,*nullDetector=NULL */ ;
	COMPLEX8FrequencySeries *transfer = NULL;

	INITSTATUS(stat, "LALBurstInjectSignals", GENERATEBURSTC);
	ATTATCHSTATUSPTR(stat);

	/* set up start and end of injection zone TODO: fix this hardwired 10 */
	injStartTime = series->epoch.gpsSeconds - 10;
	injStopTime = series->epoch.gpsSeconds + 10 + (INT4) (series->data->length * series->deltaT);

	/* 
	 *compute the transfer function 
	 */

	/* allocate memory and copy the parameters describing the freq series */
	memset(&detector, 0, sizeof(DetectorResponse));
	transfer = (COMPLEX8FrequencySeries *)
	    LALCalloc(1, sizeof(COMPLEX8FrequencySeries));
	if(!transfer) {
		ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
	}
	memcpy(&(transfer->epoch), &(resp->epoch), sizeof(LIGOTimeGPS));
	transfer->f0 = resp->f0;
	transfer->deltaF = resp->deltaF;

	tmpDetector = detector.site = (LALDetector *) LALMalloc(sizeof(LALDetector));
	/* set the detector site */
	switch (series->name[0]) {
	case 'H':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
		LALWarning(stat, "computing waveform for Hanford.");
		break;
	case 'L':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
		LALWarning(stat, "computing waveform for Livingston.");
		break;
	default:
		LALFree(detector.site);
		detector.site = NULL;
		tmpDetector = NULL;
		LALWarning(stat, "Unknown detector site, computing plus mode " "waveform with no time delay");
		break;
	}

	/* set up units for the transfer function */
	{
		RAT4 negOne = { -1, 0 };
		LALUnit unit;
		LALUnitPair pair;
		pair.unitOne = &lalADCCountUnit;
		pair.unitTwo = &lalStrainUnit;
		LALUnitRaise(stat->statusPtr, &unit, pair.unitTwo, &negOne);
		CHECKSTATUSPTR(stat);
		pair.unitTwo = &unit;
		LALUnitMultiply(stat->statusPtr, &(transfer->sampleUnits), &pair);
		CHECKSTATUSPTR(stat);
	}

	/* invert the response function to get the transfer function */
	LALCCreateVector(stat->statusPtr, &(transfer->data), resp->data->length);
	CHECKSTATUSPTR(stat);

	LALCCreateVector(stat->statusPtr, &unity, resp->data->length);
	CHECKSTATUSPTR(stat);
	for(k = 0; k < resp->data->length; ++k) {
		unity->data[k].re = 1.0;
		unity->data[k].im = 0.0;
	}

	LALCCVectorDivide(stat->statusPtr, transfer->data, unity, resp->data);
	CHECKSTATUSPTR(stat);

	LALCDestroyVector(stat->statusPtr, &unity);
	CHECKSTATUSPTR(stat);

	/* Set up a time series to hold signal in ADC counts */
	signal.deltaT = series->deltaT;
	if((signal.f0 = series->f0) != 0) {
		ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
	}
	signal.sampleUnits = lalADCCountUnit;

	signal.data = NULL;
	LALSCreateVector(stat->statusPtr, &(signal.data), series->data->length);
	CHECKSTATUSPTR(stat);

	/* loop over list of waveforms and inject into data stream */
	for(simBurst = injections; simBurst; simBurst = simBurst->next) {
		/* only do the work if the burst is in injection zone */
		if((injStartTime - simBurst->geocent_peak_time.gpsSeconds) * (injStopTime - simBurst->geocent_peak_time.gpsSeconds) > 0)
			continue;

		/* set the burt params */
		burstParam.deltaT = series->deltaT;
		if(!(strcmp(simBurst->coordinates, "HORIZON"))) {
			burstParam.system = COORDINATESYSTEM_HORIZON;
		} else if(!(strcmp(simBurst->coordinates, "ZENITH"))) {
			/* set coordinate system for completeness */
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
			detector.site = NULL;
		} else if(!(strcmp(simBurst->coordinates, "GEOGRAPHIC"))) {
			burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
		} else if(!(strcmp(simBurst->coordinates, "EQUATORIAL"))) {
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
		} else if(!(strcmp(simBurst->coordinates, "ECLIPTIC"))) {
			burstParam.system = COORDINATESYSTEM_ECLIPTIC;
		} else if(!(strcmp(simBurst->coordinates, "GALACTIC"))) {
			burstParam.system = COORDINATESYSTEM_GALACTIC;
		} else
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;

		/* generate the burst */
		memset(&waveform, 0, sizeof(CoherentGW));
		if(XLALGenerateBurst(&waveform, simBurst, &burstParam)) {
			ABORT(stat, GENERATEBURSTH_EMEM, GENERATEBURSTH_MSGEMEM);
		}

		/* must set the epoch of signal since it's used by coherent GW */
		signal.epoch = waveform.a->epoch;
		memset(signal.data->data, 0, signal.data->length * sizeof(REAL4));

		/* decide which way to calibrate the data; defaul to old way */
		if(calType)
			detector.transfer = NULL;
		else
			detector.transfer = transfer;

		/* convert this into an ADC signal */
		LALSimulateCoherentGW(stat->statusPtr, &signal, &waveform, &detector);
		CHECKSTATUSPTR(stat);

		/* if calibration using RespFilt */
		if(calType == 1)
			XLALRespFilt(&signal, transfer);

		/* inject the signal into the data channel */
		LALSSInjectTimeSeries(stat->statusPtr, series, &signal);
		CHECKSTATUSPTR(stat);

		/* free memory in coherent GW structure.  TODO:  fix this */
		LALSDestroyVectorSequence(stat->statusPtr, &(waveform.a->data));
		CHECKSTATUSPTR(stat);
		LALSDestroyVector(stat->statusPtr, &(waveform.f->data));
		CHECKSTATUSPTR(stat);
		LALDDestroyVector(stat->statusPtr, &(waveform.phi->data));
		CHECKSTATUSPTR(stat);
		LALFree(waveform.a);
		waveform.a = NULL;
		LALFree(waveform.f);
		waveform.f = NULL;
		LALFree(waveform.phi);
		waveform.phi = NULL;

		/* reset the detector site information in case it changed */
		detector.site = tmpDetector;
	}

	/* destroy the signal */
	LALSDestroyVector(stat->statusPtr, &(signal.data));
	CHECKSTATUSPTR(stat);

	LALCDestroyVector(stat->statusPtr, &(transfer->data));
	CHECKSTATUSPTR(stat);

	if(detector.site)
		LALFree(detector.site);
	LALFree(transfer);

	DETATCHSTATUSPTR(stat);
	RETURN(stat);
}
