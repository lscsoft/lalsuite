/**
 * @file LALSQTPNWaveformInterface.c
 * Contains function definitions to integrate the SpinQuadTaylor code into the other parts of the LALSuit.
 * If you want to run the program use the \ref LALSQTPNWaveformTest.c file int the
 * test directory.
 * @author László Veréb
 * @date 2010.06.27.
 */

#include <lal/LALSQTPNWaveformInterface.h>
#include <lal/LALSQTPNWaveform.h>

void LALSQTPNWaveformTemplates (LALStatus *status, REAL4Vector *signalvec1, 
		REAL4Vector *signalvec2, InspiralTemplate *params) {

	XLALPrintDeprecationWarning("LALSQTPNWaveformTemplates", 
		"XLALSQTPNWaveformTemplates");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if(XLALSQTPNWaveformTemplates(signalvec1, signalvec2, params))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNWaveformTemplates (REAL4Vector *signalvec1, 
		REAL4Vector *signalvec2, InspiralTemplate *params) {

	// Check the relevant pointers
	if( !signalvec1 || !signalvec1->data || !signalvec2 
			|| !signalvec2->data || !params )
		XLAL_ERROR(XLAL_EFAULT);
	
	// Check the parameters are sane
	if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0 			|| params->tSampling <= 0 || params->totalMass <= 0.)
		XLAL_ERROR(XLAL_EINVAL);

	InspiralInit paramsInit;
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;

	XLALInspiralInit(params, &paramsInit);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);

	memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
	memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));
	XLALSQTPNFillParams(&wave_Params, params);

	wave.waveform = NULL;
	wave.h = NULL;
	wave.hp = signalvec1;
	wave.hc = signalvec2;

	/* Call the engine function */
	if(XLALSQTPNGenerator(&wave, &wave_Params))
		XLAL_ERROR(XLAL_EFUNC);

	return XLAL_SUCCESS;
}

void LALSQTPNWaveform (LALStatus *status, REAL4Vector *signalvec, InspiralTemplate *params){

	XLALPrintDeprecationWarning("LALSQTPNWaveform", "XLALSQTPNWaveform");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if(XLALSQTPNWaveform(signalvec, params))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}

int XLALSQTPNWaveform (REAL4Vector *signalvec, InspiralTemplate *params){

	// Check the relevant pointers
	if( !signalvec || !signalvec->data || !params )
		XLAL_ERROR(XLAL_EFAULT);
	
	// Check the parameters are sane
	if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0 			|| params->tSampling <= 0 || params->totalMass <= 0.)
		XLAL_ERROR(XLAL_EINVAL);

	InspiralInit paramsInit;
	LALSQTPNWaveformParams wave_Params;
	LALSQTPNWave wave;
	memset(&wave, 0, sizeof(LALSQTPNWave));
	wave.h = NULL;
	wave.hp = signalvec;
	wave.hc = NULL;
	wave.waveform = NULL;

	XLALInspiralSetup (&(paramsInit.ak), params);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);
	if(XLALInspiralChooseModel(&(paramsInit.func),&(paramsInit.ak), params))
		XLAL_ERROR(XLAL_EFUNC);

	XLALSQTPNFillParams(&wave_Params, params);
	wave_Params.distance *= LAL_PC_SI * 1.e6;
	wave_Params.signalAmp /= LAL_PC_SI * 1.e6;

	/* Call the engine function */
	if(XLALSQTPNGenerator(&wave, &wave_Params))
		XLAL_ERROR(XLAL_EFUNC);
	params->tC = wave_Params.coalescenceTime;

	return XLAL_SUCCESS;
}

void LALSQTPNWaveformForInjection(LALStatus *status, CoherentGW *waveform,
		InspiralTemplate *params, PPNParamStruc *ppnParams) {

	XLALPrintDeprecationWarning("LALSQTPNWaveformForInjection", 
		"XLALSQTPNWaveformForInjection");
	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if(XLALSQTPNWaveformForInjection(waveform, params, ppnParams))
		ABORTXLAL(status);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}



int XLALSQTPNWaveformForInjection(CoherentGW *waveform, 
		InspiralTemplate *params, PPNParamStruc *ppnParams) {

	// Check the relevant pointers
	if( !waveform || !params || waveform->a || waveform->f
			|| waveform->phi || waveform->shift )
		XLAL_ERROR(XLAL_EFAULT);

	// variable declaration and initialization
	UINT4 i;
	InspiralInit paramsInit;

	// Compute some parameters
	XLALInspiralInit(params, &paramsInit);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);
	if (paramsInit.nbins == 0) {
		XLALPrintWarning("Warning! Waveform of zero length requested in %s. Returning empty waveform.\n ", __func__);
		return XLAL_SUCCESS;
	}

	// Allocate the waveform structures.
	xlalErrno = 0;
	XLALSQTPNAllocateCoherentGW(waveform, paramsInit.nbins);

	LALSQTPNWave wave;
	wave.h = NULL;
	wave.hp = wave.hc = NULL;
	wave.waveform = waveform;
	LALSQTPNWaveformParams wave_Params;

	// filling the parameters
	XLALSQTPNFillParams(&wave_Params, params);
	if (xlalErrno)
		XLAL_ERROR(XLAL_EFUNC);

	// calling the engine function
	if(XLALSQTPNGenerator(&wave, &wave_Params)) {
		XLALSQTPNDestroyCoherentGW(waveform);
		XLAL_ERROR(XLAL_EFUNC);
	}
	params->fFinal = wave_Params.finalFreq;
	for (i = 0; i < wave.length; i++) {
		if (waveform->phi->data->data[i] != 0.) {
			break;
		}
		if (i == wave.length - 1) {
			XLALSQTPNDestroyCoherentGW(waveform);
			XLAL_ERROR(XLAL_EFUNC);
		}
	}

	{
		if (waveform->a != NULL) {
			waveform->f->data->length = waveform->phi->data->length 				= waveform->shift->data->length = wave.length;
			waveform->a->data->length = 2 * wave.length;
			for (i = 0; i < wave.length; i++) {
				// (PPNParamStruct)ppnParams->phi === (InspiralTemplate)params->startPhase === (SimInspiralTable)injparams/this_event->coa_phase it is set to 0 in LALSQTPNWaveformTest.c at line 83.
				waveform->phi->data->data[i] 
					= waveform->phi->data->data[i] 
					- waveform->phi->data->data[wave.length-1] 
					+ ppnParams->phi;
			}
			waveform->a->deltaT = waveform->f->deltaT 
					= waveform->phi->deltaT
					= waveform->shift->deltaT 
					= 1. / params->tSampling;

			waveform->a->sampleUnits = lalStrainUnit;
			waveform->f->sampleUnits = lalHertzUnit;
			waveform->phi->sampleUnits = lalDimensionlessUnit;
			waveform->shift->sampleUnits = lalDimensionlessUnit;

			waveform->position = ppnParams->position;
			waveform->psi = ppnParams->psi;

			snprintf(waveform->a->name, LALNameLength,
					"STPN inspiral amplitudes");
			snprintf(waveform->f->name, LALNameLength,
					"STPN inspiral frequency");
			snprintf(waveform->phi->name, LALNameLength, 
					"STPN inspiral phase");
			snprintf(waveform->shift->name, LALNameLength,
					"STPN inspiral polshift");
		}
		// --- fill some output ---
		ppnParams->tc = (REAL8) (wave.length - 1) / params->tSampling;
		ppnParams->length = wave.length;
		ppnParams->dfdt 
			= ((REAL4) (waveform->f->data->data[wave.length - 1]
			- waveform->f->data->data[wave.length - 2])) 
			* ppnParams->deltaT;
		ppnParams->fStop = params->fFinal;
		ppnParams->termCode = GENERATEPPNINSPIRALH_EFSTOP;
		ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

		ppnParams->fStart = ppnParams->fStartIn;
	} // end phase condition

	return XLAL_SUCCESS;
}

int XLALSQTPNAllocateCoherentGW(CoherentGW *wave, UINT4 length) {

	if (!wave) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	if (length <= 0) {
		XLAL_ERROR(XLAL_EBADLEN);
	}
	if (wave->a || wave->f || wave->phi || wave->shift) {
		XLAL_ERROR(XLAL_EFAULT);
	}
	wave->a = (REAL4TimeVectorSeries *)LALMalloc(sizeof(REAL4TimeVectorSeries));
	wave->f = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	wave->phi = (REAL8TimeSeries *)LALMalloc(sizeof(REAL8TimeSeries));
	wave->shift = (REAL4TimeSeries *)LALMalloc(sizeof(REAL4TimeSeries));
	if (!(wave->a && wave->f && wave->phi && wave->shift)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	xlalErrno = 0;
	wave->a->data = XLALCreateREAL4VectorSequence(length, 2);
	wave->f->data = XLALCreateREAL4Vector(length);
	wave->phi->data = XLALCreateREAL8Vector(length);
	wave->shift->data = XLALCreateREAL4Vector(length);
	if (!(wave->a->data && wave->f->data && wave->phi->data && wave->shift->data)) {
		XLALSQTPNDestroyCoherentGW(wave);
		XLAL_ERROR(XLAL_ENOMEM);
	}
	return XLAL_SUCCESS;
}

void XLALSQTPNDestroyCoherentGW(CoherentGW *wave) {
	//static const char *func = "LALSQTPNDestroyCoherentGW";
	if (wave->a) {
		if (wave->a->data) {
			XLALDestroyREAL4VectorSequence(wave->a->data);
		}
		XLALFree(wave->a);

	}
	if (wave->f) {
		if (wave->f->data) {
			XLALDestroyREAL4Vector(wave->f->data);
		}
		XLALFree(wave->f);
	}
	if (wave->phi) {
		if (wave->phi->data) {
			XLALDestroyREAL8Vector(wave->phi->data);
		}
		XLALFree(wave->phi);
	}
	if (wave->shift) {
		if (wave->shift->data) {
			XLALDestroyREAL4Vector(wave->shift->data);
		}
		XLALFree(wave->shift);
	}
}

void XLALSQTPNFillParams(LALSQTPNWaveformParams *wave, InspiralTemplate *params) {
	wave->mass[0] = params->mass1;
	wave->mass[1] = params->mass2;
	wave->totalMass = wave->mass[0] + wave->mass[1];
	wave->mu = wave->mass[0] * wave->mass[1] / wave->totalMass;
	wave->eta = wave->mu / wave->totalMass;
	wave->chirpMass = wave->totalMass * pow(wave->eta, 3. / 5.);
	wave->chiAmp[0] = wave->chiAmp[1] = 0.;
	INT2 i;
	for (i = 0; i < 3; i++) {
		wave->chi[0][i] = params->spin1[i];
		wave->chi[1][i] = params->spin2[i];
		wave->chiAmp[0] += SQT_SQR(wave->chi[0][i]);
		wave->chiAmp[1] += SQT_SQR(wave->chi[1][i]);
	}
	wave->chiAmp[0] = sqrt(wave->chiAmp[0]);
	wave->chiAmp[1] = sqrt(wave->chiAmp[1]);
	for (i = 0; i < 3; i++) {
		if (wave->chiAmp[0] != 0.) {
			wave->chih[0][i] = wave->chi[0][i] / wave->chiAmp[0];
		} else {
			wave->chih[0][i] = 0.;
		}
		if (wave->chiAmp[1] != 0.) {
			wave->chih[1][i] = wave->chi[1][i] / wave->chiAmp[1];
		} else {
			wave->chih[1][i] = 0.;
		}
	}
	wave->qmParameter[0] = params->qmParameter[0];
	wave->qmParameter[1] = params->qmParameter[1];
	wave->distance = params->distance;
	wave->inclination = params->inclination;
	wave->lowerFreq = params->fLower;
	wave->finalFreq = (params->fFinal < params->fLower ? params->fCutoff : (params->fCutoff < params->fFinal ? params->fCutoff : params->fFinal));
	wave->samplingFreq = params->tSampling;
	wave->samplingTime = 1. / wave->samplingFreq;
	wave->phi = 0.;
	wave->signalAmp = 4. * wave->totalMass * wave->eta * LAL_MRSUN_SI / wave->distance;
	wave->order = params->order;
	wave->interaction = params->interaction;
	if (wave->interaction) {
		wave->interaction = (LALSimInspiralInteraction) ( wave->interaction | LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN );
	}
	/*printf("masses: %lg %lg\n", wave->mass[0], wave->mass[1]);
	printf("chis1: %lg %lg %lg\n", wave->chi[0][0], wave->chi[0][1], wave->chi[0][2]);
	printf("chis2: %lg %lg %lg\n", wave->chi[1][0], wave->chi[1][1], wave->chi[1][2]);
	printf("qmParams: %lg %lg\n", wave->qmParameter[0], wave->qmParameter[1]);
	printf("dist: %lg\n", wave->distance);
	printf("incl: %lg\n", wave->inclination);
	printf("Freq: %lg\n", wave->lowerFreq);
	printf("sF: %lg\n", wave->samplingFreq);
	printf("sT: %lg\n", wave->samplingTime);
	printf("amp: %lg\n", wave->signalAmp);
	printf("order: %d\n", wave->order);
	printf("spin: %d\n", wave->interaction);*/
}

