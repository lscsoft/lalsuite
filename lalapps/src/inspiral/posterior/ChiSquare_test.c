/* 
 *  InferenceTest.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
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

/* example command line: */
/* 
./InferenceTest --IFO [H1] --cache [/Users/john/data/triple/H1/frames.cache] --PSDstart 864162143.0 --PSDlength 1000 --srate 1024 --seglen 10 --trigtime 864162943.0
*/

#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>

//#include <mpi.h>
//#include "mpi.h"


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initVariables(LALInferenceRunState *state);

LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr, *ifoListStart;
	//ProcessParamsTable *ppt=NULL;

	
	irs = XLALCalloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " readData(): started.\n");
	irs->commandLine=commandLine;
	irs->data = LALInferenceReadData(commandLine);
	/* (this will already initialise each LALInferenceIFOData's following elements:  */
	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " readData(): finished.\n");
	if (irs->data != NULL) {
		fprintf(stdout, " initialize(): successfully read data.\n");
		
		fprintf(stdout, " injectSignal(): started.\n");
		LALInferenceInjectInspiralSignal(irs->data,commandLine);
		fprintf(stdout, " injectSignal(): finished.\n");
		
		ifoPtr = irs->data;
		ifoListStart = irs->data;
		while (ifoPtr != NULL) {
			/*If two IFOs have the same sampling rate, they should have the same timeModelh*,
			 freqModelh*, and modelParams variables to avoid excess computation 
			 in model waveform generation in the future*/
			LALInferenceIFOData * ifoPtrCompare=ifoListStart;
			int foundIFOwithSameSampleRate=0;
			while (ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
                          if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
                            ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
                            ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
                            ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;				
                            ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;				
                            ifoPtr->modelParams=ifoPtrCompare->modelParams;	
                            foundIFOwithSameSampleRate=1;	
                            break;
                          }
				ifoPtrCompare = ifoPtrCompare->next;
			}
			if(!foundIFOwithSameSampleRate){
				ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",
																	&(ifoPtr->timeData->epoch),
																	0.0,
																	ifoPtr->timeData->deltaT,
																	&lalDimensionlessUnit,
																	ifoPtr->timeData->data->length);
				ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",
																	&(ifoPtr->timeData->epoch),
																	0.0,
																	ifoPtr->timeData->deltaT,
																	&lalDimensionlessUnit,
																	ifoPtr->timeData->data->length);
				ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",
																			&(ifoPtr->freqData->epoch),
																			0.0,
																			ifoPtr->freqData->deltaF,
																			&lalDimensionlessUnit,
																			ifoPtr->freqData->data->length);
				ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",
																			 &(ifoPtr->freqData->epoch),
																			 0.0,
																			 ifoPtr->freqData->deltaF,
																			 &lalDimensionlessUnit,
																			 ifoPtr->freqData->data->length);
				ifoPtr->modelParams = XLALCalloc(1, sizeof(LALInferenceVariables));
			}
			ifoPtr = ifoPtr->next;
		}
		irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
		printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else
		fprintf(stdout, " initialize(): no data read.\n");
	
	
	return(irs);
}



/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void initVariables(LALInferenceRunState *state)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimInspiralTable *injTable=NULL;
	//LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	ProcessParamsTable *ppt=NULL;
	//INT4 AmpOrder=0;
	LALPNOrder PhaseOrder=LAL_PNORDER_TWO;
	Approximant approx=TaylorF2;
	//INT4 numberI4 = TaylorF2;
	//INT4 numberI4 = TaylorT3;
	//INT4 approx=TaylorF2;
	LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
	//REAL8 logDmin=log(1.0);
	//REAL8 logDmax=log(100.0);
	//REAL8 Dmin=1.0;
	//REAL8 Dmax=100.0;
	//REAL8 mcMin=1.0;
	//REAL8 mcMax=15.3;
	//REAL8 logmcMax,logmcMin;
	//REAL8 mMin=1.0,mMax=30.0;
	//REAL8 MTotMax=35.0;
	REAL8 etaMin=0.0312;
	REAL8 etaMax=0.25;
	//REAL8 dt=0.1;            /* Width of time prior */
	//REAL8 tmpMin,tmpMax,tmpVal;
  
  FILE *devrandom;
	struct timeval tv;
  unsigned int randomseed=0;
  
  gsl_rng_env_setup();
	state->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
  
		if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
				gettimeofday(&tv, 0);
				randomseed = tv.tv_sec + tv.tv_usec;
		} 
		else {
				size_t fread_ret = fread(&randomseed, sizeof(randomseed), 1, devrandom);
				if (fread_ret != 1) {
					fprintf(stdout, "could not read a random seed from /dev/urandom, using time of day.");
					gettimeofday(&tv, 0);
					randomseed = tv.tv_sec + tv.tv_usec;
				}
				fclose(devrandom);
		}
	

	fprintf(stdout, " initialize(): random seed: %u\n", randomseed);
	//LALInferenceAddVariable(state->algorithmParams,"random_seed",&randomseed, UINT4_t,LALINFERENCE_PARAM_FIXED);
	gsl_rng_set(state->GSLrandom, randomseed);
	gsl_rng * GSLrandom=state->GSLrandom;
	REAL8 endtime=0.0, timeParam=0.0;
	REAL8 start_mc			=4.82+gsl_ran_gaussian(GSLrandom,0.025);
	REAL8 start_eta			=etaMin+gsl_rng_uniform(GSLrandom)*(etaMax-etaMin);
	REAL8 start_phase		=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
	REAL8 start_dist		=8.07955+gsl_ran_gaussian(GSLrandom,1.1);
	REAL8 start_ra			=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
	REAL8 start_dec			=-LAL_PI/2.0+gsl_rng_uniform(GSLrandom)*(LAL_PI_2-(-LAL_PI_2));
	REAL8 start_psi			=0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
	REAL8 start_iota		=0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
	REAL8 start_a_spin1		=0.0+gsl_rng_uniform(GSLrandom)*(1.0-0.0);
	REAL8 start_theta_spin1 =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
	REAL8 start_phi_spin1	=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
	REAL8 start_a_spin2		=0.0+gsl_rng_uniform(GSLrandom)*(1.0-0.0);
	REAL8 start_theta_spin2 =0.0+gsl_rng_uniform(GSLrandom)*(LAL_PI-0.0);
	REAL8 start_phi_spin2	=0.0+gsl_rng_uniform(GSLrandom)*(LAL_TWOPI-0.0);
        INT4  numBins           = 16;
	
	memset(currentParams,0,sizeof(LALInferenceVariables));
	
	char help[]="\
	[--inj injections.xml]\tInjection XML file to use\
	[--Mmin mchirp]\tMinimum chirp mass\
	[--Mmax mchirp]\tMaximum chirp mass\
	[--dt time]\tWidth of time prior, centred around trigger (0.1s)\
	[--trigtime time]\tTrigger time to use\
	[--mc mchirp]\tTrigger chirpmass to use\
	[--eta eta]\tTrigger eta to use\
	[--phi phase]\tTrigger phase to use\
	[--iota inclination]\tTrigger inclination to use\
        [--dist dist]\tTrigger distance\
        [--ra ra]\tTrigger RA\
        [--dec dec]\tTrigger declination\
        [--psi psi]\tTrigger psi\
        [--a1 a1]\tTrigger a1\
        [--theta1 theta1]\tTrigger theta1\
        [--phi1 phi1]\tTrigger phi1\
        [--a2 a2]\tTrigger a2\
        [--theta2 theta2]\tTrigger theta2\
        [--phi2 phi2]\tTrigger phi2\
        [--time time]\tWaveform time (overrides random about trigtime)\
        [--Dmin dist]\tMinimum distance in Mpc (1)\
        [--Dmax dist]\tMaximum distance in Mpc (100)\
        [--approx ApproximantorderPN]\tSpecify a waveform to use, (default TaylorF2twoPN)\
        [--mincomp min]\tMinimum component mass (1.0)\
        [--maxcomp max]\tMaximum component mass (30.0)\
        [--MTotMax] \t Maximum total mass (35.0)\
        [--covarianceMatrix file]\tFind the Cholesky decomposition of the covariance matrix for jumps in file\
        [--num-bins\tNumber of frequency bins to use to compute Chisq statistics, (default 16)";
	
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}
	
	/* Read injection XML file for parameters if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
	if(ppt){
		SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
		if(!injTable){
			fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
			exit(1);
		}
		endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
		//AmpOrder=injTable->amp_order;
		PhaseOrder = XLALGetOrderFromString(injTable->waveform);
		if( (int) PhaseOrder == XLAL_FAILURE)
			ABORTXLAL(&status);
		approx = XLALGetApproximantFromString(injTable->waveform);
		if( (int) approx == XLAL_FAILURE)
			ABORTXLAL(&status);
	}	
	
	/* Over-ride approximant if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
	if(ppt){
		PhaseOrder = XLALGetOrderFromString(ppt->value);
		if( (int) PhaseOrder == XLAL_FAILURE)
			ABORTXLAL(&status);
		approx = XLALGetApproximantFromString(ppt->value);
		if( (int) approx == XLAL_FAILURE)
			ABORTXLAL(&status);
		//printf("%d\n",approx);
		if(strstr(ppt->value,"TaylorF2")) {approx=TaylorF2;}//numberI4 = TaylorF2;}		LALGetApproximantFromString DOES NOT HAVE TaylorF2 !!!!!!
		//if(strstr(ppt->value,"TaylorT3")) {approx=TaylorT3;}//numberI4 = TaylorT3;}
		//if(strstr(ppt->value,"SpinTaylor")) {approx=SpinTaylor;}//numberI4 = SpinTaylor;}
		fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",approx,PhaseOrder);
		//fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",numberI4,PhaseOrder);
	}
  
	/* Over-ride taper if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
	if(ppt){
		if(strstr(ppt->value,"STARTEND")) bookends=LALINFERENCE_TAPER_STARTEND;
		if(strstr(ppt->value,"STARTONLY")) bookends=LALINFERENCE_TAPER_START;
		if(strstr(ppt->value,"ENDONLY")) bookends=LALINFERENCE_TAPER_END;
		if(strstr(ppt->value,"RING")) bookends=LALINFERENCE_RING;
		if(strstr(ppt->value,"SMOOTH")) bookends=LALINFERENCE_SMOOTH;
	}
  
	/* Over-ride end time if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
	if(ppt){
		endtime=atof(ppt->value);
	}
  
	/* Over-ride chirp mass if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--mc");
	if(ppt){
		start_mc=atof(ppt->value);
	}
	
	/* Over-ride eta if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--eta");
	if(ppt){
		start_eta=atof(ppt->value);
	}
	
	/* Over-ride phase if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--phi");
	if(ppt){
		start_phase=atof(ppt->value);
	}
	
	/* Over-ride inclination if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--iota");
	if(ppt){
		start_iota=atof(ppt->value);
	}
  
  /* Over-ride distance if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dist");
  if (ppt) {
    start_dist = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--ra");
  if (ppt) {
    start_ra = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--dec");
  if (ppt) {
    start_dec = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--psi");
  if (ppt) {
    start_psi = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--a1");
  if (ppt) {
    start_a_spin1 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--theta1");
  if (ppt) {
    start_theta_spin1 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--phi1");
  if (ppt) {
    start_phi_spin1 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--a2");
  if (ppt) {
    start_a_spin2 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--theta2");
  if (ppt) {
    start_theta_spin2 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--phi2");
  if (ppt) {
    start_phi_spin2 = atof(ppt->value);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--num-bins");
  if (ppt) {
    numBins = atof(ppt->value);
  }
	
	/* Over-ride time prior if specified */
//	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
//	if(ppt){
//		dt=atof(ppt->value);
//	}
	
	/* Over-ride Distance min if specified */
//	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
//	if(ppt){
//		logDmin=log(atof(ppt->value));
//		Dmin=atof(ppt->value);
//	}
	
	/* Over-ride Distance max if specified */
//	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
//	if(ppt){
//		logDmax=log(atof(ppt->value));
//		Dmax=atof(ppt->value);
//	}
	
	/* Over-ride Mass prior if specified */
//	ppt=LALInferenceGetProcParamVal(commandLine,"--Mmin");
//	if(ppt){
//		mcMin=atof(ppt->value);
//	}
//	ppt=LALInferenceGetProcParamVal(commandLine,"--Mmax");
//	if(ppt)	mcMax=atof(ppt->value);
	
	/* Over-ride component masses */
//	ppt=LALInferenceGetProcParamVal(commandLine,"--compmin");
//	if(ppt)	mMin=atof(ppt->value);
//	LALInferenceAddVariable(priorArgs,"component_min",&mMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
//	ppt=LALInferenceGetProcParamVal(commandLine,"--compmax");
//	if(ppt)	mMax=atof(ppt->value);
//	LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
//	ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMax");
//	if(ppt)	MTotMax=atof(ppt->value);
//	LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	
	printf("Read end time %f\n",endtime);
  
        LALInferenceAddVariable(currentParams, "numbins",       &numBins,             LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	//LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	//numberI4 = LAL_PNORDER_TWO;
  LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);	
	//LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &numberI4,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
	if(ppt){
		LALInferenceAddVariable(currentParams, "INFERENCE_TAPER",     &bookends,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	}
  ppt=LALInferenceGetProcParamVal(commandLine,"--newswitch");
  int newswitch=0;
  if(ppt){
    newswitch=1;
    LALInferenceAddVariable(currentParams, "newswitch", &newswitch, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }
	/* Set up the variable parameters */
	//tmpVal=4.82+gsl_ran_gaussian(GSLrandom,0.025);//log(mcMin+(mcMax-mcMin)/2.0);
	//tmpVal=7.86508;
	LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
	//LALInferenceAddVariable(currentParams, "chirpmass",    &tmpVal,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);
	//LALInferenceAddVariable(currentParams,"logmc",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	//logmcMin=log(mcMin); logmcMax=log(mcMax);
	//addMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		LALINFERENCE_REAL8_t);
  
	//tmpVal=0.244;
	//tmpVal=0.03+gsl_rng_uniform(GSLrandom)*(0.25-0.03);
	//tmpVal=0.18957;
	LALInferenceAddVariable(currentParams, "massratio",       &start_eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	//LALInferenceAddVariable(currentParams, "massratio",       &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  //addMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	LALINFERENCE_REAL8_t);
	
	//tmpMin=endtime-dt; tmpMax=endtime+dt;
  
  /* Set up start time. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--time");
  if (ppt) {
    /* User has specified start time. */
    timeParam = atof(ppt->value);
  } else {
    timeParam = endtime+gsl_ran_gaussian(GSLrandom,0.01);          
  }
  
  LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	//addMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	
  
	//tmpVal=1.5;
	//tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.89954;
  LALInferenceAddVariable(currentParams, "phase",           &start_phase,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	//LALInferenceAddVariable(currentParams, "phase",           &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	//tmpVal=5.8287;
	//tmpVal=8.07955+gsl_ran_gaussian(GSLrandom,1.1);
	//Dmin+(Dmax-Dmin)/2.0;
	//tmpVal=46.92314;
	LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	//LALInferenceAddVariable(currentParams,"distance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);
  
	
	//tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=4.5500;//1.0;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.34650;
	LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	//LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	//tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	//tmpVal=1.0759;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=-0.90547;
	LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	//LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
  
	//tmpMin=0.0; tmpMax=LAL_PI;
	//tmpVal=0.2000;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=0.64546;
	LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	//LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	//tmpMin=0.0; tmpMax=LAL_PI;
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--max-iota");
  //if (ppt) {
  //  tmpMax = atof(ppt->value);
  //}
	//tmpVal=0.9207;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=2.86094;
 	LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	//LALInferenceAddVariable(currentParams, "inclination",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
	//addMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
	if((approx==SpinTaylor ||approx==SpinTaylorT3) && !ppt){
		
		ppt=LALInferenceGetProcParamVal(commandLine, "--singleSpin");
		if(ppt) fprintf(stdout,"Running with first spin set to 0\n");
		else {
      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      //if(ppt) tmpMin=-1.0;
      //else tmpMin=0.0;
      //tmpMax=1.0;
			LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			//addMinMaxPrior(priorArgs, "a_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
      
			ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
			if(ppt) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
			else {
				//tmpMin=0.0; tmpMax=LAL_PI;
				LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
				//addMinMaxPrior(priorArgs, "theta_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
        
				//tmpMin=0.0; tmpMax=LAL_TWOPI;
				LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
				//addMinMaxPrior(priorArgs, "phi_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
			}
		}
		
    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    //if(ppt) tmpMin=-1.0;
    //else tmpMin=0.0;
    //tmpMax=1.0;
		LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		//addMinMaxPrior(priorArgs, "a_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
		ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
		if(ppt) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
		else {
			//tmpMin=0.0; tmpMax=LAL_PI;
			LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			//addMinMaxPrior(priorArgs, "theta_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
      
			//tmpMin=0.0; tmpMax=LAL_TWOPI;
			LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			//addMinMaxPrior(priorArgs, "phi_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		}
	}
	
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
	if(approx==TaylorF2 && ppt){
		
    //tmpMin=-1.0; tmpMax=1.0;
    LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    //addMinMaxPrior(priorArgs, "spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		
		//tmpMin=-1.0; tmpMax=1.0;
		LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		//addMinMaxPrior(priorArgs, "spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
	}
  
  
    state->differentialPoints = NULL;
    state->differentialPointsLength = 0;
  
  
	
	return;
}



int main(int argc, char *argv[]){


	LALInferenceRunState *runState;
	ProcessParamsTable *procParams=NULL;
	
	/* Read command line and parse */
	procParams=LALInferenceParseCommandLine(argc,argv);
	
	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	runState = initialize(procParams);
 	initVariables(runState);
	
    REAL8 loglikelihood=0.0;
	loglikelihood = LALInferenceFreqDomainLogLikelihood(runState->currentParams, runState->data, &LALInferenceTemplateLALGenerateInspiral) - LALInferenceNullLogLikelihood(runState->data);
	printf("network LogLikelihood=%f\tSNR=%f\n",loglikelihood,sqrt(2*loglikelihood));
    LALInferenceChiSquareTest(runState->currentParams, runState->data, &LALInferenceTemplateLALGenerateInspiral);

  return 0;
}

