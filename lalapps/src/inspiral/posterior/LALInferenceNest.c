/* 
 *  InferenceNest.c:  Nested Sampling using LALInference
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


#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include "LALInference.h"
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include "LALInferenceNestedSampler.h"
#include "LALInferencePrior.h"

void initialiseNS(LALInferenceRunState *state);
LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeNS(LALInferenceRunState *runState);
void setupLivePointsArray(LALInferenceRunState *runState);
void initVariables(LALInferenceRunState *state);


LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	LALInferenceRunState *irs=NULL;
	LALIFOData *ifoPtr, *ifoListStart;
	ProcessParamsTable *ppt=NULL;
	unsigned long int randomseed;
	struct timeval tv;
	FILE *devrandom;
	
	irs = calloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " readData(): started.\n");
	irs->commandLine=commandLine;
	irs->data = readData(commandLine);
	/* (this will already initialise each LALIFOData's following elements:  */
	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " readData(): finished.\n");
	if (irs->data != NULL) {
		fprintf(stdout, " initialize(): successfully read data.\n");
		
		fprintf(stdout, " injectSignal(): started.\n");
		injectSignal(irs->data,commandLine);
		fprintf(stdout, " injectSignal(): finished.\n");
		
		ifoPtr = irs->data;
		ifoListStart = irs->data;
		while (ifoPtr != NULL) {
			/*If two IFOs have the same sampling rate, they should have the same timeModelh*,
			 freqModelh*, and modelParams variables to avoid excess computation 
			 in model waveform generation in the future*/
			LALIFOData * ifoPtrCompare=ifoListStart;
			int foundIFOwithSameSampleRate=0;
			while(ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
				if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
					ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
					ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
					ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;				
					ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;				
					ifoPtr->modelParams=ifoPtrCompare->modelParams;	
					foundIFOwithSameSampleRate=1;	
					break;
				}
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
				ifoPtr->modelParams = calloc(1, sizeof(LALVariables));
			}
			ifoPtr = ifoPtr->next;
		}
		irs->currentLikelihood=NullLogLikelihood(irs->data);
		printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else
		fprintf(stdout, " initialize(): no data read.\n");
	
	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	irs->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = getProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	else { /* otherwise generate "random" random seed: */
		if ((devrandom = fopen("/dev/random","r")) == NULL) {
			gettimeofday(&tv, 0);
			randomseed = tv.tv_sec + tv.tv_usec;
		} 
		else {
			if(1!=fread(&randomseed, sizeof(randomseed), 1, devrandom)){
			  fprintf(stderr,"Error: Unable to read random seed from /dev/random\n");
			  exit(1);
			}
			fclose(devrandom);
		}
	}
	fprintf(stdout, " initialize(): random seed: %lu\n", randomseed);
	gsl_rng_set(irs->GSLrandom, randomseed);
	
	return(irs);
}


/***** Initialise Nested Sampling structures ****/
/* Fill in samples from the prior distribution */
/* runState->algorithmParams must contain a variable "logLikelihoods" */
/* which contains a REAL8 array of likelihood values for the live */
/* points. */
/************************************************/
void initializeNS(LALInferenceRunState *runState)
{
	char help[]="\
	--Nlive N\tNumber of live points to use\n\
	--Nmcmc M\tNumber of MCMC point to use when evolving live points\n\
	[--Nruns R]\tNumber of parallel samples from logt to use(1)\n\
	[--tolerance dZ]\tTolerance of nested sampling algorithm (0.1)\n\
	[--randomseed seed]\tRandom seed of sampling distribution\n\
	[--verbose]\n";
	
	INT4 verbose=0,tmpi=0,randomseed=0;
	REAL8 tmp=0;
	ProcessParamsTable *commandLine=runState->commandLine;
	ProcessParamsTable *ppt=NULL;

	/* Print command line arguments if help requested */
	ppt=getProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}
	
	/* Initialise parameters structure */
	runState->algorithmParams=XLALCalloc(1,sizeof(LALVariables));
	runState->priorArgs=XLALCalloc(1,sizeof(LALVariables));
	runState->proposalArgs=XLALCalloc(1,sizeof(LALVariables));
	
	/* Set up the appropriate functions for the nested sampling algorithm */
	runState->algorithm=&NestedSamplingAlgorithm;
	runState->evolve=&NestedSamplingOneStep;
	runState->proposal=&LALInferenceProposalNS;

	/* This is the LAL template generator for inspiral signals */
	runState->template=&templateLAL;
	runState->likelihood=&FreqDomainLogLikelihood;
	runState->prior = &LALInferenceInspiralPrior;
	
	ppt=getProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		addVariable(runState->algorithmParams,"verbose", &verbose , INT4_t,
					PARAM_FIXED);
	}
		
	printf("set number of live points.\n");
	/* Number of live points */
	ppt=getProcParamVal(commandLine,"--Nlive");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		fprintf(stderr,"Error, must specify number of live points\n");
		exit(1);
	}
	addVariable(runState->algorithmParams,"Nlive",&tmpi, INT4_t,PARAM_FIXED);
	
	printf("set number of MCMC points.\n");
	/* Number of points in MCMC chain */
	ppt=getProcParamVal(commandLine,"--Nmcmc");
	if(ppt)
	  tmpi=atoi(ppt->value);
	else {
	  fprintf(stderr,"Error, must specify number of MCMC points\n");
	  exit(1);
	}
	addVariable(runState->algorithmParams,"Nmcmc",&tmpi,
				INT4_t,PARAM_FIXED);
	
	printf("set number of parallel runs.\n");
	/* Optionally specify number of parallel runs */
	ppt=getProcParamVal(commandLine,"--Nruns");
	if(ppt) {
		tmpi=atoi(ppt->value);
		addVariable(runState->algorithmParams,"Nruns",&tmpi,INT4_t,PARAM_FIXED);
	}
	
	printf("set tolerance.\n");
	/* Tolerance of the Nested sampling integrator */
	ppt=getProcParamVal(commandLine,"--tolerance");
	if(ppt){
		tmp=strtod(ppt->value,(char **)NULL);
		addVariable(runState->algorithmParams,"tolerance",&tmp, REAL8_t,
					PARAM_FIXED);
	}
	
	printf("set random seed.\n");
	/* Set up the random number generator */
	gsl_rng_env_setup();
	runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	
	/* (try to) get random seed from command line: */
	ppt = getProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	fprintf(stdout, " initialize(): random seed: %u\n", randomseed);
	addVariable(runState->algorithmParams,"random_seed",&randomseed, INT4_t,PARAM_FIXED);
	gsl_rng_set(runState->GSLrandom, randomseed);
	
	return;
	
}

void setupLivePointsArray(LALInferenceRunState *runState){
	/* Set up initial basket of live points, drawn from prior,
	 by copying runState->currentParams to all entries in the array*/
	
	UINT4 Nlive=(UINT4)*(INT4 *)getVariable(runState->algorithmParams,"Nlive");
	UINT4 i;
	REAL8Vector *logLs;
	
	LALVariableItem *current;
		
	/* Allocate the array */
	/* runState->livePoints=XLALCalloc(Nlive,sizeof(LALVariables *)); */
	runState->livePoints=XLALCalloc(Nlive,sizeof(LALVariables *));
	if(runState->livePoints==NULL)
	{
		fprintf(stderr,"Unable to allocate memory for %i live points\n",Nlive);
		exit(1);
	}
	
	logLs=XLALCreateREAL8Vector(Nlive);
	addVariable(runState->algorithmParams,"logLikelihoods",&logLs,REAL8Vector_t,PARAM_FIXED);
	fprintf(stdout,"Sprinkling %i live points, may take some time\n",Nlive);
	for(i=0;i<Nlive;i++)
	{
		runState->livePoints[i]=XLALCalloc(1,sizeof(LALVariables));
		
		/* Copy the param structure */
		copyVariables(runState->currentParams,runState->livePoints[i]);
		
		/* Sprinkle the varying points among prior */
		do{
			for(current=runState->livePoints[i]->head ;current!=NULL;
				current=current->next){
				if(current->vary==PARAM_CIRCULAR || current->vary==PARAM_LINEAR)
				{
					switch (current->type){
					case REAL4_t:
					{
						REAL4 tmp;
						REAL4 min,max;
						getMinMaxPrior(runState->priorArgs,current->name, 
									   (void *)&min,(void *)&max);
						tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						setVariable(runState->livePoints[i],current->name,&tmp);
						break;
					}
						
					case REAL8_t:
					{
						REAL8 tmp;
						REAL8 min,max;
						getMinMaxPrior(runState->priorArgs,current->name, 
									   (void *)&min,(void *)&max);
						tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						setVariable(runState->livePoints[i],current->name,&tmp);
						break;
					}
					case INT4_t:
					{
						INT4 tmp;
						INT4 min,max;
						getMinMaxPrior(runState->priorArgs,current->name,
									   (void *)&min,(void *)&max);
						tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						setVariable(runState->livePoints[i],current->name,&tmp);
						break;
					}
					case INT8_t:
					{
						INT8 tmp;
						INT8 min,max;
						getMinMaxPrior(runState->priorArgs,current->name,
									   (void *)&min,(void *)&max);
						tmp=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
						setVariable(runState->livePoints[i],current->name,&tmp);
						break;
					}
					default:
						fprintf(stderr,"Trying to randomise a non-numeric parameter!");
					}
				}
			}
		}while(runState->prior(runState,runState->livePoints[i])==-DBL_MAX);
		/* Populate log likelihood */
		logLs->data[i]=runState->likelihood(runState->livePoints[i],runState->data,runState->template);
	}
	
}

/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void initVariables(LALInferenceRunState *state)
{
	LALStatus status;
	SimInspiralTable *injTable=NULL;
	LALVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALVariables));
	LALVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	INT4 AmpOrder=0;
	LALPNOrder PhaseOrder=LAL_PNORDER_TWO;
	Approximant approx=TaylorF2;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
	REAL8 mcMin=1.0;
	REAL8 mcMax=20.5;
	REAL8 logmcMax,logmcMin,mMin=1.0,mMax=30.0;
	REAL8 etaMin=0.01;
	REAL8 etaMax=0.25;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;
	
	memset(currentParams,0,sizeof(LALVariables));
	
	char help[]="\
	[--injXML injections.xml]\tInjection XML file to use\n\
	[--Mmin mchirp]\tMinimum chirp mass\n\
	[--Mmax mchirp]\tMaximum chirp mass\n\
	[--dt time]\tWidth of time prior, centred around trigger (0.1s)\n\
	[--trigtime time]\tTrigger time to use\n\
	[--Dmin dist]\tMinimum distance in Mpc (1)\n\
	[--Dmax dist]\tMaximum distance in Mpc (100)\n\
	[--approx ApproximantorderPN]\tSpecify a waveform to use, (default TaylorF2twoPN)\n\
	[--mincomp min]\tMinimum component mass (1.0)\n\
	[--maxcomp max]\tMaximum component mass (30.0)\n";
	
	/* Print command line arguments if help requested */
	ppt=getProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}
	
	/* Read injection XML file for parameters if specified */
	ppt=getProcParamVal(commandLine,"--injXML");
	if(ppt){
		SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
		if(!injTable){
			fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
			exit(1);
		}
		endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
		AmpOrder=injTable->amp_order;
		LALGetOrderFromString(&status,injTable->waveform,&PhaseOrder);
		LALGetApproximantFromString(&status,injTable->waveform,&approx);
	}	

	/* Over-ride approximant if user specifies */
	ppt=getProcParamVal(commandLine,"--approx");
	if(ppt){
		LALGetOrderFromString(&status,ppt->value,&PhaseOrder);
		LALGetApproximantFromString(&status,ppt->value,&approx);
		if(strstr(ppt->value,"TaylorF2")) approx=TaylorF2;
		fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",approx,PhaseOrder);
	}
	
	/* Over-ride end time if specified */
	ppt=getProcParamVal(commandLine,"--trigtime");
	if(ppt){
		endtime=atof(ppt->value);
	}
	
	/* Over-ride time prior if specified */
	ppt=getProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}
	
	/* Over-ride Distance min if specified */
	ppt=getProcParamVal(commandLine,"--Dmin");
	if(ppt){
		logDmin=log(atof(ppt->value));
	}
	
	/* Over-ride Distance max if specified */
	ppt=getProcParamVal(commandLine,"--Dmax");
	if(ppt){
		logDmax=log(atof(ppt->value));
	}
	
	/* Over-ride Mass prior if specified */
	ppt=getProcParamVal(commandLine,"--Mmin");
	if(ppt){
		mcMin=atof(ppt->value);
	}
	ppt=getProcParamVal(commandLine,"--Mmax");
	if(ppt)	mcMax=atof(ppt->value);
	
	/* Over-ride component masses */
	ppt=getProcParamVal(commandLine,"--compmin");
	if(ppt)	mMin=atof(ppt->value);
	addVariable(priorArgs,"component_min",&mMin,REAL8_t,PARAM_FIXED);
	ppt=getProcParamVal(commandLine,"--compmax");
	if(ppt)	mMax=atof(ppt->value);
	addVariable(priorArgs,"component_max",&mMax,REAL8_t,PARAM_FIXED);
	
	
	printf("Read end time %f\n",endtime);
	
	addVariable(currentParams, "LAL_APPROXIMANT", &approx,        INT4_t, PARAM_FIXED);
    addVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        INT4_t, PARAM_FIXED);
	
	/* Set up the variable parameters */
	tmpVal=log(mcMin+(mcMax-mcMin)/2.0);
	/*addVariable(currentParams, "chirpmass",    &tmpVal,    REAL8_t,	PARAM_LINEAR);
    addMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		REAL8_t); */
	addVariable(currentParams,"logmc",&tmpVal, REAL8_t, PARAM_LINEAR);
	logmcMin=log(mcMin); logmcMax=log(mcMax);
	addMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		REAL8_t);

	tmpVal=0.24;
	addVariable(currentParams, "massratio",       &tmpVal,             REAL8_t, PARAM_LINEAR);
    addMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	REAL8_t);
	
    addVariable(currentParams, "time",            &endtime   ,           REAL8_t, PARAM_LINEAR); 
	tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
	addMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   REAL8_t);	

	tmpVal=1.0;
    addVariable(currentParams, "phase",           &tmpVal,             REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	addMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpVal=logDmin+(logDmax-logDmin)/2.0;
	addVariable(currentParams,"logdistance", &tmpVal, REAL8_t, PARAM_LINEAR);
	addMinMaxPrior(priorArgs, "logdistance",     &logDmin, &logDmax,   REAL8_t);
	
	tmpVal=1.0;
	addVariable(currentParams, "rightascension",  &tmpVal,      REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	addMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   REAL8_t);

	addVariable(currentParams, "declination",     &tmpVal,     REAL8_t, PARAM_CIRCULAR);
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	addMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   REAL8_t);
    
	addVariable(currentParams, "polarisation",    &tmpVal,     REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	addMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   REAL8_t);
	
 	addVariable(currentParams, "inclination",     &tmpVal,            REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	addMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   REAL8_t);
	
	return;
}

/*************** MAIN **********************/


int main(int argc, char *argv[]){

	LALInferenceRunState *state;
	ProcessParamsTable *procParams=NULL;
	
	/* Read command line and parse */
	procParams=parseCommandLine(argc,argv);
	
	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	state = initialize(procParams);
	
	/* Set up structures for nested sampling */
	initializeNS(state);
	
	/* Set up currentParams with variables to be used */
	initVariables(state);
	
	/* Call setupLivePointsArray() to populate live points structures */
	setupLivePointsArray(state);
	
	/* Call nested sampling algorithm */
	state->algorithm(state);
	
	
	/* end */
	return(0);
}
