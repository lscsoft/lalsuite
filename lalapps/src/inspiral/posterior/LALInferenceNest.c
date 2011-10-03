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
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lalapps.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeNS(LALInferenceRunState *runState);
void initVariables(LALInferenceRunState *state);
void initStudentt(LALInferenceRunState *state);
void initializeTemplate(LALInferenceRunState *runState);
static void mc2masses(double mc, double eta, double *m1, double *m2);

static void mc2masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m1 >= m2).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m2 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m1 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
  return;
}


LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	char help[]="\
Initialisation arguments:\n\
(--randomseed seed           Random seed for Nested Sampling)\n\n";
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr, *ifoListStart;
	ProcessParamsTable *ppt=NULL;
	unsigned long int randomseed;
	struct timeval tv;
	FILE *devrandom;
	
	irs = calloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " readData(): started.\n");
	irs->commandLine=commandLine;
	irs->data = LALInferenceReadData(commandLine);
	/* (this will already initialise each LALIFOData's following elements:  */
        ppt=LALInferenceGetProcParamVal(commandLine,"--help");
        if(ppt)
        {
                fprintf(stdout,"%s",help);
                return(irs);
        }

	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " LALInferenceReadData(): finished.\n");
	if (irs->data != NULL) {
		fprintf(stdout, " initialize(): successfully read data.\n");
		
		fprintf(stdout, " LALInferenceInjectInspiralSignal(): started.\n");
		LALInferenceInjectInspiralSignal(irs->data,commandLine);
		fprintf(stdout, " LALInferenceInjectInspiralSignal(): finished.\n");
		
		ifoPtr = irs->data;
		ifoListStart = irs->data;
		while (ifoPtr != NULL) {
			/*If two IFOs have the same sampling rate, they should have the same timeModelh*,
			 freqModelh*, and modelParams variables to avoid excess computation 
			 in model waveform generation in the future*/
			LALInferenceIFOData * ifoPtrCompare=ifoListStart;
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
				ifoPtr->modelParams = calloc(1, sizeof(LALInferenceVariables));
			}
			ifoPtr = ifoPtr->next;
		}
		irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
		printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else
		fprintf(stdout, " initialize(): no data read.\n");
	
	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	irs->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
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

void initializeTemplate(LALInferenceRunState *runState)
{
	char help[]="\
(--template [LAL,PhenSpin,LALGenerateInspiral]\tSpecify template (default LAL)\n";
	ProcessParamsTable *ppt=NULL;
	ProcessParamsTable *commandLine=runState->commandLine;
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}
	/* This is the LAL template generator for inspiral signals */
	runState->template=&LALInferenceTemplateLAL;
	ppt=LALInferenceGetProcParamVal(commandLine,"--template");
	if(ppt) {
		if(!strcmp("LALSTPN",ppt->value)){
			fprintf(stderr,"ERROR: --template LALSTPN is deprecated. Try LALGenerateInspiral instead\n");
			exit(1);
		}
		else if(!strcmp("PhenSpin",ppt->value))
			runState->template=&LALInferenceTemplatePSTRD;
		else if(!strcmp("LALGenerateInspiral",ppt->value))
			runState->template=&LALInferenceTemplateLALGenerateInspiral;
		else if(!strcmp("LAL",ppt->value))
			runState->template=&LALInferenceTemplateLAL;
		else {
			XLALPrintError("Error: unknown template %s\n",ppt->value);
			XLALPrintError(help);
			XLAL_ERROR_VOID(__func__,XLAL_EINVAL);
		}

	}
	return;
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
Nested sampling arguments:\n\
 --Nlive N\tNumber of live points to use\n\
 --Nmcmc M\tNumber of MCMC point to use when evolving live points\n\
(--Nruns R)\tNumber of parallel samples from logt to use(1)\n\
(--tolerance dZ)\tTolerance of nested sampling algorithm (0.1)\n\
(--randomseed seed)\tRandom seed of sampling distribution\n\
(--verbose)\tProduce progress information\n\
(--mcmcprop)\tUse PTMCMC proposal engine\n\
\t(--iotaDistance FRAC)\tPTMCMC: Use iota-distance jump FRAC of the time\n\
\t(--covarianceMatrix)\tPTMCMC: Propose jumps from covariance matrix of current live points\n\
\t(--differential-evolution)\tPTMCMC:Use differential evolution jumps\n";

	ProcessParamsTable *ppt=NULL;
	ProcessParamsTable *commandLine=runState->commandLine;
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}

	INT4 verbose=0,tmpi=0,randomseed=0;
	REAL8 tmp=0;
	
	/* Initialise parameters structure */
	runState->algorithmParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	runState->priorArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
	runState->proposalArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
	
	/* Set up the appropriate functions for the nested sampling algorithm */
	runState->algorithm=&LALInferenceNestedSamplingAlgorithm;
	runState->evolve=&LALInferenceNestedSamplingOneStep;
	if(LALInferenceGetProcParamVal(commandLine,"--mcmcprop")){
	  /* Use the PTMCMC proposal to sample prior */
	  runState->proposal=&NSWrapMCMCLALProposal;
	  REAL8 temp=1.0;
	  UINT4 dummy=0;
	  LALInferenceAddVariable(runState->proposalArgs, "adaptableStep", &dummy, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);
	  LALInferenceAddVariable(runState->proposalArgs, "proposedVariableNumber", &dummy, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);
	  LALInferenceAddVariable(runState->proposalArgs, "proposedArrayNumber", &dummy, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);
	  LALInferenceAddVariable(runState->proposalArgs,"temperature",&temp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	}
	else
	  runState->proposal=&LALInferenceProposalNS;

	runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
	runState->prior = &LALInferenceInspiralPrior;
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose , LALINFERENCE_INT4_t,
					LALINFERENCE_PARAM_FIXED);		
	}
	if(verbose) set_debug_level("ERROR|INFO");
	else set_debug_level("NDEBUG");
		
	printf("set number of live points.\n");
	/* Number of live points */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nlive");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		fprintf(stderr,"Error, must specify number of live points\n");
		exit(1);
	}
	LALInferenceAddVariable(runState->algorithmParams,"Nlive",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	
	printf("set number of MCMC points.\n");
	/* Number of points in MCMC chain */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nmcmc");
	if(ppt)
	  tmpi=atoi(ppt->value);
	else {
	  fprintf(stderr,"Error, must specify number of MCMC points\n");
	  exit(1);
	}
	LALInferenceAddVariable(runState->algorithmParams,"Nmcmc",&tmpi,
				LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	
	printf("set number of parallel runs.\n");
	/* Optionally specify number of parallel runs */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nruns");
	if(ppt) {
		tmpi=atoi(ppt->value);
		LALInferenceAddVariable(runState->algorithmParams,"Nruns",&tmpi,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	}
	
	printf("set tolerance.\n");
	/* Tolerance of the Nested sampling integrator */
	ppt=LALInferenceGetProcParamVal(commandLine,"--tolerance");
	if(ppt){
		tmp=strtod(ppt->value,(char **)NULL);
		LALInferenceAddVariable(runState->algorithmParams,"tolerance",&tmp, LALINFERENCE_REAL8_t,
					LALINFERENCE_PARAM_FIXED);
	}
	
	printf("set random seed.\n");
	/* Set up the random number generator */
	gsl_rng_env_setup();
	runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	
	/* (try to) get random seed from command line: */
	ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	fprintf(stdout, " initialize(): random seed: %u\n", randomseed);
	LALInferenceAddVariable(runState->algorithmParams,"random_seed",&randomseed, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	gsl_rng_set(runState->GSLrandom, randomseed);
	
	return;
	
}

/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void initVariables(LALInferenceRunState *state)
{
	LALStatus status;
	SimInspiralTable *injTable=NULL;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	INT4 AmpOrder=0;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	Approximant approx=TaylorF2;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
	REAL8 mcMin=1.0;
	REAL8 mcMax=20.5;
	REAL8 logmcMax,logmcMin,mMin=1.0,mMax=30.0;
	REAL8 a_spin2_max=1.0, a_spin1_max=1.0;
	REAL8 a_spin2_min=0.0, a_spin1_min=0.0;
	REAL8 phi_spin1_min=-LAL_PI;
	REAL8 phi_spin1_max=LAL_PI;
	REAL8 theta_spin1_min=-LAL_PI/2.0;
	REAL8 theta_spin1_max=LAL_PI/2.0;	
	REAL8 etaMin=0.01;
	REAL8 etaMax=0.25;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;
	REAL8 m1_min=0.;	
	REAL8 m1_max=0.;
	REAL8 m2_min=0.;
	REAL8 m2_max=0.;
	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
	INT4 enable_spin=0;
	INT4 aligned_spin=0;
	char help[]="\
Parameter arguments:\n\
(--injXML injections.xml)\tInjection XML file to use\n\
(--Mmin mchirp)\tMinimum chirp mass\n\
(--Mmax mchirp)\tMaximum chirp mass\n\
(--etamin eta)\tMinimum eta\n\
(--etamax eta)\tMaximum eta\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--Dmin dist)\tMinimum distance in Mpc (1)\n\
(--Dmax dist)\tMaximum distance in Mpc (100)\n\
(--approx ApproximantorderPN)\tSpecify a waveform to use, (default TaylorF2threePointFivePN)\n\
(--compmin min)\tMinimum component mass (1.0)\n\
(--compmax max)\tMaximum component mass (30.0)\n\
(--enable-spin)\tEnable spin parameters\n\
(--aligned-spin)\tUse only aligned spin parameters (uses spins between -1 and 1)\n\
(--approx ApproximantphaseOrderPN)\tSet approximant (PhenSpin implicitly enables spin)\n\
(--s1max SPIN)\tMax magnitude of spin (on both bodies!)\n\
(--s1min SPIN)\tMin magnitude of spin (on both bodies!)\n";

	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return;
	}

	
	/* Read injection XML file for parameters if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--injXML");
	if(ppt){
		SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
		if(!injTable){
			fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
			exit(1);
		}
	}
	//Select event
	ppt=LALInferenceGetProcParamVal(commandLine,"--event");
	if(ppt){
		event = atoi(ppt->value);
		while(i<event) {i++; injTable = injTable->next;}
		endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));}
		AmpOrder=injTable->amp_order;
		LALGetOrderFromString(&status,injTable->waveform,&PhaseOrder);
		LALGetApproximantFromString(&status,injTable->waveform,&approx);
	

	/* Over-ride approximant if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
	if(ppt){
		if(strstr(ppt->value,"TaylorF2")) approx=TaylorF2;
		else
		    LALGetApproximantFromString(&status,ppt->value,&approx);
        LALGetOrderFromString(&status,ppt->value,&PhaseOrder);
	}
	fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",approx,PhaseOrder);

	/* Over-ride end time if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
	if(ppt){
		endtime=atof(ppt->value);
	}
	
	/* Over-ride time prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}
	
	/* Over-ride Distance min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
	if(ppt){
		logDmin=log(atof(ppt->value));
	}
	
	/* Over-ride Distance max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
	if(ppt){
		logDmax=log(atof(ppt->value));
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--etamin");
        if(ppt)
                etaMin=atof(ppt->value);

        ppt=LALInferenceGetProcParamVal(commandLine,"--etamax");
	if(ppt)
                etaMax=atof(ppt->value);
	/* Over-ride Mass prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Mmin");
	if(ppt){
		mcMin=atof(ppt->value);
		mc2masses( mcMin,  etaMin,  &m1_min,  &m2_min);
		mMin=m2_min;
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--Mmax");
	if(ppt){	
		mcMax=atof(ppt->value);
		mc2masses(mcMax, etaMax, &m1_max, &m2_max);
		mMax=m1_max;
	}
	/* Over-ride Spin prior if specified*/

	ppt=LALInferenceGetProcParamVal(commandLine,"--s1max");
	if(ppt){
		a_spin2_max=atof(ppt->value);
		a_spin1_max=atof(ppt->value);
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--s1min");
	if(ppt){
		a_spin2_min=atof(ppt->value);
		a_spin1_min=atof(ppt->value);
	}
	/* Over-ride component masses */
	ppt=LALInferenceGetProcParamVal(commandLine,"--compmin");
	if(ppt)	mMin=atof(ppt->value);
	//fprintf(stderr,"Mmin %f, Mmax %f\n",mMin,mMax);
	LALInferenceAddVariable(priorArgs,"component_min",&mMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	ppt=LALInferenceGetProcParamVal(commandLine,"--compmax");
	if(ppt)	mMax=atof(ppt->value);
	LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	
	printf("Read end time %f\n",endtime);
	
	LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    	LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	
	/* Set up the variable parameters */
	tmpVal=log(mcMin+(mcMax-mcMin)/2.0);
	/*LALInferenceAddVariable(currentParams, "chirpmass",    &tmpVal,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t); */
	LALInferenceAddVariable(currentParams,"logmc",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	logmcMin=log(mcMin); logmcMax=log(mcMax);
	LALInferenceAddMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		LALINFERENCE_REAL8_t);

	tmpVal=0.24;
	LALInferenceAddVariable(currentParams, "massratio",       &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    	LALInferenceAddMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	LALINFERENCE_REAL8_t);
	
    	LALInferenceAddVariable(currentParams, "time",            &endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
	LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	

	tmpVal=1.0;
    	LALInferenceAddVariable(currentParams, "phase",           &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	tmpVal=logDmin+(logDmax-logDmin)/2.0;
	LALInferenceAddVariable(currentParams,"logdistance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddMinMaxPrior(priorArgs, "logdistance",     &logDmin, &logDmax,   LALINFERENCE_REAL8_t);
	
	tmpVal=1.0;
	LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

	LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
	LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
 	LALInferenceAddVariable(currentParams, "inclination",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	/* Additional parameters for spinning waveforms */
	ppt=LALInferenceGetProcParamVal(commandLine,"--template");
	if(ppt) if(!strcmp("PhenSpin",ppt->value)){ enable_spin=1;}

	if(LALInferenceGetProcParamVal(commandLine,"--enable-spin")) enable_spin=1;
	
	/* If aligned spins use magnitude in (-1,1) */
	ppt=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
	if(ppt) {enable_spin=1; aligned_spin=1; a_spin1_min=-1; a_spin2_min=-1;}
		
	if(enable_spin){
		tmpVal=a_spin1_min+(a_spin1_max-a_spin1_min)/2.0;
		LALInferenceAddVariable(currentParams, "a_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &a_spin1_min, &a_spin1_max,   LALINFERENCE_REAL8_t); 
	        
		tmpVal=a_spin2_min+(a_spin2_max-a_spin2_min)/2.0;
		LALInferenceAddVariable(currentParams, "a_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &a_spin2_min, &a_spin2_max,   LALINFERENCE_REAL8_t); 
	
		
		if(aligned_spin){ /* Set the spin angles to be parallel to orbital */
			tmpVal=LAL_PI/2;
			LALInferenceAddVariable(currentParams,"theta_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			LALInferenceAddVariable(currentParams,"theta_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			tmpVal=0;
			LALInferenceAddVariable(currentParams,"phi_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			LALInferenceAddVariable(currentParams,"phi_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		}
		else{ /* Use full spinning parameters */
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)/2.0;

			LALInferenceAddVariable(currentParams,"theta_spin1",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)/2.0;
			LALInferenceAddVariable(currentParams,"theta_spin2",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)/2.0;
	
			LALInferenceAddVariable(currentParams,"phi_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &phi_spin1_min, &phi_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)/2.0;
			LALInferenceAddVariable(currentParams,"phi_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &phi_spin1_min, &phi_spin1_max,   LALINFERENCE_REAL8_t);
		}
	}
	
	return;
}

/** Initialise student-t extra variables, set likelihood */
void initStudentt(LALInferenceRunState *state)
{
        char help[]="\
Student T Likelihood Arguments:\n\
(--studentt)\tUse student-t likelihood function\n";

        ProcessParamsTable *ppt=NULL;
	LALInferenceIFOData *ifo=state->data;

	/* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
        {
                fprintf(stdout,"%s",help);
		while(ifo) {
			fprintf(stdout,"(--dof-%s DoF)\tDegrees of freedom for %s\n",ifo->name,ifo->name);
			ifo=ifo->next;
		}
		return;
        }
	/* Don't do anything unless asked */
	if(!LALInferenceGetProcParamVal(state->commandLine,"--studentt")) return;

	/* initialise degrees of freedom parameters for each IFO */
	while(ifo){
		CHAR df_argument_name[128];
		CHAR df_variable_name[64];
		REAL8 dof=10.0; /* Degrees of freedom parameter */
		
		sprintf(df_argument_name,"--dof-%s",ifo->name);
		if((ppt=LALInferenceGetProcParamVal(state->commandLine,df_argument_name)))
			dof=atof(ppt->value);
    		sprintf(df_variable_name,"df_%s",ifo->name);
    		LALInferenceAddVariable(state->currentParams,df_variable_name,&dof,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
		fprintf(stdout,"Setting %lf degrees of freedom for %s\n",dof,ifo->name);
		ifo=ifo->next;
	}

	/* Set likelihood to student-t */
	state->likelihood = &LALInferenceFreqDomainStudentTLogLikelihood;
	
	/* Set the noise model evidence to the student t model value */
	LALInferenceTemplateNullFreqdomain(state->data);
	REAL8 noiseZ=LALInferenceFreqDomainStudentTLogLikelihood(state->currentParams,state->data,&LALInferenceTemplateNullFreqdomain);
	LALInferenceAddVariable(state->algorithmParams,"logZnoise",&noiseZ,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	fprintf(stdout,"Student-t Noise evidence %lf\n",noiseZ);

	return;
}

/*************** MAIN **********************/


int main(int argc, char *argv[]){
        char help[]="\
LALInferenceNest:\n\
Bayesian analysis tool using Nested Sampling algorithm\n\
for CBC analysis. Uses LALInference library for back-end.\n\n\
Arguments for each section follow:\n\n";

	LALInferenceRunState *state;
	ProcessParamsTable *procParams=NULL;

	/* Read command line and parse */
	procParams=LALInferenceParseCommandLine(argc,argv);
	
	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	state = initialize(procParams);
	
	/* Set template function */
	initializeTemplate(state);
	
	/* Set up structures for nested sampling */
	initializeNS(state);
	
	/* Set up currentParams with variables to be used */
	initVariables(state);
	
	/* Check for student-t and apply */
	initStudentt(state);

       /* Print command line arguments if help requested */
        if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
        {
                fprintf(stdout,"%s",help);
		exit(0);
        }

	/* Call setupLivePointsArray() to populate live points structures */
	LALInferenceSetupLivePointsArray(state);

	/* Call nested sampling algorithm */
	state->algorithm(state);

	/* end */
	return(0);
}
