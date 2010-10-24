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
#include "LALInference.h"
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
//#include "LALInferenceMCMCMPISampler.h"
#include "LALInferencePrior.h"

#include <mpi.h>
#include "mpi.h"


int MPIrank, MPIsize;


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMCMC(LALInferenceRunState *runState);
void initVariables(LALInferenceRunState *state);




LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	LALInferenceRunState *irs=NULL;
	LALIFOData *ifoPtr, *ifoListStart;
	//ProcessParamsTable *ppt=NULL;

	//int MPIrank;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	
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
	
	
	return(irs);
}

/***** Initialise Nested Sampling structures ****/
/* Fill in samples from the prior distribution */
/* runState->algorithmParams must contain a variable "logLikelihoods" */
/* which contains a REAL8 array of likelihood values for the live */
/* points. */
/************************************************/
void initializeMCMC(LALInferenceRunState *runState)
{
	char help[]="\
	--Niter N\tNumber of iteration\n\
	[--tempMax T]\tHighest temperature for parallel tempering(40.0)\n\
	[--randomseed seed]\tRandom seed of sampling distribution\n";
	
	INT4 verbose=0,tmpi=0;
	unsigned long int randomseed=0;
	REAL8 tempMax = 40.0;
	//REAL8 tmp=0;
	ProcessParamsTable *commandLine=runState->commandLine;
	ProcessParamsTable *ppt=NULL;
	FILE *devrandom;
	struct timeval tv;
	
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
	
	/* Set up the appropriate functions for the MCMC algorithm */
	runState->algorithm=PTMCMCAlgorithm;
	runState->evolve=PTMCMCOneStep;
	runState->proposal=PTMCMCLALProposal;
	//runState->proposal=PTMCMCLALAdaptationProposal;
	//runState->proposal=PTMCMCGaussianProposal;
	
	/* This is the LAL template generator for inspiral signals */
	runState->template=templateLAL;
	runState->likelihood=FreqDomainLogLikelihood;
	//runState->likelihood=GaussianLikelihood;
	runState->prior=PTUniformLALPrior;//LALInferenceInspiralPriorNonSpinning;
	//runState->prior=PTUniformGaussianPrior;

	
	
	ppt=getProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		addVariable(runState->algorithmParams,"verbose", &verbose , INT4_t,
					PARAM_FIXED);
	}
	
	printf("set iteration number.\n");
	/* Number of live points */
	ppt=getProcParamVal(commandLine,"--Niter");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		fprintf(stderr,"Error, must specify iteration number\n");
		MPI_Finalize();
		exit(1);
	}
	addVariable(runState->algorithmParams,"Niter",&tmpi, INT4_t,PARAM_FIXED);
	
	printf("set highest temperature.\n");
	/* Tolerance of the Nested sampling integrator */
	ppt=getProcParamVal(commandLine,"--tempMax");
	if(ppt){
		tempMax=strtod(ppt->value,(char **)NULL);
	}	
	addVariable(runState->algorithmParams,"tempMax",&tempMax, REAL8_t,PARAM_FIXED);
	
	printf("set random seed.\n");
	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = getProcParamVal(commandLine, "--randomseed");
	if (ppt != NULL)
		randomseed = atoi(ppt->value);
	else { /* otherwise generate "random" random seed: */
		if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
			if (MPIrank == 0) {
				gettimeofday(&tv, 0);
				randomseed = tv.tv_sec + tv.tv_usec;
			}
			MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		} 
		else {
			if (MPIrank == 0) {
				fread(&randomseed, sizeof(randomseed), 1, devrandom);
				fclose(devrandom);
			}
			MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stdout, " initialize(): random seed: %lu\n", randomseed);
	addVariable(runState->algorithmParams,"random_seed",&randomseed, INT4_t,PARAM_FIXED);
	gsl_rng_set(runState->GSLrandom, randomseed);
	
	return;
	
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
	REAL8 Dmin=1.0;
	REAL8 Dmax=100.0;
	REAL8 mcMin=1.0;
	REAL8 mcMax=20.5;
	//REAL8 logmcMax,logmcMin,mMin=1.0,mMax=30.0;
	REAL8 etaMin=0.01;
	REAL8 etaMax=0.25;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;
	
	memset(currentParams,0,sizeof(LALVariables));
	
	char help[]="\
	[--injXML injections.xml]\tInjection XML file to use\
	[--Mmin mchirp]\tMinimum chirp mass\
	[--Mmax mchirp]\tMaximum chirp mass\
	[--dt time]\tWidth of time prior, centred around trigger (0.1s)\
	[--trigtime time]\tTrigger time to use\
	[--Dmin dist]\tMinimum distance in Mpc (1)\
	[--Dmax dist]\tMaximum distance in Mpc (100)\
	[--approx ApproximantorderPN]\tSpecify a waveform to use, (default TaylorF2twoPN)\
	[--mincomp min]\tMinimum component mass (1.0)\
	[--maxcomp max]\tMaximum component mass (30.0)";
	
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
			MPI_Finalize();
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
		Dmin=atof(ppt->value);
	}
	
	/* Over-ride Distance max if specified */
	ppt=getProcParamVal(commandLine,"--Dmax");
	if(ppt){
		logDmax=log(atof(ppt->value));
		Dmax=atof(ppt->value);
	}
	
	/* Over-ride Mass prior if specified */
	ppt=getProcParamVal(commandLine,"--Mmin");
	if(ppt){
		mcMin=atof(ppt->value);
	}
	ppt=getProcParamVal(commandLine,"--Mmax");
	if(ppt)	mcMax=atof(ppt->value);
	
	/* Over-ride component masses */
	//ppt=getProcParamVal(commandLine,"--compmin");
	//if(ppt)	mMin=atof(ppt->value);
	//addVariable(priorArgs,"component_min",&mMin,REAL8_t,PARAM_FIXED);
	//ppt=getProcParamVal(commandLine,"--compmax");
	//if(ppt)	mMax=atof(ppt->value);
	//addVariable(priorArgs,"component_max",&mMax,REAL8_t,PARAM_FIXED);
	
	
	printf("Read end time %f\n",endtime);
	INT4 numberI4 = TaylorF2;
	//INT4 numberI4 = TaylorT4;
	//addVariable(currentParams, "LAL_APPROXIMANT", &approx,        INT4_t, PARAM_FIXED);
	addVariable(currentParams, "LAL_APPROXIMANT", &numberI4,        INT4_t, PARAM_FIXED);
	numberI4 = LAL_PNORDER_TWO;
    //addVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        INT4_t, PARAM_FIXED);	
	addVariable(currentParams, "LAL_PNORDER",     &numberI4,        INT4_t, PARAM_FIXED);
	
	/* Set up the variable parameters */
	tmpVal=7.0;//log(mcMin+(mcMax-mcMin)/2.0);
	tmpVal=7.86508;
	addVariable(currentParams, "chirpmass",    &tmpVal,    REAL8_t,	PARAM_LINEAR);
	addMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		REAL8_t);
	//addVariable(currentParams,"logmc",&tmpVal, REAL8_t, PARAM_LINEAR);
	//logmcMin=log(mcMin); logmcMax=log(mcMax);
	//addMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		REAL8_t);
	
	tmpVal=0.2451;
	tmpVal=0.18957;
	addVariable(currentParams, "massratio",       &tmpVal,             REAL8_t, PARAM_LINEAR);
    addMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	REAL8_t);
	
    addVariable(currentParams, "time",            &endtime   ,           REAL8_t, PARAM_LINEAR); 
	tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
	addMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   REAL8_t);	
	
	tmpVal=3.974985;
	tmpVal=3.89954;
    addVariable(currentParams, "phase",           &tmpVal,             REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	addMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpVal=30.00;//Dmin+(Dmax-Dmin)/2.0;
	tmpVal=46.92314;
	addVariable(currentParams,"distance", &tmpVal, REAL8_t, PARAM_LINEAR);
	addMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   REAL8_t);
	
	tmpVal=3.001677;//1.0;
	tmpVal=3.34650;
	addVariable(currentParams, "rightascension",  &tmpVal,      REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	addMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpVal=-1.2;
	tmpVal=-0.90547;
	addVariable(currentParams, "declination",     &tmpVal,     REAL8_t, PARAM_CIRCULAR);
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	addMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   REAL8_t);
    
	tmpVal=5.558779;
	tmpVal=0.64546;
	addVariable(currentParams, "polarisation",    &tmpVal,     REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	addMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpVal=0.2;
	tmpVal=2.86094;
 	addVariable(currentParams, "inclination",     &tmpVal,            REAL8_t, PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	addMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   REAL8_t);
	
	
	//REAL8 x0 = 0.9;
	//addVariable(currentParams, "x0", &x0,  REAL8_t, PARAM_LINEAR);

	
	return;
}




int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

  if (MPIrank == 0) fprintf(stdout," ========== LALInference_MCMCMPI ==========\n");

	LALInferenceRunState *runState;
	ProcessParamsTable *procParams=NULL;
	
	/* Read command line and parse */
	procParams=parseCommandLine(argc,argv);

	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	runState = initialize(procParams);
  
	/* Set up structures for nested sampling */
	initializeMCMC(runState);
	
	/* Set up currentParams with variables to be used */
	initVariables(runState);
	
	printf(" ==== This is thread %d of %d ====\n ", MPIrank, MPIsize);
	MPI_Barrier(MPI_COMM_WORLD);
	/* Call MCMC algorithm */
	runState->algorithm(runState);
	
  if (MPIrank == 0) printf(" ========== main(): finished. ==========\n");
  MPI_Finalize();
  return 0;
}




//void PTMCMCTest(void)
//{
//	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
//	MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
//
//	fprintf(stdout, "PTMCMC test\n");
//
//	runstate->algorithm=PTMCMCAlgorithm;
//	runstate->evolve=PTMCMCOneStep;
//	runstate->prior=PTUniformLALPrior;
//	//runstate->prior=PTUniformGaussianPrior;
//	runstate->proposal=PTMCMCLALProposal;
//	//runstate->proposal=PTMCMCLALAdaptationProposal;
//	//runstate->proposal=PTMCMCGaussianProposal;
//	runstate->proposalArgs = malloc(sizeof(LALVariables));
//	runstate->proposalArgs->head=NULL;
//	runstate->proposalArgs->dimension=0;
//	runstate->likelihood=FreqDomainLogLikelihood;
//	//runstate->likelihood=GaussianLikelihood;
//	runstate->template=templateLAL;
//	
//	
//	SimInspiralTable *injTable=NULL;
//	printf("Ninj: %d\n", SimInspiralTableFromLIGOLw(&injTable,getProcParamVal(ppt,"--injXML")->value,0,0));
//	
//	REAL8 mc = injTable->mchirp;
//	REAL8 eta = injTable->eta;
//    REAL8 iota = injTable->inclination;
//    REAL8 phi = injTable->coa_phase;
//	LIGOTimeGPS trigger_time=injTable->geocent_end_time;
//	REAL8 tc = XLALGPSGetREAL8(&trigger_time);
//	REAL8 ra_current = injTable->longitude;
//	REAL8 dec_current = injTable->latitude;
//	REAL8 psi_current = injTable->polarization;
//	REAL8 distMpc_current = injTable->distance;
//	
//    numberI4 = TaylorF2;
//    addVariable(&currentParams, "LAL_APPROXIMANT", &numberI4,        INT4_t, PARAM_FIXED);
//    numberI4 = LAL_PNORDER_TWO;
//    addVariable(&currentParams, "LAL_PNORDER",     &numberI4,        INT4_t, PARAM_FIXED);
//	
//	addVariable(&currentParams, "chirpmass",       &mc,              REAL8_t, PARAM_LINEAR);
//    addVariable(&currentParams, "massratio",       &eta,             REAL8_t, PARAM_LINEAR);
//    addVariable(&currentParams, "inclination",     &iota,            REAL8_t, PARAM_CIRCULAR);
//    addVariable(&currentParams, "phase",           &phi,             REAL8_t, PARAM_CIRCULAR);
//    addVariable(&currentParams, "time",            &tc   ,           REAL8_t, PARAM_LINEAR); 
//    addVariable(&currentParams, "rightascension",  &ra_current,      REAL8_t, PARAM_CIRCULAR);
//    addVariable(&currentParams, "declination",     &dec_current,     REAL8_t, PARAM_CIRCULAR);
//    addVariable(&currentParams, "polarisation",    &psi_current,     REAL8_t, PARAM_CIRCULAR);
//    addVariable(&currentParams, "distance",        &distMpc_current, REAL8_t, PARAM_LINEAR);
//	
//	
////	REAL8 x0 = 0.9;
////	addVariable(&currentParams, "x0", &x0,  REAL8_t, PARAM_LINEAR);
//	
//	
//	
//	
//	runstate->currentParams=&currentParams;
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	PTMCMCAlgorithm(runstate);
//	if (MPIrank == 0) fprintf(stdout, "End of PTMCMC test\n");
//}
//
