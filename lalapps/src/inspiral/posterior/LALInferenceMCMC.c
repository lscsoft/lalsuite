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
#include "LALInferenceMCMCSampler.h"
#include "LALInferencePrior.h"


#include <mpi.h>
//#include "mpi.h"


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

/********** Initialise MCMC structures *********/
/* Fill in samples from the prior distribution */
/* runState->algorithmParams must contain a variable "logLikelihoods" */
/* which contains a REAL8 array of likelihood values for the live */
/* points. */
/************************************************/
void initializeMCMC(LALInferenceRunState *runState)
{
	char help[]="\
	[--Niter] N\tNumber of iterations(2*10^6)\n\
	[--Nskip] n\tNumber of iterations between disk save(100)\n\
	[--tempMax T]\tHighest temperature for parallel tempering(40.0)\n\
	[--randomseed seed]\tRandom seed of sampling distribution\n\
        [--tdlike]\tCompute likelihood in the time domain\n";
	
	INT4 verbose=0,tmpi=0;
	unsigned int randomseed=0;
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
	runState->algorithm=&PTMCMCAlgorithm;
	runState->evolve=PTMCMCOneStep;
	//runState->evolve=&PTMCMCAdaptationOneStep;
	runState->proposal=&PTMCMCLALProposal;
	//runState->proposal=&PTMCMCLALSingleProposal;
	//runState->proposal=&PTMCMCLALAdaptationProposal;
	//runState->proposal=PTMCMCGaussianProposal;
	
	/* This is the LAL template generator for inspiral signals */
	
	ppt=getProcParamVal(commandLine,"--approx");
	if(ppt){
		/*if(strstr(ppt->value,"SpinTaylor")) {
			runState->template=&templateLALSTPN;
			fprintf(stdout,"Template function called is \"templateLALSTPN\"\n");
		}
		else {
			runState->template=&templateLAL;
			fprintf(stdout,"Template function called is \"templateLAL\"\n");
		}*/
		if(strstr(ppt->value,"TaylorF2")) {
			runState->template=&templateLAL;
			fprintf(stdout,"Template function called is \"templateLAL\"\n");
		}
    else if(strstr(ppt->value,"35phase_25amp")) {
      runState->template=&template3525TD;
			fprintf(stdout,"Template function called is \"template3525TD\"\n");
    }
		else {
			runState->template=&templateLALGenerateInspiral;
			fprintf(stdout,"Template function called is \"templateLALGenerateInspiral\"\n");
		}
		
	}
	else {runState->template=&templateLAL;}

        if (getProcParamVal(commandLine,"--tdlike")) {
          fprintf(stderr, "Computing likelihood in the time domain.\n");
          runState->likelihood=&TimeDomainLogLikelihood;
        } else if (getProcParamVal(commandLine, "--zeroLogLike")) {
          /* Use zero log(L) */
          runState->likelihood=&ZeroLogLikelihood;
        } else {
          runState->likelihood=&UndecomposedFreqDomainLogLikelihood;
        }

	/* runState->likelihood=&FreqDomainLogLikelihood; */
	/* runState->likelihood=&UndecomposedFreqDomainLogLikelihood; */
        /* runState->likelihood=&TimeDomainLogLikelihood; */
	//runState->likelihood=&UnityLikelihood;
	//runState->likelihood=GaussianLikelihood;
	//runState->prior=&PTUniformLALPrior;
	//runState->prior=&LALInferenceInspiralPrior;
	runState->prior=&LALInferenceInspiralPriorNormalised;
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
		//fprintf(stderr,"Error, must specify iteration number\n");
		//MPI_Finalize();
		//exit(1);
		tmpi=20000000;
	}
	addVariable(runState->algorithmParams,"Niter",&tmpi, INT4_t,PARAM_FIXED);
	
	printf("set iteration number between disk save.\n");
	/* Number of live points */
	ppt=getProcParamVal(commandLine,"--Nskip");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		//fprintf(stderr,"Error, must specify iteration number\n");
		//MPI_Finalize();
		//exit(1);
		tmpi=100;
	}
	addVariable(runState->algorithmParams,"Nskip",&tmpi, INT4_t,PARAM_FIXED);
	
	printf("set highest temperature.\n");
	/* Maximum temperature of the temperature ladder */
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
	fprintf(stdout, " initialize(): random seed: %u\n", randomseed);
	addVariable(runState->algorithmParams,"random_seed",&randomseed, UINT4_t,PARAM_FIXED);
	gsl_rng_set(runState->GSLrandom, randomseed);

        
        /* Now make sure that everyone is running with un-correlated
           jumps!  We re-seed rank i process with the ith output of
           the RNG stream from the rank 0 process. Otherwise the
           random stream is the same across all processes. */
        INT4 i;
        for (i = 0; i < MPIrank; i++) {
          randomseed = gsl_rng_get(runState->GSLrandom);
        }
        gsl_rng_set(runState->GSLrandom, randomseed);
	
	return;
	
}

static INT4 readSquareMatrix(gsl_matrix *m, UINT4 N, FILE *inp) {
  UINT4 i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      REAL8 value;
      INT4 nread;

      nread = fscanf(inp, " %lg ", &value);

      if (nread != 1) {
        fprintf(stderr, "Cannot read from matrix file (in %s, line %d)\n",
                __FILE__, __LINE__);
        exit(1);
      }

      gsl_matrix_set(m, i, j, value);
    }
  }

  return 0;
}


/* Setup the variables to control template generation */
/* Includes specification of prior ranges */

void initVariables(LALInferenceRunState *state)
{
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimInspiralTable *injTable=NULL;
	LALVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALVariables));
	LALVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	ProcessParamsTable *ppt=NULL;
	INT4 AmpOrder=0;
	LALPNOrder PhaseOrder=LAL_PNORDER_TWO;
	Approximant approx=TaylorF2;
	//INT4 numberI4 = TaylorF2;
	//INT4 numberI4 = TaylorT3;
	//INT4 approx=TaylorF2;
	LALInferenceApplyTaper bookends = INFERENCE_TAPER_NONE;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
	REAL8 Dmin=1.0;
	REAL8 Dmax=100.0;
	REAL8 mcMin=1.0;
	REAL8 mcMax=15.3;
	//REAL8 logmcMax,logmcMin;
	REAL8 mMin=1.0,mMax=30.0;
	REAL8 MTotMax=35.0;
	REAL8 etaMin=0.0312;
	REAL8 etaMax=0.25;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax;//,tmpVal;
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
	
	memset(currentParams,0,sizeof(LALVariables));
	
	char help[]="\
	[--injXML injections.xml]\tInjection XML file to use\
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
        [--covarianceMatrix file]\tFind the Cholesky decomposition of the covariance matrix for jumps in file";
	
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
		//printf("%d\n",approx);
		if(strstr(ppt->value,"TaylorF2")) {approx=TaylorF2;}//numberI4 = TaylorF2;}		LALGetApproximantFromString DOES NOT HAVE TaylorF2 !!!!!!
		//if(strstr(ppt->value,"TaylorT3")) {approx=TaylorT3;}//numberI4 = TaylorT3;}
		//if(strstr(ppt->value,"SpinTaylor")) {approx=SpinTaylor;}//numberI4 = SpinTaylor;}
		fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",approx,PhaseOrder);
		//fprintf(stdout,"Templates will run using Approximant %i, phase order %i\n",numberI4,PhaseOrder);
	}

        /* This flag was added to account for the broken Big Dog
           injection, which had the opposite sign in H and L compared
           to Virgo. */
        if (getProcParamVal(commandLine, "--crazyInjectionHLSign")) {
          INT4 flag = 1;
          addVariable(currentParams, "crazyInjectionHLSign", &flag, INT4_t, PARAM_FIXED);
        } else {
          INT4 flag = 0;
          addVariable(currentParams, "crazyInjectionHLSign", &flag, INT4_t, PARAM_FIXED);
        }

	/* Over-ride taper if specified */
	ppt=getProcParamVal(commandLine,"--taper");
	if(ppt){
		if(strstr(ppt->value,"STARTEND")) bookends=INFERENCE_TAPER_STARTEND;
		if(strstr(ppt->value,"STARTONLY")) bookends=INFERENCE_TAPER_START;
		if(strstr(ppt->value,"ENDONLY")) bookends=INFERENCE_TAPER_END;
		if(strstr(ppt->value,"RING")) bookends=INFERENCE_RING;
		if(strstr(ppt->value,"SMOOTH")) bookends=INFERENCE_SMOOTH;
	}

	/* Over-ride end time if specified */
	ppt=getProcParamVal(commandLine,"--trigtime");
	if(ppt){
		endtime=atof(ppt->value);
	}

	/* Over-ride chirp mass if specified */
	ppt=getProcParamVal(commandLine,"--mc");
	if(ppt){
		start_mc=atof(ppt->value);
	}
	
	/* Over-ride eta if specified */
	ppt=getProcParamVal(commandLine,"--eta");
	if(ppt){
		start_eta=atof(ppt->value);
	}
	
	/* Over-ride phase if specified */
	ppt=getProcParamVal(commandLine,"--phi");
	if(ppt){
		start_phase=atof(ppt->value);
	}
	
	/* Over-ride inclination if specified */
	ppt=getProcParamVal(commandLine,"--iota");
	if(ppt){
		start_iota=atof(ppt->value);
	}

        /* Over-ride distance if specified */
        ppt=getProcParamVal(commandLine,"--dist");
        if (ppt) {
          start_dist = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--ra");
        if (ppt) {
          start_ra = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--dec");
        if (ppt) {
          start_dec = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--psi");
        if (ppt) {
          start_psi = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--a1");
        if (ppt) {
          start_a_spin1 = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--theta1");
        if (ppt) {
          start_theta_spin1 = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--phi1");
        if (ppt) {
          start_phi_spin1 = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--a2");
        if (ppt) {
          start_a_spin2 = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--theta2");
        if (ppt) {
          start_theta_spin2 = atof(ppt->value);
        }

        ppt=getProcParamVal(commandLine,"--phi2");
        if (ppt) {
          start_phi_spin2 = atof(ppt->value);
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
	ppt=getProcParamVal(commandLine,"--compmin");
	if(ppt)	mMin=atof(ppt->value);
	addVariable(priorArgs,"component_min",&mMin,REAL8_t,PARAM_FIXED);
	ppt=getProcParamVal(commandLine,"--compmax");
	if(ppt)	mMax=atof(ppt->value);
	addVariable(priorArgs,"component_max",&mMax,REAL8_t,PARAM_FIXED);
	ppt=getProcParamVal(commandLine,"--MTotMax");
	if(ppt)	MTotMax=atof(ppt->value);
	addVariable(priorArgs,"MTotMax",&MTotMax,REAL8_t,PARAM_FIXED);
	
	
	printf("Read end time %f\n",endtime);

	addVariable(currentParams, "LAL_APPROXIMANT", &approx,        INT4_t, PARAM_FIXED);
	//addVariable(currentParams, "LAL_APPROXIMANT", &numberI4,        INT4_t, PARAM_FIXED);
	//numberI4 = LAL_PNORDER_TWO;
    addVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        INT4_t, PARAM_FIXED);	
	//addVariable(currentParams, "LAL_PNORDER",     &numberI4,        INT4_t, PARAM_FIXED);
	
	ppt=getProcParamVal(commandLine,"--taper");
	if(ppt){
		addVariable(currentParams, "INFERENCE_TAPER",     &bookends,        INT4_t, PARAM_FIXED);
	}
  ppt=getProcParamVal(commandLine,"--newswitch");
  int newswitch=0;
  if(ppt){
    newswitch=1;
    addVariable(currentParams, "newswitch", &newswitch, INT4_t, PARAM_FIXED);
  }
	/* Set up the variable parameters */
	//tmpVal=4.82+gsl_ran_gaussian(GSLrandom,0.025);//log(mcMin+(mcMax-mcMin)/2.0);
	//tmpVal=7.86508;
	ppt=getProcParamVal(commandLine,"--fixMc");
	if(ppt){
		addVariable(currentParams, "chirpmass",    &start_mc,    REAL8_t,	PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"chirpmass fixed and set to %f\n",start_mc);
	}else{
	    addVariable(currentParams, "chirpmass",    &start_mc,    REAL8_t,	PARAM_LINEAR);
    }
	addMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		REAL8_t);
	//addVariable(currentParams,"logmc",&tmpVal, REAL8_t, PARAM_LINEAR);
	//logmcMin=log(mcMin); logmcMax=log(mcMax);
	//addMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		REAL8_t);

	//tmpVal=0.244;
	//tmpVal=0.03+gsl_rng_uniform(GSLrandom)*(0.25-0.03);
	//tmpVal=0.18957;
	ppt=getProcParamVal(commandLine,"--fixEta");
	if(ppt){
	    addVariable(currentParams, "massratio",       &start_eta,             REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"eta fixed and set to %f\n",start_eta);
	}else{
	    addVariable(currentParams, "massratio",       &start_eta,             REAL8_t, PARAM_LINEAR);
	}
    addMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	REAL8_t);
	
	tmpMin=endtime-dt; tmpMax=endtime+dt;

        /* Set up start time. */
        ppt=getProcParamVal(commandLine, "--time");
        if (ppt) {
          /* User has specified start time. */
          timeParam = atof(ppt->value);
        } else {
          timeParam = endtime+gsl_ran_gaussian(GSLrandom,0.01);          
        }

	ppt=getProcParamVal(commandLine,"--fixTime");
	if(ppt){
	    addVariable(currentParams, "time",            &timeParam   ,           REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"time fixed and set to %f\n",timeParam);
	}else{
	    addVariable(currentParams, "time",            &timeParam   ,           REAL8_t, PARAM_LINEAR); 
	}
	addMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   REAL8_t);	

	//tmpVal=1.5;
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.89954;
	ppt=getProcParamVal(commandLine,"--fixPhi");
	if(ppt){
		addVariable(currentParams, "phase",           &start_phase,        REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"phase fixed and set to %f\n",start_phase);
	}else{
	    addVariable(currentParams, "phase",           &start_phase,        REAL8_t, PARAM_CIRCULAR);
	}
	addMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   REAL8_t);
	
	//tmpVal=5.8287;
	//tmpVal=8.07955+gsl_ran_gaussian(GSLrandom,1.1);
	//Dmin+(Dmax-Dmin)/2.0;
	//tmpVal=46.92314;
	ppt=getProcParamVal(commandLine,"--fixDist");
	if(ppt){
	     addVariable(currentParams,"distance", &start_dist, REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"distance fixed and set to %f\n",start_dist);
	}else{	
	     addVariable(currentParams,"distance", &start_dist, REAL8_t, PARAM_LINEAR);
	}
	addMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   REAL8_t);

	
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=4.5500;//1.0;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.34650;
	ppt=getProcParamVal(commandLine,"--fixRa");
	if(ppt){
		 addVariable(currentParams, "rightascension",  &start_ra,      REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"R.A. fixed and set to %f\n",start_ra);
    }else{
	     addVariable(currentParams, "rightascension",  &start_ra,      REAL8_t, PARAM_CIRCULAR);
	}
	addMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	//tmpVal=1.0759;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=-0.90547;
	ppt=getProcParamVal(commandLine,"--fixDec");
	if(ppt){
		addVariable(currentParams, "declination",     &start_dec,     REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"declination fixed and set to %f\n",start_dec);
	}else{
	    addVariable(currentParams, "declination",     &start_dec,     REAL8_t, PARAM_CIRCULAR);
	}
	addMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   REAL8_t);
    
	tmpMin=0.0; tmpMax=LAL_PI;
	//tmpVal=0.2000;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=0.64546;
	ppt=getProcParamVal(commandLine,"--fixPsi");
	if(ppt){
	     addVariable(currentParams, "polarisation",    &start_psi,     REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"polarisation fixed and set to %f\n",start_psi);
	}else{	
	     addVariable(currentParams, "polarisation",    &start_psi,     REAL8_t, PARAM_CIRCULAR);
	}
	addMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   REAL8_t);
	
	tmpMin=0.0; tmpMax=LAL_PI;
  
  ppt=getProcParamVal(commandLine,"--max-iota");
  if (ppt) {
    tmpMax = atof(ppt->value);
  }
	//tmpVal=0.9207;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=2.86094;

	ppt=getProcParamVal(commandLine,"--fixIota");
	if(ppt){
		addVariable(currentParams, "inclination",     &start_iota,            REAL8_t, PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"iota fixed and set to %f\n",start_iota);
	}else{
 	    addVariable(currentParams, "inclination",     &start_iota,            REAL8_t, PARAM_CIRCULAR);
	}
	addMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   REAL8_t);
	
	ppt=getProcParamVal(commandLine, "--noSpin");
	if((approx==SpinTaylor ||approx==SpinTaylorFrameless) && !ppt){
		

      ppt=getProcParamVal(commandLine, "--spinAligned");
      if(ppt) tmpMin=-1.0;
      else tmpMin=0.0;
      tmpMax=1.0;
			ppt=getProcParamVal(commandLine,"--fixA1");
			if(ppt){
			    addVariable(currentParams, "a_spin1",     &start_a_spin1,            REAL8_t, PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
			}else{
				addVariable(currentParams, "a_spin1",     &start_a_spin1,            REAL8_t, PARAM_LINEAR);
			}
			addMinMaxPrior(priorArgs, "a_spin1",     &tmpMin, &tmpMax,   REAL8_t);
				
			ppt=getProcParamVal(commandLine, "--spinAligned");
			if(ppt) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
			else {
				tmpMin=0.0; tmpMax=LAL_PI;
				ppt=getProcParamVal(commandLine,"--fixTheta1");
				if(ppt){
				    addVariable(currentParams, "theta_spin1",     &start_theta_spin1,            REAL8_t, PARAM_FIXED);
					if(MPIrank==0) fprintf(stdout,"theta 1 fixed and set to %f\n",start_theta_spin1);
				}else{
				    addVariable(currentParams, "theta_spin1",     &start_theta_spin1,            REAL8_t, PARAM_CIRCULAR);
				}
				addMinMaxPrior(priorArgs, "theta_spin1",     &tmpMin, &tmpMax,   REAL8_t);
		
				tmpMin=0.0; tmpMax=LAL_TWOPI;
				ppt=getProcParamVal(commandLine,"--fixPhi1");
				if(ppt){
					addVariable(currentParams, "phi_spin1",     &start_phi_spin1,            REAL8_t, PARAM_FIXED);
					if(MPIrank==0) fprintf(stdout,"phi 1 fixed and set to %f\n",start_phi_spin1);
				}else{
					addVariable(currentParams, "phi_spin1",     &start_phi_spin1,            REAL8_t, PARAM_CIRCULAR);
				}
				addMinMaxPrior(priorArgs, "phi_spin1",     &tmpMin, &tmpMax,   REAL8_t);
			}
		ppt=getProcParamVal(commandLine, "--singleSpin");
		if(ppt) fprintf(stdout,"Running with first spin set to 0\n");
		else {
    ppt=getProcParamVal(commandLine, "--spinAligned");
    if(ppt) tmpMin=-1.0;
    else tmpMin=0.0;
    tmpMax=1.0;
		ppt=getProcParamVal(commandLine,"--fixA2");
		if(ppt){
			addVariable(currentParams, "a_spin2",     &start_a_spin2,            REAL8_t, PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
		}else{
			addVariable(currentParams, "a_spin2",     &start_a_spin2,            REAL8_t, PARAM_LINEAR);
		}
		addMinMaxPrior(priorArgs, "a_spin2",     &tmpMin, &tmpMax,   REAL8_t);
	
		ppt=getProcParamVal(commandLine, "--spinAligned");
		if(ppt) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
		else {
			tmpMin=0.0; tmpMax=LAL_PI;
			ppt=getProcParamVal(commandLine,"--fixTheta2");
			if(ppt){
				addVariable(currentParams, "theta_spin2",     &start_theta_spin2,            REAL8_t, PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"theta spin 2 fixed and set to %f\n",start_theta_spin2);
			}else{
				addVariable(currentParams, "theta_spin2",     &start_theta_spin2,            REAL8_t, PARAM_CIRCULAR);
			}
			addMinMaxPrior(priorArgs, "theta_spin2",     &tmpMin, &tmpMax,   REAL8_t);
		
			tmpMin=0.0; tmpMax=LAL_TWOPI;
			ppt=getProcParamVal(commandLine,"--fixPhi2");
			if(ppt){
				addVariable(currentParams, "phi_spin2",     &start_phi_spin2,            REAL8_t, PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"phi 2 fixed and set to %f\n",start_phi_spin2);
			}else{
				addVariable(currentParams, "phi_spin2",     &start_phi_spin2,            REAL8_t, PARAM_CIRCULAR);
			}
			addMinMaxPrior(priorArgs, "phi_spin2",     &tmpMin, &tmpMax,   REAL8_t);
		}
	}
	}
  ppt=getProcParamVal(commandLine, "--spinAligned");
	if(approx==TaylorF2 && ppt){
		
    tmpMin=-1.0; tmpMax=1.0;
		ppt=getProcParamVal(commandLine,"--fixA1");
		if(ppt){
			addVariable(currentParams, "spin1",     &start_a_spin1,            REAL8_t, PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
		}else{
			addVariable(currentParams, "spin1",     &start_a_spin1,            REAL8_t, PARAM_LINEAR);
		}
    addMinMaxPrior(priorArgs, "spin1",     &tmpMin, &tmpMax,   REAL8_t);
		
		tmpMin=-1.0; tmpMax=1.0;
		ppt=getProcParamVal(commandLine,"--fixA2");
		if(ppt){
			addVariable(currentParams, "spin2",     &start_a_spin2,            REAL8_t, PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
		}else{
			addVariable(currentParams, "spin2",     &start_a_spin2,            REAL8_t, PARAM_LINEAR);
		}
		addMinMaxPrior(priorArgs, "spin2",     &tmpMin, &tmpMax,   REAL8_t);
    
	}
  
        /* Make sure that our initial value is within the
           prior-supported volume. */
        LALInferenceCyclicReflectiveBound(currentParams, priorArgs);

        /* Init covariance matrix, if specified.  The given file
           should contain the desired covariance matrix for the jump
           proposal, in row-major (i.e. C) order. */
        ppt=getProcParamVal(commandLine, "--covarianceMatrix");
        if (ppt) {
          FILE *inp = fopen(ppt->value, "r");
          UINT4 N = getVariableDimensionNonFixed(currentParams);
          gsl_matrix *covM = gsl_matrix_alloc(N,N);
          gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
          REAL8Vector *sigmaVec = XLALCreateREAL8Vector(N);
          UINT4 i;


          if (readSquareMatrix(covM, N, inp)) {
            fprintf(stderr, "Error reading covariance matrix (in %s, line %d)\n",
                    __FILE__, __LINE__);
            exit(1);
          }

          gsl_matrix_memcpy(covCopy, covM);

          for (i = 0; i < N; i++) {
            sigmaVec->data[i] = sqrt(gsl_matrix_get(covM, i, i)); /* Single-parameter sigma. */
          }

          addVariable(state->proposalArgs, SIGMAVECTORNAME, &sigmaVec, REAL8Vector_t, PARAM_FIXED);

          /* Set up eigenvectors and eigenvalues. */
          gsl_matrix *eVectors = gsl_matrix_alloc(N,N);
          gsl_vector *eValues = gsl_vector_alloc(N);
          REAL8Vector *eigenValues = XLALCreateREAL8Vector(N);
          gsl_eigen_symmv_workspace *ws = gsl_eigen_symmv_alloc(N);
          int gsl_status;

          if ((gsl_status = gsl_eigen_symmv(covCopy, eValues, eVectors, ws)) != GSL_SUCCESS) {
            fprintf(stderr, "Error in gsl_eigen_symmv (in %s, line %d): %d: %s\n",
                    __FILE__, __LINE__, gsl_status, gsl_strerror(gsl_status));
            exit(1);
          }

          for (i = 0; i < N; i++) {
            eigenValues->data[i] = gsl_vector_get(eValues,i);
          }

          addVariable(state->proposalArgs, "covarianceEigenvectors", &eVectors, gslMatrix_t, PARAM_FIXED);
          addVariable(state->proposalArgs, "covarianceEigenvalues", &eigenValues, REAL8Vector_t, PARAM_FIXED);

          fprintf(stdout, "Jumping with correlated jumps in %d dimensions from file %s.\n",
                  N, ppt->value);

          fclose(inp);
          gsl_eigen_symmv_free(ws);
          gsl_matrix_free(covCopy);
          gsl_vector_free(eValues);
        }

        /* Differential Evolution? */
        ppt=getProcParamVal(commandLine, "--differential-evolution");
        if (ppt) {
          FILE *dePtsFile = fopen(ppt->value, "r");
          
          if (!dePtsFile) {
            fprintf(stderr, "Could not open differential evolution file (%s, line %d).\n",
                    __FILE__, __LINE__);
            exit(1);
          } else {
            printf("Using differential evolution jumps from file %s\n", ppt->value);
          }
          
          char **headers = getHeaderLine(dePtsFile);
          size_t maxDePtsLen = 1;
          size_t dePtsLen = 1;
          LALVariables **dePts = malloc(sizeof(LALVariables *));
          
          while (!feof(dePtsFile)) {
            dePts[dePtsLen-1] = malloc(sizeof(LALVariables));
            dePts[dePtsLen-1]->head = NULL;
            dePts[dePtsLen-1]->dimension = 0;
            
            processParamLine(dePtsFile, headers, dePts[dePtsLen-1]);
            
            dePtsLen++;
            if (dePtsLen > maxDePtsLen) {
              /* Extend. */
              maxDePtsLen *= 2;
              dePts = realloc(dePts, maxDePtsLen*sizeof(LALVariables *));
            }
          }
          
          dePts = realloc(dePts, dePtsLen*sizeof(LALVariables *));
          
          state->differentialPoints = dePts;
          state->differentialPointsLength = dePtsLen;
          
          fclose(dePtsFile);
          free(headers); /* Reclaim some (but not all) the memory from
                            header.  (The individual names must stick
                            around to be keys in the LALVariables
                            structure.) */
        } else {
          state->differentialPoints = NULL;
          state->differentialPointsLength = 0;
        }

        ppt=getProcParamVal(commandLine, "--adapt");
        if (ppt) {
          fprintf(stdout, "Adapting single-param step sizes.\n");
          UINT4 N = (approx == SpinTaylor ? 15 : 9);
          if (!checkVariable(state->proposalArgs, SIGMAVECTORNAME)) {
            /* We need a sigma vector for adaptable jumps. */
            REAL8Vector *sigmaVec = XLALCreateREAL8Vector(N);
            UINT4 i = 0;
            
            for (i = 0; i < N; i++) {
              sigmaVec->data[i] = 1e-4;
            }
            
            addVariable(state->proposalArgs, SIGMAVECTORNAME, &sigmaVec, REAL8Vector_t, PARAM_FIXED);
          }
          REAL8Vector *avgPaccept = XLALCreateREAL8Vector(N);
          UINT4 i;

          for (i = 0; i < N; i++) {
            avgPaccept->data[i] = 0.0;
          }
          
          addVariable(state->proposalArgs, "adaptPacceptAvg", &avgPaccept, REAL8Vector_t, PARAM_FIXED);
        }

        INT4 adaptableStep = 0;
        addVariable(state->proposalArgs, "adaptableStep", &adaptableStep, INT4_t, PARAM_OUTPUT);

        INT4 varNumber = 0;
        addVariable(state->proposalArgs, "proposedVariableNumber", &varNumber, INT4_t, PARAM_OUTPUT);

        INT4 sigmaVecNumber = 0;
        addVariable(state->proposalArgs, "proposedSigmaNumber", &sigmaVecNumber, INT4_t, PARAM_OUTPUT);

        REAL8 tau = 1e3;
        addVariable(state->proposalArgs, "adaptTau", &tau, REAL8_t, PARAM_OUTPUT);

        ppt = getProcParamVal(commandLine, "--adaptTau");
        if (ppt) {
          tau = atof(ppt->value);
          fprintf(stdout, "Setting adapt tau = %g.\n", tau);
          setVariable(state->proposalArgs, "adaptTau", &tau);
        }
	
	return;
}




int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

  if (MPIrank == 0) fprintf(stdout," ========== LALInference_MCMCMPI ==========\n");

	LALInferenceRunState *runState;
	ProcessParamsTable *procParams=NULL;
	ProcessParamsTable *ppt=NULL;
	char *infileName;
	infileName = (char*)calloc(99,sizeof(char*));
	char str [999];
	FILE * infile;
	int n;
	char * pch;
	int fileargc = 1;
	char *fileargv[99];
	//char **fileargv[99][999] = NULL;
	char buffer [99];
	
//	for (i=0; i<argc; i++) {
//		printf("%s\n",argv[i]);
//	}
//	printf("%d\n",argc);
	
	/* Read command line and parse */
	procParams=parseCommandLine(argc,argv);
	
	ppt=getProcParamVal(procParams,"--continue-run");
	if (ppt) {
		infileName = ppt->value;
		infile = fopen(infileName,"r");
		if (infile==NULL) {fprintf(stderr,"Cannot read %s/n",infileName); exit (1);}
		n=sprintf(buffer,"lalinference_mcmcmpi_from_file_%s",infileName);
		fileargv[0] = (char*)calloc((n+1),sizeof(char*));
		fileargv[0] = buffer;
		fgets(str, 999, infile);
		fgets(str, 999, infile);
		fclose(infile);
		pch = strtok (str," ");
		while (pch != NULL)
		{
			if(strcmp(pch,"Command")!=0 && strcmp(pch,"line:")!=0)
			{
				n = strlen(pch);
				fileargv[fileargc] = (char*)calloc((n+1),sizeof(char*));
				fileargv[fileargc] = pch;
				fileargc++;
				if(fileargc>=99) {fprintf(stderr,"Too many arguments in file %s\n",infileName); exit (1);}
			}
			pch = strtok (NULL, " ");

		}
		//pch = strstr (fileargv[fileargc-1],"\n");
		//strncpy (pch,"",1);
		fileargv[fileargc-1][strlen(fileargv[fileargc-1])-1]='\0'; //in order to get rid of the '\n' than fgets returns when reading the command line.

		//for (i=0; i<fileargc; i++) {
		//	printf("%s\n",fileargv[i]);
		//}
		//printf("%d\n",fileargc);

		procParams=parseCommandLine(fileargc,fileargv);
		//ppt = getProcParamVal(procParams, "--randomseed");
		//if (ppt == NULL){
		//	ProcessParamsTable *this = procParams;
		//	ProcessParamsTable *previous = procParams;
		//	while (this != NULL) {
		//		previous = this;
		//		this = this->next;
		//	}
		//	previous->next = (ProcessParamsTable*) calloc(1, sizeof(ProcessParamsTable));
		//	previous = previous->next;
		//	strcpy(previous->program, fileargv[0]);
		//	strcpy(previous->param, "--randomseed");
		//	strcpy(previous->type, "string");
		//	strcpy(previous->value, "11111111");
		//}
	}
	
	

	/* initialise runstate based on command line */
	/* This includes reading in the data */
	/* And performing any injections specified */
	/* And allocating memory */
	runState = initialize(procParams);
  
	/* Set up structures for MCMC */
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
