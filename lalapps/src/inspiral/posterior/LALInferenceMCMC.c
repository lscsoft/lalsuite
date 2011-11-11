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


#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>

#include <mpi.h>


int MPIrank, MPIsize;

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMCMC(LALInferenceRunState *runState);
void initVariables(LALInferenceRunState *state);




LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
	LALInferenceRunState *irs=NULL;
	LALInferenceIFOData *ifoPtr, *ifoListStart;

	MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
	
	irs = calloc(1, sizeof(LALInferenceRunState));
	/* read data from files: */
	fprintf(stdout, " ==== LALInferenceReadData(): started. ====\n");
	irs->commandLine=commandLine;
	irs->data = LALInferenceReadData(commandLine);
	/* (this will already initialise each LALInferenceIFOData's following elements:  */
	/*     fLow, fHigh, detector, timeToFreqFFTPlan, freqToTimeFFTPlan,     */
	/*     window, oneSidedNoisePowerSpectrum, timeDate, freqData         ) */
	fprintf(stdout, " ==== LALInferenceReadData(): finished. ====\n");
	if (irs->data != NULL) {
		fprintf(stdout, " ==== initialize(): successfully read data. ====\n");
		
		fprintf(stdout, " ==== LALInferenceInjectInspiralSignal(): started. ====\n");
		LALInferenceInjectInspiralSignal(irs->data,commandLine);
		fprintf(stdout, " ==== LALInferenceInjectInspiralSignal(): finished. ====\n");
		
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
				ifoPtr->modelParams = calloc(1, sizeof(LALInferenceVariables));
			}
			ifoPtr = ifoPtr->next;
		}
		irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
		printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else{
		fprintf(stdout, " initialize(): no data read.\n");
    irs = NULL;
    return(irs);
	}
	
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
(--Niter N)                     Number of iterations(2*10^6)\n\
(--Nskip n)                     Number of iterations between disk save(100)\n\
(--tempMax T)                   Highest temperature for parallel tempering(40.0)\n\
(--randomseed seed)             Random seed of sampling distribution(random)\n\
(--tdlike)                      Compute likelihood in the time domain\n\
(--rapidSkyLoc)                 Use rapid sky localization jump proposals\n";
	
  /* Print command line arguments if runState was not allocated */
	if(runState==NULL)
	{
		fprintf(stdout,"%s",help);
		return;
	}
  
  
	INT4 verbose=0,tmpi=0;
	unsigned int randomseed=0;
	REAL8 tempMax = 10.0;
	//REAL8 tmp=0;
	ProcessParamsTable *commandLine=runState->commandLine;
	ProcessParamsTable *ppt=NULL;
	FILE *devrandom;
	struct timeval tv;
	
	/* Print command line arguments if help requested */
	if(LALInferenceGetProcParamVal(runState->commandLine,"--help"))
	{
		fprintf(stdout,"%s",help);
    	runState->algorithm=&PTMCMCAlgorithm;
		return;
	}
	
	/* Initialise parameters structure */
	runState->algorithmParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	runState->priorArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
	runState->proposalArgs=XLALCalloc(1,sizeof(LALInferenceVariables));
	
	/* Set up the appropriate functions for the MCMC algorithm */
	runState->algorithm=&PTMCMCAlgorithm;
	runState->evolve=PTMCMCOneStep;

    ppt=LALInferenceGetProcParamVal(commandLine,"--rapidSkyLoc");
    if(ppt)
        runState->proposal=&LALInferenceRapidSkyLocProposal;
    else
        runState->proposal=&LALInferenceDefaultProposal;
	
	/* This is the LAL template generator for inspiral signals */
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
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
			runState->template=&LALInferenceTemplateLAL;
			fprintf(stdout,"Template function called is \"templateLAL\"\n");
		}
    else if(strstr(ppt->value,"35phase_25amp")) {
      runState->template=&LALInferenceTemplate3525TD;
			fprintf(stdout,"Template function called is \"template3525TD\"\n");
    }
    //else if(strstr(ppt->value,"PhenSpinTaylorRD")) {
    //    runState->template=&LALInferenceTemplatePSTRD;
    //        fprintf(stdout,"Template function called is \"templatePSTRD\"\n");
    //}
		else {
			runState->template=&LALInferenceTemplateLALGenerateInspiral;
      //runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
			fprintf(stdout,"Template function called is \"templateLALGenerateInspiral\"\n");
		}
		
	}
	else {runState->template=&LALInferenceTemplateLAL;}

        if (LALInferenceGetProcParamVal(commandLine,"--tdlike")) {
          fprintf(stderr, "Computing likelihood in the time domain.\n");
          runState->likelihood=&LALInferenceTimeDomainLogLikelihood;
        } else if (LALInferenceGetProcParamVal(commandLine, "--zeroLogLike")) {
          /* Use zero log(L) */
          runState->likelihood=&LALInferenceZeroLogLikelihood;
        } else if (LALInferenceGetProcParamVal(commandLine, "--analyticLogLike")) {
          /* Use zero log(L) */
          runState->likelihood=&LALInferenceAnalyticLogLikelihood;
        } else {
          runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
        }

	/* runState->likelihood=&FreqDomainLogLikelihood; */
	/* runState->likelihood=&UndecomposedFreqDomainLogLikelihood; */
        /* runState->likelihood=&TimeDomainLogLikelihood; */
	//runState->likelihood=&UnityLikelihood;
	//runState->likelihood=GaussianLikelihood;
	//runState->prior=&PTUniformLALPrior;
	//runState->prior=&LALInferenceInspiralPrior;
  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    runState->prior=&LALInferenceInspiralSkyLocPrior;
    
  }else{
    runState->prior=&LALInferenceInspiralPriorNormalised;
  }
	//runState->prior=PTUniformGaussianPrior;


	
	ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
	if(ppt) {
		verbose=1;
		LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose , LALINFERENCE_UINT4_t,
					LALINFERENCE_PARAM_FIXED);
	}
	
	printf("set iteration number.\n");
	/* Number of live points */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Niter");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		//fprintf(stderr,"Error, must specify iteration number\n");
		//MPI_Finalize();
		//exit(1);
		tmpi=20000000;
	}
	LALInferenceAddVariable(runState->algorithmParams,"Niter",&tmpi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
	
	printf("set iteration number between disk save.\n");
	/* Number of live points */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nskip");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		//fprintf(stderr,"Error, must specify iteration number\n");
		//MPI_Finalize();
		//exit(1);
		tmpi=100;
	}
	LALInferenceAddVariable(runState->algorithmParams,"Nskip",&tmpi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
	
	printf("set highest temperature.\n");
	/* Maximum temperature of the temperature ladder */
	ppt=LALInferenceGetProcParamVal(commandLine,"--tempMax");
	if(ppt){
		tempMax=strtod(ppt->value,(char **)NULL);
	}	
	LALInferenceAddVariable(runState->algorithmParams,"tempMax",&tempMax, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	printf("set random seed.\n");
	/* set up GSL random number generator: */
	gsl_rng_env_setup();
	runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
	/* (try to) get random seed from command line: */
	ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
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
	LALInferenceAddVariable(runState->algorithmParams,"random_seed",&randomseed, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
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
  
  char help[]="\
(--injXML injections.xml)       Injection XML file to use\n\
(--symMassRatio)                Run with symmetric mass ratio eta, instead of q=m2/m1\n\
(--mc-min mchirp)               Minimum chirp mass\n\
(--mc-max mchirp)               Maximum chirp mass\n\
(--eta-min etaMin)              Minimum eta\n\
(--eta-max etaMax)              Maximum eta\n\
(--q-min qMin)                  Minimum q\n\
(--q-max qMax)                  Maximum q\n\
(--dt time)                     Width of time prior, centred around trigger (0.1s)\n\
(--trigtime time)               Trigger time to use\n\
(--mc mchirp)                   Trigger chirpmass to use\n\
(--fixMc)                       Do not allow chirpmass to vary\n\
(--eta eta)                     Trigger eta to use\n\
(--q q)                         Trigger q to use\n\
(--fixEta)                      Do not allow mass ratio to vary\n\
(--fixQ)                        Do not allow mass ratio to vary\n\
(--phi phase)                   Trigger phase to use\n\
(--fixPhi)                      Do not allow phase to vary\n\
(--iota inclination)            Trigger inclination to use\n\
(--fixIota)                     Do not allow inclination to vary\n\
(--dist dist)                   Trigger distance\n\
(--fixDist)                     Do not allow distance to vary\n\
(--ra ra)                       Trigger RA\n\
(--fixRa)                       Do not allow RA to vary\n\
(--dec dec)                     Trigger declination\n\
(--fixDec)                      Do not allow declination to vary\n\
(--psi psi)                     Trigger psi\n\
(--fixPsi)                      Do not allow polarization to vary\n\
(--a1 a1)                       Trigger a1\n\
(--fixA1)                       Do not allow spin to vary\n\
(--theta1 theta1)               Trigger theta1\n\
(--fixTheta1)                   Do not allow spin 1 colatitude to vary\n\
(--phi1 phi1)                   Trigger phi1\n\
(--fixPhi1)                     Do not allow spin 1 longitude to vary\n\
(--a2 a2)                       Trigger a2\n\
(--fixA2)                       Do not allow spin 2 to vary\n\
(--theta2 theta2)               Trigger theta2\n\
(--fixTheta2)                   Do not allow spin 2 colatitude to vary\n\
(--phi2 phi2)                   Trigger phi2\n\
(--fixPhi2)                     Do not allow spin 2 longitude to vary\n\
(--time time)                   Waveform time (overrides random about trigtime)\n\
(--fixTime)                     Do not allow coalescence time to vary\n\
(--Dmin dist)                   Minimum distance in Mpc (1)\n\
(--Dmax dist)                   Maximum distance in Mpc (100)\n\
(--approximant Approximant)     Specify a template approximant to use, (default TaylorF2)\n\
(--order PNorder)               Specify a PN order to use, (default threePointFivePN)\n\
(--comp-min min)                Minimum component mass (1.0)\n\
(--comp-max max)                Maximum component mass (30.0)\n\
(--MTotMax max)                 Maximum total mass (35.0)\n\
(--covarianceMatrix file)       Find the Cholesky decomposition of the covariance matrix for jumps in file\n\
(--appendOutput fname)          Basename of the file to append outputs to\n";


  /* Print command line arguments if state was not allocated */
  if(state==NULL)
  {
    fprintf(stdout,"%s",help);
		return;
	}
  
  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
  {
	  fprintf(stdout,"%s",help);
		return;
	}
  
  
	LALStatus status;
	memset(&status,0,sizeof(status));
	SimInspiralTable *injTable=NULL;
	LALInferenceVariables *priorArgs=state->priorArgs;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	LALInferenceVariables *currentParams=state->currentParams;
	ProcessParamsTable *commandLine=state->commandLine;
	ProcessParamsTable *ppt=NULL;
	INT4 AmpOrder=0;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	Approximant approx=TaylorF2;
	//INT4 numberI4 = TaylorF2;
	//INT4 numberI4 = TaylorT3;
	//INT4 approx=TaylorF2;
	LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
    UINT4 event=0;
    UINT4 i=0;
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
    REAL8 qMin=mMin/mMax;
    REAL8 qMax=1.0;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax;//,tmpVal;
	gsl_rng * GSLrandom=state->GSLrandom;
	REAL8 endtime=0.0, timeParam=0.0;
	REAL8 start_mc			=4.82+gsl_ran_gaussian(GSLrandom,0.025);
	REAL8 start_eta			=etaMin+gsl_rng_uniform(GSLrandom)*(etaMax-etaMin);
	REAL8 start_q           =qMin+gsl_rng_uniform(GSLrandom)*(qMin-qMax);
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
	
	memset(currentParams,0,sizeof(LALInferenceVariables));
	
  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    MTotMax=20.0;
    mMin=1.0;
    mMax=15.0;
    qMin=mMin/mMax;
    Dmin=10.0;
    Dmax=40.0;
    REAL8 densityVNR=1000.0;
    LALInferenceAddVariable(state->priorArgs,"densityVNR", &densityVNR , LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }
  
	/* Read injection XML file for parameters if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--injXML");
	if(ppt){
		SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
		if(!injTable){
			fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
			MPI_Finalize();
			exit(1);
		}
        ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
            event= atoi(ppt->value);
            fprintf(stderr,"Reading event %d from file\n",event);
            i=0;
            while(i<event) {i++; injTable=injTable->next;} /* select event */
		
            endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
            AmpOrder=injTable->amp_order;
            LALGetOrderFromString(&status,injTable->waveform,&PhaseOrder);
            LALGetApproximantFromString(&status,injTable->waveform,&approx);
        }
	}	
	
	/* Over-ride approximant if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
	if(ppt){
		if ( ! strcmp( "GeneratePPN", ppt->value ) )
		{
			approx = GeneratePPN;
		}
		else if ( ! strcmp( "TaylorT1", ppt->value ) )
		{
			approx = TaylorT1;
		}
		else if ( ! strcmp( "TaylorT2", ppt->value ) )
		{
			approx = TaylorT2;
		}
		else if ( ! strcmp( "TaylorT3", ppt->value ) )
		{
			approx = TaylorT3;
		}
		else if ( ! strcmp( "TaylorT4", ppt->value ) )
		{
			approx = TaylorT4;
		}
		else if ( ! strcmp( "TaylorF1", ppt->value ) )
		{
			approx = TaylorF1;
		}
		else if ( ! strcmp( "TaylorF2", ppt->value ) )
		{
			approx = TaylorF2;
		}
		else if ( ! strcmp( "EOB", ppt->value ) )
		{
			approx = EOB;
		}
		else if ( ! strcmp( "EOBNR", ppt->value ) )
		{
			approx = EOBNR;
		}
		else if ( ! strcmp( "EOBNRv2", ppt->value ) )
		{
			approx = EOBNRv2;
		}
		else if ( ! strcmp( "EOBNRv2HM", ppt->value ) )
		{
			approx = EOBNRv2HM;
		}
		else if ( ! strcmp( "SpinTaylor", ppt->value ) )
		{
			approx = SpinTaylor;
		}
		else if ( ! strcmp( "SpinTaylorT3", ppt->value ) )
		{
			approx = SpinTaylorT3;
		}
		else if ( ! strcmp( "SpinQuadTaylor", ppt->value ) )
		{
			approx = SpinQuadTaylor;
		}
		else if ( ! strcmp( "SpinTaylorFrameless", ppt->value ) )
		{
			approx = SpinTaylorFrameless;
		}    
		else if ( ! strcmp( "PhenSpinTaylorRD", ppt->value ) )
		{
			approx = PhenSpinTaylorRD;
		}
		else if ( ! strcmp( "NumRel", ppt->value ) )
		{
			approx = NumRel;
		}
		else if ( ! strcmp( "IMRPhenomA", ppt->value ) )
		{
			approx = IMRPhenomA;
		}
		else if ( ! strcmp( "IMRPhenomB", ppt->value ) )
		{
			approx = IMRPhenomB;
		}
		else
		{
			fprintf( stderr, "invalid argument to --approximant\n"
					"unknown approximant %s specified: "
					"Approximant must be one of: GeneratePPN, TaylorT1, TaylorT2,\n"
					"TaylorT3, TaylorT4, TaylorF1, TaylorF2,  EOB, EOBNR, EOBNRv2, \n"
					"EOBNRv2HM, SpinTaylor, SpinTaylorT3, SpinQuadTaylor, SpinTaylorFrameless,\n"
					"PhenSpinTaylorRD, NumRel, IMRPhenomA, IMRPhenomB \n", ppt->value);
			exit( 1 );
		}
		fprintf(stdout,"Templates will run using Approximant %s (%u)\n",ppt->value,approx);
	}
	
        /* Over-ride PN order if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--order");
	if(ppt){
        if ( ! strcmp( "newtonian", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_NEWTONIAN;
        }
        else if ( ! strcmp( "oneHalfPN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_HALF;
        }
        else if ( ! strcmp( "onePN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_ONE;
        }
        else if ( ! strcmp( "onePointFivePN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_ONE_POINT_FIVE;
        }
        else if ( ! strcmp( "twoPN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_TWO;
        }
        else if ( ! strcmp( "twoPointFive", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_TWO_POINT_FIVE;
        }
        else if ( ! strcmp( "threePN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_THREE;
        }
        else if ( ! strcmp( "threePointFivePN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_THREE_POINT_FIVE;
        }
        else if ( ! strcmp( "pseudoFourPN", ppt->value ) )
        {
          PhaseOrder = LAL_PNORDER_PSEUDO_FOUR;
        }
        else
        {
          fprintf( stderr, "invalid argument to --order:\n"
              "unknown order specified: "
              "PN order must be one of: newtonian, oneHalfPN, onePN,\n"
              "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
              "threePointFivePN\n");
          exit( 1 );
        }
	fprintf(stdout,"Templates will be generated at %.1f PN order\n",((float)(PhaseOrder))/2.0);
	}

        /* This flag was added to account for the broken Big Dog
           injection, which had the opposite sign in H and L compared
           to Virgo. */
        if (LALInferenceGetProcParamVal(commandLine, "--crazyInjectionHLSign")) {
          INT4 flag = 1;
          LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
        } else {
          INT4 flag = 0;
          LALInferenceAddVariable(currentParams, "crazyInjectionHLSign", &flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
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

    /* Over-ride q if specified */
    ppt=LALInferenceGetProcParamVal(commandLine,"--q");
    if(ppt){
            start_q=atof(ppt->value);
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
	
	/* Over-ride time prior if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
	if(ppt){
		dt=atof(ppt->value);
	}
	
	/* Over-ride Distance min if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmin");
	if(ppt){
		logDmin=log(atof(ppt->value));
		Dmin=atof(ppt->value);
	}
	
	/* Over-ride Distance max if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Dmax");
	if(ppt){
		logDmax=log(atof(ppt->value));
		Dmax=atof(ppt->value);
	}
	
	/* Over-ride Mass priors if specified */
	ppt=LALInferenceGetProcParamVal(commandLine,"--mc-min");
	if(ppt) 
        {
		mcMin=atof(ppt->value);
		if (mcMin < 0) 
		{
			fprintf(stderr,"ERROR: Minimum value of mchirp must be > 0");
			exit(1);
		}
	}

	ppt=LALInferenceGetProcParamVal(commandLine,"--mc-max");
	if(ppt)	
        {
                mcMax=atof(ppt->value);
		if (mcMax <= 0) 
		{
			fprintf(stderr,"ERROR: Maximum value of mchirp must be > 0");
			exit(1);
		}
	}

	ppt=LALInferenceGetProcParamVal(commandLine,"--eta-min");
	if(ppt)	
        {
		etaMin=atof(ppt->value);
                if (etaMin < 0.0) 
		{
			fprintf(stderr,"ERROR: Minimum value of eta must be > 0");
			exit(1);
		}
	}

	ppt=LALInferenceGetProcParamVal(commandLine,"--eta-max");
	if(ppt)	
        {
		etaMax=atof(ppt->value);
                if (etaMax > 0.25 || etaMax <= 0.0) 
		{
			fprintf(stderr,"ERROR: Maximum value of eta must be between 0 and 0.25\n");
			exit(1);
		}
	}


	/* Over-ride component masses */
	ppt=LALInferenceGetProcParamVal(commandLine,"--comp-min");
	if(ppt)	mMin=atof(ppt->value);
	LALInferenceAddVariable(priorArgs,"component_min",&mMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

	ppt=LALInferenceGetProcParamVal(commandLine,"--comp-max");
	if(ppt)	mMax=atof(ppt->value);
	LALInferenceAddVariable(priorArgs,"component_max",&mMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

	ppt=LALInferenceGetProcParamVal(commandLine,"--MTotMax");
	if(ppt)	MTotMax=atof(ppt->value);
	LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

	ppt=LALInferenceGetProcParamVal(commandLine,"--q-min");
    if(ppt) 
     {
            qMin=atof(ppt->value);
            if (qMin <= 0.0 || qMin < mMin/mMax || qMin < mMin/(MTotMax-mMin) || qMin > 1.0) 
            {
                    fprintf(stderr,"ERROR: invalid qMin ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
                    exit(1);
            }
    }

    ppt=LALInferenceGetProcParamVal(commandLine,"--q-max");
    if(ppt) 
     {
            qMax=atof(ppt->value);
            if (qMax > 1.0 || qMax <= 0.0 || qMax < mMin/mMax || qMax < mMin/(MTotMax-mMin)) 
            {
                    fprintf(stderr,"ERROR: invalid qMax ( max{0,mMin/mMax,mMin/(MTotMax-mMin) < q < 1.0} )");
                    exit(1);
            }
    }


	
	printf("Read end time %f\n",endtime);

	LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
	//LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &numberI4,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
	//numberI4 = LAL_PNORDER_TWO;
    LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);	
	//LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &numberI4,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
	
	ppt=LALInferenceGetProcParamVal(commandLine,"--taper");
	if(ppt){
		LALInferenceAddVariable(currentParams, "LALINFERENCE_TAPER",     &bookends,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
	}
  ppt=LALInferenceGetProcParamVal(commandLine,"--newswitch");
  int newswitch=0;
  if(ppt){
    newswitch=1;
    LALInferenceAddVariable(currentParams, "newswitch", &newswitch, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  }
	/* Set up the variable parameters */
	//tmpVal=4.82+gsl_ran_gaussian(GSLrandom,0.025);//log(mcMin+(mcMax-mcMin)/2.0);
	//tmpVal=7.86508;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixMc");
	if(ppt){
		LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"chirpmass fixed and set to %f\n",start_mc);
	}else{
	    LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    }
	LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);
	//LALInferenceAddVariable(currentParams,"logmc",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	//logmcMin=log(mcMin); logmcMax=log(mcMax);
	//LALInferenceAddMinMaxPrior(priorArgs,	"logmc",	&logmcMin,	&logmcMax,		LALINFERENCE_REAL8_t);

	//tmpVal=0.244;
	//tmpVal=0.03+gsl_rng_uniform(GSLrandom)*(0.25-0.03);
	//tmpVal=0.18957;

    /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
	ppt=LALInferenceGetProcParamVal(commandLine,"--symMassRatio");
	if(ppt){
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixEta");
        if(ppt){
            LALInferenceAddVariable(currentParams, "massratio",       &start_eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
            if(MPIrank==0) fprintf(stdout,"eta fixed and set to %f\n",start_eta);
        }else{
            LALInferenceAddVariable(currentParams, "massratio",       &start_eta,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs,	"massratio",	&etaMin,	&etaMax,	LALINFERENCE_REAL8_t);
    }else{
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixQ");
        if(ppt){
            LALInferenceAddVariable(currentParams, "asym_massratio",       &start_q,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
            if(MPIrank==0) fprintf(stdout,"q fixed and set to %f\n",start_q);
        }else{
            LALInferenceAddVariable(currentParams, "asym_massratio",       &start_q,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
    }

	tmpMin=endtime-dt; tmpMax=endtime+dt;

        /* Set up start time. */
        ppt=LALInferenceGetProcParamVal(commandLine, "--time");
        if (ppt) {
          /* User has specified start time. */
          timeParam = atof(ppt->value);
        } else {
          timeParam = endtime+gsl_ran_gaussian(GSLrandom,0.01);          
        }

	ppt=LALInferenceGetProcParamVal(commandLine,"--fixTime");
	if(ppt){
	    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"time fixed and set to %f\n",timeParam);
	}else{
	    LALInferenceAddVariable(currentParams, "time",            &timeParam   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	}
	LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	

	//tmpVal=1.5;
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.89954;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi");
	if(ppt){
		LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"phase fixed and set to %f\n",start_phase);
	}else{
	    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	//tmpVal=5.8287;
	//tmpVal=8.07955+gsl_ran_gaussian(GSLrandom,1.1);
	//Dmin+(Dmax-Dmin)/2.0;
	//tmpVal=46.92314;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixDist");
	if(ppt){
	     LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"distance fixed and set to %f\n",start_dist);
	}else{	
	     LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);

	
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	//tmpVal=4.5500;//1.0;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=3.34650;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixRa");
	if(ppt){
		 LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"R.A. fixed and set to %f\n",start_ra);
    }else{
	     LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	//tmpVal=1.0759;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=-0.90547;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixDec");
	if(ppt){
		LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"declination fixed and set to %f\n",start_dec);
	}else{
	    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
	tmpMin=0.0; tmpMax=LAL_PI;
	//tmpVal=0.2000;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=0.64546;
	ppt=LALInferenceGetProcParamVal(commandLine,"--fixPsi");
	if(ppt){
	     LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"polarisation fixed and set to %f\n",start_psi);
	}else{	
	     LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	tmpMin=0.0; tmpMax=LAL_PI;
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--max-iota");
  if (ppt) {
    tmpMax = atof(ppt->value);
  }
	//tmpVal=0.9207;
	//tmpVal=tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
	//tmpVal=2.86094;

	ppt=LALInferenceGetProcParamVal(commandLine,"--fixIota");
	if(ppt){
		LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		if(MPIrank==0) fprintf(stdout,"iota fixed and set to %f\n",start_iota);
	}else{
 	    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	}
	LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
	if((approx==SpinTaylor || approx==SpinTaylorFrameless || approx==PhenSpinTaylorRD) && !ppt){
		

      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      if(ppt) tmpMin=-1.0;
      else tmpMin=0.0;
      tmpMax=1.0;
			ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
			if(ppt){
			    LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
			}else{
				LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			}
			LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
				
			ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
			if(ppt) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
			else {
				tmpMin=0.0; tmpMax=LAL_PI;
				ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta1");
				if(ppt){
				    LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					if(MPIrank==0) fprintf(stdout,"theta 1 fixed and set to %f\n",start_theta_spin1);
				}else{
				    LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
				}
				LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		
				tmpMin=0.0; tmpMax=LAL_TWOPI;
				ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi1");
				if(ppt){
					LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					if(MPIrank==0) fprintf(stdout,"phi 1 fixed and set to %f\n",start_phi_spin1);
				}else{
					LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
				}
				LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
			}
		ppt=LALInferenceGetProcParamVal(commandLine, "--singleSpin");
		if(ppt) fprintf(stdout,"Running with first spin set to 0\n");
		else {
    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    if(ppt) tmpMin=-1.0;
    else tmpMin=0.0;
    tmpMax=1.0;
		ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
		if(ppt){
			LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
		}else{
			LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
		ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
		if(ppt) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
		else {
			tmpMin=0.0; tmpMax=LAL_PI;
			ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta2");
			if(ppt){
				LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"theta spin 2 fixed and set to %f\n",start_theta_spin2);
			}else{
				LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			}
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		
			tmpMin=0.0; tmpMax=LAL_TWOPI;
			ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi2");
			if(ppt){
				LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
				if(MPIrank==0) fprintf(stdout,"phi 2 fixed and set to %f\n",start_phi_spin2);
			}else{
				LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			}
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		}
	}
	}
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
	if(approx==TaylorF2 && ppt){
		
    tmpMin=-1.0; tmpMax=1.0;
		ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
		if(ppt){
			LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
		}else{
			LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "spin1",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
		
		tmpMin=-1.0; tmpMax=1.0;
		ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
		if(ppt){
			LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
		}else{
			LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
		LALInferenceAddMinMaxPrior(priorArgs, "spin2",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
	}
  
  ppt=LALInferenceGetProcParamVal(commandLine, "--TaylorF2ppE");
	if(approx==TaylorF2 && ppt){
    
    REAL8 start_alpha, start_A, start_a, start_beta, start_B, start_b;
    
    tmpMin = -1000;
    tmpMax = 1000;
    start_alpha = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppealpha");
    if (ppt) {
      start_alpha = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppealpha");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE alpha fixed and set to %f\n",start_alpha);
		}else{
			LALInferenceAddVariable(currentParams, "ppealpha",     &start_alpha,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppealpha",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
    start_beta = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppebeta");
    if (ppt) {
      start_beta = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppebeta");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE beta fixed and set to %f\n",start_beta);
		}else{
			LALInferenceAddVariable(currentParams, "ppebeta",     &start_beta,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppebeta",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
    tmpMin = -3;
    tmpMax = 3;
    start_A = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeuppera");
    if (ppt) {
      start_A = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeuppera");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE A fixed and set to %f\n",start_A);
		}else{
			LALInferenceAddVariable(currentParams, "ppeuppera",     &start_A,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppeuppera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
    start_B = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppeupperb");
    if (ppt) {
      start_B = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppeupperb");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE B fixed and set to %f\n",start_B);
		}else{
			LALInferenceAddVariable(currentParams, "ppeupperb",     &start_B,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppeupperb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
    tmpMin = -3.0;
    tmpMax = 2.0/3.0;
    start_a = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowera");
    if (ppt) {
      start_a = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowera");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE a fixed and set to %f\n",start_a);
		}else{
			LALInferenceAddVariable(currentParams, "ppelowera",     &start_a,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowera",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

    tmpMin = -4.5;
    tmpMax = 1.0;
    start_b = tmpMin+gsl_rng_uniform(GSLrandom)*(tmpMax-tmpMin);
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppelowerb");
    if (ppt) {
      start_b = atof(ppt->value);
    }    
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixppelowerb");
		if(ppt){
			LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(MPIrank==0) fprintf(stdout,"ppE b fixed and set to %f\n",start_b);
		}else{
			LALInferenceAddVariable(currentParams, "ppelowerb",     &start_b,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		}
    LALInferenceAddMinMaxPrior(priorArgs, "ppelowerb",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
  }

  /* If the currentParams are not in the prior, overwrite and pick paramaters from the priors. OVERWRITE EVEN USER CHOICES. 
  (necessary for complicated prior shapes where LALInferenceCyclicReflectiveBound() is not enought */
  while(state->prior(state, currentParams)<=-DBL_MAX){
    fprintf(stderr, "Warning initial parameter randlomy drawn from prior. (in %s, line %d)\n",__FILE__, __LINE__);
    LALInferenceVariables *temp; //
    temp=XLALCalloc(1,sizeof(LALInferenceVariables));
    memset(temp,0,sizeof(LALInferenceVariables));
    LALInferenceDrawApproxPrior(state, temp);
    LALInferenceCopyVariables(temp, currentParams);
  }
        /* Make sure that our initial value is within the
           prior-supported volume. */
        LALInferenceCyclicReflectiveBound(currentParams, priorArgs);

        /* Init covariance matrix, if specified.  The given file
           should contain the desired covariance matrix for the jump
           proposal, in row-major (i.e. C) order. */
        ppt=LALInferenceGetProcParamVal(commandLine, "--covarianceMatrix");
        if (ppt) {
          FILE *inp = fopen(ppt->value, "r");
          UINT4 N = LALInferenceGetVariableDimensionNonFixed(currentParams);
          gsl_matrix *covM = gsl_matrix_alloc(N,N);
          gsl_matrix *covCopy = gsl_matrix_alloc(N,N);
          REAL8Vector *sigmaVec = XLALCreateREAL8Vector(N);


          if (readSquareMatrix(covM, N, inp)) {
            fprintf(stderr, "Error reading covariance matrix (in %s, line %d)\n",
                    __FILE__, __LINE__);
            exit(1);
          }

          gsl_matrix_memcpy(covCopy, covM);

          for (i = 0; i < N; i++) {
            sigmaVec->data[i] = sqrt(gsl_matrix_get(covM, i, i)); /* Single-parameter sigma. */
          }

          LALInferenceAddVariable(state->proposalArgs, LALInferenceSigmaJumpName, &sigmaVec, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

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

          LALInferenceAddVariable(state->proposalArgs, "covarianceEigenvectors", &eVectors, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
          LALInferenceAddVariable(state->proposalArgs, "covarianceEigenvalues", &eigenValues, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

          fprintf(stdout, "Jumping with correlated jumps in %d dimensions from file %s.\n",
                  N, ppt->value);

          fclose(inp);
          gsl_eigen_symmv_free(ws);
          gsl_matrix_free(covCopy);
          gsl_vector_free(eValues);
        }

        /* Differential Evolution? */
        ppt=LALInferenceGetProcParamVal(commandLine, "--differential-evolution");
        if (ppt) {
          fprintf(stderr, "Using differential evolution.\nEvery Nskip parameters will be stored for use in the d.e. jump proposal.\n");

          state->differentialPoints = XLALCalloc(1, sizeof(LALInferenceVariables *));
          state->differentialPointsLength = 0;
          state->differentialPointsSize = 1;
        } else {
          state->differentialPoints = NULL;
          state->differentialPointsLength = 0;
          state->differentialPointsSize = 0;
        }

        UINT4 N = LALInferenceGetVariableDimensionNonFixed(currentParams);
  
        ppt=LALInferenceGetProcParamVal(commandLine, "--adapt");
        if (ppt) {
          fprintf(stdout, "Adapting single-param step sizes.\n");
          if (!LALInferenceCheckVariable(state->proposalArgs, SIGMAVECTORNAME)) {
            /* We need a sigma vector for adaptable jumps. */
            REAL8Vector *sigmas = XLALCreateREAL8Vector(N);
            
            for (i = 0; i < N; i++) {
              sigmas->data[i] = 1e-4;
            }
            

            LALInferenceAddVariable(state->proposalArgs, SIGMAVECTORNAME, &sigmas, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

          }
        }
        ppt=LALInferenceGetProcParamVal(commandLine, "--acceptanceRatio");
        if (ppt) {
  
          REAL8Vector *PacceptCount = XLALCreateREAL8Vector(N);
          REAL8Vector *PproposeCount = XLALCreateREAL8Vector(N);

          for (i = 0; i < N; i++) {
            PacceptCount->data[i] = 0.0;
            PproposeCount->data[i] = 0.0;
          }
          
          LALInferenceAddVariable(state->proposalArgs, "PacceptCount", &PacceptCount, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
          LALInferenceAddVariable(state->proposalArgs, "PproposeCount", &PproposeCount, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
        }

        INT4 adaptableStep = 0;
        LALInferenceAddVariable(state->proposalArgs, "adaptableStep", &adaptableStep, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

        INT4 varNumber = 0;
        LALInferenceAddVariable(state->proposalArgs, "proposedVariableNumber", &varNumber, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

        INT4 sigmasNumber = 0;
        LALInferenceAddVariable(state->proposalArgs, "proposedArrayNumber", &sigmasNumber, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

        REAL8 tau = 1e3;
        LALInferenceAddVariable(state->proposalArgs, "adaptTau", &tau, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

        ppt = LALInferenceGetProcParamVal(commandLine, "--adaptTau");
        if (ppt) {
          tau = atof(ppt->value);
          fprintf(stdout, "Setting adapt tau = %g.\n", tau);
          LALInferenceSetVariable(state->proposalArgs, "adaptTau", &tau);
        }
	
	return;
}




int main(int argc, char *argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

  if (MPIrank == 0) fprintf(stdout," ========== LALInference_MCMC ==========\n");

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
	procParams=LALInferenceParseCommandLine(argc,argv);
  
  if (LALInferenceGetProcParamVal(procParams, "--analyticLogLike")) {
    
    runState = calloc(1, sizeof(LALInferenceRunState));
    runState->commandLine=procParams;
    runState->data = NULL;
    initializeMCMC(runState);
    runState->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
    
    REAL8 Min=-1.0, Max=1.0;

    ppt=LALInferenceGetProcParamVal(procParams,"--xmin");
    if(ppt) Min=atof(ppt->value);
    ppt=LALInferenceGetProcParamVal(procParams,"--xmax");
    if(ppt)	Max=atof(ppt->value);
    
    REAL8 start_x1 = Min + gsl_rng_uniform(runState->GSLrandom)*(Max - (Min));
    ppt=LALInferenceGetProcParamVal(procParams,"--x1");
    if(ppt){start_x1=atof(ppt->value);}
    ppt=LALInferenceGetProcParamVal(procParams,"--fixX1");
    if(ppt){
      LALInferenceAddVariable(runState->currentParams, "x1",    &start_x1,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"x1 fixed and set to %f\n",start_x1);
    }else{
	    LALInferenceAddVariable(runState->currentParams, "x1",    &start_x1,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(runState->priorArgs,	"x1",	&Min,	&Max,		LALINFERENCE_REAL8_t);

    REAL8 start_x2 = Min + gsl_rng_uniform(runState->GSLrandom)*(Max - (Min));
    ppt=LALInferenceGetProcParamVal(procParams,"--x2");
    if(ppt){start_x2=atof(ppt->value);}
    ppt=LALInferenceGetProcParamVal(procParams,"--fixX2");
    if(ppt){
      LALInferenceAddVariable(runState->currentParams, "x2",    &start_x2,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"x2 fixed and set to %f\n",start_x2);
    }else{
	    LALInferenceAddVariable(runState->currentParams, "x2",    &start_x2,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(runState->priorArgs,	"x2",	&Min,	&Max,		LALINFERENCE_REAL8_t);
    
    UINT4 N = 2;
    //UINT4 N = 1;
    
    ppt=LALInferenceGetProcParamVal(procParams, "--adapt");
    if (ppt) {
      fprintf(stdout, "Adapting single-param step sizes.\n");
      if (!LALInferenceCheckVariable(runState->proposalArgs, SIGMAVECTORNAME)) {
        /* We need a sigma vector for adaptable jumps. */
        REAL8Vector *sigmas = XLALCreateREAL8Vector(N);
        UINT4 i = 0;
        
        for (i = 0; i < N; i++) {
          sigmas->data[i] = 1e-4;
        }
        
        
        LALInferenceAddVariable(runState->proposalArgs, SIGMAVECTORNAME, &sigmas, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
        
      }
    }
    ppt=LALInferenceGetProcParamVal(procParams, "--acceptanceRatio");
    if (ppt) {

      REAL8Vector *PacceptCount = XLALCreateREAL8Vector(N);
      REAL8Vector *PproposeCount = XLALCreateREAL8Vector(N);
      UINT4 i;
      
      for (i = 0; i < N; i++) {
        PacceptCount->data[i] = 0.0;
        PproposeCount->data[i] = 0.0;
      }
      
      LALInferenceAddVariable(runState->proposalArgs, "PacceptCount", &PacceptCount, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
      LALInferenceAddVariable(runState->proposalArgs, "PproposeCount", &PproposeCount, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);
    }
    
    INT4 adaptableStep = 0;
    LALInferenceAddVariable(runState->proposalArgs, "adaptableStep", &adaptableStep, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
    
    INT4 varNumber = 0;
    LALInferenceAddVariable(runState->proposalArgs, "proposedVariableNumber", &varNumber, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
    
    INT4 sigmasNumber = 0;
    LALInferenceAddVariable(runState->proposalArgs, "proposedArrayNumber", &sigmasNumber, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
    
    
  }else{
	
	ppt=LALInferenceGetProcParamVal(procParams,"--continue-run");
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

		procParams=LALInferenceParseCommandLine(fileargc,fileargv);
		//ppt = LALInferenceGetProcParamVal(procParams, "--randomseed");
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
	}//NOT analyticLogLike
  if(runState==NULL) {
    fprintf(stderr, "runState not allocated (%s, line %d).\n",
            __FILE__, __LINE__);
    exit(1);
  }
	printf(" ==== This is thread %d of %d ====\n", MPIrank, MPIsize);
	MPI_Barrier(MPI_COMM_WORLD);
	/* Call MCMC algorithm */
	runState->algorithm(runState);
	
  if (MPIrank == 0) printf(" ========== main(): finished. ==========\n");
  MPI_Finalize();
  return 0;
}




