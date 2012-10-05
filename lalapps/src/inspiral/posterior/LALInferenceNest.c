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
void initVariablesReviewEvidence(LALInferenceRunState *state);
void initVariablesReviewEvidence_bimod(LALInferenceRunState *state);
void initVariablesReviewEvidence_banana(LALInferenceRunState *state);
void initStudentt(LALInferenceRunState *state);
void initializeTemplate(LALInferenceRunState *runState);
// static void mc2masses(double mc, double eta, double *m1, double *m2);
void LogNSSampleAsMCMCSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars);                             
void LogNSSampleAsMCMCSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars);                              
 


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


void LogNSSampleAsMCMCSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  NSFillMCMCVariables(vars,state->priorArgs);
  LALInferenceLogSampleToArray(state, vars);
  return;
}

void LogNSSampleAsMCMCSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  NSFillMCMCVariables(vars,state->priorArgs);
  LALInferenceLogSampleToFile(state, vars);
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
		printf("Null Log Likelihood: %g\n", irs->currentLikelihood);
	}
	else
	{
		fprintf(stdout, " initialize(): no data read.\n");
		exit(1);
	}
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
(--template [LAL,PhenSpin,LALGenerateInspiral,LALSim]\tSpecify template (default LAL)\n";
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
			fprintf(stderr,"ERROR: --template LALSTPN is deprecated. Try LALGenerateInspiral instead...\n");
			exit(1);
		}
		else if(!strcmp("PhenSpin",ppt->value))
			runState->template=&LALInferenceTemplatePSTRD;
		else if(!strcmp("LALGenerateInspiral",ppt->value))
			runState->template=&LALInferenceTemplateLALGenerateInspiral;
		else if(!strcmp("SpinTaylor",ppt->value))
			runState->template=&LALInferenceTemplateLALGenerateInspiral;
		else if(!strcmp("LAL",ppt->value))
			runState->template=&LALInferenceTemplateLAL;
        else if(!strcmp("LALSim",ppt->value))
            runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
		else {
			XLALPrintError("Error: unknown template %s\n",ppt->value);
			XLALPrintError(help);
			XLAL_ERROR_VOID(XLAL_EINVAL);
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
(--Nmcmc M)\tOver-ride auto chain length determination and use this number of MCMC samples.\n\
(--sloppyratio S)\tNumber of sub-samples of the prior for every sample from the limited prior\n\
(--Nruns R)\tNumber of parallel samples from logt to use(1)\n\
(--tolerance dZ)\tTolerance of nested sampling algorithm (0.1)\n\
(--randomseed seed)\tRandom seed of sampling distribution\n\
(--verbose)\tProduce progress information\n\
(--iotaDistance FRAC)\tPTMCMC: Use iota-distance jump FRAC of the time\n\
(--covarianceMatrix)\tPTMCMC: Propose jumps from covariance matrix of current live points\n\
(--differential-evolution)\tPTMCMC:Use differential evolution jumps\n\
(--prior_distr )\t Set the prior to use (for the moment the only possible choice is SkyLoc which will use the sky localization project prior. All other values or skipping this option select LALInferenceInspiralPriorNormalised)\n\
(--correlatedgaussianlikelihood)\tUse analytic, correlated Gaussian for Likelihood.\n\
(--bimodalgaussianlikelihood)\tUse analytic, bimodal correlated Gaussian for Likelihood.\n\
(--rosenbrocklikelihood \tUse analytic, Rosenbrock banana for Likelihood.\n";
//(--tdlike)\tUse time domain likelihood.\n";

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
	
    /* use the ptmcmc proposal to sample prior */
    runState->proposal=&NSWrapMCMCLALProposal;
    REAL8 temp=1.0;
    LALInferenceAddVariable(runState->proposalArgs,"temperature",&temp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
	
	/* Default likelihood is the frequency domain one */
	runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;

        /* Check whether to use the SkyLocalization prior. Otherwise uses the default LALInferenceInspiralPriorNormalised. That should probably be replaced with a swhich over the possible priors. */
        ppt=LALInferenceGetProcParamVal(commandLine,"--prior_distr");
        if(ppt){
            if (!strcmp(ppt->value,"SkyLoc")) runState->prior = &LALInferenceInspiralSkyLocPrior;
        }
        else{
            runState->prior = &LALInferenceInspiralPriorNormalised;
        }
	

	if(LALInferenceGetProcParamVal(commandLine,"--correlatedgaussianlikelihood")){
        	runState->likelihood=&LALInferenceCorrelatedAnalyticLogLikelihood;
		runState->prior=LALInferenceAnalyticNullPrior;
	}
    	if(LALInferenceGetProcParamVal(commandLine,"--bimodalgaussianlikelihood")){
        	runState->likelihood=&LALInferenceBimodalCorrelatedAnalyticLogLikelihood;
		runState->prior=LALInferenceAnalyticNullPrior;
	}
        if(LALInferenceGetProcParamVal(commandLine,"--rosenbrocklikelihood")){
                runState->likelihood=&LALInferenceRosenbrockLogLikelihood;
                runState->prior=LALInferenceAnalyticNullPrior;
        }

//	if(LALInferenceGetProcParamVal(commandLine,"--tdlike")){
//		fprintf(stderr, "Computing likelihood in the time domain.\n");
//		runState->likelihood=&LALInferenceTimeDomainLogLikelihood;
//    	}
    
	#ifdef HAVE_LIBLALXML
	runState->logsample=LogNSSampleAsMCMCSampleToArray;
	#else
	runState->logsample=LogNSSampleAsMCMCSampleToFile;
	#endif
	
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
    if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--nlive");
	if(ppt)
		tmpi=atoi(ppt->value);
	else {
		fprintf(stderr,"Error, must specify number of live points\n");
		exit(1);
	}
	LALInferenceAddVariable(runState->algorithmParams,"Nlive",&tmpi, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
	
	/* Number of points in MCMC chain */
	ppt=LALInferenceGetProcParamVal(commandLine,"--Nmcmc");
    	if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--nmcmc");
	if(ppt){
	  tmpi=atoi(ppt->value);
	LALInferenceAddVariable(runState->algorithmParams,"Nmcmc",&tmpi,
				LALINFERENCE_INT4_t,LALINFERENCE_PARAM_OUTPUT);
	printf("set number of MCMC points, over-riding auto-determination!\n");
	}
	if((ppt=LALInferenceGetProcParamVal(commandLine,"--sloppyfraction")))
        	tmp=atof(ppt->value);
    	else tmp=0.0;
    	LALInferenceAddVariable(runState->algorithmParams,"sloppyfraction",&tmp,
                    LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

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

/* Setup the variable for the evidence calculation test for review */
/* 5-sigma ranges for analytic likeliood function */
/* https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/LALInferenceReviewAnalyticGaussianLikelihood */
void initVariablesReviewEvidence(LALInferenceRunState *state)
{
    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;
    char **strings=NULL;
    char *pinned_params=NULL;
    UINT4 N=0,i,j;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    }
	LALInferenceVariables *priorArgs=state->priorArgs;
        state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
        LALInferenceVariables *currentParams=state->currentParams;
	i=0;

	struct varSettings {const char *name; REAL8 val, min, max;};
	
	struct varSettings setup[]=
	{
		{.name="time", .val=0.0, .min=-0.1073625, .max=0.1073625},
		{.name="m1", .val=16., .min=14.927715, .max=17.072285},
		{.name="m2", .val=7., .min=5.829675, .max=8.170325},
		{.name="distance", .val=50., .min=37.986000000000004, .max=62.013999999999996},
		{.name="inclination", .val=LAL_PI/2., .min=1.4054428267948966, .max=1.7361498267948965},
		{.name="phase", .val=LAL_PI, .min=2.8701521535897934, .max=3.413033153589793},
		{.name="polarisation", .val=LAL_PI/2., .min=1.3885563267948966, .max=1.7530363267948965},
		{.name="rightascension", .val=LAL_PI, .min=2.813050153589793, .max=3.4701351535897933},
		{.name="declination", .val=0., .min=-0.300699, .max=0.300699},
		{.name="a_spin1", .val=0.5, .min=0.3784565, .max=0.6215435},
		{.name="a_spin2", .val=0.5, .min=0.421869, .max=0.578131},
		{.name="theta_spin1", .val=LAL_PI/2., .min=1.3993998267948966, .max=1.7421928267948965},
		{.name="theta_spin2", .val=LAL_PI/2., .min=1.4086158267948965, .max=1.7329768267948966},
		{.name="phi_spin1", .val=LAL_PI, .min=2.781852653589793, .max=3.501332653589793},
		{.name="phi_spin2", .val=LAL_PI, .min=2.777215653589793, .max=3.5059696535897933},
		{.name="END", .val=0., .min=0., .max=0.}
	};

	while(strcmp("END",setup[i].name))
	{
        LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
        /* Check if it is to be fixed */
        for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
		LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
	    LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
		i++;
	}
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
  LALInferenceIFOData *dataPtr;
  LALInferenceDomain modelDomain;
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime;
	ProcessParamsTable *ppt=NULL;
	LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
	INT4 AmpOrder=0;
	Approximant approx=TaylorF2;
	REAL8 logDmin=log(1.0);
	REAL8 logDmax=log(100.0);
	REAL8 mcMin=1.0;
	REAL8 mcMax=20.5;
	REAL8 logmcMax,logmcMin,mMin=1.0,mMax=30.0;
	REAL8 a_spin2_max=1.0, a_spin1_max=1.0;
	REAL8 a_spin2_min=0.0, a_spin1_min=0.0;
	REAL8 phi_spin1_min=0.0;
	REAL8 phi_spin1_max=2.0*LAL_PI;
	REAL8 theta_spin1_min=0.0;
	REAL8 theta_spin1_max=LAL_PI;
	REAL8 fRef=0.; /* freq. at which precessing "initial" cond. specified */
	REAL8 etaMin=0.01;
	REAL8 etaMax=0.25;
	REAL8 dt=0.1;            /* Width of time prior */
	REAL8 tmpMin,tmpMax,tmpVal;
	REAL8 m1_min=0.;	
	REAL8 m1_max=0.;
	REAL8 m2_min=0.;
	REAL8 m2_max=0.;
    REAL8 mtot_min=0.0;
    REAL8 mtot_max=0.0;
	memset(currentParams,0,sizeof(LALInferenceVariables));
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
	INT4 enable_spin=0;
	INT4 aligned_spin=0;
	char *pinned_params=NULL;
	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tInjection XML file to use\n\
(--Mmin mchirp)\tMinimum chirp mass\n\
(--Mmax mchirp)\tMaximum chirp mass\n\
(--etamin eta)\tMinimum eta\n\
(--etamax eta)\tMaximum eta\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--Dmin dist)\tMinimum distance in Mpc (1)\n\
(--Dmax dist)\tMaximum distance in Mpc (100)\n\
(--approx ApproximantorderPN)\tSpecify a waveform to use, (default TaylorF2threePointFivePN)\n\
(--amporder INT)\tSpecify post-Newtonian amplitude order to use (defaults to 0. -1 will use highest available)\n\
(--compmin min)\tMinimum component mass (1.0)\n\
(--compmax max)\tMaximum component mass (30.0)\n\
(--mtotalmin)\tMinimum total mass (2*compmin)\n\
(--mtotalmax)\tMaximum total mass (2*compmax)\n\
(--enable-spin)\tEnable spin parameters\n\
(--aligned-spin)\tUse only aligned spin parameters (uses spins between -1 and 1)\n\
(--approx ApproximantphaseOrderPN)\tSet approximant (PhenSpin implicitly enables spin)\n\
(--s1max SPIN)\tMax magnitude of spin (on both bodies!)\n\
(--s1min SPIN)\tMin magnitude of spin (on both bodies!)\n\
(--fref fRef)\tSpecify a reference frequency at which parameters are defined (default 0).\n\
(--mcq)\tUse chirp mass and asymmetric mass ratio (m1/m2) as variables\n\
(--crazyinjectionhlsign)\tFlip the sign of HL signal in likelihood function\n\
(--pinparams [mchirp,asym_massratio,etc])\n\tList of parameters to set to injected values\n\
(--no-logdistance)\tUse distance, not logdistance, as the sampling variable\n";

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
		//Select event
		ppt=LALInferenceGetProcParamVal(commandLine,"--event");
		if(ppt){
		  event = atoi(ppt->value);
		  while(i<event) {i++; injTable = injTable->next;}
		}
		endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
        fprintf(stderr,"Read trig time %lf from injection XML file\n",endtime);
		AmpOrder=injTable->amp_order;
		PhaseOrder = XLALGetOrderFromString(injTable->waveform);
		if( (int) PhaseOrder == XLAL_FAILURE)
		  ABORTXLAL(&status);
		approx = XLALGetApproximantFromString(injTable->waveform);
		if( (int) approx == XLAL_FAILURE)
		  ABORTXLAL(&status);
		/* See if there are any parameters pinned to injection values */
		if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
			pinned_params=ppt->value;
			LALInferenceVariables tempParams;
			memset(&tempParams,0,sizeof(tempParams));
			char **strings=NULL;
			UINT4 N;
			LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
			LALInferenceInjectionToVariables(injTable,&tempParams);
			LALInferenceVariableItem *node=NULL;
			while(N>0){
				N--;
				char *name=strings[N];
				node=LALInferenceGetItem(&tempParams,name);
				if(node) LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary);
				else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",name);}
			}
		}
	}

	/* Over-ride approximant if user specifies */
	ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
	if(!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
	if(ppt){
		approx = XLALGetApproximantFromString(ppt->value);
		if( (int) approx == XLAL_FAILURE)
			ABORTXLAL(&status);
        	PhaseOrder = XLALGetOrderFromString(ppt->value);
	        if( (int) PhaseOrder == XLAL_FAILURE)
	                ABORTXLAL(&status);
	}
	ppt=LALInferenceGetProcParamVal(commandLine,"--amporder");
	if(ppt) AmpOrder=atoi(ppt->value);
	fprintf(stdout,"Templates will run using Approximant %i, phase order %i, amp order %i\n",approx,PhaseOrder,AmpOrder);

	/* Set the modeldomain appropriately */
	switch(approx)
	{
		case GeneratePPN:
		case TaylorT1:
		case TaylorT2:
		case TaylorT3:
		case TaylorT4:
		case EOB:
		case EOBNR:
		case EOBNRv2:
		case EOBNRv2HM:
		case SpinTaylor:
		case SpinTaylorT4:
		case SpinQuadTaylor:
		case SpinTaylorFrameless:
		case PhenSpinTaylorRD:
		case NumRel:
			modelDomain=LALINFERENCE_DOMAIN_TIME;
			break;
		case TaylorF1:
		case TaylorF2:
		case TaylorF2RedSpin:
		case TaylorF2RedSpinTidal:
		case IMRPhenomA:
		case IMRPhenomB:
			modelDomain=LALINFERENCE_DOMAIN_FREQUENCY;
			break;
		default:
			fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
			exit(1);
			break;
	}

  /* Set model domain for all IFOs */
  dataPtr = state->data;
  while (dataPtr != NULL) {
    dataPtr->modelDomain = modelDomain;
    dataPtr = dataPtr->next;
  }

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
	
    /* Set the minimum and maximum total mass, using user values if specified */
    ppt=LALInferenceGetProcParamVal(commandLine,"--mtotalmin");
    if(ppt) mtot_min=atof(ppt->value);
    else mtot_min=2.*mMin;
    LALInferenceAddVariable(priorArgs,"MTotMin",&mtot_min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

    ppt=LALInferenceGetProcParamVal(commandLine,"--mtotalmax");
    if(ppt) mtot_max=atof(ppt->value);
    else mtot_max=2.*(mMax-mMin);
    LALInferenceAddVariable(priorArgs,"MTotMax",&mtot_max,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

    /* Set the minimum and maximum chirp mass, using user values if specified */
    ppt=LALInferenceGetProcParamVal(commandLine,"--Mmin");
    if(ppt)
        mcMin=atof(ppt->value);
    else mcMin=pow(mMin*mMin,0.6)/pow(2.0*mMin,0.2);
    ppt=LALInferenceGetProcParamVal(commandLine,"--Mmax");
    if(ppt)
        mcMax=atof(ppt->value);
    else mcMax=pow(mMax*mMax,0.6)/pow(2.0*mMax,0.2);

    INT4 tempint=1;
	if(LALInferenceGetProcParamVal(commandLine,"--crazyinjectionhlsign") || LALInferenceGetProcParamVal(commandLine,"--crazyInjectionHLSign"))
    {
        printf("Using signal sign flip in Hanford and Livingston");
        LALInferenceAddVariable(currentParams,"crazyInjectionHLSign",&tempint,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }
	printf("Read end time %f\n",endtime);
	
	LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
    	LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(currentParams, "LAL_AMPORDER", &AmpOrder, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
	
    ppt=LALInferenceGetProcParamVal(commandLine,"--mcq");
    if(ppt) /* Use MC and Q as sampling variables */
    {
        /* Set up the variable parameters */
        tmpVal=mcMin+(mcMax-mcMin)/2.0;
        
        LALInferenceAddMinMaxPrior(priorArgs,   "chirpmass",    &mcMin, &mcMax,     LALINFERENCE_REAL8_t);
        if(!LALInferenceCheckVariable(currentParams,"chirpmass")) LALInferenceAddVariable(currentParams,"chirpmass",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        tmpVal=1.5;
        REAL8 qMax=1.0;
        REAL8 qMin=mMin/mMax;
        if(!LALInferenceCheckVariable(currentParams,"asym_massratio")) LALInferenceAddVariable(currentParams, "asym_massratio",       &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddMinMaxPrior(priorArgs,   "asym_massratio",    &qMin,    &qMax,    LALINFERENCE_REAL8_t);

    }
    else /* Use log chirp mass and eta (default) */
    {
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
	}

    	if(!LALInferenceCheckVariable(currentParams,"time")) LALInferenceAddVariable(currentParams, "time",            &endtime   ,           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
	tmpMin=endtime-0.5*dt; tmpMax=endtime+0.5*dt;
	LALInferenceAddMinMaxPrior(priorArgs, "time",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);	

	tmpVal=1.0;
    	if(!LALInferenceCheckVariable(currentParams,"phase")) LALInferenceAddVariable(currentParams, "phase",           &tmpVal,             LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "phase",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
	if(LALInferenceGetProcParamVal(commandLine,"--no-logdistance"))
	{
		REAL8 Dmin=exp(logDmin);
		REAL8 Dmax=exp(logDmax);
		tmpVal=Dmin+(Dmax-Dmin)/2.;
		LALInferenceAddVariable(currentParams,"distance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);		
	}
	else 
	{
		tmpVal=logDmin+(logDmax-logDmin)/2.0;
		if(!LALInferenceCheckVariable(currentParams,"logdistance")) LALInferenceAddVariable(currentParams,"logdistance", &tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		LALInferenceAddMinMaxPrior(priorArgs, "logdistance",     &logDmin, &logDmax,   LALINFERENCE_REAL8_t);
	}
	tmpVal=1.0;
	if(!LALInferenceCheckVariable(currentParams,"rightascension")) LALInferenceAddVariable(currentParams, "rightascension",  &tmpVal,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

	if(!LALInferenceCheckVariable(currentParams,"declination")) LALInferenceAddVariable(currentParams, "declination",     &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=-LAL_PI/2.0; tmpMax=LAL_PI/2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "declination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
    
	if(!LALInferenceCheckVariable(currentParams,"polarisation")) LALInferenceAddVariable(currentParams, "polarisation",    &tmpVal,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);
	
 	if(!LALInferenceCheckVariable(currentParams,"inclination")) LALInferenceAddVariable(currentParams, "inclination",     &tmpVal,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	tmpMin=0.0; tmpMax=LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &tmpMin, &tmpMax,   LALINFERENCE_REAL8_t);

	/* 
	 * fRef used by SpinTaylorT4 to determine which frequency the reference
	 * phase and "initial" values of spin components refer to.
	 * fRef=0 is the standard behavior consistent with other approximants.
	 * it means the spin components are at the initial frequency and phiRef
	 * is the "phase at coalescence" (the last sample)
	 * fRef > 0 means the provided phiRef and spin components will be the
	 * values when the binary has GW frequency fRef.
	 */
	ppt=LALInferenceGetProcParamVal(commandLine,"--fref");
	if(ppt) fRef = atof(ppt->value);
	else fRef = 0.;
	LALInferenceAddVariable(currentParams, "fRef", &fRef, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

	/* Additional parameters for spinning waveforms */
	ppt=LALInferenceGetProcParamVal(commandLine,"--template");
	if(ppt) if(!strcmp("PhenSpin",ppt->value)){ enable_spin=1;}

	if(LALInferenceGetProcParamVal(commandLine,"--enable-spin")) enable_spin=1;
	
	/* If aligned spins use magnitude in (-1,1) */
	ppt=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
	if(ppt) {enable_spin=1; aligned_spin=1; a_spin1_min=-1; a_spin2_min=-1;}
	
	if(enable_spin){
		tmpVal=a_spin1_min+(a_spin1_max-a_spin1_min)/2.0;
		if(!LALInferenceCheckVariable(currentParams,"a_spin1")) LALInferenceAddVariable(currentParams, "a_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &a_spin1_min, &a_spin1_max,   LALINFERENCE_REAL8_t); 
	        
		tmpVal=a_spin2_min+(a_spin2_max-a_spin2_min)/2.0;
		if(!LALInferenceCheckVariable(currentParams,"a_spin2")) LALInferenceAddVariable(currentParams, "a_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR); 
		LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &a_spin2_min, &a_spin2_max,   LALINFERENCE_REAL8_t); 
	
		
		if(aligned_spin){ /* Set the spin angles to be parallel to orbital */
			tmpVal=0;
			if(!LALInferenceCheckVariable(currentParams,"theta_spin1")) LALInferenceAddVariable(currentParams,"theta_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(!LALInferenceCheckVariable(currentParams,"theta_spin2")) LALInferenceAddVariable(currentParams,"theta_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(!LALInferenceCheckVariable(currentParams,"phi_spin1")) LALInferenceAddVariable(currentParams,"phi_spin1",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
			if(!LALInferenceCheckVariable(currentParams,"phi_spin2")) LALInferenceAddVariable(currentParams,"phi_spin2",&tmpVal, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
		}
		else{ /* Use full spinning parameters */
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)/2.0;

			if(!LALInferenceCheckVariable(currentParams,"theta_spin1")) LALInferenceAddVariable(currentParams,"theta_spin1",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=theta_spin1_min+(theta_spin1_max - theta_spin1_min)/2.0;
			if(!LALInferenceCheckVariable(currentParams,"theta_spin2")) LALInferenceAddVariable(currentParams,"theta_spin2",	&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
			LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &theta_spin1_min, &theta_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)/2.0;
	
			if(!LALInferenceCheckVariable(currentParams,"phi_spin1")) LALInferenceAddVariable(currentParams,"phi_spin1",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
			LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &phi_spin1_min, &phi_spin1_max,   LALINFERENCE_REAL8_t); 
	
			tmpVal=phi_spin1_min+(phi_spin1_max - phi_spin1_min)/2.0;
			if(!LALInferenceCheckVariable(currentParams,"phi_spin2")) LALInferenceAddVariable(currentParams,"phi_spin2",		&tmpVal,	LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
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
	/* Review task needs special priors */
	if(LALInferenceGetProcParamVal(procParams,"--correlatedgaussianlikelihood"))
		initVariablesReviewEvidence(state);
        else if(LALInferenceGetProcParamVal(procParams,"--bimodalgaussianlikelihood"))
                initVariablesReviewEvidence_bimod(state);
        else if(LALInferenceGetProcParamVal(procParams,"--rosenbrocklikelihood"))
                initVariablesReviewEvidence_banana(state);
	else
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

	/* write injection with noise evidence information from algorithm */
    LALInferencePrintInjectionSample(state);

	/* end */
	return(0);
}

void initVariablesReviewEvidence_bimod(LALInferenceRunState *state)
{
    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;
    char **strings=NULL;
    char *pinned_params=NULL;
    UINT4 N=0,i,j;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    }
        LALInferenceVariables *priorArgs=state->priorArgs;
        state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
        LALInferenceVariables *currentParams=state->currentParams;
        i=0;

        struct varSettings {const char *name; REAL8 val, min, max;};

        struct varSettings setup[]=
        {
                {.name="time", .val=0.05589, .min=-0.1373625, .max=0.2491425},
                {.name="m1", .val=16.857828, .min=14.927715, .max=18.787941},
                {.name="m2", .val=7.93626, .min=5.829675, .max=10.042845},
                {.name="distance", .val=34.6112, .min=12.986, .max=56.2364},
                {.name="inclination", .val=0.9176809634, .min=0.6200446634, .max=1.2153172634},
                {.name="phase", .val=1.7879487268, .min=1.2993558268, .max=2.2765416268},
                {.name="polarisation", .val=0.9311901634, .min=0.6031581634, .max=1.2592221634},
                {.name="rightascension", .val=1.8336303268, .min=1.2422538268, .max=2.4250068268},
                {.name="declination", .val=-0.5448389634, .min=-1.0860971634, .max=-0.0035807634},
                {.name="a_spin1", .val=0.2972348, .min=0.0784565, .max=0.5160131},
                {.name="a_spin2", .val=0.2625048, .min=0.121869, .max=0.4031406},
                {.name="theta_spin1", .val=0.9225153634, .min=0.6140016634, .max=1.2310290634},
                {.name="theta_spin2", .val=0.9151425634, .min=0.6232176634, .max=1.2070674634},
                {.name="phi_spin1", .val=1.8585883268, .min=1.2110563268, .max=2.5061203268},
                {.name="phi_spin2", .val=1.8622979268, .min=1.2064193268, .max=2.5181765268},
                {.name="END", .val=0., .min=0., .max=0.}
        };

        while(strcmp("END",setup[i].name))
        {
        LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
        /* Check if it is to be fixed */
        for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
                LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
            LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
                i++;
        }
        return;
}

void initVariablesReviewEvidence_banana(LALInferenceRunState *state)
{
    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;
    char **strings=NULL;
    char *pinned_params=NULL;
    UINT4 N=0,i,j;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    }
        LALInferenceVariables *priorArgs=state->priorArgs;
        state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
        LALInferenceVariables *currentParams=state->currentParams;
        i=0;

        struct varSettings {const char *name; REAL8 val, min, max;};

        struct varSettings setup[]=
        {
                {.name="time", .val=0.0, .min=-2., .max=2.},
                {.name="m1", .val=16., .min=14., .max=18.},
                {.name="m2", .val=7., .min=5., .max=9.},
                {.name="distance", .val=50., .min=48., .max=52.},
                {.name="inclination", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
                {.name="phase", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
                {.name="polarisation", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
                {.name="rightascension", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
                {.name="declination", .val=0., .min=-2., .max=2.},
                {.name="a_spin1", .val=0.5, .min=-1.5, .max=2.5},
                {.name="a_spin2", .val=0.5, .min=-1.5, .max=2.5},
                {.name="theta_spin1", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
                {.name="theta_spin2", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
                {.name="phi_spin1", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
                {.name="phi_spin2", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
                {.name="END", .val=0., .min=0., .max=0.}
        };

        while(strcmp("END",setup[i].name))
        {
        LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
        /* Check if it is to be fixed */
        for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
                LALInferenceAddVariable(currentParams,setup[i].name, &(setup[i].val) ,LALINFERENCE_REAL8_t, type);
            LALInferenceAddMinMaxPrior(priorArgs, setup[i].name,    &(setup[i].min),    &(setup[i].max),    LALINFERENCE_REAL8_t);
                i++;
        }
        return;
}

