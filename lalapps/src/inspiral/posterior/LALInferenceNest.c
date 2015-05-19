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
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceCalibrationErrors.h>

LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeNS(LALInferenceRunState *runState);
void initializeMalmquistPrior(LALInferenceRunState *runState);

void LogNSSampleAsMCMCSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars);                             
void LogNSSampleAsMCMCSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars);                              

void LogNSSampleAsMCMCSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  LALInferenceLogSampleToArray(state, vars);
  return;
}

void LogNSSampleAsMCMCSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars)
{
  LALInferenceLogSampleToFile(state, vars);
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
  (--maxmcmc M)\tUse at most this number of MCMC points when autodetermining the chain (5000).\n\
  (--Nmcmcinitial M)\tUse this number of MCMC points when initially resampling from the prior (otherwise default is to use maxmcmc).\n\
  (--sloppyratio S)\tNumber of sub-samples of the prior for every sample from the limited prior\n\
  (--Nruns R)\tNumber of parallel samples from logt to use(1)\n\
  (--tolerance dZ)\tTolerance of nested sampling algorithm (0.1)\n\
  (--randomseed seed)\tRandom seed of sampling distribution\n\
  (--prior )\t Set the prior to use (InspiralNormalised,SkyLoc,malmquist) default: InspiralNormalised\n\n\
  (--sampleprior N)\t For Testing: Draw N samples from the prior, will not perform the nested sampling integral\n\
  (--progress)\tOutput some progress information at each iteration\n\
  (--verbose)\tOutput more info. N=1: errors, N=2 (default): warnings, N=3: info \n\
  (--resume)\tAllow non-condor checkpointing every 4 hours. If give will check for OUTFILE_resume and continue if possible\n\
  ---------------------------------------------------------------------------------------------------\n\
  --- Noise Model -----------------------------------------------------------------------------------\n\
  ---------------------------------------------------------------------------------------------------\n\
  (--glitchFit)                    Run with glitch fitting\n\
  (--glitchNmax)                   Max number of glitches to fit\n\
  (--psdFit)                       Run with PSD fitting\n\
  (--psdNblock)                    Number of noise parameters per IFO channel (8)\n\
  (--psdFlatPrior)                 Use flat prior on psd parameters (Gaussian)\n\
  (--removeLines)                  Do include persistent PSD lines in fourier-domain integration\n\
  (--KSlines)                      Run with the KS test line removal\n\
  (--KSlinesWidth)                 Width of the lines removed by the KS test (deltaF)\n\
  (--chisquaredlines)              Run with the Chi squared test line removal\n\
  (--chisquaredlinesWidth)         Width of the lines removed by the Chi squared test (deltaF)\n\
  (--powerlawlines)                Run with the power law line removal\n\
  (--powerlawlinesWidth)           Width of the lines removed by the power law test (deltaF)\n\
  (--xcorrbands)                   Run PSD fitting with correlated frequency bands\n\
  \n";
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
  INT4 verbose=0;
  INT4 x=0;
  ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
  if(ppt) {
    if(ppt->value[0]){
      x=atoi(ppt->value);
      switch(x){
        case 0:
          verbose=LALNDEBUG; /* Nothing */
          break;
        case 1:
          verbose=LALMSGLVL1; /* Only errors */
          break;
        case 2:
          verbose=LALMSGLVL2; /* Errors and warnings */
          break;
        case 3:
          verbose=LALMSGLVL3; /* Errors, warnings and info */
          break;
        default:
          verbose=LALMSGLVL2;
          break;
      }
    }
    else verbose=LALMSGLVL2; /* Errors and warnings */
    LALInferenceAddVariable(irs->algorithmParams,"verbose", &verbose , LALINFERENCE_INT4_t,
                            LALINFERENCE_PARAM_FIXED);
  }
  
  INT4 tmpi=0,randomseed=0;
  REAL8 tmp=0;
  
  /* Single thread only */
  runState->threads=LALInferenceInitThreads(1);
  LALInferenceThreadState *threadState = runState->threads[0];
  
  /* Set up the appropriate functions for the nested sampling algorithm */
  runState->algorithm=&LALInferenceNestedSamplingAlgorithm;
  runState->evolve=&LALInferenceNestedSamplingOneStep;
  
  /* use the ptmcmc proposal to sample prior */
  threadState->proposal=&NSWrapMCMCLALProposal;
  REAL8 temp=1.0;
  LALInferenceAddVariable(threadState->proposalArgs,"temperature",&temp,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  
  /* Default likelihood is the frequency domain one */
  runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
  
  /* Check whether to use the SkyLocalization prior. Otherwise uses the default LALInferenceInspiralPrior. That should probably be replaced with a swhich over the possible priors. */
  ppt=LALInferenceGetProcParamVal(commandLine,"--prior");
  if(ppt){
    if(!strcmp(ppt->value,"SkyLoc")) runState->prior = &LALInferenceInspiralSkyLocPrior;
    if(!strcmp(ppt->value,"malmquist")) initializeMalmquistPrior(runState);
  }
  else{
    runState->prior = &LALInferenceInspiralPrior;
  }
  
  /* Set up the prior for analytic tests if needed */
  if(LALInferenceGetProcParamVal(commandLine,"--correlatedGaussianLikelihood")){
    runState->prior=LALInferenceAnalyticNullPrior;
  }
  if(LALInferenceGetProcParamVal(commandLine,"--bimodalGaussianLikelihood")){
    runState->prior=LALInferenceAnalyticNullPrior;
  }
  if(LALInferenceGetProcParamVal(commandLine,"--rosenbrockLikelihood")){
    runState->prior=LALInferenceAnalyticNullPrior;
  }
  
#ifdef HAVE_LIBLALXML
  threadState->logsample=LogNSSampleAsMCMCSampleToArray;
#else
  threadState->logsample=LogNSSampleAsMCMCSampleToFile;
#endif
  
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
  
  /* Maximum number of points in MCMC chain */
  ppt=LALInferenceGetProcParamVal(commandLine,"--maxmcmc");
  if(ppt){
    tmpi=atoi(ppt->value);
    LALInferenceAddVariable(runState->algorithmParams,"maxmcmc",&tmpi,
                            LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  }
  
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
  
  /* Get initial noise model likelihood */
  if (runState->data != NULL) {
    irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
    printf("Null Log Likelihood: %g\n", irs->currentLikelihood);
  }

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
  if(LALInferenceGetProcParamVal(procParams,"--help"))
  {
    fprintf(stdout,"%s",help);
  }
  
  /* initialise runstate based on command line */
  /* This includes reading in the data */
  /* And performing any injections specified */
  /* And allocating memory */
  state = LALInferenceInitCBCRunState(procParams);
  
  /* Set up structures for nested sampling */
  initializeNS(state);
  
  /* Set the CBC model */
  state->model = LALInferenceInitCBCModel(state);
  state->templt = LALInferenceInitCBCTemplate(state);
  /* Choose the likelihood */
  LALInferenceInitLikelihood(state);
  
  /* Set up the threads */
  LALInferenceInitCBCThreads(state,1);
  LALInferenceVariables *propArgs = LALInferenceParseProposalArgs(state);
  for(UINT4 i=0;i<state->nthreads;i++)
  {
    state->threads[i]->cycle=LALInferenceSetupDefaultInspiralProposalCycle(propArgs);
    LALInferenceRandomizeProposalCycle(state->threads[i]->cycle,state->GSLrandom);
  }
  
  /* Exit since we printed all command line arguments */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
  {
    exit(0);
  }
  
  /* Call setupLivePointsArray() to populate live points structures */
  LALInferenceSetupLivePointsArray(state);
  
  /* write injection with noise evidence information from algorithm */
  LALInferencePrintInjectionSample(state);
  
  /* Call nested sampling algorithm */
  state->algorithm(state);
  
  /* end */
  return(0);
}


