/*
 *  LALInferenceMCMC.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch, Will M. Farr, and Ben Farr
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
#include <lal/LALInferenceInit.h>
#include <lalapps.h>

#include <mpi.h>


int MPIrank, MPIsize;

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


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMCMC(LALInferenceRunState *runState);
REAL8 **parseMCMCoutput(char ***params, UINT4 *nInPar, UINT4 *nInSamps, char *infilename, UINT4 burnin);


/* This contains code chopped from LALInferenceInitCBC that wasn't
 * related to the template intialisation but to the guts of the MCMC
 * algorithm */
void LALInferenceInitMCMCState(LALInferenceRunState *state);
void LALInferenceInitMCMCState(LALInferenceRunState *state)
{
  
  if(state==NULL)
  {
    return;
  }
  LALInferenceVariables *currentParams=state->currentParams;
  LALInferenceVariables *priorArgs=state->priorArgs;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  UINT4 i=0;
  
  /* Initialize variable that will store the name of the last proposal function used */
  const char *initPropName = "INITNAME";
  LALInferenceAddVariable(state->proposalArgs, LALInferenceCurrentProposalName, &initPropName, LALINFERENCE_string_t, LALINFERENCE_PARAM_LINEAR);

  /* If using a malmquist prior, force a strict prior window on distance for starting point, otherwise
   * the approximate prior draws are very unlikely to be within the malmquist prior */
  REAL8 dist_low, dist_high;
  REAL8 restricted_dist_low = 10.0;
  REAL8 restricted_dist_high = 100.0;
  INT4 changed_dist = 0;
  if (LALInferenceCheckVariable(state->priorArgs, "malmquist") && LALInferenceCheckVariableNonFixed(currentParams, "distance")) {
      changed_dist = 1;
      LALInferenceGetMinMaxPrior(state->priorArgs, "distance", &dist_low, &dist_high);
      LALInferenceRemoveMinMaxPrior(state->priorArgs, "distance");
      LALInferenceAddMinMaxPrior(state->priorArgs, "distance", &restricted_dist_low, &restricted_dist_high, LALINFERENCE_REAL8_t);
  }

  /* If the currentParams are not in the prior, overwrite and pick paramaters from the priors. OVERWRITE EVEN USER CHOICES.
   *     (necessary for complicated prior shapes where LALInferenceCyclicReflectiveBound() is not enough */
  LALInferenceVariables *temp=XLALCalloc(1,sizeof(LALInferenceVariables));
  while(state->prior(state, currentParams, state->model)<=-DBL_MAX){
    fprintf(stderr, "Warning initial parameter randlomy drawn from prior. (in %s, line %d)\n",__FILE__, __LINE__);
    LALInferenceDrawApproxPrior(state, currentParams, temp);
    LALInferenceCopyVariables(temp, currentParams);
  }
  LALInferenceClearVariables(temp);
  XLALFree(temp);

  /* Make sure that our initial value is within the
   *     prior-supported volume. */
  LALInferenceCyclicReflectiveBound(currentParams, priorArgs);
  
  /* Replace distance prior if changed for initial sample draw */
  if (changed_dist) {
      LALInferenceRemoveMinMaxPrior(state->priorArgs, "distance");
      LALInferenceAddMinMaxPrior(state->priorArgs, "distance", &dist_low, &dist_high, LALINFERENCE_REAL8_t);
  }

  /* Init covariance matrix, if specified.  The given file
   *     should contain the desired covariance matrix for the jump
   *     proposal, in row-major (i.e. C) order. */
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
  ppt=LALInferenceGetProcParamVal(commandLine, "--noDifferentialEvolution");
  if (!ppt) {
    fprintf(stderr, "Using differential evolution.\nEvery Nskip parameters will be stored for use in the d.e. jump proposal.\n");
    
    state->differentialPoints = XLALCalloc(1, sizeof(LALInferenceVariables *));
    state->differentialPointsLength = 0;
    state->differentialPointsSize = 1;
  } else {
    fprintf(stderr, "Differential evolution disabled (--noDifferentialEvolution).\n");
    state->differentialPoints = NULL;
    state->differentialPointsLength = 0;
    state->differentialPointsSize = 0;
  }
  
  /* kD Tree NCell parameter. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--kDNCell");
  if (ppt) {
    INT4 NCell = atoi(ppt->value);
    LALInferenceAddVariable(state->proposalArgs, "KDNCell", &NCell, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }

  /* KD Tree propsal. */
  ppt=LALInferenceGetProcParamVal(commandLine, "--kDTree");
  if (!ppt) {
    ppt = LALInferenceGetProcParamVal(commandLine, "--kdtree");
  }
  if (ppt) {
    LALInferenceKDTree *tree;
    REAL8 *low, *high;
    currentParams = state->currentParams;
    LALInferenceVariables *template = XLALCalloc(1,sizeof(LALInferenceVariables));
    size_t ndim = LALInferenceGetVariableDimensionNonFixed(currentParams);
    LALInferenceVariableItem *currentItem;
    
    low = XLALMalloc(ndim*sizeof(REAL8));
    high = XLALMalloc(ndim*sizeof(REAL8));
    
    currentItem = currentParams->head;
    i = 0;
    while (currentItem != NULL) {
      if (currentItem->vary != LALINFERENCE_PARAM_FIXED) {
	LALInferenceGetMinMaxPrior(state->priorArgs, currentItem->name, &(low[i]), &(high[i]));
	i++;
      }
      currentItem = currentItem->next;
    }
    
    tree = LALInferenceKDEmpty(low, high, ndim);
    LALInferenceCopyVariables(currentParams, template);
    
    LALInferenceAddVariable(state->proposalArgs, "kDTree", &tree, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(state->proposalArgs, "kDTreeVariableTemplate", &template, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
  }
  
  INT4 Neff = 0;
  ppt = LALInferenceGetProcParamVal(commandLine, "--Neff");
  if (ppt)
    Neff = atoi(ppt->value);
  LALInferenceAddVariable(state->algorithmParams, "Neff", &Neff, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
  
  return;
}



LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
  LALInferenceRunState *irs=NULL;
  LALInferenceIFOData *ifoPtr;
  UINT4 nifo = 0;
  unsigned int n_basis, n_samples, time_steps;
  n_basis = 965;//TODO: have it read from file or from command line.
  
  ProcessParamsTable *ppt=NULL;
  FILE *tempfp;
  if(LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")){
    ppt=LALInferenceGetProcParamVal(commandLine,"--roqtime_steps");
    tempfp = fopen (ppt->value,"r");
    fscanf (tempfp, "%u", &time_steps);
    fscanf (tempfp, "%u", &n_basis);
    fscanf (tempfp, "%u", &n_samples);
  }  

  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  irs = XLALCalloc(1, sizeof(LALInferenceRunState));
  /* read data from files: */
  fprintf(stdout, " ==== LALInferenceReadData(): started. ====\n");
  irs->commandLine=commandLine;
  LALInferenceCheckOptionsConsistency(commandLine);
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
    while (ifoPtr) {
        nifo++;
        ifoPtr = ifoPtr->next;
    }

    irs->currentIFOLikelihoods = XLALMalloc(nifo * sizeof(REAL8));
    irs->currentIFOSNRs = XLALMalloc(nifo * sizeof(REAL8));

    irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
    printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
  }
  else{
    fprintf(stdout, " initialize(): no data read.\n");
    irs = NULL;
  }

  return(irs);
}

/********** Initialise MCMC structures *********/

/************************************************/
void initializeMCMC(LALInferenceRunState *runState)
{
  char help[]="\
               ---------------------------------------------------------------------------------------------------\n\
               --- General Algorithm Parameters ------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--Niter N)                      Number of iterations (2*10^7).\n\
               (--Neff N)                       Number of effective samples. (ends if chain surpasses Niter)\n\
               (--Nskip N)                      Number of iterations between disk save (100).\n\
               (--trigSNR SNR)                  Network SNR from trigger, used to calculate tempMax (injection SNR).\n\
               (--randomseed seed)              Random seed of sampling distribution (random).\n\
               (--adaptTau)                     Adaptation decay power, results in adapt length of 10^tau (5).\n\
               (--noAdapt)                      Do not adapt run.\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Likelihood Functions --------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--zeroLogLike)                  Use flat, null likelihood.\n\
               (--studentTLikelihood)           Use the Student-T Likelihood that marginalizes over noise.\n\
               (--correlatedGaussianLikelihood) Use analytic, correlated Gaussian for Likelihood.\n\
               (--bimodalGaussianLikelihood)    Use analytic, bimodal correlated Gaussian for Likelihood.\n\
               (--rosenbrockLikelihood)         Use analytic, Rosenbrock banana for Likelihood.\n\
               (--analyticnullprior)            Use analytic null prior.\n\
               (--nullprior)                    Use null prior in the sampled parameters.\n\
               (--noiseonly)                    Use signal-free log likelihood (noise model only).\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Noise Model -----------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
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
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Proposals  ------------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--rapidSkyLoc)                  Use rapid sky localization jump proposals.\n\
               (--kDTree)                       Use a kDTree proposal.\n\
               (--kDNCell N)                    Number of points per kD cell in proposal.\n\
               (--covarianceMatrix file)        Find the Cholesky decomposition of the covariance matrix for jumps in file.\n\
               (--noProposalSkyRing)              Disable the proposal that rotates sky position\n\
                                                  around vector connecting any two IFOs in network.\n\
               (--noProposalCorrPsiPhi)           Disable the proponal that jumps along psi-phi \n\
                                                  correlation\n\
               (--noDifferentialEvolution)      Disable the differential-evolution proposal\n\
               (--differential-buffer-limit)    Limit the number of stored differential-evolution points\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Parallel Tempering Algorithm Parameters -------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--inverseLadder)                Space temperature uniform in 1/T, rather than geometric.\n\
               (--tempLadderBottomUp)           Construct the a geometric temperature ladder with tempDelta=1+sqrt(2/nPar).\n\
               (--tempSkip N)                   Number of iterations between temperature swap proposals (100).\n\
               (--tempKill N)                   Iteration number to stop temperature swapping (Niter).\n\
               (--tempMin T)                    Lowest temperature for parallel tempering (1.0).\n\
               (--tempMax T)                    Highest temperature for parallel tempering (50.0).\n\
               (--anneal)                       Anneal hot temperature linearly to T=1.0.\n\
               (--annealStart N)                Iteration number to start annealing (5*10^5).\n\
               (--annealLength N)               Number of iterations to anneal all chains to T=1.0 (1*10^5).\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Output ----------------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--data-dump)                    Output waveforms to file.\n\
               (--adaptVerbose)                 Output parameter jump sizes and acceptance rate stats to file.\n\
               (--tempVerbose)                  Output temperature swapping stats to file.\n\
               (--propVerbose)                  Output proposal stats to file.\n\
               (--propTrack)                    Output useful info for track proposal behavior.\n\
               (--outfile file)                 Write output files <file>.<chain_number> (PTMCMC.output.<random_seed>.<chain_number>).\n";

  /* Print command line arguments if runState was not allocated */
  if(runState==NULL)
    {
      fprintf(stdout,"%s",help);
      return;
    }

  INT4 verbose=0,tmpi=0;
  unsigned int randomseed=0;
  REAL8 trigSNR = 0.0;
  REAL8 tempMin = 1.0;
  REAL8 tempMax = 50.0;
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
  if(LALInferenceGetProcParamVal(commandLine,"--propVerbose"))
    runState->proposalStats=XLALCalloc(1,sizeof(LALInferenceVariables));

  /* Set up the appropriate functions for the MCMC algorithm */
  runState->algorithm=&PTMCMCAlgorithm;
  runState->evolve=PTMCMCOneStep;

  ppt=LALInferenceGetProcParamVal(commandLine,"--rapidSkyLoc");
  if(ppt)
    runState->proposal=&LALInferenceRapidSkyLocProposal;
  else
    runState->proposal=&LALInferenceDefaultProposal;
    //runState->proposal=&LALInferencetempProposal;

 /* runState->template=&LALInferenceTemplateLAL;
  if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
    runState->template=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
    fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveform\"\n");
  }else{
    ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
    if(ppt){
      if(strstr(ppt->value,"TaylorF2") || strstr(ppt->value,"TaylorF2RedSpin")) {
        runState->template=&LALInferenceTemplateLAL;
        fprintf(stdout,"Template function called is \"LALInferenceTemplateLAL\"\n");
      }else if(strstr(ppt->value,"35phase_25amp")) {
        runState->template=&LALInferenceTemplate3525TD;
        fprintf(stdout,"Template function called is \"LALInferenceTemplate3525TD\"\n");
      }else{
        runState->template=&LALInferenceTemplateLALGenerateInspiral;
        fprintf(stdout,"Template function called is \"LALInferenceTemplateLALGenerateInspiral\"\n");
      }
    }
  }*/

//  if (LALInferenceGetProcParamVal(commandLine,"--tdlike")) {
//    fprintf(stderr, "Computing likelihood in the time domain.\n");
//    runState->likelihood=&LALInferenceTimeDomainLogLikelihood;
//  } else

  UINT4 malmquist = 0;
  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    runState->prior=&LALInferenceInspiralSkyLocPrior;
  } else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") || 
             LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
             LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
             LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
    runState->prior=&LALInferenceAnalyticNullPrior;
  } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
    runState->prior=&LALInferenceNullPrior;
  } else if (LALInferenceGetProcParamVal(commandLine, "--malmquistprior")) {
    printf("Using malmquist prior.\n");
    malmquist = 1;
    LALInferenceAddVariable(runState->priorArgs, "malmquist", &malmquist, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
    runState->prior=&LALInferenceInspiralPrior;
  } else {
    runState->prior=&LALInferenceInspiralPriorNormalised;
  }
  //runState->prior=PTUniformGaussianPrior;

  if (malmquist) {
      REAL8 malmquist_loudest = 0.0;
      REAL8 malmquist_second_loudest = 5.0;
      REAL8 malmquist_network = 0.0;

      ppt=LALInferenceGetProcParamVal(commandLine,"--malmquist-loudest-snr");
      if(ppt)
          malmquist_loudest = atof(ppt->value);

      ppt=LALInferenceGetProcParamVal(commandLine,"--malmquist-second-loudest-snr");
      if(ppt)
          malmquist_second_loudest = atof(ppt->value);

      ppt=LALInferenceGetProcParamVal(commandLine,"--malmquist-network-snr");
      if(ppt)
          malmquist_network = atof(ppt->value);

      LALInferenceAddVariable(runState->priorArgs, "malmquist_loudest_snr", &malmquist_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      LALInferenceAddVariable(runState->priorArgs, "malmquist_second_loudest_snr", &malmquist_second_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      LALInferenceAddVariable(runState->priorArgs, "malmquist_network_snr", &malmquist_network, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }

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
    tmpi=20000000;
  }
  LALInferenceAddVariable(runState->algorithmParams,"Niter",&tmpi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);

  printf("set iteration number between disk save.\n");
  /* Number of live points */
  ppt=LALInferenceGetProcParamVal(commandLine,"--Nskip");
  if(ppt)
    tmpi=atoi(ppt->value);
  else {
    tmpi=100;
  }
  LALInferenceAddVariable(runState->algorithmParams,"Nskip",&tmpi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);

 printf("set trigger SNR.\n");
  /* Network SNR of trigger */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigSNR");
  if(ppt){
    trigSNR=strtod(ppt->value,(char **)NULL);
  }
  LALInferenceAddVariable(runState->algorithmParams,"trigSNR",&trigSNR,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  printf("set lowest temperature.\n");
  /* Minimum temperature of the temperature ladder */
  ppt=LALInferenceGetProcParamVal(commandLine,"--tempMin");
  if(ppt){
    tempMin=strtod(ppt->value,(char **)NULL);
  }
  LALInferenceAddVariable(runState->algorithmParams,"tempMin",&tempMin, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

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

  /* Differential Evolution? */
  ppt=LALInferenceGetProcParamVal(commandLine, "--noDifferentialEvolution");
  if (!ppt) {
    fprintf(stderr, "Using differential evolution.\nEvery Nskip parameters will be stored for use in the d.e. jump proposal.\n");
    
    runState->differentialPoints = XLALCalloc(1, sizeof(LALInferenceVariables *));
    runState->differentialPointsLength = 0;
    runState->differentialPointsSize = 1;
    runState->differentialPointsSkip = 1;
  } else {
    fprintf(stderr, "Differential evolution disabled (--noDifferentialEvolution).\n");
    runState->differentialPoints = NULL;
    runState->differentialPointsLength = 0;
    runState->differentialPointsSize = 0;
    runState->differentialPointsSkip = 0;
  }
  
  INT4 Neff = 0;
  ppt = LALInferenceGetProcParamVal(commandLine, "--Neff");
  if (ppt)
    Neff = atoi(ppt->value);
  LALInferenceAddVariable(runState->algorithmParams, "Neff", &Neff, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
  
  return;

}

REAL8 **parseMCMCoutput(char ***params, UINT4 *nInPar, UINT4 *nInSamps, char *infileName, UINT4 burnin) {
    char str[999];
    char header[999];
    char *word;
    UINT4 nread;
    UINT4 i=0, j=0, nCols=0, nPar=0, par=0, col=0;
    UINT4 cycle=0;
    REAL8 val=0;

    const char *non_params[] = {"cycle","logpost","logprior","logl","loglH1","loglL1","loglV1","",NULL};

    FILE *infile = fopen(infileName,"r");

    fgets(str, 999, infile);
    strcpy(header, str);
    word = strtok(header, " \t");
    // Find column headers
    while (strcmp(word,"cycle") && str != NULL) {
        fgets(str, 999, infile);
        strcpy(header, str);
        word = strtok(header, " \t");
    }

    if (str == NULL) {
        fprintf(stderr, "Couldn't find column headers in file %s\n",infileName);
        exit(1);
    }

    // Read in column names and check if they are parameters
    strcpy(header, str);
    word = strtok(header, " \t");
    while (word != NULL) {
        nCols++;
        word = strtok(NULL, " \t");
    }
    // FIXME Remove a false column due to trailing whitespace
    nCols--;

    UINT4 is_param[nCols];

    strcpy(header, str);
    word = strtok(header, " \t");
    for (i=0; i<nCols; i++) {
        j=0;
        is_param[i] = 1;
        nPar++;
        while (non_params[j] != NULL) {
            if (!strcmp(non_params[j],word)) {
                is_param[i] = 0;
                nPar--;
                break;
            }
            j++;
        }
        word = strtok(NULL, " \t");
    }

    char** in_params = XLALMalloc((nPar)*sizeof(char *));

    word = strtok(str, " \t");
    // Already assumed cycle is the first column, so skip it
    par=0;
    for (i=1; i<nCols; i++) {
        char *param_name = strtok(NULL, " \t");
        if (is_param[i]) {
            in_params[par] = param_name;
            par++;
        }
    }

    printf("Reading the following params from %s:\n", infileName);
    for (par=0; par<nPar; par++)
        printf("\t%s\n",in_params[par]);

    // Move past burnin
    INT4 ch;
    if (burnin > 0) {
        while (cycle <= burnin) {
            fscanf(infile, "%i", &cycle);
            for (j=1;j<nCols;j++)
                fscanf(infile, "%lg", &val);
        }

        // Make sure at end of line
        ch = getc(infile);
        while (ch != '\n') ch = getc(infile);
    }

    // Determine number of samples after burnin
    unsigned long startPostBurnin = ftell(infile);
    UINT4 nSamples=0;

    while ( (ch = getc(infile)) != EOF) {
        if (ch=='\n')
            ++nSamples;
    }
    fseek(infile,startPostBurnin,SEEK_SET);
    printf("%i samples read from %s.\n", nSamples, infileName);

    // Read in samples
    REAL8 **sampleArray;
    sampleArray = (REAL8**) XLALMalloc(nSamples * sizeof(REAL8*));
    
    for (i = 0; i < nSamples; i++) {
        sampleArray[i] = XLALMalloc(nPar * sizeof(REAL8));

        nread = fscanf(infile, "%i", &cycle);
        if (nread != 1) {
            fprintf(stderr, "Cannot read sample from file (in %s, line %d)\n",
            __FILE__, __LINE__);
            exit(1);
        }

        par=0;
        for (col = 1; col < nCols; col++) {
            nread = fscanf(infile, "%lg", &val);
            if (nread != 1) {
                fprintf(stderr, "Cannot read sample from file (in %s, line %d)\n",
                __FILE__, __LINE__);
                exit(1);
            }

            if (is_param[col]) {
                sampleArray[i][par] = val;
                par++;
            }
        }
    }

    *params = in_params;
    *nInPar = nPar;
    *nInSamps = nSamples;
    return sampleArray;
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
  infileName = (char*)XLALCalloc(99,sizeof(char*));
  char str [999];
  FILE * infile;
  int n;
  char * pch;
  int fileargc = 1;
  char *fileargv[99];
  char buffer [99];

  /* Read command line and parse */
  procParams=LALInferenceParseCommandLine(argc,argv);

  ppt=LALInferenceGetProcParamVal(procParams,"--continue-run");
  if (ppt) {
    infileName = ppt->value;
    infile = fopen(infileName,"r");
    if (infile==NULL) {fprintf(stderr,"Cannot read %s/n",infileName); exit (1);}
    n=sprintf(buffer,"lalinference_mcmcmpi_from_file_%s",infileName);
    fileargv[0] = (char*)XLALCalloc((n+1),sizeof(char*));
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
            fileargv[fileargc] = (char*)XLALCalloc((n+1),sizeof(char*));
            fileargv[fileargc] = pch;
            fileargc++;
            if(fileargc>=99) {fprintf(stderr,"Too many arguments in file %s\n",infileName); exit (1);}
          }
        pch = strtok (NULL, " ");

      }
    fileargv[fileargc-1][strlen(fileargv[fileargc-1])-1]='\0'; //in order to get rid of the '\n' than fgets returns when reading the command line.

    procParams=LALInferenceParseCommandLine(fileargc,fileargv);
  }


  /* initialise runstate based on command line */
  /* This includes reading in the data */
  /* And performing any injections specified */
  /* And allocating memory */
  runState = initialize(procParams);

  /* Set up structures for MCMC */
  initializeMCMC(runState);
  if (runState)
    LALInferenceAddVariable(runState->algorithmParams,"MPIrank", &MPIrank, LALINFERENCE_UINT4_t,
                          LALINFERENCE_PARAM_FIXED);

  /* Set up model struct and set currentVariables to match the initialized model params */
  runState->model = LALInferenceInitCBCModel(runState);
  runState->currentParams = XLALMalloc(sizeof(LALInferenceVariables));
  memset(runState->currentParams, 0, sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(runState->model->params, runState->currentParams);

  /* Setup ROQ */
  fprintf(stdout, " ==== LALInferenceSetupROQ(): started. ====\n");
  LALInferenceSetupROQ(runState->data, runState->model, procParams);
  fprintf(stdout, " ==== LALInferenceSetupROQ(): finished. ====\n");

  /* Set template function in runState, since it's sometime used */
  runState->templt = runState->model->templt;

  /* Choose the likelihood */
  LALInferenceInitLikelihood(runState);
 
  /* Call the extra code that was removed from previous function */
  LALInferenceInitMCMCState(runState);
 
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
