/*
 *  LALInferenceMCMC.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch and Will M. Farr
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
#include <lalapps.h>

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
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Proposals  ------------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--rapidSkyLoc)                  Use rapid sky localization jump proposals.\n\
               (--kDTree)                       Use a kDTree proposal.\n\
               (--kDNCell N)                    Number of points per kD cell in proposal.\n\
               (--covarianceMatrix file)        Find the Cholesky decomposition of the covariance matrix for jumps in file.\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Parallel Tempering Algorithm Parameters -------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--inverseLadder)                Space temperature uniform in 1/T, rather than geometric.\n\
               (--tempSkip N)                   Number of iterations between temperature swap proposals (100).\n\
               (--tempSwaps N)                  Number of random swaps proposed every <tempSkip> iterations ((nTemps-1)nTemps/2).\n\
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

  /* Choose the template generator for inspiral signals */
  runState->template=&LALInferenceTemplateLAL;
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
  }

//  if (LALInferenceGetProcParamVal(commandLine,"--tdlike")) {
//    fprintf(stderr, "Computing likelihood in the time domain.\n");
//    runState->likelihood=&LALInferenceTimeDomainLogLikelihood;
//  } else 
  if (LALInferenceGetProcParamVal(commandLine, "--zeroLogLike")) {
    /* Use zero log(L) */
    runState->likelihood=&LALInferenceZeroLogLikelihood;
  } else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood")) {
    runState->likelihood=&LALInferenceCorrelatedAnalyticLogLikelihood;
  } else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood")) {
    runState->likelihood=&LALInferenceBimodalCorrelatedAnalyticLogLikelihood;
  } else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood")) {
    runState->likelihood=&LALInferenceRosenbrockLogLikelihood;
  } else if (LALInferenceGetProcParamVal(commandLine, "--studentTLikelihood")) {
    fprintf(stderr, "Using Student's T Likelihood.\n");
    runState->likelihood=&LALInferenceFreqDomainStudentTLogLikelihood;
  } else {
    runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
  }

  if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
    runState->prior=&LALInferenceInspiralSkyLocPrior;
  } else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") || 
             LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
             LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
             LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
    runState->prior=&LALInferenceAnalyticNullPrior;
  } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
    runState->prior=&LALInferenceNullPrior;
  } else {
    runState->prior=&LALInferenceInspiralPriorNormalised;
  }
  //runState->prior=PTUniformGaussianPrior;

  ppt=LALInferenceGetProcParamVal(commandLine,"--verbose");
  if(ppt) {
    verbose=1;
    LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose , LALINFERENCE_UINT4_t,
                            LALINFERENCE_PARAM_FIXED);
    set_debug_level("ERROR|INFO");
  }
  else set_debug_level("NDEBUG");

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
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Injection Arguments ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--inj injections.xml)          Injection XML file to use.\n\
               (--event N)                     Event number from Injection XML file to use.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Template Arguments -------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--symMassRatio)                Jump in symmetric mass ratio eta, instead of q=m2/m1.\n\
               (--LALSimulation)               Interface with the LALSimulation package for template generation.\n\
               (--approximant Approximant)     Specify a template approximant to use (default TaylorF2). Possible values are: \n\
                                               default modeldomain=\"time\": GeneratePPN, TaylorT1, TaylorT2, TaylorT3, TaylorT4, \n\
                                                                           EOB, EOBNR, EOBNRv2, EOBNRv2HM, SpinTaylor, \n\
                                                                           SpinQuadTaylor, SpinTaylorFrameless, SpinTaylorT4, \n\
                                                                           PhenSpinTaylorRD, NumRel.\n\
                                               default modeldomain=\"frequency\": TaylorF1, TaylorF2, TaylorF2RedSpin, \n\
                                                                                TaylorF2RedSpinTidal, IMRPhenomA, IMRPhenomB.\n\
               (--order PNorder)               Specify a PN order in phase to use (default threePointFivePN).\n\
               (--ampOrder PNorder)            Specify a PN order in amplitude to use (defaults: LALSimulation: max available; LALInspiral: newtownian).\n\
               (--fref fRef)                   Specify a reference frequency at which parameters are defined (default 0).\n\
               (--tidal)                       Enables tidal corrections, only with LALSimulation.\n\
               (--interactionFlags)            intercation flags, only with LALSimuation (LAL_SIM_INSPIRAL_INTERACTION_ALL).\n\
               (--modeldomain)                 domain the waveform template will be computed in (\"time\" or \"frequency\").\n\
               (--spinAligned)                 template will assume spins aligned with the orbital angular momentum.\n\
                                               *Enables* spins for TaylorF2, TaylorF2RedSpin, TaylorF2RedSpinTidal, IMRPhenomB.\n\
               (--singleSpin)                  template will assume only the spin of the most massive binary component exists.\n\
               (--noSpin)                      template will assume no spins.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Starting Parameters ------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--trigtime time)               Trigger time to use.\n\
               (--time time)                   Waveform time (overrides random about trigtime).\n\
               (--mc mchirp)                   Trigger chirpmass to use.\n\
               (--eta eta)                     Trigger eta (symmetric mass ratio) to use.\n\
               (--q q)                         Trigger q (asymmetric mass ratio) to use.\n\
               (--phi phase)                   Trigger phase to use.\n\
               (--iota inclination)            Trigger inclination to use.\n\
               (--dist dist)                   Trigger distance.\n\
               (--ra ra)                       Trigger RA.\n\
               (--dec dec)                     Trigger declination.\n\
               (--psi psi)                     Trigger psi.\n\
               (--a1 a1)                       Trigger a1.\n\
               (--theta1 theta1)               Trigger theta1.\n\
               (--phi1 phi1)                   Trigger phi1.\n\
               (--a2 a2)                       Trigger a2.\n\
               (--theta2 theta2)               Trigger theta2.\n\
               (--phi2 phi2)                   Trigger phi2.\n\
               (--lambda1)                     Trigger lambda1.\n\
               (--lambda2)                     Trigger lambda2.\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Prior Arguments ----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--mc-min mchirp)               Minimum chirp mass.\n\
               (--mc-max mchirp)               Maximum chirp mass.\n\
               (--eta-min etaMin)              Minimum eta.\n\
               (--eta-max etaMax)              Maximum eta.\n\
               (--q-min qMin)                  Minimum q.\n\
               (--q-max qMax)                  Maximum q.\n\
               (--comp-min min)                Minimum component mass (1.0).\n\
               (--comp-max max)                Maximum component mass (30.0).\n\
               (--MTotMax max)                 Maximum total mass (35.0).\n\
               (--iota-max max)                Maximum inclination angle (pi).\n\
               (--Dmin dist)                   Minimum distance in Mpc (1).\n\
               (--Dmax dist)                   Maximum distance in Mpc (100).\n\
               (--lambda1-min)                 Minimum lambda1 (0).\n\
               (--lambda1-max)                 Maximum lambda1 (3000).\n\
               (--lambda2-min)                 Minimum lambda2 (0).\n\
               (--lambda2-max)                 Maximum lambda2 (3000).\n\
               (--dt time)                     Width of time prior, centred around trigger (0.1s).\n\
               \n\
               ------------------------------------------------------------------------------------------------------------------\n\
               --- Fix Parameters -----------------------------------------------------------------------------------------------\n\
               ------------------------------------------------------------------------------------------------------------------\n\
               (--fixMc)                       Do not allow chirpmass to vary.\n\
               (--fixEta)                      Do not allow mass ratio to vary.\n\
               (--fixQ)                        Do not allow mass ratio to vary.\n\
               (--fixPhi)                      Do not allow phase to vary.\n\
               (--fixIota)                     Do not allow inclination to vary.\n\
               (--fixDist)                     Do not allow distance to vary.\n\
               (--fixRa)                       Do not allow RA to vary.\n\
               (--fixDec)                      Do not allow declination to vary.\n\
               (--fixPsi)                      Do not allow polarization to vary.\n\
               (--fixA1)                       Do not allow spin to vary.\n\
               (--fixTheta1)                   Do not allow spin 1 colatitude to vary.\n\
               (--fixPhi1)                     Do not allow spin 1 longitude to vary.\n\
               (--fixA2)                       Do not allow spin 2 to vary.\n\
               (--fixTheta2)                   Do not allow spin 2 colatitude to vary.\n\
               (--fixPhi2)                     Do not allow spin 2 longitude to vary.\n\
               (--fixTime)                     Do not allow coalescence time to vary.\n\
               (--fixLambda1)                  Do not allow lambda1 to vary.\n\
               (--fixLambda2)                  Do not allow lambda2 to vary.\n";

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
  //INT4 AmpOrder=0;
  LALPNOrder PhaseOrder=LAL_PNORDER_THREE_POINT_FIVE;
  LALPNOrder AmpOrder=-1;//LAL_PNORDER_THREE_POINT_FIVE;//LAL_PNORDER_NEWTONIAN;
  Approximant approx=TaylorF2;
  REAL8 fRef = 0.0;
  LALInferenceApplyTaper bookends = LALINFERENCE_TAPER_NONE;
  UINT4 analytic=0;
  LALInferenceIFOData *dataPtr;
  LALInferenceDomain modelDomain;
  UINT4 event=0;
  UINT4 i=0;
  REAL8 m1=0;
  REAL8 m2=0;
  REAL8 logDmin=log(1.0);
  REAL8 logDmax=log(100.0);
  REAL8 Dmin=1.0;
  REAL8 Dmax=100.0;
  REAL8 mcMin=1.0;
  REAL8 mcMax=15.3;
  REAL8 mMin=1.0,mMax=30.0;
  REAL8 MTotMax=35.0;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 qMin=mMin/mMax;
  REAL8 qMax=1.0;
  REAL8 m1min=mMin,m1max=mMax;
  REAL8 m2min=mMin,m2max=mMax;
  REAL8 iotaMin=0.0,iotaMax=LAL_PI;
  REAL8 psiMin=0.0,psiMax=LAL_PI;
  REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0;
  REAL8 raMin=0.0,raMax=LAL_TWOPI;
  REAL8 phiMin=0.0,phiMax=LAL_TWOPI;
  REAL8 a1min=0.0,a1max=1.0;
  REAL8 a2min=0.0,a2max=1.0;
  REAL8 theta1min=0.0,theta1max=LAL_PI;
  REAL8 theta2min=0.0,theta2max=LAL_PI;
  REAL8 phi1min=0.0,phi1max=LAL_TWOPI;
  REAL8 phi2min=0.0,phi2max=LAL_TWOPI;
  REAL8 dt=0.1;            /* Width of time prior */
  REAL8 lambda1Min=0.0;
  REAL8 lambda1Max=3000.0;
  REAL8 lambda2Min=0.0;
  REAL8 lambda2Max=3000.0;  
  REAL8 tmpMin,tmpMax;//,tmpVal;
  gsl_rng * GSLrandom=state->GSLrandom;
  REAL8 endtime=0.0, timeParam=0.0;
  REAL8 timeMin=endtime-dt,timeMax=endtime+dt;
  REAL8 start_mc			=4.82+gsl_ran_gaussian(GSLrandom,0.025);
  REAL8 start_eta			=etaMin+gsl_rng_uniform(GSLrandom)*(etaMax-etaMin);
  REAL8 start_q           =qMin+gsl_rng_uniform(GSLrandom)*(qMax-qMin);
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
  REAL8 start_lambda1 =lambda1Min+gsl_rng_uniform(GSLrandom)*(lambda1Max-lambda1Min);
  REAL8 start_lambda2 =lambda2Min+gsl_rng_uniform(GSLrandom)*(lambda2Max-lambda2Min);
  
  memset(currentParams,0,sizeof(LALInferenceVariables));

  /* Over-ride prior bounds if analytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood")) {
    analytic  = 1;
    m1min     = 14.927715;
    m1max     = 17.072285;
    m2min     = 5.829675;
    m2max     = 8.170325;
    iotaMin   = 1.4054428267948966;
    iotaMax   = 1.7361498267948965;
    phiMin    = 2.8701521535897934;
    phiMax    = 3.413033153589793;
    psiMin    = 1.3885563267948966;
    psiMax    = 1.7530363267948965;
    raMin     = 2.813050153589793;
    raMax     = 3.4701351535897933;
    decMin    = -0.300699;
    decMax    = 0.300699;
    Dmin      = 37.986000000000004;
    Dmax      = 62.013999999999996;
    timeMin   = -0.1073625;
    timeMax   = 0.1073625;
    a1min     = 0.3784565;
    a1max     = 0.6215435;
    a2min     = 0.421869;
    a2max     = 0.578131;
    theta1min = 1.3993998267948966;
    theta1max = 1.7421928267948965;
    theta2min = 1.4086158267948965;
    theta2max = 1.7329768267948966;
    phi1min   = 2.781852653589793;
    phi1max   = 3.501332653589793;
    phi2min   = 2.777215653589793;
    phi2max   = 3.5059696535897933;
  } else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood")) {
    analytic  = 1;
    m1min     = 14.927715;
    m1max     = 18.787941;
    m2min     = 5.829675;
    m2max     = 10.042845;
    iotaMin   = 0.6200446634;
    iotaMax   = 1.2153172634;
    phiMin    = 1.2993558268;
    phiMax    = 2.2765416268;
    psiMin    = 0.6031581634;
    psiMax    = 1.2592221634;
    raMin     = 1.2422538268;
    raMax     = 2.4250068268;
    decMin    = -1.0860971634;
    decMax    = -0.0035807634;
    Dmin      = 12.986;
    Dmax      = 56.2364;
    timeMin   = -0.1373625;
    timeMax   = 0.2491425;
    a1min     = 0.0784565;
    a1max     = 0.5160131;
    a2min     = 0.121869;
    a2max     = 0.4031406;
    theta1min = 0.6140016634;
    theta1max = 1.2310290634;
    theta2min = 0.6232176634;
    theta2max = 1.2070674634;
    phi1min   = 1.2110563268;
    phi1max   = 2.5061203268;
    phi2min   = 1.2064193268;
    phi2max   = 2.5181765268;
  } else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood")) {
    analytic  = 1;
    m1min     = 14.0;
    m1max     = 18.0;
    m2min     = 5.0;
    m2max     = 9.0;
    iotaMin   = -0.429203673;
    iotaMax   = 3.570796327;
    phiMin    = 1.141592654;
    phiMax    = 5.141592654;
    psiMin    = -0.429203673;
    psiMax    = 3.570796327;
    raMin     = 1.141592654;
    raMax     = 5.141592654;
    decMin    = -2.0;
    decMax    = 2.0;
    Dmin      = 48.0;
    Dmax      = 52.0;
    timeMin   = -2.0;
    timeMax   = 2.0;
    a1min     = -1.5;
    a1max     = 2.5;
    a2min     = -1.5;
    a2max     = 2.5;
    theta1min = -0.429203673;
    theta1max = 3.570796327;
    theta2min = -0.429203673;
    theta2max = 3.570796327;
    phi1min   = 1.141592654;
    phi1max   = 5.141592654;
    phi2min   = 1.141592654;
    phi2max   = 5.141592654;
  }


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
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
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
    }
  }

  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    if ( ! strcmp( "GeneratePPN", ppt->value ) )
      {
        approx = GeneratePPN;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "TaylorT1", ppt->value ) )
      {
        approx = TaylorT1;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "TaylorT2", ppt->value ) )
      {
        approx = TaylorT2;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "TaylorT3", ppt->value ) )
      {
        approx = TaylorT3;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "TaylorT4", ppt->value ) )
      {
        approx = TaylorT4;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "TaylorF1", ppt->value ) )
      {
        approx = TaylorF1;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
    else if ( ! strcmp( "TaylorF2", ppt->value ) )
      {
        approx = TaylorF2;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
    else if ( ! strcmp( "TaylorF2RedSpin", ppt->value ) )
      {
        approx = TaylorF2RedSpin;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
      else if ( ! strcmp( "TaylorF2RedSpinTidal", ppt->value ) )
      {
        approx = TaylorF2RedSpinTidal;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
    else if ( ! strcmp( "EOB", ppt->value ) )
      {
        approx = EOB;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "EOBNR", ppt->value ) )
      {
        approx = EOBNR;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "EOBNRv2", ppt->value ) )
      {
        approx = EOBNRv2;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "EOBNRv2HM", ppt->value ) )
      {
        approx = EOBNRv2HM;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "SpinTaylor", ppt->value ) )
      {
        approx = SpinTaylor;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "SpinTaylorT4", ppt->value ) )
      {
        approx = SpinTaylorT4;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "SpinQuadTaylor", ppt->value ) )
      {
        approx = SpinQuadTaylor;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "SpinTaylorFrameless", ppt->value ) )
      {
        approx = SpinTaylorFrameless;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "PhenSpinTaylorRD", ppt->value ) )
      {
        approx = PhenSpinTaylorRD;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "NumRel", ppt->value ) )
      {
        approx = NumRel;
        modelDomain = LALINFERENCE_DOMAIN_TIME;
      }
    else if ( ! strcmp( "IMRPhenomA", ppt->value ) )
      {
        approx = IMRPhenomA;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
    else if ( ! strcmp( "IMRPhenomB", ppt->value ) )
      {
        approx = IMRPhenomB;
        modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
      }
    else
      {
        fprintf( stderr, "invalid argument to --approximant\n"
                 "unknown approximant %s specified: "
                 "Approximant must be one of: GeneratePPN, TaylorT1, TaylorT2,\n"
                 "TaylorT3, TaylorT4, TaylorF1, TaylorF2, TaylorF2RedSpin, TaylorF2RedSpinTidal, EOB, EOBNR, EOBNRv2, \n"
                 "EOBNRv2HM, SpinTaylor, SpinQuadTaylor, SpinTaylorFrameless, SpinTaylorT4\n"
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
    else if ( ! strcmp( "twoPointFivePN", ppt->value ) )
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
    fprintf(stdout,"Templates will be generated at %.1f PN order in phase\n",((float)(PhaseOrder))/2.0);
  }
  
  /* Over-ride amplitude PN order if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--ampOrder");
  if(ppt){
    if ( ! strcmp( "newtonian", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_NEWTONIAN;
      }
    else if ( ! strcmp( "oneHalfPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_HALF;
      }
    else if ( ! strcmp( "onePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_ONE;
      }
    else if ( ! strcmp( "onePointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_ONE_POINT_FIVE;
      }
    else if ( ! strcmp( "twoPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_TWO;
      }
    else if ( ! strcmp( "twoPointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_TWO_POINT_FIVE;
      }
    else if ( ! strcmp( "threePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_THREE;
      }
    else if ( ! strcmp( "threePointFivePN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_THREE_POINT_FIVE;
      }
    else if ( ! strcmp( "pseudoFourPN", ppt->value ) )
      {
        AmpOrder = LAL_PNORDER_PSEUDO_FOUR;
      }
    else
      {
        fprintf( stderr, "invalid argument to --ampOrder:\n"
                 "unknown order specified: "
                 "PN order must be one of: newtonian, oneHalfPN, onePN,\n"
                 "onePointFivePN, twoPN, twoPointFivePN, threePN or\n"
                 "threePointFivePN\n");
        exit( 1 );
      }
    fprintf(stdout,"Templates will be generated at %.1f PN order in amplitude\n",((float)(AmpOrder))/2.0);
  }
  
  ppt=LALInferenceGetProcParamVal(commandLine, "--fref");
  if (ppt) fRef = atof(ppt->value);

  ppt=LALInferenceGetProcParamVal(commandLine,"--modeldomain");
  if(ppt){
    if ( ! strcmp( "time", ppt->value ) )
    {
      modelDomain = LALINFERENCE_DOMAIN_TIME;
    }
    else if ( ! strcmp( "frequency", ppt->value ) )
    {
      modelDomain = LALINFERENCE_DOMAIN_FREQUENCY;
    }
    else
    {
      fprintf( stderr, "invalid argument to --modeldomain:\n"
              "unknown domain specified: "
              "domain must be one of: time, frequency\n");
      exit( 1 );
    }
  }
  
  dataPtr = state->data;
  while (dataPtr != NULL) {
    dataPtr->modelDomain = modelDomain;
    dataPtr = dataPtr->next;
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
  if(ppt && !analytic){
    endtime=atof(ppt->value);
    timeMin=endtime-dt; timeMax=endtime+dt;
    printf("Read end time %f\n",endtime);

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
    LALInferenceMcEta2Masses(start_mc, start_eta, &m1, &m2);
    start_q=m2/m1;
  }

  /* Over-ride q if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--q");
  if(ppt){
    start_q=atof(ppt->value);
    LALInferenceQ2Eta(start_q, &start_eta);
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

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1");
  if (ppt) {
    start_lambda1 = atof(ppt->value);
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2");
  if (ppt) {
    start_lambda2 = atof(ppt->value);
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

  /* Over-ride lambda1 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-min");
  if(ppt){
    lambda1Min=atof(ppt->value);
  }
  
  /* Over-ride lambda1 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda1-max");
  if(ppt){
    lambda1Max=atof(ppt->value);
  }
  
  /* Over-ride lambda2 min if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-min");
  if(ppt){
    lambda2Min=atof(ppt->value);
  }
  
  /* Over-ride lambda2 max if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--lambda2-max");
  if(ppt){
    lambda2Max=atof(ppt->value);
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

  ppt=LALInferenceGetProcParamVal(commandLine,"--iota-max");
  if (ppt) {
    iotaMax = atof(ppt->value);
  }


  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(currentParams, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  if(LALInferenceGetProcParamVal(commandLine,"--ampOrder")) 
    LALInferenceAddVariable(currentParams, "LAL_AMPORDER",     &AmpOrder,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  if (LALInferenceGetProcParamVal(commandLine, "--fref"))
    LALInferenceAddVariable(currentParams, "fRef", &fRef, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

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
  /* Jump in component masses if anaytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") || 
      LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
      LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
      LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
    LALInferenceMcEta2Masses(start_mc, start_eta, &m1, &m2);
    LALInferenceAddVariable(currentParams, "m1",    &m1,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs,	"m1",	&m1min,	&m1max,		LALINFERENCE_REAL8_t);
    LALInferenceAddVariable(currentParams, "m2",    &m2,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddMinMaxPrior(priorArgs,	"m2",	&m2min,	&m2max,		LALINFERENCE_REAL8_t);
  } else {

    ppt=LALInferenceGetProcParamVal(commandLine,"--fixMc");
    if(ppt){
      LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"chirpmass fixed and set to %f\n",start_mc);
    }else{
      LALInferenceAddVariable(currentParams, "chirpmass",    &start_mc,    LALINFERENCE_REAL8_t,	LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs,	"chirpmass",	&mcMin,	&mcMax,		LALINFERENCE_REAL8_t);

    /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixQ");
    if(ppt){
      LALInferenceAddVariable(currentParams, "asym_massratio", &start_q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
      if(MPIrank==0) fprintf(stdout,"q fixed and set to %f\n",start_q);
    }else{
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixEta");
      if(ppt){
        LALInferenceAddVariable(currentParams, "massratio", &start_eta, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddMinMaxPrior(priorArgs,	"massratio", &etaMin, &etaMax, LALINFERENCE_REAL8_t);
        if(MPIrank==0) fprintf(stdout,"eta fixed and set to %f\n",start_eta);
      }else{
        ppt=LALInferenceGetProcParamVal(commandLine,"--symMassRatio");
        if(ppt){
          LALInferenceAddVariable(currentParams, "massratio", &start_eta, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
          LALInferenceAddMinMaxPrior(priorArgs,	"massratio", &etaMin, &etaMax, LALINFERENCE_REAL8_t);
        }else{
          LALInferenceAddVariable(currentParams, "asym_massratio", &start_q, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
          LALInferenceAddMinMaxPrior(priorArgs,	"asym_massratio",	&qMin,	&qMax,	LALINFERENCE_REAL8_t);
        }
      }
    }
  }

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
  LALInferenceAddMinMaxPrior(priorArgs, "time",     &timeMin, &timeMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"phase fixed and set to %f\n",start_phase);
  }else{
    LALInferenceAddVariable(currentParams, "phase",           &start_phase,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "phase",     &phiMin, &phiMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDist");
  if(ppt){
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"distance fixed and set to %f\n",start_dist);
  }else{
    LALInferenceAddVariable(currentParams,"distance", &start_dist, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "distance",     &Dmin, &Dmax,   LALINFERENCE_REAL8_t);


  ppt=LALInferenceGetProcParamVal(commandLine,"--fixRa");
  if(ppt){
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"R.A. fixed and set to %f\n",start_ra);
  }else{
    LALInferenceAddVariable(currentParams, "rightascension",  &start_ra,      LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "rightascension",     &raMin, &raMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixDec");
  if(ppt){
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"declination fixed and set to %f\n",start_dec);
  }else{
    LALInferenceAddVariable(currentParams, "declination",     &start_dec,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "declination",     &decMin, &decMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixPsi");
  if(ppt){
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"polarisation fixed and set to %f\n",start_psi);
  }else{
    LALInferenceAddVariable(currentParams, "polarisation",    &start_psi,     LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "polarisation",     &psiMin, &psiMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine,"--fixIota");
  if(ppt){
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    if(MPIrank==0) fprintf(stdout,"iota fixed and set to %f\n",start_iota);
  }else{
    LALInferenceAddVariable(currentParams, "inclination",     &start_iota,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
  LALInferenceAddMinMaxPrior(priorArgs, "inclination",     &iotaMin, &iotaMax,   LALINFERENCE_REAL8_t);

  ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
  if((approx==SpinTaylor || approx==SpinTaylorFrameless || approx==PhenSpinTaylorRD || approx==SpinTaylorT4) && !ppt){


    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    if(ppt) a1min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "a_spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "a_spin1",     &a1min, &a1max,   LALINFERENCE_REAL8_t);

    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    if(ppt) fprintf(stdout,"Running with spin1 aligned to the orbital angular momentum.\n");
    else {
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(MPIrank==0) fprintf(stdout,"theta 1 fixed and set to %f\n",start_theta_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "theta_spin1",     &start_theta_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1",     &theta1min, &theta1max,   LALINFERENCE_REAL8_t);

      ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi1");
      if(ppt){
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(MPIrank==0) fprintf(stdout,"phi 1 fixed and set to %f\n",start_phi_spin1);
      }else{
        LALInferenceAddVariable(currentParams, "phi_spin1",     &start_phi_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "phi_spin1",     &phi1min, &phi1max,   LALINFERENCE_REAL8_t);
    }
    ppt=LALInferenceGetProcParamVal(commandLine, "--singleSpin");
    if(ppt) fprintf(stdout,"Running with first spin set to 0\n");
    else {
      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      if(ppt) a2min=-1.0;
      ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
      if(ppt){
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
      }else{
        LALInferenceAddVariable(currentParams, "a_spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      LALInferenceAddMinMaxPrior(priorArgs, "a_spin2",     &a2min, &a2max,   LALINFERENCE_REAL8_t);

      ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
      if(ppt) fprintf(stdout,"Running with spin2 aligned to the orbital angular momentum.\n");
      else {
        ppt=LALInferenceGetProcParamVal(commandLine,"--fixTheta2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          if(MPIrank==0) fprintf(stdout,"theta spin 2 fixed and set to %f\n",start_theta_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "theta_spin2",     &start_theta_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2",     &theta2min, &theta2max,   LALINFERENCE_REAL8_t);

        ppt=LALInferenceGetProcParamVal(commandLine,"--fixPhi2");
        if(ppt){
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
          if(MPIrank==0) fprintf(stdout,"phi 2 fixed and set to %f\n",start_phi_spin2);
        }else{
          LALInferenceAddVariable(currentParams, "phi_spin2",     &start_phi_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
        }
        LALInferenceAddMinMaxPrior(priorArgs, "phi_spin2",     &phi2min, &phi2max,   LALINFERENCE_REAL8_t);
      }
    }
  }
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
  if((approx==TaylorF2 || approx==TaylorF2RedSpin || approx==TaylorF2RedSpinTidal || approx==IMRPhenomB) && ppt){

    a1min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"spin 1 fixed and set to %f\n",start_a_spin1);
    }else{
      LALInferenceAddVariable(currentParams, "spin1",     &start_a_spin1,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin1",     &a1min, &a1max,   LALINFERENCE_REAL8_t);

    a2min=-1.0;
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixA2");
    if(ppt){
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"spin 2 fixed and set to %f\n",start_a_spin2);
    }else{
      LALInferenceAddVariable(currentParams, "spin2",     &start_a_spin2,            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "spin2",     &a2min, &a2max,   LALINFERENCE_REAL8_t);

  }

  if (LALInferenceGetProcParamVal(commandLine, "--studentTLikelihood")) {
    
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
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--tidal");
  if(ppt){
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda1");
    if(ppt){
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"phase fixed and set to %f\n",start_lambda1);
    }else{
      LALInferenceAddVariable(currentParams, "lambda1",           &start_lambda1,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda1",     &lambda1Min, &lambda1Max,   LALINFERENCE_REAL8_t);
  
    ppt=LALInferenceGetProcParamVal(commandLine,"--fixLambda2");
    if(ppt){
    LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      if(MPIrank==0) fprintf(stdout,"phase fixed and set to %f\n",start_lambda2);
    }else{
      LALInferenceAddVariable(currentParams, "lambda2",           &start_lambda2,        LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    }
    LALInferenceAddMinMaxPrior(priorArgs, "lambda2",     &lambda2Min, &lambda2Max,   LALINFERENCE_REAL8_t);
  }
  
  LALSimInspiralInteraction interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine,"--interactionFlags");
  if(ppt){
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_NONE")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_NONE;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_15PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_SPIN_SELF_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_QUAD_MONO_2PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_SPIN_ORBIT_25PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_5PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_TIDAL_6PN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL_SPIN;
    if(strstr(ppt->value,"LAL_SIM_INSPIRAL_INTERACTION_ALL")) interactionFlags=LAL_SIM_INSPIRAL_INTERACTION_ALL;
    LALInferenceAddVariable(currentParams, "interactionFlags", &interactionFlags,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED); 
 }
 
  

  /* Initialize variable that will store the name of the last proposal function used */
  const char *initPropName = "INITNAME";
  LALInferenceAddVariable(state->proposalArgs, LALInferenceCurrentProposalName, &initPropName, LALINFERENCE_string_t, LALINFERENCE_PARAM_LINEAR);

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
  char buffer [99];

  /* Read command line and parse */
  procParams=LALInferenceParseCommandLine(argc,argv);

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

  /* Set up currentParams with variables to be used */
  initVariables(runState);
  
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
