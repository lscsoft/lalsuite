/*
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
 *
 *  Copyright (C) 2009, 2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch, Will M. Farr, and Ben Farr
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
#include <stdlib.h>
#include <lal/LALInspiral.h>
#include <lal/DetResponse.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <mpi.h>
#include <lal/LALInference.h>
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <sys/time.h>

#include <LALAppsVCSInfo.h>
#include <lal/LALStdlib.h>

#define PROGRAM_NAME "LALInferenceMCMCSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

const char *const parallelSwapProposalName = "ParallelSwap";
const char *const clusteredKDEProposalName = "ClusteredKDEProposal";

static void
thinDifferentialEvolutionPoints(LALInferenceRunState *runState) {
  size_t i;
  size_t newSize;
  
  /* Delete all the even-index points. */
  for (i = 0; i < runState->differentialPointsLength; i += 2) {
    LALInferenceClearVariables(runState->differentialPoints[i]);
    XLALFree(runState->differentialPoints[i]);
    runState->differentialPoints[i] = NULL;
  }
  
  /* Copy the odd points into the first part of the array. */
  for (i = 1; i < runState->differentialPointsLength; i += 2) {
    runState->differentialPoints[i/2] = runState->differentialPoints[i];
    runState->differentialPoints[i] = NULL;
  }

  newSize = runState->differentialPointsLength / 2;

  /* Now shrink the buffer down. */
  runState->differentialPoints = XLALRealloc(runState->differentialPoints, 2*newSize*sizeof(LALInferenceVariables *));
  runState->differentialPointsSize = 2*newSize;
  runState->differentialPointsLength = newSize;
  runState->differentialPointsSkip *= 2;
}

static void
accumulateDifferentialEvolutionSample(LALInferenceRunState *runState) {
  if (runState->differentialPointsSize == runState->differentialPointsLength) {
    size_t newSize = runState->differentialPointsSize*2;
    ProcessParamsTable *ppt = LALInferenceGetProcParamVal(runState->commandLine, "--differential-buffer-limit");

    if (ppt && (size_t)atoi(ppt->value) < newSize) {
      /* Then thin, and record sample. */
      thinDifferentialEvolutionPoints(runState);
      return accumulateDifferentialEvolutionSample(runState);
    } else {
      runState->differentialPoints = XLALRealloc(runState->differentialPoints, newSize*sizeof(LALInferenceVariables *));
      runState->differentialPointsSize = newSize;
    }
  }

  runState->differentialPoints[runState->differentialPointsLength] = XLALCalloc(1, sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(runState->currentParams, runState->differentialPoints[runState->differentialPointsLength]);
  runState->differentialPointsLength += 1;
}

static void
resetDifferentialEvolutionBuffer(LALInferenceRunState *runState) {
  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  size_t i;

  for (i = 0; i < runState->differentialPointsLength; i++) {
    LALInferenceClearVariables(runState->differentialPoints[i]);
    XLALFree(runState->differentialPoints[i]);
    runState->differentialPoints[i] = NULL;
  }

  runState->differentialPoints = XLALRealloc(runState->differentialPoints, 1*sizeof(LALInferenceVariables *));
  runState->differentialPointsLength = 0;
  runState->differentialPointsSize = 1;
  runState->differentialPointsSkip = Nskip;
}

static void
accumulateKDTreeSample(LALInferenceRunState *runState) {
  LALInferenceVariables *proposalParams = runState->proposalArgs;

  if (!LALInferenceCheckVariable(proposalParams, "kDTree") || !LALInferenceCheckVariable(proposalParams, "kDTreeVariableTemplate")) {
    /* Not setup correctly. */
    return;
  }

  LALInferenceKDTree *tree = *(LALInferenceKDTree **)LALInferenceGetVariable(proposalParams, "kDTree");
  LALInferenceVariables *template = *(LALInferenceVariables **)LALInferenceGetVariable(proposalParams, "kDTreeVariableTemplate");
  size_t ndim = LALInferenceGetVariableDimensionNonFixed(template);
  REAL8 *pt = XLALMalloc(ndim*sizeof(REAL8));

  LALInferenceKDVariablesToREAL8(runState->currentParams, pt, template);

  LALInferenceKDAddPoint(tree, pt);

  /* Don't free pt, since it is now stored in tree. */
}

static void
BcastDifferentialEvolutionPoints(LALInferenceRunState *runState, INT4 sourceTemp) {
  INT4 MPIrank;
  UINT4 i=0;
  REAL8** packedDEsamples;
  REAL8*  temp;

  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  UINT4 nPar = (UINT4)LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  UINT4 nPoints = (UINT4)runState->differentialPointsLength;

  /* Prepare a DE buffer of the proper size */
  MPI_Bcast(&nPoints, 1, MPI_INT, sourceTemp, MPI_COMM_WORLD);

  /* Prepare 2D array for DE points */
  packedDEsamples = (REAL8**) XLALMalloc(nPoints * sizeof(REAL8*));
  temp = (REAL8*) XLALMalloc(nPoints * nPar * sizeof(REAL8));
  for (i=0; i < nPoints; i++)
    packedDEsamples[i] = temp + (i*nPar);

  /* Pack it up */
  if (MPIrank == sourceTemp)
    LALInferenceBufferToArray(runState, packedDEsamples);

  /* Send it out */
  MPI_Bcast(packedDEsamples[0], nPoints*nPar, MPI_DOUBLE, sourceTemp, MPI_COMM_WORLD);

  /* Unpack it */
  if (MPIrank != sourceTemp)
    LALInferenceArrayToBuffer(runState, packedDEsamples, nPoints);

  /* Clean up */
  XLALFree(temp);
  XLALFree(packedDEsamples);
  MPI_Barrier(MPI_COMM_WORLD);
}

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  INT4 i,t; //indexes for for() loops
  INT4 nChain;
  INT4 MPIrank, MPIsize;
  MPI_Status MPIstatus;
  LALStatus status;
  memset(&status,0,sizeof(status));
  INT4 runComplete=0;
  INT4 acceptanceCount = 0;
  REAL8 nullLikelihood;
  REAL8 trigSNR = 0.0;
  REAL8 *ladder = NULL;			//the ladder
  REAL8 *annealDecay = NULL;
  INT4 parameter=0;
  INT4 annealStartIter = 0;
  UINT4 hotChain = 0;                 // Affects proposal setup
  REAL8 temp = 1;
  REAL8 tempDelta = 0.0;
  MPI_Request MPIrequest;

  LALInferenceMPIswapAcceptance swapReturn;
  LALInferenceVariables *propStats = runState->proposalStats; // Print proposal acceptance rates to file
  LALInferenceProposalStatistics *propStat;
  UINT4 propTrack = 0;
  UINT4 (*parallelSwap)(LALInferenceRunState *, REAL8 *, INT4, FILE *);
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 Neff = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Neff");
  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");
  INT4 acl=(INT4)INFINITY, PTacl=(INT4)INFINITY;
  INT4 iEff=0;
  REAL8 ladderMin,ladderMax;
  REAL8 timestamp,timestamp_epoch=0.0;
  struct timeval tv;

  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  nChain = MPIsize;		//number of parallel chain

  LALInferenceMCMCRunPhase *runPhase_p = XLALCalloc(sizeof(LALInferenceMCMCRunPhase), 1);
  *runPhase_p = LALINFERENCE_ONLY_PT;
  LALInferenceAddVariable(runState->algorithmParams, "runPhase", &runPhase_p,  LALINFERENCE_MCMCrunphase_ptr_t, LALINFERENCE_PARAM_FIXED);

  //
  /* Command line flags (avoid repeated checks of runState->commandLine) */
  //
  UINT4 annealingOn = 0; // Chains will be annealed
  UINT4 annealingStarted = 0; // Flag indicating proper bookkeeping complete for annealing
  if (LALInferenceGetProcParamVal(runState->commandLine, "--anneal")) {
    annealingOn = 1;
    *runPhase_p = LALINFERENCE_TEMP_PT;

    /* Non-blocking receive to accept change in state.
     * Run phase changes are triggered by the coldest chain,
     * then propogated from the top of the ladder down to avoid
     * lockup from PT swaps */
    if (MPIrank != 0) {
        if (MPIrank == nChain-1)
            MPI_Irecv(runPhase_p, 1, MPI_INT, 0, RUN_PHASE_COM, MPI_COMM_WORLD, &MPIrequest);
        else
            MPI_Irecv(runPhase_p, 1, MPI_INT, MPIrank+1, RUN_PHASE_COM, MPI_COMM_WORLD, &MPIrequest);
    }
  }

  /* Clustered-KDE proposal updates */
  INT4 kde_update_start    = 200;  // rough number of effective samples to start KDE updates
  INT4 kde_update_interval = 0;    // proposal will be updated 5 times per decade, so this interval will change
  INT4 last_kde_update = 0;        // effective sample size at last KDE update


  UINT4 diffEvo = 1; // Differential evolution
  if (LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution") || LALInferenceGetProcParamVal(runState->commandLine, "--nodifferentialevolution")) {
    diffEvo = 0;
  }

  UINT4 adapting = 0;

  if (LALInferenceGetProcParamVal(runState->commandLine,"--varyFlow")) {
    /* Metropolis-coupled MCMC Swap (assumes likelihood function differs between chains).*/
    parallelSwap = &LALInferenceMCMCMCswap;
  } else {
    /* Standard parallel tempering swap. */
    parallelSwap = &LALInferencePTswap;
  }

  UINT4 tempVerbose = 0; // Print temperature swaps to file
  if (LALInferenceGetProcParamVal(runState->commandLine, "--tempVerbose")) {
    tempVerbose = 1;
  }

  UINT4 adaptVerbose = 0; // Print adaptation info to file
  if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")) {
    adaptVerbose = 1;
  }

  UINT4 benchmark = 0; // Print timestamps to chain outputs
  if (LALInferenceGetProcParamVal(runState->commandLine, "--benchmark")) {
    benchmark = 1;
  }

  ProcessParamsTable *ppt;

  ladder = malloc(nChain * sizeof(REAL8));                  // Array of temperatures for parallel tempering.
  annealDecay = malloc(nChain * sizeof(REAL8));           			// Used by annealing scheme

  /* If not specified otherwise, set effective sample size to total number of iterations */
  if (!Neff) {
    Neff = Niter;
    LALInferenceSetVariable(runState->algorithmParams, "Neff", &Neff);
  }

  /* Determine network SNR if injection was done */
  REAL8 networkSNRsqrd = 0.0;
  LALInferenceIFOData *IFO = runState->data;
  while (IFO != NULL) {
    networkSNRsqrd  += IFO->SNR * IFO->SNR;
    IFO = IFO->next;
  }

  /* Adaptation settings */
  LALInferenceSetupAdaptiveProposals(runState);
  INT4  adaptationOn = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptationOn")); // Run adapts
  INT4  adaptTau     = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptTau"));     // Sets decay of adaption function
  INT4  adaptLength       = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptLength"));// Number of iterations to adapt before turning off
  REAL8 s_gamma           = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "s_gamma"));                // Sets the size of changes to jump size during adaptation
  INT4  endOfPhase        = Neff;


  LALInferenceAddVariable(runState->algorithmParams, "acl", &acl,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);

  /* burnin settings */
  INT4 burnin       = 1;                           // Burnin phase where proposal ratio is ignored
  INT4 burninLength = (INT4)(0.25 * adaptLength);  // Number of iterations to turn off proposal ratio

  LALInferenceAddVariable(runState->proposalArgs, "burnin", &burnin,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "burninLength", &burninLength,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);

  /* Standard parallel tempering */
  ladderMin = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMin");   // Min temp in ladder
  ladderMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");   // Max temp in ladder

  REAL8 targetHotLike       = nPar/2.;          // Targeted max 'experienced' log(likelihood) of hottest chain
  INT4 hotThreshold         = nChain/2-1;        // If MPIrank > hotThreshold, use different proposal set

  /*  If running hot chains for Thermodynamic Integration make all draw from prior */
  if(ladderMin>0.0) hotThreshold=-1;

  /* Set maximum temperature (command line value take precidence) */
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempMax")) {
    if(MPIrank==0)
      fprintf(stdout,"Using tempMax specified by commandline: %f.\n", ladderMax);
  } else if (LALInferenceGetProcParamVal(runState->commandLine,"--trigSNR")) {        //--trigSNR given, choose ladderMax to get targetHotLike
    trigSNR = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "trigSNR");
    networkSNRsqrd = trigSNR * trigSNR;
    ladderMax = networkSNRsqrd/(2*targetHotLike);
    if(MPIrank==0)
      fprintf(stdout,"Trigger SNR of %f specified, setting tempMax to %f.\n", trigSNR, ladderMax);
  } else if (networkSNRsqrd > 0.0) {                                                  //injection, choose ladderMax to get targetHotLike
    ladderMax = networkSNRsqrd/(2*targetHotLike);
    if(MPIrank==0)
      fprintf(stdout,"Injecting SNR of %f, setting tempMax to %f.\n", sqrt(networkSNRsqrd), ladderMax);
  } else {                                                                            //If all else fails, use the default
    ladderMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");
    if(MPIrank==0)
      fprintf(stdout,"No --trigSNR or --tempMax specified, and not injecting a signal. Setting tempMax to default of %f.\n", ladderMax);
  }
  LALInferenceSetVariable(runState->algorithmParams, "tempMax", &ladderMax);

  if (ladderMin > ladderMax) {
    fprintf(stdout,"WARNING: tempMin > tempMax.  Forcing tempMin=1.0.\n");
    ladderMin = 1.0;
    LALInferenceSetVariable(runState->algorithmParams, "tempMin", &ladderMin);
  }


  /* Parallel tempering settings */
  INT4 Tskip            = 100;                                // Number of iterations between proposed temperature swaps 
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip"))
    Tskip = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip")->value);
  LALInferenceAddVariable(runState->algorithmParams, "tempSkip", &Tskip,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);

  INT4  annealStart     = 500;                                // # of autocorrelation lengths after adaptation before annealing
  INT4  annealLength    = 100;                                // # of autocorrelation lenghts to cool temperatures to ~1.0

  if (annealingOn) {
    if (LALInferenceGetProcParamVal(runState->commandLine,"--annealStart"))
      annealStart = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--annealStart")->value);

    if (LALInferenceGetProcParamVal(runState->commandLine,"--annealLength"))
      annealLength = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--annealLength")->value);

    endOfPhase=annealStart;
  }

  for (t=0; t<nChain; ++t) {
    ladder[t] = 0.0;
  }

  if (runState->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood ||
      runState->likelihood==&LALInferenceFreqDomainLogLikelihood){
    nullLikelihood = LALInferenceNullLogLikelihood(runState->data);
  } else if (runState->likelihood==&LALInferenceFreqDomainStudentTLogLikelihood || 
	     (runState->likelihood==&LALInferenceMarginalisedTimeLogLikelihood &&
        !LALInferenceGetProcParamVal(runState->commandLine, "--malmquistPrior")) ) {
    LALInferenceIFOData *headData = runState->data;
    REAL8 d = *(REAL8 *)LALInferenceGetVariable(runState->currentParams, "distance");
    REAL8 bigD = INFINITY;

    /* Don't store to cache, since distance scaling won't work */
    LALSimInspiralWaveformCache *cache = runState->model->waveformCache;
    runState->model->waveformCache = NULL;

    LALInferenceSetVariable(runState->currentParams, "distance", &bigD);
    nullLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->model);

    runState->model->waveformCache = cache;
    while (headData != NULL) {
      headData->nullloglikelihood = runState->model->loglikelihood;
      headData = headData->next;
    }

    LALInferenceSetVariable(runState->currentParams, "distance", &d);
  } else {
    nullLikelihood = 0.0;
    LALInferenceIFOData *headData = runState->data;
    while (headData != NULL) {
      headData->nullloglikelihood = 0.0;
      headData = headData->next;
    }
  }

  //null log likelihood logic doesn't work with noise parameters
  if(runState->likelihood==&LALInferenceMarginalisedTimeLogLikelihood &&
     (LALInferenceGetProcParamVal(runState->commandLine,"--psdFit") ||
      LALInferenceGetProcParamVal(runState->commandLine,"--glitchFit") ) )
    nullLikelihood = 0.0;

  LALInferenceAddVariable(runState->algorithmParams, "nChain", &nChain,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->algorithmParams, "nPar", &nPar,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs, "parameter",&parameter, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "nullLikelihood", &nullLikelihood, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);

  /* Construct temperature ladder */
  if(nChain > 1){
    if(LALInferenceGetProcParamVal(runState->commandLine, "--inverseLadder")) {     //Spacing uniform in 1/T
      tempDelta = (1.0/ladderMin - 1.0/ladderMax)/(REAL8)(nChain-1);
      for (t=0; t<nChain; ++t) {
        ladder[t]=1.0/(REAL8)(1.0/ladderMin-t*tempDelta);
      }
    } else {                                                                        //Geometric spacing
      if (LALInferenceGetProcParamVal(runState->commandLine, "--tempLadderBottomUp") && nPar != 1)
        tempDelta=1.+sqrt(2./(REAL8)nPar);
      else
        tempDelta=pow(ladderMax/ladderMin,1.0/(REAL8)(nChain-1));
      for (t=0;t<nChain; ++t) {
        ladder[t]=ladderMin*pow(tempDelta,t);
      }
    }
  } else {
    if(LALInferenceGetProcParamVal(runState->commandLine,"--tempMax")){   //assume --tempMax specified intentionally
      ladder[0]=ladderMax;
    }else{
      ladder[0]=1.0;
      ladderMax=1.0;
    }
  }
  temp=ladder[MPIrank];

  if (MPIrank > hotThreshold) {
    hotChain = 1;
  }

  LALInferenceAddVariable(runState->proposalArgs, "temperature", &temp,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "hotChain", &hotChain, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

  FILE * chainoutput = NULL;
  FILE *statfile = NULL;
  FILE *propstatfile = NULL;
  FILE *proptrackfile = NULL;
  FILE *swapfile = NULL;
  char statfilename[256];
  char propstatfilename[256];
  char proptrackfilename[256];
  char swapfilename[256];
  if (tempVerbose && MPIrank > 0) {
    sprintf(swapfilename,"PTMCMC.tempswaps.%u.%2.2d",randomseed,MPIrank);
    swapfile = fopen(swapfilename, "w");
    fprintf(swapfile, "cycle\tlog(chain_swap)\tlow_temp_likelihood\thigh_temp_likelihood\tswap_accepted\n");  // Print header for temp stat file
  }

  if (adaptationOn && adaptVerbose) {
    sprintf(statfilename,"PTMCMC.statistics.%u.%2.2d",randomseed,MPIrank);
    statfile = fopen(statfilename, "w");

    /* Print header for adaptation stats */
    fprintf(statfile,"cycle\ts_gamma");
    LALInferenceVariableItem *ptr=runState->currentParams->head;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        fprintf(statfile, "\tsigma_%s", LALInferenceTranslateInternalToExternalParamName(ptr->name));
      }
      ptr=ptr->next;
    }
    ptr=runState->currentParams->head;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        fprintf(statfile, "\tPaccept_%s", LALInferenceTranslateInternalToExternalParamName(ptr->name));
      }
      ptr=ptr->next;
    }
    fprintf(statfile,"\n");
  }

  if (propStats) {
    sprintf(propstatfilename,"PTMCMC.propstats.%u.%2.2d",randomseed,MPIrank);
    propstatfile = fopen(propstatfilename, "w");
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--propTrack")){
    propTrack=1;
    runState->preProposalParams = XLALCalloc(1, sizeof(LALInferenceVariableItem));
    runState->proposedParams = XLALCalloc(1, sizeof(LALInferenceVariableItem));

    UINT4 a=0;
    LALInferenceAddVariable(runState->proposalArgs, "accepted", &a, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
    sprintf(proptrackfilename,"PTMCMC.proptrack.%u.%2.2d",randomseed,MPIrank);
    proptrackfile = fopen(proptrackfilename, "w");
  }

  // initialize starting likelihood value:
  runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->model);
  runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);

  REAL8 logLAtAdaptStart = runState->currentLikelihood;
  LALInferenceSetVariable(runState->proposalArgs, "logLAtAdaptStart", &(logLAtAdaptStart));

  chainoutput = LALInferencePrintPTMCMCHeaderOrResume(runState);
  if (MPIrank == 0) {
    LALInferencePrintPTMCMCInjectionSample(runState);
  }

  if (MPIrank == 0){
    printf("\nTemperature ladder:\n");
    for (t=0; t<nChain; ++t) {
      printf(" ladder[%d]=%f\n",t,ladder[t]);
    }
  }

  /* Add parallel swaps to proposal tracking structure */
  if (propStats && MPIrank < nChain-1) {
    if(!LALInferenceCheckVariable(propStats, parallelSwapProposalName)) {
        LALInferenceProposalStatistics newPropStat = {
        .weight = 0,
        .proposed = 0,
        .accepted = 0};
      LALInferenceAddVariable(propStats, parallelSwapProposalName, (void *)&newPropStat, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_LINEAR);
    }
  }


  if (benchmark)
    timestamp_epoch = *((REAL8 *)LALInferenceGetVariable(runState->algorithmParams, "timestamp_epoch"));

  /* Print run details */
  if (MPIrank == 0) {
    printf("\nParallel Behavior:\n");
    if (adaptationOn)
      printf(" Adapting with decay power %i for %i iterations after max log(L) increases by nParams/2 (%1.2f).\n", adaptTau, adaptLength, (double)nPar/2.0);
    else
      printf(" Adaptation off.\n");
    if (annealingOn)
      printf(" Annealing linearly for %i effective samples.\n", annealLength);
    else
      printf(" Annealing off.\n");
    if (Neff != Niter)
      printf(" Collecting %i effective samples.\n", Neff);
  }


  if (MPIrank == 0) {
    printf("\nPTMCMCAlgorithm(); starting parameter values:\n");
    LALInferencePrintVariables(runState->currentParams);
    printf(" MCMC iteration: 0\t");
    printf("%f\t", runState->currentLikelihood - nullLikelihood);
    printf("\n");

    /* Print to file the contents of model->freqhPlus. */
    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--data-dump");
    if (ppt) {
      LALInferenceDataDump(runState->data, runState->model);
    }
  }

  /* Setup non-blocking recieve that will trigger from chain 0 and propigate down the ladder */
  if (MPIrank == nChain-1)
    MPI_Irecv(&runComplete, 1, MPI_INT, 0, RUN_COMPLETE, MPI_COMM_WORLD, &MPIrequest);
  else if (MPIrank != 0)
    MPI_Irecv(&runComplete, 1, MPI_INT, MPIrank+1, RUN_COMPLETE, MPI_COMM_WORLD, &MPIrequest);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  // iterate:
  i=0;
  while (!runComplete) {
    /* Increment iteration counter */
    i++;

    if (*runPhase_p == LALINFERENCE_LADDER_UPDATE)
      LALInferenceLadderUpdate(runState, 0, i);

    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));

    if (adaptationOn)
      LALInferenceAdaptation(runState, i);

    // Parallel tempering phase
    if (*runPhase_p == LALINFERENCE_ONLY_PT || *runPhase_p == LALINFERENCE_TEMP_PT) {

      //ACL calculation during parallel tempering
      if (i % (100*Nskip) == 0) {
        adapting = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));

        if (adapting)
          iEff = 0;
        else
          iEff = LALInferenceComputeEffectiveSampleSize(runState);
      }

      if (MPIrank==0 && iEff > endOfPhase) {
        if (*runPhase_p==LALINFERENCE_ONLY_PT) {
          fprintf(stdout,"Chain %i has %i effective samples. Stopping...\n", MPIrank, iEff);
          runComplete = 1;          // Sampling is done!
        } else if (*runPhase_p==LALINFERENCE_TEMP_PT) {
          /* Broadcast new run phase to other chains */
          *runPhase_p=LALINFERENCE_ANNEALING;
          MPI_Send(runPhase_p, 1, MPI_INT, nChain-1, RUN_PHASE_COM, MPI_COMM_WORLD);
        }
      }
    }

    // Annealing phase
    if (*runPhase_p==LALINFERENCE_ANNEALING) {
      if (!annealingStarted) {
        annealingStarted = 1;
        annealStartIter=i;

        /* Propogate the change in runPhase down the ladder, waiting for the next attempted PT swap */
        if (MPIrank > 1) {
            MPI_Probe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);  // Completes when PT swap attempted
            MPI_Send(runPhase_p, 1, MPI_INT, MPIrank-1, RUN_PHASE_COM, MPI_COMM_WORLD);  // Update runPhase
            LALInferenceFlushPTswap();                                // Rejects the swap
        }

        /* Broadcast the cold chain ACL from parallel tempering */
        PTacl = *((INT4*) LALInferenceGetVariable(runState->algorithmParams, "acl"));
        MPI_Bcast(&PTacl, 1, MPI_INT, 0, MPI_COMM_WORLD);

        /* Switch to a new set of jump proposals */
        runState->proposal = &LALInferencePostPTProposal;

        if(MPIrank==0)
          printf("Starting to anneal at iteration %i.\n",i);

        /* Share DE buffer from cold chain */
        if (diffEvo)
          BcastDifferentialEvolutionPoints(runState, 0);

        /* Build KDE proposal */
        if (!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-kde")) {
            if (MPIrank!=0)
              LALInferenceSetupClusteredKDEProposalFromDEBuffer(runState);
          last_kde_update = 0;
        }

        /* Force chains to re-adapt */
        if (adaptationOn)
          LALInferenceAdaptationRestart(runState, i);

        /* Calculate speed of annealing based on ACL */
        for (t=0; t<nChain; ++t)
          annealDecay[t] = (ladder[t]-1.0)/(REAL8)(annealLength*PTacl);

        /* Reset effective sample size and ACL */
        iEff=0;
        acl = PTacl;
        LALInferenceSetVariable(runState->algorithmParams, "acl", &acl);
      }

      if (i-annealStartIter < PTacl*annealLength) {
        for (t=0;t<nChain; ++t) {
          ladder[t] = ladder[t] - annealDecay[t];
          LALInferenceSetVariable(runState->proposalArgs, "temperature", &(ladder[MPIrank]));
        }
      } else {
        *runPhase_p = LALINFERENCE_SINGLE_CHAIN;
        if (MPIrank==0)
          printf(" Single-chain sampling starting at iteration %i.\n", i);
      }
    } //if (runState==ANNEALING)


    //Post-annealing single-chain sampling
    if (*runPhase_p==LALINFERENCE_SINGLE_CHAIN && i % (100*Nskip) == 0) {
      adapting = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));
      if (!adapting) {
        iEff = LALInferenceComputeEffectiveSampleSize(runState);
        if (iEff >= Neff/nChain) {
          fprintf(stdout,"Chain %i has %i effective samples. Stopping...\n", MPIrank, iEff);
          break;                                 // Sampling is done for this chain!
        }
      }
    } //if (*runPhase_p==SINGLE_CHAIN)

    INT4 accepted = runState->evolve(runState); //evolve the chain at temperature ladder[t]
    acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

    LALInferenceTrackProposalAcceptance(runState, accepted);

    /* Print proposal tracking headers now that the proposal cycle should be built. */
    if (i==1){
      if (propStats) {
          fprintf(propstatfile, "cycle\t");
          LALInferencePrintProposalStatsHeader(propstatfile, runState->proposalStats);
          fflush(propstatfile);
      }

      if (propTrack) {
          fprintf(proptrackfile, "cycle\t");
          LALInferencePrintProposalTrackingHeader(proptrackfile, runState->currentParams);
          fflush(proptrackfile);
      }
    }

    if ((i % Nskip) == 0) {

      /* Update clustered-KDE proposal every time the buffer is expanded */
      if (!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-kde")
          && (iEff > kde_update_start)
          && (((iEff - last_kde_update) > kde_update_interval) ||
              ((last_kde_update - iEff) > kde_update_interval))) {
        LALInferenceSetupClusteredKDEProposalFromDEBuffer(runState);

        /* Update 5 times each decade.  This keeps hot chains (with lower ACLs) under control */
        kde_update_interval = 2 * ((INT4) pow(10.0, floor(log10((REAL8) iEff))));

        /* Reset proposal counting */
        if (propStats && LALInferenceCheckVariable(propStats, clusteredKDEProposalName)) {
            propStat = (LALInferenceProposalStatistics *)LALInferenceGetVariable(propStats, clusteredKDEProposalName);
            propStat->proposed = 0;
            propStat->accepted = 0;
        }

        last_kde_update = iEff;
      }

      if (diffEvo) {
	    if (i % (runState->differentialPointsSkip) == 0)
	      accumulateDifferentialEvolutionSample(runState);
      }

      if (LALInferenceGetProcParamVal(runState->commandLine, "--kDTree") || LALInferenceGetProcParamVal(runState->commandLine, "--kdtree")) {
        accumulateKDTreeSample(runState);
      }

      fseek(chainoutput, 0L, SEEK_END);
      fprintf(chainoutput, "%d\t%f\t%f\t", i,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior,runState->currentPrior);
      LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
      fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);

      UINT4 ifo = 0;
      LALInferenceIFOData *headIFO = runState->data;
      while (headIFO != NULL) {
        fprintf(chainoutput, "%f\t", runState->currentIFOLikelihoods[ifo] - headIFO->nullloglikelihood);
        ifo++;
        headIFO = headIFO->next;
      }

      if (LALInferenceGetProcParamVal(runState->commandLine, "--output-SNRs")) {
          headIFO = runState->data;
          ifo = 0;
          while (headIFO != NULL) {
            fprintf(chainoutput, "%f\t", runState->currentIFOSNRs[ifo]);
            headIFO = headIFO->next;
            ifo++;
          }
          fprintf(chainoutput, "%f\t", runState->model->SNR);
      }

      if (benchmark) {
        gettimeofday(&tv, NULL);
        timestamp = tv.tv_sec + tv.tv_usec/1E6 - timestamp_epoch;
        fprintf(chainoutput, "%f\t", timestamp);
      }
      fprintf(chainoutput,"\n");
      fflush(chainoutput);

      if (adaptationOn == 1 && adaptVerbose) {
        fseek(statfile, 0L, SEEK_END);
        fprintf(statfile,"%d\t",i);
        if (LALInferenceCheckVariable(runState->proposalArgs, "s_gamma")) {
          s_gamma = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "s_gamma");
        } else {
          s_gamma = 0.0;
        }
        fprintf(statfile,"%f\t",s_gamma);
        
        for(LALInferenceVariableItem *item=runState->currentParams->head;item;item=item->next){
            char tmpname[MAX_STRLEN]="";
            sprintf(tmpname,"%s_%s",item->name,ADAPTSUFFIX);
            REAL8 *sigma=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
              fprintf(statfile,"%g\t",*sigma);
        }
        for(LALInferenceVariableItem *item=runState->currentParams->head;item;item=item->next){
            char tmpname[MAX_STRLEN]=""; 
            sprintf(tmpname,"%s_%s",item->name,ACCEPTSUFFIX);
            REAL8 *naccepted=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname); 
            sprintf(tmpname,"%s_%s",item->name,PROPOSEDSUFFIX);
            REAL8 *nproposed=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname); 
              fprintf(statfile,"%f\t",*naccepted/( *nproposed==0 ? 1.0 : *nproposed ));
        }
        fprintf(statfile,"\n");
        fflush(statfile);
      }

      if (propStats){
        fprintf(propstatfile, "%d\t", i);
        LALInferencePrintProposalStats(propstatfile, propStats);
        fflush(propstatfile);
      }

      if (propTrack) {
        REAL8 logProposalRatio = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logProposalRatio");
        fprintf(proptrackfile, "%d\t", i);
        LALInferencePrintProposalTracking(proptrackfile, runState->proposalArgs, runState->preProposalParams, runState->proposedParams, logProposalRatio, accepted);
      }
    }

    /* Excute swap proposal. */
    if (*runPhase_p == LALINFERENCE_ONLY_PT || *runPhase_p == LALINFERENCE_TEMP_PT) {
      swapReturn = parallelSwap(runState, ladder, i, swapfile); 
      if (propStats) {
        if (swapReturn != NO_SWAP_PROPOSED) {
          propStat = ((LALInferenceProposalStatistics *)LALInferenceGetVariable(propStats, parallelSwapProposalName));
          propStat->proposed++;
          if (swapReturn == ACCEPTED_SWAP)
            propStat->accepted++;
        }
      }
    }

    if (MPIrank==0 && i > Niter)
      runComplete=1;

  }// while (!runComplete)
  
  /* Send complete message to hottest chain, and it will propogate down to ensure no
   * hanging swap proposals happen */
  LALInferenceShutdownLadder();

  fclose(chainoutput);

  if(MPIrank == 0){
    if (adaptVerbose) {
      fclose(statfile);
    }
    if (tempVerbose) {
      fclose(swapfile);
    }
    if (propStats) {
      fclose(propstatfile);
    }
    if (propTrack) {
      fclose(proptrackfile);
    }
  }

  XLALFree(ladder);
  XLALFree(annealDecay);
}

INT4 PTMCMCOneStep(LALInferenceRunState *runState)
  // Metropolis-Hastings sampler.
{
  int MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  REAL8 logPriorCurrent, logPriorProposed;
  REAL8 logLikelihoodCurrent, logLikelihoodProposed;
  LALInferenceVariables proposedParams;
  REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
  REAL8 logAcceptanceProbability;
  REAL8 temperature;
  REAL8 targetAcceptance = 0.234;
  INT4 acceptanceCount;
  INT4 accepted = 0;

  // current values:
  logPriorCurrent      = runState->currentPrior;
  logLikelihoodCurrent = runState->currentLikelihood;

  temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");
  acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

  // generate proposal:
  proposedParams.head = NULL;
  proposedParams.dimension = 0;

  logProposalRatio = runState->proposal(runState, runState->currentParams, &proposedParams);

  // compute prior & likelihood:
  logPriorProposed = runState->prior(runState, &proposedParams, runState->model);
  if (logPriorProposed > -DBL_MAX)
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->model);
  else
    logLikelihoodProposed = -DBL_MAX;

  if (runState->preProposalParams != NULL)
    LALInferenceCopyVariables(runState->currentParams, runState->preProposalParams);

  if (runState->proposedParams != NULL)
    LALInferenceCopyVariables(&proposedParams, runState->proposedParams);

  //REAL8 nullLikelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
  //printf("%10.10f\t%10.10f\t%10.10f\n", logPriorProposed-logPriorCurrent, logLikelihoodProposed-nullLikelihood, logProposalRatio);
  //LALInferencePrintVariables(&proposedParams);
  
  // determine acceptance probability:
  logAcceptanceProbability = (1.0/temperature)*(logLikelihoodProposed - logLikelihoodCurrent)
    + (logPriorProposed - logPriorCurrent)
    + logProposalRatio;

  // accept/reject:
  if ((logAcceptanceProbability > 0)
      || (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
    LALInferenceCopyVariables(&proposedParams, runState->currentParams);
    runState->currentLikelihood = logLikelihoodProposed;
    runState->currentPrior = logPriorProposed;

    /* Calculate SNR if requested, and not already calculated by prior */
    LALInferenceIFOData *ifoPtr = runState->data;
    UINT4 ifo = 0;
    if (LALInferenceGetProcParamVal(runState->commandLine, "--output-SNRs")) {
      if (runState->model->SNR == 0.0)
        LALInferenceNetworkSNR(runState->currentParams, runState->data, runState->model);
      runState->currentSNR = runState->model->SNR;
      while (ifoPtr) {
        runState->currentIFOSNRs[ifo] = runState->model->ifo_SNRs[ifo];
        ifo++;
        ifoPtr = ifoPtr->next;
      }
    }

    ifoPtr = runState->data;
    ifo = 0;
    while (ifoPtr) {
      runState->currentIFOLikelihoods[ifo] = runState->model->ifo_loglikelihoods[ifo];
      ifo++;
      ifoPtr = ifoPtr->next;
    }

    acceptanceCount++;
    accepted = 1;
    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount);

  }

  if (LALInferenceCheckVariable(runState->proposalArgs, "accepted"))
    LALInferenceSetVariable(runState->proposalArgs, "accepted", &accepted);

  if (!LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio"))
    LALInferenceAddVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceSetVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);

  LALInferenceUpdateAdaptiveJumps(runState, accepted, targetAcceptance);
  LALInferenceClearVariables(&proposedParams);

  return accepted;
}


//-----------------------------------------
// Swap routines:
//-----------------------------------------
UINT4 LALInferencePTswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile)
{
  INT4 MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Status MPIstatus;
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  INT4 Tskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "tempSkip");
  REAL8 adjCurrentLikelihood, adjCurrentPrior;
  REAL8 logChainSwap;
  INT4 readyToSwap = 0;
  UINT4 swapProposed=0;
  INT4 swapAccepted=0;
  LALInferenceMPIswapAcceptance swapReturn;

  /* If Tskip reached, then block until next chain in ladder is prepared to accept swap proposal */
  if (((i % Tskip) == 0) && MPIrank < nChain-1) {
    swapProposed = 1;
    /* Send current likelihood for swap proposal */
    MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

    /* Determine if swap was accepted */
    MPI_Recv(&swapAccepted, 1, MPI_INT, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);

    /* Perform Swap */
    if (swapAccepted) {
      /* Set new likelihood */
      MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      runState->currentLikelihood = adjCurrentLikelihood;

      /* Exchange current prior values */
      MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
      MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
      runState->currentPrior = adjCurrentPrior;

      /* Package and send parameters */
      REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
      LALInferenceCopyVariablesToArray(runState->currentParams, parameters);
      MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

      /* Recieve and unpack parameters */
      REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
      MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      LALInferenceCopyArrayToVariables(adjParameters, runState->currentParams);

      XLALFree(parameters);
      XLALFree(adjParameters);
    }
  }

  /* Check if next lower temperature is ready to swap */
  else if (MPIrank > 0) {
    MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &readyToSwap, &MPIstatus);

    /* Hotter chain decides acceptance */
    if (readyToSwap) {
      /* Receive adjacent likelilhood */
      MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);

      /* Determine if swap is accepted and tell the other chain */
      logChainSwap = (1.0/ladder[MPIrank-1]-1.0/ladder[MPIrank])*(runState->currentLikelihood-adjCurrentLikelihood);
      if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) {
        swapAccepted = 1;
      } else {
        swapAccepted = 0;
      }

      MPI_Send(&swapAccepted, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);

      /* Print to file if verbose is chosen */
      if (swapfile != NULL) {
        fprintf(swapfile,"%d\t%f\t%f\t%f\t%i\n",i,logChainSwap,adjCurrentLikelihood,runState->currentLikelihood,swapAccepted);
        fflush(swapfile);
      }

      /* Perform Swap */
      if (swapAccepted) {
        /* Swap likelihoods */
        MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        runState->currentLikelihood=adjCurrentLikelihood;

        /* Exchange current prior values */
        MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        runState->currentPrior = adjCurrentPrior;

        /* Package parameters */
        REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
        LALInferenceCopyVariablesToArray(runState->currentParams, parameters);

        /* Swap parameters */
        REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
        MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);

        /* Unpack parameters */
        LALInferenceCopyArrayToVariables(adjParameters, runState->currentParams);

        XLALFree(parameters);
        XLALFree(adjParameters);
      }
    }
  }

  /* Return values for colder chain: 0=nothing happened; 1=swap proposed, not accepted; 2=swap proposed & accepted */
  if (swapProposed) {
      if (swapAccepted)
          swapReturn = ACCEPTED_SWAP;
      else
          swapReturn = REJECTED_SWAP;
  } else {
      swapReturn = NO_SWAP_PROPOSED;
  }
  return swapReturn;
}

/* Utility function that checks if the next coldest chain is attempting a PT swap.  If so, it is rejected. */
void LALInferenceFlushPTswap() {
    INT4 MPIrank;
    INT4 attemptingSwap=0;
    INT4 swapRejection = 0;
    REAL8 dummyLikelihood;
    MPI_Status MPIstatus;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    if (MPIrank==0) {
        return;
    } else {
        printf("Chain %i waiting for PT comm.\n", MPIrank);
        while (!attemptingSwap)
            MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &attemptingSwap, &MPIstatus);
        printf("Chain %i flushing swap proposed by %i.\n", MPIrank, MPIrank-1);
        MPI_Recv(&dummyLikelihood, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&swapRejection, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);
    }
    return;
}

/* Function to shut down the ladder */
void LALInferenceShutdownLadder() {
    INT4 MPIrank, MPIsize;
    INT4 attemptingSwap = 0;
    INT4 swapRejection = 0;
    INT4 runComplete = 1;
    REAL8 dummyLikelihood;
    MPI_Status MPIstatus;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

    if (MPIrank == 0)
        MPI_Send(&runComplete, 1, MPI_INT, MPIsize-1, RUN_COMPLETE, MPI_COMM_WORLD);
    else if (MPIrank != 1){
        printf("Chain %i waiting for PT comm.\n", MPIrank);
        while (!attemptingSwap)
            MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &attemptingSwap, &MPIstatus);
        MPI_Send(&runComplete, 1, MPI_INT, MPIrank-1, RUN_COMPLETE, MPI_COMM_WORLD);
        printf("Chain %i flushing swap proposed by %i.\n", MPIrank, MPIrank-1);
        MPI_Recv(&dummyLikelihood, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&swapRejection, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return;
}


// TODO: Fix after the deletion of "acknowledgedRunPhase"
void LALInferenceLadderUpdate(LALInferenceRunState *runState, INT4 sourceChainFlag, INT4 cycle)
{
  INT4 MPIrank, chain;
  INT4 readyToSend=0;
  MPI_Status MPIstatus;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  LALInferenceMCMCRunPhase *runPhase = *(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "runPhase");

  REAL8 *params = XLALMalloc(nPar * sizeof(REAL8));

  if (sourceChainFlag == 1) {
    /* Inform hotter chains of ladder update */
    *runPhase=LALINFERENCE_LADDER_UPDATE;
    for (chain=MPIrank+1; chain<nChain; chain++)
      MPI_Send(runPhase, 1, MPI_INT, chain, RUN_PHASE_COM, MPI_COMM_WORLD);

    /* Package and send current parameters */
    LALInferenceCopyVariablesToArray(runState->currentParams, params);

    /* Start from the top down since colder chains may still be getting PT swaps rejected */
    for (chain=nChain-1; chain>MPIrank; chain--) {
      MPI_Send(params, nPar, MPI_DOUBLE, chain, LADDER_UPDATE_COM, MPI_COMM_WORLD);
    }

  } else {
    /* Flush out any lingering swap proposals from colder chains */
    LALInferenceFlushPTswap();

    /* Wait until source chain is ready to send parameters */
    while (!readyToSend) {
      MPI_Iprobe(MPI_ANY_SOURCE, LADDER_UPDATE_COM, MPI_COMM_WORLD, &readyToSend, &MPIstatus);
    }

    /* Recieve new parameters and unpack into current params */
    MPI_Recv(params, nPar, MPI_DOUBLE, MPIstatus.MPI_SOURCE, LADDER_UPDATE_COM, MPI_COMM_WORLD, &MPIstatus);

    LALInferenceCopyArrayToVariables(params, runState->currentParams);

    /* Update prior and likelihood */
    runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->model);

    /* Restart adaptation to tune for the current location */
    LALInferenceAdaptationRestart(runState, cycle);
  }

  /* Reset runPhase to the last phase each chain was in */
  // TODO: Fix after the deletion of "acknowledgedRunPhase"

  XLALFree(params);
}

UINT4 LALInferenceMCMCMCswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile)
{
  INT4 MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Status MPIstatus;
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  INT4 Tskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "tempSkip");
  REAL8 adjCurrentPrior;
  REAL8 logChainSwap;
  INT4 readyToSwap = 0;
  UINT4 swapProposed=0;
  INT4 swapAccepted=0;

  REAL8 lowLikeLowParams = 0;
  REAL8 highLikeHighParams = 0;
  REAL8 lowLikeHighParams = 0;
  REAL8 highLikeLowParams = 0;

  LALInferenceParamVaryType fLowVary = LALINFERENCE_PARAM_FIXED;
  LALInferenceVariables *adjCurrentParams = NULL;
  adjCurrentParams = (LALInferenceVariables *)XLALCalloc(sizeof(LALInferenceVariables), 1);

  /* If Tskip reached, then block until next chain in ladder is prepared to accept swap proposal */
  if (((i % Tskip) == 0) && MPIrank < nChain-1) {
    swapProposed = 1;
    /* Send current likelihood for swap proposal */
    MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

    /* Fix fLow if it is being varied, to avoid swapping it */
    if(LALInferenceCheckVariable(runState->currentParams,"fLow")) {
      fLowVary = LALInferenceGetVariableVaryType(runState->currentParams,"fLow");
      if(fLowVary != LALINFERENCE_PARAM_FIXED) {
        LALInferenceSetParamVaryType(runState->currentParams,"fLow",LALINFERENCE_PARAM_FIXED);
        nPar -= 1;
      }
    }

    /* Prepare Variables structure to recieve adjacent parameters */
    LALInferenceCopyVariables(runState->currentParams, adjCurrentParams);

    /* Package, swap, and unpack parameters */
    REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
    REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
    LALInferenceCopyVariablesToArray(runState->currentParams, parameters);
    MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
    MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
    LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
    XLALFree(parameters);
    XLALFree(adjParameters);

    /* Calculate likelihood at adjacent parameters and send */
    lowLikeHighParams = runState->likelihood(adjCurrentParams, runState->data, runState->model);
    MPI_Send(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

    /* Determine if swap was accepted */
    MPI_Recv(&swapAccepted, 1, MPI_INT, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);

    if (swapAccepted) {
      /* Set new likelihood */
      runState->currentLikelihood = lowLikeHighParams;

      MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
      MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      runState->currentPrior = adjCurrentPrior;

      /* Set new parameters */
      LALInferenceCopyVariables(adjCurrentParams,runState->currentParams);
    }

    /* Unfix fLow if it was originally unfixed */
    if(fLowVary!=LALINFERENCE_PARAM_FIXED) {
      LALInferenceSetParamVaryType(runState->currentParams,"fLow",fLowVary);
      nPar += 1;
    }

    if (swapAccepted) {
      /* Calculate prior at new values */
      runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
    }
  }

  /* Check if next lower temperature is ready to swap */
  else if (MPIrank > 0) {
    MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &readyToSwap, &MPIstatus);

    /* Hotter chain decides acceptance */
    if (readyToSwap) {
      /* Receive adjacent likelilhood */
      MPI_Recv(&lowLikeLowParams, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      highLikeHighParams = runState->currentLikelihood;

      /* Fix fLow if it is being varied, to avoid swapping it */
      if(LALInferenceCheckVariable(runState->currentParams,"fLow")) {
        fLowVary = LALInferenceGetVariableVaryType(runState->currentParams,"fLow");
        if(fLowVary != LALINFERENCE_PARAM_FIXED) {
          LALInferenceSetParamVaryType(runState->currentParams,"fLow",LALINFERENCE_PARAM_FIXED);
          nPar -= 1;
        }
      }

      /* Prepare Variables structure to recieve adjacent parameters */
      LALInferenceCopyVariables(runState->currentParams, adjCurrentParams);

      /* Package, swap, and unpack parameters */
      REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
      REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
      LALInferenceCopyVariablesToArray(runState->currentParams, parameters);
      MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
      LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
      XLALFree(parameters);
      XLALFree(adjParameters);

      /* Calculate likelihood at adjacent parameters */
      highLikeLowParams = runState->likelihood(adjCurrentParams, runState->data, runState->model);

      /* Recieve likelihood from adjacent chain */
      MPI_Recv(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);

      /* Propose swap */
      logChainSwap = (1./ladder[MPIrank-1])*(lowLikeHighParams-lowLikeLowParams)+(1./ladder[MPIrank])*(highLikeLowParams-highLikeHighParams);

      if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) {
        swapAccepted = 1;
      } else {
        swapAccepted = 0;
      }
      MPI_Send(&swapAccepted, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);

      /* Print to file if verbose is chosen */
      if (swapfile != NULL) {
        fprintf(swapfile,"%d\t%f\t%f\t%f\t%i\n",i,logChainSwap,lowLikeLowParams,highLikeHighParams,swapAccepted);
        fflush(swapfile);
      }

      if (swapAccepted) {
        /* Swap likelihoods */
        runState->currentLikelihood = highLikeLowParams;

        /* Exchange current prior values (assumes prior functions are identical) */
        MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        runState->currentPrior = adjCurrentPrior;

        /* Set new parameters */
        LALInferenceCopyVariables(adjCurrentParams,runState->currentParams);
      }

      /* Unfix fLow if it was originally unfixed */
      if(fLowVary!=LALINFERENCE_PARAM_FIXED) {
        LALInferenceSetParamVaryType(runState->currentParams,"fLow",fLowVary);
        nPar += 1;
      }

      if (swapAccepted) {
        /* Calculate prior at new values */
        runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
      }
    }
  }

  /* Return values for colder chain: 0=nothing happened; 1=swap proposed, not accepted; 2=swap proposed & accepted */
  if (swapProposed && swapAccepted)
    swapProposed++;

  LALInferenceClearVariables(adjCurrentParams);
  XLALFree(adjCurrentParams);
  return swapProposed;
}


//-----------------------------------------
// Adaptation:
//-----------------------------------------
void LALInferenceAdaptation(LALInferenceRunState *runState, INT4 cycle)
{
  INT4 MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 adapting = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adapting");
  INT4 burnin = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "burnin");
  INT4 burninLength = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "burninLength");
  INT4 adaptStart = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart");
  INT4 adaptLength = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptLength");
  REAL8 logLAtAdaptStart = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logLAtAdaptStart");

  /* Only burnin at the beginning of the run */
  if (burnin && (cycle > burninLength)) {
      burnin = 0;
      LALInferenceSetVariable(runState->proposalArgs, "burnin", &burnin);
  }

  /* if maximum logL has increased by more than nParam/2, restart it */
  if (runState->currentLikelihood > logLAtAdaptStart+(REAL8)nPar/2) {
    if (!adapting)
      fprintf(stdout,"Turning on adaptation for chain %u at iteration %u.\n",MPIrank,cycle);
    LALInferenceAdaptationRestart(runState, cycle);
  } else if (adapting) {
    /* Turn off adaption after adaptLength steps without restarting */
    if ((cycle-adaptStart) > adaptLength) {
      adapting = 0;  //turn off adaptation
      LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
      LALInferenceRemoveVariable(runState->proposalArgs,"s_gamma");

      /* Clear differential buffer so that it contains only post-burnin samples */
      resetDifferentialEvolutionBuffer(runState);

      fprintf(stdout,"Ending adaptation for chain %u at iteration %u.\n",MPIrank,cycle);

    /* Else set adaptation envelope */
    } else {
      LALInferenceAdaptationEnvelope(runState, cycle);
    }
  }
}


//-----------------------------------------
// Restart adaptation:
//-----------------------------------------
void LALInferenceAdaptationRestart(LALInferenceRunState *runState, INT4 cycle)
{
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 adapting=1;

  INT4 MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  for(LALInferenceVariableItem *item=runState->currentParams->head;item;item=item->next){
    if (item->vary != LALINFERENCE_PARAM_FIXED && item->vary != LALINFERENCE_PARAM_OUTPUT) {
      char tmpname[MAX_STRLEN]="";

      sprintf(tmpname,"%s_%s",item->name,ACCEPTSUFFIX);
      REAL8 *accepted=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
      *accepted = 0;

      sprintf(tmpname,"%s_%s",item->name,PROPOSEDSUFFIX);
      REAL8 *proposed=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
      *proposed = 0;
    }
  }

  LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
  LALInferenceSetVariable(runState->proposalArgs, "adaptStart", &cycle);
  LALInferenceSetVariable(runState->proposalArgs, "logLAtAdaptStart", &(runState->currentLikelihood));
  LALInferenceSetVariable(runState->algorithmParams, "acl", &Niter);
  LALInferenceAdaptationEnvelope(runState, cycle);
}

//-----------------------------------------
// Adaptation envelope function:
//-----------------------------------------
void LALInferenceAdaptationEnvelope(LALInferenceRunState *runState, INT4 cycle)
{
  INT4 adaptStart = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart");
  INT4 adaptTau = *(INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptTau");
  INT4 adaptLength = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptLength");
  INT4 adaptResetBuffer = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptResetBuffer");
  REAL8 s_gamma = 0.0;

  if (cycle-adaptStart <= adaptResetBuffer) {
    s_gamma=(((REAL8)cycle-(REAL8)adaptStart)/(REAL8)adaptResetBuffer)*(((REAL8)cycle-(REAL8)adaptStart)/(REAL8)(adaptResetBuffer));
  } else if (cycle-adaptStart < adaptLength) {
    s_gamma=10.0*exp(-(1.0/adaptTau)*log((REAL8)(cycle-adaptStart)))-1;
  } else {
    s_gamma=0.0;
  }

  if (LALInferenceCheckVariable(runState->proposalArgs, "s_gamma"))
    LALInferenceSetVariable(runState->proposalArgs, "s_gamma", &s_gamma);
  else
    LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
}


//-----------------------------------------
// file output routines:
//-----------------------------------------
FILE *LALInferencePrintPTMCMCHeaderOrResume(LALInferenceRunState *runState) {
  ProcessParamsTable *ppt;
  char *outFileName = NULL;
  FILE *chainoutput = NULL;
  UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");
  int MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
  if (ppt) {
    outFileName = (char*)XLALCalloc(strlen(ppt->value)+255,sizeof(char*));
    sprintf(outFileName,"%s.%2.2d",ppt->value,MPIrank);
  } else {
    outFileName = (char*)XLALCalloc(255,sizeof(char*));
    sprintf(outFileName,"PTMCMC.output.%u.%2.2d",randomseed,MPIrank);
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--resume") && access(outFileName, R_OK) == 0) {
    /* Then file already exists for reading, and we're going to resume
       from it, so don't write the header. */

    chainoutput = fopen(outFileName, "r");
    if (chainoutput == NULL) {
      XLALErrorHandler = XLALExitErrorHandler;
      XLALPrintError("Error reading resume file (in %s, line %d)\n", __FILE__, __LINE__);
      XLAL_ERROR_NULL(XLAL_EIO);
    }
    
    LALInferenceMCMCResumeRead(runState, chainoutput);

    fclose(chainoutput);
  } else {
    chainoutput = fopen(outFileName,"w");
    if(chainoutput == NULL){
      XLALErrorHandler = XLALExitErrorHandler;
      XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
      XLAL_ERROR_NULL(XLAL_EIO);
    }
    
    LALInferencePrintPTMCMCHeaderFile(runState, chainoutput);

    fclose(chainoutput);
  }

  chainoutput = fopen(outFileName, "a");
  if (chainoutput == NULL) {
    XLALErrorHandler = XLALExitErrorHandler;
    XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
    XLAL_ERROR_NULL(XLAL_EIO);
  }
  
  XLALFree(outFileName);

  return chainoutput;
}

void LALInferencePrintPTMCMCHeaderFile(LALInferenceRunState *runState, FILE *chainoutput)
{
  int MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  REAL8 nullLikelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  REAL8 temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");
  REAL8 SampleRate=4096.0; //default value of the sample rate from LALInferenceReadData()
  UINT4 nIFO=0;
  LALInferenceIFOData *ifodata1=runState->data;
  REAL8 timestamp;
  struct timeval tv;

  REAL8 fRef= 0.0;
  if(LALInferenceCheckVariable(runState->currentParams,"fRef")) fRef = *(REAL8*)LALInferenceGetVariable(runState->currentParams, "fRef");

  while(ifodata1){
    nIFO++;
    ifodata1=ifodata1->next;
  }

  int waveform = 0;
  if(LALInferenceCheckVariable(runState->currentParams,"LAL_APPROXIMANT")) waveform= *(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_APPROXIMANT");
  double pnorder = 0.0;
  if(LALInferenceCheckVariable(runState->currentParams,"LAL_PNORDER")) pnorder = ((double)(*(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_PNORDER")))/2.0;

  char *str;
  str = LALInferencePrintCommandLine(runState->commandLine);

  REAL8 networkSNR=0.0;
  ifodata1=runState->data;
  while(ifodata1){
    networkSNR+=ifodata1->SNR*ifodata1->SNR;
    ifodata1=ifodata1->next;
  }
  networkSNR=sqrt(networkSNR);

  if(LALInferenceGetProcParamVal(runState->commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(runState->commandLine,"--srate")->value);

  UINT4 benchmark=0;
  if(LALInferenceGetProcParamVal(runState->commandLine,"--benchmark")) benchmark=1;

    fprintf(chainoutput, "  LALInference version:%s,%s,%s,%s,%s\n", lalAppsVCSId,lalAppsVCSDate,lalAppsVCSBranch,lalAppsVCSAuthor,lalAppsVCSStatus);
    fprintf(chainoutput,"  %s\n",str);
    fprintf(chainoutput, "%10s  %10s  %6s  %20s  %6s %8s   %6s  %10s  %12s  %9s  %9s  %8s %8s\n",
        "nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tchain","Network SNR","Waveform","pN order","Npar","fRef");
    fprintf(chainoutput, "%10d  %10d  %u  %20.10lf  %6d %8d   %6d%12.1f%14.6f  %9i  %9.1f  %8i %12.1f\n",
        Niter,0,randomseed,nullLikelihood,nIFO,0,nChain,temperature,networkSNR,waveform,(double)pnorder,nPar,fRef);
    fprintf(chainoutput, "\n%16s  %16s  %10s  %10s  %10s  %10s  %20s  %15s  %12s  %12s  %12s\n",
        "Detector","SNR","f_low","f_high","before tc","after tc","Sample start (GPS)","Sample length","Sample rate","Sample size","FT size");
    ifodata1=runState->data;
    while(ifodata1){
      fprintf(chainoutput, "%16s  %16.8lf  %10.2lf  %10.2lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %.1f  %12d  %12d\n",
          ifodata1->detector->frDetector.name,ifodata1->SNR,ifodata1->fLow,ifodata1->fHigh,atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)-2.0,2.00,
          XLALGPSGetREAL8(&(ifodata1->epoch)),atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value),SampleRate,
          (int)(atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)*SampleRate),
          (int)(atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value)*SampleRate));
      ifodata1=ifodata1->next;
    }
    fprintf(chainoutput, "\n\n%31s\n","");
    fprintf(chainoutput, "cycle\tlogpost\tlogprior\t");
    LALInferenceFprintParameterNonFixedHeaders(chainoutput, runState->currentParams);
    fprintf(chainoutput, "logl\t");
    LALInferenceIFOData *headIFO = runState->data;
    while (headIFO != NULL) {
      fprintf(chainoutput, "logl");
      fprintf(chainoutput, "%s",headIFO->name);
      fprintf(chainoutput, "\t");
      headIFO = headIFO->next;
    }
    if (LALInferenceGetProcParamVal(runState->commandLine, "--output-SNRs")) {
        headIFO = runState->data;
        while (headIFO != NULL) {
          fprintf(chainoutput, "SNR");
          fprintf(chainoutput, "%s",headIFO->name);
          fprintf(chainoutput, "\t");
          headIFO = headIFO->next;
        }
        fprintf(chainoutput, "SNR\t");
    }

    if (benchmark)
      fprintf(chainoutput, "timestamp\t");
    fprintf(chainoutput,"\n");
    fprintf(chainoutput, "%d\t%f\t%f\t", 0, (runState->currentLikelihood - nullLikelihood)+runState->currentPrior, runState->currentPrior);
    LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
    fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);
    headIFO = runState->data;
    while (headIFO != NULL) {
      fprintf(chainoutput, "%f\t", runState->currentLikelihood - headIFO->nullloglikelihood);
      headIFO = headIFO->next;
    }
    if(benchmark) {
      gettimeofday(&tv, NULL);
      timestamp = tv.tv_sec + tv.tv_usec/1E6;
      LALInferenceAddVariable(runState->algorithmParams, "timestamp_epoch", &timestamp,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
      fprintf(chainoutput, "%f\t", 0.0);
    }
    fprintf(chainoutput,"\n");
}

void LALInferencePrintPTMCMCInjectionSample(LALInferenceRunState *runState) {
  ProcessParamsTable *ppt;

  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--inj");
  if (ppt) {
    ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
    UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");
    FILE *out = NULL;
    char *fname = NULL;
    LALInferenceVariables *saveParams = NULL;

    saveParams = (LALInferenceVariables *)XLALCalloc(sizeof(LALInferenceVariables), 1);

    if (ppt2) {
      fname = (char *) XLALCalloc((strlen(ppt2->value)+255)*sizeof(char), 1);
      sprintf(fname, "%s.injection", ppt2->value);
    } else {
      fname = (char *) XLALCalloc(255*sizeof(char), 1);
      sprintf(fname, "PTMCMC.output.%u.injection", randomseed);
    }
    out = fopen(fname, "w");

    LALInferenceCopyVariables(runState->currentParams, saveParams);

    SimInspiralTable *injTable = NULL;
    SimInspiralTable *theEventTable = NULL;

    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
    
    ppt2 = LALInferenceGetProcParamVal(runState->commandLine, "--event");
    if (ppt2) {
      UINT4 event = atoi(ppt2->value);
      UINT4 i;
      theEventTable = injTable;
      for (i = 0; i < event; i++) {
        theEventTable = theEventTable->next;
      }
      theEventTable->next = NULL;
    } else {
      theEventTable=injTable;
      theEventTable->next = NULL;
    }

    REAL8 m1 = theEventTable->mass1;
    REAL8 m2 = theEventTable->mass2;
    REAL8 q = m2/m1;
    REAL8 eta = m1*m2/(m1+m2)/(m1+m2);

    if (q > 1.0) q = 1.0/q;

    REAL8 sx = theEventTable->spin1x;
    REAL8 sy = theEventTable->spin1y;
    REAL8 sz = theEventTable->spin1z;

    REAL8 a_spin1 = sqrt(sx*sx + sy*sy + sz*sz);
    
    REAL8 theta_spin1, phi_spin1;
    if (a_spin1 == 0.0) {
      theta_spin1 = 0.0;
      phi_spin1 = 0.0;
    } else {
      theta_spin1 = acos(sz / a_spin1);
      phi_spin1 = atan2(sy, sx);
      if (phi_spin1 < 0.0) phi_spin1 += 2.0*M_PI;
    }

    sx = theEventTable->spin2x;
    sy = theEventTable->spin2y;
    sz = theEventTable->spin2z;
    
    REAL8 a_spin2 = sqrt(sx*sx + sy*sy + sz*sz), theta_spin2, phi_spin2;
    if (a_spin2 == 0.0) {
      theta_spin2 = 0.0;
      phi_spin2 = 0.0;
    } else {
      theta_spin2 = acos(sz / a_spin2);
      phi_spin2 = atan2(sy, sx);
      if (phi_spin2 < 0.0) phi_spin2 += 2.0*M_PI;
    }

    REAL8 psi = theEventTable->polarization;
    if (psi>=M_PI) psi -= M_PI;

    REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->geocent_end_time));

    REAL8 chirpmass = theEventTable->mchirp;

    REAL8 dist = theEventTable->distance;
    REAL8 inclination = theEventTable->inclination;
    REAL8 phase = theEventTable->coa_phase;
    REAL8 dec = theEventTable->latitude;
    REAL8 ra = theEventTable->longitude;

    LALInferenceSetVariable(runState->currentParams, "chirpmass", &chirpmass);
    if (LALInferenceCheckVariable(runState->currentParams, "asym_massratio")) {
      LALInferenceSetVariable(runState->currentParams, "asym_massratio", &q);
    } else if (LALInferenceCheckVariable(runState->currentParams, "massratio")) {
      LALInferenceSetVariable(runState->currentParams, "massratio", &eta);
    } else {
      /* Restore state, cleanup, and throw error */
      LALInferenceCopyVariables(saveParams, runState->currentParams);
      XLALFree(fname);
      LALInferenceClearVariables(saveParams);
      XLALFree(saveParams);
      XLAL_ERROR_VOID(XLAL_EINVAL, "unknown mass ratio parameter name (allowed are 'massratio' or 'asym_massratio')");
    }

    UINT4 added_time_param = 0;
    if (!LALInferenceCheckVariable(runState->currentParams, "time")) {
        added_time_param = 1;
        LALInferenceAddVariable(runState->currentParams, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    } else {
        LALInferenceSetVariable(runState->currentParams, "time", &injGPSTime);
    }

    UINT4 added_phase_param = 0;
    if (!LALInferenceCheckVariable(runState->currentParams, "phase")) {
      added_phase_param = 1;
      LALInferenceAddVariable(runState->currentParams, "phase", &phase, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    } else {
      LALInferenceSetVariable(runState->currentParams, "phase", &phase);
    }

    LALInferenceSetVariable(runState->currentParams, "distance", &dist);
    LALInferenceSetVariable(runState->currentParams, "theta_JN", &inclination);
    LALInferenceSetVariable(runState->currentParams, "polarisation", &(psi));
    LALInferenceSetVariable(runState->currentParams, "declination", &dec);
    LALInferenceSetVariable(runState->currentParams, "rightascension", &ra);
    if (LALInferenceCheckVariable(runState->currentParams, "a_spin1")) {
      LALInferenceSetVariable(runState->currentParams, "a_spin1", &a_spin1);
    }
    if (LALInferenceCheckVariable(runState->currentParams, "theta_spin1")) {
      LALInferenceSetVariable(runState->currentParams, "theta_spin1", &theta_spin1);
    }
    if (LALInferenceCheckVariable(runState->currentParams, "phi_spin1")) {
      LALInferenceSetVariable(runState->currentParams, "phi_spin1", &phi_spin1);
    }
    if (LALInferenceCheckVariable(runState->currentParams, "a_spin2")) {
      LALInferenceSetVariable(runState->currentParams, "a_spin2", &a_spin2);
    }
    if (LALInferenceCheckVariable(runState->currentParams, "theta_spin2")) {
      LALInferenceSetVariable(runState->currentParams, "theta_spin2", &theta_spin2);
    }
    if (LALInferenceCheckVariable(runState->currentParams, "phi_spin2")) {
      LALInferenceSetVariable(runState->currentParams, "phi_spin2", &phi_spin2);
    }

    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->model);
    runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
    LALInferencePrintPTMCMCHeaderFile(runState, out);
    fclose(out);

    if (added_time_param) {
        LALInferenceRemoveVariable(runState->currentParams, "time");
        LALInferenceRemoveMinMaxPrior(runState->priorArgs, "time");
    }

    if (added_phase_param) {
      LALInferenceRemoveVariable(runState->currentParams, "phase");
      LALInferenceRemoveMinMaxPrior(runState->priorArgs, "phase");
    }

    LALInferenceCopyVariables(saveParams, runState->currentParams);
    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->model);
    runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);

    XLALFree(fname);
    LALInferenceClearVariables(saveParams);
    XLALFree(saveParams);
  }
}

void LALInferenceDataDump(LALInferenceIFOData *data, LALInferenceModel *model){

  const UINT4 nameLength=256;
  char filename[nameLength];
  FILE *out;
  UINT4 ui;

  snprintf(filename, nameLength, "freqTemplatehPlus.dat");
  out = fopen(filename, "w");
  for (ui = 0; ui < model->freqhPlus->data->length; ui++) {
    REAL8 f = model->freqhPlus->deltaF * ui;
    COMPLEX16 d = model->freqhPlus->data->data[ui];

    fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
  }
  fclose(out);

  snprintf(filename, nameLength, "freqTemplatehCross.dat");
  out = fopen(filename, "w");
  for (ui = 0; ui < model->freqhCross->data->length; ui++) {
    REAL8 f = model->freqhCross->deltaF * ui;
    COMPLEX16 d = model->freqhCross->data->data[ui];

    fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
  }
  fclose(out);

  while (data != NULL) {
    snprintf(filename, nameLength, "%s-freqTemplateStrain.dat", data->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < model->freqhCross->data->length; ui++) {
      REAL8 f = model->freqhCross->deltaF * ui;
      COMPLEX16 d;
      d = data->fPlus * model->freqhPlus->data->data[ui] +
             data->fCross * model->freqhCross->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d) );
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplateStrain.dat", data->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < model->timehCross->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(model->timehCross->epoch)) +
        data->timeshift + ui*model->timehCross->deltaT;
      REAL8 d = data->fPlus*model->timehPlus->data->data[ui] +
        data->fCross*model->timehCross->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplatehPlus.dat", data->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < model->timehPlus->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(model->timehPlus->epoch)) +
        ui * model->timehPlus->deltaT;
      REAL8 d = model->timehPlus->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplatehCross.dat", data->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < model->timehCross->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(model->timehCross->epoch)) +
        ui * model->timehCross->deltaT;
      REAL8 d = model->timehCross->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    data = data->next;
  }
}

void LALInferenceMCMCResumeRead(LALInferenceRunState *runState, FILE *resumeFile) {
  /* Hope that the line is shorter than 16K! */
  const long len = 16384;
  char linebuf[len];
  char *last_line = NULL;
  long flen, line_length;
  int cycle;
  float loglike, logprior;  

  fseek(resumeFile, 0L, SEEK_END);
  flen = ftell(resumeFile);

  if (flen < len) {
    fseek(resumeFile, 0L, SEEK_SET);
    fread(linebuf, flen, 1, resumeFile);
    linebuf[flen-1] = '\0'; /* Strip off trailing newline. */
  } else {
    fseek(resumeFile, -len, SEEK_END);
    fread(linebuf, len, 1, resumeFile);
    linebuf[len-1] = '\0'; /* Strip off trailing newline.... */
  }

  last_line = strrchr(linebuf, '\n'); /* Last occurence of '\n' */
  last_line += 1;

  line_length = strlen(last_line);

  /* Go to the beginning of the last line. */
  fseek(resumeFile, -line_length, SEEK_END);

  fscanf(resumeFile, "%d %f %f", &cycle, &loglike, &logprior);

  LALInferenceReadSampleNonFixed(resumeFile, runState->currentParams);
}
