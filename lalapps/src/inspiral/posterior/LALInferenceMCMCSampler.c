/*
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
 *
 *  Copyright (C) 2009, 2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch and Will M. Farr
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
  size_t i;

  for (i = 0; i < runState->differentialPointsLength; i++) {
    LALInferenceClearVariables(runState->differentialPoints[i]);
    XLALFree(runState->differentialPoints[i]);
    runState->differentialPoints[i] = NULL;
  }

  runState->differentialPoints = XLALRealloc(runState->differentialPoints, 1*sizeof(LALInferenceVariables *));
  runState->differentialPointsLength = 0;
  runState->differentialPointsSize = 1;
  runState->differentialPointsSkip = 1;
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
  INT4 i=0;
  REAL8** packedDEsamples;
  REAL8*  temp;

  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  INT4 nPoints = runState->differentialPointsLength;

  /* Prepare a DE buffer of the proper size */
  MPI_Bcast(&nPoints, 1, MPI_INT, sourceTemp, MPI_COMM_WORLD);
  INT4 startCycle=0, endCycle=nPoints*Nskip;

  /* Prepare 2D array for DE points */
  packedDEsamples = (REAL8**) XLALMalloc(nPoints * sizeof(REAL8*));
  temp = (REAL8*) XLALMalloc(nPoints * nPar * sizeof(REAL8));
  for (i=0; i < nPoints; i++) {
    packedDEsamples[i] = temp + (i*nPar);
  }

  /* Pack it up */
  if (MPIrank == sourceTemp)
    LALInferenceBufferToArray(runState, startCycle, endCycle, packedDEsamples);

  /* Send it out */
  MPI_Bcast(packedDEsamples[0], nPoints*nPar, MPI_DOUBLE, sourceTemp, MPI_COMM_WORLD);

  /* Unpack it */
  if (MPIrank != sourceTemp)
    LALInferenceArrayToBuffer(runState, packedDEsamples);

  /* Clean up */
  XLALFree(temp);
  MPI_Barrier(MPI_COMM_WORLD);
}

static void
computeMaxAutoCorrLen(LALInferenceRunState *runState, INT4 startCycle, INT4 endCycle, INT4* maxACL) {
/* Define ACL as the smallest s such that
 *
 * 1 + 2*ACF(1) + 2*ACF(2) + ... + 2*ACF(M*s) < s,
 *
 * the short length so that the sum of the ACF function
 * is smaller than that length over a window of M times
 * that length.
 *
 * The maximum window length is restricted to be N/K as
 * a safety precaution against relying on data near the
 * extreme of the lags in the ACF, where there is a lot
 * of noise.
*/
  INT4 M=5, K=2;
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  INT4 totalPoints = runState->differentialPointsLength;
  INT4 start = (INT4)ceil((REAL8)startCycle/(REAL8)Nskip);
  INT4 end = (INT4)floor((REAL8)endCycle/(REAL8)Nskip);
  /* Include last point */
  if (end > totalPoints-1)
    end = totalPoints-1;
  INT4 nPoints = end - start + 1;
  REAL8** DEarray;
  REAL8*  temp;
  REAL8 mean, ACL, ACF, max=0;
  INT4 par=0, lag=0, i=0, imax;
  REAL8 cumACF, s;

  if (nPoints > 1) {
    imax = nPoints/K;
    /* Prepare 2D array for DE points */
    DEarray = (REAL8**) XLALMalloc(nPoints * sizeof(REAL8*));
    temp = (REAL8*) XLALMalloc(nPoints * nPar * sizeof(REAL8));
    for (i=0; i < nPoints; i++) {
      DEarray[i] = temp + (i*nPar);
    }

    LALInferenceBufferToArray(runState, startCycle, endCycle, DEarray);

    for (par=0; par<nPar; par++) {
      mean = gsl_stats_mean(&DEarray[0][par], nPar, nPoints);
      for (i=0; i<nPoints; i++)
        DEarray[i][par] -= mean;

      lag=1;
      ACL=1.0;
      ACF=1.0;
      s=1.0;
      cumACF=1.0;
      while (cumACF >= s) {
        ACF = gsl_stats_correlation(&DEarray[0][par], nPar, &DEarray[lag][par], nPar, nPoints-lag);
        cumACF += 2.0 * ACF;
        lag++;
        s = (REAL8)lag/(REAL8)M;
        if (lag > imax) {
          ACL=(REAL8)Niter/(REAL8)Nskip;
          break;
        }
      }
      ACL = s*Nskip;
      if (ACL>max)
        max=ACL;
    }
    XLALFree(temp);
    //XLALFree(DEarray);
  } else {
    max = Niter;
  }

  /* Account for any thinning of the DE buffer that has happend */
  max *= runState->differentialPointsSkip;

  *maxACL = (INT4)max;
}

static void
updateMaxAutoCorrLen(LALInferenceRunState *runState, INT4 currentCycle) {
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 proposedACL=0;
  INT4 adaptStart = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart");
  INT4 adaptLength = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptLength");
  INT4 iEffStart = adaptStart+adaptLength;
  // Calculate ACL with latter half of data to avoid ACL overestimation from chain equilibrating after adaptation
  INT4 aclStart = (currentCycle+iEffStart)/2;
  INT4 acl=Niter;
  INT4 goodACL=0;

  if (aclStart<currentCycle)
    computeMaxAutoCorrLen(runState, aclStart, currentCycle, &proposedACL);

  if (proposedACL < Niter && proposedACL > 0) {
    goodACL=1;
    acl = proposedACL;
  }
  LALInferenceSetVariable(runState->algorithmParams, "goodACL", &goodACL);

  LALInferenceSetVariable(runState->algorithmParams, "acl", &acl);
}


void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  INT4 i,t,c; //indexes for for() loops
  INT4 nChain;
  INT4 MPIrank, MPIsize;
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
  INT4 iEffStart = 0;
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
  INT4 acl=Niter, PTacl=Niter, goodACL=0;
  INT4 numACLchecks=10;             // Number of times to recalculate ACL
  INT4 aclCheckCounter=1;
  INT4 iEff=0;
  REAL8 ladderMin,ladderMax;
  REAL8 timestamp,timestamp_epoch=0.0;
  struct timeval tv;

  LALInferenceMCMCRunPhase *runPhase_p = XLALCalloc(sizeof(LALInferenceMCMCRunPhase), 1);
  *runPhase_p = LALINFERENCE_ONLY_PT;
  LALInferenceAddVariable(runState->algorithmParams, "runPhase", &runPhase_p,  LALINFERENCE_MCMCrunphase_ptr_t, LALINFERENCE_PARAM_FIXED);

  //
  /* Command line flags (avoid repeated checks of runState->commandLine) */
  //
  UINT4 annealingOn = 0; // Chains will be annealed
  if (LALInferenceGetProcParamVal(runState->commandLine, "--anneal")) {
    annealingOn = 1;
    *runPhase_p = LALINFERENCE_TEMP_PT;
  }

  /* Setup non-blocking recieve that will allow other chains to update the runPhase */
  acknowledgePhase(runState);

  LALInferenceMCMCRunPhase *phaseAcknowledged = *(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "acknowledgedRunPhase");

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

  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  nChain = MPIsize;		//number of parallel chain
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
  INT4  adaptStart        = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptStart"));                  // Keeps track of last iteration adaptation was restarted
  INT4  endOfPhase        = Neff;


  LALInferenceAddVariable(runState->algorithmParams, "acl", &acl,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->algorithmParams, "goodACL", &goodACL,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);


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
    LALSimInspiralWaveformCache *cache = headData->waveformCache;
    while (headData != NULL) {
      headData->waveformCache = NULL;
      headData = headData->next;
    }
    headData = runState->data;

    LALInferenceSetVariable(runState->currentParams, "distance", &bigD);
    nullLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->templt);

    while (headData != NULL) {
      headData->nullloglikelihood = headData->loglikelihood;
      headData->waveformCache = cache;
      headData = headData->next;
    }

    LALInferenceSetVariable(runState->currentParams, "distance", &d);
  } else {
    nullLikelihood = 0.0;
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
  runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->templt);
  LALInferenceIFOData *headData = runState->data;
  while (headData != NULL) {
    headData->acceptedloglikelihood = headData->loglikelihood;
    headData = headData->next;
  }
  runState->currentPrior = runState->prior(runState, runState->currentParams);

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

    /* Print to file the contents of ->freqModelhPlus->data->length. */
    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--data-dump");
    if (ppt) {
      LALInferenceDataDump(runState);
    }
  }

  /* Setup non-blocking recieve that will succeed when chain 0 is complete. */
  if (MPIrank!=0)
    MPI_Irecv(&runComplete, 1, MPI_INT, 0, RUN_COMPLETE, MPI_COMM_WORLD, &MPIrequest);

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
      if (i % (100*Nskip) == 0 && MPIrank == 0) {
        adapting = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));

        /* Check if cold chain ACL has been calculated */
        if (!adapting) {
          goodACL = *((INT4*) LALInferenceGetVariable(runState->algorithmParams, "goodACL"));
          if (!goodACL)
            updateMaxAutoCorrLen(runState, i);
          acl = *((INT4*) LALInferenceGetVariable(runState->algorithmParams, "acl"));

          if (goodACL) {
            adaptStart = *((INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart"));
            iEffStart = adaptStart+adaptLength;
            iEff = (i - iEffStart)/acl;

            /* Periodically recalculate ACL */
            if (*runPhase_p != LALINFERENCE_SINGLE_CHAIN && iEff > floor((REAL8)(aclCheckCounter*endOfPhase)/(REAL8)numACLchecks)) {
              updateMaxAutoCorrLen(runState, i);
              acl = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "acl");
              iEff = (i - iEffStart)/acl;
              aclCheckCounter = ceil((REAL8)iEff / ((REAL8)endOfPhase/(REAL8)numACLchecks));
            }
          } //if (goodACL)
        } //if (!adapting)
      } //if (i % (100*Nskip) == 0)

      if (MPIrank==0) {
        if (aclCheckCounter > numACLchecks) {
          if (*runPhase_p==LALINFERENCE_ONLY_PT) {
            fprintf(stdout,"Chain %i has %i effective samples. Stopping...\n", MPIrank, iEff);
            runComplete = 1;          // Sampling is done!
          } else if (*runPhase_p==LALINFERENCE_TEMP_PT) {
            /* Broadcast new run phase to other chains */
            *runPhase_p=LALINFERENCE_ANNEALING;
            for (c=1; c<nChain; c++) {
              MPI_Send(&*runPhase_p, 1, MPI_INT, c, RUN_PHASE_COM, MPI_COMM_WORLD);
            }
          }
        }
      }
    }

    // Annealing phase
    if (*runPhase_p==LALINFERENCE_ANNEALING) {
      if (*phaseAcknowledged!=LALINFERENCE_ANNEALING) {
        acknowledgePhase(runState);

        /* Broadcast the cold chain ACL from parallel tempering */
        PTacl = acl;
        MPI_Bcast(&PTacl, 1, MPI_INT, 0, MPI_COMM_WORLD);
        annealStartIter=i;
        runState->proposal = &LALInferencePostPTProposal;
        if(MPIrank==0)
          printf("Starting to anneal at iteration %i.\n",i);

        /* Share DE buffer from cold chain */
        if (diffEvo)
          BcastDifferentialEvolutionPoints(runState, 0);

        /* Force chains to re-adapt */
        if (adaptationOn)
          LALInferenceAdaptationRestart(runState, i);

        /* Calculate speed of annealing based on ACL */
        for (t=0; t<nChain; ++t) {
          annealDecay[t] = (ladder[t]-1.0)/(REAL8)(annealLength*PTacl);
        }

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
        *phaseAcknowledged=*runPhase_p;
        if (MPIrank==0)
          printf(" Single-chain sampling starting at iteration %i.\n", i);
      }
    } //if (runState==ANNEALING)


    //Post-annealing single-chain sampling
    if (*runPhase_p==LALINFERENCE_SINGLE_CHAIN) {
      adapting = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));
      adaptStart = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart");
      iEffStart = adaptStart+adaptLength;
      if (!adapting) {
        iEff = (i-iEffStart)/acl;
        if (iEff >= Neff/nChain) {
          /* Double check ACL before ending */
          updateMaxAutoCorrLen(runState, i);
          acl = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "acl");
          iEff = (i-iEffStart)/acl;
          if (iEff >= Neff/nChain) {
            fprintf(stdout,"Chain %i has %i effective samples. Stopping...\n", MPIrank, iEff);
            break;                                 // Sampling is done for this chain!
          }
        }
      }
    } //if (*runPhase_p==SINGLE_CHAIN)

    runState->evolve(runState); //evolve the chain at temperature ladder[t]
    acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

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
      if (diffEvo) {
	if (i % (Nskip*runState->differentialPointsSkip) == 0) 
	  accumulateDifferentialEvolutionSample(runState);
      }

      if (LALInferenceGetProcParamVal(runState->commandLine, "--kDTree") || LALInferenceGetProcParamVal(runState->commandLine, "--kdtree")) {
        accumulateKDTreeSample(runState);
      }

      fseek(chainoutput, 0L, SEEK_END);
      fprintf(chainoutput, "%d\t%f\t%f\t", i,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior,runState->currentPrior);
      LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
      fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);

      LALInferenceIFOData *headIFO = runState->data;
      while (headIFO != NULL) {
        fprintf(chainoutput, "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
        headIFO = headIFO->next;
      }

      REAL8 networkSNR = 0.0;
      headIFO = runState->data;
      while (headIFO != NULL) {
        fprintf(chainoutput, "%f\t", headIFO->acceptedSNR);
        networkSNR += headIFO->acceptedSNR * headIFO->acceptedSNR;
        headIFO = headIFO->next;
      }
      networkSNR = sqrt(networkSNR);
      fprintf(chainoutput, "%f\t", networkSNR);

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
            REAL8 *accepted=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname); 
            sprintf(tmpname,"%s_%s",item->name,PROPOSEDSUFFIX);
            REAL8 *proposed=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname); 
              fprintf(statfile,"%f\t",*accepted/( *proposed==0 ? 1.0 : *proposed ));
        }
        fprintf(statfile,"\n");
        fflush(statfile);
      }

      if (propStats){
        fprintf(propstatfile, "%d\t", i);
        LALInferencePrintProposalStats(propstatfile,runState->proposalStats);
        fflush(propstatfile);
      }

      if (propTrack) {
        fprintf(proptrackfile, "%d\t", i);
        LALInferencePrintProposalTracking(proptrackfile, runState->proposalArgs, runState->preProposalParams, runState->proposedParams);
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

    if (MPIrank==0 && runComplete==1) {
      for (c=1; c<nChain; c++) {
        MPI_Send(&runComplete, 1, MPI_INT, c, RUN_COMPLETE, MPI_COMM_WORLD);
      }
    }
  }// while (!runComplete)
  
  /* Flush any remaining PT swap attempts before moving on */
  LALInferenceFlushPTswap();
  MPI_Barrier(MPI_COMM_WORLD);

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

void acknowledgePhase(LALInferenceRunState *runState)
/* Handle MPI communications and runState params for a change in run phase */
{
  INT4 MPIrank;
  MPI_Request MPIrequest;

  LALInferenceMCMCRunPhase *runPhase = *(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "runPhase");
  if (!LALInferenceCheckVariable(runState->algorithmParams, "acknowledgedRunPhase")) {
    LALInferenceMCMCRunPhase acknowledgedRunPhase = *runPhase;
    LALInferenceMCMCRunPhase *acknowledgedRunPhase_p = &acknowledgedRunPhase;
    LALInferenceAddVariable(runState->algorithmParams, "acknowledgedRunPhase", &acknowledgedRunPhase_p,  LALINFERENCE_MCMCrunphase_ptr_t, LALINFERENCE_PARAM_FIXED);
  } else {
    LALInferenceMCMCRunPhase *acknowledgedRunPhase = *(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "acknowledgedRunPhase");
    *acknowledgedRunPhase = *runPhase;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  if (MPIrank!=0)
    MPI_Irecv(runPhase, 1, MPI_INT, MPI_ANY_SOURCE, RUN_PHASE_COM, MPI_COMM_WORLD, &MPIrequest);
}

void PTMCMCOneStep(LALInferenceRunState *runState)
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

  runState->proposal(runState, &proposedParams);
  if (LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio"))
    logProposalRatio = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logProposalRatio");

  // compute prior & likelihood:
  logPriorProposed = runState->prior(runState, &proposedParams);
  if (logPriorProposed > -DBL_MAX)
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->templt);
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

  //printf("%g=(1/%g)*(%g-%g)+(%g-%g)+%g\n",logAcceptanceProbability,temperature,logLikelihoodProposed,logLikelihoodCurrent,logPriorProposed,logPriorCurrent,logProposalRatio);

  // accept/reject:
  if ((logAcceptanceProbability > 0)
      || (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
    LALInferenceCopyVariables(&proposedParams, runState->currentParams);
    gsl_matrix_memcpy (runState->data->glitch_x, runState->data->glitch_y);
    runState->currentLikelihood = logLikelihoodProposed;
    LALInferenceIFOData *headData = runState->data;
    while (headData != NULL) {
      headData->acceptedloglikelihood = headData->loglikelihood;
      headData->acceptedSNR = headData->currentSNR;
      headData = headData->next;
    }
    runState->currentPrior = logPriorProposed;
    acceptanceCount++;
    accepted = 1;
    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount);

  }

  if (LALInferenceCheckVariable(runState->proposalArgs, "accepted"))
    LALInferenceSetVariable(runState->proposalArgs, "accepted", &accepted);

  LALInferenceUpdateAdaptiveJumps(runState, accepted, targetAcceptance);
  LALInferenceClearVariables(&proposedParams);
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

  REAL8Vector * parameters = NULL;
  REAL8Vector * adjParameters = NULL;

  /* If Tskip reached, then block until next chain in ladder is prepared to accept swap proposal */
  if (((i % Tskip) == 0) && MPIrank < nChain-1) {
    swapProposed = 1;
    /* Send current likelihood for swap proposal */
    MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

    /* Determine if swap was accepted */
    MPI_Recv(&swapAccepted, 1, MPI_INT, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);

    /* Perform Swap */
    if (swapAccepted) {
      adjParameters = XLALCreateREAL8Vector(nPar);

      /* Set new likelihood */
      MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      runState->currentLikelihood = adjCurrentLikelihood;

      /* Exchange current prior values */
      MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
      MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
      runState->currentPrior = adjCurrentPrior;

      /* Package and send parameters */
      parameters = LALInferenceCopyVariablesToArray(runState->currentParams);
      MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);

      /* Recieve and unpack parameters */
      MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      LALInferenceCopyArrayToVariables(adjParameters, runState->currentParams);

      XLALDestroyREAL8Vector(parameters);
      XLALDestroyREAL8Vector(adjParameters);
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
        adjParameters = XLALCreateREAL8Vector(nPar);

        /* Swap likelihoods */
        MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        runState->currentLikelihood=adjCurrentLikelihood;

        /* Exchange current prior values */
        MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        runState->currentPrior = adjCurrentPrior;

        /* Package parameters */
        parameters = LALInferenceCopyVariablesToArray(runState->currentParams);

        /* Swap parameters */
        MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);

        /* Unpack parameters */
        LALInferenceCopyArrayToVariables(adjParameters, runState->currentParams);

        /* Recompute glitch model */
        UINT4 glitchFlag = 0;
        if(LALInferenceCheckVariable(runState->currentParams,"glitchFitFlag"))
          glitchFlag = *((UINT4 *)LALInferenceGetVariable(runState->currentParams, "glitchFitFlag"));

        if(glitchFlag)
        {
          UINT4Vector *gsize = *(UINT4Vector **) LALInferenceGetVariable(runState->currentParams, "glitch_size");
          UINT4 ifo = 0;
          UINT4 n   = 0;

          /* Remove wavlet form linear combination */
          for(ifo=0; ifo<gsize->length; ifo++) for(n=0; n<runState->data->glitch_x->size2; n++) gsl_matrix_set(runState->data->glitch_x, ifo, n, 0.0);
          for(ifo=0; ifo<gsize->length; ifo++) for(n=0; n<gsize->data[ifo]; n++) UpdateWaveletSum(runState, runState->currentParams, runState->data->glitch_x, ifo, n, 1);
        }

        XLALDestroyREAL8Vector(parameters);
        XLALDestroyREAL8Vector(adjParameters);
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
    INT4 swapRejection=0;
    REAL8 dummyLikelihood;
    MPI_Status MPIstatus;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    if (MPIrank==0) {
        return;
    } else {
        MPI_Send(&swapRejection, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &attemptingSwap, &MPIstatus);
        if (attemptingSwap) {
          MPI_Recv(&dummyLikelihood, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
          //MPI_Send(&swapRejection, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);
        }
    }
    return;
}


void LALInferenceLadderUpdate(LALInferenceRunState *runState, INT4 sourceChainFlag, INT4 cycle)
{
  INT4 MPIrank, chain;
  INT4 readyToSend=0;
  MPI_Status MPIstatus;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  LALInferenceMCMCRunPhase *runPhase = *(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "runPhase");

  REAL8Vector *params = NULL;

  if (sourceChainFlag == 1) {
    /* Inform hotter chains of ladder update */
    *runPhase=LALINFERENCE_LADDER_UPDATE;
    for (chain=MPIrank+1; chain<nChain; chain++)
      MPI_Send(runPhase, 1, MPI_INT, chain, RUN_PHASE_COM, MPI_COMM_WORLD);

    /* Package and send current parameters */
    params = LALInferenceCopyVariablesToArray(runState->currentParams);

    /* Start from the top down since colder chains may still be getting PT swaps rejected */
    for (chain=nChain-1; chain>MPIrank; chain--) {
      MPI_Send(params->data, nPar, MPI_DOUBLE, chain, LADDER_UPDATE_COM, MPI_COMM_WORLD);
    }

  } else {
    /* Flush out any lingering swap proposals from colder chains */
    LALInferenceFlushPTswap();

    params = XLALCreateREAL8Vector(nPar);

    /* Wait until source chain is ready to send parameters */
    while (!readyToSend) {
      MPI_Iprobe(MPI_ANY_SOURCE, LADDER_UPDATE_COM, MPI_COMM_WORLD, &readyToSend, &MPIstatus);
    }

    /* Recieve new parameters and unpack into current params */
    MPI_Recv(params->data, nPar, MPI_DOUBLE, MPIstatus.MPI_SOURCE, LADDER_UPDATE_COM, MPI_COMM_WORLD, &MPIstatus);

    LALInferenceCopyArrayToVariables(params, runState->currentParams);

    /* Update prior and likelihood */
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->templt);

    /* Restart adaptation to tune for the current location */
    LALInferenceAdaptationRestart(runState, cycle);
  }

  /* Reset runPhase to the last phase each chain was in */
  LALInferenceMCMCRunPhase acknowledgedRunPhase = **(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "acknowledgedRunPhase");
  *runPhase = acknowledgedRunPhase;
  acknowledgePhase(runState);

  XLALDestroyREAL8Vector(params);
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
  REAL8Vector * parameters = NULL;
  REAL8Vector * adjParameters = NULL;
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
    adjParameters = XLALCreateREAL8Vector(nPar);
    parameters = LALInferenceCopyVariablesToArray(runState->currentParams);
    MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
    MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
    LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
    XLALDestroyREAL8Vector(parameters);
    XLALDestroyREAL8Vector(adjParameters);

    /* Calculate likelihood at adjacent parameters and send */
    lowLikeHighParams = runState->likelihood(adjCurrentParams, runState->data, runState->templt);
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
      runState->currentPrior = runState->prior(runState, runState->currentParams);
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
      adjParameters = XLALCreateREAL8Vector(nPar);
      parameters = LALInferenceCopyVariablesToArray(runState->currentParams);
      MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
      MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
      LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
      XLALDestroyREAL8Vector(parameters);
      XLALDestroyREAL8Vector(adjParameters);

      /* Calculate likelihood at adjacent parameters */
      highLikeLowParams = runState->likelihood(adjCurrentParams, runState->data, runState->templt);

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
        runState->currentPrior = runState->prior(runState, runState->currentParams);
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
  INT4 adaptStart = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptStart");
  INT4 adaptLength = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adaptLength");
  REAL8 logLAtAdaptStart = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logLAtAdaptStart");

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
  //LALInferenceMCMCRunPhase runPhase = **(LALInferenceMCMCRunPhase **) LALInferenceGetVariable(runState->algorithmParams, "runPhase");
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  //INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  INT4 adapting=1;
  INT4 goodACL=0;

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

  /* Also clear differential buffer when adaptation restarts. */
  resetDifferentialEvolutionBuffer(runState);

  LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
  LALInferenceSetVariable(runState->proposalArgs, "adaptStart", &cycle);
  LALInferenceSetVariable(runState->proposalArgs, "logLAtAdaptStart", &(runState->currentLikelihood));
  LALInferenceSetVariable(runState->algorithmParams, "acl", &Niter);
  LALInferenceSetVariable(runState->algorithmParams, "goodACL", &goodACL);
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

    fprintf(chainoutput, "  LALInference version:%s,%s,%s,%s,%s\n", LALAPPS_VCS_ID,LALAPPS_VCS_DATE,LALAPPS_VCS_BRANCH,LALAPPS_VCS_AUTHOR,LALAPPS_VCS_STATUS);
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
    headIFO = runState->data;
    while (headIFO != NULL) {
      fprintf(chainoutput, "SNR");
      fprintf(chainoutput, "%s",headIFO->name);
      fprintf(chainoutput, "\t");
      headIFO = headIFO->next;
    }
    fprintf(chainoutput, "SNR\t");

    if (benchmark)
      fprintf(chainoutput, "timestamp\t");
    fprintf(chainoutput,"\n");
    fprintf(chainoutput, "%d\t%f\t%f\t", 0, (runState->currentLikelihood - nullLikelihood)+runState->currentPrior, runState->currentPrior);
    LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
    fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);
    headIFO = runState->data;
    while (headIFO != NULL) {
      fprintf(chainoutput, "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
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

static void setIFOAcceptedLikelihoods(LALInferenceRunState *runState) {
  LALInferenceIFOData *data = runState->data;
  LALInferenceIFOData *ifo = NULL;

  for (ifo = data; ifo != NULL; ifo = ifo->next) {
    ifo->acceptedloglikelihood = ifo->loglikelihood;
    ifo->acceptedSNR = ifo->currentSNR;
  }
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

    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->templt);
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    setIFOAcceptedLikelihoods(runState);
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
    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->templt);
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    setIFOAcceptedLikelihoods(runState);    

    XLALFree(fname);
    LALInferenceClearVariables(saveParams);
    XLALFree(saveParams);
  }
}

void LALInferenceDataDump(LALInferenceRunState *runState){

  const UINT4 nameLength=256;
  char filename[nameLength];
  FILE *out;
  LALInferenceIFOData *headData = runState->data;
  UINT4 ui;

  while (headData != NULL) {

    snprintf(filename, nameLength, "%s-freqTemplatehPlus.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->freqModelhPlus->data->length; ui++) {
      REAL8 f = headData->freqModelhPlus->deltaF * ui;
      COMPLEX16 d = headData->freqModelhPlus->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-freqTemplatehCross.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->freqModelhCross->data->length; ui++) {
      REAL8 f = headData->freqModelhCross->deltaF * ui;
      COMPLEX16 d = headData->freqModelhCross->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-freqTemplateStrain.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->freqModelhCross->data->length; ui++) {
      REAL8 f = headData->freqModelhCross->deltaF * ui;
      COMPLEX16 d;
      d = headData->fPlus * headData->freqModelhPlus->data->data[ui] +
             headData->fCross * headData->freqModelhCross->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d) );
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplatehPlus.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->timeModelhPlus->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhPlus->epoch)) +
        ui * headData->timeModelhPlus->deltaT;
      REAL8 d = headData->timeModelhPlus->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplatehCross.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->timeModelhCross->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhCross->epoch)) +
        ui * headData->timeModelhCross->deltaT;
      REAL8 d = headData->timeModelhCross->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeTemplateStrain.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->timeModelhCross->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhCross->epoch)) +
        headData->timeshift + ui*headData->timeModelhCross->deltaT;
      REAL8 d = headData->fPlus*headData->timeModelhPlus->data->data[ui] +
        headData->fCross*headData->timeModelhCross->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    headData = headData->next;
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
