/*
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
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

#include <LALAppsVCSInfo.h>
#include <lal/LALStdlib.h>

#define PROGRAM_NAME "LALInferenceMCMCSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

void LALInferencePTswap(LALInferenceRunState *runState, double *TcurrentLikelihood, REAL8 *parametersVec,
                        REAL8 *tempLadder, int lowerRank, int upperRank, int i, FILE *tempfile);
FILE* LALInferencePrintPTMCMCHeader(LALInferenceRunState *runState);
void LALInferenceDataDump(LALInferenceRunState *runState);

static void
accumulateDifferentialEvolutionSample(LALInferenceRunState *runState) {
  if (runState->differentialPointsSize == runState->differentialPointsLength) {
    size_t newSize = runState->differentialPointsSize*2;
    runState->differentialPoints = XLALRealloc(runState->differentialPoints, newSize*sizeof(LALInferenceVariables *));
    runState->differentialPointsSize = newSize;
  }

  runState->differentialPoints[runState->differentialPointsLength] = XLALCalloc(1, sizeof(LALInferenceVariables));
  LALInferenceCopyVariables(runState->currentParams, runState->differentialPoints[runState->differentialPointsLength]);
  runState->differentialPointsLength += 1;
}

static void
accumulateKDTreeSample(LALInferenceRunState *runState) {
  LALInferenceVariables *proposalParams = runState->proposalArgs;

  if (!LALInferenceCheckVariable(proposalParams, "kDTree") || !LALInferenceCheckVariable(proposalParams, "kDTreeVariableTemplate")) {
    /* Improper setup---bail! */
    return;
  }

  LALInferenceKDTree *tree = *(LALInferenceKDTree **)LALInferenceGetVariable(proposalParams, "kDTree");
  LALInferenceVariables *template = *(LALInferenceVariables **)LALInferenceGetVariable(proposalParams, "kDTreeVariableTemplate");
  size_t ndim = LALInferenceGetVariableDimensionNonFixed(template);
  REAL8 *pt = XLALMalloc(ndim*sizeof(REAL8));

  LALInferenceKDVariablesToREAL8(runState->currentParams, pt, template);

  LALInferenceKDAddPoint(tree, pt);

  XLALFree(pt);
}

static void
BcastDifferentialEvolutionPoints(LALInferenceRunState *runState, int sourceTemp) {
  int MPIrank;
  double** packedDEpoints;
  double*  temp;
  LALInferenceVariableItem *ptr;
  int i=0,p=0;

  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nPoints = runState->differentialPointsLength;

  /* Prepare 2D array for DE points */
  packedDEpoints = (double**) XLALMalloc(nPoints * sizeof(double*));
  temp = (double*) XLALMalloc(nPoints * nPar * sizeof(double));
  for (i=0; i < nPoints; i++) {
    packedDEpoints[i] = temp + (i*nPar);
  }

  /* Pack it up */
  if (MPIrank==sourceTemp) {
    for (i=0; i < nPoints; i++) {
      ptr=runState->differentialPoints[i]->head;
      p=0;
      while(ptr!=NULL) {
        if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
          packedDEpoints[i][p]=*(double *)ptr->value;
          p++;
        }
        ptr=ptr->next;
      }
    }
  }

  /* Send it out */
  MPI_Bcast(packedDEpoints[0], nPoints*nPar, MPI_DOUBLE, sourceTemp, MPI_COMM_WORLD);

  /* Unpack it */
  if (MPIrank != sourceTemp) {
    for (i=0; i < nPoints; i++) {
      ptr=runState->differentialPoints[i]->head;
      p=0;
      while(ptr!=NULL) {
        if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
          *((REAL8 *)ptr->value) = (REAL8)packedDEpoints[i][p];
          p++;
        }
        ptr=ptr->next;
      }
    }
  }
  /* Clean up */
  XLALFree(temp);
  MPI_Barrier(MPI_COMM_WORLD);
}

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  int i,t,p,lowerRank,upperRank; //indexes for for() loops
  int nChain;
  int MPIrank, MPIsize;
  LALStatus status;
  memset(&status,0,sizeof(status));
  INT4 acceptanceCount = 0;
  INT4 swapAttempt=0;
  REAL8 nullLikelihood;
  REAL8 *tempLadder = NULL;			//the temperature ladder
  REAL8 *annealDecay = NULL;
  INT4 *acceptanceCountLadder = NULL;	//array of acceptance counts to compute the acceptance ratios.
  double *TcurrentLikelihood = NULL; //the current likelihood for each chain
  INT4 parameter=0;
  UINT4 hotChain = 0;                 // Affects proposal setup
  REAL8Vector *sigmas = NULL;
  REAL8Vector *PacceptCount = NULL;
  REAL8Vector *PproposeCount = NULL;
  REAL8 *parametersVec = NULL;
  REAL8 tempDelta = 0.0;
  REAL8Vector * parameters = NULL;

  INT4 annealingOn = 0;
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");

  ProcessParamsTable *ppt;

  MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  nChain = MPIsize;		//number of parallel chain
  tempLadder = malloc(nChain * sizeof(REAL8));                  // Array of temperatures for parallel tempering.
  acceptanceCountLadder = (int*) malloc(sizeof(int)*nChain);		// Array of acceptance counts to compute the acceptance ratios.
  annealDecay = malloc(nChain * sizeof(REAL8));           			// Used by annealing scheme


  if(MPIrank == 0){
    parametersVec = (REAL8 *)malloc(MPIsize*nPar*sizeof(REAL8));
    for (p=0;p<(nChain*nPar);++p){
      parametersVec[p] = 0.0;
    }
  }

  parameters = XLALCreateREAL8Vector(nPar);

  LALInferenceVariableItem *ptr=runState->currentParams->head;
  p=0;
  while(ptr!=NULL) {
    if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
      parameters->data[p]=*(REAL8 *)ptr->value;
      p++;
    }
    ptr=ptr->next;
  }

  /* Determine network SNR if injection was done */
  REAL8 networkSNRsqrd = 0.0;
  LALInferenceIFOData *IFO = runState->data;
  while (IFO != NULL) {
    networkSNRsqrd  += IFO->SNR * IFO->SNR;
    IFO = IFO->next;
  }

  /* Adaptation settings */
  INT4  adaptationOn = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptationOn")); // Run adapts
  INT4  adapting     = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));     // Current step being adapted
  REAL8 s_gamma           = 1.0;                // Sets the size of changes to jump size during adaptation
  INT4  adaptStart        = 0;                  // Keeps track of last iteration adaptation was restarted
  INT4  adaptResetBuffer  = 100;                // Number of iterations before adapting after a restart
  INT4  adaptTau          = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptTau")); // Sets the length and slope of adaption function
  INT4  adaptationLength  = pow(10,adaptTau);   // Number of iterations to adapt before turning off

  /* Temperature ladder settings */
  REAL8 tempMin = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMin");   // Min temperature in ladder
  REAL8 tempMax = 0.0;
  REAL8 trigSNR = 0.0;
  REAL8 targetHotLike       = 15;               // Targeted max 'experienced' log(likelihood) of hottest chain
  INT4  hotThreshold        = nChain/2-1;         // If MPIrank > hotThreshold, use proposals with higher acceptance rates for hot chains

  /* Set maximum temperature (command line value take precidence) */
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempMax")) {
    tempMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");
  } else if (LALInferenceGetProcParamVal(runState->commandLine,"--trigSNR")) {
    trigSNR = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "trigSNR");
    networkSNRsqrd = trigSNR * trigSNR;
    tempMax = networkSNRsqrd/(2*targetHotLike); // If trigSNR specified, choose max temp so targetHotLike is achieved
    if(MPIrank==0)
      fprintf(stdout,"Trigger SNR of %f specified, setting tempMax to %f.\n", trigSNR, tempMax);
  } else if (networkSNRsqrd > 0.0) {
    tempMax = networkSNRsqrd/(2*targetHotLike); // If injection, choose max temp so targetHotLike is achieved
    if(MPIrank==0)
      fprintf(stdout,"Injecting SNR of %f, setting tempMax to %f.\n", sqrt(networkSNRsqrd), tempMax);
  } else {
    tempMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax"); // Otherwise use default
    if(MPIrank==0)
      fprintf(stdout,"No --trigSNR or --tempMax specified, and not injecting a signal. Setting tempMax to default of %f.\n", tempMax);
  }
  LALInferenceSetVariable(runState->algorithmParams, "tempMax", &tempMax);

  if (tempMin > tempMax) {
    fprintf(stdout,"WARNING: tempMin > tempMax.  Forcing tempMin=1.0.\n");
    tempMin = 1.0;
    LALInferenceSetVariable(runState->algorithmParams, "tempMin", &tempMin);
  }


  /* Parallel tempering settings */
  INT4 tempSwaps        = (nChain-1)*nChain/2;                // Number of proposed swaps between temperatures in one swap iteration
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempSwaps"))
    tempSwaps = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempSwaps")->value);

  INT4  annealStart     = 500000;                             // Iteration where annealing starts
  if (LALInferenceGetProcParamVal(runState->commandLine,"--annealStart"))
    annealStart = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--annealStart")->value);

  INT4  annealLength    = 100000;                             // Number of iterations to cool temperatures to ~1.0
  if (LALInferenceGetProcParamVal(runState->commandLine,"--annealLength"))
    annealLength = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--annealLength")->value);

  INT4 Tskip            = 100;                                // Number of iterations between proposed temperature swaps 
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip"))
    Tskip = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip")->value);

  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--anneal");
  if (ppt) {
    annealingOn = 1;                                          // Flag to indicate annealing is being used during the run
  }

  INT4 Tkill            = Niter;                              // Iteration where parallel tempering ends 
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempKill"))
    Tkill = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempKill")->value);
  else if (annealingOn)
    Tkill = annealStart+annealLength;

  for (t=0; t<nChain; ++t) {
    tempLadder[t] = 0.0;
    acceptanceCountLadder[t] = 0;
  }

  if (MPIrank == 0) {
    TcurrentLikelihood = (double*) malloc(sizeof(double)*nChain);
  }


  if (runState->likelihood==&LALInferenceTimeDomainLogLikelihood) {
    fprintf(stderr, "Computing null likelihood in time domain.\n");
    nullLikelihood = LALInferenceTimeDomainNullLogLikelihood(runState->data);
  } else if (runState->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood ||
      runState->likelihood==&LALInferenceFreqDomainLogLikelihood) {
    nullLikelihood = LALInferenceNullLogLikelihood(runState->data);
  } else if (runState->likelihood==&LALInferenceFreqDomainStudentTLogLikelihood) {
    REAL8 d = *(REAL8 *)LALInferenceGetVariable(runState->currentParams, "distance");
    REAL8 bigD = 1.0 / 0.0;

    LALInferenceSetVariable(runState->currentParams, "distance", &bigD);
    nullLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
    LALInferenceSetVariable(runState->currentParams, "distance", &d);
  } else if (runState->likelihood==&LALInferenceZeroLogLikelihood) {
    nullLikelihood = 0.0;
  } else if (runState->likelihood==&LALInferenceCorrelatedAnalyticLogLikelihood) {
    nullLikelihood = 0.0;
  } else if (runState->likelihood==&LALInferenceBimodalCorrelatedAnalyticLogLikelihood) {
    nullLikelihood = 0.0;
  } else {
    fprintf(stderr, "Unrecognized log(L) function (in %s, line %d)\n",
        __FILE__, __LINE__);
    exit(1);
  }

  // initialize starting likelihood value:
  runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
  LALInferenceIFOData *headData = runState->data;
  while (headData != NULL) {
    headData->acceptedloglikelihood = headData->loglikelihood;
    headData = headData->next;
  }
  runState->currentPrior = runState->prior(runState, runState->currentParams);

  LALInferenceAddVariable(runState->algorithmParams, "nChain", &nChain,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->algorithmParams, "nPar", &nPar,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs, "parameter",&parameter, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "nullLikelihood", &nullLikelihood, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
  REAL8 logLAtAdaptStart = runState->currentLikelihood;

  /* Construct temperature ladder */
  if(nChain > 1){
    if(LALInferenceGetProcParamVal(runState->commandLine, "--inverseLadder")) {     //temperature spacing uniform in 1/T
      tempDelta = (1.0/tempMin - 1.0/tempMax)/(REAL8)(nChain-1);
      for (t=0; t<nChain; ++t) {
        tempLadder[t]=1.0/(REAL8)(1.0/tempMin-t*tempDelta);
      }
    } else {                                                                        //Geometric spacing
      tempDelta=pow(tempMax-tempMin+1,1.0/(REAL8)(nChain-1));
      for (t=0;t<nChain; ++t) {
        tempLadder[t]=tempMin + pow(tempDelta,t) - 1.0;
        annealDecay[t] = (tempLadder[t]-1.0)/(REAL8)annealLength;
      }
    }
  } else {
    if(LALInferenceGetProcParamVal(runState->commandLine,"--tempMax")){
      tempLadder[0]=tempMax;
    }else{
      tempLadder[0]=1.0;
      tempMax=1.0;
    }
  }

  if (MPIrank > hotThreshold) {
    hotChain = 1;
  }

  LALInferenceAddVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]),  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "hotChain", &hotChain, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

  if (MPIrank == 0){
    printf("\nTemperature ladder:\n");
    for (t=0; t<nChain; ++t) {
      printf(" tempLadder[%d]=%f\n",t,tempLadder[t]);
    }
  }


  FILE * chainoutput = NULL;

  FILE *statfile = NULL;
  FILE *propstatfile = NULL;
  FILE *tempfile = NULL;
  char statfilename[256];
  char propstatfilename[256];
  char tempfilename[256];
  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--tempVerbose")) {
      sprintf(tempfilename,"PTMCMC.tempswaps.%u",randomseed);
      tempfile = fopen(tempfilename, "w");
      /* Print header */
      fprintf(tempfile, "cycle\tlog(chain_swap)\ttemp_low\ttemp_high\n");
    }
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")) {
    sprintf(statfilename,"PTMCMC.statistics.%u.%2.2d",randomseed,MPIrank);
    statfile = fopen(statfilename, "w");
    /* Print header */
    fprintf(statfile,"cycle\ts_gamma");
    ptr=runState->currentParams->head;
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

  if (LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose")) {
    sprintf(propstatfilename,"PTMCMC.propstats.%u.%2.2d",randomseed,MPIrank);
    propstatfile = fopen(propstatfilename, "w");
  }

  chainoutput = LALInferencePrintPTMCMCHeader(runState);


  UINT4 adaptStart = 0;
  REAL8 s_gamma = 1.0;
  REAL8 logLAtAdaptStart = runState->currentLikelihood;

  if (adaptationOn == 1) {
    LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    sigmas = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, SIGMAVECTORNAME));
    PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
    PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
  }

  /* Print run details */
  if (MPIrank == 0) {
    printf("\nParallel Behavior:\n");
    if (adaptationOn)
      printf(" Adapting with decay power %i for %i iterations after max log(L) increases by nParams/2 (%1.2f).\n", adaptTau, adaptationLength, (double)nPar/2.0);
    else
      printf(" Adaptation off.\n");
    if (annealingOn)
      printf(" Annealing linearly starting at iteration %i for %i iterations.\n", annealStart, annealLength);
    else
      printf(" Annealing off.\n");
    if (Tkill != Niter)
      printf(" Parallel tempering for %i iterations.\n", Tkill);
    else
      printf(" Parallel tempering for the entire run.\n");
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

  // iterate:
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=1; i<=Niter; i++) {


    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));

    if (adaptationOn) {
      /* if maximum logL has increased by more than nParam/2, restart it */
      if (runState->currentLikelihood > logLAtAdaptStart+nPar/2) {
        for (p=0; p<nPar; ++p) {
          PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
          PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
          PacceptCount->data[p] =0;
          PproposeCount->data[p]=0;
        }
        adaptStart = i;
        logLAtAdaptStart = runState->currentLikelihood;
        if (!adapting) {
          adapting = 1;
          LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
          LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
      }
      if (adapting) {
        /* Turn off adaption after after adaptationLength step without restarting */
        if ((i-adaptStart) > adaptationLength) {
          adapting = 0;  //turn off adaptation
          LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
          LALInferenceRemoveVariable(runState->proposalArgs,"s_gamma");
          fprintf(stdout,"Ending adaptation for temperature %u at iteration %u.\n",MPIrank,i);

        /* Else set adaptation envelope */
        } else if (i-adaptStart < adaptResetBuffer) {
          s_gamma=(((double)i-(double)adaptStart)/(double)adaptResetBuffer)*(((double)i-(double)adaptStart)/(double)(adaptResetBuffer));
        } else if (i-adaptStart == adaptResetBuffer) {
          s_gamma=1;
        } else {
          s_gamma=10.0*exp(-(1.0/adaptTau)*log((double)(i-adaptStart)))-1;
        }

        /* If adaptation hasn't been turned off, set envelope variable */
        if (adapting)
          LALInferenceSetVariable(runState->proposalArgs, "s_gamma", &s_gamma);
      }
    }

    if (annealingOn) {
      /* Annealing */
      if (i == annealStart) {
        runState->proposal = &LALInferencePostPTProposal;
        if (!LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution"))
          BcastDifferentialEvolutionPoints(runState, 0);
        if (adaptationOn && tempLadder[MPIrank] != 1.0) {
          /* Force hot chains to re-adapt */
          if (!adapting) {
            adapting = 1;
            LALInferenceSetVariable(runState->proposalArgs, "adapting", &adapting);
            LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
          }
          for (p=0; p<nPar; ++p) {
            PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
            PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
            PacceptCount->data[p] =0;
            PproposeCount->data[p]=0;
          }
          s_gamma = 0.0;
          LALInferenceSetVariable(runState->proposalArgs, "s_gamma", &s_gamma);
          adaptStart = i;
        }
      } else if (i > annealStart) {
        for (t=0;t<nChain; ++t) {
          tempLadder[t] = tempLadder[t] - annealDecay[t];
          LALInferenceSetVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]));
        }
        if (annealLength == i - annealStart)
          annealingOn = 0;
      }
    }

    runState->evolve(runState); //evolve the chain with the parameters TcurrentParams[t] at temperature tempLadder[t]
    acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

    if (i==1){
      ppt = LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose");
      if (ppt) {
        // Make sure numbers are initialized!!!
        LALInferenceProposalStatistics *propStat;
        LALInferenceVariableItem *this;
        this = runState->proposalStats->head;
        while(this){
          propStat = (LALInferenceProposalStatistics *)this->value;
          propStat->accepted = 0;
          propStat->proposed = 0;
          this = this->next;
        }
        fprintf(propstatfile, "cycle\t");
        LALInferencePrintProposalStatsHeader(propstatfile, runState->proposalStats);
        fflush(propstatfile);
      }
    }

    if ((i % Nskip) == 0) {
      if (!LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution")) {
        accumulateDifferentialEvolutionSample(runState);
      }

      if (LALInferenceGetProcParamVal(runState->commandLine, "--kDTree")) {
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

      fprintf(chainoutput,"\n");
      fflush(chainoutput);

      if (adaptationOn == 1) {
        if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")) {
          fseek(statfile, 0L, SEEK_END);
          fprintf(statfile,"%d\t",i);

          if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")){
            fprintf(statfile,"%f\t",s_gamma);
            for (p=0; p<nPar; ++p) {
              fprintf(statfile,"%f\t",sigmas->data[p]);
            }
            for (p=0; p<nPar; ++p) {
              fprintf(statfile,"%f\t",PacceptCount->data[p]/( PproposeCount->data[p]==0 ? 1.0 : PproposeCount->data[p] ));
            }
          }
          fprintf(statfile,"\n");
          fflush(statfile);
        }
      }

      if (LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose")){
        fprintf(propstatfile, "%d\t", i);
        LALInferencePrintProposalStats(propstatfile,runState->proposalStats);
        fflush(propstatfile);
      }
    }

    if ((i % Tskip) == 0) {
      ptr=runState->currentParams->head;
      p=0;
      while(ptr!=NULL) {
        if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
          parameters->data[p]=*(REAL8 *)ptr->value;
          p++;
        }
        ptr=ptr->next;
      }

      if (i <= Tkill) {
        ptr=runState->currentParams->head;
        p=0;
        while(ptr!=NULL) {
          if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
            parameters->data[p]=*(REAL8 *)ptr->value;

            p++;
          }
          ptr=ptr->next;
        }

        MPI_Gather(&(runState->currentLikelihood), 1, MPI_DOUBLE, TcurrentLikelihood, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&acceptanceCount, 1, MPI_INT, acceptanceCountLadder, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(parameters->data,nPar,MPI_DOUBLE,parametersVec,nPar,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (MPIrank == 0) { //swap parameters and likelihood between chains
          if(LALInferenceGetProcParamVal(runState->commandLine, "--oldPT")) {
            for(lowerRank=0;lowerRank<nChain-1;lowerRank++) { //swap parameters and likelihood between chains
              for(upperRank=lowerRank+1;upperRank<nChain;upperRank++) {
                LALInferencePTswap(runState, TcurrentLikelihood, parametersVec, tempLadder, lowerRank, upperRank, i, tempfile);
              } //for(upperRank=lowerRank+1;upperRank<nChain;upperRank++)
            } //for(lowerRank=0;lowerRank<nChain-1;lowerRank++)
          } else {
            for(swapAttempt=0; swapAttempt<tempSwaps; ++swapAttempt) {
              lowerRank = gsl_rng_uniform_int(runState->GSLrandom, nChain-1);
              upperRank = lowerRank+1;
              LALInferencePTswap(runState, TcurrentLikelihood, parametersVec, tempLadder, lowerRank, upperRank, i, tempfile);
            } //for(swapAttempt=0; swapAttempt<50; ++swapAttempt)
          } //else
        } //if (MPIrank == 0)

        MPI_Scatter(parametersVec,nPar,MPI_DOUBLE,parameters->data,nPar,MPI_DOUBLE,0, MPI_COMM_WORLD);
        MPI_Scatter(TcurrentLikelihood, 1, MPI_DOUBLE, &(runState->currentLikelihood), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(acceptanceCountLadder, 1, MPI_INT, &acceptanceCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        ptr=runState->currentParams->head;
        p=0;
        while(ptr!=NULL) {
          if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
            memcpy(ptr->value,&(parameters->data[p]),LALInferenceTypeSize[ptr->type]);
            p++;
          }
          ptr=ptr->next;
        }

        MPI_Barrier(MPI_COMM_WORLD);
      }// if (i <= Tkill)
    }// if ((i % Tskip) == 0)
  }// for (i=1; i<=Niter; i++)

  MPI_Barrier(MPI_COMM_WORLD);

  fclose(chainoutput);

  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")) {
      fclose(statfile);
    }
    if (LALInferenceGetProcParamVal(runState->commandLine, "--tempVerbose")) {
      fclose(tempfile);
    }
    if (LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose")) {
      fclose(propstatfile);
    }
  }

  free(tempLadder);
  free(acceptanceCountLadder);
  free(annealDecay);
  free(parametersVec);

  if (MPIrank == 0) {
    free(TcurrentLikelihood);
  }
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
  REAL8 acceptanceRate = 0.0;
  INT4 acceptanceCount;
  INT4 accepted = 0;
  const char *currentProposalName;
  LALInferenceProposalStatistics *propStat;

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
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
  else
    logLikelihoodProposed = -DBL_MAX;

  // determine acceptance probability:
  //printf("%f\t%f\n",logPriorProposed, logLikelihoodProposed);
  logAcceptanceProbability = (1.0/temperature)*(logLikelihoodProposed - logLikelihoodCurrent)
    + (logPriorProposed - logPriorCurrent)
    + logProposalRatio;

  // accept/reject:
  if ((logAcceptanceProbability > 0)
      || (log(gsl_rng_uniform(runState->GSLrandom)) < logAcceptanceProbability)) {   //accept
    LALInferenceCopyVariables(&proposedParams, runState->currentParams);
    runState->currentLikelihood = logLikelihoodProposed;
    LALInferenceIFOData *headData = runState->data;
    while (headData != NULL) {
      headData->acceptedloglikelihood = headData->loglikelihood;
      headData = headData->next;
    }
    runState->currentPrior = logPriorProposed;
    acceptanceCount++;
    accepted = 1;
    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &acceptanceCount);

  }

  INT4 adaptableStep = 0;
  INT4 adapting = 0;
  INT4 i = 0;
  adaptableStep = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptableStep"));
  adapting = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));
  if (adaptableStep && adapting) {
    i = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedArrayNumber"));
    REAL8Vector *PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
    REAL8Vector *PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
    PproposeCount->data[i]+=1;
    if(accepted == 1){
      PacceptCount->data[i]+=1;
    }
    acceptanceRate = PacceptCount->data[i] / PproposeCount->data[i];
  }
  /* Update proposal statistics */
  if (runState->proposalStats){
    currentProposalName = *((const char **)LALInferenceGetVariable(runState->proposalArgs, LALInferenceCurrentProposalName));
    propStat = ((LALInferenceProposalStatistics *)LALInferenceGetVariable(runState->proposalStats, currentProposalName));
    propStat->proposed++;
    if (accepted == 1){
      propStat->accepted++;
    }
  }

  /* Adapt if desired. */
  if (LALInferenceCheckVariable(runState->proposalArgs, "proposedArrayNumber") &&
      LALInferenceCheckVariable(runState->proposalArgs, "proposedVariableNumber") &&
      LALInferenceCheckVariable(runState->proposalArgs, "s_gamma") &&
      LALInferenceCheckVariable(runState->proposalArgs, SIGMAVECTORNAME) &&
      LALInferenceCheckVariable(runState->proposalArgs, "adapting") &&
      LALInferenceCheckVariable(runState->proposalArgs, "adaptableStep")) {

    adaptableStep = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptableStep"));
    if (adaptableStep) {

      i = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedArrayNumber"));
      INT4 varNr = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedVariableNumber"));
      REAL8 s_gamma = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "s_gamma");
      REAL8Vector *sigmas = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, SIGMAVECTORNAME));

      REAL8 sigma = sigmas->data[i];
      char *name = LALInferenceGetVariableName(&proposedParams, varNr);

      char nameMin[VARNAME_MAX], nameMax[VARNAME_MAX];
      REAL8 priorMin, priorMax, dprior;

      sprintf(nameMin, "%s_min", name);
      sprintf(nameMax, "%s_max", name);

      priorMin = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMin));
      priorMax = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMax));

      dprior = priorMax - priorMin;

      if(accepted == 1){
        sigma=sigma+s_gamma*(dprior/100.0)*(1.0-targetAcceptance);
      }else{
        sigma=sigma-s_gamma*(dprior/100.0)*(targetAcceptance);
      }

      sigma = (sigma > dprior ? dprior : sigma);
      sigma = (sigma < 0 ? 0 : sigma);

      sigmas->data[i] = sigma;

      /* Make sure we don't do this again until we take another adaptable step.*/
    }
  }
  adaptableStep = 0;
  LALInferenceSetVariable(runState->proposalArgs, "adaptableStep", &adaptableStep);
  LALInferenceDestroyVariables(&proposedParams);
}


//-----------------------------------------
// temperature swap routine:
//-----------------------------------------
void LALInferencePTswap(LALInferenceRunState *runState,
                        double *TcurrentLikelihood,
                        REAL8 *parametersVec,
                        REAL8 *tempLadder,
                        int lowerRank,
                        int upperRank,
                        int i,
                        FILE *tempfile)
{
  REAL8 logChainSwap;
  REAL8 dummyR8;
  INT4 p;
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);

  logChainSwap = (1.0/tempLadder[lowerRank]-1.0/tempLadder[upperRank]) * (TcurrentLikelihood[upperRank]-TcurrentLikelihood[lowerRank]);

  if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) { //Then swap...
    // Check if --tempVerbose was specified
    if (tempfile != NULL) {
      fprintf(tempfile,"%d\t%f\t%f\t%f\n",i,logChainSwap,tempLadder[lowerRank],tempLadder[upperRank]);
      fflush(tempfile);
    }
    for (p=0; p<(nPar); ++p){
      dummyR8=parametersVec[p+nPar*upperRank];
      parametersVec[p+nPar*upperRank]=parametersVec[p+nPar*lowerRank];
      parametersVec[p+nPar*lowerRank]=dummyR8;
    }
  dummyR8 = TcurrentLikelihood[upperRank];
  TcurrentLikelihood[upperRank] = TcurrentLikelihood[lowerRank];
  TcurrentLikelihood[lowerRank] = dummyR8;
  }
}


//-----------------------------------------
// file output routines:
//-----------------------------------------

FILE* LALInferencePrintPTMCMCHeader(LALInferenceRunState *runState)
{

  int MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  REAL8 tempMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");
  UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");
  REAL8 nullLikelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  REAL8 temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");
  REAL8 SampleRate=4096.0; //default value of the sample rate from LALInferenceReadData()
  UINT4 nIFO=0;
  LALInferenceIFOData *ifodata1=runState->data;
  ProcessParamsTable *ppt;

  while(ifodata1){
    nIFO++;
    ifodata1=ifodata1->next;
  }

  FILE * chainoutput = NULL;
  char *outfileName = NULL;

  int waveform = 0;
  if(LALInferenceCheckVariable(runState->currentParams,"LAL_APPROXIMANT")) waveform= *(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_APPROXIMANT");
  double pnorder = 0.0;
  if(LALInferenceCheckVariable(runState->currentParams,"LAL_PNORDER")) pnorder = ((double)(*(INT4 *)LALInferenceGetVariable(runState->currentParams,"LAL_PNORDER")))/2.0;

  char str[999];
  LALInferencePrintCommandLine(runState->commandLine, str);

  REAL8 networkSNR=0.0;
  ifodata1=runState->data;
  while(ifodata1){
    networkSNR+=ifodata1->SNR*ifodata1->SNR;
    ifodata1=ifodata1->next;
  }
  networkSNR=sqrt(networkSNR);

  if(LALInferenceGetProcParamVal(runState->commandLine,"--srate")) SampleRate=atof(LALInferenceGetProcParamVal(runState->commandLine,"--srate")->value);

  outfileName = (char*)calloc(99,sizeof(char*));
  ppt = LALInferenceGetProcParamVal(runState->commandLine, "--appendOutput");
  if (ppt) {
    sprintf(outfileName, "%s.%2.2d", ppt->value, MPIrank);
  } else {
    sprintf(outfileName,"PTMCMC.output.%u.%2.2d",randomseed,MPIrank);
  }
  if (!ppt) { /* Skip header output if we are appending. */
    chainoutput = fopen(outfileName,"w");
    fprintf(chainoutput, "  LALInference version:%s,%s,%s,%s,%s\n", LALAPPS_VCS_ID,LALAPPS_VCS_DATE,LALAPPS_VCS_BRANCH,LALAPPS_VCS_AUTHOR,LALAPPS_VCS_STATUS);
    fprintf(chainoutput,"  %s\n",str);
    fprintf(chainoutput, "%10s  %10s  %6s  %20s  %6s %8s   %6s  %8s  %10s  %12s  %9s  %9s  %8s\n",
        "nIter","Nburn","seed","null likelihood","Ndet","nCorr","nTemps","Tmax","Tchain","Network SNR","Waveform","pN order","Npar");
    fprintf(chainoutput, "%10d  %10d  %u  %20.10lf  %6d %8d   %6d%10d%12.1f%14.6f  %9i  %9.1f  %8i\n",
        Niter,0,randomseed,nullLikelihood,nIFO,0,nChain,(int)tempMax,temperature,networkSNR,waveform,(double)pnorder,nPar);
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
    fprintf(chainoutput,"\n");
    fprintf(chainoutput, "%d\t%f\t%f\t", 0,(runState->currentLikelihood - nullLikelihood)+runState->currentPrior, runState->currentPrior);
    LALInferencePrintSampleNonFixed(chainoutput,runState->currentParams);
    fprintf(chainoutput,"%f\t",runState->currentLikelihood - nullLikelihood);
    headIFO = runState->data;
    while (headIFO != NULL) {
      fprintf(chainoutput, "%f\t", headIFO->acceptedloglikelihood - headIFO->nullloglikelihood);
      headIFO = headIFO->next;
    }

    fprintf(chainoutput,"\n");
    fclose(chainoutput);
  }

  chainoutput = fopen(outfileName,"a");
  free(outfileName);

  return chainoutput;
}

void LALInferenceDataDump(LALInferenceRunState *runState){

  const UINT4 nameLength=256;
  char filename[nameLength];
  FILE *out;
  LALInferenceIFOData *headData = runState->data;
  UINT4 ui;

  while (headData != NULL) {

    snprintf(filename, nameLength, "%s-freqModelhPlus.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->freqModelhPlus->data->length; ui++) {
      REAL8 f = headData->freqModelhPlus->deltaF * ui;
      COMPLEX16 d = headData->freqModelhPlus->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-freqModelhCross.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->freqModelhCross->data->length; ui++) {
      REAL8 f = headData->freqModelhCross->deltaF * ui;
      COMPLEX16 d = headData->freqModelhCross->data->data[ui];

      fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeModelhPlus.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->timeModelhPlus->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhPlus->epoch)) +
        ui * headData->timeModelhPlus->deltaT;
      REAL8 d = headData->timeModelhPlus->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeModelhCross.dat", headData->name);
    out = fopen(filename, "w");
    for (ui = 0; ui < headData->timeModelhCross->data->length; ui++) {
      REAL8 tt = XLALGPSGetREAL8(&(headData->timeModelhCross->epoch)) +
        ui * headData->timeModelhCross->deltaT;
      REAL8 d = headData->timeModelhCross->data->data[ui];

      fprintf(out, "%.6f %g\n", tt, d);
    }
    fclose(out);

    snprintf(filename, nameLength, "%s-timeModel.dat", headData->name);
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
