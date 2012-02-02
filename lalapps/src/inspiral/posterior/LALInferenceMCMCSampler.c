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

void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState)
{
  int i,t,p,lowerRank,upperRank,x; //indexes for for() loops
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
  INT4 **pdf = NULL;
  INT4 pdf_count = 0;
  INT4 param_count = 0;
  INT4 parameter=0;
  UINT4 tMaxSearch = 0;
  UINT4 hotChain = 0;                 // Affects proposal setup
  REAL8Vector *sigmas = NULL;
  REAL8Vector *PacceptCount = NULL;
  REAL8Vector *PproposeCount = NULL;
  REAL8 *parametersVec = NULL;
  REAL8 flatPriorTestVal = 0.0;
  REAL8 randVal = 0.0;
  REAL8 paramVal = 0.0;
  REAL8 tempCurrentPrior = 0.0;
  REAL8 tempCurrentLikelihood = 0.0;
  REAL8 priorMin, priorMax, dprior;
  REAL8Vector * parameters = NULL;
  LALInferenceVariables tempCurrentParams;
  LALInferenceVariables flatPriorTestParams;
  LALInferenceVariables flatPriorParams;
  LALInferenceProposalFunction *tempProposal;
  LALInferencePriorFunction *tempPrior;
  char *name = NULL;
  char nameMin[VARNAME_MAX], nameMax[VARNAME_MAX];

  INT4 adaptationOn = 0;
  INT4 annealingOn= 0;
  INT4 acceptanceRatioOn = 0;
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


  /* Adaptation settings */
  INT4 adaptTau           = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptTau"));
  INT4 adaptResetBuffer   = 100;              // Number of iterations before adapting after a restart
  UINT4 adaptationLength  = pow(10,adaptTau);   // Number of iterations to adapt before turning off
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatio");
  if(ppt){
    acceptanceRatioOn = 1;
  }
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--adapt");
  if (ppt) {
    adaptationOn = 1;
  }

  /* Temperature ladder settings */
  UINT4 hotThreshold        = 2;                // If MPIrank > hotThreshold, use different "hot" jump proposals
  REAL8 tempMin = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMin");   //min temperature in the temperature ladder
  REAL8 tempMax = *(REAL8*) LALInferenceGetVariable(runState->algorithmParams, "tempMax");   //max temperature in the temperature ladder
  INT4  TmaxSearchLen       = 10000;            // Lentgh of mini-MCMC testing for ideal maximum temperature
  REAL8 tempSearchLow       = tempMin;          // Lower bound for maximum temperature search
  REAL8 tempSearchHigh      = tempMax;          // Upper bound for maximum temperature search
  INT4  nFlatPriorBins      = 10;               // Number of bins for params w/ flat priors for max temp test
  REAL8 flatPriorTolerance  = .3;               // Percentage of tolerance allowed in each bin when testing for flatness
  REAL8 tempDelta           = (tempSearchHigh-tempSearchLow)/(REAL8)(nChain-1);
  INT4 flatBinLow  = (INT4) (TmaxSearchLen / nFlatPriorBins * (1.0 - flatPriorTolerance));
  INT4 flatBinHigh = (INT4) (TmaxSearchLen / nFlatPriorBins * (1.0 + flatPriorTolerance));

  /* Annealing settings */
  INT4 startAnnealing       = 200000;           // Iteration where annealing starts
  INT4 annealLength         = 20000;            // Number of iterations to cool temperatures to ~1.0
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--anneal");
  if (ppt) {
    adaptationOn = 1;
  }

  /* Parallel tempering settings */
  INT4 nSwaps = (nChain-1)*nChain/2;            // Number of proposed swaps between temperatures
  INT4 Tskip = 100;                             // Number of iterations between temperature swaps proposals
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip"))
    Tskip = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempSkip")->value);
  INT4 Tkill = Niter;                             // Iteration number at which temperature swaps are no longer proposed
  if (LALInferenceGetProcParamVal(runState->commandLine,"--tempKill"))
    Tkill = atoi(LALInferenceGetProcParamVal(runState->commandLine,"--tempKill")->value);

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
  LALInferenceAddVariable(runState->proposalArgs, "tMaxSearch", &tMaxSearch, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
  LALInferenceAddVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]),  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "hotChain", &hotChain, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

  if (!LALInferenceGetProcParamVal(runState->commandLine, "--noTempSearch") && nChain > 1) {
    /*
     * Determine how high of a temperature is needed to recover the prior.
     *
     * A linear temperature ladder is constructed, and the parameters with
     * flat priors are binned and checked for flatness.
     */
    tMaxSearch = 1;
    LALInferenceSetVariable(runState->proposalArgs, "tMaxSearch", &(tMaxSearch));

    /* Save values for after temperature testing */
    tempCurrentParams.head      = NULL;
    tempCurrentParams.dimension = 0;
    LALInferenceCopyVariables(runState->currentParams, &tempCurrentParams);
    tempCurrentPrior = runState->currentPrior;
    tempCurrentLikelihood = runState->currentLikelihood;
    tempPrior = runState->prior;
    tempProposal = runState->proposal;

    /* Find parameters with flat prior */
    INT4 nFlatPar = nPar;
    flatPriorTestParams.head      = NULL;
    flatPriorTestParams.dimension = 0;
    flatPriorParams.head          = NULL;
    flatPriorParams.dimension     = 0;
    runState->prior = &LALInferenceInspiralPriorNormalised;
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    flatPriorTestVal = runState->currentPrior;
    LALInferenceCopyVariables(runState->currentParams, &flatPriorParams);

    for(p=0;p<nPar;++p){
      LALInferenceCopyVariables(runState->currentParams, &flatPriorTestParams);

      name = LALInferenceGetVariableName(runState->currentParams, (p+1));
      sprintf(nameMin, "%s_min", name);
      sprintf(nameMax, "%s_max", name);

      priorMin = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMin));
      priorMax = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMax));
      dprior = priorMax - priorMin;

      for(x=0;x<10000;++x){
        randVal = gsl_rng_uniform(runState->GSLrandom);
        paramVal = priorMin + randVal * (priorMax - priorMin);
        LALInferenceSetVariable(&flatPriorTestParams, name, &paramVal);
        flatPriorTestVal = runState->prior(runState, &flatPriorTestParams);
        if(flatPriorTestVal != runState->currentPrior){
          LALInferenceRemoveVariable(&flatPriorParams, name);
          nFlatPar -= 1;
          break;
        }
      }
    }

    /* Construct temporary linear temperature ladder to probe for best max temp */
    for(t=0; t<nChain; ++t){
      tempLadder[t]=tempSearchLow+t*tempDelta;
    }

    LALInferenceSetVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]));

    /* Use specialized jump proposal and run a short MCMC */
    if(MPIrank==0)
      fprintf(stdout,"Running exploratory MCMC to determine best temperature ladder.\n");

    runState->proposal = &LALInferencePTTempTestProposal;
    runState->prior = &LALInferenceInspiralPriorNormalised;
    pdf=(UINT4**)calloc(nPar,sizeof(INT4 *));

    while (tMaxSearch == 1) {
      for(p=0;p<nPar;++p){
        pdf[p]=calloc(10,sizeof(INT4));
        for(x=0;x<10;++x){
          pdf[p][x]=0;
        }
      }


      for(i=0;i<TmaxSearchLen;++i){
        PTMCMCOneStep(runState);

        ptr=runState->currentParams->head;
        p=0;
        while(ptr!=NULL) {
          if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
            parameters->data[p]=*(REAL8 *)ptr->value;
            p++;
          }
          ptr=ptr->next;
        }

        /* Bin parameterm values */
        for (p=0;p<nPar;++p){
          name = LALInferenceGetVariableName(runState->currentParams, (p+1));
          sprintf(nameMin, "%s_min", name);
          sprintf(nameMax, "%s_max", name);
          priorMin = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMin));
          priorMax = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, nameMax));
          dprior = priorMax - priorMin;
          x=(int)(((parameters->data[p] - priorMin)/dprior)*nFlatPriorBins);
          if(x<0) x=0;
          if(x>nFlatPriorBins-1) x=nFlatPriorBins-1;
          pdf[p][x]++;
        }
      }//for(i=0;i<TmaxSearchLen;++i)

      /* Check for flat PDFs in parameters w/ flat priors */
      param_count=0;
      for (p=0;p<nPar;++p){
        name = LALInferenceGetVariableName(runState->currentParams, (p+1));
        if(LALInferenceCheckVariable(&flatPriorParams, name)){
          pdf_count=0;
          for(x=0;x<nFlatPriorBins;++x){
            if(pdf[p][x]<flatBinLow || pdf[p][x]>flatBinHigh) pdf_count++;
          }
          if(pdf_count==0) param_count++;
        }
      }

      if (param_count == nFlatPar) {
        acceptanceCount = 1;
      } else {
        acceptanceCount = 0;
      }

      MPI_Allgather(&acceptanceCount, 1, MPI_INT, acceptanceCountLadder, 1, MPI_INT, MPI_COMM_WORLD);

      UINT4 recoveredPrior = 0;
      tempMax = tempLadder[nChain-1];
      for (i=0;i<nChain;++i) {
        if (acceptanceCountLadder[i]) {
          if (!recoveredPrior) {
            recoveredPrior = 1;
            tempMax = tempLadder[i];
            tMaxSearch = 0;
          }
        } else {
          if (recoveredPrior) {
            if(MPIrank==0)
              fprintf(stdout,"Inconsistent temperature performance, possibly due to stuck chain.  Re-running exploritory MCMC");
            recoveredPrior = 0;
            tempMax = tempLadder[nChain-1];
            tMaxSearch = 1;
            break;
          }
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    } //while (tMaxSearch == 1)

    if (tempMax == tempLadder[nChain-1])
      fprintf(stdout,"WARNING: The search set max temperature to the maximum allowed temperature (%f). \
              This may be insufficient. Recommend allowing higher temperatures using --Tmax=<Tmax>.\n",tempMax);

    LALInferenceSetVariable(runState->algorithmParams, "tempMax", &(tempMax));
    LALInferenceSetVariable(runState->proposalArgs, "tMaxSearch", &(tMaxSearch));

    /* Reset values to those before temperature search */
    LALInferenceCopyVariables(&tempCurrentParams, runState->currentParams);
    runState->currentLikelihood = tempCurrentLikelihood;
    runState->prior = tempPrior;
    runState->currentPrior = tempCurrentPrior;
    runState->proposal = tempProposal;
    acceptanceCount = 0;
    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));
    LALInferenceDeleteProposalCycle(runState);
    for (t=0; t<nChain; ++t) {
      acceptanceCountLadder[t] = 0;
    }
  }//if(!LALInferenceGetProcParamVal(runState->commandLine, "--noTempSearch") && nChain>1)

  /* Construct temperature ladder */
  if(nChain > 1){
    if(LALInferenceGetProcParamVal(runState->commandLine, "--inverseLadder")){ //temperature spacing uniform in 1/T
      tempDelta = (1.0 - 1.0/tempMax)/(REAL8)(nChain-1);
      for (t=0; t<nChain; ++t) {
        tempLadder[t]=1.0/(REAL8)(1.0-t*tempDelta);
      }
    }
    else if(LALInferenceGetProcParamVal(runState->commandLine, "--geomLadder")){ //Geometric spacing (most efficient so far. Should become default?
      tempDelta=pow(tempMax,1.0/(REAL8)(nChain-1));
      //tempDelta=pow(tempMax,1.0/(REAL8)(nChain+4-1));
      for (t=0;t<nChain; ++t) {
        tempLadder[t]=pow(tempDelta,t);
        //tempLadder[t]=pow(tempDelta,t+4);
        annealDecay[t] = t*(REAL8)annealLength / (((t-1)/(nChain-1))*log(tempMax - 1)+1);
        //annealDecay[t] = t*(REAL8)annealLength / (((t+4-1)/(nChain+4-1))*log(tempMax - 1)+1);
      }
    }
    else{ //epxonential spacing
      tempDelta = log(tempMax)/(REAL8)(nChain-1);
      for (t=0; t<nChain; ++t) {
        tempLadder[t]=exp(t*tempDelta);
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

  LALInferenceSetVariable(runState->proposalArgs, "temperature", &(tempLadder[MPIrank]));
  LALInferenceSetVariable(runState->proposalArgs, "hotChain", &(hotChain));

  if (MPIrank == 0){
    for (t=0; t<nChain; ++t) {
      printf("tempLadder[%d]=%f\n",t,tempLadder[t]);
    }
  }


  FILE * chainoutput = NULL;

  FILE *stat = NULL;
  FILE *propstatfile = NULL;
  FILE *tempfile = NULL;
  char statfilename[256];
  char propstatfilename[256];
  char tempfilename[256];
  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--tempVerbose")) {
      sprintf(tempfilename,"PTMCMC.tempswaps.%u",randomseed);
      tempfile = fopen(tempfilename, "a");
    }
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
    sprintf(statfilename,"PTMCMC.statistics.%u.%2.2d",randomseed,MPIrank);
    stat = fopen(statfilename, "a");
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose")) {
    sprintf(propstatfilename,"PTMCMC.propstats.%u.%2.2d",randomseed,MPIrank);
    propstatfile = fopen(propstatfilename, "a");
  }

  chainoutput = LALInferencePrintPTMCMCHeader(runState);


  UINT4 adaptStart = 0;
  REAL8 s_gamma = 1.0;
  REAL8 logLAtAdaptStart = runState->currentLikelihood;

  if (adaptationOn == 1) {
    LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
    sigmas = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, SIGMAVECTORNAME));
  }

  if (acceptanceRatioOn == 1){
    PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
    PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
  }

  if (MPIrank == 0) {
    printf(" PTMCMCAlgorithm(); starting parameter values:\n");
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
  MPI_Barrier(MPI_COMM_WORLD);
  for (i=1; i<=Niter; i++) {


    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));

    /* Turn off adaptation if necessary */
    if (adaptationOn) {
      if ((i-adaptStart) > adaptationLength) {
        fprintf(stdout,"Ending adaptation for temperature %u at iteration %u.\n",MPIrank,i);
        LALInferenceRemoveVariable(runState->proposalArgs,"s_gamma");
        adaptationOn = 0;  //turn off adaptation
      }
    }

    /* Annealing */
    if (annealingOn) {
      if (i == startAnnealing) {
        if (!adaptationOn) {
          adaptationOn = 1;
          //s_gamma = 1.0;
          s_gamma = 0.0;
          LALInferenceAddVariable(runState->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        }
        for (p=0; p<nPar; ++p) {
          PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
          PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
          PacceptCount->data[p] =0;
          PproposeCount->data[p]=0;
        }
        adaptStart = i;
      } 

      if (i >= startAnnealing) {
        for (t=0;t<nChain; ++t) {
          tempLadder[t]= pow(tempMax-1, t/(nChain-1))*exp(-(i-startAnnealing)/annealDecay[t]) + 1;
        }
      }
    }

    if (adaptationOn == 1) {
      /* if logL has increased by more than nParam/2 since adaptation began, restart it */
      if (runState->currentLikelihood > logLAtAdaptStart+nPar/2) {
        for (p=0; p<nPar; ++p) {
          PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
          PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
          PacceptCount->data[p] =0;
          PproposeCount->data[p]=0;
        }
        adaptStart = i;
        logLAtAdaptStart = runState->currentLikelihood;
      }

      if (i-adaptStart < adaptResetBuffer) {
        //s_gamma=(((double)i-(double)adaptStart)/(double)adaptResetBuffer)*(((double)i-(double)adaptStart)/(double)(adaptResetBuffer));
        s_gamma=0;
      } else if (i-adaptStart == adaptResetBuffer) {
        s_gamma=1;
      }
      //s_gamma=s_gamma * (10*exp(-(1.0/adaptTau)*log((double)(i-adaptStart)))-1);
      LALInferenceSetVariable(runState->proposalArgs, "s_gamma", &s_gamma);
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
        if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
          fseek(stat, 0L, SEEK_END);
          fprintf(stat,"%d\t",i);

          if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")){
            fprintf(stat,"%f\t",s_gamma);
            for (p=0; p<nPar; ++p) {
              fprintf(stat,"%f\t",sigmas->data[p]);
            }
          }
          if(LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")){
            for (p=0; p<nPar; ++p) {
              fprintf(stat,"%f\t",PacceptCount->data[p]/( PproposeCount->data[p]==0 ? 1.0 : PproposeCount->data[p] ));
            }
          }
          fprintf(stat,"\n");
          fflush(stat);
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

        if (i < startAnnealing) {
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
              for(swapAttempt=0; swapAttempt<nSwaps; ++swapAttempt) {
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
        }
      }// if (i <= Tkill)
    }// if ((i % Tskip) == 0)
  }// for (i=1; i<=Niter; i++)

  MPI_Barrier(MPI_COMM_WORLD);

  fclose(chainoutput);

  if(MPIrank == 0){
    if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose") || LALInferenceGetProcParamVal(runState->commandLine, "--acceptanceRatioVerbose")) {
      fclose(stat);
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
  REAL8 diff = 0.0;
  INT4 acceptanceCount;
  INT4 accepted = 0;
  UINT4 tMaxSearch = 0;
  const char *currentProposalName;
  LALInferenceProposalStatistics *propStat;
  ProcessParamsTable *ppt, *commandLine = runState->commandLine;

  // current values:
  logPriorCurrent      = runState->currentPrior;
  logLikelihoodCurrent = runState->currentLikelihood;

  temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");
  acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");
  tMaxSearch = *(UINT4*) LALInferenceGetVariable(runState->proposalArgs, "tMaxSearch");

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
  INT4 i = 0;
  adaptableStep = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adaptableStep"));
  if (adaptableStep) {
    i = *((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "proposedArrayNumber"));
    REAL8Vector *PacceptCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PacceptCount"));
    REAL8Vector *PproposeCount = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "PproposeCount"));
    PproposeCount->data[i]+=1;
    if(accepted == 1){
      PacceptCount->data[i]+=1;
    }
    acceptanceRate = PacceptCount->data[i] / PproposeCount->data[i];
  }
  /* Update proposal statistics unless we are searching for max temperature */
  if (runState->proposalStats && !tMaxSearch){
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

      diff = acceptanceRate - targetAcceptance;
      if (diff > 0){
        sigma=sigma*(1+s_gamma*diff*diff);
      } else {
        sigma=sigma*(1-s_gamma*diff*diff);
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
