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

void LALInferencePTswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile);
void LALInferenceMCMCMCswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile);
void LALInferenceAdaptation(LALInferenceRunState *runState, INT4 cycle);
void LALInferenceAdaptationRestart(LALInferenceRunState *runState, INT4 cycle);
void LALInferenceAdaptationEnvelope(LALInferenceRunState *runState, INT4 cycle);
FILE* LALInferencePrintPTMCMCHeader(LALInferenceRunState *runState);
void LALInferencePrintPTMCMCHeaderFile(LALInferenceRunState *runState, FILE *file);
void LALInferencePrintPTMCMCInjectionSample(LALInferenceRunState *runState);
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

static void DEbuffer2array(LALInferenceRunState *runState, INT4 startCycle, INT4 endCycle, REAL8** DEarray) {
  LALInferenceVariableItem *ptr;
  INT4 i=0,p=0;

  INT4 Nskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Nskip");
  INT4 totalPoints = runState->differentialPointsLength;
  INT4 start = (INT4)ceil((REAL8)startCycle/(REAL8)Nskip);
  INT4 end = (INT4)floor((REAL8)endCycle/(REAL8)Nskip);
  /* Include last point */
  if (end > totalPoints-1)
    end = totalPoints-1;

  for (i = start; i <= end; i++) {
    ptr=runState->differentialPoints[i]->head;
    p=0;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        DEarray[i-start][p]=*(REAL8 *)ptr->value;
        p++;
      }
      ptr=ptr->next;
    }
  }
}

static void
replaceDEbuffer(LALInferenceRunState *runState, REAL8** DEarray) {
  LALInferenceVariableItem *ptr;
  UINT4 i=0,p=0;
  UINT4 nPoints = sizeof(DEarray) / sizeof(REAL8*);

  /* Save last LALInferenceVariables item from buffer to keep fixed params consistent for chain */
  LALInferenceVariables templateParamSet;
  LALInferenceCopyVariables(runState->differentialPoints[runState->differentialPointsLength-1], &templateParamSet);

  /* Free old DE buffer */
  XLALFree(runState->differentialPoints);

  /* Expand DE buffer */
  size_t newSize = runState->differentialPointsSize;
  while (nPoints > newSize) {
    newSize = newSize*2;
  }

  runState->differentialPoints = XLALCalloc(newSize, sizeof(LALInferenceVariables *));
  runState->differentialPointsLength = nPoints;
  runState->differentialPointsSize = newSize;

  for (i=0; i<nPoints; i++) {
    runState->differentialPoints[i] = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(&templateParamSet, runState->differentialPoints[i]);
    ptr = runState->differentialPoints[i]->head;
    while(ptr!=NULL) {
      if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
        *((REAL8 *)ptr->value) = (REAL8)DEarray[i][p];
        p++;
      }
      ptr=ptr->next;
    }
  }
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
    DEbuffer2array(runState, startCycle, endCycle, packedDEsamples);

  /* Send it out */
  MPI_Bcast(packedDEsamples[0], nPoints*nPar, MPI_DOUBLE, sourceTemp, MPI_COMM_WORLD);

  /* Unpack it */
  if (MPIrank != sourceTemp)
    replaceDEbuffer(runState, packedDEsamples);

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

    DEbuffer2array(runState, startCycle, endCycle, DEarray);

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
  } else {
    max = Niter;
  }

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

static REAL8Vector*
CopyLALInferenceVariablesToArray(LALInferenceVariables *origin) {
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(origin);
  REAL8Vector * parameters = NULL;
  gsl_matrix *m = NULL; //for dealing with noise parameters
  UINT4 j,k;

  parameters = XLALCreateREAL8Vector(nPar);

  LALInferenceVariableItem *ptr=origin->head;
  INT4 p=0;
  while(ptr!=NULL) {
    if (ptr->vary != LALINFERENCE_PARAM_FIXED) {
      //Generalized to allow for parameters stored in gsl_matrix
      if(ptr->type == LALINFERENCE_gslMatrix_t)
      {
        m = *((gsl_matrix **)ptr->value);
        for(j=0; j<m->size1; j++)
        {
          for(k=0; k<m->size2; k++)
          {
            parameters->data[p]=gsl_matrix_get(m,j,k);
            p++;
          }
        }
      }
      else
      {
        parameters->data[p]=*(REAL8 *)ptr->value;
        p++;
      }
    }
    ptr=ptr->next;
  }

  return parameters;
}

static void
CopyArrayToLALInferenceVariables(REAL8Vector *origin, LALInferenceVariables *target) {
  gsl_matrix *m = NULL; //for dealing with noise parameters
  UINT4 j,k;

  LALInferenceVariableItem *ptr = target->head;
  INT4 p=0;
  while(ptr!=NULL) {
    if (ptr->vary != LALINFERENCE_PARAM_FIXED)
    {
      //Generalized to allow for parameters stored in gsl_matrix
      if(ptr->type == LALINFERENCE_gslMatrix_t)
      {
        m = *((gsl_matrix **)ptr->value);
        for(j=0; j<m->size1; j++)
        {
          for(k=0; k<m->size2; k++)
          {
            gsl_matrix_set(m,j,k,origin->data[p]);
            p++;
          }
        }
      }
      else
      {
        memcpy(ptr->value,&(origin->data[p]),LALInferenceTypeSize[ptr->type]);
        p++;
      }
    }
    ptr=ptr->next;
  }

  /* update stored noise parameters */
  if(LALInferenceCheckVariable(target,"psdscale"))
  {
    gsl_matrix_memcpy(*((gsl_matrix **)LALInferenceGetVariable(target, "psdstore")),
                      *((gsl_matrix **)LALInferenceGetVariable(target, "psdscale")));
  }
  return;
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
  INT4 *intVec = NULL;
  INT4 annealStartIter = 0;
  INT4 phaseAcknowledged = 0;
  INT4 iEffStart = 0;
  UINT4 hotChain = 0;                 // Affects proposal setup
  REAL8 temp = 1;
  REAL8 tempDelta = 0.0;
  MPI_Request MPIrequest;

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
  REAL8 timestamp,timestamp_epoch;
  struct timeval tv;

  //
  /* Command line flags (avoid repeated checks of runState->commandLine) */
  //
  UINT4 runPhase = 0;    // Phase of run. (0=PT-only run, 1=temporary PT, 2=annealing, 3=single-chain sampling)
  UINT4 annealingOn = 0; // Chains will be annealed
  if (LALInferenceGetProcParamVal(runState->commandLine, "--anneal")) {
    annealingOn = 1;
    runPhase=1;
  }

  UINT4 diffEvo = 1; // Differential evolution
  if (LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution") || LALInferenceGetProcParamVal(runState->commandLine, "--nodifferentialevolution")) {
    diffEvo = 1;
  }

  UINT4 adapting = 0;

  UINT4 MCMCMC=0; // Metropolis-Coupled MCMC
  if (LALInferenceGetProcParamVal(runState->commandLine,"--flowSearch")) {
    MCMCMC=1;
  }

  UINT4 tempVerbose = 0; // Print temperature swaps to file
  if (LALInferenceGetProcParamVal(runState->commandLine, "--tempVerbose")) {
    tempVerbose = 1;
  }

  UINT4 adaptVerbose = 0; // Print adaptation info to file
  if (LALInferenceGetProcParamVal(runState->commandLine, "--adaptVerbose")) {
    adaptVerbose = 1;
  }

  UINT4 propVerbose = 0; // Print proposal acceptance rates to file
  if (LALInferenceGetProcParamVal(runState->commandLine, "--propVerbose")) {
    propVerbose = 1;
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
  intVec = malloc(nChain * sizeof(INT4));

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

  REAL8 tempMaxMin          = 10.0;             // Don't let tempMax go too low
  REAL8 targetHotLike       = 15;               // Targeted max 'experienced' log(likelihood) of hottest chain
  INT4 hotThreshold        = nChain/2-1;       // If MPIrank > hotThreshold, use different proposal set

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
    if (ladderMax < tempMaxMin)
      ladderMax = tempMaxMin;
    if(MPIrank==0)
      fprintf(stdout,"Trigger SNR of %f specified, setting tempMax to %f.\n", trigSNR, ladderMax);
  } else if (networkSNRsqrd > 0.0) {                                                  //injection, choose ladderMax to get targetHotLike
    ladderMax = networkSNRsqrd/(2*targetHotLike);
    if (ladderMax < tempMaxMin)
      ladderMax = tempMaxMin;
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
      runState->likelihood==&LALInferenceFreqDomainLogLikelihood ||
      runState->likelihood==&LALInferenceNoiseOnlyLogLikelihood){
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
  } else if (runState->likelihood==&LALInferenceRosenbrockLogLikelihood) {
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
  LALInferenceSetVariable(runState->proposalArgs, "logLAtAdaptStart", &(logLAtAdaptStart));

  /* Construct temperature ladder */
  if(nChain > 1){
    if(LALInferenceGetProcParamVal(runState->commandLine, "--inverseLadder")) {     //Spacing uniform in 1/T
      tempDelta = (1.0/ladderMin - 1.0/ladderMax)/(REAL8)(nChain-1);
      for (t=0; t<nChain; ++t) {
        ladder[t]=1.0/(REAL8)(1.0/ladderMin-t*tempDelta);
      }
    } else {                                                                        //Geometric spacing
      if (LALInferenceGetProcParamVal(runState->commandLine, "--bottomUp") && nPar != 1)
        tempDelta=(nPar+1.)/(nPar-1.);
      else
        tempDelta=pow(ladderMax-ladderMin+1,1.0/(REAL8)(nChain-1));
      for (t=0;t<nChain; ++t) {
        ladder[t]=ladderMin + pow(tempDelta,t) - 1.0;
      }
    }
  } else {                                                                        //Geometric spacing
    tempDelta=pow(ladderMax/ladderMin,1.0/(REAL8)(nChain-1));
    for (t=0;t<nChain; ++t) {
      ladder[t]=ladderMin*pow(tempDelta,t);
    }
  }
  temp=ladder[MPIrank];

  if (MPIrank > hotThreshold) {
    hotChain = 1;
  }

  LALInferenceAddVariable(runState->proposalArgs, "temperature", &temp,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(runState->proposalArgs, "hotChain", &hotChain, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

  if (MPIrank == 0){
    printf("\nTemperature ladder:\n");
    for (t=0; t<nChain; ++t) {
      printf(" ladder[%d]=%f\n",t,ladder[t]);
    }
  }


  FILE * chainoutput = NULL;

  FILE *statfile = NULL;
  FILE *propstatfile = NULL;
  FILE *swapfile = NULL;
  char statfilename[256];
  char propstatfilename[256];
  char swapfilename[256];
  if(MPIrank == 0){
    if (tempVerbose) {
      sprintf(swapfilename,"PTMCMC.tempswaps.%u",randomseed);
      swapfile = fopen(swapfilename, "w");
      fprintf(swapfile, "cycle\tlog(chain_swap)\ttemp_low\ttemp_high\n");  // Print header for temp stat file
    }
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

  if (propVerbose) {
    sprintf(propstatfilename,"PTMCMC.propstats.%u.%2.2d",randomseed,MPIrank);
    propstatfile = fopen(propstatfilename, "w");
  }

  chainoutput = LALInferencePrintPTMCMCHeader(runState);
  if (MPIrank == 0) {
    LALInferencePrintPTMCMCInjectionSample(runState);
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


  /* MPI tags used:
   *   0: Parallel tempering communications
   *   1: runPhase passing
   *   2: runComplete flag
   */

  /* Setup non-blocking recieve that will update the runPhase when chain 0 does */
  if (annealingOn && MPIrank!=0)
    MPI_Irecv(&runPhase, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &MPIrequest);

  /* Setup non-blocking recieve that will succeed when chain 0 is complete. */
  if (MPIrank!=0)
    MPI_Irecv(&runComplete, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &MPIrequest);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  // iterate:
  i=0;
  while (!runComplete) {
    /* Increment iteration counter */
    i++;

    LALInferenceSetVariable(runState->proposalArgs, "acceptanceCount", &(acceptanceCount));

    if (adaptationOn)
      LALInferenceAdaptation(runState, i);

    // Parallel tempering phase
    if (runPhase < 2) {

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
            if (runPhase<2 && iEff > floor((REAL8)(aclCheckCounter*endOfPhase)/(REAL8)numACLchecks)) {
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
          if (runPhase==0) {
            fprintf(stdout,"Chain %i has %i effective samples. Stopping...\n", MPIrank, iEff);
            runComplete = 1;          // Sampling is done!
          } else if (runPhase==1) {
            /* Broadcast new run phase to other chains */
            runPhase++;
            for (c=1; c<nChain; c++) {
              MPI_Send(&runPhase, 1, MPI_INT, c, 1, MPI_COMM_WORLD);
            }
          }
        }
      }
    }

    // Annealing phase
    if (runPhase==2) {
      if (phaseAcknowledged!=2) {
        phaseAcknowledged=2;

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
        runPhase++;
        phaseAcknowledged=runPhase;
        if (MPIrank==0)
          printf(" Single-chain sampling starting at iteration %i.\n", i);
      }
    } //if (runState==2)


    //Post-annealing single-chain sampling
    if (runPhase==3) {
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
    } //if (runPhase==3)

    runState->evolve(runState); //evolve the chain at temperature ladder[t]
    acceptanceCount = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "acceptanceCount");

    if (i==1){
      if (propVerbose) {
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
      if (diffEvo) {
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

      if (propVerbose){
        fprintf(propstatfile, "%d\t", i);
        LALInferencePrintProposalStats(propstatfile,runState->proposalStats);
        fflush(propstatfile);
      }
    }

    /* Propose parellel tempering swap */ 
    if (runPhase < 2) {
      if (MCMCMC) {
        /* Metropolis-coupled MCMC Swap (assumes likelihood function differs between chains).*/
        LALInferenceMCMCMCswap(runState, ladder, i, swapfile);
      } else {
        /* Standard parallel tempering swap. */
        LALInferencePTswap(runState, ladder, i, swapfile);
      }
    }// if (runPhase < 2)

    if (MPIrank==0 && i > Niter)
      runComplete=1;
    if (MPIrank==0 && runComplete==1) {
      for (c=1; c<nChain; c++) {
        MPI_Send(&runComplete, 1, MPI_INT, c, 2, MPI_COMM_WORLD);
      }
    }
  }// while (!runComplete)
  

  MPI_Barrier(MPI_COMM_WORLD);

  fclose(chainoutput);

  if(MPIrank == 0){
    if (adaptVerbose) {
      fclose(statfile);
    }
    if (tempVerbose) {
      fclose(swapfile);
    }
    if (propVerbose) {
      fclose(propstatfile);
    }
  }

  free(ladder);
  free(annealDecay);
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
    logLikelihoodProposed = runState->likelihood(&proposedParams, runState->data, runState->template);
  else
    logLikelihoodProposed = -DBL_MAX;

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

  /* special (clumsy) treatment for noise parameters */
  //???: Is there a better way to deal with accept/reject of noise parameters?
  if(LALInferenceCheckVariable(runState->currentParams,"psdscale"))
  {
    gsl_matrix *nx = *((gsl_matrix **)LALInferenceGetVariable(runState->currentParams, "psdstore"));
    gsl_matrix *ny = *((gsl_matrix **)LALInferenceGetVariable(runState->currentParams, "psdscale"));
    if(accepted == 1) gsl_matrix_memcpy(nx,ny);
    else              gsl_matrix_memcpy(ny,nx);
  }

  LALInferenceUpdateAdaptiveJumps(runState, accepted, targetAcceptance);
  LALInferenceDestroyVariables(&proposedParams);
}


//-----------------------------------------
// Swap routines:
//-----------------------------------------
void LALInferencePTswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile)
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
  INT4 swapAccepted=0;

  REAL8Vector * parameters = NULL;
  REAL8Vector * adjParameters = NULL;
  adjParameters = XLALCreateREAL8Vector(nPar);

  /* If Tskip reached, then block until next chain in ladder is prepared to accept swap proposal */
  if (((i % Tskip) == 0) && MPIrank < nChain-1) {
    /* Send current likelihood for swap proposal */
    MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);

    /* Determine if swap was accepted */
    MPI_Recv(&swapAccepted, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);

    /* Perform Swap */
    if (swapAccepted) {
      /* Set new likelihood */
      MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
      runState->currentLikelihood = adjCurrentLikelihood;

      /* Exchange current prior values */
      MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);
      MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
      runState->currentPrior = adjCurrentPrior;

      /* Package and send parameters */
      parameters = CopyLALInferenceVariablesToArray(runState->currentParams);
      MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);

      /* Recieve and unpack parameters */
      MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
      CopyArrayToLALInferenceVariables(adjParameters, runState->currentParams);
    }
  }

  /* Check if next lower temperature is ready to swap */
  if (MPIrank > 0) {
    MPI_Iprobe(MPIrank-1, 0, MPI_COMM_WORLD, &readyToSwap, &MPIstatus);

    /* Hotter chain decides acceptance */
    if (readyToSwap) {
      /* Receive adjacent likelilhood */
      MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);

      /* Determine if swap is accepted and tell the other chain */
      logChainSwap = (1.0/ladder[MPIrank-1]-1.0/ladder[MPIrank])*(runState->currentLikelihood-adjCurrentLikelihood);
      if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) {
        swapAccepted = 1;
      }

      MPI_Send(&swapAccepted, 1, MPI_INT, MPIrank-1, 0, MPI_COMM_WORLD);

      /* Perform Swap */
      if (swapAccepted) {
        /* Swap likelihoods */
        MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD);
        runState->currentLikelihood=adjCurrentLikelihood;

        /* Exchange current prior values */
        MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD);
        runState->currentPrior = adjCurrentPrior;

        /* Package parameters */
        parameters = CopyLALInferenceVariablesToArray(runState->currentParams);

        /* Swap parameters */
        MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD);

        /* Unpack parameters */
        CopyArrayToLALInferenceVariables(adjParameters, runState->currentParams);

        /* Print to file if verbose is chosen */
        if (swapfile != NULL) {
          fprintf(swapfile,"%d\t%f\t%f\t%f\n",i,logChainSwap,ladder[MPIrank-1],ladder[MPIrank]);
          fflush(swapfile);
        }
      }
    }
  }

  return;
}

void LALInferenceMCMCMCswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile)
{
  REAL8 nullLikelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
  INT4 MPIrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
  MPI_Status MPIstatus;
  INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
  INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
  INT4 Tskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "tempSkip");
  REAL8 adjCurrentPrior;
  REAL8 logChainSwap;
  INT4 readyToSwap = 0;
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
    /* Send current likelihood for swap proposal */
    MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);

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
    parameters = CopyLALInferenceVariablesToArray(runState->currentParams);
    MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);
    MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
    CopyArrayToLALInferenceVariables(adjParameters, adjCurrentParams);

    /* Calculate likelihood at adjacent parameters and send */
    lowLikeHighParams = runState->likelihood(adjCurrentParams, runState->data, runState->template);
    MPI_Send(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);

    /* Determine if swap was accepted */
    if(MPIrank==0){
      REAL8 flow = *(REAL8*) LALInferenceGetVariable(runState->currentParams, "fLow");
      fprintf(stdout,"%i:%f\n",MPIrank,flow);
    }
    MPI_Recv(&swapAccepted, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);

    if (swapAccepted) {
      /* Set new likelihood */
      runState->currentLikelihood = lowLikeHighParams;

      MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD);
      MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, 0, MPI_COMM_WORLD, &MPIstatus);
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
  if (MPIrank > 0) {
    MPI_Iprobe(MPIrank-1, 0, MPI_COMM_WORLD, &readyToSwap, &MPIstatus);

    /* Hotter chain decides acceptance */
    if (readyToSwap) {
      /* Receive adjacent likelilhood */
      MPI_Recv(&lowLikeLowParams, 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);
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
      parameters = CopyLALInferenceVariablesToArray(runState->currentParams);
      MPI_Recv(adjParameters->data, nPar, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);
      MPI_Send(parameters->data, nPar, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD);
      CopyArrayToLALInferenceVariables(adjParameters, adjCurrentParams);

      /* Calculate likelihood at adjacent parameters */
      highLikeLowParams = runState->likelihood(adjCurrentParams, runState->data, runState->template);

      /* Recieve likelihood from adjacent chain */
      MPI_Recv(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);

      /* Propose swap */
      if(MPIrank==1){
        REAL8 flow = *(REAL8*) LALInferenceGetVariable(runState->currentParams, "fLow");
        fprintf(stdout,"%i:%f\n",MPIrank,flow);
        fprintf(stdout,"ll: %f\nlh: %f\nhl %f\nhh:%f\n",lowLikeLowParams-nullLikelihood,lowLikeHighParams-nullLikelihood,highLikeLowParams-nullLikelihood,highLikeHighParams-nullLikelihood);
      }
      logChainSwap = (1./ladder[MPIrank-1])*(lowLikeHighParams-lowLikeLowParams)+(1./ladder[MPIrank])*(highLikeLowParams-highLikeHighParams);

      if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) {
        if(MPIrank==1) fprintf(stdout,"accepted\n\n");
        swapAccepted = 1;
      } else 
        if(MPIrank==1) fprintf(stdout,"rejected\n\n");
      MPI_Send(&swapAccepted, 1, MPI_INT, MPIrank-1, 0, MPI_COMM_WORLD);

      if (swapAccepted) {
        /* Swap likelihoods */
        runState->currentLikelihood = highLikeLowParams;

        /* Exchange current prior values (assumes prior functions are identical) */
        MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD, &MPIstatus);
        MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, 0, MPI_COMM_WORLD);
        runState->currentPrior = adjCurrentPrior;

        /* Set new parameters */
        LALInferenceCopyVariables(adjCurrentParams,runState->currentParams);

        /* Print to file if verbose is chosen */
        if (swapfile != NULL) {
          fprintf(swapfile,"%d\t%f\t%f\t%f\n",i,logChainSwap,ladder[MPIrank-1],ladder[MPIrank]);
          fflush(swapfile);
        }
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

  LALInferenceDestroyVariables(adjCurrentParams);
  XLALFree(adjCurrentParams);
  return;
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
  INT4 Niter = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "Niter");
  INT4 adapting=1;
  INT4 goodACL=0;

  for(LALInferenceVariableItem *item=runState->currentParams->head;item;item=item->next){
    if (item->vary != LALINFERENCE_PARAM_FIXED) {
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
FILE *LALInferencePrintPTMCMCHeader(LALInferenceRunState *runState) {
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

  chainoutput = fopen(outFileName,"w");
  if(chainoutput == NULL){
    XLALErrorHandler = XLALExitErrorHandler;
    XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
    XLAL_ERROR_NULL(XLAL_EIO);
  }
  
  LALInferencePrintPTMCMCHeaderFile(runState, chainoutput);

  fclose(chainoutput);

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
      LALInferenceDestroyVariables(saveParams);
      XLALFree(saveParams);
      XLAL_ERROR_VOID(XLAL_EINVAL, "unknown mass ratio parameter name (allowed are 'massratio' or 'asym_massratio')");
    }
    LALInferenceSetVariable(runState->currentParams, "time", &injGPSTime);
    LALInferenceSetVariable(runState->currentParams, "distance", &dist);
    LALInferenceSetVariable(runState->currentParams, "inclination", &inclination);
    LALInferenceSetVariable(runState->currentParams, "polarisation", &(psi));
    LALInferenceSetVariable(runState->currentParams, "phase", &phase);
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

    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    setIFOAcceptedLikelihoods(runState);
    LALInferencePrintPTMCMCHeaderFile(runState, out);
    fclose(out);
    
    LALInferenceCopyVariables(saveParams, runState->currentParams);
    runState->currentLikelihood = runState->likelihood(runState->currentParams, runState->data, runState->template);
    runState->currentPrior = runState->prior(runState, runState->currentParams);
    setIFOAcceptedLikelihoods(runState);    

    XLALFree(fname);
    LALInferenceDestroyVariables(saveParams);
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
