/* 
 *  LALInferenceProposal.c:  Bayesian Followup, jump proposals.
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch, Will M. Farr
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
#include <lal/SkyCoordinates.h>
#include <lal/LALInference.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/XLALError.h>

#include <lal/LALStdlib.h>

static const char *cycleArrayName = "Proposal Cycle";
static const char *cycleArrayLengthName = "Proposal Cycle Length";
static const char *cycleArrayCounterName = "Proposal Cycle Counter";

const char *LALInferenceSigmaJumpName = "sigmaJump";

/* Mode hopping fraction for the differential evoultion proposals. */
static const REAL8 modeHoppingFrac = 0.1;

static int
same_detector_location(LALInferenceIFOData *d1, LALInferenceIFOData *d2) {
  UINT4 i;

  for (i = 0; i < 3; i++) {
    if (d1->detector->location[i] != d2->detector->location[i]) return 0;
  }

  return 1;
}

static UINT4 
numDetectorsUniquePositions(LALInferenceRunState *runState) {
  UINT4 nIFO = 0;
  UINT4 nCollision = 0;
  LALInferenceIFOData *currentIFO = NULL;

  for (currentIFO = runState->data; currentIFO; currentIFO = currentIFO->next) {
    LALInferenceIFOData *subsequentIFO = NULL;
    nIFO++;
    for (subsequentIFO = currentIFO->next; subsequentIFO; subsequentIFO = subsequentIFO->next) {
      if (same_detector_location(subsequentIFO, currentIFO)) {
        nCollision++;
        break;
      }
    }
  }

  return nIFO - nCollision;
}


static void
LALInferenceSetLogProposalRatio(LALInferenceRunState *runState, REAL8 logP) {
  if (LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio")) {
    LALInferenceSetVariable(runState->proposalArgs, "logProposalRatio", &logP);
  } else {
    LALInferenceAddVariable(runState->proposalArgs, "logProposalRatio", &logP, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
  }
}

void
LALInferenceAddProposalToCycle(LALInferenceRunState *runState, LALInferenceProposalFunction *prop, UINT4 weight) {
  const char *fname = "LALInferenceAddProposalToCycle";

  UINT4 length = 0;
  LALInferenceProposalFunction **cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;

  /* Quit without doing anything if weight = 0. */
  if (weight == 0) {
    return;
  }

  if (LALInferenceCheckVariable(propArgs, cycleArrayName) && LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* Have all the data in proposal args. */
    UINT4 i;

    length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));
    cycle = *((LALInferenceProposalFunction ***)LALInferenceGetVariable(propArgs, cycleArrayName));

    cycle = XLALRealloc(cycle, (length+weight)*sizeof(LALInferenceProposalFunction *));
    if (cycle == NULL) {
      XLALError(fname, __FILE__, __LINE__, XLAL_ENOMEM);
      exit(1);
    }

    for (i = length; i < length + weight; i++) {
      cycle[i] = prop;
    }

    length += weight;

    LALInferenceSetVariable(propArgs, cycleArrayLengthName, &length);
    LALInferenceSetVariable(propArgs, cycleArrayName, &cycle);
  } else {
    /* There are no data in proposal args.  Set some. */
    UINT4 i;
    
    length = weight;
    cycle = XLALMalloc(length*sizeof(LALInferenceProposalFunction *));
    if (cycle == NULL) {
      XLALError(fname, __FILE__, __LINE__, XLAL_ENOMEM);
      exit(1);
    }

    for (i = 0; i < length; i++) {
      cycle[i] = prop;
    }

    LALInferenceAddVariable(propArgs, cycleArrayLengthName, &length, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(propArgs, cycleArrayName, &cycle, LALINFERENCE_PROPOSAL_ARRAY_t, LALINFERENCE_PARAM_LINEAR);
  }
}

void
LALInferenceRandomizeProposalCycle(LALInferenceRunState *runState) {
  const char *fname = "LALInferenceRandomizeProposalCycle";
  UINT4 length = 0;
  LALInferenceProposalFunction **cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;

  UINT4 i;

  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  }

  cycle = *((LALInferenceProposalFunction ***)LALInferenceGetVariable(propArgs, cycleArrayName));
  length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));

  for (i = length - 1; i > 0; i--) {
    /* Fill in array from right to left, chosen randomly from remaining proposals. */
    UINT4 j;
    LALInferenceProposalFunction *prop;

    j = gsl_rng_uniform_int(runState->GSLrandom, i+1);
    prop = cycle[j];
    cycle[j] = cycle[i];
    cycle[i] = prop;
  }
}

void NSFillMCMCVariables(LALInferenceVariables *proposedParams)
{
  REAL8 distance=0.0,mc=0.0;
  if(LALInferenceCheckVariable(proposedParams,"logdistance"))
  {
    distance=exp(*(REAL8*)LALInferenceGetVariable(proposedParams,"logdistance"));
    LALInferenceAddVariable(proposedParams,"distance",&distance,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  if(LALInferenceCheckVariable(proposedParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(proposedParams,"logmc"));
    LALInferenceAddVariable(proposedParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  return;
}

void NSWrapMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{ /* PTMCMCLALProposal needs a few params converted */
 
  /* PTMCMC likes to read this directly so we have to plug our mangled values in*/
  LALInferenceVariables *currentParamsBackup=runState->currentParams;
  NSFillMCMCVariables(proposedParams);

  runState->currentParams=proposedParams; 
  LALInferenceDefaultProposal(runState,proposedParams);
  /* Restore currentParams */
  runState->currentParams=currentParamsBackup;
}

void 
LALInferenceCyclicProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *fname = "LALInferenceCyclicProposal";
  UINT4 length = 0;
  UINT4 i = 0;
  LALInferenceProposalFunction **cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;

  /* Must have cycle array and cycle array length in propArgs. */
  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  }

  length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));
  cycle = *((LALInferenceProposalFunction ***)LALInferenceGetVariable(propArgs, cycleArrayName));

  /* If there is not a proposal counter, put one into the variables, initialized to zero. */
  if (!LALInferenceCheckVariable(propArgs, cycleArrayCounterName)) {
    i = 0;
    LALInferenceAddVariable(propArgs, cycleArrayCounterName, &i, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_CIRCULAR);
  }

  i = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayCounterName));

  if (i >= length) {
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  }

  /* Call proposal. */
  (cycle[i])(runState, proposedParams);

  /* Increment counter for the next time around. */
  i = (i+1) % length;
  LALInferenceSetVariable(propArgs, cycleArrayCounterName, &i);
}

void
LALInferenceDeleteProposalCycle(LALInferenceRunState *runState) {
  LALInferenceVariables *propArgs = runState->proposalArgs;
  
  if (LALInferenceCheckVariable(propArgs, cycleArrayName)) {
    LALInferenceProposalFunction **cycle = *((LALInferenceProposalFunction ***)LALInferenceGetVariable(propArgs, cycleArrayName));
    XLALFree(cycle);
    LALInferenceRemoveVariable(propArgs, cycleArrayName);
  }

  if (LALInferenceCheckVariable(propArgs, cycleArrayCounterName)) {
    LALInferenceRemoveVariable(propArgs, cycleArrayCounterName);
  }

  if (LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    LALInferenceRemoveVariable(propArgs, cycleArrayLengthName);
  }
}

static void
SetupDefaultProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const UINT4 BIGWEIGHT = 20;
  const UINT4 SMALLWEIGHT = 5;
  const UINT4 TINYWEIGHT = 1;

  ProcessParamsTable *ppt;
        
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  /* The default, single-parameter updates. */
  LALInferenceAddProposalToCycle(runState, &LALInferenceSingleAdaptProposal, BIGWEIGHT);

  LALInferenceAddProposalToCycle(runState, &LALInferenceSkyLocWanderJump, SMALLWEIGHT);

  UINT4 nDet = numDetectorsUniquePositions(runState);
  if (nDet == 3) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceSkyReflectDetPlane, TINYWEIGHT);
  }

  LALInferenceAddProposalToCycle(runState, &LALInferenceDrawApproxPrior, TINYWEIGHT);

  /* Now add various special proposals that are conditional on
     command-line arguments or variables in the params. */
  if(LALInferenceCheckVariable(proposedParams,"inclination")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceInclinationFlip, TINYWEIGHT);
  }

  if(LALInferenceCheckVariable(proposedParams,"phase")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceOrbitalPhaseJump, TINYWEIGHT);
  }

  if (LALInferenceCheckVariable(proposedParams, "theta_spin1")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceRotateSpins, SMALLWEIGHT);
  }

  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--covarianceMatrix");
  if (ppt) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceCovarianceEigenvectorJump, BIGWEIGHT);
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--differential-evolution")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionFull, BIGWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionMasses, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionAmp, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionSpins, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionSky, SMALLWEIGHT);
  } 

  LALInferenceRandomizeProposalCycle(runState);
}

void LALInferenceDefaultProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  LALInferenceVariables *propArgs = runState->proposalArgs;

  /* If the cyclic proposal is not yet set up, set it up.  Note that
     this means that you can set up your own proposal cycle and it
     will be used in this function. */
  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* In case there is a partial cycle set up already, delete it. */
    LALInferenceDeleteProposalCycle(runState);
    SetupDefaultProposal(runState, proposedParams);
  }

  LALInferenceCyclicProposal(runState, proposedParams);
}

void LALInferenceSingleAdaptProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *args = runState->proposalArgs;
  ProcessParamsTable *ppt = LALInferenceGetProcParamVal(runState->commandLine, "--adapt");
  
  if (!LALInferenceCheckVariable(args, LALInferenceSigmaJumpName) || !ppt) {
    /* We are not adaptive, or for some reason don't have a sigma
       vector---fall back on old proposal. */
    LALInferenceSingleProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    LALInferenceVariableItem *param = NULL, *dummyParam = NULL;
    REAL8 T = *(REAL8 *)LALInferenceGetVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    UINT4 dim;
    UINT4 i;
    UINT4 varNr;
    REAL8Vector *sigmas = *(REAL8Vector **) LALInferenceGetVariable(args, LALInferenceSigmaJumpName);

    LALInferenceCopyVariables(runState->currentParams, proposedParams);

    dim = proposedParams->dimension;

    do {
      varNr = 1+gsl_rng_uniform_int(rng, dim);
      param = LALInferenceGetItemNr(proposedParams, varNr);
    } while (param->vary == LALINFERENCE_PARAM_FIXED || param->vary == LALINFERENCE_PARAM_OUTPUT);

    for (dummyParam = proposedParams->head, i = 0; dummyParam != NULL; dummyParam = dummyParam->next) {
      if (!strcmp(dummyParam->name, param->name)) {
        /* Found it; i = index into sigma vector. */
        break;
      } else if (dummyParam->vary == LALINFERENCE_PARAM_FIXED || dummyParam->vary == LALINFERENCE_PARAM_OUTPUT) {
        /* Don't increment i, since we're not dealing with a "real" parameter. */
        continue;
      } else {
        i++;
        continue;
      }
    }

    if (param->type != LALINFERENCE_REAL8_t) {
      fprintf(stderr, "Attempting to set non-REAL8 parameter with numerical sigma (in %s, %d)\n",
              __FILE__, __LINE__);
      exit(1);
    } 

    if (i >= sigmas->length) {
      fprintf(stderr, "Attempting to draw single-parameter jump %d past the end of sigma array %d.\n(Maybe you used a non-spinning correlation matrix for a spinning run?)\nError in %s, line %d.\n",
              i,sigmas->length,__FILE__, __LINE__);
      exit(1);
    }

    *((REAL8 *)param->value) += gsl_ran_ugaussian(rng)*sigmas->data[i]*sqrtT;

    LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);

    /* Set the log of the proposal ratio to zero, since this is a
       symmetric proposal. */
    LALInferenceSetLogProposalRatio(runState, 0.0);

    INT4 as = 1;
    LALInferenceSetVariable(args, "adaptableStep", &as);

    LALInferenceSetVariable(args, "proposedVariableNumber", &varNr);
    
    LALInferenceSetVariable(args, "proposedArrayNumber", &i);
  }
}

void LALInferenceSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  gsl_rng * GSLrandom=runState->GSLrandom;
  LALInferenceVariableItem *param=NULL, *dummyParam=NULL;
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 T = *(REAL8 *)LALInferenceGetVariable(runState->proposalArgs, "temperature");
	
  REAL8 sigma = 0.1*sqrt(T); /* Adapt step to temperature. */
  REAL8 big_sigma = 1.0;
  UINT4 dim;
  UINT4 i;
  UINT4 varNr;
  
  if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in a parameter
  if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in a parameter

  dim = proposedParams->dimension;
  
  do {
    varNr = 1+gsl_rng_uniform_int(GSLrandom, dim);
    param = LALInferenceGetItemNr(proposedParams, varNr);
  } while (param->vary == LALINFERENCE_PARAM_FIXED || param->vary == LALINFERENCE_PARAM_OUTPUT);
  
  for (dummyParam = proposedParams->head, i = 0; dummyParam != NULL; dummyParam = dummyParam->next) {
    if (!strcmp(dummyParam->name, param->name)) {
      /* Found it; i = index into sigma vector. */
      break;
    } else if (dummyParam->vary == LALINFERENCE_PARAM_FIXED || dummyParam->vary == LALINFERENCE_PARAM_OUTPUT) {
      /* Don't increment i, since we're not dealing with a "real" parameter. */
      continue;
    } else {
      i++;
      continue;
    }
  }	//printf("%s\n",param->name);
		
  if (LALInferenceGetProcParamVal(runState->commandLine, "--zeroLogLike")) {
    if (!strcmp(param->name, "massratio")) {
      sigma = 0.02;
    } else if (!strcmp(param->name, "chirpmass")) {
      sigma = 1.0;
    } else if (!strcmp(param->name, "time")) {
      sigma = 0.02;
    } else if (!strcmp(param->name, "phase")) {
      sigma = 0.6;
    } else if (!strcmp(param->name, "distance")) {
      sigma = 10.0;
    } else if (!strcmp(param->name, "declination")) {
      sigma = 0.3;
    } else if (!strcmp(param->name, "rightascension")) {
      sigma = 0.6;
    } else if (!strcmp(param->name, "polarisation")) {
      sigma = 0.6;
    } else if (!strcmp(param->name, "inclination")) {
      sigma = 0.3;
    } else if (!strcmp(param->name, "a_spin1")) {
      sigma = 0.1;
    } else if (!strcmp(param->name, "theta_spin1")) {
      sigma = 0.3;
    } else if (!strcmp(param->name, "phi_spin1")) {
      sigma = 0.6;
    } else if (!strcmp(param->name, "a_spin2")) {
      sigma = 0.1;
    } else if (!strcmp(param->name, "theta_spin2")) {
      sigma = 0.3;
    } else if (!strcmp(param->name, "phi_spin2")) {
      sigma = 0.6;
    } else {
      fprintf(stderr, "Could not find parameter %s!", param->name);
      exit(1);
    }
    *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*sigma;
  } else {
    if (!strcmp(param->name,"massratio") || !strcmp(param->name,"time") || !strcmp(param->name,"a_spin2") || !strcmp(param->name,"a_spin1")){
      *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
    } else if (!strcmp(param->name,"polarisation") || !strcmp(param->name,"phase") || !strcmp(param->name,"inclination")){
      *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
    } else {
      *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
    }
  }
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
  
  /* Symmetric Proposal. */
  LALInferenceSetLogProposalRatio(runState, 0.0);

  INT4 as = 1;
  LALInferenceSetVariable(runState->proposalArgs, "adaptableStep", &as);
  
  LALInferenceSetVariable(runState->proposalArgs, "proposedVariableNumber", &varNr);
  
  LALInferenceSetVariable(runState->proposalArgs, "proposedArrayNumber", &i);
  
}

void LALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  REAL8 phi;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  
  phi = *((REAL8 *) LALInferenceGetVariable(proposedParams, "phase"));

  phi = fmod(phi+M_PI, 2.0*M_PI);

  LALInferenceSetVariable(proposedParams, "phase", &phi);

  LALInferenceSetLogProposalRatio(runState, 0.0);

  /* Probably not needed, but play it safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void LALInferenceInclinationFlip(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  REAL8 iota;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  iota = *((REAL8 *) LALInferenceGetVariable(proposedParams, "inclination"));
  
  iota = M_PI - iota;

  LALInferenceSetVariable(proposedParams, "inclination", &iota);

  /* Symmetric proposal, because Cos(I) = Cos(Pi - I). */
  LALInferenceSetLogProposalRatio(runState, 0.0);

  /* Not really needed (probably), but let's be safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs); 
}

void LALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *proposalArgs = runState->proposalArgs;
  gsl_matrix *eigenvectors = *((gsl_matrix **)LALInferenceGetVariable(proposalArgs, "covarianceEigenvectors"));
  REAL8Vector *eigenvalues = *((REAL8Vector **)LALInferenceGetVariable(proposalArgs, "covarianceEigenvalues"));
  REAL8 temp = *((REAL8 *)LALInferenceGetVariable(proposalArgs, "temperature"));
  UINT4 N = eigenvalues->length;
  gsl_rng *rng = runState->GSLrandom;
  UINT4 i = gsl_rng_uniform_int(rng, N);
  REAL8 jumpSize = sqrt(temp*eigenvalues->data[i])*gsl_ran_ugaussian(rng);
  UINT4 j;
  LALInferenceVariableItem *proposeIterator;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  j = 0;
  proposeIterator = proposedParams->head;
  if (proposeIterator == NULL) {
    fprintf(stderr, "Bad proposed params in %s, line %d\n",
            __FILE__, __LINE__);
    exit(1);
  }
  do {
    if (proposeIterator->vary != LALINFERENCE_PARAM_FIXED && proposeIterator->vary != LALINFERENCE_PARAM_OUTPUT) {
      REAL8 tmp = *((REAL8 *)proposeIterator->value);
      REAL8 inc = jumpSize*gsl_matrix_get(eigenvectors, j, i);

      tmp += inc;

      memcpy(proposeIterator->value, &tmp, sizeof(REAL8));

      j++;
    }
  } while ((proposeIterator = proposeIterator->next) != NULL && j < N);

  LALInferenceSetLogProposalRatio(runState, 0.0);
}

void LALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  gsl_rng *rng = runState->GSLrandom;
  LALInferenceVariables *proposalArgs = runState->proposalArgs;
  REAL8 temp = *((REAL8 *)LALInferenceGetVariable(proposalArgs, "temperature"));
  REAL8 one_deg = 1.0 / (2.0*M_PI);
  REAL8 sigma = sqrt(temp)*one_deg;
  REAL8 XU = gsl_ran_ugaussian(rng);
  REAL8 YU = gsl_ran_ugaussian(rng);
  REAL8 jumpX = sigma*XU;
  REAL8 jumpY = sigma*YU;
  REAL8 RA, DEC;
  REAL8 newRA, newDEC;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  RA = *((REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension"));
  DEC = *((REAL8 *)LALInferenceGetVariable(proposedParams, "declination"));

  newRA = RA + jumpX;
  newDEC = DEC + jumpY;

  LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
  LALInferenceSetVariable(proposedParams, "declination", &newDEC);

  LALInferenceSetLogProposalRatio(runState, 0.0);
}

void LALInferenceDifferentialEvolutionFull(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceDifferentialEvolutionNames(runState, proposedParams, NULL);
}

void LALInferenceDifferentialEvolutionNames(LALInferenceRunState *runState, 
                                            LALInferenceVariables *proposedParams,
                                            const char **names) {
  if (names == NULL) {
    size_t i;
    size_t N = LALInferenceGetVariableDimension(runState->currentParams) + 1;
    names = alloca(N*sizeof(char *)); /* Hope we have alloca---saves
                                         having to deallocate after
                                         proposal. */
    for (i = 1; i <= N-1; i++) {
      names[i-1] = LALInferenceGetVariableName(runState->currentParams, i);
    }

    names[N-1]=NULL; /* Terminate */
  }


  LALInferenceVariables **dePts = runState->differentialPoints;
  size_t nPts = runState->differentialPointsLength;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  if (dePts == NULL || nPts <= 1) {
    return; /* Quit now, since we don't have any points to use. */
  }

  size_t i,j;

  i = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  do {
    j = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  } while (j == i);

  LALInferenceVariables *ptI = dePts[i];
  LALInferenceVariables *ptJ = dePts[j];
  REAL8 scale;

  /* Some small fraction of the time, we do a "mode hopping" jump,
     where we jump exactly along the difference vector. */
  if (gsl_rng_uniform(runState->GSLrandom) < modeHoppingFrac) {
    scale = 1.0;
  } else {      
    scale = 1.66511*gsl_ran_ugaussian(runState->GSLrandom); 
  }

  for (i = 0; names[i] != NULL; i++) {
    if (!LALInferenceCheckVariable(proposedParams, names[i]) || !LALInferenceCheckVariable(ptJ, names[i]) || !LALInferenceCheckVariable(ptI, names[i])) {
      /* Ignore variable if it's not in each of the params. */
    } else {
      REAL8 x = *((REAL8 *)LALInferenceGetVariable(proposedParams, names[i]));
      x += scale * (*((REAL8 *) LALInferenceGetVariable(ptJ, names[i])));
      x -= scale * (*((REAL8 *) LALInferenceGetVariable(ptI, names[i])));
      
      LALInferenceSetVariable(proposedParams, names[i], &x);
    }
  }
  
  LALInferenceSetLogProposalRatio(runState, 0.0); /* Symmetric proposal. */
}

void LALInferenceDifferentialEvolutionMasses(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"chirpmass", "massratio", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionAmp(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", "time", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionSpins(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"a_spin1", "a_spin2", "phi_spin1", "phi_spin2", "theta_spin1", "theta_spin2", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionSky(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"rightascension", "declination", "time", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

static REAL8 
draw_distance(LALInferenceRunState *runState) {
  REAL8 dmin, dmax;

  LALInferenceGetMinMaxPrior(runState->priorArgs, "distance", &dmin, &dmax);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);

  return pow(x*(dmax*dmax*dmax - dmin*dmin*dmin) + dmin*dmin*dmin, 1.0/3.0);
}

static REAL8 
draw_colatitude(LALInferenceRunState *runState, const char *name) {
  REAL8 min, max;

  LALInferenceGetMinMaxPrior(runState->priorArgs, name, &min, &max);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);

  return acos(cos(min) - x*(cos(min) - cos(max)));
}

static REAL8 
draw_dec(LALInferenceRunState *runState) {
  REAL8 min, max;
  
  LALInferenceGetMinMaxPrior(runState->priorArgs, "declination", &min, &max);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);
  
  return asin(x*(sin(max) - sin(min)) + sin(min));
}

static REAL8 
draw_flat(LALInferenceRunState *runState, const char *name) {
  REAL8 min, max;

  LALInferenceGetMinMaxPrior(runState->priorArgs, name, &min, &max);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);

  return min + x*(max - min);
}

static REAL8 
draw_chirp(LALInferenceRunState *runState) {
  REAL8 min, max;

  LALInferenceGetMinMaxPrior(runState->priorArgs, "chirpmass", &min, &max);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);

  return pow(pow(min, -5.0/6.0) - x*(pow(min, -5.0/6.0) - pow(max, -5.0/6.0)), -6.0/5.0);
}

static REAL8
approxLogPrior(LALInferenceVariables *params) {
  REAL8 logP = 0.0;

  REAL8 Mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
  logP += -11.0/6.0*log(Mc);

  /* Flat in eta. */

  REAL8 iota = *(REAL8 *)LALInferenceGetVariable(params, "inclination");
  logP += log(sin(iota));

  /* Flat in time, ra, psi, phi. */

  REAL8 dist = *(REAL8 *)LALInferenceGetVariable(params, "distance");
  logP += 2.0*log(dist);

  REAL8 dec = *(REAL8 *)LALInferenceGetVariable(params, "declination");
  logP += log(cos(dec));

  if (LALInferenceCheckVariable(params, "theta_spin1")) {
    REAL8 theta1 = *(REAL8 *)LALInferenceGetVariable(params, "theta_spin1");
    logP += log(sin(theta1));
  }

  if (LALInferenceCheckVariable(params, "theta_spin2")) {
    REAL8 theta2 = *(REAL8 *)LALInferenceGetVariable(params, "theta_spin2");
    logP += log(sin(theta2));
  }

  return logP;
}

void 
LALInferenceDrawApproxPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 Mc = draw_chirp(runState);
  LALInferenceSetVariable(proposedParams, "chirpmass", &Mc);

  REAL8 eta = draw_flat(runState, "massratio");
  LALInferenceSetVariable(proposedParams, "massratio", &eta);

  REAL8 theTime = draw_flat(runState, "time");
  LALInferenceSetVariable(proposedParams, "time", &theTime);

  REAL8 phase = draw_flat(runState, "phase");
  LALInferenceSetVariable(proposedParams, "phase", &phase);

  REAL8 inc = draw_colatitude(runState, "inclination");
  LALInferenceSetVariable(proposedParams, "inclination", &inc);

  REAL8 pol = draw_flat(runState, "polarisation");
  LALInferenceSetVariable(proposedParams, "polarisation", &pol);

  REAL8 dist = draw_distance(runState);
  LALInferenceSetVariable(proposedParams, "distance", &dist);

  REAL8 ra = draw_flat(runState, "rightascension");
  LALInferenceSetVariable(proposedParams, "rightascension", &ra);

  REAL8 dec = draw_dec(runState);
  LALInferenceSetVariable(proposedParams, "declination", &dec);

  if (LALInferenceCheckVariable(proposedParams, "a_spin1")) {
    REAL8 a1 = draw_flat(runState, "a_spin1");
    LALInferenceSetVariable(proposedParams, "a_spin1", &a1);
  }

  if (LALInferenceCheckVariable(proposedParams, "a_spin2")) {
    REAL8 a2 = draw_flat(runState, "a_spin2");
    LALInferenceSetVariable(proposedParams, "a_spin2", &a2);
  }

  if (LALInferenceCheckVariable(proposedParams, "phi_spin1")) {
    REAL8 phi1 = draw_flat(runState, "phi_spin1");
    LALInferenceSetVariable(proposedParams, "phi_spin1", &phi1);
  }

  if (LALInferenceCheckVariable(proposedParams, "phi_spin2")) {
    REAL8 phi2 = draw_flat(runState, "phi_spin2");
    LALInferenceSetVariable(proposedParams, "phi_spin2", &phi2);
  }

  if (LALInferenceCheckVariable(proposedParams, "theta_spin1")) {
    REAL8 theta1 = draw_colatitude(runState, "theta_spin1");
    LALInferenceSetVariable(proposedParams, "theta_spin1", &theta1);
  }

  if (LALInferenceCheckVariable(proposedParams, "theta_spin2")) {
    REAL8 theta2 = draw_colatitude(runState, "theta_spin2");
    LALInferenceSetVariable(proposedParams, "theta_spin2", &theta2);
  }

  LALInferenceSetLogProposalRatio(runState, approxLogPrior(runState->currentParams) - approxLogPrior(proposedParams));
}

static void
cross_product(REAL8 x[3], const REAL8 y[3], const REAL8 z[3]) {
  x[0] = y[1]*z[2]-y[2]*z[1];
  x[1] = y[2]*z[0]-y[0]*z[2];
  x[2] = y[0]*z[1]-y[1]*z[0];
}

static REAL8
norm(const REAL8 x[3]) {
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static void 
unit_vector(REAL8 v[3], const REAL8 w[3]) {
  REAL8 n = norm(w);

  if (n == 0.0) { 
    XLALError("unit_vector", __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  } else {
    v[0] = w[0] / n;
    v[1] = w[1] / n;
    v[2] = w[2] / n;
  }
}

static REAL8 
dot(const REAL8 v[3], const REAL8 w[3]) {
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

static void
project_along(REAL8 vproj[3], const REAL8 v[3], const REAL8 w[3]) {
  REAL8 what[3];
  REAL8 vdotw;

  unit_vector(what, w);
  vdotw = dot(v, w);

  vproj[0] = what[0]*vdotw;
  vproj[1] = what[1]*vdotw;
  vproj[2] = what[2]*vdotw;
}

static void
vsub(REAL8 diff[3], const REAL8 w[3], const REAL8 v[3]) {
  diff[0] = w[0] - v[0];
  diff[1] = w[1] - v[1];
  diff[2] = w[2] - v[2];
}

static void
vadd(REAL8 sum[3], const REAL8 w[3], const REAL8 v[3]) {
  sum[0] = w[0] + v[0];
  sum[1] = w[1] + v[1];
  sum[2] = w[2] + v[2];
}

static void
reflect_plane(REAL8 pref[3], const REAL8 p[3], 
              const REAL8 x[3], const REAL8 y[3], const REAL8 z[3]) {
  REAL8 n[3], nhat[3], xy[3], xz[3], pn[3], pnperp[3];

  vsub(xy, y, x);
  vsub(xz, z, x);

  cross_product(n, xy, xz);
  unit_vector(nhat, n);

  project_along(pn, p, nhat);
  vsub(pnperp, p, pn);

  vsub(pref, pnperp, pn);
}

static void 
sph_to_cart(REAL8 cart[3], const REAL8 lat, const REAL8 longi) {
  cart[0] = cos(longi)*sin(lat);
  cart[1] = sin(longi)*sin(lat);
  cart[2] = cos(lat);
}

static void
cart_to_sph(const REAL8 cart[3], REAL8 *lat, REAL8 *longi) {
  *longi = atan2(cart[1], cart[0]);
  *lat = acos(cart[2] / sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]));
}

static void
reflected_position_and_time(LALInferenceRunState *runState, const REAL8 ra, const REAL8 dec, const REAL8 oldTime,
                            REAL8 *newRA, REAL8 *newDec, REAL8 *newTime) {
  LALStatus status;
  SkyPosition currentEqu, currentGeo, newEqu, newGeo;
  currentEqu.latitude = dec;
  currentEqu.longitude = ra;
  currentEqu.system = COORDINATESYSTEM_EQUATORIAL;
  currentGeo.system = COORDINATESYSTEM_GEOGRAPHIC;
  LALEquatorialToGeographic(&status, &currentGeo, &currentEqu, &(runState->data->epoch));

  /* This function should only be called when we know that we have
     three detectors, or the following will crash. */
  REAL8 x[3], y[3], z[3];
  LALInferenceIFOData *xD = runState->data;
  memcpy(x, xD->detector->location, 3*sizeof(REAL8));

  LALInferenceIFOData *yD = xD->next;
  while (same_detector_location(yD, xD)) {
    yD = yD->next;
  }
  memcpy(y, yD->detector->location, 3*sizeof(REAL8));

  LALInferenceIFOData *zD = yD->next;
  while (same_detector_location(zD, yD) || same_detector_location(zD, xD)) {
    zD = zD->next;
  }
  memcpy(z, zD->detector->location, 3*sizeof(REAL8));

  REAL8 currentLoc[3];
  sph_to_cart(currentLoc, currentGeo.latitude, currentGeo.longitude);

  REAL8 newLoc[3];
  reflect_plane(newLoc, currentLoc, x, y, z);

  REAL8 newGeoLat, newGeoLongi;
  cart_to_sph(newLoc, &newGeoLat, &newGeoLongi);

  newGeo.latitude = newGeoLat;
  newGeo.longitude = newGeoLongi;
  newGeo.system = COORDINATESYSTEM_GEOGRAPHIC;
  newEqu.system = COORDINATESYSTEM_EQUATORIAL;
  LALGeographicToEquatorial(&status, &newEqu, &newGeo, &(runState->data->epoch));

  REAL8 oldDt, newDt;
  oldDt = XLALTimeDelayFromEarthCenter(runState->data->detector->location, currentEqu.longitude,
                                       currentEqu.latitude, &(runState->data->epoch));
  newDt = XLALTimeDelayFromEarthCenter(runState->data->detector->location, newEqu.longitude,
                                       newEqu.latitude, &(runState->data->epoch));

  *newRA = newEqu.longitude;
  *newDec = newEqu.latitude;
  *newTime = oldTime + oldDt - newDt;
}

void LALInferenceSkyReflectDetPlane(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  /* Find the number of distinct-position detectors. */
  /* Exit with same parameters (with a warning the first time) if
     there are not three detectors. */
  static UINT4 warningDelivered = 0;
  if (numDetectorsUniquePositions(runState) != 3) {
    if (warningDelivered) {
      /* Do nothing. */
    } else {
      fprintf(stderr, "WARNING: trying to reflect through the decector plane with %d\n", numDetectorsUniquePositions(runState));
      fprintf(stderr, "WARNING: geometrically independent locations,\n");
      fprintf(stderr, "WARNING: but this proposal should only be used with exactly 3 independent detectors.\n");
      fprintf(stderr, "WARNING: %s, line %d\n", __FILE__, __LINE__);
      warningDelivered = 1;
    }

    return; 
  }

  REAL8 ra = *(REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension");
  REAL8 dec = *(REAL8 *)LALInferenceGetVariable(proposedParams, "declination");
  REAL8 baryTime = *(REAL8 *)LALInferenceGetVariable(proposedParams, "time");

  REAL8 newRA, newDec, newTime;
  reflected_position_and_time(runState, ra, dec, baryTime, &newRA, &newDec, &newTime);

  /* Unit normal deviates, used to "fuzz" the state. */
  REAL8 nRA, nDec, nTime;
  const REAL8 epsTime = 6e-6; /* 1e-1 / (16 kHz) */
  const REAL8 epsAngle = 3e-4; /* epsTime*c/R_Earth */
  
  nRA = gsl_ran_ugaussian(runState->GSLrandom);
  nDec = gsl_ran_ugaussian(runState->GSLrandom);
  nTime = gsl_ran_ugaussian(runState->GSLrandom);

  newRA += epsAngle*nRA;
  newDec += epsAngle*nDec;
  newTime += epsTime*nTime;

  /* And the doubly-reflected position (near the original, but not
     exactly due to the fuzzing). */
  REAL8 refRA, refDec, refTime;
  reflected_position_and_time(runState, newRA, newDec, newTime, &refRA, &refDec, &refTime);

  /* The Gaussian increments required to shift us back to the original
     position from the doubly-reflected position. */
  REAL8 nRefRA, nRefDec, nRefTime;
  nRefRA = (ra - refRA)/epsAngle;
  nRefDec = (dec - refDec)/epsAngle;
  nRefTime = (baryTime - refTime)/epsTime;

  REAL8 pForward, pReverse;
  pForward = gsl_ran_ugaussian_pdf(nRA)*gsl_ran_ugaussian_pdf(nDec)*gsl_ran_ugaussian_pdf(nTime);
  pReverse = gsl_ran_ugaussian_pdf(nRefRA)*gsl_ran_ugaussian_pdf(nRefDec)*gsl_ran_ugaussian_pdf(nRefTime);

  LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
  LALInferenceSetVariable(proposedParams, "declination", &newDec);
  LALInferenceSetVariable(proposedParams, "time", &newTime);
  LALInferenceSetLogProposalRatio(runState, log(pReverse/pForward));
}

static void
rotateVectorAboutAxis(REAL8 vrot[3], REAL8 dvRot[3], 
                      const REAL8 v[3],
                      const REAL8 axis[3],
                      const REAL8 theta) {
  REAL8 vperp[3], vpar[3], vperprot[3];
  REAL8 xhat[3], yhat[3], zhat[3];
  REAL8 vp;
  UINT4 i;

  project_along(vpar, v, axis);
  vsub(vperp, v, vpar);

  vp = norm(vperp);

  unit_vector(zhat, axis);
  unit_vector(xhat, vperp);
  cross_product(yhat, zhat, xhat);

  for (i = 0; i < 3; i++) {
    vperprot[i] = vp*(cos(theta)*xhat[i] + sin(theta)*yhat[i]);
    dvRot[i] = vp*(cos(theta)*yhat[i] - sin(theta)*xhat[i]);
  }

  vadd(vrot, vpar, vperprot);
}

static void
vectorToColatLong(const REAL8 v[3], const REAL8 dv[3], 
                  REAL8 *colat, REAL8 *longi,
                  REAL8 *dcolat, REAL8 *dlongi) {
  REAL8 r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  REAL8 rho2 = v[0]*v[0] + v[1]*v[1];
  REAL8 r = sqrt(r2);
  REAL8 rho = sqrt(rho2);

  *longi = atan2(v[1], v[0]);
  if (*longi < 0.0) {
    *longi += 2.0*M_PI;
  }

  *colat = acos(v[2] / r);

  *dlongi = (v[0]*dv[1] - v[1]*dv[0]) / rho2;
  *dcolat = (v[2]*(v[0]*dv[0] + v[1]*dv[1]) - dv[2]*rho2) / (rho*r2);
}

void 
LALInferenceRotateSpins(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 theta1 = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);
  REAL8 theta2 = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);

  REAL8 logPr = 0.0;

  if (LALInferenceCheckVariable(proposedParams, "theta_spin1")) {
    REAL8 theta, phi, iota;
    REAL8 s1[3], L[3], newS[3], dS[3];
    
    theta = *(REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin1");
    phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin1");

    iota = *(REAL8 *)LALInferenceGetVariable(proposedParams, "inclination");

    s1[0] = cos(phi)*sin(theta);
    s1[1] = sin(phi)*sin(theta);
    s1[2] = cos(theta);

    L[0] = sin(iota);
    L[1] = 0.0;
    L[2] = cos(iota);

    rotateVectorAboutAxis(newS, dS, s1, L, theta1);

    REAL8 newPhi, newTheta, dPhi, dTheta;

    vectorToColatLong(newS, dS, &newTheta, &newPhi, &dTheta, &dPhi);

    REAL8 LxS[3], dOldTheta, dOldPhi;

    cross_product(LxS, L, s1);
    vectorToColatLong(s1, LxS, &theta, &phi, &dOldTheta, &dOldPhi);

    logPr += 0.5*log((dTheta*dTheta + dPhi*dPhi)/(dOldTheta*dOldTheta + dOldPhi*dOldPhi));

    LALInferenceSetVariable(proposedParams, "phi_spin1", &newPhi);
    LALInferenceSetVariable(proposedParams, "theta_spin1", &newTheta);
  }

  if (LALInferenceCheckVariable(proposedParams, "theta_spin2")) {
    REAL8 theta, phi, iota;
    REAL8 s2[3], L[3], newS[3], dS[3];
    
    theta = *(REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin2");
    phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin2");

    iota = *(REAL8 *)LALInferenceGetVariable(proposedParams, "inclination");

    s2[0] = cos(phi)*sin(theta);
    s2[1] = sin(phi)*sin(theta);
    s2[2] = cos(theta);

    L[0] = sin(iota);
    L[1] = 0.0;
    L[2] = cos(iota);

    rotateVectorAboutAxis(newS, dS, s2, L, theta2);

    REAL8 newPhi, newTheta, dPhi, dTheta;

    vectorToColatLong(newS, dS, &newTheta, &newPhi, &dTheta, &dPhi);

    REAL8 LxS[3], dOldTheta, dOldPhi;

    cross_product(LxS, L, s2);
    vectorToColatLong(s2, LxS, &theta, &phi, &dOldTheta, &dOldPhi);

    logPr += 0.5*log((dTheta*dTheta + dPhi*dPhi)/(dOldTheta*dOldTheta + dOldPhi*dOldPhi));

    LALInferenceSetVariable(proposedParams, "phi_spin2", &newPhi);
    LALInferenceSetVariable(proposedParams, "theta_spin2", &newTheta);
  }

  LALInferenceSetLogProposalRatio(runState, logPr);
}
