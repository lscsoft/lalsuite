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
  } if (LALInferenceCheckVariable(propArgs, cycleArrayName) && LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* There is one or the other of array or length in propArgs, which is an error. */
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
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
        
  /* The default, single-parameter updates. */
  LALInferenceAddProposalToCycle(runState, &LALInferenceSingleAdaptProposal, BIGWEIGHT);

  LALInferenceAddProposalToCycle(runState, &LALInferenceSkyLocWanderJump, SMALLWEIGHT);

  LALInferenceAddProposalToCycle(runState, &LALInferenceDrawUniformlyFromPrior, TINYWEIGHT);

  /* Now add various special proposals that are conditional on
     command-line arguments or variables in the params. */
  if(LALInferenceCheckVariable(proposedParams,"inclination")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceInclinationFlip, TINYWEIGHT);
  }

  if(LALInferenceCheckVariable(proposedParams,"phase")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceOrbitalPhaseJump, TINYWEIGHT);
  }

  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--covarianceMatrix");
  if (ppt) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceCovarianceEigenvectorJump, BIGWEIGHT);
  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--differential-evolution")) {
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionFull, BIGWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionMasses, TINYWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionAmp, TINYWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionSpins, TINYWEIGHT);
    LALInferenceAddProposalToCycle(runState, &LALInferenceDifferentialEvolutionSky, TINYWEIGHT);
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

  newRA = RA + jumpX/cos(DEC);
  newDEC = DEC + jumpY;

  LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
  LALInferenceSetVariable(proposedParams, "declination", &newDEC);

  LALInferenceSetLogProposalRatio(runState, log(cos(DEC)/cos(newDEC)));
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

  if (dePts == NULL || nPts <= 1) {
    fprintf(stderr, "Trying to differentially evolve without enough points (in %s, %d)\n",
            __FILE__, __LINE__);
    exit(1);
  }

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

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
  const char *names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionSpins(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"a_spin1", "a_spin2", "phi_spin1", "phi_spin2", "theta_spin1", "theta_spin2", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionSky(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"rightascension", "declination", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

/*draws a value from the prior, uniformly in individual parameters used for jumps.*/
void LALInferenceDrawFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {

  REAL8 value=0, min=0, max=0;
  //printf("%s\n",runState->currentParams->head->name);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  LALInferenceVariableItem *item=proposedParams->head;
  
  do {
    item=proposedParams->head;
    for(;item;item=item->next){
      if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
        continue;
      else
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, item->name, (void *)&min, (void *)&max);
          value=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
          LALInferenceSetVariable(proposedParams, item->name, &(value));
          //printf("%s\t%f\t%f\t%f\t%f\n",item->name, *(REAL8 *)item->value, value, min, max);
        }
    }
    //LALInferencePrintVariables(proposedParams);
    //printf("%f\n",runState->prior(runState, proposedParams));
  } while(runState->prior(runState, proposedParams)<=-DBL_MAX);
}

/*draws a value from the prior, using Von Neumann rejection sampling.*/
void LALInferenceDrawUniformlyFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  
  REAL8 value=0, min=0, max=0, b=0.01, alpha=0;
  //printf("%s\n",runState->currentParams->head->name);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  LALInferenceVariableItem *item=proposedParams->head;
  if(LALInferenceCheckVariable(runState->priorArgs, "densityVNR")){
    b = *((REAL8 *)LALInferenceGetVariable(runState->priorArgs, "densityVNR"));
  }
  
  do {
    item=proposedParams->head;
    for(;item;item=item->next){
      if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
        continue;
      else
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, item->name, (void *)&min, (void *)&max);
          value=min+(max-min)*gsl_rng_uniform(runState->GSLrandom);
          LALInferenceSetVariable(proposedParams, item->name, &(value));
          //printf("%s\t%f\t%f\t%f\t%f\n",item->name, *(REAL8 *)item->value, value, min, max);
        }
    }
    alpha=gsl_rng_uniform(runState->GSLrandom);
    //LALInferencePrintVariables(proposedParams);
  } while(exp(runState->prior(runState, proposedParams))<=alpha*b);
  //printf("%f\t%f\t%f\n",exp(runState->prior(runState, proposedParams)),alpha*b,b);

  REAL8 logRatio = runState->prior(runState, runState->currentParams) - runState->prior(runState, proposedParams);
  LALInferenceSetLogProposalRatio(runState, logRatio);
}
