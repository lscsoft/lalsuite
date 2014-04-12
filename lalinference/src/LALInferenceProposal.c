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

#define LAL_USE_OLD_COMPLEX_STRUCTS

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
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimNoise.h>
#include <lal/XLALError.h>

#include <lal/LALStdlib.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

const char *const cycleArrayName = "Proposal Cycle";
const char *const cycleArrayLengthName = "Proposal Cycle Length";
const char *const cycleArrayCounterName = "Proposal Cycle Counter";

const char *const LALInferenceCurrentProposalName = "Current Proposal";

/* Proposal Names */
const char *const singleAdaptProposalName = "Single";
const char *const singleProposalName = "Single";
const char *const orbitalPhaseJumpName = "OrbitalPhase";
const char *const covarianceEigenvectorJumpName = "CovarianceEigenvector";
const char *const skyLocWanderJumpName = "SkyLocWander";
const char *const differentialEvolutionFullName = "DifferentialEvolutionFull";
const char *const differentialEvolutionIntrinsicName = "DifferentialEvolutionIntrinsic";
const char *const differentialEvolutionExtrinsicName = "DifferentialEvolutionExtrinsic";
const char *const drawApproxPriorName = "DrawApproxPrior";
const char *const skyReflectDetPlaneName = "SkyReflectDetPlane";
const char *const skyRingProposalName = "SkyRingProposal";
const char *const PSDFitJumpName = "PSDFitJump";
const char *const rotateSpinsName = "RotateSpins";
const char *const polarizationPhaseJumpName = "PolarizationPhase";
const char *const polarizationCorrPhaseJumpName = "CorrPolarizationPhase";
const char *const extrinsicParamProposalName = "ExtrinsicParamProposal";
const char *const KDNeighborhoodProposalName = "KDNeighborhood";
const char *const frequencyBinJumpName = "FrequencyBin";
const char *const GlitchMorletJumpName = "glitchMorletJump";
const char *const GlitchMorletReverseJumpName = "glitchMorletReverseJump";

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
LALInferenceAddProposalToCycle(LALInferenceRunState *runState, const char *propName, LALInferenceProposalFunction prop, UINT4 weight) {
  const char *fname = "LALInferenceAddProposalToCycle";

  UINT4 length = 0;
  LALInferenceProposalFunction *cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;
  LALInferenceVariables *propStats = runState->proposalStats;
  UINT4 i;

  /* Quit without doing anything if weight = 0. */
  if (weight == 0) {
    return;
  }

  if (LALInferenceCheckVariable(propArgs, cycleArrayName) && LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* Have all the data in proposal args. */

    length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));
    cycle = *((LALInferenceProposalFunction **)LALInferenceGetVariable(propArgs, cycleArrayName));
  }
  
  cycle = XLALRealloc(cycle, (length+weight)*sizeof(LALInferenceProposalFunction));
  if (cycle == NULL) {
    XLALError(fname, __FILE__, __LINE__, XLAL_ENOMEM);
    exit(1);
  }
  for (i = length; i < length + weight; i++) {
    cycle[i] = prop;
  }

  length += weight;

  LALInferenceAddVariable(propArgs, cycleArrayLengthName, &length, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(propArgs, cycleArrayName, (void *)&cycle, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_LINEAR);

  /* If propStats is not NULL, add counters for proposal function if they aren't already there */
  if(propStats){
    if(!LALInferenceCheckVariable(propStats, propName)){
      LALInferenceProposalStatistics propStat = {
        .weight = weight,
        .proposed = 0,
        .accepted = 0};
      LALInferenceAddVariable(propStats, propName, (void *)&propStat, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_LINEAR);
    }
  }
}

void
LALInferenceRandomizeProposalCycle(LALInferenceRunState *runState) {
  const char *fname = "LALInferenceRandomizeProposalCycle";
  UINT4 length = 0;
  LALInferenceProposalFunction *cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;

  UINT4 i;

  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  }

  cycle = *((LALInferenceProposalFunction **)LALInferenceGetVariable(propArgs, cycleArrayName));
  length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));

  for (i = length - 1; i > 0; i--) {
    /* Fill in array from right to left, chosen randomly from remaining proposals. */
    UINT4 j;
    LALInferenceProposalFunction prop;

    j = gsl_rng_uniform_int(runState->GSLrandom, i+1);
    prop = cycle[j];
    cycle[j] = cycle[i];
    cycle[i] = prop;
  }
}

/* Convert NS to MCMC variables (call before calling MCMC proposal from NS) */
void NSFillMCMCVariables(LALInferenceVariables *proposedParams, LALInferenceVariables *priorArgs)
{
  REAL8 distance=0.0,mc=0.0,dmin,dmax,mmin,mmax;
  if(LALInferenceCheckVariable(proposedParams,"logdistance"))
  {
    distance=exp(*(REAL8*)LALInferenceGetVariable(proposedParams,"logdistance"));
    LALInferenceAddVariable(proposedParams,"distance",&distance,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  if(!LALInferenceCheckMinMaxPrior(priorArgs,"distance") &&
     LALInferenceCheckMinMaxPrior(priorArgs,"logdistance"))
  {
    LALInferenceGetMinMaxPrior(priorArgs,"logdistance",&dmin,&dmax);
    dmin=exp(dmin); dmax=exp(dmax);
    LALInferenceAddMinMaxPrior(priorArgs,"distance",&dmin,&dmax,LALINFERENCE_REAL8_t);
  }
  if(LALInferenceCheckVariable(proposedParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(proposedParams,"logmc"));
    LALInferenceAddVariable(proposedParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  if(!LALInferenceCheckMinMaxPrior(priorArgs,"chirpmass") &&
     LALInferenceCheckMinMaxPrior(priorArgs,"logmc"))
  {
    LALInferenceGetMinMaxPrior(priorArgs,"logmc",&mmin,&mmax);
    mmin=exp(mmin); mmax=exp(mmax);
    LALInferenceAddMinMaxPrior(priorArgs,"chirpmass",&mmin,&mmax,LALINFERENCE_REAL8_t);
  }
  return;
}


void
LALInferenceCyclicProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *fname = "LALInferenceCyclicProposal";
  UINT4 length = 0;
  UINT4 i = 0;
  LALInferenceProposalFunction *cycle = NULL;
  LALInferenceVariables *propArgs = runState->proposalArgs;

  /* Must have cycle array and cycle array length in propArgs. */
  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    XLALError(fname, __FILE__, __LINE__, XLAL_FAILURE);
    exit(1);
  }

  length = *((UINT4 *)LALInferenceGetVariable(propArgs, cycleArrayLengthName));
  cycle = *((LALInferenceProposalFunction **)LALInferenceGetVariable(propArgs, cycleArrayName));

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
    LALInferenceProposalFunction *cycle = *((LALInferenceProposalFunction **)LALInferenceGetVariable(propArgs, cycleArrayName));
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

void LALInferenceSetupDefaultNSProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const UINT4 BIGWEIGHT = 20;
  const UINT4 SMALLWEIGHT = 5;
  const UINT4 TINYWEIGHT = 1;
  const char defaultPropName[]="none";
  UINT4 fullProp = 1;
  UINT4 nDet = numDetectorsUniquePositions(runState);

  if(!runState->proposalStats) runState->proposalStats = XLALCalloc(1,sizeof(LALInferenceVariables));

  if(!LALInferenceCheckVariable(runState->proposalArgs,LALInferenceCurrentProposalName))
      LALInferenceAddVariable(runState->proposalArgs,LALInferenceCurrentProposalName, (void*)&defaultPropName, LALINFERENCE_string_t, LALINFERENCE_PARAM_OUTPUT);

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  /* The default, single-parameter updates. */
  if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-singleadapt"))
  {
    LALInferenceSetupAdaptiveProposals(runState);
    LALInferenceAddProposalToCycle(runState, singleAdaptProposalName, &LALInferenceSingleAdaptProposal, TINYWEIGHT);
  }

  if(!LALInferenceGetProcParamVal(runState->commandLine,"--margphi") && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-psiphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi"))
    LALInferenceAddProposalToCycle(runState, polarizationPhaseJumpName, &LALInferencePolarizationPhaseJump, TINYWEIGHT);

  LALInferenceFrame frame=LALINFERENCE_FRAME_RADIATION;
  if (LALInferenceCheckVariable(runState->currentParams, "LALINFERENCE_FRAME"))
    frame = *(LALInferenceFrame*) LALInferenceGetVariable(runState->currentParams, "LALINFERENCE_FRAME");

  if (nDet >= 3 && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-extrinsicparam") && frame == LALINFERENCE_FRAME_RADIATION ) {
    LALInferenceAddProposalToCycle(runState, extrinsicParamProposalName, &LALInferenceExtrinsicParamProposal, SMALLWEIGHT);
  }

  if (fullProp) {
    if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-skywander"))
    {   /* If there are not 3 detectors, the other sky jumps are not used, so increase the % of wandering jumps */
        if(nDet<3) LALInferenceAddProposalToCycle(runState, skyLocWanderJumpName, &LALInferenceSkyLocWanderJump, BIGWEIGHT);
        else LALInferenceAddProposalToCycle(runState, skyLocWanderJumpName, &LALInferenceSkyLocWanderJump, 3.0*SMALLWEIGHT);
    }
    if (nDet >= 3 && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-skyreflect")) {
      LALInferenceAddProposalToCycle(runState, skyReflectDetPlaneName, &LALInferenceSkyReflectDetPlane, TINYWEIGHT);
    }
    if (nDet>=2 && !LALInferenceGetProcParamVal(runState->commandLine,"--noProposalSkyRing")) {
      LALInferenceAddProposalToCycle(runState, skyRingProposalName, &LALInferenceSkyRingProposal, SMALLWEIGHT);
    }
    if(LALInferenceGetProcParamVal(runState->commandLine,"--proposal-drawprior"))
      LALInferenceAddProposalToCycle(runState, drawApproxPriorName, &LALInferenceDrawApproxPrior, TINYWEIGHT);

    if(!LALInferenceGetProcParamVal(runState->commandLine,"--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi"))
        if(LALInferenceCheckVariableNonFixed(proposedParams,"phase")) {
          LALInferenceAddProposalToCycle(runState, orbitalPhaseJumpName, &LALInferenceOrbitalPhaseJump, TINYWEIGHT);
          if (!LALInferenceGetProcParamVal(runState->commandLine,"--noProposalCorrPsiPhi"))
            LALInferenceAddProposalToCycle(runState, polarizationCorrPhaseJumpName, &LALInferenceCorrPolarizationPhaseJump, SMALLWEIGHT);
        }
  }

  /* Now add various special proposals that are conditional on
     command-line arguments or variables in the params. */

  if (LALInferenceCheckVariable(proposedParams, "theta_spin1")&&!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-rotate-spins")) {
  	if(LALInferenceGetVariableVaryType(proposedParams,"theta_spin1")==LALINFERENCE_PARAM_CIRCULAR
  		|| LALInferenceGetVariableVaryType(proposedParams,"theta_spin1")==LALINFERENCE_PARAM_LINEAR )
	    LALInferenceAddProposalToCycle(runState, rotateSpinsName, &LALInferenceRotateSpins, SMALLWEIGHT);
  }

  /* Always use the covariance method */
    LALInferenceAddProposalToCycle(runState, covarianceEigenvectorJumpName, &LALInferenceCovarianceEigenvectorJump, BIGWEIGHT);

  /* Use differential evolution unless turned off */
  if (!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-differentialevolution")) {
    LALInferenceAddProposalToCycle(runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull, BIGWEIGHT);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionIntrinsicName, &LALInferenceDifferentialEvolutionIntrinsic, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionExtrinsicName, &LALInferenceDifferentialEvolutionExtrinsic, SMALLWEIGHT);
  }

  //Add LALInferencePSDFitJump to the cycle
  if(LALInferenceGetProcParamVal(runState->commandLine, "--psdFit"))
  {
    LALInferenceAddProposalToCycle (runState, PSDFitJumpName, *LALInferencePSDFitJump, SMALLWEIGHT);
  }



  /********** TURNED OFF - very small acceptance with nested sampling, slows everything down ****************/
  /*
  if (!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-kdtree")) {
    LALInferenceAddProposalToCycle(runState, KDNeighborhoodProposalName, &LALInferenceKDNeighborhoodProposal, SMALLWEIGHT);
  }
  */

  LALInferenceRandomizeProposalCycle(runState);
  LALInferenceZeroProposalStats(runState);
}



static void
SetupDefaultProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const UINT4 BIGWEIGHT = 20;
  const UINT4 SMALLWEIGHT = 5;
  const UINT4 TINYWEIGHT = 1;
  const char defaultPropName[]="none";
  UINT4 fullProp = 1;
  UINT4 nDet = numDetectorsUniquePositions(runState);

  /* If MCMC w/ parallel tempering, use reduced set of proposal functions for cold chains */
  if(LALInferenceCheckVariable(runState->proposalArgs, "hotChain")) {
    fullProp = *(UINT4 *)LALInferenceGetVariable(runState->proposalArgs, "hotChain");
  }

  ProcessParamsTable *ppt;
  if(!LALInferenceCheckVariable(runState->proposalArgs,LALInferenceCurrentProposalName))
      LALInferenceAddVariable(runState->proposalArgs,LALInferenceCurrentProposalName, (void*)&defaultPropName, LALINFERENCE_string_t, LALINFERENCE_PARAM_OUTPUT);

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  /* Only add signal proposals if signal is part of the model */
  if(!LALInferenceGetProcParamVal(runState->commandLine,"--noiseonly"))
  {
  /* The default, single-parameter updates. */
  if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-singleadapt"))
    LALInferenceAddProposalToCycle(runState, singleAdaptProposalName, &LALInferenceSingleAdaptProposal, BIGWEIGHT);

  if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-psiphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi")) {
    LALInferenceAddProposalToCycle(runState, polarizationPhaseJumpName, &LALInferencePolarizationPhaseJump, TINYWEIGHT);
  }

  LALInferenceFrame frame=LALINFERENCE_FRAME_RADIATION;
  if (LALInferenceCheckVariable(runState->currentParams, "LALINFERENCE_FRAME"))
    frame = *(LALInferenceFrame*) LALInferenceGetVariable(runState->currentParams, "LALINFERENCE_FRAME");

  if (nDet >= 3 && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-extrinsicparam") && frame == LALINFERENCE_FRAME_RADIATION ) {
    LALInferenceAddProposalToCycle(runState, extrinsicParamProposalName, &LALInferenceExtrinsicParamProposal, SMALLWEIGHT);
  }

  if (fullProp) {
    if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-skywander"))
      LALInferenceAddProposalToCycle(runState, skyLocWanderJumpName, &LALInferenceSkyLocWanderJump, SMALLWEIGHT);

    if (nDet == 3 && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-skyreflect") ) {
      LALInferenceAddProposalToCycle(runState, skyReflectDetPlaneName, &LALInferenceSkyReflectDetPlane, TINYWEIGHT);
    }

    if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-drawprior"))
      LALInferenceAddProposalToCycle(runState, drawApproxPriorName, &LALInferenceDrawApproxPrior, TINYWEIGHT);

    if(LALInferenceCheckVariableNonFixed(proposedParams,"phase") && !LALInferenceGetProcParamVal(runState->commandLine, "--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi")) {
      LALInferenceAddProposalToCycle(runState, orbitalPhaseJumpName, &LALInferenceOrbitalPhaseJump, TINYWEIGHT);
    }
  }

  /* Now add various special proposals that are conditional on
     command-line arguments or variables in the params. */

  if (LALInferenceCheckVariable(proposedParams, "theta_spin1")) {
    LALInferenceAddProposalToCycle(runState, rotateSpinsName, &LALInferenceRotateSpins, SMALLWEIGHT);
  }

  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--covariancematrix");
  if(!ppt){
    ppt=LALInferenceGetProcParamVal(runState->commandLine, "--covarianceMatrix");
    if(ppt) XLALPrintWarning("WARNING: Deprecated --covarianceMatrix option will be removed, please change to --covariancematrix");
  }
  if (ppt) {
    LALInferenceAddProposalToCycle(runState, covarianceEigenvectorJumpName, &LALInferenceCovarianceEigenvectorJump, BIGWEIGHT);
  }

  if (!LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution")
      && !LALInferenceGetProcParamVal(runState->commandLine, "--nodifferentialevolution") && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-differentialevolution")) {
    LALInferenceAddProposalToCycle(runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull, BIGWEIGHT);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionIntrinsicName, &LALInferenceDifferentialEvolutionIntrinsic, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionExtrinsicName, &LALInferenceDifferentialEvolutionExtrinsic, SMALLWEIGHT);

  }

  if (LALInferenceGetProcParamVal(runState->commandLine, "--kDTree") || LALInferenceGetProcParamVal(runState->commandLine,"--kdtree")) {
    LALInferenceAddProposalToCycle(runState, KDNeighborhoodProposalName, &LALInferenceKDNeighborhoodProposal, SMALLWEIGHT);
  }

  if (nDet >= 2 && !LALInferenceGetProcParamVal(runState->commandLine,"--noProposalSkyRing")) {
    LALInferenceAddProposalToCycle(runState, skyRingProposalName, &LALInferenceSkyRingProposal, SMALLWEIGHT);
  }

  if (!LALInferenceGetProcParamVal(runState->commandLine,"--noProposalCorrPsiPhi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi")) {
    LALInferenceAddProposalToCycle(runState, polarizationCorrPhaseJumpName, &LALInferenceCorrPolarizationPhaseJump, SMALLWEIGHT);
  }
  }//End noise-only conditional

  /* Add LALInferencePSDFitJump to the cycle */
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--psdFit");
  if(ppt)
  {
    LALInferenceAddProposalToCycle(runState, PSDFitJumpName, *LALInferencePSDFitJump, SMALLWEIGHT);
  }

  /* Add glitch-fitting proposals to cycle */
  ppt=LALInferenceGetProcParamVal(runState->commandLine, "--glitchFit");
  if(ppt)
  {
    //Morlet wavelet propposals
    LALInferenceAddProposalToCycle(runState, GlitchMorletJumpName, *LALInferenceGlitchMorletProposal, SMALLWEIGHT);
    LALInferenceAddProposalToCycle(runState, GlitchMorletReverseJumpName, *LALInferenceGlitchMorletReverseJump, SMALLWEIGHT);
  }

  LALInferenceRandomizeProposalCycle(runState);
  LALInferenceZeroProposalStats(runState);
}

static void
SetupRapidSkyLocProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  LALInferenceAddProposalToCycle(runState, singleAdaptProposalName, &LALInferenceSingleAdaptProposal, 100);
  //LALInferenceAddProposalToCycle(runState, skyLocWanderJumpName, &LALInferenceSkyLocWanderJump, 1);
  if (!LALInferenceGetProcParamVal(runState->commandLine, "--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi"))
    LALInferenceAddProposalToCycle(runState, polarizationPhaseJumpName, &LALInferencePolarizationPhaseJump, 1);

  LALInferenceFrame frame=LALINFERENCE_FRAME_RADIATION;
  if (LALInferenceCheckVariable(runState->currentParams, "LALINFERENCE_FRAME"))
    frame = *(LALInferenceFrame*) LALInferenceGetVariable(runState->currentParams, "LALINFERENCE_FRAME");

  UINT4 nDet = numDetectorsUniquePositions(runState);
  if (nDet >= 3 && !LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-extrinsicparam") && frame == LALINFERENCE_FRAME_RADIATION && !LALInferenceGetProcParamVal(runState->commandLine,"--margtime") && !LALInferenceGetProcParamVal(runState->commandLine,"--margtimephi")) {
    LALInferenceAddProposalToCycle(runState, extrinsicParamProposalName, &LALInferenceExtrinsicParamProposal, 20);
  }

  /*
  UINT4 nDet = numDetectorsUniquePositions(runState);
  if (nDet == 3) {
    LALInferenceAddProposalToCycle(runState, skyReflectDetPlaneName, &LALInferenceSkyReflectDetPlane, 1);
  }
  */

  if (!LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution")) {
    LALInferenceAddProposalToCycle(runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull, 10);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionExtrinsicName, &LALInferenceDifferentialEvolutionExtrinsic, 5);
  }

  LALInferenceRandomizeProposalCycle(runState);
}

static void
SetupPostPTProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  if (!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-singleadapt"))
    LALInferenceAddProposalToCycle(runState, singleAdaptProposalName, &LALInferenceSingleAdaptProposal, 5);

  if(!LALInferenceGetProcParamVal(runState->commandLine,"--proposal-no-psiphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margphi") && !LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi"))
    LALInferenceAddProposalToCycle(runState, polarizationPhaseJumpName, &LALInferencePolarizationPhaseJump, 1);

  if (!LALInferenceGetProcParamVal(runState->commandLine, "--noDifferentialEvolution")) {
    LALInferenceAddProposalToCycle(runState, differentialEvolutionFullName, &LALInferenceDifferentialEvolutionFull, 2);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionIntrinsicName, &LALInferenceDifferentialEvolutionIntrinsic, 4);
    LALInferenceAddProposalToCycle(runState, differentialEvolutionExtrinsicName, &LALInferenceDifferentialEvolutionExtrinsic, 4);
  }

  LALInferenceRandomizeProposalCycle(runState);
}

void LALInferencePostPTProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *propArgs = runState->proposalArgs;
  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* In case there is a partial cycle set up already, delete it. */
    LALInferenceDeleteProposalCycle(runState);
    SetupPostPTProposal(runState, proposedParams);
  }

  LALInferenceCyclicProposal(runState, proposedParams);
}


void LALInferenceRapidSkyLocProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *propArgs = runState->proposalArgs;

  if (!LALInferenceCheckVariable(propArgs, cycleArrayName) || !LALInferenceCheckVariable(propArgs, cycleArrayLengthName)) {
    /* In case there is a partial cycle set up already, delete it. */
    LALInferenceDeleteProposalCycle(runState);
    SetupRapidSkyLocProposal(runState, proposedParams);
  }

  LALInferenceCyclicProposal(runState, proposedParams);
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
  /* Set adapting flag to 0, it will be set to 1 if an adapting step is used*/
  LALInferenceCyclicProposal(runState, proposedParams);
}

void LALInferenceSingleAdaptProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = singleAdaptProposalName;
  LALInferenceVariables *args = runState->proposalArgs;
  LALInferenceSetVariable(args, LALInferenceCurrentProposalName, &propName);
  ProcessParamsTable *ppt = LALInferenceGetProcParamVal(runState->commandLine, "--noAdapt");
  gsl_matrix *m=NULL;
  UINT4Vector *v=NULL;

  if (ppt) {
    /* We are not adaptive, or for some reason don't have a sigma
       vector---fall back on old proposal. */
    LALInferenceSingleProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    LALInferenceVariableItem *param = NULL, *dummyParam = NULL;
    REAL8 T = 1.0;
    if(LALInferenceCheckVariable(args,"temperature"))
      T=*(REAL8 *)LALInferenceGetVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    UINT4 dim;
    UINT4 i;
    UINT4 varNr;
    char tmpname[MAX_STRLEN]="";

    LALInferenceCopyVariables(runState->currentParams, proposedParams);

    dim = proposedParams->dimension;

    do {
    varNr = 1+gsl_rng_uniform_int(rng, dim);
    param = LALInferenceGetItemNr(proposedParams, varNr);

    } while (param->vary == LALINFERENCE_PARAM_FIXED || param->vary == LALINFERENCE_PARAM_OUTPUT || param->type != LALINFERENCE_REAL8_t);

    for (dummyParam = proposedParams->head, i = 0; dummyParam != NULL; dummyParam = dummyParam->next) {
          if (!strcmp(dummyParam->name, param->name)) {
            /* Found it; i = index into sigma vector. */
            break;
          } else if (dummyParam->vary == LALINFERENCE_PARAM_FIXED || dummyParam->vary == LALINFERENCE_PARAM_OUTPUT) {
            /* Don't increment i, since we're not dealing with a "real" parameter. */
            continue;
          } else if (param->type == LALINFERENCE_gslMatrix_t)
          {
            /*increment i by number of noise parameters, since they aren't included in adaptive jumps*/
            m = *((gsl_matrix **)dummyParam->value);
            i += (int)( m->size1*m->size2 );
          } else if (param->type == LALINFERENCE_UINT4Vector_t)
          {
            /*
             increment i by number of size of vectors --
             number of wavelets in glitch model is not
             part of adaptive proposal
             */
            v = *((UINT4Vector **)dummyParam->value);
            i += (int)( v->length );
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

    sprintf(tmpname,"%s_%s",param->name,ADAPTSUFFIX);
    if (!LALInferenceCheckVariable(runState->proposalArgs,tmpname))
    {
      fprintf(stderr, "Attempting to draw single-parameter jump for %s but cannot find sigma!\nError in %s, line %d.\n",
              param->name,__FILE__, __LINE__);
      exit(1);
    }
    REAL8 *sigma=LALInferenceGetVariable(runState->proposalArgs,tmpname);

    /* Save the name of the proposed variable */
    if(LALInferenceCheckVariable(args,"proposedVariableName")){
      char *nameBuffer=*(char **)LALInferenceGetVariable(args,"proposedVariableName");
      strncpy(nameBuffer, param->name, MAX_STRLEN-1);
    }

    *((REAL8 *)param->value) += gsl_ran_ugaussian(rng) * *sigma * sqrtT;

    LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);

    /* Set the log of the proposal ratio to zero, since this is a
       symmetric proposal. */
    LALInferenceSetLogProposalRatio(runState, 0.0);

    INT4 as = 1;
    LALInferenceSetVariable(args, "adaptableStep", &as);

  }
}

void LALInferenceSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  const char *propName = singleProposalName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  gsl_rng * GSLrandom=runState->GSLrandom;
  LALInferenceVariableItem *param=NULL, *dummyParam=NULL;
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 T = 1.0;
  if(LALInferenceCheckVariable(runState->proposalArgs,"temperature")) T=*(REAL8 *)LALInferenceGetVariable(runState->proposalArgs, "temperature");

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
  } while (param->vary == LALINFERENCE_PARAM_FIXED || param->vary == LALINFERENCE_PARAM_OUTPUT || param->type != LALINFERENCE_REAL8_t);

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

  if (LALInferenceGetProcParamVal(runState->commandLine, "--zeroLogLike") || LALInferenceGetProcParamVal(runState->commandLine,"--zerologlike")) {
    if (!strcmp(param->name, "massratio")) {
      sigma = 0.02;
    } else if (!strcmp(param->name, "asym_massratio")) {
      sigma = 0.08;
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
    } else if (!strcmp(param->name, "inclination")||!strcmp(param->name,"theta_JN")) {
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
    if (!strcmp(param->name,"massratio") || !strcmp(param->name,"asym_massratio") || !strcmp(param->name,"time") || !strcmp(param->name,"a_spin2") || !strcmp(param->name,"a_spin1")){
      *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
    } else if (!strcmp(param->name,"polarisation") || !strcmp(param->name,"phase") || !strcmp(param->name,"inclination") ||!strcmp(param->name,"theta_JN")){
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
  const char *propName = orbitalPhaseJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  REAL8 phi;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  phi = *((REAL8 *) LALInferenceGetVariable(proposedParams, "phase"));

  phi = fmod(phi+M_PI, 2.0*M_PI);

  LALInferenceSetVariable(proposedParams, "phase", &phi);

  LALInferenceSetLogProposalRatio(runState, 0.0);

  /* Probably not needed, but play it safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void LALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = covarianceEigenvectorJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceVariables *proposalArgs = runState->proposalArgs;
  gsl_matrix *eigenvectors = *((gsl_matrix **)LALInferenceGetVariable(proposalArgs, "covarianceEigenvectors"));
  REAL8Vector *eigenvalues = *((REAL8Vector **)LALInferenceGetVariable(proposalArgs, "covarianceEigenvalues"));
  REAL8 temp = 1.0;
  if(LALInferenceCheckVariable(proposalArgs,"temperature"))
    temp=*((REAL8 *)LALInferenceGetVariable(proposalArgs, "temperature"));
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
    if (proposeIterator->vary != LALINFERENCE_PARAM_FIXED && proposeIterator->vary != LALINFERENCE_PARAM_OUTPUT && strcmp(proposeIterator->name,"psdscale")) {
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
  const char *propName = skyLocWanderJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  gsl_rng *rng = runState->GSLrandom;
  LALInferenceVariables *proposalArgs = runState->proposalArgs;
  REAL8 temp = 1.0;
  if(LALInferenceCheckVariable(proposalArgs,"temperature"))
    temp=*((REAL8 *)LALInferenceGetVariable(proposalArgs, "temperature"));
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
  const char *propName = differentialEvolutionFullName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceDifferentialEvolutionNames(runState, proposedParams, NULL);
}


void LALInferenceDifferentialEvolutionNames(LALInferenceRunState *runState,
                                            LALInferenceVariables *proposedParams,
                                            const char **names) {
  size_t i,j;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  if (names == NULL) {

    size_t N = LALInferenceGetVariableDimension(runState->currentParams) + 1; /* More names than we need. */
    names = alloca(N*sizeof(char *)); /* Hope we have alloca---saves
                                         having to deallocate after
                                         proposal. */

    LALInferenceVariableItem *item = runState->currentParams->head;
    i = 0;
    while (item != NULL) {
      if (item->vary != LALINFERENCE_PARAM_FIXED && item->vary != LALINFERENCE_PARAM_OUTPUT && item->type==LALINFERENCE_REAL8_t ) {
        names[i] = item->name;
        i++;
      }

      item = item->next;
    }
    names[i]=NULL; /* Terminate */
  }


  size_t Ndim = 0;
  for(Ndim=0,i=0; names[i] != NULL; i++ ) {
    if(LALInferenceCheckVariableNonFixed(proposedParams,names[i]))
      Ndim++;
  }
  LALInferenceVariables **dePts = runState->differentialPoints;
  size_t nPts = runState->differentialPointsLength;

  if (dePts == NULL || nPts <= 1) {
    LALInferenceSetLogProposalRatio(runState, 0.0);
    return; /* Quit now, since we don't have any points to use. */
  }

  i = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  do {
    j = gsl_rng_uniform_int(runState->GSLrandom, nPts);
  } while (j == i);

  LALInferenceVariables *ptI = dePts[i];
  LALInferenceVariables *ptJ = dePts[j];


  REAL8 scale;

  const REAL8 modeHoppingFrac = 0.5;
  /* Some fraction of the time, we do a "mode hopping" jump,
     where we jump exactly along the difference vector. */
  if (gsl_rng_uniform(runState->GSLrandom) < modeHoppingFrac) {
      scale = 1.0;
  } else {
  /* Otherwise scale is chosen uniform in log between 0.1 and 10 times the
     desired jump size. */
      scale = 2.38/sqrt(Ndim) * exp(log(0.1) + log(100.0)*gsl_rng_uniform(runState->GSLrandom));
  }


  for (i = 0; names[i] != NULL; i++) {
    if (!LALInferenceCheckVariableNonFixed(proposedParams, names[i]) || !LALInferenceCheckVariable(ptJ, names[i]) || !LALInferenceCheckVariable(ptI, names[i])) {
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

void LALInferenceDifferentialEvolutionIntrinsic(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *propName = differentialEvolutionIntrinsicName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  const char *names[] = {"chirpmass", "asym_massratio", "massratio", "m1", "m2", "a_spin1", "a_spin2",
      "tilt_spin1", "tilt_spin2", "phi12", "phi_spin1", "phi_spin2", "theta_spin1", "theta_spin2", NULL};
  LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void LALInferenceDifferentialEvolutionExtrinsic(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *propName = differentialEvolutionExtrinsicName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);

  const char *names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", "phase", "time", "theta_JN", NULL};
  const char *marg_time_names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", "phase", "theta_JN", NULL};
  const char *marg_time_phase_names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", "theta_JN", NULL};

  if (LALInferenceGetProcParamVal(runState->commandLine, "--margtimephi"))
    LALInferenceDifferentialEvolutionNames(runState, pp, marg_time_phase_names);
  else if (LALInferenceGetProcParamVal(runState->commandLine, "--margtime"))
    LALInferenceDifferentialEvolutionNames(runState, pp, marg_time_names);
  else
    LALInferenceDifferentialEvolutionNames(runState, pp, names);
}

static REAL8
draw_distance(LALInferenceRunState *runState) {
  REAL8 dmin, dmax;

  LALInferenceGetMinMaxPrior(runState->priorArgs, "distance", &dmin, &dmax);

  REAL8 x = gsl_rng_uniform(runState->GSLrandom);

  return cbrt(x*(dmax*dmax*dmax - dmin*dmin*dmin) + dmin*dmin*dmin);
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

  REAL8 mMin56 = pow(min, 5.0/6.0);
  REAL8 mMax56 = pow(max, 5.0/6.0);

  REAL8 delta = 1.0/mMin56 - 1.0/mMax56;

  REAL8 u = delta*gsl_rng_uniform(runState->GSLrandom);

  return pow(1.0/(1.0/mMin56 - u), 6.0/5.0);
}

static REAL8
approxLogPrior(LALInferenceVariables *params) {
  REAL8 logP = 0.0;

  REAL8 Mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
  logP += -11.0/6.0*log(Mc);

  /* Flat in eta. */

  if (LALInferenceCheckVariable(params, "inclination")) {
    REAL8 iota = *(REAL8 *)LALInferenceGetVariable(params, "inclination");
    logP += log(sin(iota));
  }else if (LALInferenceCheckVariable(params, "theta_JN")) {
    REAL8 thetaJN = *(REAL8 *)LALInferenceGetVariable(params, "theta_JN");
    logP += log(sin(thetaJN));
  }

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
  const char *propName = drawApproxPriorName;

  REAL8 tmp = 0.0;
  UINT4 analyticTest = 0;
  REAL8 logBackwardJump;

  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  if (runState->likelihood==&LALInferenceCorrelatedAnalyticLogLikelihood ||
      runState->likelihood==&LALInferenceBimodalCorrelatedAnalyticLogLikelihood ||
      runState->likelihood==&LALInferenceRosenbrockLogLikelihood) {
    analyticTest = 1;
  }

  if (analyticTest) {
    LALInferenceVariableItem *ptr = runState->currentParams->head;
    while(ptr!=NULL) {
      if(ptr->vary != LALINFERENCE_PARAM_FIXED) {
        tmp = draw_flat(runState, ptr->name);
        LALInferenceSetVariable(proposedParams, ptr->name, &tmp);
      }
      ptr=ptr->next;
    }
  } else {
    logBackwardJump = approxLogPrior(runState->currentParams);

    REAL8 Mc = draw_chirp(runState);
    LALInferenceSetVariable(proposedParams, "chirpmass", &Mc);

    if (LALInferenceCheckVariableNonFixed(runState->currentParams, "asym_massratio")) {
      REAL8 q = draw_flat(runState, "asym_massratio");
      LALInferenceSetVariable(proposedParams, "asym_massratio", &q);
    }
    else if (LALInferenceCheckVariableNonFixed(runState->currentParams, "massratio")) {
      REAL8 eta = draw_flat(runState, "massratio");
      LALInferenceSetVariable(proposedParams, "massratio", &eta);
    }

    if (LALInferenceCheckVariableNonFixed(runState->currentParams, "time")) {
      REAL8 theTime = draw_flat(runState, "time");
      LALInferenceSetVariable(proposedParams, "time", &theTime);
    }

    if (LALInferenceCheckVariableNonFixed(runState->currentParams, "phase")) {
      REAL8 phase = draw_flat(runState, "phase");
      LALInferenceSetVariable(proposedParams, "phase", &phase);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "inclination")) {
      REAL8 inc = draw_colatitude(runState, "inclination");
      LALInferenceSetVariable(proposedParams, "inclination", &inc);
    }

    REAL8 pol = draw_flat(runState, "polarisation");
    LALInferenceSetVariable(proposedParams, "polarisation", &pol);

    REAL8 dist = draw_distance(runState);
    LALInferenceSetVariable(proposedParams, "distance", &dist);

    REAL8 ra = draw_flat(runState, "rightascension");
    LALInferenceSetVariable(proposedParams, "rightascension", &ra);

    REAL8 dec = draw_dec(runState);
    LALInferenceSetVariable(proposedParams, "declination", &dec);

    if (LALInferenceCheckVariableNonFixed(proposedParams, "theta_JN")) {
      REAL8 thetaJN = draw_colatitude(runState, "theta_JN");
      LALInferenceSetVariable(proposedParams, "theta_JN", &thetaJN);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "phi_JL")) {
      REAL8 phiJL = draw_flat(runState, "phi_JL");
      LALInferenceSetVariable(proposedParams, "phi_JL", &phiJL);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "phi12")) {
      REAL8 phi12 = draw_flat(runState, "phi12");
      LALInferenceSetVariable(proposedParams, "phi12", &phi12);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "tilt_spin1")) {
      REAL8 tilt1 = draw_colatitude(runState, "tilt_spin1");
      LALInferenceSetVariable(proposedParams, "tilt_spin1", &tilt1);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "tilt_spin2")) {
      REAL8 tilt2 = draw_colatitude(runState, "tilt_spin2");
      LALInferenceSetVariable(proposedParams, "tilt_spin2", &tilt2);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "a_spin1")) {
      REAL8 a1 = draw_flat(runState, "a_spin1");
      LALInferenceSetVariable(proposedParams, "a_spin1", &a1);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "a_spin2")) {
      REAL8 a2 = draw_flat(runState, "a_spin2");
      LALInferenceSetVariable(proposedParams, "a_spin2", &a2);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "spin1")) {
      REAL8 a1 = draw_flat(runState, "spin1");
      LALInferenceSetVariable(proposedParams, "spin1", &a1);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "spin2")) {
      REAL8 a2 = draw_flat(runState, "spin2");
      LALInferenceSetVariable(proposedParams, "spin2", &a2);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "phi_spin1")) {
      REAL8 phi1 = draw_flat(runState, "phi_spin1");
      LALInferenceSetVariable(proposedParams, "phi_spin1", &phi1);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "phi_spin2")) {
      REAL8 phi2 = draw_flat(runState, "phi_spin2");
      LALInferenceSetVariable(proposedParams, "phi_spin2", &phi2);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "theta_spin1")) {
      REAL8 theta1 = draw_colatitude(runState, "theta_spin1");
      LALInferenceSetVariable(proposedParams, "theta_spin1", &theta1);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "theta_spin2")) {
      REAL8 theta2 = draw_colatitude(runState, "theta_spin2");
      LALInferenceSetVariable(proposedParams, "theta_spin2", &theta2);
    }

    if (LALInferenceCheckVariableNonFixed(proposedParams, "psdscale")) {
      REAL8 x, min, max;
      UINT4 i,j;
      min=0.10;
      max=10.0;
      gsl_matrix *eta = *((gsl_matrix **)LALInferenceGetVariable(proposedParams, "psdscale"));

      for(i=0;i<(UINT8)eta->size1;i++)
      {
        for(j=0;j<(UINT8)eta->size2;j++)
        {
          x = min + gsl_rng_uniform(runState->GSLrandom)*(max - min);
          gsl_matrix_set(eta,i,j,x);
        }
      }

    }//end if(psdscale)
  }

  if (analyticTest) {
    /* Flat in every variable means uniform jump probability. */
    LALInferenceSetLogProposalRatio(runState, 0.0);
  } else {
    LALInferenceSetLogProposalRatio(runState, logBackwardJump - approxLogPrior(proposedParams));
  }
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
  cart[0] = cos(longi)*cos(lat);
  cart[1] = sin(longi)*cos(lat);
  cart[2] = sin(lat);
}

static void
cart_to_sph(const REAL8 cart[3], REAL8 *lat, REAL8 *longi) {
  *longi = atan2(cart[1], cart[0]);
  *lat = asin(cart[2] / sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]));
}

static void
reflected_position_and_time(LALInferenceRunState *runState, const REAL8 ra, const REAL8 dec, const REAL8 oldTime,
                            REAL8 *newRA, REAL8 *newDec, REAL8 *newTime) {
  LALStatus status;
  memset(&status,0,sizeof(status));
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

/*
static REAL8 evaluate_morlet_proposal(LALInferenceRunState *runState, UNUSED LALInferenceVariables *proposedParams, UNUSED int ifo, UNUSED int k)
*/
static REAL8 evaluate_morlet_proposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams, int ifo, int k)
{
  REAL8 prior = 0.0;

  REAL8 component_min,component_max;


  //Uniform priors on f,Q,t,phi

  component_min = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_f0_prior_min"));
  component_max = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_f0_prior_max"));
  prior -= log(component_max-component_min);

  component_min = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_Q_prior_min"));
  component_max = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_Q_prior_max"));
  prior -= log(component_max-component_min);


  component_min = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_Amp_prior_min"));
  component_max = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_Amp_prior_max"));
  prior -= log(component_max-component_min);


  component_min = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_t0_prior_min"));
  component_max = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_t0_prior_max"));
  prior -= log(component_max-component_min);

  component_min = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_phi_prior_min"));
  component_max = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"morlet_phi_prior_max"));
  prior -= log(component_max-component_min);

  //"Malmquist" prior on A
  /*
  gsl_matrix *glitch_f = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
  gsl_matrix *glitch_Q = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
  gsl_matrix *glitch_A = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");

  REAL8 A,f,Q;
  A = gsl_matrix_get(glitch_A,ifo,k);
  Q = gsl_matrix_get(glitch_Q,ifo,k);
  f = gsl_matrix_get(glitch_f,ifo,k);


  prior += logGlitchAmplitudeDensity(A,Q,f);


  REAL8 Anorm = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"glitch_norm");

  prior += logGlitchAmplitudeDensity(A*Anorm,Q,f);
   */
  return prior;
}

/*
static REAL8 glitchAmplitudeDraw(REAL8 Q, REAL8 f, gsl_rng *r)
{
  REAL8 SNR;
  REAL8 PIterm = 0.5*LAL_2_SQRTPI*LAL_SQRT1_2;
  REAL8 SNRPEAK = 5.0;

  UINT4 k;
  REAL8 den, alpha;
  REAL8 max;

  // x/a^2 exp(-x/a) prior on SNR. Peaks at x = a. Good choice is a=5

  max = 1.0/(SNRPEAK*LAL_E);


  // rejection sample. Envelope function on runs out to ten times the peak
  // don't even bother putting in this minute correction to the normalization
  // (it is a 5e-4 correction).

  k = 0;
  do
  {

    SNR = 20.0*SNRPEAK*gsl_rng_uniform(r);

    den = SNR/(SNRPEAK*SNRPEAK)*exp(-SNR/SNRPEAK);

    den /= max;

    alpha = gsl_rng_uniform(r);

    k++;

  } while(alpha > den);

  return SNR/sqrt( (PIterm*Q/f) );
}
*/
void LALInferenceSkyRingProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  UINT4 i,j,l,ifo,nifo,timeflag=0;
  const char *propName = skyRingProposalName;
  REAL8 baryTime;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  LIGOTimeGPS GPSlal;

  LALInferenceIFOData *dataPtr;
  dataPtr = runState->data;

  REAL8 dL       = *(REAL8 *)LALInferenceGetVariable(proposedParams, "distance");
  REAL8 ra       = *(REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension");
  REAL8 dec      = *(REAL8 *)LALInferenceGetVariable(proposedParams, "declination");
  REAL8 psi      = *(REAL8 *)LALInferenceGetVariable(proposedParams, "polarisation");
  if(LALInferenceCheckVariable(proposedParams,"time")){
    baryTime = *(REAL8 *)LALInferenceGetVariable(proposedParams, "time");
    timeflag=1;
  }
  else
  {
    baryTime = XLALGPSGetREAL8(&(runState->data->epoch));
  }

  REAL8 newRA, newDec, newTime, newPsi, newDL;

  XLALGPSSetREAL8(&GPSlal, baryTime);
  REAL8 gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  //remap gmst back to [0:2pi]
  REAL8 intpart;
  REAL8 decpart;
  gmst /= LAL_TWOPI;
  intpart = (int)( gmst );
  decpart = gmst - (REAL8)intpart;
  gmst = decpart*LAL_TWOPI;
  gmst = gmst < 0. ? gmst + LAL_TWOPI : gmst;

  /*
   line-of-sight vector
   */
  REAL8 k[3];
  k[0] = cos(gmst-ra)*cos(dec);
  k[1] =-sin(gmst-ra)*cos(dec);
  k[2] = sin(dec);

  REAL8 IFO1[3],IFO2[3];
  REAL8 IFOX[3];

  /*
   Store location for each detector
   */
  nifo=0;
  while(dataPtr != NULL)
  {
    dataPtr = dataPtr->next;
    nifo++;
  }

  gsl_matrix *IFO = gsl_matrix_alloc(nifo,3);

  dataPtr = runState->data;
  for(ifo=0; ifo<nifo; ifo++)
  {
    memcpy(IFOX, dataPtr->detector->location, 3*sizeof(REAL8));
    for(i=0; i<3; i++) gsl_matrix_set(IFO,ifo,i,IFOX[i]);
    dataPtr=dataPtr->next;
  }

  /*
   Randomly select two detectors from the network
   -this assumes there are no co-located detectors
   */
  i=j=0;
  while(i==j)
  {
    i=gsl_rng_uniform_int(runState->GSLrandom, nifo);
    j=gsl_rng_uniform_int(runState->GSLrandom, nifo);
  }

  for(l=0; l<3; l++)
  {
    IFO1[l]=gsl_matrix_get(IFO,i,l);
    IFO2[l]=gsl_matrix_get(IFO,j,l);
  }

  /*
   detector axis
   */
  REAL8 normalize;
  REAL8 n[3];

  normalize=0.0;
  for(i=0; i<3; i++)
  {
    n[i]  = IFO1[i]-IFO2[i];
    normalize += n[i]*n[i];
  }
  normalize = 1./sqrt(normalize);
  for(i=0; i<3; i++) n[i] *= normalize;

  /*
   rotation angle
   */
  REAL8 omega    = LAL_TWOPI*gsl_rng_uniform(runState->GSLrandom);
  REAL8 cosomega = cos(omega);
  REAL8 sinomega = sin(omega);
  REAL8 c1momega = 1.0 - cosomega;

  /*
   rotate k' = Rk
   */
  REAL8 kp[3];
  kp[0] = (c1momega*n[0]*n[0] + cosomega)     *k[0] + (c1momega*n[0]*n[1] - sinomega*n[2])*k[1] + (c1momega*n[0]*n[2] + sinomega*n[1])*k[2];
  kp[1] = (c1momega*n[0]*n[1] + sinomega*n[2])*k[0] + (c1momega*n[1]*n[1] + cosomega)     *k[1] + (c1momega*n[1]*n[2] - sinomega*n[0])*k[2];
  kp[2] = (c1momega*n[0]*n[2] - sinomega*n[1])*k[0] + (c1momega*n[1]*n[2] + sinomega*n[0])*k[1] + (c1momega*n[2]*n[2] + cosomega)     *k[2];

  /*
   convert k' back to ra' and dec'
   */
  newDec = asin(kp[2]);
  newRA  = atan2(kp[1],kp[0]) + gmst;
  if (newRA < 0.0)
    newRA += LAL_TWOPI;
  else if (newRA >= LAL_TWOPI)
    newRA -= LAL_TWOPI;
  /*
   compute new geocenter time using
   fixed arrival time at IFO1 (arbitrary)
   */
  REAL8 tx; //old time shift = k * n
  REAL8 ty; //new time shift = k'* n
  tx=ty=0;
  for(i=0; i<3; i++)
  {
    tx += -IFO1[i]*k[i] /LAL_C_SI;
    ty += -IFO1[i]*kp[i]/LAL_C_SI;
  }
  newTime = tx + baryTime - ty;

  XLALGPSSetREAL8(&GPSlal, newTime);
  REAL8 newGmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  /*
   draw new polarisation angle uniformally
   for now
   MARK: Need to be smarter about psi in sky-ring jump
   */
  newPsi = LAL_PI*gsl_rng_uniform(runState->GSLrandom);

  /*
   compute new luminosity distance,
   maintaining F+^2 + Fx^2 across the network
   */
  REAL8 Fx,Fy;
  REAL8 Fp,Fc;
  Fx=0;Fy=0;

  dataPtr = runState->data;
  while(dataPtr != NULL)
  {
    XLALComputeDetAMResponse(&Fp, &Fc, (const REAL4(*)[3])dataPtr->detector->response, ra, dec, psi, gmst);
    Fx += Fp*Fp+Fc*Fc;

    XLALComputeDetAMResponse(&Fp, &Fc, (const REAL4(*)[3])dataPtr->detector->response, newRA, newDec, newPsi, newGmst);
    Fy += Fp*Fp+Fc*Fc;

    dataPtr = dataPtr->next;
  }
  newDL = dL*sqrt(Fy/Fx);

  /*
   update new parameters and exit.  woo!
   */
  LALInferenceSetVariable(proposedParams, "distance",       &newDL);
  LALInferenceSetVariable(proposedParams, "polarisation",   &newPsi);
  LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
  LALInferenceSetVariable(proposedParams, "declination",    &newDec);
  if(timeflag) LALInferenceSetVariable(proposedParams, "time",           &newTime);

  REAL8 pForward, pReverse;
  pForward = cos(newDec);
  pReverse = cos(dec);
  gsl_matrix_free(IFO);
  LALInferenceSetLogProposalRatio(runState, log(pReverse/pForward));
}

void LALInferenceSkyReflectDetPlane(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = skyReflectDetPlaneName;
  int timeflag=0;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
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
  REAL8 baryTime;
  if(LALInferenceCheckVariable(proposedParams,"time")){
    baryTime = *(REAL8 *)LALInferenceGetVariable(proposedParams, "time");
    timeflag=1;
  }
  else
  {
    baryTime = XLALGPSGetREAL8(&(runState->data->epoch));
  }

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
  if(timeflag) LALInferenceSetVariable(proposedParams, "time", &newTime);
  LALInferenceSetLogProposalRatio(runState, log(pReverse/pForward));
}

static void
rotateVectorAboutAxis(REAL8 vrot[3],
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
  }

  vadd(vrot, vpar, vperprot);
}

static void
vectorToColatLong(const REAL8 v[3],
                  REAL8 *colat, REAL8 *longi) {
  *longi = atan2(v[1], v[0]);
  if (*longi < 0.0) {
    *longi += 2.0*M_PI;
  }

  *colat = acos(v[2] / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));
}

void LALInferencePSDFitJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  const char *propName = PSDFitJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);

  INT4 i,j;
  INT4 nifo;
  INT4 N;

  REAL8 draw=0.0;
  //REAL8 var = *(REAL8 *)(LALInferenceGetVariable(runState->proposalArgs, "psdsigma"));
  REAL8Vector *var = *((REAL8Vector **)LALInferenceGetVariable(runState->proposalArgs, "psdsigma"));

  //Get current state of chain into workable form
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  gsl_matrix *ny = *((gsl_matrix **)LALInferenceGetVariable(proposedParams, "psdscale"));

  //Get size of noise parameter array
  nifo = (int)ny->size1;
  N    = (int)ny->size2;

  //perturb noise parameter
  for(i=0; i<nifo; i++)
  {
    for(j=0; j<N; j++)
    {
      draw = gsl_matrix_get(ny,i,j) + gsl_ran_ugaussian(runState->GSLrandom)*var->data[j];
      gsl_matrix_set(ny,i,j,draw);
    }
  }
  LALInferenceSetLogProposalRatio(runState, 0.0);
}

static void UpdateWaveletSum(LALInferenceRunState *runState, LALInferenceVariables *proposedParams, gsl_matrix *glitchFD, UINT4 ifo, UINT4 n, UINT4 flag)
{
  UINT4 i=0,j=0;
  LALInferenceIFOData *dataPtr = runState->data;
  REAL8FrequencySeries *noiseASD = NULL;

  /* get dataPtr pointing to correct IFO */
  while(dataPtr!=NULL)
  {
    if(ifo==i) noiseASD = dataPtr->noiseASD;
    dataPtr=dataPtr->next;
    i++;
  }
  dataPtr = runState->data;

  REAL8 deltaT = dataPtr->timeData->deltaT;
  REAL8 Tobs   = (((double)dataPtr->timeData->data->length) * deltaT);
  REAL8 deltaF = 1.0 / Tobs;

  UINT4 lower  = (UINT4)ceil(dataPtr->fLow / deltaF);
  UINT4 upper  = (UINT4)floor(dataPtr->fHigh / deltaF);

  gsl_matrix *glitch_f = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
  gsl_matrix *glitch_Q = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
  gsl_matrix *glitch_A = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");
  gsl_matrix *glitch_t = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_t0");
  gsl_matrix *glitch_p = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_phi");

  REAL8 Q,Amp,t0,ph0,f0; //sine-Gaussian parameters
  REAL8 amparg,phiarg,Ai;//helpers for computing sineGaussian
  REAL8 gRe,gIm;         //real and imaginary parts of current glitch model
  Q    = gsl_matrix_get(glitch_Q,ifo,n);
  Amp  = gsl_matrix_get(glitch_A,ifo,n);
  t0   = gsl_matrix_get(glitch_t,ifo,n);
  ph0  = gsl_matrix_get(glitch_p,ifo,n);
  f0   = gsl_matrix_get(glitch_f,ifo,n);

  //6 x decay time of sine Gaussian (truncate how much glitch we compute)
  REAL8 tau  = Q/LAL_TWOPI/f0;
  UINT4 glitchLower = (int)floor((f0 - 1./tau)/deltaF);
  UINT4 glitchUpper = (int)floor((f0 + 1./tau)/deltaF);

  //set glitch model to zero
  if(flag==0)
  {
    for(j=lower; j<=upper; j++)
    {
      gsl_matrix_set(glitchFD,ifo,2*j,0.0);
      gsl_matrix_set(glitchFD,ifo,2*j+1,0.0);
    }
  }

  for(j=glitchLower; j<glitchUpper; j++)
  {
    if(j>=lower && j<=upper)
    {
      gRe = gsl_matrix_get(glitchFD,ifo,2*j);
      gIm = gsl_matrix_get(glitchFD,ifo,2*j+1);
      amparg = ((REAL8)j*deltaF - f0)*LAL_PI*tau;
      phiarg = LAL_PI*(REAL8)j + ph0 - LAL_TWOPI*(REAL8)j*deltaF*(t0-Tobs/2.);//TODO: SIMPLIFY PHASE FOR SINEGAUSSIAN
      Ai = Amp*tau*0.5*sqrt(LAL_PI)*exp(-amparg*amparg)*noiseASD->data->data[j]/sqrt(Tobs);

      switch(flag)
      {
          // Remove wavelet from model
        case -1:
          gRe -= Ai*cos(phiarg);
          gIm -= Ai*sin(phiarg);
          break;
          // Add wavelet to model
        case  1:
          gRe += Ai*cos(phiarg);
          gIm += Ai*sin(phiarg);
          break;
          // Replace model with wavelet
        case 0:
          gRe = Ai*cos(phiarg);
          gIm = Ai*sin(phiarg);
          break;
          //Do nothing
        default:
          break;
      }//end switch

      //update glitch model
      gsl_matrix_set(glitchFD,ifo,2*j,gRe);
      gsl_matrix_set(glitchFD,ifo,2*j+1,gIm);

    }//end upper/lower check
  }//end loop over glitch samples
}

static void phase_blind_time_shift(REAL8 *corr, REAL8 *corrf, COMPLEX16Vector *data1, COMPLEX16Vector *data2, LALInferenceIFOData *IFOdata)
{
  UINT4 i;

  UINT4 N  = IFOdata->timeData->data->length;   // Number of data points
  UINT4 N2 = IFOdata->freqData->data->length-1; // 1/2 number of data points (plus 1)

  REAL8 deltaF = IFOdata->freqData->deltaF;
  REAL8 deltaT = IFOdata->timeData->deltaT;

  UINT4 lower  = (UINT4)ceil(  IFOdata->fLow  / deltaF );
  UINT4 upper  = (UINT4)floor( IFOdata->fHigh / deltaF );

  COMPLEX16FrequencySeries *corrFD  = XLALCreateCOMPLEX16FrequencySeries("cf1", &(IFOdata->freqData->epoch), 0.0, deltaF, &lalDimensionlessUnit, N2+1);
  COMPLEX16FrequencySeries *corrfFD = XLALCreateCOMPLEX16FrequencySeries("cf2", &(IFOdata->freqData->epoch), 0.0, deltaF, &lalDimensionlessUnit, N2+1);

  REAL8TimeSeries *corrTD  = XLALCreateREAL8TimeSeries("ct1", &(IFOdata->timeData->epoch), 0.0, deltaT, &lalDimensionlessUnit, N);
  REAL8TimeSeries *corrfTD = XLALCreateREAL8TimeSeries("ct2", &(IFOdata->timeData->epoch), 0.0, deltaT, &lalDimensionlessUnit, N);

  REAL8Vector *psd = IFOdata->oneSidedNoisePowerSpectrum->data;

  //convolution of signal & template
  for (i=0; i < N2; i++)
  {
    /*
		corrFD->data->data[i].real_FIXME	= 0.0;
		corrFD->data->data[i].imag_FIXME	= 0.0;
    corrfFD->data->data[i].real_FIXME = 0.0;
    corrfFD->data->data[i].imag_FIXME = 0.0;
    */
    corrFD->data->data[i]  = crect(0.0,0.0);
    corrfFD->data->data[i] = crect(0.0,0.0);

    if(i>lower && i<upper)
    {
      /*
      corrFD->data->data[i].real_FIXME	= ( creal(data1->data[i])*creal(data2->data[i]) + cimag(data1->data[i])*cimag(data2->data[i])) / psd->data[i];
      corrFD->data->data[i].imag_FIXME	= ( cimag(data1->data[i])*creal(data2->data[i]) - creal(data1->data[i])*cimag(data2->data[i])) / psd->data[i];
      corrfFD->data->data[i].real_FIXME = ( creal(data1->data[i])*cimag(data2->data[i]) - cimag(data1->data[i])*creal(data2->data[i])) / psd->data[i];
      corrfFD->data->data[i].imag_FIXME = ( cimag(data1->data[i])*cimag(data2->data[i]) + creal(data1->data[i])*creal(data2->data[i])) / psd->data[i];
      */
      corrFD->data->data[i]	 = crect( ( creal(data1->data[i])*creal(data2->data[i]) + cimag(data1->data[i])*cimag(data2->data[i])) / psd->data[i], ( cimag(data1->data[i])*creal(data2->data[i]) - creal(data1->data[i])*cimag(data2->data[i])) / psd->data[i] );
      corrfFD->data->data[i] = crect( ( creal(data1->data[i])*cimag(data2->data[i]) - cimag(data1->data[i])*creal(data2->data[i])) / psd->data[i], ( cimag(data1->data[i])*cimag(data2->data[i]) + creal(data1->data[i])*creal(data2->data[i])) / psd->data[i] );
    }

  }

  //invFFT convolutions to find time offset
  XLALREAL8FreqTimeFFT(corrTD, corrFD, IFOdata->freqToTimeFFTPlan);
  XLALREAL8FreqTimeFFT(corrfTD, corrfFD, IFOdata->freqToTimeFFTPlan);

  for (i=0; i < N; i++)
  {
		corr[i]  = corrTD->data->data[i];
		corrf[i] = corrfTD->data->data[i];
  }

  XLALDestroyREAL8TimeSeries(corrTD);
  XLALDestroyREAL8TimeSeries(corrfTD);
  XLALDestroyCOMPLEX16FrequencySeries(corrFD);
  XLALDestroyCOMPLEX16FrequencySeries(corrfFD);
}

static void MaximizeGlitchParameters(LALInferenceVariables *currentParams, LALInferenceRunState *runState, int ifo, int n)//int N, int NI, double Tobs, COMPLEX16Vector **s, COMPLEX16Vector **h, double *params, REAL8Vector **Sn)
{
  UINT4 i;

  LALInferenceIFOData *dataPtr = runState->data;
  i=0;
  while(i<(UINT4)ifo)
  {
    i++;
    dataPtr = dataPtr->next;
  }


  UINT4 N = (UINT4)dataPtr->timeData->data->length;
  REAL8 deltaT = (REAL8)(dataPtr->timeData->deltaT);
  REAL8 Tobs = (REAL8)(deltaT*N);
  REAL8 sqTwoDeltaToverN = sqrt(2.0 * deltaT / ((double) N) );

  REAL8 deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
  UINT4 lower  = (UINT4)ceil(dataPtr->fLow / deltaF);
  UINT4 upper  = (UINT4)floor(dataPtr->fHigh / deltaF);

  COMPLEX16Vector *s = dataPtr->freqData->data;
  COMPLEX16Vector *h = XLALCreateCOMPLEX16Vector(N/2);
  COMPLEX16Vector *r = XLALCreateCOMPLEX16Vector(N/2);
  REAL8Vector *Sn = dataPtr->oneSidedNoisePowerSpectrum->data;

  /* Get parameters for new wavelet */
  UINT4Vector *gsize   = *(UINT4Vector **) LALInferenceGetVariable(currentParams, "glitch_size");

  //gsl_matrix *glitchFD = runState->data->glitch_y;//*(gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_FD");
  gsl_matrix *glitchFD = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_FD");
  gsl_matrix *glitch_A = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_Amp");
  gsl_matrix *glitch_t = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_t0");
  gsl_matrix *glitch_p = *(gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_phi");

  REAL8 Amp,t0,ph0; //sine-Gaussian parameters
  Amp  = gsl_matrix_get(glitch_A,ifo,n);
  t0   = gsl_matrix_get(glitch_t,ifo,n);
  ph0  = gsl_matrix_get(glitch_p,ifo,n);

  REAL8 *corr;
  REAL8 *AC,*AF;

  /* Make new wavelet */
  gsl_matrix *hmatrix = gsl_matrix_alloc(ifo+1,N);
  gsl_matrix_set_all(hmatrix, 0.0);

  UpdateWaveletSum(runState, currentParams, hmatrix, ifo, n, 1);

  /* Copy to appropriate template array*/
  REAL8 rho=0.0;
  REAL8 hRe,hIm;
  REAL8 gRe,gIm;
  for(i=0; i<N/2; i++)
  {
    hRe = 0.0;
    hIm = 0.0;
    gRe = 0.0;
    gIm = 0.0;
    /*
    r->data[i].real_FIXME = 0.0;
    r->data[i].imag_FIXME = 0.0;
    */
    r->data[i] = crect(0.0,0.0);

    if(i>lower && i<upper)
    {
      hRe = sqTwoDeltaToverN * gsl_matrix_get(hmatrix,ifo,2*i);
      hIm = sqTwoDeltaToverN * gsl_matrix_get(hmatrix,ifo,2*i+1);
      /*
      h->data[i].real_FIXME = hRe;
      h->data[i].imag_FIXME = hIm;
      */
      h->data[i] = crect(hRe,hIm);
      //compute SNR of new wavelet
      rho += (hRe*hRe + hIm*hIm) / Sn->data[i];

      //form up residual while we're in here (w/out new template)
      if(gsize->data[ifo]>0)
      {
        gRe = gsl_matrix_get(glitchFD,ifo,2*i);
        gIm = gsl_matrix_get(glitchFD,ifo,2*i+1);
      }
      /*
      r->data[i].real_FIXME = sqTwoDeltaToverN * (creal(s->data[i])/deltaT-gRe);
      r->data[i].imag_FIXME = sqTwoDeltaToverN * (cimag(s->data[i])/deltaT-gIm);
      */
      r->data[i] = crect(sqTwoDeltaToverN * (creal(s->data[i])/deltaT-gRe),sqTwoDeltaToverN * (cimag(s->data[i])/deltaT-gIm));
    }
  }
  rho*=4.0;

  /* Compute correlation of data & template */
  corr = XLALMalloc(sizeof(REAL8) * N);
  AF   = XLALMalloc(sizeof(REAL8) * N);
  AC   = XLALMalloc(sizeof(REAL8) * N);

  for(i=0; i<N; i++)
  {
    corr[i] = 0.0;
  }


  /* Cross-correlate template & residual */
  phase_blind_time_shift(AC, AF, r, h, dataPtr);

  for(i=0; i<N; i++)
  {
    corr[i] += sqrt(AC[i]*AC[i] + AF[i]*AF[i]);
  }

  /* Find element where correlation is maximized */
  REAL8 max = corr[0];
  UINT4 imax = 0;

  for(i=1; i<N; i++)
  {
    if(corr[i] > max)
    {
      max  = corr[i];
      imax = i;
    }
  }
  max *= 4.0;

  /* Use max correlation to rescale amplitude */
  //REAL8 dAmplitude = max/rho;

  /* Get phase shift at max correlation */
  REAL8 dPhase = atan2(AF[imax],AC[imax]);

  /* Compute time shift needed for propsed template */
  REAL8 dTime;
  if(imax < (N/2)-1) dTime = ( ( (double)imax           )/(double)N )*Tobs;
  else               dTime = ( ( (double)imax-(double)N )/(double)N )*Tobs;

  /* Shift template parameters accordingly */
  //printf("%g %g %g %g %g %g\n",dTime,dPhase,dAmplitude,max/rho, max, rho);
  //printf("max amp: %g->%g  max time:  %g->%g, %g\n",Amp,Amp*dAmplitude,t0,t0+dTime,sqTwoDeltaToverN);
  t0  += dTime;
  //Amp *= dAmplitude;
  ph0 -= dPhase;

  /* Map time & phase back in range if necessary */
  if(ph0 < 0.0)            ph0 += LAL_TWOPI;
  else if(ph0 > LAL_TWOPI) ph0 -= LAL_TWOPI;

  if(t0 < 0.0)       t0 += Tobs;
  else if(t0 > Tobs) t0 -= Tobs;

  gsl_matrix_set(glitch_t,ifo,n,t0);
  gsl_matrix_set(glitch_A,ifo,n,Amp);
  gsl_matrix_set(glitch_p,ifo,n,ph0);


  gsl_matrix_free(hmatrix);

  XLALDestroyCOMPLEX16Vector(h);
  XLALDestroyCOMPLEX16Vector(r);

  XLALFree(corr);
  XLALFree(AF);
  XLALFree(AC);

}

static void MorletDiagonalFisherMatrix(REAL8Vector *params, REAL8Vector *sigmas)
{
  REAL8 f0;
  REAL8 Q;
  REAL8 Amp;

  REAL8 sigma_t0;
  REAL8 sigma_f0;
  REAL8 sigma_Q;
  REAL8 sigma_Amp;
  REAL8 sigma_phi0;

  REAL8 SNR   = 0.0;
  REAL8 sqrt3 = 1.7320508;

  f0   = params->data[1];
  Q    = params->data[2];
  Amp  = params->data[3];

  SNR = Amp*sqrt(Q/(2.0*sqrt(LAL_TWOPI)*f0));

  // this caps the size of the proposed jumps
  if(SNR < 5.0) SNR = 5.0;

  sigma_t0   = 1.0/(LAL_TWOPI*f0*SNR);
  sigma_f0   = 2.0*f0/(Q*SNR);
  sigma_Q    = 2.0*Q/(sqrt3*SNR);
  sigma_Amp  = Amp/SNR;
  sigma_phi0 = 1.0/SNR;

  // Map diagonal Fisher elements to sigmas vector
  sigmas->data[0] = sigma_t0;
  sigmas->data[1] = sigma_f0;
  sigmas->data[2] = sigma_Q;
  sigmas->data[3] = sigma_Amp;
  sigmas->data[4] = sigma_phi0;

}

void LALInferenceGlitchMorletProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  const char *propName = GlitchMorletJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);

  UINT4 i,ifo;
  UINT4 n;

  REAL8 t0;
  REAL8 f0;
  REAL8 Q;
  REAL8 Amp;
  REAL8 phi0;

  REAL8 scale;

  REAL8 qyx;
  REAL8 qxy;

  /*
   Vectors to store wavelet parameters.
   Order:
   [0] t0
   [1] f0
   [2] Q
   [3] Amp
   [4] phi0
 */
  REAL8Vector *params_x = XLALCreateREAL8Vector(5);
  REAL8Vector *params_y = XLALCreateREAL8Vector(5);

  REAL8Vector *sigmas_x = XLALCreateREAL8Vector(5);
  REAL8Vector *sigmas_y = XLALCreateREAL8Vector(5);

  /* Copy parameter structures and get local pointers to glitch parameters */
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  /* Get glitch meta paramters (dimnsion, proposal) */
  UINT4Vector *gsize = *(UINT4Vector **) LALInferenceGetVariable(proposedParams, "glitch_size");

  //gsl_matrix *glitchFD = runState->data->glitch_y;
  gsl_matrix *glitchFD = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_FD");
  gsl_matrix *glitch_f = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
  gsl_matrix *glitch_Q = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
  gsl_matrix *glitch_A = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");
  gsl_matrix *glitch_t = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_t0");
  gsl_matrix *glitch_p = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_phi");

  REAL8 Anorm = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"glitch_norm");

  /* Choose which IFO */
  ifo = (UINT4)floor( gsl_rng_uniform(runState->GSLrandom)*(REAL8)(gsize->length) );

  /* Bail out of proposal if no wavelets */
  if(gsize->data[ifo]==0)
  {
    LALInferenceSetLogProposalRatio(runState, 0.0);

    XLALDestroyREAL8Vector(params_x);
    XLALDestroyREAL8Vector(params_y);

    XLALDestroyREAL8Vector(sigmas_x);
    XLALDestroyREAL8Vector(sigmas_y);

    return;
  }

  /* Choose which glitch */
  n   = (UINT4)floor( gsl_rng_uniform(runState->GSLrandom)*(REAL8)(gsize->data[ifo]));

  /* Remove wavlet form linear combination */
  UpdateWaveletSum(runState, proposedParams, glitchFD, ifo, n, -1);

  /* Get parameters of n'th glitch int params vector */
  t0   = gsl_matrix_get(glitch_t,ifo,n); //Centroid time
  f0   = gsl_matrix_get(glitch_f,ifo,n); //Frequency
  Q    = gsl_matrix_get(glitch_Q,ifo,n); //Quality
  Amp  = gsl_matrix_get(glitch_A,ifo,n); //Amplitude
  phi0 = gsl_matrix_get(glitch_p,ifo,n); //Centroid phase

  /* Map to params Vector and compute Fisher */
  params_x->data[0] = t0;
  params_x->data[1] = f0;
  params_x->data[2] = Q;
  params_x->data[3] = Amp*(0.25*Anorm);
  params_x->data[4] = phi0;

  MorletDiagonalFisherMatrix(params_x, sigmas_x);

  /* Jump from x -> y:  y = x + N[0,sigmas_x]*scale */
  scale = 0.4082482; // 1/sqrt(6)

  for(i=0; i<5; i++)
    params_y->data[i] = params_x->data[i] + gsl_ran_ugaussian(runState->GSLrandom)*sigmas_x->data[i]*scale;

  /* Set parameters of n'th glitch int params vector */
  /* Map to params Vector and compute Fisher */
  t0   = params_y->data[0];
  f0   = params_y->data[1];
  Q    = params_y->data[2];
  Amp  = params_y->data[3]/(0.25*Anorm);
  phi0 = params_y->data[4];

  gsl_matrix_set(glitch_t,ifo,n,t0);
  gsl_matrix_set(glitch_f,ifo,n,f0);
  gsl_matrix_set(glitch_Q,ifo,n,Q);
  gsl_matrix_set(glitch_A,ifo,n,Amp);
  gsl_matrix_set(glitch_p,ifo,n,phi0);

  /* Add wavlet to linear combination */
  UpdateWaveletSum(runState, proposedParams, glitchFD, ifo, n, 1);

  /* Now compute proposal ratio using Fisher at y */
  MorletDiagonalFisherMatrix(params_y, sigmas_y);


  REAL8 sx  = 1.0; // sigma
  REAL8 sy  = 1.0;
  REAL8 dx  = 1.0; // (params_x - params_y)/sigma
  REAL8 dy  = 1.0;
  REAL8 exy = 0.0; // argument of exponential part of q
  REAL8 eyx = 0.0;
  REAL8 nxy = 1.0; // argument of normalization part of q
  REAL8 nyx = 1.0; // (computed as product to avoid too many log()s )
  for(i=0; i<5; i++)
  {
    sx = scale*sigmas_x->data[i];
    sy = scale*sigmas_y->data[i];

    dx = (params_x->data[i] - params_y->data[i])/sx;
    dy = (params_x->data[i] - params_y->data[i])/sy;

    nxy *= sy;
    nyx *= sx;

    exy += -dy*dy/2.0;
    eyx += -dx*dx/2.0;
  }

  qyx = eyx - log(nyx); //probabiltiy of proposing y given x
  qxy = exy - log(nxy); //probability of proposing x given y


  LALInferenceSetLogProposalRatio(runState, qxy-qyx);

  XLALDestroyREAL8Vector(params_x);
  XLALDestroyREAL8Vector(params_y);

  XLALDestroyREAL8Vector(sigmas_x);
  XLALDestroyREAL8Vector(sigmas_y);

}

void LALInferenceGlitchMorletReverseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{

  const char *propName = GlitchMorletReverseJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);

  UINT4 i,n;
  UINT4 ifo;
  UINT4 rj,nx,ny;

  REAL8 draw;
  REAL8 val;
  /*
  REAL8 t,f;
  */
  REAL8 t=0,f=0,Q=0,A=0;

  REAL8 qx       = 0.0; //log amp proposals
  REAL8 qy       = 0.0;
  REAL8 qyx      = 0.0; //log pixel proposals
  REAL8 qxy      = 0.0;

  REAL8 pForward; //combined p() & q() probabilities for ...
  REAL8 pReverse; //...RJMCMC hastings ratio

  REAL8 temperature = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "temperature");

  gsl_matrix *params = NULL;

  /* Copy parameter structures and get local pointers to glitch parameters */
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  UINT4Vector *gsize = *(UINT4Vector **)LALInferenceGetVariable(proposedParams, "glitch_size");
  gsl_matrix *glitchFD = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_FD");
  //gsl_matrix *glitchFD = runState->data->glitch_y;//*(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_FD");

  //REAL8 Anorm = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"glitch_norm");

  UINT4 nmin = (UINT4)(*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"glitch_dim_min"));
  UINT4 nmax = (UINT4)(*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"glitch_dim_max"));

  INT4 adapting = *(INT4*) LALInferenceGetVariable(runState->proposalArgs, "adapting");

  /* Choose which IFO */
  ifo = (UINT4)floor( gsl_rng_uniform(runState->GSLrandom)*(REAL8)(gsize->length) );
  nx  = gsize->data[ifo];

  /* Choose birth or death move */
  draw = gsl_rng_uniform(runState->GSLrandom);
  if( (draw < 0.5 && nx < nmax) || nx == nmin ) rj = 1;
  else rj = -1;

  //find dimension of proposed model
  ny = nx + rj;

  switch(rj)
  {
    /* Birth */
    case 1:

      //Add new wavelet to glitch model
//      flag=0;
//      while(flag==0)
//      {
//        /* Rejection sample on time and frequency */
//
//        t = draw_flat(runState, "morlet_t0_prior");
//        f = draw_flat(runState, "morlet_f0_prior");
//
//        //find TF pixel
//        j = floor( log(f*T)/log(2.) );
//        i = floor( int_pow(2,j)*t/T );
//        k = int_pow(2,j) + i;
//        //draw U[0,max_power]
//        draw = -10.0;//gsl_rng_uniform(runState->GSLrandom)*max->data[ifo];
//
//        //if draw<max_power[ifo][k] exit loop
//        printf("%g %g\n",draw,gsl_matrix_get(power,ifo,k));
//        if(draw < gsl_matrix_get(power,ifo,k)) flag=1;
//        printf("here\n");
//      }
      t = draw_flat(runState, "morlet_t0_prior");
      f = draw_flat(runState, "morlet_f0_prior");

      //Centroid time
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_t0");
      gsl_matrix_set(params,ifo,nx,t);

      //Frequency
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
      gsl_matrix_set(params,ifo,nx,f);

      //Quality
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
      val    = draw_flat(runState, "morlet_Q_prior");
      gsl_matrix_set(params,ifo,nx,val);
      Q=val;

      //Amplitude
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");

      val    = draw_flat(runState, "morlet_Amp_prior");
      gsl_matrix_set(params,ifo,nx,val);
      /*
      val    = glitchAmplitudeDraw(Q, f, runState->GSLrandom)/Anorm;
      gsl_matrix_set(params,ifo,nx,val);
      A=val;
       */
      //Centroid phase
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_phi");
      val    = draw_flat(runState, "morlet_phi_prior");
      gsl_matrix_set(params,ifo,nx,val);

      //Maximize phase, time, and amplitude using cross-correlation of data & wavelet
      if(adapting)MaximizeGlitchParameters(proposedParams, runState, ifo, nx);

      //Add wavlet to linear combination
      UpdateWaveletSum(runState, proposedParams, glitchFD, ifo, nx, 1);

      //Compute probability of drawing parameters
      qy = evaluate_morlet_proposal(runState,proposedParams,ifo,nx);// + log(gsl_matrix_get(power,ifo,k));

      //Compute reverse probability of dismissing k
      qxy = 0.0;//-log((double)runState->data->freqData->data->length);//-log( (REAL8)ny );

      if(adapting) qy += 10.0/temperature;

      break;

    /* Death */
    case -1:

      //Choose wavelet to remove from glitch model
      draw = gsl_rng_uniform(runState->GSLrandom);
      n    = (UINT4)( floor( draw*(REAL8)nx ) );     //choose which hot pixel

      // Remove wavlet from linear combination
      UpdateWaveletSum(runState, proposedParams, glitchFD, ifo, n, -1);

      //Get t and f of removed wavelet
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
      f = gsl_matrix_get(params,ifo,n);
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_t0");
      t = gsl_matrix_get(params,ifo,n);
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
      Q = gsl_matrix_get(params,ifo,n);
      params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");
      A = gsl_matrix_get(params,ifo,n);

      //Shift morlet parameters to fill in array
      for(i=n; i<ny; i++)
      {
        params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_f0");
        gsl_matrix_set(params,ifo,i,gsl_matrix_get(params,ifo,i+1) );
        params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Q");
        gsl_matrix_set(params,ifo,i,gsl_matrix_get(params,ifo,i+1) );
        params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_Amp");
        gsl_matrix_set(params,ifo,i,gsl_matrix_get(params,ifo,i+1) );
        params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_t0");
        gsl_matrix_set(params,ifo,i,gsl_matrix_get(params,ifo,i+1) );
        params = *(gsl_matrix **)LALInferenceGetVariable(proposedParams, "morlet_phi");
        gsl_matrix_set(params,ifo,i,gsl_matrix_get(params,ifo,i+1) );
      }

      //Compute reverse probability of drawing parameters
      //find TF pixel
//      j = floor( log(f*T)/log(2.) );
//      i = floor( int_pow(2,j)*t/T );
//      k = int_pow(2,j) + i;

      qx = evaluate_morlet_proposal(runState,runState->currentParams,ifo,n);// + log(gsl_matrix_get(power,ifo,k));

      //Compute forward probability of dismissing k
      qyx = 0.0;//-log((double)runState->data->freqData->data->length);//0.0;//-log( (REAL8)nx );

      if(adapting) qx += 10.0/temperature;

      break;

    default:
      break;
  }

  /* Update proposal structure for return to MCMC */

  //Update model meta-date
  gsize->data[ifo]=ny;

  //Re-package prior and proposal ratios into runState
  pForward = qxy + qx;
  pReverse = qyx + qy;

  /*
   fprintf(stdout," Trying %i->%i\n",nx,ny);
   fprintf(stdout,"  {%f,%f,%f,%f}\n",f,t,Q,A);
   fprintf(stdout,"  pForward = %g = %g + %g\n",pForward,qx,qxy);
   fprintf(stdout,"  pReverse = %g = %g + %g\n",pReverse,qy,qyx);
   fprintf(stdout,"  qxyRatio = %g\n",pForward-pReverse);
  */

  LALInferenceSetLogProposalRatio(runState, pForward-pReverse);
}

void
LALInferenceRotateSpins(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = rotateSpinsName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  REAL8 iota=0.;

  if(LALInferenceCheckVariable(proposedParams,"theta_JN"))
    iota=*(REAL8 *)LALInferenceGetVariable(proposedParams,"theta_JN");
  else
    iota=*(REAL8 *)LALInferenceGetVariable(proposedParams,"inclination");

  REAL8 theta1 = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);
  REAL8 theta2 = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);

  REAL8 logPr = 0.0;

  if (LALInferenceCheckVariableNonFixed(proposedParams, "theta_spin1")) {
    REAL8 theta, phi;
    REAL8 s1[3], L[3], newS[3];

    theta = *(REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin1");
    phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin1");

    s1[0] = cos(phi)*sin(theta);
    s1[1] = sin(phi)*sin(theta);
    s1[2] = cos(theta);

    L[0] = sin(iota);
    L[1] = 0.0;
    L[2] = cos(iota);

    rotateVectorAboutAxis(newS, s1, L, theta1);

    REAL8 newPhi, newTheta;

    vectorToColatLong(newS, &newTheta, &newPhi);

    /* Since the proposal is inherently uniform on the surface of the
       sphere, we only need to account for the volume factors between
       cos(theta) and theta. */
    logPr += log(sin(theta)/sin(newTheta));

    LALInferenceSetVariable(proposedParams, "phi_spin1", &newPhi);
    LALInferenceSetVariable(proposedParams, "theta_spin1", &newTheta);
  }

  if (LALInferenceCheckVariableNonFixed(proposedParams, "theta_spin2")) {
    REAL8 theta, phi;
    REAL8 s2[3], L[3], newS[3];

    theta = *(REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin2");
    phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin2");

    s2[0] = cos(phi)*sin(theta);
    s2[1] = sin(phi)*sin(theta);
    s2[2] = cos(theta);

    L[0] = sin(iota);
    L[1] = 0.0;
    L[2] = cos(iota);

    rotateVectorAboutAxis(newS, s2, L, theta2);

    REAL8 newPhi, newTheta;

    vectorToColatLong(newS, &newTheta, &newPhi);

    /* Since the proposal is inherently uniform on the surface of the
       sphere, we only need to account for the volume factors between
       cos(theta) and theta. */
    logPr += log(sin(theta)/sin(newTheta));

    LALInferenceSetVariable(proposedParams, "phi_spin2", &newPhi);
    LALInferenceSetVariable(proposedParams, "theta_spin2", &newTheta);
  }

  LALInferenceSetLogProposalRatio(runState, logPr);
}

void
LALInferencePolarizationPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = polarizationPhaseJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 psi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "polarisation");
  REAL8 phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phase");

  phi += M_PI;
  psi += M_PI/2;

  phi = fmod(phi, 2.0*M_PI);
  psi = fmod(psi, M_PI);

  LALInferenceSetVariable(proposedParams, "polarisation", &psi);
  LALInferenceSetVariable(proposedParams, "phase", &phi);

  LALInferenceSetLogProposalRatio(runState, 0.0);
}

void LALInferenceCorrPolarizationPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  const char *propName = polarizationCorrPhaseJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

	REAL8 alpha,beta;
	REAL8 draw;

  REAL8 psi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "polarisation");
  REAL8 phi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "phase");

	alpha = psi + phi;
	beta  = psi - phi;

  //alpha =>   0:3pi
	//beta  => -2pi:pi

	//big jump in either alpha (beta) or beta (alpha)
  draw=gsl_rng_uniform(runState->GSLrandom);
	if(draw < 0.5) alpha = gsl_rng_uniform(runState->GSLrandom)*3.0*LAL_PI;
	else           beta  = -LAL_TWOPI+gsl_rng_uniform(runState->GSLrandom)*3.0*LAL_PI;

	//transform back to psi,phi space
	psi =  (alpha + beta)*0.5;
	phi =  (alpha - beta)*0.5;

  //map back in range
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);

  LALInferenceSetVariable(proposedParams, "polarisation", &psi);
  LALInferenceSetVariable(proposedParams, "phase", &phi);

  LALInferenceSetLogProposalRatio(runState, 0.0);
}
typedef enum {
  USES_DISTANCE_VARIABLE,
  USES_LOG_DISTANCE_VARIABLE
} DistanceParam;


void LALInferenceKDNeighborhoodProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  size_t NCell;
  LALInferenceVariables *proposalArgs = runState->proposalArgs;

  if (LALInferenceCheckVariable(runState->proposalArgs, "KDNCell")) {
    NCell = *(INT4 *)LALInferenceGetVariable(runState->proposalArgs, "KDNCell");
  } else if (LALInferenceCheckVariable(runState->proposalArgs, "kdncell")) {
    NCell = *(INT4 *)LALInferenceGetVariable(runState->proposalArgs, "kdncell");
  } else {
    /* NCell default value. */
    NCell = 64;
  }

  const char *propName = KDNeighborhoodProposalName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  if (!LALInferenceCheckVariable(proposalArgs, "kDTree") || !LALInferenceCheckVariable(proposalArgs, "kDTreeVariableTemplate")) {
    /* For whatever reason, the appropriate data are not set up in the
       proposalArgs, so just propose the current point again and
       bail. */
    LALInferenceSetLogProposalRatio(runState, 0.0);
    return;
  }

  LALInferenceKDTree *tree = *(LALInferenceKDTree **)LALInferenceGetVariable(proposalArgs, "kDTree");
  LALInferenceVariables *templt = *(LALInferenceVariables **)LALInferenceGetVariable(proposalArgs, "kDTreeVariableTemplate");
  /* If tree has zero points, bail. */
  if (tree->npts == 0) {
    LALInferenceSetLogProposalRatio(runState, 0.0);
    return;
  }

  REAL8 *currentPt = XLALCalloc(tree->dim, sizeof(REAL8));
  REAL8 *proposedPt = XLALCalloc(tree->dim, sizeof(REAL8));

  /* Get the coordinates of the current point. */
  LALInferenceKDVariablesToREAL8(runState->currentParams, currentPt, templt);

  /* A randomly-chosen point from those in the tree. */
  LALInferenceKDDrawEigenFrame(runState->GSLrandom, tree, proposedPt, NCell);
  LALInferenceKDREAL8ToVariables(proposedParams, proposedPt, templt);

  REAL8 logPropRatio = LALInferenceKDLogProposalRatio(tree, currentPt, proposedPt, NCell);

  LALInferenceSetLogProposalRatio(runState, logPropRatio);

  /* Cleanup the allocated storage for currentPt. */
  XLALFree(currentPt);
  XLALFree(proposedPt);
}

void
LALInferenceFrequencyBinJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = frequencyBinJumpName;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  REAL8 f0 = *(REAL8 *)LALInferenceGetVariable(proposedParams, "f0");
  REAL8 df = *(REAL8 *)LALInferenceGetVariable(proposedParams, "df");

  REAL8 plusminus = gsl_rng_uniform(runState->GSLrandom);
  if ( plusminus < 0.5 ) { f0 -= df; }
  else { f0 += df; }

  LALInferenceSetVariable(proposedParams, "f0", &f0);

  LALInferenceSetLogProposalRatio(runState, 0.0);
}


static void
reflected_extrinsic_parameters(LALInferenceRunState *runState, const REAL8 ra, const REAL8 dec, const REAL8 baryTime,
                               const REAL8 dist, const REAL8 iota, const REAL8 psi,
                               REAL8 *newRA, REAL8 *newDec, REAL8 *newTime,
                               REAL8 *newDist, REAL8 *newIota, REAL8 *newPsi) {

//This proposal needs to be called with exactly 3 independent detector locations.

  LIGOTimeGPS GPSlal;
  REAL8 R2[4];
  REAL8 newGmst;
  REAL8 dist2;

  XLALGPSSetREAL8(&GPSlal, baryTime);
  REAL8 gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  reflected_position_and_time(runState, ra, dec, baryTime, newRA, newDec, newTime);

  XLALGPSSetREAL8(&GPSlal, *newTime);
  newGmst = XLALGreenwichMeanSiderealTime(&GPSlal);

  dist2=dist*dist;

  REAL8 cosIota = cos(iota);
  REAL8 cosIota2 = cosIota*cosIota;

  double Fplus, Fcross, psi_temp;
  double x[4],y[4],x2[4],y2[4];
  int i=1,j=0;
  LALInferenceIFOData *dataPtr;

  dataPtr = runState->data;

  /* Loop over interferometers */
  while (dataPtr != NULL) {

    psi_temp = 0.0;
    XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, *newRA, *newDec, psi_temp, newGmst);
    j=i-1;
    while (j>0){
      if(Fplus==x[j]){
        dataPtr = dataPtr->next;
        XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, *newRA, *newDec, psi_temp, newGmst);
      }
      j--;
    }
    x[i]=Fplus;
    x2[i]=Fplus*Fplus;
    y[i]=Fcross;
    y2[i]=Fcross*Fcross;

    XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, ra, dec, psi, gmst);
    R2[i] = (((1.0+cosIota2)*(1.0+cosIota2))/(4.0*dist2))*Fplus*Fplus
    + ((cosIota2)/(dist2))*Fcross*Fcross;

    dataPtr = dataPtr->next;
    i++;
  }

  REAL8 a,a2,b;

  a=(R2[3]*x2[2]*y2[1] - R2[2]*x2[3]*y2[1] - R2[3]*x2[1]*y2[2] + R2[1]*x2[3]*y2[2] + R2[2]*x2[1]*y2[3] -
     R2[1]*x2[2]*y2[3]);
  a2=a*a;
  b=(-(R2[3]*x[1]*x2[2]*y[1]) + R2[2]*x[1]*x2[3]*y[1] + R2[3]*x2[1]*x[2]*y[2] - R2[1]*x[2]*x2[3]*y[2] +
     R2[3]*x[2]*y2[1]*y[2] - R2[3]*x[1]*y[1]*y2[2] - R2[2]*x2[1]*x[3]*y[3] + R2[1]*x2[2]*x[3]*y[3] - R2[2]*x[3]*y2[1]*y[3] + R2[1]*x[3]*y2[2]*y[3] +
     R2[2]*x[1]*y[1]*y2[3] - R2[1]*x[2]*y[2]*y2[3]);

  (*newPsi)=(2.*atan((b - a*sqrt((a2 + b*b)/(a2)))/a))/4.;

  while((*newPsi)<0){
    (*newPsi)=(*newPsi)+LAL_PI/4.0;
  }
  while((*newPsi)>LAL_PI/4.0){
    (*newPsi)=(*newPsi)-LAL_PI/4.0;
  }

  REAL8 newFplus[4], newFplus2[4], newFcross[4], newFcross2[4];

  for (i = 1; i < 4; i++){

    newFplus[i]=x[i]*cos(2.0*(*newPsi))+y[i]*sin(2.0*(*newPsi));
    newFplus2[i]=newFplus[i]*newFplus[i];

    newFcross[i]=y[i]*cos(2.0*(*newPsi))-x[i]*sin(2.0*(*newPsi));
    newFcross2[i]=newFcross[i]*newFcross[i];

  }

  REAL8 c12;

  c12 = -2.0*((R2[1]*(newFcross2[2])-R2[2]*(newFcross2[1]))
              /(R2[1]*(newFplus2[2])-R2[2]*(newFplus2[1])))-1.0;

  if(c12<1.0){
    c12 = (3.0-c12)/(1.0+c12);
    (*newPsi)=(*newPsi)+LAL_PI/4.0;

    for (i = 1; i < 4; i++){

      newFplus[i]=x[i]*cos(2.0*(*newPsi))+y[i]*sin(2.0*(*newPsi));
      newFplus2[i]=newFplus[i]*newFplus[i];

      newFcross[i]=y[i]*cos(2.0*(*newPsi))-x[i]*sin(2.0*(*newPsi));
      newFcross2[i]=newFcross[i]*newFcross[i];

    }
  }

  if(c12<1){
    *newIota=iota;
    *newDist=dist;
    return;
  }

  REAL8 cosnewIota, cosnewIota2;
  cosnewIota2 = c12-sqrt(c12*c12-1.0);
  cosnewIota = sqrt(cosnewIota2);
  *newIota = acos(cosnewIota);

  *newDist = sqrt((
                  ((((1.0+cosnewIota2)*(1.0+cosnewIota2))/(4.0))*newFplus2[1]
                   + (cosnewIota2)*newFcross2[1])
                  )/ R2[1]);

  if(Fplus*newFplus[3]<0){
    (*newPsi)=(*newPsi)+LAL_PI/2.;
    newFcross[3]=-newFcross[3];
  }

  if(Fcross*cosIota*cosnewIota*newFcross[3]<0){
    (*newIota)=LAL_PI-(*newIota);
  }

}


void LALInferenceExtrinsicParamProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  const char *propName = extrinsicParamProposalName;
  int timeflag=0;
  REAL8 baryTime;
  LALInferenceSetVariable(runState->proposalArgs, LALInferenceCurrentProposalName, &propName);
  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  int USES_THETA_JN=0;
  /* Find the number of distinct-position detectors. */
  /* Exit with same parameters (with a warning the first time) if
   there are not EXACTLY three unique detector locations. */
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

  DistanceParam distParam;

  if (LALInferenceCheckVariable(proposedParams, "distance")) {
    distParam = USES_DISTANCE_VARIABLE;
  } else if (LALInferenceCheckVariable(proposedParams, "logdistance")) {
    distParam = USES_LOG_DISTANCE_VARIABLE;
  } else {
    XLAL_ERROR_VOID(XLAL_FAILURE, "could not find 'distance' or 'logdistance' in current params");
  }

  REAL8 ra = *(REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension");
  REAL8 dec = *(REAL8 *)LALInferenceGetVariable(proposedParams, "declination");
  if(LALInferenceCheckVariable(proposedParams,"time")){
    baryTime = *(REAL8 *)LALInferenceGetVariable(proposedParams, "time");
    timeflag=1;
  }
  else
  {
    baryTime = XLALGPSGetREAL8(&(runState->data->epoch));
  }
  REAL8 iota=0.;
  
  if(LALInferenceCheckVariable(proposedParams,"inclination"))
    iota = *(REAL8 *)LALInferenceGetVariable(proposedParams, "inclination");
  else if(LALInferenceCheckVariable(proposedParams,"theta_JN"))
  {
    iota = *(REAL8 *)LALInferenceGetVariable(proposedParams, "theta_JN");
    USES_THETA_JN=1;
  }
  else fprintf(stderr,"LALInferenceExtrinsicParamProposal: No inclination or theta_JN parameter!\n");
  
  REAL8 psi = *(REAL8 *)LALInferenceGetVariable(proposedParams, "polarisation");
  REAL8 dist;
  if (distParam == USES_DISTANCE_VARIABLE) {
    dist = *(REAL8 *)LALInferenceGetVariable(proposedParams, "distance");
  } else {
    dist = exp(*(REAL8 *)LALInferenceGetVariable(proposedParams, "logdistance"));
  }

  REAL8 newRA, newDec, newTime, newDist, newIota, newPsi;

  reflected_extrinsic_parameters(runState, ra, dec, baryTime, dist, iota, psi, &newRA, &newDec, &newTime, &newDist, &newIota, &newPsi);

  /* Unit normal deviates, used to "fuzz" the state. */
  REAL8 nRA, nDec, nTime, nDist, nIota, nPsi;
  const REAL8 epsDist = 1e-8;
  const REAL8 epsTime = 1e-8;
  const REAL8 epsAngle = 1e-8;

  nRA = gsl_ran_ugaussian(runState->GSLrandom);
  nDec = gsl_ran_ugaussian(runState->GSLrandom);
  nTime = gsl_ran_ugaussian(runState->GSLrandom);
  nDist = gsl_ran_ugaussian(runState->GSLrandom);
  nIota = gsl_ran_ugaussian(runState->GSLrandom);
  nPsi = gsl_ran_ugaussian(runState->GSLrandom);

  newRA += epsAngle*nRA;
  newDec += epsAngle*nDec;
  newTime += epsTime*nTime;
  newDist += epsDist*nDist;
  newIota += epsAngle*nIota;
  newPsi += epsAngle*nPsi;

  /* And the doubly-reflected position (near the original, but not
   exactly due to the fuzzing). */
  REAL8 refRA, refDec, refTime, refDist, refIota, refPsi;
  reflected_extrinsic_parameters(runState, newRA, newDec, newTime, newDist, newIota, newPsi, &refRA, &refDec, &refTime, &refDist, &refIota, &refPsi);

  /* The Gaussian increments required to shift us back to the original
   position from the doubly-reflected position. */
  REAL8 nRefRA, nRefDec, nRefTime, nRefDist, nRefIota, nRefPsi;
  nRefRA = (ra - refRA)/epsAngle;
  nRefDec = (dec - refDec)/epsAngle;
  nRefTime = (baryTime - refTime)/epsTime;
  nRefDist = (dist - refDist)/epsDist;
  nRefIota = (iota - refIota)/epsAngle;
  nRefPsi = (psi - refPsi)/epsAngle;

  REAL8 pForward, pReverse;
  REAL8 cst = log(1./(sqrt(2.*LAL_PI)));
  pReverse = 6*cst-0.5*(nRefRA*nRefRA+nRefDec*nRefDec+nRefTime*nRefTime+nRefDist*nRefDist+nRefIota*nRefIota+nRefPsi*nRefPsi);
  pForward = 6*cst-0.5*(nRA*nRA+nDec*nDec+nTime*nTime+nDist*nDist+nIota*nIota+nPsi*nPsi);

  LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
  LALInferenceSetVariable(proposedParams, "declination", &newDec);
  if(timeflag) LALInferenceSetVariable(proposedParams, "time", &newTime);
  if (distParam == USES_DISTANCE_VARIABLE) {
    LALInferenceSetVariable(proposedParams, "distance", &newDist);
  } else {
    REAL8 logNewDist = log(newDist);
    LALInferenceSetVariable(proposedParams, "logdistance", &logNewDist);
  }
  if(USES_THETA_JN) LALInferenceSetVariable(proposedParams,"theta_JN",&newIota);
  else LALInferenceSetVariable(proposedParams, "inclination", &newIota);
  LALInferenceSetVariable(proposedParams, "polarisation", &newPsi);

  LALInferenceSetLogProposalRatio(runState, pReverse-pForward);

  return;
}


void NSWrapMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
  /* PTMCMC likes to read currentParams directly, whereas NS expects proposedParams
   to be modified by the proposal. Back up currentParams and then restore it after
   calling the MCMC proposal function. */
  REAL8 oldlogdist=-1.0,oldlogmc=-1.0;
  REAL8 newdist,newmc;
  LALInferenceVariables *currentParamsBackup=runState->currentParams;
  /* Create the proposal if none exists */
  if (!LALInferenceCheckVariable(runState->proposalArgs, cycleArrayName) || !LALInferenceCheckVariable(runState->proposalArgs, cycleArrayLengthName))
   {
    /* In case there is a partial cycle set up already, delete it. */
    LALInferenceDeleteProposalCycle(runState);
    if(LALInferenceGetProcParamVal(runState->commandLine,"--mcmcprop"))
	   { SetupDefaultProposal(runState, proposedParams); }
	 else {
	 	LALInferenceSetupDefaultNSProposal(runState,proposedParams);
	 }
  }

  /* PTMCMC expects some variables that NS doesn't use by default, so create them */

  if(LALInferenceCheckVariable(proposedParams,"logdistance"))
    oldlogdist=*(REAL8 *)LALInferenceGetVariable(proposedParams,"logdistance");
  if(LALInferenceCheckVariable(proposedParams,"logmc"))
    oldlogmc=*(REAL8*)LALInferenceGetVariable(proposedParams,"logmc");

  NSFillMCMCVariables(proposedParams,runState->priorArgs);

  runState->currentParams=proposedParams;
	  LALInferenceCyclicProposal(runState,proposedParams);
  /* Restore currentParams */
  runState->currentParams=currentParamsBackup;

  /* If the remapped variables are not updated do it here */
  if(oldlogdist!=-1.0)
    if(oldlogdist==*(REAL8*)LALInferenceGetVariable(proposedParams,"logdistance"))
      {
		newdist=*(REAL8*)LALInferenceGetVariable(proposedParams,"distance");
		newdist=log(newdist);
		LALInferenceSetVariable(proposedParams,"logdistance",&newdist);
      }
  if(oldlogmc!=-1.0)
    if(oldlogmc==*(REAL8*)LALInferenceGetVariable(proposedParams,"logmc"))
    {
      newmc=*(REAL8*)LALInferenceGetVariable(proposedParams,"chirpmass");
      newmc=log(newmc);
      LALInferenceSetVariable(proposedParams,"logmc",&newmc);
    }

}

/** Setup adaptive proposals. Should be called when state->currentParams is already filled with an initial sample */
void LALInferenceSetupAdaptiveProposals(LALInferenceRunState *state)
{
        INT4 adaptationOn=1;
        LALInferenceVariableItem *this=state->currentParams->head;
        if (LALInferenceGetProcParamVal(state->commandLine, "--noAdapt"))
                adaptationOn=0;

        for(this=state->currentParams->head;this;this=this->next)
        {
                char *name=this->name;
                REAL8 sigma=0.01;
                if (!strcmp(name,"massratio") || !strcmp(name,"asym_massratio") || !strcmp(name,"time") || !strcmp(name,"a_spin2") || !strcmp(name,"a_spin1")){
                        sigma = 0.001;
                } else if (!strcmp(name,"polarisation") || !strcmp(name,"phase") || !strcmp(name,"inclination")||!strcmp(name,"theta_JN")){
                        sigma = 0.1;
                }
                /* Set up variables to store current sigma, proposed and accepted */
                char varname[MAX_STRLEN]="";
                sprintf(varname,"%s_%s",name,ADAPTSUFFIX);
                LALInferenceAddVariable(state->proposalArgs,varname,&sigma,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
                sigma=0.0;
                sprintf(varname,"%s_%s",name,ACCEPTSUFFIX);
                LALInferenceAddVariable(state->proposalArgs,varname,&sigma,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
                sprintf(varname,"%s_%s",name,PROPOSEDSUFFIX);
                LALInferenceAddVariable(state->proposalArgs,varname,&sigma,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
                if(LALInferenceCheckVariable(state->algorithmParams,"verbose")) printf("Setup adaptive proposal for %s\n",name);
        }


        INT4 adapting = adaptationOn;      // Indicates if current iteration is being adapted
        LALInferenceAddVariable(state->proposalArgs, "adaptationOn", &adaptationOn, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
        LALInferenceAddVariable(state->proposalArgs, "adapting", &adapting, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

        INT4 adaptableStep = 0;
        LALInferenceAddVariable(state->proposalArgs, "adaptableStep", &adaptableStep, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
        char *nameBuffer=XLALCalloc(MAX_STRLEN,sizeof(char));
        sprintf(nameBuffer,"none");
        LALInferenceAddVariable(state->proposalArgs, "proposedVariableName", &nameBuffer, LALINFERENCE_string_t, LALINFERENCE_PARAM_OUTPUT);

        INT4 tau = 5;
        LALInferenceAddVariable(state->proposalArgs, "adaptTau", &tau, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

        ProcessParamsTable *ppt = LALInferenceGetProcParamVal(state->commandLine, "--adaptTau");
        if (ppt) {
                tau = atof(ppt->value);
                fprintf(stdout, "Setting adapt tau = %i.\n", tau);
                LALInferenceSetVariable(state->proposalArgs, "adaptTau", &tau);
        }
        INT4  adaptTau     = *((INT4 *)LALInferenceGetVariable(state->proposalArgs, "adaptTau"));     // Sets decay of adaption function
        INT4  adaptLength       = pow(10,adaptTau);   // Number of iterations to adapt before turning off
        INT4  adaptResetBuffer  = 100;                // Number of iterations before adapting after a restart
        REAL8 s_gamma           = 1.0;                // Sets the size of changes to jump size during adaptation
        INT4  adaptStart        = 0;                  // Keeps track of last iteration adaptation was restarted
        REAL8 logLAtAdaptStart  = 0.0;                // max log likelihood as of last adaptation restart
        LALInferenceAddVariable(state->proposalArgs, "adaptLength", &adaptLength,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddVariable(state->proposalArgs, "adaptResetBuffer", &adaptResetBuffer,  LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddVariable(state->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddVariable(state->proposalArgs, "adaptStart", &adaptStart, LALINFERENCE_INT4_t, LALINFERENCE_PARAM_LINEAR);
        LALInferenceAddVariable(state->proposalArgs, "logLAtAdaptStart", &logLAtAdaptStart, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
        return;
}


/** Update proposal statistics if tracking */
void LALInferenceTrackProposalAcceptance(LALInferenceRunState *runState, INT4 accepted){
    const char *currentProposalName;
    LALInferenceProposalStatistics *propStat;

    /* Update proposal statistics */
    if (runState->proposalStats){
        currentProposalName = *((const char **)LALInferenceGetVariable(runState->proposalArgs, LALInferenceCurrentProposalName));
        propStat = ((LALInferenceProposalStatistics *)LALInferenceGetVariable(runState->proposalStats, currentProposalName));
        propStat->proposed++;
        if (accepted == 1){
            propStat->accepted++;
        }
    }
    return;
}

/* Zero out proposal statistics counters */
void LALInferenceZeroProposalStats(LALInferenceRunState *runState){
    LALInferenceVariables *propStats = runState->proposalStats;
    if(propStats==NULL) return;

    LALInferenceVariableItem *ptr=propStats->head;
    while(ptr!=NULL) {
      ((LALInferenceProposalStatistics *) ptr->value)->proposed = 0;
      ((LALInferenceProposalStatistics *) ptr->value)->accepted = 0;
      ptr=ptr->next;
    }
    return;
}

/** Update the adaptive proposal. Whether or not a jump was accepted is passed with accepted */
void LALInferenceUpdateAdaptiveJumps(LALInferenceRunState *runState, INT4 accepted, REAL8 targetAcceptance){
        INT4 *adaptableStep = NULL;
        INT4 *adapting = NULL;

        if( LALInferenceCheckVariable(runState->proposalArgs, "adaptableStep" ) &&
                        LALInferenceCheckVariable(runState->proposalArgs, "adapting" ) ){
                adaptableStep = ((INT4 *)LALInferenceGetVariable(runState->proposalArgs,
                                        "adaptableStep"));
                adapting = ((INT4 *)LALInferenceGetVariable(runState->proposalArgs, "adapting"));
        }
        /* Don't do anything if these are not found */
        else return;

        if (*adaptableStep && *adapting) {
                char *name=*(char **)LALInferenceGetVariable(runState->proposalArgs,"proposedVariableName");
                char tmpname[MAX_STRLEN]="";

                sprintf(tmpname,"%s_%s",name,PROPOSEDSUFFIX);
                REAL8 *propose=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
                *propose+=1;
                sprintf(tmpname,"%s_%s",name,ACCEPTSUFFIX);
                REAL8 *accept=(REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);
                if(accepted == 1){
                        *accept+=1;
                }
        }
        /* Adapt if desired. */
        if (LALInferenceCheckVariable(runState->proposalArgs, "proposedVariableName") &&
            LALInferenceCheckVariable(runState->proposalArgs, "s_gamma") &&
            LALInferenceCheckVariable(runState->proposalArgs, "adapting") &&
            LALInferenceCheckVariable(runState->proposalArgs, "adaptableStep")) {

                if (*adaptableStep) {
                        char *name=*(char **)LALInferenceGetVariable(runState->proposalArgs,"proposedVariableName");
                        char tmpname[MAX_STRLEN]="";

                        REAL8 s_gamma = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "s_gamma");
                        sprintf(tmpname,"%s_%s",name,ADAPTSUFFIX);
                        REAL8 *sigma = (REAL8 *)LALInferenceGetVariable(runState->proposalArgs,tmpname);

                        REAL8 priorMin, priorMax, dprior;

                        LALInferenceGetMinMaxPrior(runState->priorArgs, name, &priorMin, &priorMax);
                        dprior = priorMax - priorMin;

                        if(accepted == 1){
                                *sigma=*sigma+s_gamma*(dprior/100.0)*(1.0-targetAcceptance);
                        }else{
                                *sigma=*sigma-s_gamma*(dprior/100.0)*(targetAcceptance);
                        }

                        *sigma = (*sigma > dprior ? dprior : *sigma);
                        *sigma = (*sigma < DBL_MIN ? DBL_MIN : *sigma);

                        //printf("Adapting step size for %s to %lf\n", name, *sigma);

                        /* Make sure we don't do this again until we take another adaptable step.*/
                }
        }
        *adaptableStep = 0;
}

INT4 LALInferencePrintProposalTrackingHeader(FILE *fp,LALInferenceVariables *params) {
      fprintf(fp, "proposal\t");
      LALInferenceFprintParameterNonFixedHeaders(fp, params);
      LALInferenceFprintParameterNonFixedHeadersWithSuffix(fp, params, "p");
      fprintf(fp, "prop_ratio\taccepted\t");
      fprintf(fp, "\n");
      return 0;
}

void LALInferencePrintProposalTracking(FILE *fp, LALInferenceVariables *propArgs, LALInferenceVariables *theta, LALInferenceVariables *theta_prime){
  const char *currentProposalName;
  currentProposalName = *((const char **)LALInferenceGetVariable(propArgs, LALInferenceCurrentProposalName));
  REAL8 logProposalRatio = *(REAL8*) LALInferenceGetVariable(propArgs, "logProposalRatio");
  UINT4 accepted = *(UINT4 *) LALInferenceGetVariable(propArgs, "accepted");

  fprintf(fp, "%s\t", currentProposalName);
  LALInferencePrintSampleNonFixed(fp, theta);
  LALInferencePrintSampleNonFixed(fp, theta_prime);
  fprintf(fp, "%9.5f\t", exp(logProposalRatio));
  fprintf(fp, "%d\t", accepted);
  fprintf(fp, "\n");
  return;
}
