/*
 *  LALInferenceProposal.c:  Bayesian Followup, jump proposals.
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch, Will M. Farr, Ben Farr
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
#include <lal/LALInferenceInit.h>
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
#include <lal/LALInferenceClusteredKDE.h>
#include <lal/LALInferenceNestedSampler.h>
#include <alloca.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef enum {
  USES_DISTANCE_VARIABLE,
  USES_LOG_DISTANCE_VARIABLE
} DistanceParam;

const char *const cycleArrayName = "Proposal Cycle";
const char *const cycleArrayLengthName = "Proposal Cycle Length";
const char *const cycleArrayCounterName = "Proposal Cycle Counter";

/* Proposal Names */
const char *const nullProposalName = "NULL";
const char *const singleAdaptProposalName = "Single";
const char *const singleProposalName = "Single";
const char *const orbitalPhaseJumpName = "OrbitalPhase";
const char *const covarianceEigenvectorJumpName = "CovarianceEigenvector";
const char *const skyLocWanderJumpName = "SkyLocWander";
const char *const differentialEvolutionFullName = "DifferentialEvolutionFull";
const char *const differentialEvolutionIntrinsicName = "DifferentialEvolutionIntrinsic";
const char *const differentialEvolutionExtrinsicName = "DifferentialEvolutionExtrinsic";
const char *const ensembleStretchFullName = "EnsembleStretchFull";
const char *const ensembleStretchIntrinsicName = "EnsembleStretchIntrinsic";
const char *const ensembleStretchExtrinsicName = "EnsembleStretchExtrinsic";
const char *const drawApproxPriorName = "DrawApproxPrior";
const char *const drawFlatPriorName = "DrawFlatPrior";
const char *const skyReflectDetPlaneName = "SkyReflectDetPlane";
const char *const skyRingProposalName = "SkyRingProposal";
const char *const PSDFitJumpName = "PSDFitJump";
const char *const polarizationPhaseJumpName = "PolarizationPhase";
const char *const polarizationCorrPhaseJumpName = "CorrPolarizationPhase";
const char *const extrinsicParamProposalName = "ExtrinsicParamProposal";
const char *const frequencyBinJumpName = "FrequencyBin";
const char *const GlitchMorletJumpName = "glitchMorletJump";
const char *const GlitchMorletReverseJumpName = "glitchMorletReverseJump";
const char *const ensembleWalkFullName = "EnsembleWalkFull";
const char *const ensembleWalkIntrinsicName = "EnsembleWalkIntrinsic";
const char *const ensembleWalkExtrinsicName = "EnsembleWalkExtrinsic";
const char *const clusteredKDEProposalName = "ClusteredKDEProposal";
const char *const splineCalibrationProposalName = "SplineCalibration";
const char *const distanceLikelihoodProposalName = "DistanceLikelihood";

static const char *intrinsicNames[] = {"chirpmass", "q", "eta", "mass1", "mass2", "a_spin1", "a_spin2",
  "tilt_spin1", "tilt_spin2", "phi12", "phi_jl", "frequency", "quality", "duration","polar_angle", "phase", "polar_eccentricity","dchi0","dchi1","dchi2","dchi3","dchi4","dchi5","dchi5l","dchi6","dchi6l","dchi7","aPPE","alphaPPE","bPPE","betaPPE","betaStep","fStep","dxi1","dxi2","dxi3","dxi4","dxi5","dxi6","dalpha1","dalpha2","dalpha3","dalpha4","dalpha5","dbeta1","dbeta2","dbeta3","dsigma1","dsigma2","dsigma3","dsigma4",NULL};

static const char *extrinsicNames[] = {"rightascension", "declination", "cosalpha", "azimuth", "polarisation", "distance",
  "logdistance", "time", "costheta_jn", "t0", "theta","hrss", "loghrss", NULL};

static INT4 same_detector_location(LALDetector *d1, LALDetector *d2) {
    INT4 i;

    for (i = 0; i < 3; i++) {
        if (d1->location[i] != d2->location[i])
            return 0;
    }

    return 1;
}

static INT4 numDetectorsUniquePositions(LALInferenceIFOData *data) {
    INT4 nIFO = 0;
    INT4 nCollision = 0;
    LALInferenceIFOData *currentIFO = NULL;

    for (currentIFO = data; currentIFO; currentIFO = currentIFO->next) {
        LALInferenceIFOData *subsequentIFO = NULL;
        nIFO++;
        for (subsequentIFO = currentIFO->next; subsequentIFO; subsequentIFO = subsequentIFO->next) {
            if (same_detector_location(subsequentIFO->detector, currentIFO->detector)) {
                nCollision++;
                break;
            }
        }
    }

    return nIFO - nCollision;
}

LALInferenceProposal *LALInferenceInitProposal(LALInferenceProposalFunction func, const char *name)
{
  LALInferenceProposal *proposal = XLALCalloc(1,sizeof(LALInferenceProposal));
  proposal->func = func;
  proposal->proposed = 0;
  proposal->accepted = 0;
  strcpy(proposal->name, name);
  return proposal;
}


void LALInferenceRegisterProposal(LALInferenceVariables *propArgs, const char *name, INT4 *flag, ProcessParamsTable *command_line) {
    char offopt[VARNAME_MAX+15];
    char onopt[VARNAME_MAX+12];

    sprintf(offopt, "--proposal-no-%s", name);
    sprintf(onopt, "--proposal-%s", name);

    if (LALInferenceGetProcParamVal(command_line, offopt))
        *flag = 0;
    else if (LALInferenceGetProcParamVal(command_line, onopt))
        *flag = 1;

    LALInferenceAddINT4Variable(propArgs, name, *flag, LALINFERENCE_PARAM_FIXED);
}

void LALInferenceAddProposalToCycle(LALInferenceProposalCycle *cycle, LALInferenceProposal *prop, INT4 weight) {
    const char *fname = "LALInferenceAddProposalToCycle";
    INT4 i;

    /* Quit without doing anything if weight = 0. */
    if (weight == 0)
        return;

    cycle->order = XLALRealloc(cycle->order, (cycle->length + weight)*sizeof(INT4));
    if (cycle->order == NULL) {
        XLALError(fname, __FILE__, __LINE__, XLAL_ENOMEM);
        exit(1);
    }

    for (i = cycle->length; i < cycle->length + weight; i++) {
        cycle->order[i] = cycle->nProposals;
    }

    cycle->nProposals += 1;
    cycle->proposals = XLALRealloc(cycle->proposals, (cycle->nProposals)*sizeof(LALInferenceProposal));
    if (cycle->proposals == NULL) {
        XLALError(fname, __FILE__, __LINE__, XLAL_ENOMEM);
        exit(1);
    }
    cycle->proposals[cycle->nProposals-1] = prop;

    cycle->length += weight;
}



void LALInferenceRandomizeProposalCycle(LALInferenceProposalCycle *cycle, gsl_rng *rng) {
    INT4 i, j, temp;

    for (i = cycle->length - 1; i > 0; i--) {
        /* Fill in array from right to left, chosen randomly from remaining proposals. */
        j = gsl_rng_uniform_int(rng, i+1);

        temp = cycle->order[j];
        cycle->order[j] = cycle->order[i];
        cycle->order[i] = temp;
    }
}


REAL8 LALInferenceCyclicProposal(LALInferenceThreadState *thread,
                                 LALInferenceVariables *currentParams,
                                 LALInferenceVariables *proposedParams) {
    INT4 i = 0;
    LALInferenceProposalCycle *cycle=NULL;

    /* Must have cycle array and cycle array length in propArgs. */
    cycle = thread->cycle;
    if (cycle == NULL) {
        XLALError("LALInferenceCyclicProposal()",__FILE__,__LINE__,XLAL_FAILURE);
        exit(1);
    }

    if (cycle->counter >= cycle->length) {
        XLALError("LALInferenceCyclicProposal()",__FILE__,__LINE__,XLAL_FAILURE);
        exit(1);
    }

    /* One instance of each proposal object is stored in cycle->proposals.
        cycle->order is a list of elements to call from the proposals */

    REAL8 logPropRatio=-INFINITY;
    do
    {
      i = cycle->order[cycle->counter];
      logPropRatio = cycle->proposals[i]->func(thread, currentParams, proposedParams);
      strcpy(cycle->last_proposal_name, cycle->proposals[i]->name);
      cycle->counter = (cycle->counter + 1) % cycle->length;
    }
    /* Call proposals until one succeeds */
    while (proposedParams->head == NULL);

    return logPropRatio;
}

LALInferenceProposalCycle* LALInferenceInitProposalCycle(void) {
  LALInferenceProposalCycle *cycle = XLALCalloc(1,sizeof(LALInferenceProposalCycle));
  strcpy(cycle->last_proposal_name, nullProposalName);

  return cycle;
}

void LALInferenceDeleteProposalCycle(LALInferenceProposalCycle *cycle) {
    XLALFree(cycle->proposals);
    XLALFree(cycle->order);
}

LALInferenceVariables *LALInferenceParseProposalArgs(LALInferenceRunState *runState) {
    INT4 i;
    ProcessParamsTable *ppt;
    LALInferenceIFOData *ifo = runState->data;

    /* This will copy any existing arguments over from runState. I (John) don't think this should be necessary
     * as this function is used to initialise these arguments in the first place. */
    LALInferenceVariables *propArgs = XLALCalloc(1, sizeof(LALInferenceVariables));
    if(runState->proposalArgs && runState->proposalArgs->dimension>0) LALInferenceCopyVariables(runState->proposalArgs, propArgs);

    INT4 Nskip = 1;
    INT4 noise_only = 0;
    INT4 cyclic_reflective_kde = 0;

    /* Flags for proposals, initialized with the MCMC defaults */

    INT4 singleadapt = 1; /* Disabled for bug checking */
    INT4 psiphi = 1;
    INT4 ext_param = 1;
    INT4 skywander = 1;
    INT4 skyreflect = 1;
    INT4 drawprior = 1;
    INT4 covjump = 0;
    INT4 diffevo = 1;
    INT4 stretch = 0;
    INT4 walk = 0;
    INT4 skyring = 1;
    INT4 distance = 1;
    INT4 kde = 0;
    INT4 spline_cal = 0;
    INT4 psdfit = 0;
    INT4 glitchfit = 0;

    if (runState->algorithm == &LALInferenceNestedSamplingAlgorithm) {
        singleadapt = 0;
        psiphi = 0;
        ext_param = 0;
        skywander = 0;
        skyreflect = 0;
        drawprior = 0;
        covjump = 1;
        diffevo = 1;
        stretch = 1;
        walk = 1;
        skyring = 0;
        distance = 1;
        kde = 0;
        spline_cal = 0;
        psdfit = 0;
        glitchfit = 0;
    }
    if (LALInferenceCheckVariable(runState->algorithmParams, "LIB"))
      distance=0;

    ProcessParamsTable *command_line = runState->commandLine;

    INT4 verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--verbose"))
        verbose = 1;
    LALInferenceAddINT4Variable(propArgs, "verbose", verbose, LALINFERENCE_PARAM_FIXED);

    LIGOTimeGPS epoch = ifo->epoch;
    LALInferenceAddVariable(propArgs, "epoch", &(epoch), LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);

    /* Determine the number of iterations between each entry in the DE buffer */
    if (LALInferenceCheckVariable(runState->algorithmParams, "Nskip"))
        Nskip = LALInferenceGetINT4Variable(runState->algorithmParams, "Nskip");
    LALInferenceAddINT4Variable(propArgs, "Nskip", Nskip, LALINFERENCE_PARAM_FIXED);

    /* Count the number of IFOs and uniquely-located IFOs to decide which sky-related proposals to use */
    INT4 nDet = 0;
    ifo=runState->data;
    while (ifo) {
        nDet++;
        ifo = ifo->next;
    }
    LALInferenceAddINT4Variable(propArgs, "nDet", nDet, LALINFERENCE_PARAM_FIXED);

    INT4 nUniqueDet = numDetectorsUniquePositions(runState->data);
    LALInferenceAddINT4Variable(propArgs, "nUniqueDet", nUniqueDet, LALINFERENCE_PARAM_FIXED);

    LALDetector *detectors = XLALCalloc(nDet, sizeof(LALDetector));
    for (i=0,ifo=runState->data; i<nDet; i++,ifo=ifo->next)
        detectors[i] = *(ifo->detector);
    LALInferenceAddVariable(propArgs, "detectors", &detectors, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);

    char **ifo_names = XLALCalloc(nDet, sizeof(char*));
    for(ifo=runState->data,i=0;ifo;ifo=ifo->next,i++) {
        ifo_names[i] = XLALCalloc(DETNAMELEN, sizeof(char));
        strcpy(ifo_names[i], ifo->name);
    }
    LALInferenceAddVariable(propArgs, "detector_names", &ifo_names, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);

    INT4 marg_timephi = 0;
    if (LALInferenceGetProcParamVal(command_line, "--margtimephi"))
        marg_timephi = 1;

    INT4 marg_time = 0;
    if (marg_timephi || LALInferenceGetProcParamVal(command_line, "--margtime"))
        marg_time = 1;
    LALInferenceAddINT4Variable(propArgs, "marg_time", marg_time, LALINFERENCE_PARAM_FIXED);

    INT4 marg_phi = 0;
    if (marg_timephi || LALInferenceGetProcParamVal(command_line, "--margphi"))
        marg_phi = 1;
    LALInferenceAddINT4Variable(propArgs, "marg_phi", marg_phi, LALINFERENCE_PARAM_FIXED);

    INT4 analytic_test = 0;
    if (LALInferenceGetProcParamVal(command_line, "--correlatedGaussianLikelihood") ||
        LALInferenceGetProcParamVal(command_line, "--bimodalGaussianLikelihood") ||
        LALInferenceGetProcParamVal(command_line, "--rosenbrockLikelihood")) {
        analytic_test = 1;
        distance = 0;
    }
    LALInferenceAddINT4Variable(propArgs, "analytical_test", analytic_test, LALINFERENCE_PARAM_FIXED);

    INT4 skyframe = 1;
    if (LALInferenceGetProcParamVal(command_line, "--no-sky-frame"))
        skyframe = 0;

    INT4 noAdapt = 0;
    if (LALInferenceGetProcParamVal(command_line, "--no-adapt") ||
        LALInferenceGetProcParamVal(command_line, "--noiseonly"))
        noAdapt = 1;
    INT4 adapting = !noAdapt;
    LALInferenceAddINT4Variable(propArgs, "no_adapt", noAdapt, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddINT4Variable(propArgs, "adapting", adapting, LALINFERENCE_PARAM_LINEAR);

    INT4 tau = 5;
    ppt = LALInferenceGetProcParamVal(command_line, "--adapt-tau");
    if (ppt)
        tau = atof(ppt->value);
    LALInferenceAddINT4Variable(propArgs, "adaptTau", tau, LALINFERENCE_PARAM_FIXED);

    INT4 sampling_prior = 0;
    ppt = LALInferenceGetProcParamVal(command_line, "--zerologlike");
    if (ppt)
        sampling_prior = 1;
    LALInferenceAddINT4Variable(propArgs, "sampling_prior", sampling_prior, LALINFERENCE_PARAM_FIXED);

    if (LALInferenceGetProcParamVal(command_line, "--enable-spline-calibration"))
        spline_cal = 1;

    ppt = LALInferenceGetProcParamVal(command_line, "--psd-fit");
    if (!ppt)
        ppt = LALInferenceGetProcParamVal(command_line, "--psdFit");
    if (ppt)
        psdfit = 1;

    if (LALInferenceGetProcParamVal(command_line, "--glitch-fit"))
        glitchfit = 1;

    /* Convenience option for turning off ensemble moves */
    if (LALInferenceGetProcParamVal(command_line, "--proposal-no-ensemble")) {
        stretch = 0;
        walk = 0;
    }

    /* Check if imposing cyclic reflective bounds */
    if (LALInferenceGetProcParamVal(runState->commandLine, "--cyclic-reflective-kde"))
        cyclic_reflective_kde = 1;
    LALInferenceAddINT4Variable(propArgs, "cyclic_reflective_kde", cyclic_reflective_kde, LALINFERENCE_PARAM_FIXED);

    if (LALInferenceGetProcParamVal(command_line, "--noiseonly"))
        noise_only = 1;
    LALInferenceAddINT4Variable(propArgs, "noiseonly", noise_only, LALINFERENCE_PARAM_FIXED);

    /* Turn off signal proposals if no signal is in the model */
    if (noise_only) {
        singleadapt = 0;
        psiphi = 0;
        ext_param = 0;
        skywander = 0;
        skyreflect = 0;
        drawprior = 0;
        covjump = 0;
        diffevo = 0;
        stretch = 0;
        walk = 0;
        distance = 0;
        skyring = 0;
        spline_cal = 0;
    }

    /* Turn off phi-related proposals if marginalizing over phi in likelihood */
    if (marg_phi) {
        psiphi = 0;
    }

    /* Disable proposals that won't work with the current number of unique detectors */
    if (nUniqueDet < 2) {
        skyring = 0;
    }

    if (nUniqueDet != 3) {
        skyreflect = 0;
    }

    if (nUniqueDet >= 3) {
        ext_param = 0;
    }

    /* Turn off ra-dec related proposals when using the sky-frame coordinate system */
    if (skyframe) {
        ext_param = 0;
        skywander = 0;
        skyreflect = 0;
        skyring = 0;
    }

    /* Register all proposal functions, check for explicit command-line requests */
    LALInferenceRegisterProposal(propArgs, "single-adapt", &singleadapt, command_line);
    LALInferenceRegisterProposal(propArgs, "psiphi", &psiphi, command_line);
    LALInferenceRegisterProposal(propArgs, "extrinsic-param", &ext_param, command_line);
    LALInferenceRegisterProposal(propArgs, "sky-wander", &skywander, command_line);
    LALInferenceRegisterProposal(propArgs, "sky-reflect", &skyreflect, command_line);
    LALInferenceRegisterProposal(propArgs, "draw-prior", &drawprior, command_line);
    LALInferenceRegisterProposal(propArgs, "eigenvectors", &covjump, command_line);
    LALInferenceRegisterProposal(propArgs, "differential-evolution", &diffevo, command_line);
    LALInferenceRegisterProposal(propArgs, "ensemble-stretch", &stretch, command_line);
    LALInferenceRegisterProposal(propArgs, "ensemble-walk", &walk, command_line);
    LALInferenceRegisterProposal(propArgs, "sky-ring", &skyring, command_line);
    LALInferenceRegisterProposal(propArgs, "distance", &distance, command_line);
    LALInferenceRegisterProposal(propArgs, "kde", &kde, command_line);
    LALInferenceRegisterProposal(propArgs, "spline_cal", &spline_cal, command_line);
    LALInferenceRegisterProposal(propArgs, "psdfit", &psdfit, command_line);
    LALInferenceRegisterProposal(propArgs, "glitchfit", &glitchfit, command_line);

    /* Setup adaptive proposals */
    if (singleadapt){
      LALInferenceModel *model = LALInferenceInitCBCModel(runState);
      LALInferenceSetupAdaptiveProposals(propArgs, model->params);
      XLALFree(model);
    }
    /* Setup buffer now since threads aren't accessible to the main setup function */
    if (diffevo || stretch || walk) {
        for (i=0; i<runState->nthreads; i++)
            LALInferenceSetupDifferentialEvolutionProposal(runState->threads[i]);
    }

    /* Setup now since we need access to the data */
    if (glitchfit)
        LALInferenceSetupGlitchProposal(runState->data, propArgs);

    return propArgs;
}


LALInferenceProposalCycle* LALInferenceSetupDefaultInspiralProposalCycle(LALInferenceVariables *propArgs) {
    LALInferenceProposal *prop;

    const INT4 BIGWEIGHT = 20;
    const INT4 SMALLWEIGHT = 5;
    const INT4 TINYWEIGHT = 1;

    LALInferenceProposalCycle *cycle = XLALCalloc(1, sizeof(LALInferenceProposalCycle));

    if (LALInferenceGetINT4Variable(propArgs, "single-adapt")) {
        prop = LALInferenceInitProposal(&LALInferenceSingleAdaptProposal, singleAdaptProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "psiphi")) {
        prop = LALInferenceInitProposal(&LALInferencePolarizationPhaseJump, polarizationPhaseJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, TINYWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "extrinsic-param")) {
        prop = LALInferenceInitProposal(&LALInferenceExtrinsicParamProposal, extrinsicParamProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "sky-wander")) {
        prop = LALInferenceInitProposal(&LALInferenceSkyLocWanderJump, skyLocWanderJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "sky-reflect")) {
        prop = LALInferenceInitProposal(&LALInferenceSkyReflectDetPlane, skyReflectDetPlaneName);
        LALInferenceAddProposalToCycle(cycle, prop, TINYWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "draw-prior")) {
        prop = LALInferenceInitProposal(&LALInferenceDrawApproxPrior, drawApproxPriorName);
        LALInferenceAddProposalToCycle(cycle, prop, TINYWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "eigenvectors")) {
        prop = LALInferenceInitProposal(&LALInferenceCovarianceEigenvectorJump, covarianceEigenvectorJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "differential-evolution")) {
        prop = LALInferenceInitProposal(&LALInferenceDifferentialEvolutionFull, differentialEvolutionFullName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceDifferentialEvolutionIntrinsic, differentialEvolutionIntrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceDifferentialEvolutionExtrinsic, differentialEvolutionExtrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "ensemble-stretch")) {
        prop = LALInferenceInitProposal(&LALInferenceEnsembleStretchFull, ensembleStretchFullName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceEnsembleStretchIntrinsic, ensembleStretchIntrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceEnsembleStretchExtrinsic, ensembleStretchExtrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "ensemble-walk")) {
        prop = LALInferenceInitProposal(&LALInferenceEnsembleWalkFull, ensembleWalkFullName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceEnsembleWalkIntrinsic, ensembleWalkIntrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceEnsembleWalkExtrinsic, ensembleWalkExtrinsicName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "sky-ring")) {
        prop = LALInferenceInitProposal(&LALInferenceSkyRingProposal, skyRingProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    /* Distance proposal */
    if (LALInferenceGetINT4Variable(propArgs, "distance")) {
        prop = LALInferenceInitProposal(&LALInferenceDistanceLikelihoodProposal, distanceLikelihoodProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "kde")) {
        prop = LALInferenceInitProposal(&LALInferenceClusteredKDEProposal, clusteredKDEProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, BIGWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "spline_cal")) {
        prop = LALInferenceInitProposal(&LALInferenceSplineCalibrationProposal, splineCalibrationProposalName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "psdfit")) {
        prop = LALInferenceInitProposal(&LALInferencePSDFitJump, PSDFitJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    if (LALInferenceGetINT4Variable(propArgs, "glitchfit")) {
        prop = LALInferenceInitProposal(&LALInferenceGlitchMorletProposal, GlitchMorletJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);

        prop = LALInferenceInitProposal(&LALInferenceGlitchMorletReverseJump, GlitchMorletReverseJumpName);
        LALInferenceAddProposalToCycle(cycle, prop, SMALLWEIGHT);
    }

    return cycle;
}


REAL8 LALInferenceSingleAdaptProposal(LALInferenceThreadState *thread,
                                      LALInferenceVariables *currentParams,
                                      LALInferenceVariables *proposedParams) {
    INT4 dim, varNr;
    REAL8 logPropRatio, sqrttemp, sigma;
    char tmpname[MAX_STRLEN] = "";
    LALInferenceVariableItem *param = NULL;

    LALInferenceCopyVariables(currentParams, proposedParams);
    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng *rng = thread->GSLrandom;

    if (!LALInferenceGetINT4Variable(args, "no_adapt")) {
        if (!LALInferenceCheckVariable(args, "adapting"))
            LALInferenceSetupAdaptiveProposals(args, currentParams);

        sqrttemp = sqrt(thread->temperature);
        dim = proposedParams->dimension;

        do {
            varNr = 1 + gsl_rng_uniform_int(rng, dim);
            param = LALInferenceGetItemNr(proposedParams, varNr);
        } while (!LALInferenceCheckVariableNonFixed(proposedParams, param->name) || param->type != LALINFERENCE_REAL8_t);

        if (param->type != LALINFERENCE_REAL8_t) {
            fprintf(stderr, "Attempting to set non-REAL8 parameter with numerical sigma (in %s, %d)\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        sprintf(tmpname,"%s_%s",param->name,ADAPTSUFFIX);
        if (!LALInferenceCheckVariable(thread->proposalArgs, tmpname)) {
            fprintf(stderr, "Attempting to draw single-parameter jump for %s but cannot find sigma!\nError in %s, line %d.\n",
                    param->name,__FILE__, __LINE__);
            exit(1);
        }

        sigma = LALInferenceGetREAL8Variable(thread->proposalArgs, tmpname);

        /* Save the name of the proposed variable */
        LALInferenceAddstringVariable(args, "proposedVariableName", param->name, LALINFERENCE_PARAM_OUTPUT);

        *((REAL8 *)param->value) += gsl_ran_ugaussian(rng) * sigma * sqrttemp;

        LALInferenceCyclicReflectiveBound(proposedParams, thread->priorArgs);

        /* Set the log of the proposal ratio to zero, since this is a
        symmetric proposal. */
        logPropRatio = 0.0;

        INT4 as = 1;
        LALInferenceSetVariable(args, "adaptableStep", &as);

    } else {
        /* We are not adaptive, or for some reason don't have a sigma
           vector---fall back on old proposal. */
        logPropRatio = LALInferenceSingleProposal(thread, currentParams, proposedParams);
    }

    return logPropRatio;
}

REAL8 LALInferenceSingleProposal(LALInferenceThreadState *thread,
                                 LALInferenceVariables *currentParams,
                                 LALInferenceVariables *proposedParams) {
    LALInferenceVariableItem *param=NULL;
    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng * GSLrandom = thread->GSLrandom;
    REAL8 sigma, big_sigma;
    INT4 dim, varNr;

    LALInferenceCopyVariables(currentParams, proposedParams);

    sigma = 0.1 * sqrt(thread->temperature); /* Adapt step to temperature. */
    big_sigma = 1.0;

    if (gsl_ran_ugaussian(GSLrandom) < 1.0e-3)
        big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in a parameter
    if (gsl_ran_ugaussian(GSLrandom) < 1.0e-4)
        big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in a parameter

    dim = proposedParams->dimension;

    do {
        varNr = 1 + gsl_rng_uniform_int(GSLrandom, dim);
        param = LALInferenceGetItemNr(proposedParams, varNr);
    } while (!LALInferenceCheckVariableNonFixed(proposedParams, param->name) || param->type != LALINFERENCE_REAL8_t);

    /* Scale jumps proposal appropriately for prior sampling */
    if (LALInferenceGetINT4Variable(args, "sampling_prior")) {
        if (!strcmp(param->name, "eta")) {
            sigma = 0.02;
        } else if (!strcmp(param->name, "q")) {
            sigma = 0.08;
        } else if (!strcmp(param->name, "chirpmass")) {
            sigma = 1.0;
        } else if (!strcmp(param->name, "time")) {
            sigma = 0.02;
		} else if (!strcmp(param->name, "t0")) {
		    sigma = 0.02;
        } else if (!strcmp(param->name, "phase")) {
            sigma = 0.6;
        } else if (!strcmp(param->name, "distance")) {
            sigma = 10.0;
        } else if (!strcmp(param->name, "declination")) {
            sigma = 0.3;
		} else if (!strcmp(param->name, "azimuth")) {
			sigma = 0.6;
		} else if (!strcmp(param->name, "cosalpha")) {
			sigma = 0.1;
		} else if (!strcmp(param->name, "rightascension")) {
            sigma = 0.6;
        } else if (!strcmp(param->name, "polarisation")) {
            sigma = 0.6;
        } else if (!strcmp(param->name,"costheta_jn")) {
            sigma = 0.3;
        } else if (!strcmp(param->name, "a_spin1")) {
            sigma = 0.1;
        } else if (!strcmp(param->name, "a_spin2")) {
            sigma = 0.1;
        } else {
            fprintf(stderr, "Could not find parameter %s!", param->name);
            exit(1);
        }

        *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*sigma;
    } else {
        if (!strcmp(param->name,"eta") || !strcmp(param->name,"q") || !strcmp(param->name,"time") || !strcmp(param->name,"t0") || !strcmp(param->name,"a_spin2") || !strcmp(param->name,"a_spin1")){
            *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
        } else if (!strcmp(param->name,"polarisation") || !strcmp(param->name,"phase") || !strcmp(param->name,"costheta_jn")){
            *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
        } else {
            *(REAL8 *)param->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
        }
    }

    LALInferenceCyclicReflectiveBound(proposedParams, thread->priorArgs);

    /* Symmetric Proposal. */
    REAL8 logPropRatio = 0.0;

    return logPropRatio;
}


REAL8 LALInferenceCovarianceEigenvectorJump(LALInferenceThreadState *thread,
                                            LALInferenceVariables *currentParams,
                                            LALInferenceVariables *proposedParams) {
    LALInferenceVariableItem *proposeIterator;
    REAL8Vector *eigenvalues;
    gsl_matrix *eigenvectors;
    REAL8 jumpSize, tmp, inc;
    REAL8 logPropRatio = 0.0;
    INT4 N, i, j;

    LALInferenceCopyVariables(currentParams, proposedParams);

    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng *rng = thread->GSLrandom;

    eigenvalues = LALInferenceGetREAL8VectorVariable(args, "covarianceEigenvalues");
    eigenvectors = LALInferenceGetgslMatrixVariable(args, "covarianceEigenvectors");

    N = eigenvalues->length;

    i = gsl_rng_uniform_int(rng, N);
    jumpSize = sqrt(thread->temperature * eigenvalues->data[i]) * gsl_ran_ugaussian(rng);

    j = 0;
    proposeIterator = proposedParams->head;
    if (proposeIterator == NULL) {
        fprintf(stderr, "Bad proposed params in %s, line %d\n",
                __FILE__, __LINE__);
        exit(1);
    }

    do {
        if (LALInferenceCheckVariableNonFixed(proposedParams, proposeIterator->name) &&
            proposeIterator->type==LALINFERENCE_REAL8_t) {
            tmp = LALInferenceGetREAL8Variable(proposedParams, proposeIterator->name);
            inc = jumpSize * gsl_matrix_get(eigenvectors, j, i);

            tmp += inc;

            LALInferenceSetVariable(proposedParams, proposeIterator->name, &tmp);

            j++;
        }
    } while ((proposeIterator = proposeIterator->next) != NULL && j < N);

    return logPropRatio;
}

REAL8 LALInferenceSkyLocWanderJump(LALInferenceThreadState *thread,
                                   LALInferenceVariables *currentParams,
                                   LALInferenceVariables *proposedParams) {
    REAL8 sigma;
    REAL8 jumpX, jumpY;
    REAL8 RA, DEC;
    REAL8 newRA, newDEC;
    REAL8 one_deg = 1.0 / (2.0*M_PI);
    REAL8 logPropRatio = 0.0;

    LALInferenceCopyVariables(currentParams, proposedParams);

    gsl_rng *rng = thread->GSLrandom;

    sigma = sqrt(thread->temperature) * one_deg;
    jumpX = sigma * gsl_ran_ugaussian(rng);
    jumpY = sigma * gsl_ran_ugaussian(rng);

    RA = LALInferenceGetREAL8Variable(proposedParams, "rightascension");
    DEC = LALInferenceGetREAL8Variable(proposedParams, "declination");

    newRA = RA + jumpX;
    newDEC = DEC + jumpY;

    LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
    LALInferenceSetVariable(proposedParams, "declination", &newDEC);


    return logPropRatio;
}

REAL8 LALInferenceDifferentialEvolutionFull(LALInferenceThreadState *thread,
                                            LALInferenceVariables *currentParams,
                                            LALInferenceVariables *proposedParams) {
    return(LALInferenceDifferentialEvolutionNames(thread, currentParams, proposedParams, NULL));
}

REAL8 LALInferenceEnsembleStretchFull(LALInferenceThreadState *thread,
                                      LALInferenceVariables *currentParams,
                                      LALInferenceVariables *proposedParams) {
    return(LALInferenceEnsembleStretchNames(thread, currentParams, proposedParams, NULL));
}


REAL8 LALInferenceEnsembleStretchIntrinsic(LALInferenceThreadState *thread,
                                           LALInferenceVariables *currentParams,
                                           LALInferenceVariables *proposedParams) {

    return(LALInferenceEnsembleStretchNames(thread, currentParams, proposedParams, intrinsicNames));

}

REAL8 LALInferenceEnsembleStretchExtrinsic(LALInferenceThreadState *thread,
                                           LALInferenceVariables *currentParams,
                                           LALInferenceVariables *proposedParams) {
    REAL8 logPropRatio;

    logPropRatio = LALInferenceEnsembleStretchNames(thread, currentParams, proposedParams, extrinsicNames);

    return logPropRatio;
}

/* This jump uses the current sample 'A' and another randomly
 * drawn 'B' from the ensemble of live points, and proposes
 * C = B+Z(A-B) where Z is a scale factor */
REAL8 LALInferenceEnsembleStretchNames(LALInferenceThreadState *thread,
                                       LALInferenceVariables *currentParams,
                                       LALInferenceVariables *proposedParams,
                                       const char **names) {
    size_t i, N, Ndim, nPts;
    REAL8 logPropRatio;
    REAL8 maxScale, Y, logmax, X, scale;
    REAL8 cur, other, x;
    LALInferenceVariableItem *item;
    LALInferenceVariables **dePts;
    LALInferenceVariables *ptI;

    N = LALInferenceGetVariableDimension(currentParams) + 1; /* More names than we need. */
    const char* local_names[N];
    if (names == NULL) {
        names = local_names;

        item = currentParams->head;
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

    Ndim = 0;
    for(Ndim=0,i=0; names[i] != NULL; i++ ) {
        if(LALInferenceCheckVariableNonFixed(currentParams,names[i]))
        Ndim++;
    }

    LALInferenceCopyVariables(currentParams, proposedParams);

    Ndim = 0;
    for(Ndim=0, i=0; names[i] != NULL; i++ ) {
        if (LALInferenceCheckVariableNonFixed(proposedParams, names[i]))
            Ndim++;
    }

    dePts = thread->differentialPoints;
    nPts = thread->differentialPointsLength;

    if (dePts == NULL || nPts <= 1) {
        logPropRatio = 0.0;
        return logPropRatio; /* Quit now, since we don't have any points to use. */
    }

    /* Choose a different sample */
    do {
        i = gsl_rng_uniform_int(thread->GSLrandom, nPts);
    } while (!LALInferenceCompareVariables(currentParams, dePts[i]));

    ptI = dePts[i];

    /* Scale z is chosen according to be symmetric under z -> 1/z */
    /* so p(x) \propto 1/z between 1/a and a */

    /* TUNABLE PARAMETER (a), must be >1. Larger value -> smaller acceptance */
    maxScale=3.0;

    /* Draw sample between 1/max and max */
    Y = gsl_rng_uniform(thread->GSLrandom);
    logmax = log(maxScale);
    X = 2.0*logmax*Y - logmax;
    scale = exp(X);

    for (i = 0; names[i] != NULL; i++) {
        /* Ignore variable if it's not in each of the params. */
        if (LALInferenceCheckVariableNonFixed(proposedParams, names[i]) &&
            LALInferenceCheckVariableNonFixed(ptI, names[i])) {
                cur = LALInferenceGetREAL8Variable(proposedParams, names[i]);
                other= LALInferenceGetREAL8Variable(ptI, names[i]);
                x = other + scale*(cur-other);

                LALInferenceSetVariable(proposedParams, names[i], &x);
        }
    }

    if (scale < maxScale && scale > (1.0/maxScale))
        logPropRatio = log(scale)*((REAL8)Ndim);
    else
        logPropRatio = -INFINITY;

    return logPropRatio;
}


REAL8 LALInferenceEnsembleWalkFull(LALInferenceThreadState *thread,
                                   LALInferenceVariables *currentParams,
                                   LALInferenceVariables *proposedParams) {
    return(LALInferenceEnsembleWalkNames(thread, currentParams, proposedParams, NULL));
}


REAL8 LALInferenceEnsembleWalkIntrinsic(LALInferenceThreadState *thread,
                                        LALInferenceVariables *currentParams,
                                        LALInferenceVariables *proposedParams) {

    return(LALInferenceEnsembleWalkNames(thread, currentParams, proposedParams, intrinsicNames));

}

REAL8 LALInferenceEnsembleWalkExtrinsic(LALInferenceThreadState *thread,
                                        LALInferenceVariables *currentParams,
                                        LALInferenceVariables *proposedParams) {
    return(LALInferenceEnsembleWalkNames(thread, currentParams, proposedParams, extrinsicNames));
}


REAL8 LALInferenceEnsembleWalkNames(LALInferenceThreadState *thread,
                                    LALInferenceVariables *currentParams,
                                    LALInferenceVariables *proposedParams,
                                    const char **names) {
  size_t i;
  LALInferenceVariableItem *item;
  REAL8 logPropRatio = 0.0;

  size_t N = LALInferenceGetVariableDimension(currentParams) + 1; /* More names than we need. */
  const char* local_names[N];
  if (names == NULL) {
    names = local_names;

    item = currentParams->head;
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
    if(LALInferenceCheckVariableNonFixed(currentParams,names[i]))
      Ndim++;
  }

  LALInferenceVariables **pointsPool = thread->differentialPoints;
  size_t k=0;
  size_t sample_size=3;

  LALInferenceVariables **dePts = thread->differentialPoints;
  size_t nPts = thread->differentialPointsLength;

  if (dePts == NULL || nPts <= 1) {
    logPropRatio = 0.0;
    return logPropRatio; /* Quit now, since we don't have any points to use. */
  }

  LALInferenceCopyVariables(currentParams, proposedParams);

  UINT4 indices[sample_size];
  UINT4 all_indices[nPts];

  for (i=0;i<nPts;i++) all_indices[i]=i;
  gsl_ran_choose(thread->GSLrandom,indices, sample_size, all_indices, nPts, sizeof(UINT4));

  double w=0.0;
  double univariate_normals[sample_size];
  for(i=0;i<sample_size;i++) univariate_normals[i] = gsl_ran_ugaussian(thread->GSLrandom);

  /* Note: Simplified this loop on master 2015-08-12, take this version when rebasing */
  for(k=0;names[k]!=NULL;k++)
  {
    if(!LALInferenceCheckVariableNonFixed(proposedParams,names[k]) || LALInferenceGetVariableType(proposedParams,names[k])!=LALINFERENCE_REAL8_t) continue;
    REAL8 centre_of_mass=0.0;
    /* Compute centre of mass */
    for(i=0;i<sample_size;i++)
    {
      centre_of_mass+=LALInferenceGetREAL8Variable(pointsPool[indices[i]],names[k])/((REAL8)sample_size);
    }
    /* Compute offset */
    for(i=0,w=0.0;i<sample_size;i++)
    {
      w+= univariate_normals[i] * (LALInferenceGetREAL8Variable(pointsPool[indices[i]],names[k]) - centre_of_mass);
    }
    REAL8 tmp = LALInferenceGetREAL8Variable(proposedParams,names[k]) + w;
    LALInferenceSetVariable(proposedParams,names[k],&tmp);
  }

  logPropRatio = 0.0;

  return logPropRatio;
}

REAL8 LALInferenceDifferentialEvolutionNames(LALInferenceThreadState *thread,
                                    LALInferenceVariables *currentParams,
                                    LALInferenceVariables *proposedParams,
                                    const char **names) {
    size_t i, j, N, Ndim, nPts;
    LALInferenceVariableItem *item;
    LALInferenceVariables **dePts;
    LALInferenceVariables *ptI, *ptJ;
    REAL8 logPropRatio = 0.0;
    REAL8 scale, x;


    gsl_rng *rng = thread->GSLrandom;

    if (names == NULL) {
        N = LALInferenceGetVariableDimension(currentParams) + 1; /* More names than we need. */
        names = alloca(N * sizeof(char *)); /* Hope we have alloca---saves
                                               having to deallocate after
                                               proposal. */

        item = currentParams->head;
        i = 0;
        while (item != NULL) {
            if (LALInferenceCheckVariableNonFixed(currentParams, item->name) && item->type==LALINFERENCE_REAL8_t ) {
                names[i] = item->name;
                i++;
            }

            item = item->next;
        }
        names[i]=NULL; /* Terminate */
    }


    Ndim = 0;
    for (Ndim=0, i=0; names[i] != NULL; i++ ) {
        if (LALInferenceCheckVariableNonFixed(currentParams, names[i]))
            Ndim++;
    }

    dePts = thread->differentialPoints;
    nPts = thread->differentialPointsLength;

    if (dePts == NULL || nPts <= 1)
        return logPropRatio; /* Quit now, since we don't have any points to use. */

    LALInferenceCopyVariables(currentParams, proposedParams);


    i = gsl_rng_uniform_int(rng, nPts);
    do {
        j = gsl_rng_uniform_int(rng, nPts);
    } while (j == i);

    ptI = dePts[i];
    ptJ = dePts[j];

    const REAL8 modeHoppingFrac = 0.5;
    /* Some fraction of the time, we do a "mode hopping" jump,
       where we jump exactly along the difference vector. */
    if (gsl_rng_uniform(rng) < modeHoppingFrac) {
        scale = 1.0;
    } else {
        /* Otherwise scale is chosen uniform in log between 0.1 and 10 times the
        desired jump size. */
        scale = 2.38/sqrt(Ndim) * exp(log(0.1) + log(100.0) * gsl_rng_uniform(rng));
    }

    for (i = 0; names[i] != NULL; i++) {
        if (!LALInferenceCheckVariableNonFixed(currentParams, names[i]) ||
            !LALInferenceCheckVariable(ptJ, names[i]) ||
            !LALInferenceCheckVariable(ptI, names[i])) {
        /* Ignore variable if it's not in each of the params. */
        } else {
            x = LALInferenceGetREAL8Variable(currentParams, names[i]);
            x += scale * LALInferenceGetREAL8Variable(ptJ, names[i]);
            x -= scale * LALInferenceGetREAL8Variable(ptI, names[i]);
            LALInferenceSetVariable(proposedParams, names[i], &x);
        }
    }

    return logPropRatio;
}

REAL8 LALInferenceDifferentialEvolutionIntrinsic(LALInferenceThreadState *thread,
                                                 LALInferenceVariables *currentParams,
                                                 LALInferenceVariables *proposedParams) {

    return(LALInferenceDifferentialEvolutionNames(thread, currentParams, proposedParams, intrinsicNames));
}

REAL8 LALInferenceDifferentialEvolutionExtrinsic(LALInferenceThreadState *thread,
                                                 LALInferenceVariables *currentParams,
                                                 LALInferenceVariables *proposedParams) {
    return(LALInferenceDifferentialEvolutionNames(thread, currentParams, proposedParams, extrinsicNames));
}

static REAL8 draw_distance(LALInferenceThreadState *thread) {
    REAL8 dmin, dmax, x;

    LALInferenceGetMinMaxPrior(thread->priorArgs, "distance", &dmin, &dmax);

    x = gsl_rng_uniform(thread->GSLrandom);

    return cbrt(x*(dmax*dmax*dmax - dmin*dmin*dmin) + dmin*dmin*dmin);
}

static REAL8 draw_logdistance(LALInferenceThreadState *thread) {
    REAL8 logdmin, logdmax;

    LALInferenceGetMinMaxPrior(thread->priorArgs, "logdistance", &logdmin, &logdmax);

    REAL8 dmin=exp(logdmin);
    REAL8 dmax=exp(logdmax);

    REAL8 x = gsl_rng_uniform(thread->GSLrandom);

    return log(cbrt(x*(dmax*dmax*dmax - dmin*dmin*dmin) + dmin*dmin*dmin));
}

static REAL8 draw_colatitude(LALInferenceThreadState *thread, const char *name) {
    REAL8 min, max, x;

    LALInferenceGetMinMaxPrior(thread->priorArgs, name, &min, &max);

    x = gsl_rng_uniform(thread->GSLrandom);

    return acos(cos(min) - x*(cos(min) - cos(max)));
}

static REAL8 draw_dec(LALInferenceThreadState *thread) {
    REAL8 min, max, x;

    LALInferenceGetMinMaxPrior(thread->priorArgs, "declination", &min, &max);

    x = gsl_rng_uniform(thread->GSLrandom);

    return asin(x*(sin(max) - sin(min)) + sin(min));
}

static REAL8 draw_flat(LALInferenceThreadState *thread, const char *name) {
    REAL8 min, max, x;

    LALInferenceGetMinMaxPrior(thread->priorArgs, name, &min, &max);

    x = gsl_rng_uniform(thread->GSLrandom);

    return min + x*(max - min);
}

static REAL8 draw_chirp(LALInferenceThreadState *thread) {
    REAL8 min, max, delta;
    REAL8 mMin56, mMax56, u;

    LALInferenceGetMinMaxPrior(thread->priorArgs, "chirpmass", &min, &max);

    mMin56 = pow(min, 5.0/6.0);
    mMax56 = pow(max, 5.0/6.0);

    delta = 1.0/mMin56 - 1.0/mMax56;

    u = delta*gsl_rng_uniform(thread->GSLrandom);

    return pow(1.0/(1.0/mMin56 - u), 6.0/5.0);
}

static REAL8 approxLogPrior(LALInferenceVariables *params) {
    REAL8 logP = 0.0;

    if(LALInferenceCheckVariable(params, "chirpmass")) {
        REAL8 Mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
        logP += -11.0/6.0*log(Mc);
    }

    /* Flat in time, ra, psi, phi. */

    if(LALInferenceCheckVariable(params,"logdistance"))
      logP += 3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
    else if(LALInferenceCheckVariable(params,"distance"))
      logP += 2.0*log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));

    if(LALInferenceCheckVariable(params,"declination"))
      logP += log(cos(*(REAL8 *)LALInferenceGetVariable(params, "declination")));

    if (LALInferenceCheckVariable(params, "tilt_spin1")) {
      logP += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params, "tilt_spin1"))));
    }

    if (LALInferenceCheckVariable(params, "tilt_spin2")) {
      logP += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params, "tilt_spin2"))));
    }

    return logP;
}

/* WARNING: If you add any non-flat draws to this proposal, you MUST
   update the above approxLogPrior function with the corresponding
   density.  The proposal ratio is calculated using approxLogPrior, so
   non-flat proposals that do not have a corresponding density term in
   approxLogPrior will result in a failure to sample from the
   posterior density! */
REAL8 LALInferenceDrawApproxPrior(LALInferenceThreadState *thread,
                                  LALInferenceVariables *currentParams,
                                  LALInferenceVariables *proposedParams) {
    REAL8 tmp = 0.0;
    INT4 analytic_test, i;
    REAL8 logBackwardJump;
    REAL8 logPropRatio;
    LALInferenceVariableItem *ptr;

    LALInferenceCopyVariables(currentParams, proposedParams);

    const char *flat_params[] = {"q", "eta", "t0", "azimuth", "cosalpha", "time", "phase", "polarisation",
                                 "rightascension", "costheta_jn", "phi_jl",
                                 "phi12", "a_spin1", "a_spin2", NULL};

    LALInferenceVariables *args = thread->proposalArgs;

    analytic_test = LALInferenceGetINT4Variable(args, "analytical_test");

    if (analytic_test) {
        ptr = currentParams->head;
        while (ptr!=NULL) {
            if (LALInferenceCheckVariableNonFixed(currentParams, ptr->name)) {
                tmp = draw_flat(thread, ptr->name);
                LALInferenceSetVariable(proposedParams, ptr->name, &tmp);
            }
            ptr=ptr->next;
        }
    } else {
        logBackwardJump = approxLogPrior(currentParams);

        for (i = 0; flat_params[i] != NULL; i++) {
            if (LALInferenceCheckVariableNonFixed(proposedParams, flat_params[i])) {
                REAL8 val = draw_flat(thread, flat_params[i]);
                LALInferenceSetVariable(proposedParams, flat_params[i], &val);
            }
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "chirpmass")) {
            REAL8 Mc = draw_chirp(thread);
            LALInferenceSetVariable(proposedParams, "chirpmass", &Mc);
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "logdistance")) {
            REAL8 logdist = draw_logdistance(thread);
            LALInferenceSetVariable(proposedParams, "logdistance", &logdist);
        } else if (LALInferenceCheckVariableNonFixed(proposedParams, "distance")) {
            REAL8 dist = draw_distance(thread);
            LALInferenceSetVariable(proposedParams, "distance", &dist);
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "declination")) {
            REAL8 dec = draw_dec(thread);
            LALInferenceSetVariable(proposedParams, "declination", &dec);
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "tilt_spin1")) {
            REAL8 tilt1 = draw_colatitude(thread, "tilt_spin1");
            LALInferenceSetVariable(proposedParams, "tilt_spin1", &tilt1);
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "tilt_spin2")) {
            REAL8 tilt2 = draw_colatitude(thread, "tilt_spin2");
            LALInferenceSetVariable(proposedParams, "tilt_spin2", &tilt2);
        }

        if (LALInferenceCheckVariableNonFixed(proposedParams, "psdscale")) {
            REAL8 x, min, max;
            INT4 j;

            min=0.10;
            max=10.0;

            gsl_matrix *eta = LALInferenceGetgslMatrixVariable(proposedParams, "psdscale");

            for(i=0; i<(INT8)eta->size1; i++) {
                for(j=0; j<(INT8)eta->size2; j++) {
                    x = min + gsl_rng_uniform(thread->GSLrandom) * (max - min);
                    gsl_matrix_set(eta, i, j, x);
                }
            }
        }//end if(psdscale)
    }

    if (analytic_test) {
        /* Flat in every variable means uniform jump probability. */
        logPropRatio = 0.0;
    } else {
        logPropRatio = logBackwardJump - approxLogPrior(proposedParams);
    }

    return logPropRatio;
}

REAL8 LALInferenceDrawFlatPrior(LALInferenceThreadState *threadState, LALInferenceVariables *currentParams, LALInferenceVariables *proposedParams) {
  REAL8 logPropRatio = 0., tmp = 0.;

  LALInferenceCopyVariables(currentParams, proposedParams);
  LALInferenceVariableItem *ptr = currentParams->head;

  while(ptr!=NULL) {
    if(ptr->vary != LALINFERENCE_PARAM_FIXED && LALInferenceCheckMinMaxPrior(threadState->priorArgs, ptr->name ) ) {
      tmp = draw_flat(threadState, ptr->name);
      LALInferenceSetVariable(proposedParams, ptr->name, &tmp);
    }
    ptr=ptr->next;
  }

  return logPropRatio;
}

static void cross_product(REAL8 x[3], const REAL8 y[3], const REAL8 z[3]) {
    x[0] = y[1]*z[2]-y[2]*z[1];
    x[1] = y[2]*z[0]-y[0]*z[2];
    x[2] = y[0]*z[1]-y[1]*z[0];
}

static REAL8 norm(const REAL8 x[3]) {
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

static void unit_vector(REAL8 v[3], const REAL8 w[3]) {
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

static REAL8 dot(const REAL8 v[3], const REAL8 w[3]) {
    return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

static void project_along(REAL8 vproj[3], const REAL8 v[3], const REAL8 w[3]) {
    REAL8 what[3];
    REAL8 vdotw;

    unit_vector(what, w);
    vdotw = dot(v, w);

    vproj[0] = what[0]*vdotw;
    vproj[1] = what[1]*vdotw;
    vproj[2] = what[2]*vdotw;
}

static void vsub(REAL8 diff[3], const REAL8 w[3], const REAL8 v[3]) {
    diff[0] = w[0] - v[0];
    diff[1] = w[1] - v[1];
    diff[2] = w[2] - v[2];
}

static void reflect_plane(REAL8 pref[3], const REAL8 p[3],
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

static void sph_to_cart(REAL8 cart[3], const REAL8 lat, const REAL8 longi) {
    cart[0] = cos(longi)*cos(lat);
    cart[1] = sin(longi)*cos(lat);
    cart[2] = sin(lat);
}

static void cart_to_sph(const REAL8 cart[3], REAL8 *lat, REAL8 *longi) {
    *longi = atan2(cart[1], cart[0]);
    *lat = asin(cart[2] / sqrt(cart[0]*cart[0] + cart[1]*cart[1] + cart[2]*cart[2]));
}

static void reflected_position_and_time(LALInferenceThreadState *thread, const REAL8 ra, const REAL8 dec,
                                        const REAL8 oldTime, REAL8 *newRA, REAL8 *newDec, REAL8 *newTime) {
    LALStatus status;
    memset(&status, 0, sizeof(status));
    SkyPosition currentEqu, currentGeo, newEqu, newGeo;
    LALDetector *detectors;
    LIGOTimeGPS *epoch;
    REAL8 x[3], y[3], z[3];
    REAL8 currentLoc[3], newLoc[3];
    REAL8 newGeoLat, newGeoLongi;
    REAL8 oldDt, newDt;
    LALDetector xD, yD, zD;

    currentEqu.latitude = dec;
    currentEqu.longitude = ra;
    currentEqu.system = COORDINATESYSTEM_EQUATORIAL;
    currentGeo.system = COORDINATESYSTEM_GEOGRAPHIC;

    LALInferenceVariables *args = thread->proposalArgs;

    epoch = (LIGOTimeGPS *)LALInferenceGetVariable(args, "epoch");
    detectors = *(LALDetector **)LALInferenceGetVariable(args, "detectors");

    LALEquatorialToGeographic(&status, &currentGeo, &currentEqu, epoch);

    /* This function should only be called when we know that we have
     three detectors, or the following will crash. */
    xD = detectors[0];
    memcpy(x, xD.location, 3*sizeof(REAL8));

    INT4 det = 1;
    yD = detectors[det];
    while (same_detector_location(&yD, &xD)) {
        det++;
        yD = detectors[det];
    }
    memcpy(y, yD.location, 3*sizeof(REAL8));
    det++;

    zD = detectors[det];
    while (same_detector_location(&zD, &yD) || same_detector_location(&zD, &xD)) {
        det++;
        zD = detectors[det];
    }
    memcpy(z, zD.location, 3*sizeof(REAL8));

    sph_to_cart(currentLoc, currentGeo.latitude, currentGeo.longitude);

    reflect_plane(newLoc, currentLoc, x, y, z);

    cart_to_sph(newLoc, &newGeoLat, &newGeoLongi);

    newGeo.latitude = newGeoLat;
    newGeo.longitude = newGeoLongi;
    newGeo.system = COORDINATESYSTEM_GEOGRAPHIC;
    newEqu.system = COORDINATESYSTEM_EQUATORIAL;
    LALGeographicToEquatorial(&status, &newEqu, &newGeo, epoch);

    oldDt = XLALTimeDelayFromEarthCenter(detectors[0].location, currentEqu.longitude,
                                         currentEqu.latitude, epoch);
    newDt = XLALTimeDelayFromEarthCenter(detectors[0].location, newEqu.longitude,
                                         newEqu.latitude, epoch);

    *newRA = newEqu.longitude;
    *newDec = newEqu.latitude;
    *newTime = oldTime + oldDt - newDt;
}

static REAL8 evaluate_morlet_proposal(LALInferenceThreadState *thread,
                                      LALInferenceVariables *proposedParams,
                                      INT4 ifo, INT4 k) {
    REAL8 prior = 0.0;
    REAL8 component_min,component_max;
    REAL8 A, f, Q, Anorm;
    gsl_matrix *glitch_f, *glitch_Q, *glitch_A;

    component_min = LALInferenceGetREAL8Variable(thread->priorArgs,"morlet_f0_prior_min");
    component_max = LALInferenceGetREAL8Variable(thread->priorArgs,"morlet_f0_prior_max");
    prior -= log(component_max - component_min);

    component_min = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_Q_prior_min");
    component_max = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_Q_prior_max");
    prior -= log(component_max - component_min);

    component_min = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_t0_prior_min");
    component_max = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_t0_prior_max");
    prior -= log(component_max - component_min);

    component_min = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_phi_prior_min");
    component_max = LALInferenceGetREAL8Variable(thread->priorArgs, "morlet_phi_prior_max");
    prior -= log(component_max - component_min);

    //"Malmquist" prior on A
    glitch_f = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
    glitch_Q = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
    glitch_A = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");

    A = gsl_matrix_get(glitch_A, ifo, k);
    Q = gsl_matrix_get(glitch_Q, ifo, k);
    f = gsl_matrix_get(glitch_f, ifo, k);

    Anorm = LALInferenceGetREAL8Variable(thread->priorArgs, "glitch_norm");

    prior += logGlitchAmplitudeDensity(A*Anorm, Q, f);

    return prior;
}


static REAL8 glitchAmplitudeDraw(REAL8 Q, REAL8 f, gsl_rng *r) {
    REAL8 SNR;
    REAL8 PIterm = 0.5*LAL_2_SQRTPI*LAL_SQRT1_2;
    REAL8 SNRPEAK = 5.0;

    INT4 k=0;
    REAL8 den=0.0, alpha=1.0;
    REAL8 max= 1.0/(SNRPEAK*LAL_E);;

    // x/a^2 exp(-x/a) prior on SNR. Peaks at x = a. Good choice is a=5

    // rejection sample. Envelope function on runs out to ten times the peak
    // don't even bother putting in this minute correction to the normalization
    // (it is a 5e-4 correction).
    do {
        SNR = 20.0 * SNRPEAK * gsl_rng_uniform(r);

        den = SNR/(SNRPEAK*SNRPEAK) * exp(-SNR/SNRPEAK);

        den /= max;

        alpha = gsl_rng_uniform(r);

        k++;
    } while (alpha > den);

    return SNR/sqrt((PIterm*Q/f));
}

REAL8 LALInferenceSkyRingProposal(LALInferenceThreadState *thread,
                                  LALInferenceVariables *currentParams,
                                  LALInferenceVariables *proposedParams) {
    INT4 i, j, l;
    INT4 ifo, nifo, timeflag=0;
    REAL8 logPropRatio = 0.0;
    REAL8 ra, dec;
    REAL8 baryTime, gmst;
    REAL8 newRA, newDec, newTime, newPsi;
    REAL8 intpart, decpart;
    REAL8 omega, cosomega, sinomega, c1momega;
    REAL8 IFO1[3], IFO2[3];
    REAL8 IFOX[3], k[3];
    REAL8 normalize;
    REAL8 pForward, pReverse;
    REAL8 n[3];
    REAL8 kp[3];
    LIGOTimeGPS GPSlal, *epoch;
    LALDetector *detectors;
    gsl_matrix *IFO;

    LALInferenceCopyVariables(currentParams, proposedParams);

    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng *rng = thread->GSLrandom;

    epoch = (LIGOTimeGPS *)LALInferenceGetVariable(args, "epoch");
    detectors = *(LALDetector **)LALInferenceGetVariable(args, "detectors");

    ra = LALInferenceGetREAL8Variable(proposedParams, "rightascension");
    dec = LALInferenceGetREAL8Variable(proposedParams, "declination");

    if (LALInferenceCheckVariable(proposedParams, "time")){
        baryTime = LALInferenceGetREAL8Variable(proposedParams, "time");
        timeflag = 1;
    } else {
        baryTime = XLALGPSGetREAL8(epoch);
    }

    XLALGPSSetREAL8(&GPSlal, baryTime);
    gmst = XLALGreenwichMeanSiderealTime(&GPSlal);

    //remap gmst back to [0:2pi]
    gmst /= LAL_TWOPI;
    intpart = (INT4)gmst;
    decpart = gmst - (REAL8)intpart;
    gmst = decpart*LAL_TWOPI;
    gmst = gmst < 0. ? gmst + LAL_TWOPI : gmst;

    /*
    line-of-sight vector
    */
    k[0] = cos(gmst-ra) * cos(dec);
    k[1] =-sin(gmst-ra) * cos(dec);
    k[2] = sin(dec);

    /*
    Store location for each detector
    */
    nifo = LALInferenceGetINT4Variable(args, "nDet");

    IFO = gsl_matrix_alloc(nifo, 3);

    for(ifo=0; ifo<nifo; ifo++) {
        memcpy(IFOX, detectors[ifo].location, 3*sizeof(REAL8));
        for (i=0; i<3; i++)
            gsl_matrix_set(IFO, ifo, i, IFOX[i]);
    }

    /*
    Randomly select two detectors from the network
    -this assumes there are no co-located detectors
    */
    i = j = 0;
    while (i==j) {
        i=gsl_rng_uniform_int(rng, nifo);
        j=gsl_rng_uniform_int(rng, nifo);
    }

    for(l=0; l<3; l++) {
        IFO1[l] = gsl_matrix_get(IFO, i, l);
        IFO2[l] = gsl_matrix_get(IFO, j, l);
    }

    /*
    detector axis
    */
    normalize=0.0;
    for(i=0; i<3; i++) {
        n[i] = IFO1[i] - IFO2[i];
        normalize += n[i] * n[i];
    }
    normalize = 1./sqrt(normalize);
    for(i=0; i<3; i++)
        n[i] *= normalize;

    /*
    rotation angle
    */
    omega = LAL_TWOPI * gsl_rng_uniform(rng);
    cosomega = cos(omega);
    sinomega = sin(omega);
    c1momega = 1.0 - cosomega;

    /*
    rotate k' = Rk
    */
    kp[0] = (c1momega*n[0]*n[0] + cosomega) * k[0]
            + (c1momega*n[0]*n[1] - sinomega*n[2]) * k[1]
            + (c1momega*n[0]*n[2] + sinomega*n[1]) * k[2];
    kp[1] = (c1momega*n[0]*n[1] + sinomega*n[2]) * k[0]
            + (c1momega*n[1]*n[1] + cosomega) * k[1]
            + (c1momega*n[1]*n[2] - sinomega*n[0]) * k[2];
    kp[2] = (c1momega*n[0]*n[2] - sinomega*n[1]) * k[0]
            + (c1momega*n[1]*n[2] + sinomega*n[0]) * k[1]
            + (c1momega*n[2]*n[2] + cosomega) * k[2];

    /*
    convert k' back to ra' and dec'
    */
    newDec = asin(kp[2]);
    newRA = atan2(kp[1], kp[0]) + gmst;
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
    for(i=0; i<3; i++) {
        tx += -IFO1[i]*k[i] /LAL_C_SI;
        ty += -IFO1[i]*kp[i]/LAL_C_SI;
    }
    newTime = tx + baryTime - ty;

    XLALGPSSetREAL8(&GPSlal, newTime);

    /*
    draw new polarisation angle uniformally
    for now
    MARK: Need to be smarter about psi in sky-ring jump
    */
    newPsi = LAL_PI * gsl_rng_uniform(rng);

    LALInferenceSetVariable(proposedParams, "polarisation", &newPsi);
    LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
    LALInferenceSetVariable(proposedParams, "declination", &newDec);
    if (timeflag)
        LALInferenceSetVariable(proposedParams, "time", &newTime);

    pForward = cos(newDec);
    pReverse = cos(dec);

    gsl_matrix_free(IFO);

    logPropRatio = log(pReverse/pForward);

    return logPropRatio;
}

REAL8 LALInferenceSkyReflectDetPlane(LALInferenceThreadState *thread,
                                     LALInferenceVariables *currentParams,
                                     LALInferenceVariables *proposedParams) {
    INT4 timeflag=0;
    REAL8 ra, dec, baryTime;
    REAL8 newRA, newDec, newTime;
    REAL8 nRA, nDec, nTime;
    REAL8 refRA, refDec, refTime;
    REAL8 nRefRA, nRefDec, nRefTime;
    REAL8 pForward, pReverse;
    REAL8 logPropRatio = 0.0;
    INT4 nUniqueDet;
    LIGOTimeGPS *epoch;


    /* Find the number of distinct-position detectors. */
    /* Exit with same parameters (with a warning the first time) if
    there are not three detectors. */
    static INT4 warningDelivered = 0;

    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng *rng = thread->GSLrandom;

    epoch = (LIGOTimeGPS *)LALInferenceGetVariable(args, "epoch");
    nUniqueDet = LALInferenceGetINT4Variable(args, "nUniqueDet");

    if (nUniqueDet != 3) {
        if (!warningDelivered) {
            fprintf(stderr, "WARNING: trying to reflect through the decector plane with %d\n", nUniqueDet);
            fprintf(stderr, "WARNING: geometrically independent locations,\n");
            fprintf(stderr, "WARNING: but this proposal should only be used with exactly 3 independent detectors.\n");
            fprintf(stderr, "WARNING: %s, line %d\n", __FILE__, __LINE__);
            warningDelivered = 1;
        }

        return logPropRatio;
    }
    LALInferenceCopyVariables(currentParams, proposedParams);

    ra = LALInferenceGetREAL8Variable(currentParams, "rightascension");
    dec = LALInferenceGetREAL8Variable(currentParams, "declination");

    if (LALInferenceCheckVariable(currentParams, "time")){
        baryTime = LALInferenceGetREAL8Variable(currentParams, "time");
        timeflag=1;
    } else {
        baryTime = XLALGPSGetREAL8(epoch);
    }

    reflected_position_and_time(thread, ra, dec, baryTime, &newRA, &newDec, &newTime);

    /* Unit normal deviates, used to "fuzz" the state. */
    const REAL8 epsTime = 6e-6; /* 1e-1 / (16 kHz) */
    const REAL8 epsAngle = 3e-4; /* epsTime*c/R_Earth */

    nRA = gsl_ran_ugaussian(rng);
    nDec = gsl_ran_ugaussian(rng);
    nTime = gsl_ran_ugaussian(rng);

    newRA += epsAngle*nRA;
    newDec += epsAngle*nDec;
    newTime += epsTime*nTime;

    /* And the doubly-reflected position (near the original, but not
    exactly due to the fuzzing). */
    reflected_position_and_time(thread, newRA, newDec, newTime, &refRA, &refDec, &refTime);

    /* The Gaussian increments required to shift us back to the original
    position from the doubly-reflected position. */
    nRefRA = (ra - refRA)/epsAngle;
    nRefDec = (dec - refDec)/epsAngle;
    nRefTime = (baryTime - refTime)/epsTime;

    pForward = gsl_ran_ugaussian_pdf(nRA) * gsl_ran_ugaussian_pdf(nDec) * gsl_ran_ugaussian_pdf(nTime);
    pReverse = gsl_ran_ugaussian_pdf(nRefRA) * gsl_ran_ugaussian_pdf(nRefDec) * gsl_ran_ugaussian_pdf(nRefTime);

    LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
    LALInferenceSetVariable(proposedParams, "declination", &newDec);

    if (timeflag)
        LALInferenceSetVariable(proposedParams, "time", &newTime);

    logPropRatio = log(pReverse/pForward);

    return logPropRatio;
}

REAL8 LALInferencePSDFitJump(LALInferenceThreadState *thread,
                             LALInferenceVariables *currentParams,
                             LALInferenceVariables *proposedParams) {
    INT4 i,j;
    INT4 N, nifo;
    REAL8 draw=0.0;
    REAL8 logPropRatio = 0.0;
    REAL8Vector *var;
    gsl_matrix *ny;

    LALInferenceCopyVariables(currentParams, proposedParams);

    var = LALInferenceGetREAL8VectorVariable(thread->proposalArgs, "psdsigma");

    //Get current state of chain into workable form
    ny = LALInferenceGetgslMatrixVariable(proposedParams, "psdscale");

    //Get size of noise parameter array
    nifo = (INT4)ny->size1;
    N = (INT4)ny->size2;

    //perturb noise parameter
    for(i=0; i<nifo; i++) {
        for(j=0; j<N; j++) {
            draw = gsl_matrix_get(ny, i, j) + gsl_ran_ugaussian(thread->GSLrandom) * var->data[j];
            gsl_matrix_set(ny, i, j, draw);
        }
    }

    return logPropRatio;
}

static void UpdateWaveletSum(LALInferenceThreadState *thread,
                             LALInferenceVariables *proposedParams,
                             gsl_matrix *glitchFD, INT4 ifo, INT4 n, INT4 flag) {
    INT4 i=0;
    INT4 lower, upper;
    INT4 glitchLower, glitchUpper;
    REAL8FrequencySeries **asds, *asd = NULL;
    REAL8Vector *flows;
    REAL8 deltaT, Tobs, deltaF;
    REAL8 Q, Amp, t0, ph0, f0; //sine-Gaussian parameters
    REAL8 amparg, phiarg, Ai;//helpers for computing sineGaussian
    REAL8 gRe, gIm;         //real and imaginary parts of current glitch model
    REAL8 tau;
    gsl_matrix *glitch_f, *glitch_Q, *glitch_A;
    gsl_matrix *glitch_t, *glitch_p;
    REAL8TimeSeries **td_data;

    LALInferenceVariables *args = thread->proposalArgs;

    asds = *(REAL8FrequencySeries ***)LALInferenceGetVariable(args, "asds");
    flows = LALInferenceGetREAL8VectorVariable(args, "flows");
    td_data = *(REAL8TimeSeries ***)LALInferenceGetVariable(args, "td_data");

    /* get dataPtr pointing to correct IFO */
    asd = asds[ifo];

    deltaT = td_data[ifo]->deltaT;
    Tobs = (((REAL8)td_data[ifo]->data->length) * deltaT);
    deltaF = 1.0 / Tobs;

    lower = (INT4)ceil(flows->data[ifo] / deltaF);
    upper = (INT4)floor(flows->data[ifo] / deltaF);

    glitch_f = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
    glitch_Q = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
    glitch_A = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");
    glitch_t = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_t0");
    glitch_p = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_phi");

    Q = gsl_matrix_get(glitch_Q, ifo, n);
    Amp = gsl_matrix_get(glitch_A, ifo, n);
    t0 = gsl_matrix_get(glitch_t, ifo, n);
    ph0 = gsl_matrix_get(glitch_p, ifo, n);
    f0 = gsl_matrix_get(glitch_f, ifo, n);

    //6 x decay time of sine Gaussian (truncate how much glitch we compute)
    tau = Q/LAL_TWOPI/f0;
    glitchLower = (INT4)floor((f0 - 1./tau)/deltaF);
    glitchUpper = (INT4)floor((f0 + 1./tau)/deltaF);

    //set glitch model to zero
    if (flag==0) {
        for(i=lower; i<=upper; i++) {
            gsl_matrix_set(glitchFD, ifo, 2*i, 0.0);
            gsl_matrix_set(glitchFD, ifo, 2*i+1, 0.0);
        }
    }

    for (i=glitchLower; i<glitchUpper; i++) {
        if (i>=lower && i<=upper) {
            gRe = gsl_matrix_get(glitchFD, ifo, 2*i);
            gIm = gsl_matrix_get(glitchFD, ifo, 2*i+1);
            amparg = ((REAL8)i*deltaF - f0)*LAL_PI*tau;
            phiarg = LAL_PI*(REAL8)i + ph0 - LAL_TWOPI*(REAL8)i*deltaF*(t0-Tobs/2.);//TODO: SIMPLIFY PHASE FOR SINEGAUSSIAN
            Ai = Amp*tau*0.5*sqrt(LAL_PI)*exp(-amparg*amparg)*asd->data->data[i]/sqrt(Tobs);

            switch(flag) {
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
            gsl_matrix_set(glitchFD, ifo, 2*i, gRe);
            gsl_matrix_set(glitchFD, ifo, 2*i+1, gIm);

        }//end upper/lower check
    }//end loop over glitch samples
}

static void phase_blind_time_shift(REAL8 *corr, REAL8 *corrf, COMPLEX16Vector *data1,
                                   COMPLEX16Vector *data2, INT4 ifo, LALInferenceVariables *args) {
    INT4 i, N, N2;
    INT4 lower, upper;
    REAL8 deltaF, deltaT;
    REAL8FrequencySeries **psds;
    REAL8FrequencySeries *psd;
    COMPLEX16FrequencySeries *corrFD, *corrfFD;
    REAL8TimeSeries *corrTD, *corrfTD;
    REAL8TimeSeries **td_data;
    COMPLEX16FrequencySeries **fd_data;
    REAL8FFTPlan **plans;
    REAL8Vector *flows, *fhighs;

    psds = *(REAL8FrequencySeries ***)LALInferenceGetVariable(args, "psds");
    flows = LALInferenceGetREAL8VectorVariable(args, "flows");
    fhighs = LALInferenceGetREAL8VectorVariable(args, "fhighs");

    td_data = *(REAL8TimeSeries ***)LALInferenceGetVariable(args, "td_data");
    fd_data = *(COMPLEX16FrequencySeries ***)LALInferenceGetVariable(args, "fd_data");

    plans = *(REAL8FFTPlan ***)LALInferenceGetVariable(args, "f2t_plans");

    /* get dataPtr pointing to correct IFO */
    psd = psds[ifo];

    N  = td_data[ifo]->data->length;   // Number of data points
    N2 = fd_data[ifo]->data->length-1; // 1/2 number of data points (plus 1)

    deltaF = fd_data[ifo]->deltaF;
    deltaT = td_data[ifo]->deltaT;

    lower  = (INT4)ceil(flows->data[ifo]  / deltaF);
    upper  = (INT4)floor(fhighs->data[ifo] / deltaF);

    corrFD  = XLALCreateCOMPLEX16FrequencySeries("cf1", &(fd_data[ifo]->epoch), 0.0, deltaF, &lalDimensionlessUnit, N2+1);
    corrfFD = XLALCreateCOMPLEX16FrequencySeries("cf2", &(fd_data[ifo]->epoch), 0.0, deltaF, &lalDimensionlessUnit, N2+1);

    corrTD  = XLALCreateREAL8TimeSeries("ct1", &(td_data[ifo]->epoch), 0.0, deltaT, &lalDimensionlessUnit, N);
    corrfTD = XLALCreateREAL8TimeSeries("ct2", &(td_data[ifo]->epoch), 0.0, deltaT, &lalDimensionlessUnit, N);

    //convolution of signal & template
    for (i=0; i < N2; i++) {
        corrFD->data->data[i]  = crect(0.0,0.0);
        corrfFD->data->data[i] = crect(0.0,0.0);

        if(i>lower && i<upper) {
            corrFD->data->data[i] = crect( ( creal(data1->data[i])*creal(data2->data[i]) + cimag(data1->data[i])*cimag(data2->data[i])) / psd->data->data[i],
                                           ( cimag(data1->data[i])*creal(data2->data[i]) - creal(data1->data[i])*cimag(data2->data[i])) / psd->data->data[i] );
            corrfFD->data->data[i] = crect( ( creal(data1->data[i])*cimag(data2->data[i]) - cimag(data1->data[i])*creal(data2->data[i])) / psd->data->data[i],
                                            ( cimag(data1->data[i])*cimag(data2->data[i]) + creal(data1->data[i])*creal(data2->data[i])) / psd->data->data[i] );
        }
    }

    //invFFT convolutions to find time offset
    XLALREAL8FreqTimeFFT(corrTD, corrFD, plans[ifo]);
    XLALREAL8FreqTimeFFT(corrfTD, corrfFD, plans[ifo]);

    for (i=0; i < N; i++) {
        corr[i]  = corrTD->data->data[i];
        corrf[i] = corrfTD->data->data[i];
    }

    XLALDestroyREAL8TimeSeries(corrTD);
    XLALDestroyREAL8TimeSeries(corrfTD);
    XLALDestroyCOMPLEX16FrequencySeries(corrFD);
    XLALDestroyCOMPLEX16FrequencySeries(corrfFD);
}

static void MaximizeGlitchParameters(LALInferenceThreadState *thread,
                                     LALInferenceVariables *currentParams,
                                     INT4 ifo, INT4 n)
{
    INT4 i, imax, N;
    INT4 lower, upper;
    REAL8 deltaT, Tobs, deltaF, sqTwoDeltaToverN;
    REAL8 Amp, t0, ph0;
    REAL8 rho=0.0;
    REAL8 hRe, hIm;
    REAL8 gRe, gIm;
    REAL8 dPhase, dTime;
    REAL8 max;
    REAL8 *corr, *AC, *AF;
    REAL8FrequencySeries **psds;
    REAL8Vector *flows, *fhighs, *Sn;
    INT4Vector *gsize;
    COMPLEX16Sequence *s, *h, *r;
    gsl_matrix *glitchFD, *glitch_A, *glitch_t, *glitch_p, *hmatrix;
    REAL8TimeSeries **td_data;
    COMPLEX16FrequencySeries **fd_data;

    LALInferenceVariables *args = thread->proposalArgs;

    INT4 nDet = LALInferenceGetINT4Variable(args, "nDet");
    psds = *(REAL8FrequencySeries ***)LALInferenceGetVariable(args, "psds");
    flows = LALInferenceGetREAL8VectorVariable(args, "flows");
    fhighs = LALInferenceGetREAL8VectorVariable(args, "fhighs");

    td_data = XLALCalloc(nDet, sizeof(REAL8TimeSeries *));
    fd_data = XLALCalloc(nDet, sizeof(COMPLEX16FrequencySeries *));

    N = td_data[ifo]->data->length;
    deltaT = td_data[ifo]->deltaT;
    Tobs = (REAL8)(deltaT*N);
    sqTwoDeltaToverN = sqrt(2.0 * deltaT / ((REAL8) N) );

    deltaF = 1.0 / (((REAL8)N) * deltaT);
    lower = (INT4)ceil(flows->data[ifo] / deltaF);
    upper = (INT4)floor(fhighs->data[ifo] / deltaF);

    s = fd_data[ifo]->data;
    h = XLALCreateCOMPLEX16Vector(N/2);
    r = XLALCreateCOMPLEX16Vector(N/2);
    Sn = psds[ifo]->data;

    /* Get parameters for new wavelet */
    gsize = LALInferenceGetINT4VectorVariable(currentParams, "glitch_size");

    glitchFD = LALInferenceGetgslMatrixVariable(currentParams, "morlet_FD");
    glitch_A = LALInferenceGetgslMatrixVariable(currentParams, "morlet_Amp");
    glitch_t = LALInferenceGetgslMatrixVariable(currentParams, "morlet_t0");
    glitch_p = LALInferenceGetgslMatrixVariable(currentParams, "morlet_phi");

    /* sine-Gaussian parameters */
    Amp = gsl_matrix_get(glitch_A, ifo, n);
    t0 = gsl_matrix_get(glitch_t, ifo, n);
    ph0 = gsl_matrix_get(glitch_p, ifo, n);

    /* Make new wavelet */
    hmatrix = gsl_matrix_alloc(ifo+1, N);
    gsl_matrix_set_all(hmatrix, 0.0);

    UpdateWaveletSum(thread, currentParams, hmatrix, ifo, n, 1);

    /* Copy to appropriate template array*/
    for (i=0; i<N/2; i++) {
        hRe = 0.0;
        hIm = 0.0;
        gRe = 0.0;
        gIm = 0.0;
        r->data[i] = crect(0.0, 0.0);

        if(i>lower && i<upper) {
            hRe = sqTwoDeltaToverN * gsl_matrix_get(hmatrix, ifo, 2*i);
            hIm = sqTwoDeltaToverN * gsl_matrix_get(hmatrix, ifo, 2*i+1);
            h->data[i] = crect(hRe, hIm);
            //compute SNR of new wavelet
            rho += (hRe*hRe + hIm*hIm) / Sn->data[i];

            //form up residual while we're in here (w/out new template)
            if(gsize->data[ifo]>0) {
                gRe = gsl_matrix_get(glitchFD, ifo, 2*i);
                gIm = gsl_matrix_get(glitchFD, ifo, 2*i+1);
            }
            r->data[i] = crect(sqTwoDeltaToverN * (creal(s->data[i])/deltaT-gRe),
                               sqTwoDeltaToverN * (cimag(s->data[i])/deltaT-gIm));
        }
    }
    rho*=4.0;

    /* Compute correlation of data & template */
    corr = XLALMalloc(sizeof(REAL8) * N);
    AF = XLALMalloc(sizeof(REAL8) * N);
    AC = XLALMalloc(sizeof(REAL8) * N);

    for(i=0; i<N; i++)
        corr[i] = 0.0;

    /* Cross-correlate template & residual */
    phase_blind_time_shift(AC, AF, r, h, ifo, thread->proposalArgs);

    for(i=0; i<N; i++)
        corr[i] += sqrt(AC[i]*AC[i] + AF[i]*AF[i]);

    /* Find element where correlation is maximized */
    max = corr[0];
    imax = 0;
    for(i=1; i<N; i++) {
        if(corr[i] > max) {
            max  = corr[i];
            imax = i;
        }
    }
    max *= 4.0;

    /* Get phase shift at max correlation */
    dPhase = atan2(AF[imax], AC[imax]);

    /* Compute time shift needed for propsed template */
    if (imax < (N/2)-1)
        dTime = ((REAL8)imax/(REAL8)N) * Tobs;
    else
        dTime = (((REAL8)imax-(REAL8)N)/(REAL8)N) * Tobs;

    /* Shift template parameters accordingly */
    t0 += dTime;
    Amp *= 1.0;//dAmplitude;
    ph0 -= dPhase;

    /* Map time & phase back in range if necessary */
    if (ph0 < 0.0)
        ph0 += LAL_TWOPI;
    else if (ph0 > LAL_TWOPI)
        ph0 -= LAL_TWOPI;

    if (t0 < 0.0)
        t0 += Tobs;
    else if (t0 > Tobs)
        t0 -= Tobs;

    gsl_matrix_set(glitch_t, ifo, n, t0);
    gsl_matrix_set(glitch_A, ifo, n, Amp);
    gsl_matrix_set(glitch_p, ifo, n, ph0);

    gsl_matrix_free(hmatrix);

    XLALDestroyCOMPLEX16Vector(h);
    XLALDestroyCOMPLEX16Vector(r);

    XLALFree(corr);
    XLALFree(AF);
    XLALFree(AC);

}

static void MorletDiagonalFisherMatrix(REAL8Vector *params, REAL8Vector *sigmas) {
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

REAL8 LALInferenceGlitchMorletProposal(LALInferenceThreadState *thread,
                                       LALInferenceVariables *currentParams,
                                       LALInferenceVariables *proposedParams) {
    INT4 i, ifo;
    INT4 n;

    REAL8 logPropRatio = 0.0;
    REAL8 t0;
    REAL8 f0;
    REAL8 Q;
    REAL8 Amp;
    REAL8 phi0;

    REAL8 scale;

    REAL8 qyx;
    REAL8 qxy;
    REAL8 Anorm;

    LALInferenceCopyVariables(currentParams, proposedParams);

    gsl_matrix *glitchFD, *glitch_f, *glitch_Q, *glitch_A, *glitch_t, *glitch_p;

    gsl_rng *rng = thread->GSLrandom;
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

    /* Get glitch meta paramters (dimnsion, proposal) */
    INT4Vector *gsize = LALInferenceGetINT4VectorVariable(proposedParams, "glitch_size");

    glitchFD = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_FD");
    glitch_f = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
    glitch_Q = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
    glitch_A = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");
    glitch_t = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_t0");
    glitch_p = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_phi");

    Anorm = LALInferenceGetREAL8Variable(thread->priorArgs, "glitch_norm");

    /* Choose which IFO */
    ifo = (INT4)floor(gsl_rng_uniform(rng) * (REAL8)(gsize->length));

    /* Bail out of proposal if no wavelets */
    if (gsize->data[ifo]==0) {
        XLALDestroyREAL8Vector(params_x);
        XLALDestroyREAL8Vector(params_y);
        XLALDestroyREAL8Vector(sigmas_x);
        XLALDestroyREAL8Vector(sigmas_y);

        return logPropRatio;
    }

    /* Choose which glitch */
    n = (INT4)floor(gsl_rng_uniform(rng) * (REAL8)(gsize->data[ifo]));

    /* Remove wavlet form linear combination */
    UpdateWaveletSum(thread, proposedParams, glitchFD, ifo, n, -1);

    /* Get parameters of n'th glitch int params vector */
    t0 = gsl_matrix_get(glitch_t, ifo, n); //Centroid time
    f0 = gsl_matrix_get(glitch_f, ifo, n); //Frequency
    Q = gsl_matrix_get(glitch_Q, ifo, n); //Quality
    Amp = gsl_matrix_get(glitch_A, ifo, n); //Amplitude
    phi0 = gsl_matrix_get(glitch_p, ifo, n); //Centroid phase


    /* Map to params Vector and compute Fisher */
    params_x->data[0] = t0;
    params_x->data[1] = f0;
    params_x->data[2] = Q;
    params_x->data[3] = Amp * (0.25*Anorm);//TODO: What is the 0.25*Anorm about?
    params_x->data[4] = phi0;

    MorletDiagonalFisherMatrix(params_x, sigmas_x);

    /* Jump from x -> y:  y = x + N[0,sigmas_x]*scale */
    scale = 0.4082482; // 1/sqrt(6)

    for(i=0; i<5; i++)
        params_y->data[i] = params_x->data[i] + gsl_ran_ugaussian(rng)*sigmas_x->data[i]*scale;

    /* Set parameters of n'th glitch int params vector */
    /* Map to params Vector and compute Fisher */
    t0   = params_y->data[0];
    f0   = params_y->data[1];
    Q    = params_y->data[2];
    Amp  = params_y->data[3]/(0.25*Anorm);
    phi0 = params_y->data[4];

    gsl_matrix_set(glitch_t, ifo, n, t0);
    gsl_matrix_set(glitch_f, ifo, n, f0);
    gsl_matrix_set(glitch_Q, ifo, n, Q);
    gsl_matrix_set(glitch_A, ifo, n, Amp);
    gsl_matrix_set(glitch_p, ifo, n, phi0);

    /* Add wavlet to linear combination */
    UpdateWaveletSum(thread, proposedParams, glitchFD, ifo, n, 1);

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
    for(i=0; i<5; i++) {
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

    logPropRatio = qxy-qyx;

    XLALDestroyREAL8Vector(params_x);
    XLALDestroyREAL8Vector(params_y);

    XLALDestroyREAL8Vector(sigmas_x);
    XLALDestroyREAL8Vector(sigmas_y);

    return logPropRatio;
}

REAL8 LALInferenceGlitchMorletReverseJump(LALInferenceThreadState *thread,
                                          LALInferenceVariables *currentParams,
                                          LALInferenceVariables *proposedParams) {
    INT4 i,n;
    INT4 ifo;
    INT4 rj, nx, ny;
    REAL8 draw;
    REAL8 val, Anorm;
    REAL8 t=0, f=0, Q=0, A=0;
    REAL8 qx = 0.0; //log amp proposals
    REAL8 qy = 0.0;
    REAL8 qyx = 0.0; //log pixel proposals
    REAL8 qxy = 0.0;
    REAL8 pForward = 0.0; //combined p() & q() probabilities for ...
    REAL8 pReverse = 0.0; //...RJMCMC hastings ratio
    REAL8 logPropRatio = 0.0;
    INT4 adapting=1;
    gsl_matrix *params = NULL;

    gsl_rng *rng = thread->GSLrandom;
    LALInferenceVariables *propArgs = thread->proposalArgs;

    LALInferenceCopyVariables(currentParams, proposedParams);

    INT4Vector *gsize = LALInferenceGetINT4VectorVariable(proposedParams, "glitch_size");
    gsl_matrix *glitchFD = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_FD");

    INT4 nmin = (INT4)(LALInferenceGetREAL8Variable(thread->priorArgs,"glitch_dim_min"));
    INT4 nmax = (INT4)(LALInferenceGetREAL8Variable(thread->priorArgs,"glitch_dim_max"));

    if (LALInferenceCheckVariable(propArgs, "adapting"))
        adapting = LALInferenceGetINT4Variable(propArgs, "adapting");

    /* Choose which IFO */
    ifo = (INT4)floor(gsl_rng_uniform(rng) * (REAL8)(gsize->length) );
    nx = gsize->data[ifo];

    /* Choose birth or death move */
    draw = gsl_rng_uniform(rng);
    if (draw < 0.5)
        rj = 1;
    else
        rj = -1;

    /* find dimension of proposed model */
    ny = nx + rj;

    /* Check that new dimension is allowed */
    if(ny<nmin || ny>=nmax) {
        logPropRatio = -DBL_MAX;
        return logPropRatio;
    }

    switch(rj) {
        /* Birth */
        case 1:
            //Add new wavelet to glitch model
            t = draw_flat(thread, "morlet_t0_prior");
            f = draw_flat(thread, "morlet_f0_prior");

            //Centroid time
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_t0");
            gsl_matrix_set(params, ifo, nx, t);

            //Frequency
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
            gsl_matrix_set(params, ifo, nx, f);

            //Quality
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
            val = draw_flat(thread, "morlet_Q_prior");
            gsl_matrix_set(params, ifo, nx, val);
            Q = val;

            //Amplitude
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");
            val = glitchAmplitudeDraw(Q, f, rng);
            Anorm = LALInferenceGetREAL8Variable(thread->priorArgs, "glitch_norm");
            A = val/Anorm;

            gsl_matrix_set(params, ifo, nx, A);

            //Centroid phase
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_phi");
            val = draw_flat(thread, "morlet_phi_prior");
            gsl_matrix_set(params, ifo, nx, val);

            //Maximize phase, time, and amplitude using cross-correlation of data & wavelet
            if (adapting)
                MaximizeGlitchParameters(thread, proposedParams, ifo, nx);

            //Add wavlet to linear combination
            UpdateWaveletSum(thread, proposedParams, glitchFD, ifo, nx, 1);

            //Compute probability of drawing parameters
            qy = evaluate_morlet_proposal(thread, proposedParams, ifo, nx);// + log(gsl_matrix_get(power,ifo,k));

            //Compute reverse probability of dismissing k
            qxy = 0.0;//-log((double)runState->data->freqData->data->length);//-log( (REAL8)ny );

            if (adapting)
                qy += 10.0;

            break;

        /* Death */
        case -1:
            //Choose wavelet to remove from glitch model
            draw = gsl_rng_uniform(rng);
            n = (INT4)(floor(draw * (REAL8)nx));     //choose which hot pixel

            // Remove wavlet from linear combination
            UpdateWaveletSum(thread, proposedParams, glitchFD, ifo, n, -1);

            //Get t and f of removed wavelet
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
            f = gsl_matrix_get(params, ifo, n);
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_t0");
            t = gsl_matrix_get(params, ifo, n);
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
            Q = gsl_matrix_get(params, ifo, n);
            params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");
            A = gsl_matrix_get(params, ifo, n);

            //Shift morlet parameters to fill in array
            for(i=n; i<ny; i++) {
                params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_f0");
                gsl_matrix_set(params, ifo, i, gsl_matrix_get(params, ifo, i+1));
                params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Q");
                gsl_matrix_set(params, ifo, i, gsl_matrix_get(params, ifo, i+1));
                params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_Amp");
                gsl_matrix_set(params, ifo, i, gsl_matrix_get(params, ifo, i+1));
                params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_t0");
                gsl_matrix_set(params, ifo, i, gsl_matrix_get(params, ifo, i+1));
                params = LALInferenceGetgslMatrixVariable(proposedParams, "morlet_phi");
                gsl_matrix_set(params, ifo, i, gsl_matrix_get(params, ifo, i+1));
            }

            //Compute reverse probability of drawing parameters
            //find TF pixel

            qx = evaluate_morlet_proposal(thread, currentParams, ifo, n);// + log(gsl_matrix_get(power,ifo,k));

            //Compute forward probability of dismissing k
            qyx = 0.0;//-log((double)runState->data->freqData->data->length);//0.0;//-log( (REAL8)nx );

            if(adapting)
                qx += 10.0;

            break;

        default:
            break;
    }

    /* Update proposal structure for return to MCMC */

    //Update model meta-date
    gsize->data[ifo] = ny;

    //Re-package prior and proposal ratios into runState
    pForward = qxy + qx;
    pReverse = qyx + qy;

    logPropRatio = pForward-pReverse;

    return logPropRatio;
}

REAL8 LALInferencePolarizationPhaseJump(UNUSED LALInferenceThreadState *thread,
                                        LALInferenceVariables *currentParams,
                                        LALInferenceVariables *proposedParams) {
    REAL8 logPropRatio = 0.0;

    LALInferenceCopyVariables(currentParams, proposedParams);

    REAL8 psi = LALInferenceGetREAL8Variable(proposedParams, "polarisation");
    REAL8 phi = LALInferenceGetREAL8Variable(proposedParams, "phase");

    phi += M_PI;
    psi += M_PI/2;

    phi = fmod(phi, 2.0*M_PI);
    psi = fmod(psi, M_PI);

    LALInferenceSetVariable(proposedParams, "polarisation", &psi);
    LALInferenceSetVariable(proposedParams, "phase", &phi);

    return logPropRatio;
}

REAL8 LALInferenceCorrPolarizationPhaseJump(LALInferenceThreadState *thread,
                                            LALInferenceVariables *currentParams,
                                            LALInferenceVariables *proposedParams) {
    REAL8 alpha,beta;
    REAL8 draw;
    REAL8 psi, phi;
    REAL8 logPropRatio = 0.0;

    LALInferenceCopyVariables(currentParams, proposedParams);

    gsl_rng *rng = thread->GSLrandom;

    psi = LALInferenceGetREAL8Variable(proposedParams, "polarisation");
    phi = LALInferenceGetREAL8Variable(proposedParams, "phase");

    alpha = psi + phi;
    beta  = psi - phi;

    //alpha =>   0:3pi
    //beta  => -2pi:pi

    //big jump in either alpha (beta) or beta (alpha)
    draw = gsl_rng_uniform(rng);
    if (draw < 0.5)
        alpha = gsl_rng_uniform(rng)*3.0*LAL_PI;
    else
        beta = -LAL_TWOPI + gsl_rng_uniform(rng)*3.0*LAL_PI;

    //transform back to psi,phi space
    psi = (alpha + beta)*0.5;
    phi = (alpha - beta)*0.5;

    //map back in range
    LALInferenceCyclicReflectiveBound(proposedParams, thread->priorArgs);

    LALInferenceSetVariable(proposedParams, "polarisation", &psi);
    LALInferenceSetVariable(proposedParams, "phase", &phi);

    return logPropRatio;
}

REAL8 LALInferenceFrequencyBinJump(LALInferenceThreadState *thread,
                                   LALInferenceVariables *currentParams,
                                   LALInferenceVariables *proposedParams) {
    REAL8 f0, df;
    REAL8 plusminus;
    REAL8 logPropRatio = 0.0;

    LALInferenceCopyVariables(currentParams, proposedParams);

    f0 = LALInferenceGetREAL8Variable(proposedParams, "f0");
    df = LALInferenceGetREAL8Variable(proposedParams, "df");

    plusminus = gsl_rng_uniform(thread->GSLrandom);
    if ( plusminus < 0.5 )
        f0 -= df;
    else
        f0 += df;

    LALInferenceSetVariable(proposedParams, "f0", &f0);

    return logPropRatio;
}

//This proposal needs to be called with exactly 3 independent detector locations.
static void reflected_extrinsic_parameters(LALInferenceThreadState *thread, const REAL8 ra, const REAL8 dec,
                                           const REAL8 baryTime, const REAL8 dist, const REAL8 iota, const REAL8 psi,
                                           REAL8 *newRA, REAL8 *newDec, REAL8 *newTime,
                                           REAL8 *newDist, REAL8 *newIota, REAL8 *newPsi) {
    REAL8 R2[4];
    REAL8 dist2;
    REAL8 gmst, newGmst;
    REAL8 cosIota, cosIota2;
    REAL8 Fplus, Fcross, psi_temp;
    REAL8 x[4], y[4], x2[4], y2[4];
    REAL8 newFplus[4], newFplus2[4], newFcross[4], newFcross2[4];
    REAL8 a, a2, b, c12;
    REAL8 cosnewIota, cosnewIota2;
    LIGOTimeGPS GPSlal;
    INT4 nUniqueDet, det;
    LALDetector *detectors;

    detectors = (LALDetector *)LALInferenceGetVariable(thread->proposalArgs, "detectors");
    nUniqueDet = LALInferenceGetINT4Variable(thread->proposalArgs, "nUniqueDet");

    XLALGPSSetREAL8(&GPSlal, baryTime);
    gmst = XLALGreenwichMeanSiderealTime(&GPSlal);

    reflected_position_and_time(thread, ra, dec, baryTime, newRA, newDec, newTime);

    XLALGPSSetREAL8(&GPSlal, *newTime);
    newGmst = XLALGreenwichMeanSiderealTime(&GPSlal);

    dist2 = dist*dist;

    cosIota = cos(iota);
    cosIota2 = cosIota*cosIota;

    /* Loop over interferometers */
    INT4 i=1, j=0;
    for (det=0; det < nUniqueDet; det++) {
        psi_temp = 0.0;

        XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])detectors[det].response, *newRA, *newDec, psi_temp, newGmst);
        j=i-1;
        while (j>0){
            if (Fplus==x[j]){
                det++;
                XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])detectors[det].response, *newRA, *newDec, psi_temp, newGmst);
            }
            j--;
        }
        x[i]=Fplus;
        x2[i]=Fplus*Fplus;
        y[i]=Fcross;
        y2[i]=Fcross*Fcross;

        XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])detectors[det].response, ra, dec, psi, gmst);
        R2[i] = (((1.0+cosIota2)*(1.0+cosIota2))/(4.0*dist2))*Fplus*Fplus
                + ((cosIota2)/(dist2))*Fcross*Fcross;

        i++;
    }

    a = (R2[3]*x2[2]*y2[1] - R2[2]*x2[3]*y2[1] - R2[3]*x2[1]*y2[2] + R2[1]*x2[3]*y2[2] + R2[2]*x2[1]*y2[3] -
        R2[1]*x2[2]*y2[3]);
    a2 = a*a;
    b = (-(R2[3]*x[1]*x2[2]*y[1]) + R2[2]*x[1]*x2[3]*y[1] + R2[3]*x2[1]*x[2]*y[2] - R2[1]*x[2]*x2[3]*y[2] +
        R2[3]*x[2]*y2[1]*y[2] - R2[3]*x[1]*y[1]*y2[2] - R2[2]*x2[1]*x[3]*y[3] + R2[1]*x2[2]*x[3]*y[3] - R2[2]*x[3]*y2[1]*y[3] + R2[1]*x[3]*y2[2]*y[3] +
        R2[2]*x[1]*y[1]*y2[3] - R2[1]*x[2]*y[2]*y2[3]);

    (*newPsi) = (2.*atan((b - a*sqrt((a2 + b*b)/(a2)))/a))/4.;

    while ((*newPsi)<0){
        (*newPsi)=(*newPsi)+LAL_PI/4.0;
    }

    while ((*newPsi)>LAL_PI/4.0){
        (*newPsi)=(*newPsi)-LAL_PI/4.0;
    }

    for (i = 1; i < 4; i++){
        newFplus[i] = x[i]*cos(2.0*(*newPsi)) + y[i]*sin(2.0*(*newPsi));
        newFplus2[i] = newFplus[i] * newFplus[i];

        newFcross[i] = y[i]*cos(2.0*(*newPsi)) - x[i]*sin(2.0*(*newPsi));
        newFcross2[i] = newFcross[i] * newFcross[i];
    }

    c12 = -2.0*((R2[1]*(newFcross2[2])-R2[2]*(newFcross2[1]))
          /(R2[1]*(newFplus2[2])-R2[2]*(newFplus2[1])))-1.0;

    if (c12<1.0){
        c12 = (3.0-c12)/(1.0+c12);
        (*newPsi) = (*newPsi)+LAL_PI/4.0;

        for (i = 1; i < 4; i++){
            newFplus[i] = x[i]*cos(2.0*(*newPsi)) + y[i]*sin(2.0*(*newPsi));
            newFplus2[i] = newFplus[i] * newFplus[i];

            newFcross[i] = y[i]*cos(2.0*(*newPsi)) - x[i]*sin(2.0*(*newPsi));
            newFcross2[i] = newFcross[i] * newFcross[i];
        }
    }

    if (c12<1){
        *newIota = iota;
        *newDist = dist;
        return;
    }

    cosnewIota2 = c12-sqrt(c12*c12-1.0);
    cosnewIota = sqrt(cosnewIota2);
    *newIota = acos(cosnewIota);

    *newDist = sqrt((
                    ((((1.0+cosnewIota2)*(1.0+cosnewIota2))/(4.0))*newFplus2[1]
                    + (cosnewIota2)*newFcross2[1])
                    )/ R2[1]);

    if (Fplus*newFplus[3]<0){
        (*newPsi)=(*newPsi)+LAL_PI/2.;
        newFcross[3]=-newFcross[3];
    }

    if (Fcross*cosIota*cosnewIota*newFcross[3]<0){
        (*newIota)=LAL_PI-(*newIota);
    }
}

REAL8 LALInferenceDistanceLikelihoodProposal(LALInferenceThreadState *thread,
                                             LALInferenceVariables *currentParams,
                                             LALInferenceVariables *proposedParams) {
  REAL8 old_d=0.0;
  DistanceParam distParam;

  if (LALInferenceCheckVariable(currentParams, "distance")) {
    distParam = USES_DISTANCE_VARIABLE;
    old_d = LALInferenceGetREAL8Variable(currentParams,"distance");
  } else if (LALInferenceCheckVariable(currentParams, "logdistance")) {
    distParam = USES_LOG_DISTANCE_VARIABLE;
    old_d = exp(LALInferenceGetREAL8Variable(currentParams,"logdistance"));
  } else {
    XLAL_ERROR_REAL8(XLAL_FAILURE, "could not find 'distance' or 'logdistance' in current params");
  }
  if(!LALInferenceCheckVariable(currentParams,"optimal_snr") || !LALInferenceCheckVariable(currentParams,"matched_filter_snr"))
    /* Force a likelihood calculation to generate the parameters */
    thread->parent->likelihood(currentParams,thread->parent->data,thread->model);
  if(!LALInferenceCheckVariable(currentParams,"optimal_snr") || !LALInferenceCheckVariable(currentParams,"matched_filter_snr"))
  {
		  /* If still can't find the required parameters, reject this jump */
		  LALInferenceCopyVariables(currentParams,proposedParams);
		  return(-INFINITY);
  }
  REAL8 OptimalSNR = LALInferenceGetREAL8Variable(currentParams,"optimal_snr");
  REAL8 MatchedFilterSNR = LALInferenceGetREAL8Variable(currentParams,"matched_filter_snr");
  REAL8 d_inner_h = MatchedFilterSNR * OptimalSNR;

  /* Get params of sampling distribution */
  REAL8 xmean = d_inner_h/(OptimalSNR*OptimalSNR*old_d);
  REAL8 xsigma = 1.0/(OptimalSNR*old_d);
  REAL8 old_x = 1.0/old_d;

  /* Sample new x. Do not check for x<0 since that can throw the code
   * into a difficult-to-terminate loop */
  REAL8 new_x;
  new_x = xmean + gsl_ran_gaussian(thread->GSLrandom,xsigma);
  REAL8 new_d = 1.0/new_x;

  LALInferenceCopyVariables(currentParams,proposedParams);
  /* Adjust SNRs */
  OptimalSNR *= new_x / old_x;
  LALInferenceSetREAL8Variable(proposedParams,"optimal_snr",OptimalSNR);

  /* Update individual detector information */
  const char *ifonames[5]={"H1","L1","V1","I1","J1"};
  char pname[64]="";
  for(UINT4 i=0;i<5;i++)
  {
    sprintf(pname,"%s_optimal_snr",ifonames[i]);
    if(LALInferenceCheckVariable(currentParams,pname))
      LALInferenceSetREAL8Variable(proposedParams,pname,LALInferenceGetREAL8Variable(currentParams,pname) * (new_x/old_x));
  }

  REAL8 logxdjac;
  /* Set distance */
  if(distParam == USES_DISTANCE_VARIABLE)
  {
    /* The natural proposal density is in x, but we need to compute
       the proposal ratio density in d.  So, we need a factor of

       |dx_old/dd_old| / |dx_new/dd_new| = d_new^2 / d_old^2

       to correct the x proposal density.
    */
    LALInferenceSetVariable(proposedParams,"distance",&new_d);
    logxdjac = 2.0*log(new_d) - 2.0*log(old_d);
  }
  else
  {
    /* The natural proposal density is in x, but we need to compute
       the proposal ratio in log(d).  So, we need a factor of

       |dx_old/dlog(d_old)| / |dx_new/dlog(d_new)| = d_new / d_old

       to correct the x proposal density.
     */
    REAL8 new_logd = log(new_d);
    LALInferenceSetVariable(proposedParams,"logdistance",&new_logd);
    logxdjac = log(new_d) - log(old_d);
  }
  /* Proposal ratio is not symmetric */
  /* P.R. = p(xold|xmean,xsigma) / p(xnew|xmean,xsigma) * jacobian */
  REAL8 log_p_reverse = -0.5*(old_x-xmean)*(old_x-xmean)/(xsigma*xsigma);
  REAL8 log_p_forward = -0.5*(new_x-xmean)*(new_x-xmean)/(xsigma*xsigma);

  return(log_p_reverse - log_p_forward + logxdjac);
}


REAL8 LALInferenceExtrinsicParamProposal(LALInferenceThreadState *thread,
                                         LALInferenceVariables *currentParams,
                                         LALInferenceVariables *proposedParams) {
    INT4 nUniqueDet;
    INT4 timeflag=0;
    REAL8 baryTime;
    REAL8 ra, dec;
    REAL8 psi, dist;
    REAL8 newRA, newDec, newTime, newDist, newIota, newPsi;
    REAL8 nRA, nDec, nTime, nDist, nIota, nPsi;
    REAL8 refRA, refDec, refTime, refDist, refIota, refPsi;
    REAL8 nRefRA, nRefDec, nRefTime, nRefDist, nRefIota, nRefPsi;
    REAL8 pForward, pReverse;
    REAL8 cst;
    REAL8 iota=0.0;
    REAL8 logPropRatio = 0.0;
    /* Find the number of distinct-position detectors. */
    /* Exit with same parameters (with a warning the first time) if
       there are not EXACTLY three unique detector locations. */
    static INT4 warningDelivered = 0;
    LIGOTimeGPS *epoch;


    LALInferenceVariables *args = thread->proposalArgs;
    gsl_rng *rng = thread->GSLrandom;
    epoch = (LIGOTimeGPS *)LALInferenceGetVariable(args, "epoch");

    nUniqueDet = LALInferenceGetINT4Variable(args, "nUniqueDet");
    if (nUniqueDet != 3) {
        if (!warningDelivered) {
            fprintf(stderr, "WARNING: trying to reflect through the decector plane with %d\n", nUniqueDet);
            fprintf(stderr, "WARNING: geometrically independent locations,\n");
            fprintf(stderr, "WARNING: but this proposal should only be used with exactly 3 independent detectors.\n");
            fprintf(stderr, "WARNING: %s, line %d\n", __FILE__, __LINE__);
            warningDelivered = 1;
        }

        return logPropRatio;
    }
    LALInferenceCopyVariables(currentParams, proposedParams);


    ra = LALInferenceGetREAL8Variable(proposedParams, "rightascension");
    dec = LALInferenceGetREAL8Variable(proposedParams, "declination");

    if (LALInferenceCheckVariable(proposedParams,"time")){
        baryTime = LALInferenceGetREAL8Variable(proposedParams, "time");
        timeflag = 1;
    } else {
        baryTime = XLALGPSGetREAL8(epoch);
    }

    if (LALInferenceCheckVariable(proposedParams,"costheta_jn"))
        iota = acos(LALInferenceGetREAL8Variable(proposedParams, "costheta_jn"));
    else
        fprintf(stderr, "LALInferenceExtrinsicParamProposal: No  theta_jn parameter!\n");

    psi = LALInferenceGetREAL8Variable(proposedParams, "polarisation");

    dist = exp(LALInferenceGetREAL8Variable(proposedParams, "logdistance"));

    reflected_extrinsic_parameters(thread, ra, dec, baryTime, dist, iota, psi, &newRA, &newDec, &newTime, &newDist, &newIota, &newPsi);

    /* Unit normal deviates, used to "fuzz" the state. */
    const REAL8 epsDist = 1e-8;
    const REAL8 epsTime = 1e-8;
    const REAL8 epsAngle = 1e-8;

    nRA = gsl_ran_ugaussian(rng);
    nDec = gsl_ran_ugaussian(rng);
    nTime = gsl_ran_ugaussian(rng);
    nDist = gsl_ran_ugaussian(rng);
    nIota = gsl_ran_ugaussian(rng);
    nPsi = gsl_ran_ugaussian(rng);

    newRA += epsAngle*nRA;
    newDec += epsAngle*nDec;
    newTime += epsTime*nTime;
    newDist += epsDist*nDist;
    newIota += epsAngle*nIota;
    newPsi += epsAngle*nPsi;

    /* And the doubly-reflected position (near the original, but not
    exactly due to the fuzzing). */
    reflected_extrinsic_parameters(thread, newRA, newDec, newTime, newDist, newIota, newPsi, &refRA, &refDec, &refTime, &refDist, &refIota, &refPsi);

    /* The Gaussian increments required to shift us back to the original
    position from the doubly-reflected position. */
    nRefRA = (ra - refRA)/epsAngle;
    nRefDec = (dec - refDec)/epsAngle;
    nRefTime = (baryTime - refTime)/epsTime;
    nRefDist = (dist - refDist)/epsDist;
    nRefIota = (iota - refIota)/epsAngle;
    nRefPsi = (psi - refPsi)/epsAngle;

    cst = log(1./(sqrt(2.*LAL_PI)));
    pReverse = 6*cst-0.5*(nRefRA*nRefRA+nRefDec*nRefDec+nRefTime*nRefTime+nRefDist*nRefDist+nRefIota*nRefIota+nRefPsi*nRefPsi);
    pForward = 6*cst-0.5*(nRA*nRA+nDec*nDec+nTime*nTime+nDist*nDist+nIota*nIota+nPsi*nPsi);

    LALInferenceSetVariable(proposedParams, "rightascension", &newRA);
    LALInferenceSetVariable(proposedParams, "declination", &newDec);
    if (timeflag)
        LALInferenceSetVariable(proposedParams, "time", &newTime);

    REAL8 logNewDist = log(newDist);
    LALInferenceSetVariable(proposedParams, "logdistance", &logNewDist);

    REAL8 newcosIota = cos(newIota);
    LALInferenceSetVariable(proposedParams, "costheta_jn", &newcosIota);
    LALInferenceSetVariable(proposedParams, "polarisation", &newPsi);

    logPropRatio = pReverse - pForward;

    return logPropRatio;
}


void LALInferenceSetupGlitchProposal(LALInferenceIFOData *data, LALInferenceVariables *propArgs) {
    INT4 i, nDet;
    REAL8Vector *flows, *fhighs;
    REAL8FrequencySeries **asds, **psds;
    REAL8TimeSeries **td_data;
    COMPLEX16FrequencySeries **fd_data;
    REAL8FFTPlan **plans;

    nDet = LALInferenceGetINT4Variable(propArgs, "nDet");

    flows = XLALCreateREAL8Vector(nDet);
    fhighs = XLALCreateREAL8Vector(nDet);
    asds = XLALCalloc(nDet, sizeof(REAL8FrequencySeries *));
    psds = XLALCalloc(nDet, sizeof(REAL8FrequencySeries *));
    td_data = XLALCalloc(nDet, sizeof(REAL8TimeSeries *));
    fd_data = XLALCalloc(nDet, sizeof(COMPLEX16FrequencySeries *));
    plans = XLALCalloc(nDet, sizeof(REAL8FFTPlan *));

    for (i=0; i<nDet; i++) {
        flows->data[i] = data->fLow;
        fhighs->data[i] = data->fHigh;

        asds[i] = data->noiseASD;
        psds[i] = data->oneSidedNoisePowerSpectrum;

        td_data[i] = data->timeData;
        fd_data[i] = data->freqData;

        plans[i] = data->freqToTimeFFTPlan;
        data = data->next;
    }

    LALInferenceAddREAL8VectorVariable(propArgs, "flows", flows, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddREAL8VectorVariable(propArgs, "fhighs", fhighs, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(propArgs, "asds", asds, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(propArgs, "psds", psds, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(propArgs, "td_data", td_data, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(propArgs, "fd_data", fd_data, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(propArgs, "f2t_plans", plans, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED);

}

/* Initialize differential evolution proposal */
void LALInferenceSetupDifferentialEvolutionProposal(LALInferenceThreadState *thread) {
    thread->differentialPoints = XLALCalloc(1, sizeof(LALInferenceVariables *));
    thread->differentialPointsLength = 0;
    thread->differentialPointsSize = 1;
    thread->differentialPointsSkip = 1;
}


/** Setup adaptive proposals. Should be called when state->currentParams is already filled with an initial sample */
void LALInferenceSetupAdaptiveProposals(LALInferenceVariables *propArgs, LALInferenceVariables *params) {
    INT4 no_adapt, adapting;
    INT4 adaptTau, adaptableStep, adaptLength, adaptResetBuffer;
    REAL8 sigma, s_gamma;
    REAL8 logLAtAdaptStart = -DBL_MAX;

    LALInferenceVariableItem *this;

    for(this=params->head; this; this=this->next) {
        if (LALInferenceCheckVariableNonFixed(params, this->name) && this->type == LALINFERENCE_REAL8_t) {
            char *name = this->name;

            if (!strcmp(name, "eta") || !strcmp(name, "q") || !strcmp(name, "time") || !strcmp(name, "a_spin2") || !strcmp(name, "a_spin1") || !strcmp(name,"t0")){
                sigma = 0.001;
            } else if (!strcmp(name, "polarisation") || !strcmp(name, "phase") || !strcmp(name, "costheta_jn")){
                sigma = 0.1;
            } else {
                sigma = 0.01;
            }

            /* Set up variables to store current sigma, proposed and accepted */
            char varname[MAX_STRLEN] = "";
            sprintf(varname, "%s_%s", name, ADAPTSUFFIX);
            LALInferenceAddREAL8Variable(propArgs, varname, sigma, LALINFERENCE_PARAM_LINEAR);

            sigma = 0.0;
            sprintf(varname, "%s_%s", name, ACCEPTSUFFIX);
            LALInferenceAddREAL8Variable(propArgs, varname, sigma, LALINFERENCE_PARAM_LINEAR);

            sprintf(varname, "%s_%s", name, PROPOSEDSUFFIX);
            LALInferenceAddREAL8Variable(propArgs, varname, sigma, LALINFERENCE_PARAM_LINEAR);
        }
    }

    no_adapt = LALInferenceGetINT4Variable(propArgs, "no_adapt");
    adapting = !no_adapt;      // Indicates if current iteration is being adapted
    LALInferenceAddINT4Variable(propArgs, "adapting", adapting, LALINFERENCE_PARAM_LINEAR);

    char name[MAX_STRLEN] = "none";
    LALInferenceAddstringVariable(propArgs, "proposedVariableName", name, LALINFERENCE_PARAM_OUTPUT);

    adaptableStep = 0;
    adaptTau = LALInferenceGetINT4Variable(propArgs, "adaptTau");  // Sets decay of adaption function
    adaptLength = pow(10, adaptTau);  // Number of iterations to adapt before turning off
    adaptResetBuffer = 100; // Number of iterations before adapting after a restart
    s_gamma = 1.0; // Sets the size of changes to jump size during adaptation

    LALInferenceAddINT4Variable(propArgs, "adaptableStep", adaptableStep, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddINT4Variable(propArgs, "adaptLength", adaptLength, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddINT4Variable(propArgs, "adaptResetBuffer", adaptResetBuffer, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddREAL8Variable(propArgs, "s_gamma", s_gamma, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddREAL8Variable(propArgs, "logLAtAdaptStart", logLAtAdaptStart, LALINFERENCE_PARAM_LINEAR);

    return;
}


/** Update proposal statistics if tracking */
void LALInferenceTrackProposalAcceptance(LALInferenceThreadState *thread) {
    INT4 i = 0;

    LALInferenceProposal *prop = thread->cycle->proposals[i];

    /* Find the proposal that was last called (by name) */
    while (strcmp(prop->name, thread->cycle->last_proposal_name)) {
        i++;
        prop = thread->cycle->proposals[i];
    }

    /* Update proposal statistics */
    prop->proposed++;
    if (thread->accepted == 1){
        prop->accepted++;
    }

    return;
}

/* Zero out proposal statistics counters */
void LALInferenceZeroProposalStats(LALInferenceProposalCycle *cycle) {
    INT4 i=0;

    for (i=0; i<cycle->nProposals; i++) {
        LALInferenceProposal *prop = cycle->proposals[i];

        prop->proposed = 0;
        prop->accepted = 0;
    }

    return;
}

/** Update the adaptive proposal. Whether or not a jump was accepted is passed with accepted */
void LALInferenceUpdateAdaptiveJumps(LALInferenceThreadState *thread, REAL8 targetAcceptance){
    INT4 adaptableStep = 0;
    INT4 adapting = 0;
    REAL8 priorMin, priorMax, dprior, s_gamma;
    REAL8 accept, propose, sigma;
    char *name;

    LALInferenceVariables *args = thread->proposalArgs;

    if (LALInferenceCheckVariable(args, "adaptableStep" ) &&
        LALInferenceCheckVariable(args, "adapting" )){
        adaptableStep = LALInferenceGetINT4Variable(args, "adaptableStep");
        adapting = LALInferenceGetINT4Variable(args, "adapting");
    }
    /* Don't do anything if these are not found */
    else return;

    if (adaptableStep && adapting) {
        char tmpname[MAX_STRLEN] = "";
        name = LALInferenceGetstringVariable(args, "proposedVariableName");

        sprintf(tmpname, "%s_%s", name, PROPOSEDSUFFIX);
        propose = LALInferenceGetREAL8Variable(args, tmpname) + 1;
        LALInferenceSetVariable(args, tmpname, &propose);

        sprintf(tmpname, "%s_%s", name,  ACCEPTSUFFIX);
        accept = LALInferenceGetREAL8Variable(args, tmpname) + thread->accepted;
        LALInferenceSetVariable(args, tmpname, &accept);
    }

    /* Adapt if desired. */
    if (LALInferenceCheckVariable(args, "proposedVariableName") &&
        LALInferenceCheckVariable(args, "s_gamma") &&
        LALInferenceCheckVariable(args, "adapting") &&
        LALInferenceCheckVariable(args, "adaptableStep")) {

        if (adaptableStep) {
            name = LALInferenceGetstringVariable(args, "proposedVariableName");
            char tmpname[MAX_STRLEN]="";

            sprintf(tmpname, "%s_%s", name, ADAPTSUFFIX);

            s_gamma = LALInferenceGetREAL8Variable(args, "s_gamma");

            sigma = LALInferenceGetREAL8Variable(args, tmpname);

            LALInferenceGetMinMaxPrior(thread->priorArgs, name, &priorMin, &priorMax);
            dprior = priorMax - priorMin;

            if (thread->accepted == 1){
                sigma += s_gamma * (dprior/100.0) * (1.0-targetAcceptance);
            } else {
                sigma -= s_gamma * (dprior/100.0) * (targetAcceptance);
            }

            sigma = (sigma > dprior ? dprior : sigma);
            sigma = (sigma < DBL_MIN ? DBL_MIN : sigma);
            LALInferenceSetVariable(args, tmpname, &sigma);
        }
    }

    /* Make sure we don't do this again until we take another adaptable step.*/
    adaptableStep = 0;
    LALInferenceSetVariable(args, "adaptableStep", &adaptableStep);
}




/**
 * Setup all clustered-KDE proposals with samples read from file.
 *
 * Constructed clustered-KDE proposals from all sample lists provided in
 * files given on the command line.
 * @param thread The LALInferenceThreadState to get command line options from and to the proposal cycle of.
 * @param input The input file pointer
 * @param burnin The number of burn-in points
 * @param weight The weight?
 * @param ptmcmc Flag to determine if using parallel tempering MCMC
 */
void LALInferenceSetupClusteredKDEProposalsFromASCII(LALInferenceThreadState *thread, FILE *input, INT4 burnin, REAL8 weight, INT4 ptmcmc) {
    LALInferenceVariableItem *item;
    INT4 j=0, k=0;

    INT4 cyclic_reflective = LALInferenceGetINT4Variable(thread->proposalArgs, "cyclic_reflective_kde");

    LALInferenceClusteredKDE *kde = XLALCalloc(1, sizeof(LALInferenceClusteredKDE));

    INT4 nInSamps;
    INT4 nCols;
    REAL8 *sampleArray;

    if (ptmcmc)
        LALInferenceDiscardPTMCMCHeader(input);

    char params[128][VARNAME_MAX];
    LALInferenceReadAsciiHeader(input, params, &nCols);

    LALInferenceVariables *backwardClusterParams = XLALCalloc(1, sizeof(LALInferenceVariables));

    /* Only cluster parameters that are being sampled */
    INT4 nValidCols=0;
    INT4 *validCols = XLALCalloc(nCols, sizeof(INT4));
    for (j=0; j<nCols; j++)
        validCols[j] = 0;

    INT4 logl_idx = 0;
    for (j=0; j<nCols; j++) {
        if (!strcmp("logl", params[j])) {
            logl_idx = j;
            continue;
        }

        char* internal_param_name = XLALCalloc(512, sizeof(char));
        LALInferenceTranslateExternalToInternalParamName(internal_param_name, params[j]);

        for (item = thread->currentParams->head; item; item = item->next) {
            if (!strcmp(item->name, internal_param_name) &&
                LALInferenceCheckVariableNonFixed(thread->currentParams, item->name)) {
                nValidCols++;
                validCols[j] = 1;
                LALInferenceAddVariable(backwardClusterParams, item->name, item->value, item->type, item->vary);
                break;
            }
        }
    }

    /* LALInferenceAddVariable() builds the array backwards, so reverse it. */
    LALInferenceVariables *clusterParams = XLALCalloc(1, sizeof(LALInferenceVariables));

    for (item = backwardClusterParams->head; item; item = item->next)
        LALInferenceAddVariable(clusterParams, item->name, item->value, item->type, item->vary);

    /* Burn in samples and parse the remainder */
    if (ptmcmc)
        LALInferenceBurninPTMCMC(input, logl_idx, nValidCols);
    else
        LALInferenceBurninStream(input, burnin);

    sampleArray = LALInferenceParseDelimitedAscii(input, nCols, validCols, &nInSamps);

    /* Downsample PTMCMC file to have independent samples */
    if (ptmcmc) {
        REAL8 acl_real8 = LALInferenceComputeMaxAutoCorrLen(sampleArray, nInSamps, nValidCols);
        INT4 acl;
        if (acl_real8 == INFINITY)
            acl = INT_MAX;
        else if (acl_real8 < 1.0)
            acl = 1;
        else
            acl = (INT4)acl_real8;

        INT4 downsampled_size = ceil((REAL8)nInSamps/acl);
        REAL8 *downsampled_array = (REAL8 *)XLALCalloc(downsampled_size * nValidCols, sizeof(REAL8));
        printf("Downsampling to achieve %i samples.\n", downsampled_size);
        for (k=0; k < downsampled_size; k++) {
            for (j=0; j < nValidCols; j++)
                downsampled_array[k*nValidCols + j] = sampleArray[k*nValidCols*acl + j];
        }
        XLALFree(sampleArray);
        sampleArray = downsampled_array;
        nInSamps = downsampled_size;
    }

    /* Build the KDE estimate and add to the KDE proposal set */
    INT4 ntrials = 50;  // Number of trials at fixed-k to find optimal BIC
    LALInferenceInitClusteredKDEProposal(thread, kde, sampleArray, nInSamps, clusterParams, clusteredKDEProposalName, weight, LALInferenceOptimizedKmeans, cyclic_reflective, ntrials);

    /* If kmeans construction failed, halt the run */
    if (!kde->kmeans) {
        fprintf(stderr, "\nERROR: Couldn't build kmeans clustering from the file specified.\n");
        XLALFree(kde);
        exit(-1);
    }

    LALInferenceAddClusteredKDEProposalToSet(thread->proposalArgs, kde);

    LALInferenceClearVariables(backwardClusterParams);
    XLALFree(backwardClusterParams);
    XLALFree(sampleArray);
}


/**
 * Initialize a clustered-KDE proposal.
 *
 * Estimates the underlying distribution of a set of points with a clustered kernel density estimate
 * and constructs a jump proposal from the estimated distribution.
 * @param      thread   The current LALInferenceThreadState.
 * @param[out] kde      An empty proposal structure to populate with the clustered-KDE estimate.
 * @param[in]  array    The data to estimate the underlying distribution of.
 * @param[in]  nSamps   Number of samples contained in \a array.
 * @param[in]  params   The parameters contained in \a array.
 * @param[in]  name     The name of the proposal being constructed.
 * @param[in]  weight   The relative weight this proposal is to have against other KDE proposals.
 * @param[in]  cluster_method A pointer to the clustering function to be used.
 * @param[in]  cyclic_reflective Flag to check for cyclic/reflective bounds.
 * @param[in]  ntrials  Number of kmeans attempts at fixed k to find optimal BIC.
 */
void LALInferenceInitClusteredKDEProposal(LALInferenceThreadState *thread,
                                          LALInferenceClusteredKDE *kde,
                                          REAL8 *array,
                                          INT4 nSamps,
                                          LALInferenceVariables *params,
                                          const char *name,
                                          REAL8 weight,
                                          LALInferenceKmeans* (*cluster_method)(gsl_matrix*, INT4, gsl_rng*),
                                          INT4 cyclic_reflective,
                                          INT4 ntrials) {
    INT4 dim;
    INT4 ndraws = 1000;
    gsl_matrix_view mview;
    char outp_name[256];
    char outp_draws_name[256];

    strcpy(kde->name, name);
    dim = LALInferenceGetVariableDimensionNonFixed(params);

    /* If kmeans is already assigned, assume it was calculated elsewhere */
    if (!kde->kmeans) {
        mview = gsl_matrix_view_array(array, nSamps, dim);
        kde->kmeans = (*cluster_method)(&mview.matrix, ntrials, thread->GSLrandom);
    }

    /* Return if kmeans setup failed */
    if (!kde->kmeans)
        return;

    kde->dimension = kde->kmeans->dim;
    kde->params = params;

    kde->weight = weight;
    kde->next = NULL;

    /* Selectivey impose bounds on KDEs */
    LALInferenceKmeansImposeBounds(kde->kmeans, params, thread->priorArgs, cyclic_reflective);

    /* Print out clustered samples, assignments, and PDF values if requested */
    if (LALInferenceGetINT4Variable(thread->proposalArgs, "verbose")) {
        printf("Thread %i found %i clusters.\n", thread->id, kde->kmeans->k);

        sprintf(outp_name, "clustered_samples.%2.2d", thread->id);
        sprintf(outp_draws_name, "clustered_draws.%2.2d", thread->id);

        LALInferenceDumpClusteredKDE(kde, outp_name, array);
        LALInferenceDumpClusteredKDEDraws(kde, outp_draws_name, ndraws);
    }
}


/**
 * Dump draws from a KDE to file.
 *
 * Print out the samples used to estimate the distribution, along with their
 * cluster assignments, and the PDF evaluated at each sample.
 * @param[in] kde       The clustered KDE to dump the info of.
 * @param[in] outp_name The name of the output file.
 * @param[in] array     The array of samples used for the KDE (it only stores a whitened version).
 */
void LALInferenceDumpClusteredKDE(LALInferenceClusteredKDE *kde, char *outp_name, REAL8 *array) {
    FILE *outp;
    REAL8 PDF;
    INT4 i, j;

    outp = fopen(outp_name, "w");
    LALInferenceFprintParameterNonFixedHeaders(outp, kde->params);
    fprintf(outp, "cluster\tweight\tPDF\n");

    for (i=0; i<kde->kmeans->npts; i++) {
        PDF = LALInferenceKmeansPDF(kde->kmeans, array + i*kde->dimension);
        for (j=0; j<kde->dimension; j++)
            fprintf(outp, "%g\t", array[i*kde->dimension + j]);
        fprintf(outp, "%i\t%f\t%g\n", kde->kmeans->assignments[i], kde->kmeans->weights[kde->kmeans->assignments[i]], PDF);
    }
    fclose(outp);
}


/**
 * Dump clustered KDE information to file.
 *
 * Dump a requested number of draws from a clustered-KDE to file,
 * along with the value of the PDF at each point.
 * @param[in] kde        The clustered-KDE proposal to draw from.
 * @param[in] outp_name  The name of the file to write to.
 * @param[in] nSamps     The number of draws to write.
 */
void LALInferenceDumpClusteredKDEDraws(LALInferenceClusteredKDE *kde, char *outp_name, INT4 nSamps) {
    FILE *outp;
    INT4 i, j;
    REAL8 *draw, PDF;

    outp = fopen(outp_name, "w");
    LALInferenceFprintParameterNonFixedHeaders(outp, kde->params);
    fprintf(outp, "PDF\n");

    for (i=0; i<nSamps; i++) {
        draw = LALInferenceKmeansDraw(kde->kmeans);
        PDF = LALInferenceKmeansPDF(kde->kmeans, draw);
        for (j=0; j<kde->dimension; j++)
            fprintf(outp, "%g\t", draw[j]);
        fprintf(outp, "%g\n", PDF);
        XLALFree(draw);
    }
    fclose(outp);
}


/**
 * Add a KDE proposal to the KDE proposal set.
 *
 * If other KDE proposals already exist, the provided KDE is appended to the list, otherwise it is added
 * as the first of such proposals.
 * @param     propArgs The proposal arguments to be added to.
 * @param[in] kde      The proposal to be added to \a thread->cycle.
 */
void LALInferenceAddClusteredKDEProposalToSet(LALInferenceVariables *propArgs, LALInferenceClusteredKDE *kde) {
    /* If proposal doesn't already exist, add to proposal args */
    if (!LALInferenceCheckVariable(propArgs, clusteredKDEProposalName)) {
        LALInferenceAddVariable(propArgs, clusteredKDEProposalName, (void *)&kde, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_LINEAR);

    /* If proposals already exist, add to the end */
    } else {
        LALInferenceClusteredKDE *existing_kde = *((LALInferenceClusteredKDE **)LALInferenceGetVariable(propArgs, clusteredKDEProposalName));
        LALInferenceClusteredKDE *old_kde = NULL;

        /* If the first proposal has the same name, replace it */
        if (!strcmp(existing_kde->name, kde->name)) {
            old_kde = existing_kde;
            kde->next = existing_kde->next;
            LALInferenceSetVariable(propArgs, clusteredKDEProposalName, (void *)&kde);
        } else {
            while (existing_kde->next != NULL) {
                /* Replace proposal with the same name if found */
                if (!strcmp(existing_kde->next->name, kde->name)) {
                    old_kde = existing_kde->next;
                    kde->next = old_kde->next;
                    existing_kde->next = kde;
                    break;
                }
                existing_kde = existing_kde->next;
            }

            /* If a proposal was not replaced, add the proposal to the end of the list */
            existing_kde->next=kde;
        }

        LALInferenceDestroyClusteredKDEProposal(old_kde);
    }

    return;
}


/**
 * Destroy an existing clustered-KDE proposal.
 *
 * Convenience function for freeing a clustered KDE proposal that
 * already exists.  This is particularly useful for a proposal that
 * is updated during a run.
 * @param proposal The proposal to be destroyed.
 */
void LALInferenceDestroyClusteredKDEProposal(LALInferenceClusteredKDE *proposal) {
    if (proposal != NULL) {
        LALInferenceClearVariables(proposal->params);
        LALInferenceKmeansDestroy(proposal->kmeans);
        XLALFree(proposal->params);
    }
    return;
}


/**
 * Setup a clustered-KDE proposal from the differential evolution buffer.
 *
 * Reads the samples currently in the differential evolution buffer and construct a
 * jump proposal from its clustered kernel density estimate.
 * @param thread The LALInferenceThreadState to get the buffer from and add the proposal to.
 */
void LALInferenceSetupClusteredKDEProposalFromDEBuffer(LALInferenceThreadState *thread) {
    INT4 i;

    /* If ACL can be estimated, thin DE buffer to only have independent samples */
    REAL8 bufferSize = (REAL8) thread->differentialPointsLength;
    REAL8 effSampleSize = (REAL8) LALInferenceComputeEffectiveSampleSize(thread);

    /* Correlations wont effect the proposal much, so floor is taken instead of ceil
     * when determining the step size */
    INT4 step = 1;
    if (effSampleSize > 0)
        step = (INT4) floor(bufferSize/effSampleSize);

    if (step == 0)
        step = 1;
    INT4 nPoints = (INT4) ceil(bufferSize/(REAL8)step);

    /* Get points to be clustered from the differential evolution buffer. */
    INT4 nPar = LALInferenceGetVariableDimensionNonFixed(thread->currentParams);
    REAL8** DEsamples = (REAL8**) XLALCalloc(nPoints, sizeof(REAL8*));
    REAL8*  temp = (REAL8*) XLALCalloc(nPoints * nPar, sizeof(REAL8));
    for (i=0; i < nPoints; i++)
      DEsamples[i] = temp + (i*nPar);

    LALInferenceThinnedBufferToArray(thread, DEsamples, step);

    /* Check if imposing cyclic reflective bounds */
    INT4 cyclic_reflective = LALInferenceGetINT4Variable(thread->proposalArgs, "cyclic_reflective_kde");

    INT4 ntrials = 5;
    LALInferenceSetupClusteredKDEProposalFromRun(thread, DEsamples[0], nPoints, cyclic_reflective, ntrials);

    /* The proposal copies the data, so the local array can be freed */
    XLALFree(temp);
    XLALFree(DEsamples);
}

/**
 * Setup a clustered-KDE proposal from the parameters in a run.
 *
 * Reads the samples currently in the differential evolution buffer and construct a
 * jump proposal from its clustered kernel density estimate.
 * @param thread The LALInferenceThreadState to get the buffer from and add the proposal to.
 * @param samples  The samples to estimate the distribution of.  Column order expected to match
 *                     the order in \a thread->currentParams.
 * @param size     Number of samples in \a samples.
 * @param cyclic_reflective Flag to check for cyclic/reflective bounds.
 * @param ntrials  Number of tirals at fixed-k to find optimal BIC
 */
void LALInferenceSetupClusteredKDEProposalFromRun(LALInferenceThreadState *thread, REAL8 *samples, INT4 size, INT4 cyclic_reflective, INT4 ntrials) {
    REAL8 weight=2.;

    /* Keep track of clustered parameter names */
    LALInferenceVariables *backwardClusterParams = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceVariables *clusterParams = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceVariableItem *item;
    for (item = thread->currentParams->head; item; item = item->next)
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name))
            LALInferenceAddVariable(backwardClusterParams, item->name, item->value, item->type, item->vary);
    for (item = backwardClusterParams->head; item; item = item->next)
        LALInferenceAddVariable(clusterParams, item->name, item->value, item->type, item->vary);

    /* Build the proposal */
    LALInferenceClusteredKDE *proposal = XLALCalloc(1, sizeof(LALInferenceClusteredKDE));
    LALInferenceInitClusteredKDEProposal(thread, proposal, samples, size, clusterParams, clusteredKDEProposalName, weight, LALInferenceOptimizedKmeans, cyclic_reflective, ntrials);

    /* Only add the kmeans was successfully setup */
    if (proposal->kmeans)
        LALInferenceAddClusteredKDEProposalToSet(thread->proposalArgs, proposal);
    else {
        LALInferenceClearVariables(clusterParams);
        XLALFree(clusterParams);
        XLALFree(proposal);
    }

    LALInferenceClearVariables(backwardClusterParams);
    XLALFree(backwardClusterParams);
}


/**
 * A proposal based on the clustered kernal density estimate of a set of samples.
 *
 * Proposes samples from the estimated distribution of a collection of points.
 * The distribution is estimated with a clustered kernel density estimator.  This
 * proposal is added to the proposal cycle with a specified weight, and in turn
 * chooses at random a KDE-estimate from a linked list.
 * @param      thread         The current LALInferenceThreadState.
 * @param      currentParams  The current parameters.
 * @param[out] proposedParams The proposed parameters.
 * @return proposal_ratio     The (log) proposal ratio for maintaining detailed balance
 */
REAL8 LALInferenceClusteredKDEProposal(LALInferenceThreadState *thread, LALInferenceVariables *currentParams, LALInferenceVariables *proposedParams) {
    REAL8 logPropRatio;

    logPropRatio = LALInferenceStoredClusteredKDEProposal(thread, currentParams, proposedParams, NULL);

    return logPropRatio;
}

/**
 * An interface to the KDE proposal that avoids a KDE evaluation if possible.
 *
 * If the value of the KDE at the current location is known, use it.  Otherwise
 * calculate and return.
 * @param      thread         The current LALInferenceThreadState.
 * @param      currentParams  The current parameters.
 * @param[out] proposedParams The proposed parameters.
 * @param      propDensity    If input is not NULL or >-DBL_MAX, assume this is the
 *                              proposal density at \a currentParams, otherwise
 *                              calculate.  It is then replaced with the proposal
 *                              density at \a proposedParams.
 * @return proposal_ratio    The (log) proposal ratio for maintaining detailed balance
 */
REAL8 LALInferenceStoredClusteredKDEProposal(LALInferenceThreadState *thread, LALInferenceVariables *currentParams, LALInferenceVariables *proposedParams, REAL8 *propDensity) {
    REAL8 cumulativeWeight, totalWeight;
    REAL8 logPropRatio = 0.0;

    LALInferenceVariableItem *item;
    LALInferenceVariables *propArgs = thread->proposalArgs;

    if (!LALInferenceCheckVariable(propArgs, clusteredKDEProposalName)) {
        LALInferenceClearVariables(proposedParams);
        return logPropRatio; /* Quit now, since there is no proposal to call */
    }

    LALInferenceCopyVariables(currentParams, proposedParams);

    /* Clustered KDE estimates are stored in a linked list, with possibly different weights */
    LALInferenceClusteredKDE *kdes = *((LALInferenceClusteredKDE **)LALInferenceGetVariable(propArgs, clusteredKDEProposalName));

    totalWeight = 0.;
    LALInferenceClusteredKDE *kde = kdes;
    while (kde!=NULL) {
        totalWeight += kde->weight;
        kde = kde->next;
    }

    /* If multiple KDE estimates exists, draw one at random */
    REAL8 randomDraw = gsl_rng_uniform(thread->GSLrandom);

    kde = kdes;
    cumulativeWeight = kde->weight;
    while(cumulativeWeight/totalWeight < randomDraw) {
        kde = kde->next;
        cumulativeWeight += kde->weight;
    }

    /* Draw a sample and fill the proposedParams variable with the parameters described by the KDE */
    REAL8 *current = XLALCalloc(kde->dimension, sizeof(REAL8));
    REAL8 *proposed = LALInferenceKmeansDraw(kde->kmeans);

    INT4 i=0;
    for (item = kde->params->head; item; item = item->next) {
        if (LALInferenceCheckVariableNonFixed(kde->params, item->name)) {
            current[i] = *(REAL8 *) LALInferenceGetVariable(currentParams, item->name);
            LALInferenceSetVariable(proposedParams, item->name, &(proposed[i]));
            i++;
        }
    }

    /* Calculate the proposal ratio */
    REAL8 logCurrentP;
    if (propDensity == NULL || *propDensity == -DBL_MAX)
        logCurrentP = LALInferenceKmeansPDF(kde->kmeans, current);
    else
        logCurrentP = *propDensity;

    REAL8 logProposedP = LALInferenceKmeansPDF(kde->kmeans, proposed);

    logPropRatio = logCurrentP - logProposedP;

    if (propDensity != NULL)
        *propDensity = logProposedP;

    XLALFree(current);
    XLALFree(proposed);

    return logPropRatio;
}


/**
 * A wrapper for the KDE proposal that doesn't store KDE evaluations.
 *
 */


/**
 * Compute the maximum ACL from the differential evolution buffer.
 *
 * Given the current differential evolution buffer, the maximum
 * one-dimensional autocorrelation length is found.
 * @param thread The thread state containing the differential evolution buffer.
 * @param maxACL UNDOCUMENTED
*/
void LALInferenceComputeMaxAutoCorrLenFromDE(LALInferenceThreadState *thread, INT4* maxACL) {
    INT4 nPar, nPoints, nSkip;
    INT4 i;
    REAL8** DEarray;
    REAL8*  temp;
    REAL8 max_acl;

    nPar = LALInferenceGetVariableDimensionNonFixed(thread->currentParams);
    nPoints = thread->differentialPointsLength;

    /* Determine the number of iterations between each entry in the DE buffer */
    nSkip = thread->differentialPointsSkip;

    /* Prepare 2D array for DE points */
    DEarray = (REAL8**) XLALCalloc(nPoints, sizeof(REAL8*));
    temp = (REAL8*) XLALCalloc(nPoints * nPar, sizeof(REAL8));
    for (i=0; i < nPoints; i++)
        DEarray[i] = temp + (i*nPar);

    LALInferenceBufferToArray(thread, DEarray);
    max_acl = nSkip * LALInferenceComputeMaxAutoCorrLen(DEarray[nPoints/2], nPoints-nPoints/2, nPar);

    if (max_acl == INFINITY)
        max_acl = INT_MAX;
    else if (max_acl < 1.0)
        max_acl = 1.0;

    *maxACL = (INT4)max_acl;
    XLALFree(temp);
    XLALFree(DEarray);
}

/**
 * Compute the maximum single-parameter autocorrelation length.  Each
 * parameter's ACL is the smallest s such that
 *
 * 1 + 2*ACF(1) + 2*ACF(2) + ... + 2*ACF(M*s) < s,
 *
 * The parameter M controls the length of the window over which we sum
 * the ACF to estimate the ACL; there must be at least M*ACL samples
 * summed.  This ensures that we obtain a reliable estimate of the ACL
 * by incorporating lags that are much longer that the estimated ACL.
 *
 * The maximum window length is also restricted to be N/K as a safety
 * precaution against relying on data near the extreme of the lags in
 * the ACF, where there is a lot of noise.
 *
 * By default, safe parameters are M = 5, K = 2.
 *
 * If no estimate can be obtained, then return Infinity.
 *
 * @param array Array with rows containing samples.
 * @param nPoints UNDOCUMENTED
 * @param nPar UNDOCUMENTED
 * @return The maximum one-dimensional autocorrelation length
*/
REAL8 LALInferenceComputeMaxAutoCorrLen(REAL8 *array, INT4 nPoints, INT4 nPar) {
    INT4 M=5, K=2;

    REAL8 mean, ACL, ACF, maxACL=0;
    INT4 par=0, lag=0, i=0, imax;
    REAL8 cumACF, s;

    if (nPoints > 1) {
        imax = nPoints/K;

        for (par=0; par<nPar; par++) {
            mean = gsl_stats_mean(array+par, nPar, nPoints);
            for (i=0; i<nPoints; i++)
                array[i*nPar + par] -= mean;

            lag=1;
            ACL=1.0;
            ACF=1.0;
            s=1.0/(REAL8)M;
            cumACF=1.0;
            while (cumACF >= s) {
                ACF = gsl_stats_correlation(array + par, nPar, array + lag*nPar + par, nPar, nPoints-lag);
                cumACF += 2.0 * ACF;
                lag++;
                s = (REAL8)lag/(REAL8)M;
                if (lag > imax) {
                    return INFINITY; /* Short circuit: this parameter has indeterminate ACL */
                }
            }
            ACL = cumACF;
            if (ACL>maxACL)
                maxACL=ACL;

            for (i=0; i<nPoints; i++)
                array[i*nPar + par] += mean;
        }
    } else {
        maxACL = INFINITY;
    }

    return maxACL;
}

/**
 * Update the estimatate of the autocorrelation length.
 *
 * @param      thread      The current LALInferenceThreadState.
*/
void LALInferenceUpdateMaxAutoCorrLen(LALInferenceThreadState *thread) {
  // Calculate ACL with latter half of data to avoid ACL overestimation from chain equilibrating after adaptation
  INT4 acl;

  LALInferenceComputeMaxAutoCorrLenFromDE(thread, &acl);

  LALInferenceSetVariable(thread->proposalArgs, "acl", &acl);
}

/**
 * Determine the effective sample size based on the DE buffer.
 *
 * Compute the number of independent samples in the differential evolution
 * buffer.
 * @param      thread      The current LALInferenceThreadState.
 */
INT4 LALInferenceComputeEffectiveSampleSize(LALInferenceThreadState *thread) {
    /* Update the ACL estimate, assuming a thinned DE buffer if ACL isn't available */
    INT4 acl = 1;
    if (LALInferenceCheckVariable(thread->proposalArgs, "acl")) {
        LALInferenceUpdateMaxAutoCorrLen(thread);
        acl = LALInferenceGetINT4Variable(thread->proposalArgs, "acl");
    }

    /* Estimate the total number of samples post-burnin based on samples in DE buffer */
    INT4 nPoints =  thread->differentialPointsLength * thread->differentialPointsSkip;
    INT4 iEff = nPoints/acl;
    return iEff;
}


INT4 LALInferencePrintProposalTrackingHeader(FILE *fp,LALInferenceVariables *params) {
      fprintf(fp, "proposal\t");
      LALInferenceFprintParameterNonFixedHeaders(fp, params);
      LALInferenceFprintParameterNonFixedHeadersWithSuffix(fp, params, "p");
      fprintf(fp, "prop_ratio\taccepted\t");
      fprintf(fp, "\n");
      return 0;
}

void LALInferencePrintProposalTracking(FILE *fp, LALInferenceProposalCycle *cycle, LALInferenceVariables *theta, LALInferenceVariables *theta_prime, REAL8 logPropRatio, INT4 accepted){
  fprintf(fp, "%s\t", cycle->proposals[cycle->counter]->name);
  LALInferencePrintSampleNonFixed(fp, theta);
  LALInferencePrintSampleNonFixed(fp, theta_prime);
  fprintf(fp, "%9.5f\t", exp(logPropRatio));
  fprintf(fp, "%d\t", accepted);
  fprintf(fp, "\n");
  return;
}

REAL8 LALInferenceSplineCalibrationProposal(LALInferenceThreadState *thread, LALInferenceVariables *currentParams, LALInferenceVariables *proposedParams) {
  char **ifo_names;
  INT4 ifo;
  INT4 nifo = LALInferenceGetINT4Variable(thread->proposalArgs, "nDet");

  LALInferenceCopyVariables(currentParams, proposedParams);

  ifo_names = *(char ***)LALInferenceGetVariable(thread->proposalArgs, "detector_names");
  for (ifo=0; ifo<nifo; ifo++) {
    UINT4 i;

    char ampName[VARNAME_MAX];
    char phaseName[VARNAME_MAX];
    REAL8 dummy;

    REAL8 ampWidth;
    REAL8 phaseWidth;
    UINT4 nspl = LALInferenceGetUINT4Variable(proposedParams, "spcal_npts");
    for (i = 0; i < nspl; i++) {
      snprintf(ampName, VARNAME_MAX, "%s_spcal_amp_%i", ifo_names[ifo], i);
      snprintf(phaseName, VARNAME_MAX, "%s_spcal_phase_%i", ifo_names[ifo], i);

      LALInferenceGetGaussianPrior(thread->priorArgs, ampName, &dummy, &ampWidth);
      REAL8 amp = LALInferenceGetREAL8Variable(proposedParams, ampName);
      amp += ampWidth*gsl_ran_ugaussian(thread->GSLrandom)/sqrt(nifo*(REAL8)nspl);
      LALInferenceSetREAL8Variable(proposedParams, ampName, amp);

      LALInferenceGetGaussianPrior(thread->priorArgs, phaseName, &dummy, &phaseWidth);
      REAL8 ph = LALInferenceGetREAL8Variable(proposedParams, phaseName);
      ph += phaseWidth*gsl_ran_ugaussian(thread->GSLrandom)/sqrt(nifo*(REAL8)nspl);
      LALInferenceSetREAL8Variable(proposedParams, phaseName, ph);
    }
  };

  return 0.0;
}
