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
#include <lal/LALInference.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>


#include <lal/LALStdlib.h>




//Test LALPriorFunction
//REAL8 PTUniformLALPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
/****************************************/
/* Returns unnormalized (!),            */
/* logarithmic (!) prior density.      	*/
/****************************************/
/* Assumes the following parameters	*/
/* exist (e.g., for TaylorT1):		*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,		*/
/* desclination, polarisation, distance.*/
/* Prior is flat if within range	*/
/****************************************/
//{
//	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;	
//	REAL8 logdensity;
//	
//	mc   = *(REAL8*) LALInferenceGetVariable(params, "chirpmass");		/* solar masses*/
//	eta  = *(REAL8*) LALInferenceGetVariable(params, "massratio");		/* dim-less    */
//	iota = *(REAL8*) LALInferenceGetVariable(params, "inclination");		/* radian      */
//	tc   = *(REAL8*) LALInferenceGetVariable(params, "time");			/* GPS seconds */
//	phi  = *(REAL8*) LALInferenceGetVariable(params, "phase");		/* radian      */
//	ra   = *(REAL8*) LALInferenceGetVariable(params, "rightascension");	/* radian      */
//	dec  = *(REAL8*) LALInferenceGetVariable(params, "declination");		/* radian      */
//	psi  = *(REAL8*) LALInferenceGetVariable(params, "polarisation"); 	/* radian      */
//	dist = *(REAL8*) LALInferenceGetVariable(params, "distance");		/* Mpc         */
//	
//	if(mc>2.41 && mc<=9.64 && eta>0.03 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
//	   && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_PI && dist>0.0 && dist<=100.0
//	   && tc>=968654557.90 && tc<=968654558.20)	
//		logdensity = 0.0;
//	else
//		logdensity = -DBL_MAX;
//	//TODO: should be properly normalized; pass in range via priorArgs?	
//	
//	return(logdensity);
//}

static void PTMCMCCombinedProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams,
                                   LALInferenceProposalFunction *props[], REAL8 weights[]) {
  REAL8 totalWeight, uniformRand;
  INT4 NProps, i;

  NProps=0;
  while (props[NProps] != NULL) NProps++;

  totalWeight = 0.0;
  for (i = 0; i < NProps; i++) {
    totalWeight += weights[i];
  }

  uniformRand = gsl_rng_uniform(runState->GSLrandom);

  i = 0;
  while (uniformRand > weights[i] / totalWeight) {
    uniformRand -= weights[i]/totalWeight;
    i++;
  }

  (props[i])(runState, proposedParams);

  return;
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
  PTMCMCLALProposal(runState,proposedParams);
  /* Restore currentParams */
  runState->currentParams=currentParamsBackup;
}

void PTMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
	UINT4 nIFO=0;
	LALInferenceIFOData *ifo=runState->data;
	REAL8 BLOCKFRAC=0.0, /* Removed block jumps because of
                                roundoff problems in cholesky
                                decomp. */
          SINGLEFRAC=1.0,
          //SKYFRAC=0.0, /* Not symmetric! */
          INCFRAC=0.0,
          PHASEFRAC=0.0,
          SKYLOCSMALLWANDERFRAC=0.0; /* Not symmetric! Was: 0.05; */
        /* No spin rotations, because they are actually not symmetric! */
        REAL8 SPINROTFRAC = 0.0; /* (runState->template == &templateLALSTPN ? 0.05 : 0.0); */
        REAL8 COVEIGENFRAC;
        REAL8 IOTADISTANCEFRAC=0.0; /* Not symmetric! Stop! */
        REAL8 DIFFFULLFRAC;
        REAL8 DIFFPARTIALFRAC;
        REAL8 PRIORFRAC=0.05;
        ProcessParamsTable *ppt;
        
        if(LALInferenceCheckVariable(proposedParams,"inclination")) INCFRAC=0.05;
        if(LALInferenceCheckVariable(proposedParams,"phase")) PHASEFRAC=0.05;
          
        ppt=LALInferenceGetProcParamVal(runState->commandLine, "--iotaDistance");
        if (ppt) {
          IOTADISTANCEFRAC = atof(ppt->value);
        }
        ppt=LALInferenceGetProcParamVal(runState->commandLine, "--covarianceMatrix");
        if (ppt) {
          COVEIGENFRAC = 1.0;
        } else {
          COVEIGENFRAC = 0.0;
        }

        if (LALInferenceGetProcParamVal(runState->commandLine, "--differential-evolution")) {
          DIFFFULLFRAC = 1.0;
          DIFFPARTIALFRAC = 1.0 / 4.0;
        } else {
          DIFFFULLFRAC = 0.0;
          DIFFPARTIALFRAC = 0.0;
        }

        nIFO = 0;
        while(ifo){ifo=ifo->next; nIFO++;}
        
        if (nIFO < 2) {
          REAL8 weights[] = {BLOCKFRAC, SINGLEFRAC, INCFRAC, PHASEFRAC, SPINROTFRAC, COVEIGENFRAC, SKYLOCSMALLWANDERFRAC, IOTADISTANCEFRAC, DIFFFULLFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, PRIORFRAC};
          LALInferenceProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          &PTMCMCLALInferenceDrawUniformlyFromPrior,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
          return;
        } else if (nIFO < 3) {
          /* Removed the rotate sky function from proposal because it's not symmetric. */
          REAL8 weights[] = {BLOCKFRAC, SINGLEFRAC, 0.0 /* SKYFRAC */, INCFRAC, PHASEFRAC, SPINROTFRAC, COVEIGENFRAC, SKYLOCSMALLWANDERFRAC, IOTADISTANCEFRAC, DIFFFULLFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, DIFFPARTIALFRAC, PRIORFRAC};
          LALInferenceProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceRotateSky,
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          &PTMCMCLALInferenceDrawUniformlyFromPrior,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
        } else {
          /* Removed the rotate sky function because it's not symmetric. */
          REAL8 weights[] = { BLOCKFRAC, 
                              SINGLEFRAC, 
                              0.0 /* SKYFRAC */, 
                              0.0 /* SKYFRAC */, 
                              INCFRAC, 
                              PHASEFRAC, 
                              SPINROTFRAC, 
                              COVEIGENFRAC, 
                              SKYLOCSMALLWANDERFRAC, 
                              IOTADISTANCEFRAC, 
                              DIFFFULLFRAC, 
                              DIFFPARTIALFRAC, 
                              DIFFPARTIALFRAC, 
                              DIFFPARTIALFRAC, 
                              DIFFPARTIALFRAC,  
                              PRIORFRAC};
          LALInferenceProposalFunction *props[] = {&PTMCMCLALBlockCorrelatedProposal,
                                          &PTMCMCLALSingleAdaptProposal,
                                          &PTMCMCLALInferenceRotateSky,
                                          (LALInferenceProposalFunction *)(&PTMCMCLALInferenceReflectDetPlane),
                                          &PTMCMCLALInferenceInclinationFlip,
                                          &PTMCMCLALInferenceOrbitalPhaseJump,
                                          &PTMCMCLALInferenceRotateSpins,
                                          &PTMCMCLALInferenceCovarianceEigenvectorJump,
                                          &PTMCMCLALInferenceSkyLocWanderJump,
                                          &PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump,
                                          &PTMCMCLALInferenceDifferentialEvolutionFull,
                                          &PTMCMCLALInferenceDifferentialEvolutionMasses,
                                          &PTMCMCLALInferenceDifferentialEvolutionAmp,
                                          &PTMCMCLALInferenceDifferentialEvolutionSpins,
                                          &PTMCMCLALInferenceDifferentialEvolutionSky,
                                          &PTMCMCLALInferenceDrawUniformlyFromPrior,
                                          0};
          PTMCMCCombinedProposal(runState, proposedParams, props, weights);
        }
}

void PTMCMCLALBlockProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALInferenceVariableItem *paraHead=NULL;
	INT4 i;
	LALInferenceCopyVariables(runState->currentParams, proposedParams);

        REAL8 T = *(REAL8 *)LALInferenceGetVariable(runState->proposalArgs, "temperature");
	
	REAL8 sigma = 0.1*sqrt(T); /* Adapt to temperature. */
	REAL8 big_sigma = 1.0;
	
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
	
	/* loop over all parameters */
	for (paraHead=proposedParams->head,i=0; paraHead; paraHead=paraHead->next)
	{ 
		if(paraHead->vary==LALINFERENCE_PARAM_LINEAR || paraHead->vary==LALINFERENCE_PARAM_CIRCULAR){
			
			if (!strcmp(paraHead->name,"massratio") || !strcmp(paraHead->name,"time") || !strcmp(paraHead->name,"a_spin2") || !strcmp(paraHead->name,"a_spin1")){
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
			}else if (!strcmp(paraHead->name,"polarisation") || !strcmp(paraHead->name,"phase") || !strcmp(paraHead->name,"inclination")){
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
			}else{
				*(REAL8 *)paraHead->value += gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
			}
			i++;
		}
	}
	
	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
	
}

void PTMCMCLALSingleAdaptProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *args = runState->proposalArgs;
  ProcessParamsTable *ppt = LALInferenceGetProcParamVal(runState->commandLine, "--adapt");
  
  if (!LALInferenceCheckVariable(args, SIGMAVECTORNAME) || !ppt) {
    /* We are not adaptive, or for some reason don't have a sigma
       vector---fall back on old proposal. */
    PTMCMCLALSingleProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    LALInferenceVariableItem *param = NULL, *dummyParam = NULL;
    REAL8 T = *(REAL8 *)LALInferenceGetVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    UINT4 dim;
    UINT4 i;
    UINT4 varNr;
    REAL8Vector *sigmas = *(REAL8Vector **) LALInferenceGetVariable(args, SIGMAVECTORNAME);

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

    INT4 as = 1;
    LALInferenceSetVariable(args, "adaptableStep", &as);

    LALInferenceSetVariable(args, "proposedVariableNumber", &varNr);
    
    LALInferenceSetVariable(args, "proposedArrayNumber", &i);
  }
}

void PTMCMCLALSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
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
  
  INT4 as = 1;
  LALInferenceSetVariable(runState->proposalArgs, "adaptableStep", &as);
  
  LALInferenceSetVariable(runState->proposalArgs, "proposedVariableNumber", &varNr);
  
  LALInferenceSetVariable(runState->proposalArgs, "proposedArrayNumber", &i);
  
}

void PTMCMCLALBlockCorrelatedProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  LALInferenceVariables *args = runState->proposalArgs;

  if (!LALInferenceCheckVariable(args, COVMATRIXNAME) || !LALInferenceCheckVariable(args, UNCORRSAMPNAME)) {
    /* No correlation matrix! */
    PTMCMCLALBlockProposal(runState, proposedParams);
  } else {
    gsl_rng *rng = runState->GSLrandom;
    REAL8 T = *(REAL8 *)LALInferenceGetVariable(args, "temperature");
    REAL8 sqrtT = sqrt(T);
    gsl_matrix *covarianceMatrix = *(gsl_matrix **)LALInferenceGetVariable(args, COVMATRIXNAME);
    REAL8Vector *uncorrelatedSample = *(REAL8Vector **)LALInferenceGetVariable(args, UNCORRSAMPNAME);
    UINT4 i;
    UINT4 N = uncorrelatedSample->length;
    LALInferenceVariableItem *param = NULL;

    LALInferenceCopyVariables(runState->currentParams, proposedParams);

    if (covarianceMatrix->size1 != N || covarianceMatrix->size2 != N) {
      fprintf(stderr, "ERROR: covariance matrix and sample vector sizes do not agree (in %s, line %d)\n",
              __FILE__, __LINE__);
      exit(1);
    }
    
    for (i = 0; i < N; i++) {
      uncorrelatedSample->data[i] = gsl_ran_ugaussian(rng)*sqrtT; /* Normalized to magnitude sqrt(T) */
    }

    for (i = 0, param = proposedParams->head; param != NULL; param = param->next) {
      if (param->vary != LALINFERENCE_PARAM_FIXED && param->vary != LALINFERENCE_PARAM_OUTPUT) {
        /* Then it's a parameter to set. */
        UINT4 j;
        REAL8 sum;

        for (j = 0, sum = 0.0; j < N; j++) {
          sum += gsl_matrix_get(covarianceMatrix, i, j)*uncorrelatedSample->data[j];
        }

        if (param->type != LALINFERENCE_REAL8_t) {
          fprintf(stderr, "Trying to use covariance matrix to set non-REAL8 parameter (in %s, line %d)\n",
                  __FILE__, __LINE__);
          exit(1);
        } else {
          *((REAL8 *)param->value) += sum;
        }
        
        i++;
      }
    }

    LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
  }
}

/* Reflect the inclination about the observing plane, iota -> Pi - iota */
void PTMCMCLALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  REAL8 phi;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);
  
  phi = *((REAL8 *) LALInferenceGetVariable(proposedParams, "phase"));

  phi = fmod(phi+M_PI, 2.0*M_PI);

  LALInferenceSetVariable(proposedParams, "phase", &phi);

  /* Probably not needed, but play it safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}


//Test LALInferenceProposalFunction
//void PTMCMCLALInferenceProposaltemp(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
/****************************************/
/* Assumes the following parameters		*/
/* exist (e.g., for TaylorT1):			*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,			*/
/* desclination, polarisation, distance.*/
/* Simply picks a new value based on	*/
/* fixed Gaussian;						*/
/* need smarter wall bounces in future.	*/
/****************************************/
//{
//	REAL8 mc, eta, iota, phi, tc, ra, dec, psi, dist;
//	REAL8 a_spin1, a_spin2, theta_spin1, theta_spin2, phi_spin1, phi_spin2;
//	REAL8 mc_proposed, eta_proposed, iota_proposed, phi_proposed, tc_proposed, 
//	ra_proposed, dec_proposed, psi_proposed, dist_proposed;
//	REAL8 a_spin1_proposed, a_spin2_proposed, theta_spin1_proposed, 
//	theta_spin2_proposed, phi_spin1_proposed, phi_spin2_proposed;
//	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
//	gsl_rng * GSLrandom=runState->GSLrandom;
//	LALInferenceVariables * currentParams = runState->currentParams;
//	LALInferenceCopyVariables(currentParams, proposedParams);
//	
//	REAL8 sigma = 0.1;
//	REAL8 big_sigma = 1.0;
//	
//	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-3) big_sigma = 1.0e1;    //Every 1e3 iterations, take a 10x larger jump in all parameters
//	if(gsl_ran_ugaussian(GSLrandom) < 1.0e-4) big_sigma = 1.0e2;    //Every 1e4 iterations, take a 100x larger jump in all parameters
//	
//	
//	mc   = *(REAL8*) LALInferenceGetVariable(currentParams, "chirpmass");		/* solar masses*/
//	eta  = *(REAL8*) LALInferenceGetVariable(currentParams, "massratio");		/* dim-less    */
//	iota = *(REAL8*) LALInferenceGetVariable(currentParams, "inclination");		/* radian      */
//	tc   = *(REAL8*) LALInferenceGetVariable(currentParams, "time");				/* GPS seconds */
//	phi  = *(REAL8*) LALInferenceGetVariable(currentParams, "phase");			/* radian      */
//	ra   = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension");	/* radian      */
//	dec  = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");		/* radian      */
//	psi  = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");		/* radian      */
//	dist = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");			/* Mpc         */
//	
//	if (LALInferenceCheckVariable(currentParams, "a_spin1")){
//		a_spin1 = *(REAL8*) LALInferenceGetVariable(currentParams, "a_spin1");
//		a_spin1_proposed = a_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
//		LALInferenceSetVariable(proposedParams, "a_spin1",      &a_spin1_proposed);
//	}
//	if (LALInferenceCheckVariable(currentParams, "theta_spin1")){
//		theta_spin1 = *(REAL8*) LALInferenceGetVariable(currentParams, "theta_spin1");
//		theta_spin1_proposed = theta_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		LALInferenceSetVariable(proposedParams, "theta_spin1",      &theta_spin1_proposed);
//	}
//	if (LALInferenceCheckVariable(currentParams, "phi_spin1")){
//		phi_spin1 = *(REAL8*) LALInferenceGetVariable(currentParams, "phi_spin1");
//		phi_spin1_proposed = phi_spin1 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		LALInferenceSetVariable(proposedParams, "phi_spin1",      &phi_spin1_proposed);
//	}
//	if (LALInferenceCheckVariable(currentParams, "a_spin2")){
//		a_spin2 = *(REAL8*) LALInferenceGetVariable(currentParams, "a_spin2");
//		a_spin2_proposed = a_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001;
//		LALInferenceSetVariable(proposedParams, "a_spin2",      &a_spin2_proposed);
//	}
//	if (LALInferenceCheckVariable(currentParams, "theta_spin2")){
//		theta_spin2 = *(REAL8*) LALInferenceGetVariable(currentParams, "theta_spin2");
//		theta_spin2_proposed = theta_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		LALInferenceSetVariable(proposedParams, "theta_spin2",      &theta_spin2_proposed);
//	}
//	if (LALInferenceCheckVariable(currentParams, "phi_spin2")){
//		phi_spin2 = *(REAL8*) LALInferenceGetVariable(currentParams, "phi_spin2");
//		phi_spin2_proposed = phi_spin2 + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//		LALInferenceSetVariable(proposedParams, "phi_spin2",      &phi_spin2_proposed);
//	}
//
//	//mc_proposed   = mc*(1.0+gsl_ran_ugaussian(GSLrandom)*0.01);	/*mc changed by 1% */
//	// (above proposal is not symmetric!)
//	mc_proposed   = mc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;	/*mc changed by 0.0001 */
//	//mc_proposed   = mc * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.01);          /* mc changed by ~0.1% */
//	
//	eta_proposed  = eta  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*eta changed by 0.01*/
//	//TODO: if(eta_proposed>0.25) eta_proposed=0.25-(eta_proposed-0.25); etc.
//	iota_proposed = iota + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.5;
//	tc_proposed   = tc   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.001; /*time changed by 5 ms*/
//	phi_proposed  = phi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
//	ra_proposed   = ra   + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//	dec_proposed  = dec  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.01;
//	psi_proposed  = psi  + gsl_ran_ugaussian(GSLrandom)*big_sigma*sigma*0.1;
//	//dist_proposed = dist + gsl_ran_ugaussian(GSLrandom)*0.5;
//	dist_proposed = dist * exp(gsl_ran_ugaussian(GSLrandom)*sigma*0.1); // ~10% change
//	
//	
//	
//	LALInferenceSetVariable(proposedParams, "chirpmass",      &mc_proposed);		
//	LALInferenceSetVariable(proposedParams, "massratio",      &eta_proposed);
//	LALInferenceSetVariable(proposedParams, "inclination",    &iota_proposed);
//	LALInferenceSetVariable(proposedParams, "phase",          &phi_proposed);
//	LALInferenceSetVariable(proposedParams, "time",           &tc_proposed); 
//	LALInferenceSetVariable(proposedParams, "rightascension", &ra_proposed);
//	LALInferenceSetVariable(proposedParams, "declination",    &dec_proposed);
//	LALInferenceSetVariable(proposedParams, "polarisation",   &psi_proposed);
//	LALInferenceSetVariable(proposedParams, "distance",       &dist_proposed);
//	
//	LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
//	
//	dist_proposed = *(REAL8*) LALInferenceGetVariable(proposedParams, "distance");
//	logProposalRatio *= dist_proposed / dist;
//	//logProposalRatio *= mc_proposed / mc;   // (proposal ratio for above "scaled log-normal" proposal)
//	
//	// return ratio of proposal densities (for back & forth jumps) 
//	// in "runState->proposalArgs" vector:
//	if (LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio"))
//		LALInferenceSetVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
//	else
//		LALInferenceAddVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
//}


/*REAL8 GaussianLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, LALTemplateFunction *template)
{
	
	double result=0.0;
	double sumsq=0.0;
	//double norm=0.0;
	//int i=0;
	double x[20];
	double xmax=0.0;
	double deltax=0.01;
	
	x[0]=*(REAL8 *)LALInferenceGetVariable(currentParams,"x0");
	//for(i=0;i<run.nMCMCpar;i++){
	//	x[i]= par->par[run.parRevID[185+i]];
	//}
	
//	for(i=0;i<run.nMCMCpar;i++){
	//	sumsq+=(x[i]-xmax)*(x[i]-xmax)/(2*deltax);
		//norm+=-0.91893853320468-log(sqrt(deltax));
//	}
	sumsq=(x[0]-xmax)*(x[0]-xmax)/(2*deltax);
    //result=log(100*exp(-sumsq));
	//result=15/(2*deltax)-sumsq;
	result=1.0/(2.0*deltax)-sumsq;
	return result;

}*/

/*REAL8 UnityLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data, LALTemplateFunction *template)
{
	return 1.0;
}*/



/*REAL8 PTUniformGaussianPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{

	REAL8 x0;	
	REAL8 logdensity;
	
	x0   = *(REAL8*) LALInferenceGetVariable(params, "x0");

	if(x0>=-1.0 && x0<=1.0)	
		logdensity = 0.0;
	else
		logdensity = -DBL_MAX;
	//TODO: should be properly normalized; pass in range via priorArgs?	
	
	return(logdensity);
	
}*/

/*void PTMCMCGaussianProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams)
{
	
	REAL8 x0;
	REAL8 x0_proposed;
	REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
	gsl_rng * GSLrandom=runState->GSLrandom;
	LALInferenceVariables * currentParams = runState->currentParams;
	REAL8 sigma = 0.1;
	
	x0   = *(REAL8*) LALInferenceGetVariable(currentParams, "x0");

	x0_proposed   = x0 + gsl_ran_ugaussian(GSLrandom)*sigma;
	//logProposalRatio *= x0_proposed / x0;   // (proposal ratio for above "scaled log-normal" proposal)

	
	LALInferenceCopyVariables(currentParams, proposedParams);
	LALInferenceSetVariable(proposedParams, "x0",      &(x0_proposed));		

	
	// return ratio of proposal densities (for back & forth jumps) 
	// in "runState->proposalArgs" vector:
	if (LALInferenceCheckVariable(runState->proposalArgs, "logProposalRatio"))
		LALInferenceSetVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio);
	else
		LALInferenceAddVariable(runState->proposalArgs, "logProposalRatio", &logProposalRatio, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
	
	
	
}*/



void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude);
void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude)
{
	vec[0]=cos(longitude)*cos(latitude);
	vec[1]=sin(longitude)*cos(latitude);
	vec[2]=sin(latitude);
	return;
}

void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude);
void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude)
{
	REAL8 longi,lat,dist;
	dist=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	/*XLALMCMCSetParameter(parameter,"distMpc",dist);*/
	longi=atan2(pos[1]/dist,pos[0]/dist);
	if(longi<0.0) longi=LAL_TWOPI+longi;
	lat=asin(pos[2]/dist);
	*longitude=longi;
	*latitude=lat;	
	return;
}

void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);
void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3])
{
	out[0]=x[1]*y[2] - x[2]*y[1];
	out[1]=y[0]*x[2] - x[0]*y[2];
	out[2]=x[0]*y[1] - x[1]*y[0];
	return;
}

void normalise(REAL8 vec[3]);
void normalise(REAL8 vec[3]){
	REAL8 my_abs=0.0;
	my_abs=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
	vec[0]/=my_abs;
	vec[1]/=my_abs;
	vec[2]/=my_abs;
	return;
}


void PTMCMCLALInferenceRotateSky(
						   LALInferenceRunState *state,
						   LALInferenceVariables *parameter
						   )
{ /* Function to rotate the current sample around the vector between two random detectors */
	static LALStatus status;
	INT4 IFO1,IFO2;
	REAL4 randnum;
	REAL8 vec[3];
	REAL8 cur[3];
	REAL8 longi,lat;
	REAL8 vec_abs=0.0,theta,c,s;
	UINT4 i,j;
	
	LALInferenceCopyVariables(state->currentParams, parameter);
	
	UINT4 nIFO=0;
	LALInferenceIFOData *ifodata1=state->data;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALInferenceIFOData **IFOs=calloc(nIFO,sizeof(LALInferenceIFOData *));
	for(i=0,ifodata1=state->data;i<nIFO;i++){
		IFOs[i]=ifodata1;
		ifodata1=ifodata1->next;
	}
	
	
	if(nIFO<2) return;
	if(nIFO==2 && IFOs[0]==IFOs[1]) return;
	
	longi = *(REAL8 *)LALInferenceGetVariable(parameter,"rightascension");
	lat = *(REAL8 *)LALInferenceGetVariable(parameter,"declination");
	
	/* Convert the RA/dec to geodetic coordinates, as the detectors use these */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(IFOs[0]->epoch));
	longi=geodetic.longitude;
	lat=geodetic.latitude;
	cur[0]=cos(lat)*cos(longi);
	cur[1]=cos(lat)*sin(longi);
	cur[2]=sin(lat);
	
	IFO1 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	do{ /* Pick random interferometer other than the first one */
		IFO2 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO2==IFO1 || IFOs[IFO1]->detector==IFOs[IFO2]->detector);
	
	/*	fprintf(stderr,"Rotating around %s-%s vector\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2]);*/
	/* Calc normalised direction vector */
	for(i=0;i<3;i++) vec[i]=IFOs[IFO2]->detector->location[i]-IFOs[IFO1]->detector->location[i];
	for(i=0;i<3;i++) vec_abs+=vec[i]*vec[i];
	vec_abs=sqrt(vec_abs);
	for(i=0;i<3;i++) vec[i]/=vec_abs;
	
	/* Chose random rotation angle */
	randnum=gsl_rng_uniform(state->GSLrandom);
	theta=LAL_TWOPI*randnum;
	c=cos(-theta); s=sin(-theta);
	/* Set up rotation matrix */
	double R[3][3] = {{c+vec[0]*vec[0]*(1.0-c), 
		vec[0]*vec[1]*(1.0-c)-vec[2]*s,
		vec[0]*vec[2]*(1.0-c)+vec[1]*s},
		{vec[1]*vec[0]*(1.0-c)+vec[2]*s,
			c+vec[1]*vec[1]*(1.0-c),
			vec[1]*vec[2]*(1.0-c)-vec[0]*s},
		{vec[2]*vec[0]*(1.0-c)-vec[1]*s,
			vec[2]*vec[1]*(1.0-c)+vec[0]*s,
			c+vec[2]*vec[2]*(1.0-c)}};
	REAL8 new[3]={0.0,0.0,0.0};
	for (i=0; i<3; ++i)
		for (j=0; j<3; ++j)
			new[i] += R[i][j]*cur[j];
	double newlong = atan2(new[1],new[0]);
	if(newlong<0.0) newlong=LAL_TWOPI+newlong;
	
	geodetic.longitude=newlong;
	geodetic.latitude=asin(new[2]);
	/* Convert back into equatorial (sky) coordinates */
	LALGeographicToEquatorial(&status,&equatorial,&geodetic,&(IFOs[0]->epoch));
	newlong=equatorial.longitude;
	double newlat=equatorial.latitude;
	
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	dtold = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, longi, lat, &(IFOs[0]->epoch)); /* Compute time delay */
	dtnew = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, newlong, newlat, &(IFOs[0]->epoch)); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=*(REAL8 *)LALInferenceGetVariable(parameter,"time");
	LALInferenceSetVariable(parameter,"time",&deltat);	
	LALInferenceSetVariable(parameter,"declination",&newlat);
	LALInferenceSetVariable(parameter,"rightascension",&newlong);
	/*fprintf(stderr,"Skyrotate: new pos = %lf %lf %lf => %lf %lf\n",new[0],new[1],new[2],newlong,asin(new[2]));*/
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	return;
}


INT4 PTMCMCLALInferenceReflectDetPlane(
								 LALInferenceRunState *state,
								 LALInferenceVariables *parameter
								 )
{ /* Function to reflect a point on the sky about the plane of 3 detectors */
	/* Returns -1 if not possible */
	static LALStatus status;
	UINT4 i;
	int DetCollision=0;
//	REAL4 randnum; - set but not used error 
	REAL8 longi,lat;
	REAL8 dist;
	REAL8 pos[3];
	REAL8 normal[3];
	REAL8 w1[3]; /* work vectors */
	REAL8 w2[3];
	INT4 IFO1,IFO2,IFO3;
	
	LALInferenceCopyVariables(state->currentParams, parameter);
	
	UINT4 nIFO=0;
	LALInferenceIFOData *ifodata1=state->data;
	LALInferenceIFOData *ifodata2=NULL;
	while(ifodata1){
		nIFO++;
		ifodata1=ifodata1->next;
	}
	
	LALInferenceIFOData **IFOs=calloc(nIFO,sizeof(LALInferenceIFOData *));
	if(!IFOs) {
		printf("Unable to allocate memory for %i LALInferenceIFOData *s\n",nIFO);
		exit(1);
	}
	for(i=0,ifodata1=state->data;i<nIFO;i++){
		IFOs[i]=ifodata1;
		ifodata1=ifodata1->next;
	}
	
	if(nIFO<3) return(-1) ; /* not enough IFOs to construct a plane */
	for(ifodata1=state->data;ifodata1;ifodata1=ifodata1->next)
		for(ifodata2=ifodata1->next;ifodata2;ifodata2=ifodata2->next)
			if(ifodata1->detector==ifodata2->detector) DetCollision+=1;
	
	if(nIFO-DetCollision<3) return(-1); /* Not enough independent IFOs */
	
	/* Select IFOs to use */
	IFO1=gsl_rng_uniform_int(state->GSLrandom,nIFO);
	do {
		IFO2=gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO1==IFO2 || IFOs[IFO1]==IFOs[IFO2]);
	//randnum=gsl_rng_uniform(state->GSLrandom); - set but not used
	do {
		IFO3 = gsl_rng_uniform_int(state->GSLrandom,nIFO);
	}while(IFO3==IFO1
		   || IFO3==IFO2
		   || IFOs[IFO3]==IFOs[IFO1]
		   || IFOs[IFO3]==IFOs[IFO2]);
	/*fprintf(stderr,"Using %s, %s and %s for plane\n",inputMCMC->ifoID[IFO1],inputMCMC->ifoID[IFO2],inputMCMC->ifoID[IFO3]);*/
	
	longi = *(REAL8 *)LALInferenceGetVariable(parameter,"rightascension");
	lat = *(REAL8 *)LALInferenceGetVariable(parameter,"declination");
	
	double deltalong=0;
	
	/* Convert to earth coordinates */
	SkyPosition geodetic,equatorial;
	equatorial.longitude=longi;
	equatorial.latitude=lat;
	equatorial.system=COORDINATESYSTEM_EQUATORIAL;
	geodetic.system=COORDINATESYSTEM_GEOGRAPHIC;
	LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(state->data->epoch));
	//LALEquatorialToGeographic(&status,&geodetic,&equatorial,&(state->data->epoch));
        deltalong=geodetic.longitude-equatorial.longitude;
	
	/* Add offset to RA to convert to earth-fixed */
	
	/* Calculate cartesian version of earth-fixed sky position */
	GetCartesianPos(pos,geodetic.longitude,lat); /* Get sky position in cartesian coords */
	
	
	/* calculate the unit normal vector of the detector plane */
	for(i=0;i<3;i++){ /* Two vectors in the plane */
		w1[i]=IFOs[IFO2]->detector->location[i] - IFOs[IFO1]->detector->location[i];
		w2[i]=IFOs[IFO3]->detector->location[i] - IFOs[IFO1]->detector->location[i];
	}
	crossProduct(normal,w1,w2);
	normalise(normal);
	
	/* Calculate the distance between the point and the plane n.(point-IFO1) */
	for(dist=0.0,i=0;i<3;i++) dist+=normal[i]*pos[i];
	/* Reflect the point pos across the plane */
	for(i=0;i<3;i++) pos[i]=pos[i]-2.0*dist*normal[i];
	
	REAL8 newLongGeo,newLat;
	CartesianToSkyPos(pos,&newLongGeo,&newLat);
	REAL8 newLongSky=newLongGeo-deltalong;
	
	
	LALInferenceSetVariable(parameter,"rightascension",&newLongSky);
	LALInferenceSetVariable(parameter,"declination",&newLat);
	
	/* Compute change in tgeocentre for this change in sky location */
	REAL8 dtold,dtnew,deltat;
	dtold = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, longi, lat, &(IFOs[0]->epoch)); /* Compute time delay */
	dtnew = XLALTimeDelayFromEarthCenter(IFOs[0]->detector->location, newLongSky, newLat, &(IFOs[0]->epoch)); /* Compute time delay */
	deltat=dtold-dtnew; /* deltat is change in arrival time at geocentre */
	deltat+=*(REAL8 *)LALInferenceGetVariable(parameter,"time");
	LALInferenceSetVariable(parameter,"time",&deltat);
	
	LALInferenceCyclicReflectiveBound(parameter,state->priorArgs);
	free(IFOs);
	
	return(0);
}



static REAL8 dot(REAL8 v[3], REAL8 w[3]) {
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

static REAL8 norm3(REAL8 v[3]) { return sqrt(dot(v,v)); }

static void cross(REAL8 x[3], REAL8 y[3], REAL8 z[3]) {
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}

static void rotateVectorAboutVector(REAL8 v[3], REAL8 axis[3], REAL8 theta) {
  REAL8 an = norm3(axis);
  REAL8 x[3], y[3], z[3];
  REAL8 zdotv, xnorm;
  INT4 i;

  for (i = 0; i < 3; i++) {
    z[i] = axis[i]/an; /* zhat = axisHat */
  }

  zdotv = dot(z, v);

  for (i = 0; i < 3; i++) {
    x[i] = v[i] - zdotv*z[i]; /* Remove the z component from v, store in x. */
  }

  xnorm = norm3(x);

  if (xnorm == 0.0) return; /* v is along axis, rotation is done. */

  for (i = 0; i < 3; i++) {
    x[i] /= xnorm;
  }

  cross(z, x, y);  /* y = z \times x*/

  for (i = 0; i < 3; i++) {
    v[i] = zdotv*z[i] + xnorm*(cos(theta)*x[i]+sin(theta)*y[i]);
  }
}

static void thetaPhiToVector(REAL8 norm, REAL8 theta, REAL8 phi, REAL8 v[3]) {
  v[2] = norm*cos(theta);
  v[0] = norm*sin(theta)*cos(phi);
  v[1] = norm*sin(theta)*sin(phi);
}

static void vectorToThetaPhi(REAL8 *nrm, REAL8 *theta, REAL8 *phi, REAL8 v[3]) {
  *nrm = norm3(v);

  *theta = acos(v[2]/(*nrm));

  *phi = atan2(v[1], v[0]);

  if (*phi < 0.0) {
    *phi += 2.0*M_PI; /* We use 0 <= phi < 2*M_PI, while atan2 uses -M_PI < phi <= M_PI. */
  }
}

void PTMCMCLALInferenceRotateSpins(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  REAL8 inc, theta1, theta2, phi1, phi2;
  REAL8 L[3], A1[3], A2[3];
  REAL8 rotAngle;
  REAL8 selector;
  REAL8 dummy;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  inc = *((REAL8 *)LALInferenceGetVariable(proposedParams, "inclination"));
  theta1 = *((REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin1"));
  phi1 = *((REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin1"));
  theta2 = *((REAL8 *)LALInferenceGetVariable(proposedParams, "theta_spin2"));
  phi2 = *((REAL8 *)LALInferenceGetVariable(proposedParams, "phi_spin2"));

  thetaPhiToVector(1.0, inc, 0.0, L); 
  thetaPhiToVector(1.0, theta1, phi1, A1);
  thetaPhiToVector(1.0, theta2, phi2, A2);

  rotAngle = 2.0*M_PI*gsl_rng_uniform(runState->GSLrandom);

  selector = gsl_rng_uniform(runState->GSLrandom);
  
  if (selector < 1.0/3.0) {
    /* Rotate both spins about L */
    rotateVectorAboutVector(A1, L, rotAngle);
    rotateVectorAboutVector(A2, L, rotAngle);
  } else if (selector < 2.0 / 3.0) {
    /* Rotate only A1 */
    rotateVectorAboutVector(A1, L, rotAngle);
  } else {
    /* Rotate only A2 */
    rotateVectorAboutVector(A2, L, rotAngle);
  }

  vectorToThetaPhi(&dummy, &theta1, &phi1, A1);
  vectorToThetaPhi(&dummy, &theta2, &phi2, A2);

  LALInferenceSetVariable(proposedParams, "theta_spin1", &theta1);
  LALInferenceSetVariable(proposedParams, "phi_spin1", &phi1);
  LALInferenceSetVariable(proposedParams, "theta_spin2", &theta2);
  LALInferenceSetVariable(proposedParams, "phi_spin2", &phi2);
}

void PTMCMCLALInferenceInclinationFlip(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  REAL8 iota;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  iota = *((REAL8 *) LALInferenceGetVariable(proposedParams, "inclination"));
  
//   if (runState->template==&LALInferenceTemplateLALSTPN) {
//     /* Handle spins. */
//     REAL8 dummyNorm, newIota, newPhi;
//     REAL8 theta1, theta2, phi1, phi2;
//     REAL8 L[3], a1[3], a2[3], xhat[3] = {1,0,0};
// 
//     theta1 = *((REAL8 *) LALInferenceGetVariable(proposedParams, "theta_spin1"));
//     theta2 = *((REAL8 *) LALInferenceGetVariable(proposedParams, "theta_spin2"));
// 
//     phi1 = *((REAL8 *) LALInferenceGetVariable(proposedParams, "phi_spin1"));
//     phi2 = *((REAL8 *) LALInferenceGetVariable(proposedParams, "phi_spin2"));
// 
//     thetaPhiToVector(1.0, iota, 0.0, L);
//     thetaPhiToVector(1.0, theta1, phi1, a1);
//     thetaPhiToVector(1.0, theta2, phi2, a2);
// 
//     rotateVectorAboutVector(L, xhat, M_PI-2.0*iota);
//     rotateVectorAboutVector(a1, xhat, M_PI-2.0*iota);
//     rotateVectorAboutVector(a2, xhat, M_PI-2.0*iota);
// 
//     vectorToThetaPhi(&dummyNorm, &newIota, &newPhi, L);
//     vectorToThetaPhi(&dummyNorm, &theta1, &phi1, a1);
//     vectorToThetaPhi(&dummyNorm, &theta2, &phi2, a2);
// 
//     if (fabs(newIota + iota - M_PI) > 1e-8 || fabs(newPhi) > 1e-8) {
//       fprintf(stderr, "ERROR: inclination swap not implemented properly.\n");
//       fprintf(stderr, "ERROR: should have new iota = Pi - iota, phi = 0 instead have\n");
//       fprintf(stderr, "ERROR: new iota = %g, old iota = %g, phi = %g\n", newIota, iota, newPhi);
//       exit(1);
//     }
// 
//     LALInferenceSetVariable(proposedParams, "phi_spin1", &phi1);
//     LALInferenceSetVariable(proposedParams, "phi_spin2", &phi2);
//     LALInferenceSetVariable(proposedParams, "theta_spin1", &theta1);
//     LALInferenceSetVariable(proposedParams, "theta_spin2", &theta2);
//     /* Don't need to set iota because it will happen outside the if statement. */
//   }

  iota = M_PI - iota;

  LALInferenceSetVariable(proposedParams, "inclination", &iota);

  /* Not really needed (probably), but let's be safe. */
  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs); 
}

void PTMCMCLALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
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

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void PTMCMCLALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  gsl_rng *rng = runState->GSLrandom;
  LALInferenceVariables *proposalArgs = runState->proposalArgs;
  REAL8 temp = *((REAL8 *)LALInferenceGetVariable(proposalArgs, "temperature"));
  REAL8 jumpX = sqrt(temp)*0.01*gsl_ran_ugaussian(rng)/sqrt(2.0);
  REAL8 jumpY = sqrt(temp)*0.01*gsl_ran_ugaussian(rng)/sqrt(2.0);
  REAL8 RA, DEC;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  RA = *((REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension"));
  DEC = *((REAL8 *)LALInferenceGetVariable(proposedParams, "declination"));

  RA += jumpX/cos(DEC);
  DEC += jumpY;

  LALInferenceSetVariable(proposedParams, "rightascension", &RA);
  LALInferenceSetVariable(proposedParams, "declination", &DEC);

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

/* Choose a jump in inclination and distance such that the signal
   amplitude in one of the detectors remains constant. */
void PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump(LALInferenceRunState *runState, 
                                                             LALInferenceVariables *proposedParams) {
  UINT4 nIFO, iIFO;
  LALInferenceIFOData *ifoData;
  REAL8 fPlus, fCross;
  LIGOTimeGPS timeGPS;
  REAL8 ra, dec, psi, t, gmst;
  REAL8 iotaNew, cosIotaNew, iota, cosIota;
  REAL8 dNew, d;
  REAL8 norm;

  LALInferenceCopyVariables(runState->currentParams, proposedParams);

  ra = *((REAL8 *)LALInferenceGetVariable(proposedParams, "rightascension"));
  dec = *((REAL8 *)LALInferenceGetVariable(proposedParams, "declination"));
  psi = *((REAL8 *)LALInferenceGetVariable(proposedParams, "polarisation"));
  t = *((REAL8 *)LALInferenceGetVariable(proposedParams, "time"));

  XLALGPSSetREAL8(&timeGPS, t);
  gmst = XLALGreenwichMeanSiderealTime(&timeGPS);

  /* Find number of IFO's. */
  nIFO=0;
  ifoData=runState->data;
  do {
    nIFO++;
    ifoData=ifoData->next;
  } while (ifoData != NULL);

  /* Now choose one at random, and get its data. */
  iIFO=gsl_rng_uniform_int(runState->GSLrandom, nIFO);
  ifoData=runState->data;
  while (iIFO > 0) {
    ifoData=ifoData->next;
    iIFO--;
  }

  /* Now we know which IFO we want to keep the magnitude constant in.
     Now compute f+, fx for that IFO. */
  XLALComputeDetAMResponse(&fPlus, &fCross, ifoData->detector->response, ra, dec, psi, gmst);

  iotaNew = M_PI*gsl_rng_uniform(runState->GSLrandom);
  cosIotaNew = cos(iotaNew);

  d = *((REAL8 *)LALInferenceGetVariable(proposedParams, "distance"));

  iota = *((REAL8 *)LALInferenceGetVariable(proposedParams, "inclination"));
  cosIota = cos(iota);

  norm = fabs((fPlus*(0.5*(1.0 + cosIota*cosIota)) + fCross*cosIota)/d);

  dNew = fabs((fPlus*(0.5*(1.0 + cosIotaNew*cosIotaNew)) + fCross*cosIotaNew) / norm);

  LALInferenceSetVariable(proposedParams, "distance", &dNew);
  if(LALInferenceCheckVariable(proposedParams,"logdistance")){
    REAL8 logdist=log(dNew);
    LALInferenceSetVariable(proposedParams,"logdistance",&logdist);
  }
  LALInferenceSetVariable(proposedParams, "inclination", &iotaNew);

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
  
}

void PTMCMCLALInferenceDifferentialEvolutionFull(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, proposedParams, NULL);
}

void PTMCMCLALInferenceDifferentialEvolutionNames(LALInferenceRunState *runState, 
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

  for (i = 0; names[i] != NULL; i++) {
    if (!LALInferenceCheckVariable(proposedParams, names[i]) || !LALInferenceCheckVariable(ptJ, names[i]) || !LALInferenceCheckVariable(ptI, names[i])) {
      /* Ignore variable if it's not in each of the params. */
    } else {
      REAL8 x = *((REAL8 *)LALInferenceGetVariable(proposedParams, names[i]));
      REAL8 scale = 1.66511*gsl_ran_ugaussian(runState->GSLrandom); 
      /* 1.66511 = number of sigma where Gaussian PDF drops to
         0.25. */
      
      x += scale * (*((REAL8 *) LALInferenceGetVariable(ptJ, names[i])));
      x -= scale * (*((REAL8 *) LALInferenceGetVariable(ptI, names[i])));
      
      LALInferenceSetVariable(proposedParams, names[i], &x);
    }
  }

  LALInferenceCyclicReflectiveBound(proposedParams, runState->priorArgs);
}

void PTMCMCLALInferenceDifferentialEvolutionMasses(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"chirpmass", "massratio", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionAmp(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"rightascension", "declination", "polarisation", "inclination", "distance", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionSpins(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"a_spin1", "a_spin2", "phi_spin1", "phi_spin2", "theta_spin1", "theta_spin2", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

void PTMCMCLALInferenceDifferentialEvolutionSky(LALInferenceRunState *runState, LALInferenceVariables *pp) {
  const char *names[] = {"rightascension", "declination", NULL};
  PTMCMCLALInferenceDifferentialEvolutionNames(runState, pp, names);
}

/*draws a value from the prior, uniformly in individual parameters used for jumps.*/
void PTMCMCLALInferenceDrawFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {

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
void PTMCMCLALInferenceDrawUniformlyFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams) {
  
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
}


