/* 
 *  InferenceNest.c:  Nested Sampling using LALInference
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
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferenceNestedSampler.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceReadBurstData.h>
#include <lal/LALInferenceCalibrationErrors.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>

/*************** MAIN **********************/


int main(int argc, char *argv[]){
  char help[]="\
  LALInferenceNest:\n\
  Bayesian analysis tool using Nested Sampling algorithm\n\
  for Burst analysis. Uses LALInference library for back-end.\n\n\
  Arguments for each section follow:\n\n";
  
  LALInferenceRunState *state;
  ProcessParamsTable *procParams=NULL;
  
  /* Read command line and parse */
  procParams=LALInferenceParseCommandLine(argc,argv);
  if(LALInferenceGetProcParamVal(procParams,"--help"))
  {
    fprintf(stdout,"%s",help);
  }
  
  /* initialise runstate based on command line */
  /* This includes reading in the data */
  /* And performing any injections specified */
  /* And allocating memory */
  state = LALInferenceInitRunState(procParams);
  
  ProcessParamsTable *ppt=NULL;
  if ((ppt=LALInferenceGetProcParamVal(state->commandLine,"--binj"))){
    /* Perform injections if data successful read or created */
    LALInferenceInjectBurstSignal(state->data, state->commandLine);
  }
  else{
    /* Perform CBC injection if required */
    LALInferenceInjectInspiralSignal(state->data, state->commandLine);
  }
  if (LALInferenceGetProcParamVal(state->commandLine,"--inject_from_mdc")){
      fprintf(stdout,"WARNING: Injecting a signal from MDC has not been carefully tested yet! \n"); 
      LALInferenceInjectFromMDC(state->commandLine, state->data);
  }

  /* Simulate calibration errors. 
  * NOTE: this must be called after both ReadData and (if relevant) 
  * injectInspiralTD/FD are called! */
  LALInferenceApplyCalibrationErrors(state->data, state->commandLine);

  /* Set up the appropriate functions for the nested sampling algorithm */
  if (state){
    /* Set up the appropriate functions for the nested sampling algorithm */
    state->algorithm=&LALInferenceNestedSamplingAlgorithm;
    state->evolve=&LALInferenceNestedSamplingOneStep;
    INT4 one=1;
    LALInferenceAddVariable(state->algorithmParams,"LIB",&one, LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    state->proposalArgs = LALInferenceParseProposalArgs(state);
  }

  /* Check if recovery is LIB or CBC */
  if ((ppt=LALInferenceGetProcParamVal(state->commandLine,"--approx"))){
    if (XLALCheckBurstApproximantFromString(ppt->value)){
      /* Set up the threads */
      LALInferenceInitBurstThreads(state,1);
      /* Init the prior */
      LALInferenceInitLIBPrior(state);
    }
    else{
      /* Set up the threads */
      LALInferenceInitCBCThreads(state,1);
      /* Init the prior */
      LALInferenceInitCBCPrior(state);
    }
  }
  else{
    fprintf(stderr,"Must specify the approximant while using lalinference_burst\n");
    exit(1);
    }
  /* Set up structures for nested sampling */
  LALInferenceNestedSamplingAlgorithmInit(state);
  
  for(INT4 i=0;i<state->nthreads;i++)
  {
    state->threads[i]->cycle=LALInferenceSetupDefaultInspiralProposalCycle(state->threads[i]->proposalArgs);
    LALInferenceRandomizeProposalCycle(state->threads[i]->cycle,state->GSLrandom);
  }

  /* Choose the likelihood and set some auxiliary variables */
  LALInferenceInitLikelihood(state);
  
  /* Exit since we printed all command line arguments */
  if(state == NULL || LALInferenceGetProcParamVal(state->commandLine,"--help"))
  {
    exit(0);
  }
  
  /* Call setupLivePointsArray() to populate live points structures */
  LALInferenceSetupLivePointsArray(state);
  
  /* write injection with noise evidence information from algorithm */
  // SALVO FIXME
  //LALInferencePrintInjectionSample(state);
  
  /* Call nested sampling algorithm */
  state->algorithm(state);
  
  /* end */
  return(0);
}


