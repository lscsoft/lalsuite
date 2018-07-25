/*
 *  *  LALInferenceInjectedLike.c: Util bin to create an *injection file with the true parameters and the injected logL/P
 *   *
 *    *  Copyright (C) 2018 salvatore vitale 
 *     *
 *      *
 *       *  This program is free software; you can redistribute it and/or modify
 *        *  it under the terms of the GNU General Public License as published by
 *         *  the Free Software Foundation; either version 2 of the License, or
 *          *  (at your option) any later version.
 *           *
 *            *  This program is distributed in the hope that it will be useful,
 *             *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *               *  GNU General Public License for more details.
 *                *
 *                 *  You should have received a copy of the GNU General Public License
 *                  *  along with with program; see the file COPYING. If not, write to the
 *                   *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *                    *  MA  02111-1307  USA
 *                     */


#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceCalibrationErrors.h>
#include <lal/LALInferenceVCSInfo.h>

/*************** MAIN **********************/


int main(int argc, char *argv[]){
  int helpflag=0;
  char help[]="\
  LALInferenceInjectedLike:\n\
  Print the injected values of parameters and (delta)LogL/P \n";

  LALInferenceRunState *state;
  ProcessParamsTable *procParams=NULL;
  LALInferenceIFOData *data = NULL;

  /* Read command line and parse */
  procParams=LALInferenceParseCommandLine(argc,argv);
  if(LALInferenceGetProcParamVal(procParams,"--help"))
  {
    helpflag=1;
    fprintf(stdout,"%s",help);
  }
  /* write down git information */
  fprintf(stdout,"\n\nLALInference version:%s,%s,%s,%s,%s\n\n", lalInferenceVCSInfo.vcsId,lalInferenceVCSInfo.vcsDate,lalInferenceVCSInfo.vcsBranch,lalInferenceVCSInfo.vcsAuthor,lalInferenceVCSInfo.vcsStatus);
  
  /* initialise runstate based on command line */
  /* This includes reading in the data */
  /* And performing any injections specified */
  /* And allocating memory */
  state = LALInferenceInitRunState(procParams);
  /* Create header  */
  if (state!=NULL && !helpflag){
    ProcessParamsTable *ppt=NULL;
    ppt=LALInferenceGetProcParamVal(state->commandLine,"--outfile");
    if(!ppt){
    ppt=LALInferenceGetProcParamVal(state->commandLine,"--outhdf");
    if(!ppt){
      fprintf(stderr,"Must specify --outfile <filename.dat> or --outhdf <filename.h5>\n");
      exit(1);
      }
  }
    
    }
  if (state == NULL) {
      if (!helpflag) {
          fprintf(stderr, "run state not allocated (%s, line %d).\n",
                  __FILE__, __LINE__);
      }
  } else {
      data = state->data;
  }

  /* Perform injections if data successful read or created */
  if (state&&!helpflag){
    LALInferenceInjectInspiralSignal(data, state->commandLine);
  }

  /* Simulate calibration errors. 
 *   * NOTE: this must be called after both ReadData and (if relevant) 
 *     * injectInspiralTD/FD are called! */
  LALInferenceApplyCalibrationErrors(data, procParams);


  if (!helpflag && LALInferenceGetProcParamVal(state->commandLine, "--roqtime_steps")){

        LALInferenceSetupROQdata(state->data, state->commandLine);
        fprintf(stderr, "done LALInferenceSetupROQdata\n");

     }

  /* Set up the threads */
  LALInferenceInitCBCThreads(state,1);

  /* Init the prior */
  LALInferenceInitCBCPrior(state);

  /* Choose the likelihood and set some auxiliary variables */
  LALInferenceInitLikelihood(state);

  /* Exit since we printed all command line arguments */
  if(state == NULL || LALInferenceGetProcParamVal(state->commandLine,"--help"))
  {
    exit(0);
  }

  /* write injection with noise evidence information from algorithm */
  LALInferencePrintInjectionSample(state);

  /* end */
  return(0);
}

