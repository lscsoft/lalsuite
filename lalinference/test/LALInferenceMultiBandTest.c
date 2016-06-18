/*
 *  LALInferenceBench.c:  Benchmark LALInference functions
 *
 *  Copyright (C) 2015 John Veitch
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
#include <lal/LALInference.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceCalibrationErrors.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceLikelihood.h>
#include <sys/resource.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

const char HELPSTR[]=\
"LALInferenceMultiBandTest: Unit test for consistency between multiband and regular template functions.\n\
 Example (for 1.4-1.4 binary with seglen 32, srate 4096): \n\
 $ ./LALInferenceMultiBandTest --psdlength 1000 --psdstart 1 --seglen 32 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --fix-chirpmass 1.218 --fix-q 1.0 --margphi\n\n\n\
";

void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model);
void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model)
{
  return;
}

int compare_template(LALInferenceRunState *runState);
int compare_template(LALInferenceRunState *runState)
{
  REAL8 logLnormal=0.0;
  REAL8 logLmultiband=0.0;
  REAL8 tolerance = 0.01; /* Maximum allowable difference between standard and multiband */

  REAL8 oldSNR=0.0,mbSNR=0.0;
  REAL8 oldPhase=0.0,mbPhase=0.0;
  REAL8 CrossMatchPlus,CrossMatchCross;
  REAL8 h0MatchPlus,h0MatchCross;
  REAL8 hmbMatchPlus,hmbMatchCross;


  COMPLEX16FrequencySeries *oldTemplatePlus=NULL,*oldTemplateCross=NULL,*newTemplatePlus=NULL,*newTemplateCross=NULL;

  
  runState->threads[0]->model->templt=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
  LALInferenceTemplateNullFreqdomain(runState->threads[0]->model);
  /* FIXME: Have to call the function once before result is repeatable!!! */
  logLnormal = runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  oldTemplatePlus = runState->threads[0]->model->freqhPlus;
  oldTemplateCross = runState->threads[0]->model->freqhCross;

  runState->threads[0]->model->freqhPlus=XLALCreateCOMPLEX16FrequencySeries("mbtemplate",&oldTemplatePlus->epoch,oldTemplatePlus->f0,oldTemplatePlus->deltaF,&lalDimensionlessUnit,oldTemplatePlus->data->length);
  runState->threads[0]->model->freqhCross=XLALCreateCOMPLEX16FrequencySeries("mbtemplate",&oldTemplateCross->epoch,oldTemplateCross->f0,oldTemplateCross->deltaF,&lalDimensionlessUnit,oldTemplateCross->data->length);
  
  /* Clear the template */
  LALInferenceTemplateNullFreqdomain(runState->threads[0]->model);

  logLnormal = runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  oldSNR=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"optimal_snr");
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    oldPhase=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"phase_maxl");
  
  runState->threads[0]->model->templt=&LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated;
  
  logLmultiband = runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  int result = tolerance > fabs(logLmultiband-logLnormal)?0:1;
  mbSNR=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"optimal_snr");
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    mbPhase=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"phase_maxl");
  
  newTemplatePlus = runState->threads[0]->model->freqhPlus;
  newTemplateCross = runState->threads[0]->model->freqhCross;

  CrossMatchPlus = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplatePlus->data,oldTemplatePlus->data);
  CrossMatchCross = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplateCross->data,oldTemplateCross->data);

  h0MatchPlus = LALInferenceComputeFrequencyDomainOverlap(runState->data,oldTemplatePlus->data,oldTemplatePlus->data);
  h0MatchCross = LALInferenceComputeFrequencyDomainOverlap(runState->data,oldTemplateCross->data,oldTemplateCross->data);

  hmbMatchPlus = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplatePlus->data,newTemplatePlus->data);
  hmbMatchCross = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplateCross->data,newTemplateCross->data);
  /* Want to check that <h0|h0>+<h_mb|h_mb>-2<h0|h_mb> < tolerance */
  
  fprintf(stdout,"Parameter values:\n");
  LALInferencePrintVariables(runState->threads[0]->model->params);
  
  fprintf(stdout,"\n\n");
  fprintf(stdout,"matchPlus = %lf, matchCross = %lf\n",CrossMatchPlus,CrossMatchCross);
  fprintf(stdout,"h0Plus = %lf, h0Cross = %lf\n",h0MatchPlus,h0MatchCross);
  fprintf(stdout,"hmbPlus = %lf, hmbCross = %lf\n",hmbMatchPlus,hmbMatchCross);
  fprintf(stdout,"Difference <h0|h0> + <hmb|hmb> - 2<h0|hmb> = %lf\n",hmbMatchPlus+h0MatchPlus-2.0*CrossMatchPlus);

  fprintf(stdout,"Optimal SNR:\tnormal = %lf, phaseinterp = %lf\n",oldSNR,mbSNR);
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    fprintf(stdout,"max Like phase:\tnormal = %lf, phaseinterp = %lf\n",oldPhase,mbPhase);
  fprintf(stdout,"logL:\tnormal = %lf, logL multiband = %lf. Test result: %i\n",logLnormal,logLmultiband,result);
  fprintf(stdout,"Fractional difference epsilon = %lf\n",(logLnormal-logLmultiband)/logLnormal);
  return(result);
}

int main(int argc, char *argv[]){
  
  //#define default_command_line_len 18
  //char * const default_command_line[default_command_line_len];
  
  //"--psdlength","1000","--psdstart","1","--seglen","64","--srate","4096","--trigtime","0","--ifo","H1","--H1-channel","LALSimAdLIGO","--H1-cache","LALSimAdLIGO","--dataseed","1234"};

  
  ProcessParamsTable *procParams = NULL;
  LALInferenceRunState *runState=NULL;
  
  procParams=LALInferenceParseCommandLine(argc,argv);

  runState = LALInferenceInitRunState(procParams);
  
  if(runState) {
    LALInferenceInjectInspiralSignal(runState->data,runState->commandLine);

    /* Simulate calibration errors */
    LALInferenceApplyCalibrationErrors(runState->data,runState->commandLine);
  }
  
  /* Set up the template and likelihood functions */
  LALInferenceInitCBCThreads(runState,1);
  LALInferenceInitLikelihood(runState);

  /* Disable waveform caching */
  runState->threads[0]->model->waveformCache=NULL;
  
  int result = compare_template(runState);
  
  return(result);
}
