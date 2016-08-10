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

COMPLEX16 compute_mismatch(LALInferenceIFOData *data, COMPLEX16FrequencySeries *a, COMPLEX16FrequencySeries *b);


void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model);
void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model)
{
  return;
}

int compare_template(LALInferenceRunState *runState);
int compare_template(LALInferenceRunState *runState)
{
  REAL8 oldPhase=0.0,mbPhase=0.0;
  REAL8 h0MatchPlus,h0MatchCross;
  REAL8 hmbMatchPlus,hmbMatchCross;

  REAL8 target_snr = 50; /* Max SNR that we expect to handle */
  REAL8 tolerance = 0.1/(target_snr*target_snr); /* Error in likelihood of 0.1 at target SNR */

  COMPLEX16FrequencySeries *oldTemplatePlus=NULL,*oldTemplateCross=NULL,*newTemplatePlus=NULL,*newTemplateCross=NULL;

  
  runState->threads[0]->model->templt=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
  LALInferenceTemplateNullFreqdomain(runState->threads[0]->model);
  /* FIXME: Have to call the function once before result is repeatable!!! */
  runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  oldTemplatePlus = runState->threads[0]->model->freqhPlus;
  oldTemplateCross = runState->threads[0]->model->freqhCross;

  runState->threads[0]->model->freqhPlus=XLALCreateCOMPLEX16FrequencySeries("mbtemplate",&oldTemplatePlus->epoch,oldTemplatePlus->f0,oldTemplatePlus->deltaF,&lalDimensionlessUnit,oldTemplatePlus->data->length);
  runState->threads[0]->model->freqhCross=XLALCreateCOMPLEX16FrequencySeries("mbtemplate",&oldTemplateCross->epoch,oldTemplateCross->f0,oldTemplateCross->deltaF,&lalDimensionlessUnit,oldTemplateCross->data->length);
  
  /* Clear the template */
  LALInferenceTemplateNullFreqdomain(runState->threads[0]->model);

  runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  LALInferenceDumpWaveforms(runState->threads[0]->model, "normal");
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    oldPhase=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"phase_maxl");
  
  runState->threads[0]->model->templt=&LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated;
  
  runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  LALInferenceDumpWaveforms(runState->threads[0]->model, "multiband");
  
  REAL8 SNRsqPlus =LALInferenceComputeFrequencyDomainOverlap(runState->data,runState->threads[0]->model->freqhPlus->data,runState->threads[0]->model->freqhPlus->data);
  REAL8 SNRsqCross =LALInferenceComputeFrequencyDomainOverlap(runState->data,runState->threads[0]->model->freqhCross->data,runState->threads[0]->model->freqhCross->data);
  
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    mbPhase=LALInferenceGetREAL8Variable(runState->threads[0]->model->params,"phase_maxl");
  
  newTemplatePlus = runState->threads[0]->model->freqhPlus;
  newTemplateCross = runState->threads[0]->model->freqhCross;

  h0MatchPlus = LALInferenceComputeFrequencyDomainOverlap(runState->data,oldTemplatePlus->data,oldTemplatePlus->data);
  h0MatchCross = LALInferenceComputeFrequencyDomainOverlap(runState->data,oldTemplateCross->data,oldTemplateCross->data);

  hmbMatchPlus = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplatePlus->data,newTemplatePlus->data);
  hmbMatchCross = LALInferenceComputeFrequencyDomainOverlap(runState->data,newTemplateCross->data,newTemplateCross->data);
  /* Want to check that <h0|h0>+<h_mb|h_mb>-2<h0|h_mb> < tolerance */
  
  fprintf(stdout,"Parameter values:\n");
  LALInferencePrintVariables(runState->threads[0]->model->params);
  
  COMPLEX16 mismatchplus = compute_mismatch(runState->data, newTemplatePlus , oldTemplatePlus);
  COMPLEX16 mismatchcross = compute_mismatch(runState->data, newTemplateCross , oldTemplateCross);
  
  COMPLEX16 innerPlus = LALInferenceComputeFrequencyDomainComplexOverlap(runState->data, newTemplatePlus->data, oldTemplatePlus->data);
  COMPLEX16 innerCross = LALInferenceComputeFrequencyDomainComplexOverlap(runState->data, newTemplateCross->data, oldTemplateCross->data);
  
  int result = (1.0-cabs(innerPlus)/h0MatchPlus < tolerance) && (1.0 - cabs(innerCross)/h0MatchCross < tolerance);
  
  fprintf(stdout,"\n\n");
  fprintf(stdout,"|<h0|hmb>/<h0|h0>| = %lf (+) %lf (x)\n",(cabs(innerPlus)/h0MatchPlus),(cabs(innerCross)/h0MatchCross));
  fprintf(stdout,"SNR    = %lf (+), %lf (x)\n",sqrt(h0MatchPlus),sqrt(h0MatchCross));
  fprintf(stdout,"SNR mb = %lf (+), %lf (x)\n",sqrt(hmbMatchPlus),sqrt(hmbMatchCross));
  fprintf(stdout,"cplx <h|h0> (+): %lf*exp(%lfi), (x): %lf*exp(%lfi)\n",cabs(innerPlus),carg(innerPlus),cabs(innerCross),carg(innerCross));

  fprintf(stdout,"Normalised mismatch (|h0-h|/|h0|)^2: %le (+), %le (x)\n",cabs(mismatchplus)/SNRsqPlus,cabs(mismatchcross)/SNRsqCross);
  fprintf(stdout,"Tolerance = %le\n",tolerance);
  fprintf(stdout,"Test result: %s\n",result?"passed":"failed");
  fprintf(stdout,"Interpolation good up to SNR %lf (+), %lf (x)\n",
	  sqrt(tolerance*target_snr*target_snr *SNRsqPlus/ mismatchplus),
	  sqrt(tolerance*target_snr*target_snr *SNRsqCross/ mismatchplus) );
  if(LALInferenceCheckVariable(runState->threads[0]->model->params,"phase_maxl"))
    fprintf(stdout,"max Like phase:\tnormal = %lf, phaseinterp = %lf\n",oldPhase,mbPhase);
  return(result);
}

/* Computes <a-b|a-b> */
COMPLEX16 compute_mismatch(LALInferenceIFOData *data, COMPLEX16FrequencySeries *a, COMPLEX16FrequencySeries *b)
{
  REAL8 aa = LALInferenceComputeFrequencyDomainComplexOverlap(data,a->data,a->data);
  REAL8 bb = LALInferenceComputeFrequencyDomainComplexOverlap(data,b->data,b->data);
  REAL8 ab = LALInferenceComputeFrequencyDomainComplexOverlap(data,a->data,b->data);
  REAL8 ba = LALInferenceComputeFrequencyDomainComplexOverlap(data,b->data,a->data);
  return aa+bb-ab-ba;
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
