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
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <sys/resource.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

const char HELPSTR[]=\
"lalinference_bench: Benchmark template and likelihood functions.\n\
 Options:\n\
    --Niter            : Number of calls to time (delfault 1000) \n\
    --bench-template   : Only benchmark template function\n\
    --bench-likelihood : Only benchmark likelihood function\n\
                         (defaults to benchmarking both)\n\
 Example (for 1.0-1.0 binary with seglen 8, srate 4096): \n\
 $ ./lalinference_bench --psdlength 1000 --psdstart 1 --seglen 8 --srate 4096 --trigtime 0 --ifo H1 --H1-channel LALSimAdLIGO --H1-cache LALSimAdLIGO --dataseed 1324 --Niter 10000 --fix-chirpmass 1.218 --fix-q 1.0\n\n\n\
";

void fprintf_bench(FILE *fp, struct rusage start, struct rusage end, UINT4 Niter);
void fprintf_bench(FILE *fp, struct rusage start, struct rusage end, UINT4 Niter)
{
  REAL8 utime = (end.ru_utime.tv_sec - start.ru_utime.tv_sec) + 1e-6 * (end.ru_utime.tv_usec - start.ru_utime.tv_usec);
  REAL8 stime = (end.ru_stime.tv_sec - start.ru_stime.tv_sec) + 1e-6 * (end.ru_stime.tv_usec - start.ru_stime.tv_usec);
  
  fprintf(fp,"USER Total: %lf s\n",utime);
  fprintf(fp,"USER Per iteration: %e s\n",utime / (double) Niter);
  
  fprintf(fp,"SYS Total: %lf s\n",stime);
  fprintf(fp,"SYS Per iteration: %e s\n",stime / (double) Niter);
}

void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model);
void LALInferenceTemplateNoop(UNUSED LALInferenceModel *model)
{
  return;
}

void bench_likelihood(LALInferenceRunState *runState,UINT4 Niter);
void bench_likelihood(LALInferenceRunState *runState,UINT4 Niter)
{
  UINT4 i=0;
  struct rusage r_usage_start,r_usage_end;
  
  /* Clear the template */
  LALInferenceTemplateNullFreqdomain(runState->threads[0]->model);
  
  LALInferenceTemplateFunction old_templt=runState->threads[0]->model->templt;
  runState->threads[0]->model->templt=LALInferenceTemplateNoop;
  
  fprintf(stdout,"Benchmarking likelihood:\n");
  getrusage(RUSAGE_SELF, &r_usage_start);
  for(i=0;i<Niter;i++)
  {
    runState->likelihood(runState->threads[0]->model->params,runState->data, runState->threads[0]->model);
  }
  getrusage(RUSAGE_SELF, &r_usage_end);
  fprintf_bench(stdout, r_usage_start, r_usage_end, Niter);
  runState->threads[0]->model->templt=old_templt;
  
}

void bench_template(LALInferenceRunState *runState, UINT4 Niter);
void bench_template(LALInferenceRunState *runState, UINT4 Niter)
{
  UINT4 i=0;
  struct rusage r_usage_start,r_usage_end;
  fprintf(stdout,"Benchmarking template:\n");
  getrusage(RUSAGE_SELF, &r_usage_start);
  for(i=0;i<Niter;i++)
  {
    runState->threads[0]->model->templt(runState->threads[0]->model);
  }
  getrusage(RUSAGE_SELF, &r_usage_end);
  fprintf_bench(stdout, r_usage_start, r_usage_end, Niter);
}

int main(int argc, char *argv[]){
  ProcessParamsTable *procParams = NULL,*ppt=NULL;
  LALInferenceRunState *runState=NULL;
  UINT4 Niter=1000;
  UINT4 bench_L=1;
  UINT4 bench_T=1;
  
  procParams=LALInferenceParseCommandLine(argc,argv);

  if(LALInferenceGetProcParamVal(procParams,"--help"))
  {
    fprintf(stdout,"%s",HELPSTR);
    bench_T=bench_L=0;
  }
  if((ppt=LALInferenceGetProcParamVal(procParams,"--Niter")))
     Niter=atoi(ppt->value);
  if(LALInferenceGetProcParamVal(procParams,"--bench-template"))
  {
    bench_T=1; bench_L=0;
  }
  if(LALInferenceGetProcParamVal(procParams,"--bench-likelihood"))
  {
    bench_T=0; bench_L=1;
  }

  
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
  
  if(bench_T)
  {
    printf("Template test will run with parameters:\n");
    LALInferencePrintVariables(runState->threads[0]->model->params);
    printf("\n");

    bench_template(runState,Niter);
    printf("\n");
  }
  if(bench_L)
  {
    bench_likelihood(runState,Niter);
    printf("\n");
  }
  
  return(0);
}
