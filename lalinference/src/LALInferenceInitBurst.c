/*
 *  LALInferenceCBCInit.c:  Bayesian Followup initialisation routines.
 *
 *  Copyright (C) 2013 James Clark, John Veitch, Salvatore Vitale
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
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceReadBurstData.h>
#include <lal/LALInferenceInit.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/GenerateBurst.h>
#include <lal/LALSimBurst.h>

LALInferenceTemplateFunction LALInferenceInitBurstTemplate(LALInferenceRunState *runState)
{

  char help[]="(--approx [SineGaussian,SineGaussianF,Gaussian,GaussianF,RingdownF]\tSpecify approximant to use (default SineGaussianF)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  LALInferenceTemplateFunction templt = &LALInferenceTemplateXLALSimInspiralChooseWaveform;
  
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt) {
    if(!strcmp("PrincipalComp",ppt->value)){
        printf("Using LALInferenceTemplatePrincipalComp \n");
        templt=&LALInferenceTemplatePrincipalComp;
    }
    else if(!strcmp("PrincipalCompSpec",ppt->value)){
        printf("Using LALInferenceTemplatePrincipalCompSpec \n");
        templt=&LALInferenceTemplatePrincipalCompSpec;
    }
    else if(XLALSimBurstImplementedFDApproximants(XLALGetBurstApproximantFromString(ppt->value))){
        templt=&LALInferenceTemplateXLALSimBurstChooseWaveform;
    }
    else if(XLALSimBurstImplementedTDApproximants(XLALGetBurstApproximantFromString(ppt->value))){
        templt=&LALInferenceTemplateXLALSimBurstChooseWaveform;
    }
    else {
      XLALPrintError("Error: unknown template %s\n",ppt->value);
      XLALPrintError(help);
      //XLAL_ERROR_VOID(XLAL_EINVAL);
    }
  }
  return templt;
}

/* Setup the variables to control Burst template generation */
/* Includes specification of prior ranges */

LALInferenceModel * LALInferenceInitBurstModel(LALInferenceRunState *state)
{
  fprintf(stderr,"Using LALInferenceBurstVariables!\n");

  LALStatus status;
  SimBurst *BinjTable=NULL;
  SimInspiralTable *inj_table=NULL;
	state->currentParams=XLALCalloc(1,sizeof(LALInferenceVariables));
	ProcessParamsTable *commandLine=state->commandLine;
	REAL8 endtime=-1;
	REAL8 endtime_from_inj=-1;
  ProcessParamsTable *ppt=NULL;
	memset(&status,0,sizeof(LALStatus));
	INT4 event=0;	
	INT4 i=0;
  BurstApproximant approx = (BurstApproximant) 0;
  char *pinned_params=NULL;
	char help[]="\
Parameter arguments:\n\
(--inj injections.xml)\tSimInspiral or SimBurst Injection XML file to use\n\
(--dt time)\tWidth of time prior, centred around trigger (0.1s)\n\
(--trigtime time)\tTrigger time to use\n\
(--approx Approximant)\tSet approximant (SineGaussianF,SineGaussian,Gaussian,RingdownF)\n\
(--fref fRef)\tSpecify a reference frequency at which parameters are defined (default 0).\n\
(--pinparams [frequency,q,loghrss, etc])\n\tList of parameters to set to injected values\n";
	/* Print command line arguments if help requested */
	ppt=LALInferenceGetProcParamVal(commandLine,"--help");
	if(ppt)
	{
		fprintf(stdout,"%s",help);
		return 0;
	}
  
  LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  memset(model->params, 0, sizeof(LALInferenceVariables));
  LALInferenceVariables *currentParams=model->params;
  
  //state->proposal=&NSWrapMCMCSineGaussProposal;
  
  /* We may have used a CBC injection... test */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if (ppt)
      endtime=atof(ppt->value);
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if (ppt) {
    BinjTable=XLALSimBurstTableFromLIGOLw(LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
    if (BinjTable){
      //burst_inj=1;
      ppt=LALInferenceGetProcParamVal(commandLine,"--event");
      if(ppt){
        event = atoi(ppt->value);
        while(i<event) {i++; BinjTable = BinjTable->next;}
      }
      else
        fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
      endtime_from_inj=XLALGPSGetREAL8(&(BinjTable->time_geocent_gps));
    }
    else{
      SimInspiralTableFromLIGOLw(&inj_table,LALInferenceGetProcParamVal(commandLine,"--inj")->value,0,0);
      if (inj_table){
        ppt=LALInferenceGetProcParamVal(commandLine,"--event");
        if(ppt){
          event= atoi(ppt->value);
          fprintf(stderr,"Reading event %d from file\n",event);
          i =0;
          while(i<event) {i++; inj_table=inj_table->next;} /* select event */
          endtime_from_inj=XLALGPSGetREAL8(&(inj_table->geocent_end_time));
        }
        else
            fprintf(stdout,"WARNING: You did not provide an event number with you --inj. Using default event=0 which may not be what you want!!!!\n");
      }
    }
  }
  if(!(BinjTable || inj_table || endtime>=0.0 )){
    printf("Did not provide --trigtime or an xml file and event... Exiting.\n");
    exit(1);
  }
  if (endtime_from_inj!=endtime){
    if(endtime_from_inj>0 && endtime>0)
      fprintf(stderr,"WARNING!!! You set trigtime %lf with --trigtime but event %i seems to trigger at time %lf\n",endtime,event,endtime_from_inj);
      else if(endtime_from_inj>0)
    endtime=endtime_from_inj;
  }
        
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    char **strings=NULL;
    UINT4 N;
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    LALInferenceBurstInjectionToVariables(BinjTable,&tempParams);

    LALInferenceVariableItem *node=NULL;
    while(N>0){
      N--;
      char *name=strings[N];
      node=LALInferenceGetItem(&tempParams,name);
      if(node) {LALInferenceAddVariable(currentParams,node->name,node->value,node->type,node->vary); printf("pinned %s \n",node->name);}
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
    }
  }
    
  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    approx = XLALGetBurstApproximantFromString(ppt->value);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt){
    approx = XLALGetBurstApproximantFromString(ppt->value);
    }
   /* Set the model domain appropriately */
  if (XLALSimBurstImplementedFDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  } else if (XLALSimBurstImplementedTDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_TIME;
  } else {
    fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
    exit(1);
  }
  
  LALInferenceAddVariable(currentParams, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
 
    REAL8 psiMin=0.0,psiMax=LAL_PI; 
    REAL8 raMin=0.0,raMax=LAL_TWOPI; 
    REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0; 
    REAL8 qMin=3., qMax=100.0;
    REAL8 ffMin=40., ffMax=1300.0;
    REAL8 durMin=1.0e-4; // min and max value of duration for gaussian templates 
    REAL8 durMax=.5;
    REAL8 hrssMin=1.e-23, hrssMax=1.0e-21;
    REAL8 loghrssMin=log(hrssMin),loghrssMax=log(hrssMax);
    REAL8 dt=0.1;
    REAL8 timeMin=endtime-0.5*dt; 
    REAL8 timeMax=endtime+0.5*dt;
    REAL8 phiMin=0.0,phiMax=LAL_TWOPI;

    REAL8 zero=0.0;

    ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
    if (ppt) dt=atof(ppt->value);
    
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", zero, timeMin, timeMax, LALINFERENCE_PARAM_LINEAR);
    
    /* If we are marginalising over the time, remove that variable from the model (having set the prior above) */
    /* Also set the prior in model->params, since Likelihood can't access the state! (ugly hack) */
    if(LALInferenceGetProcParamVal(commandLine,"--margtime") || LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
        LALInferenceVariableItem *p=LALInferenceGetItem(state->priorArgs,"time_min");
        LALInferenceAddVariable(model->params,"time_min",p->value,p->type,p->vary);
        p=LALInferenceGetItem(state->priorArgs,"time_max");
        LALInferenceAddVariable(model->params,"time_max",p->value,p->type,p->vary);
        LALInferenceRemoveVariable(model->params,"time");
        if (LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
            UINT4 margphi = 1;
            LALInferenceAddVariable(model->params, "margtimephi", &margphi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
        }
    }

    LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, raMin, raMax, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decMin, decMax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psiMin, psiMax, LALINFERENCE_PARAM_LINEAR);

    ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
    if (!strcmp("SineGaussian",ppt->value) || !strcmp("SineGaussianF",ppt->value)|| !strcmp("DampedSinusoid",ppt->value) || !strcmp("DampedSinusoidF",ppt->value)){
     
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "frequency",  zero, ffMin, ffMax,   LALINFERENCE_PARAM_LINEAR);
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "quality",  zero,qMin, qMax,   LALINFERENCE_PARAM_LINEAR);
      if(!LALInferenceGetProcParamVal(commandLine,"--margphi") && !LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phase", zero, phiMin, phiMax, LALINFERENCE_PARAM_CIRCULAR);
      }
    }
    else if (!strcmp("Gaussian",ppt->value) || !strcmp("GaussianF",ppt->value)){
      LALInferenceRegisterUniformVariableREAL8(state, model->params,"duration", zero, durMin,durMax, LALINFERENCE_PARAM_LINEAR);
    }
    
    if (LALInferenceGetProcParamVal(commandLine,"--use-hrss")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "hrss",  zero,hrssMin, hrssMax,   LALINFERENCE_PARAM_LINEAR);
    }
    else
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "loghrss",  zero,loghrssMin, loghrssMax,   LALINFERENCE_PARAM_LINEAR);
    
    LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", zero,0.0,2*LAL_PI ,LALINFERENCE_PARAM_CIRCULAR);
    ppt=LALInferenceGetProcParamVal(commandLine,"--cross_only");
    if (ppt){
      printf("Fixing alpha to Pi/2 in template ---> only cross polarization will be used\n");
      LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", LAL_PI/2.0,0.0,2*LAL_PI ,LALINFERENCE_PARAM_FIXED);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--plus_only");
    if (ppt){
        printf("Fixing alpha to 0 in template ---> only plus polarization will be used\n");
        LALInferenceRegisterUniformVariableREAL8(state,model->params, "alpha", 0.0,0.0,2*LAL_PI, LALINFERENCE_PARAM_FIXED);
    }
    
    /* Needs two condition: must be a burst template and the burst injection must have been provided to do those checks
    if (BinjTable && burst_inj){
        
        if (log(BinjTable->hrss) > loghrssMax || log(BinjTable->hrss) < loghrssMin)
            fprintf(stdout,"WARNING: The injected value of loghrss (%.4e) lies outside the prior range. That may be ok if your template and injection belong to different WF families\n",log(BinjTable->hrss));
        if (BinjTable->q > qMax || BinjTable->q < qMin )
            fprintf(stdout,"WARNING: The injected value of q (%lf) lies outside the prior range\n",BinjTable->q);
         if (BinjTable->frequency > ffMax || BinjTable->frequency < ffMin )
            fprintf(stdout,"WARNING: The injected value of centre_frequency (%lf) lies outside the prior range\n",BinjTable->frequency);
        
        ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
        if (!strcmp("Gaussian",ppt->value) || !strcmp("GaussianF",ppt->value)){
          if (BinjTable->duration >durMax || BinjTable->duration < durMin )
            fprintf(stdout,"WARNING: The injected value of centre_frequency (%lf) lies outside the prior range\n",BinjTable->frequency);
        }
        // Check the max Nyquist frequency for this parameter range
        if ( (ffmax+ 3.0*ffmax/qmin) > state->data->fHigh){
            fprintf(stderr,"WARNING, some of the template in your parameter space will be generated at a frequency higher than Nyquist (%lf). This is bound to produce unwanted consequences. Consider increasing the sampling rate, or reducing (increasing) the max (min) value of frequency (Q). With current setting, srate must be higher than %lf\n",state->data->fHigh,2*(ffmax+ 3.0*ffmax/qmin));
            //exit(1);
        }
            
    }*/
  /* Set model sampling rates to be consistent with data */
  model->deltaT = state->data->timeData->deltaT;
  model->deltaF = state->data->freqData->deltaF;
  UINT4 nifo=0;
  LALInferenceIFOData *dataPtr = state->data;
  while (dataPtr != NULL)
  {
    dataPtr = dataPtr->next;
    nifo++;
  }
  /* Initialize waveform buffers */
  model->timehPlus  = XLALCreateREAL8TimeSeries("timehPlus",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->timehCross = XLALCreateREAL8TimeSeries("timehCross",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->freqhPlus = XLALCreateCOMPLEX16FrequencySeries("freqhPlus",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);
  model->freqhCross = XLALCreateCOMPLEX16FrequencySeries("freqhCross",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  /* Choose proper template */
  model->templt = LALInferenceInitBurstTemplate(state);

  /* Use same window and FFT plans on model as data */
  model->window = state->data->window;
  model->padding = state->data->padding;
  model->timeToFreqFFTPlan = state->data->timeToFreqFFTPlan;
  model->freqToTimeFFTPlan = state->data->freqToTimeFFTPlan;
  /* Initialize waveform cache */
  model->burstWaveformCache = XLALCreateSimBurstWaveformCache();

  return model;
}

LALInferenceModel *LALInferenceInitModelReviewBurstEvidence_unimod(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.001, .min=-0.006121, .max=0.008121},
    {.name="frequency", .val=210., .min=205.346948, .max=216.653052},
    {.name="quality", .val=6.03626, .min=5.043829, .max=6.956171},
    {.name="loghrss", .val=-46., .min=-46.985195, .max=-45.014805},
    {.name="phase", .val=1.008, .min=0.718919, .max=1.281081},
    {.name="polarisation", .val=0.73, .min=0.427564, .max=0.972436},
    {.name="rightascension", .val=LAL_PI, .min=2.837864, .max=3.445321},
    {.name="declination", .val=0.04, .min= -0.334492, .max=0.334492},
    {.name="alpha", .val=0.58, .min=0.200742, .max=0.799258},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceRegisterUniformVariableREAL8(state, model->params, setup[i].name, setup[i].val, setup[i].min, setup[i].max, type);
    i++;
  }
  return(model);
}

LALInferenceModel *LALInferenceInitModelReviewBurstEvidence_bimod(LALInferenceRunState *state)
{
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i,j;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }
  LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  i=0;
  
  struct varSettings {const char *name; REAL8 val, min, max;};
  
  struct varSettings setup[]=
  {
    {.name="time", .val=0.001, .min=-0.006121, .max=0.019514},
    {.name="frequency", .val=211., .min=205.346948, .max=225.697936},
    {.name="quality", .val=6.0, .min=5.043829, .max=8.486044},
    {.name="loghrss", .val=-46., .min=-46.985195, .max=-43.438492},
    {.name="phase", .val=1.0, .min=0.718919, .max=1.730810},
    {.name="polarisation", .val=0.73, .min=0.427564,.max=1.408335},
    {.name="rightascension", .val=LAL_PI, .min=2.837864, .max=3.931287},
    {.name="declination", .val=0.0, .min=-0.334492, .max=0.869678},
    {.name="alpha", .val=0.5, .min=0.200742, .max=1.278070},
    {.name="END", .val=0., .min=0., .max=0.}
  };
  
  while(strcmp("END",setup[i].name))
  {
    LALInferenceParamVaryType type=LALINFERENCE_PARAM_CIRCULAR;
    /* Check if it is to be fixed */
    for(j=0;j<N;j++) if(!strcmp(setup[i].name,strings[j])) {type=LALINFERENCE_PARAM_FIXED; printf("Fixing parameter %s\n",setup[i].name); break;}
    LALInferenceRegisterUniformVariableREAL8(state, model->params, setup[i].name, setup[i].val, setup[i].min, setup[i].max, type);
    i++;
  }
  return(model);
}

/* Setup the variables to control Principal component template generation */
/* Includes specification of prior ranges */

LALInferenceModel * LALInferenceInitPrincipalCompModel(LALInferenceRunState *state)
{

    printf("-----Using LALInferenceInitPrincipalCompVariables!\n");


    LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
    model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
    memset(model->params, 0, sizeof(LALInferenceVariables));

    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;

    /* PC model configuration TODO: sanity check values */
    UINT4 nPCs = 0; // number of principle components to use 
    UINT4 nPCs_max = 20; // max number of PCs
    UINT4 ncatrows = 0; // number of rows in pc matrix (frequency samples) 
    UINT4 ncatcols = 0; // number of columns in pc matrix (waveforms) 
    char *PCfilePlus = NULL; // name of PC file
    char *PCfileCross = NULL; // name of PC file

    if((ppt=LALInferenceGetProcParamVal(commandLine,"--nPCs"))){
        nPCs=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --nPCs");
        exit(1);
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--ncatrows"))){
        ncatrows=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --ncatrows");
        exit(1);
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--ncatcols"))){
        ncatcols=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --ncatcols");
        exit(1);
    }

    /* Check PC matrix is compatible with data */
    if (state->data->freqData->data->length != ncatrows){
        XLALPrintError("length of F-domain data (%d) does not match length of F-domain PCs (%d)\n",
                state->data->freqData->data->length, ncatrows);
        exit(1);
    }

    /* Retrieve PCs */
    model->pcs = XLALMalloc(sizeof(LALInferencePCsModel));
    model->pcs->nPCs = nPCs;

    if((ppt=LALInferenceGetProcParamVal(commandLine,"--PCfile_pluspol"))){
        PCfilePlus=ppt->value;
        fprintf(stdout, "Loading PCs from + file\n");
        /* Read PC matrix */
        gsl_matrix_complex *tmp = get_complex_matrix_from_file(PCfilePlus, ncatrows, ncatcols);
        /* Now put first nPCs into a ncatrows x 10 matrix for easily
         * standardised template functions */
        model->pcs->pcs_plus = copy_npcs_from_complex_matrix(tmp, nPCs_max, nPCs, ncatrows, ncatcols);
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--PCfile_crosspol"))){
        PCfileCross=ppt->value;
        fprintf(stdout, "Loading PCs from x file\n");
        /* Read PC matrix */
        gsl_matrix_complex *tmp = get_complex_matrix_from_file(PCfileCross, ncatrows, ncatcols);
        /* Now put first nPCs into a ncatrows x 10 matrix for easily
         * standardised template functions */
        model->pcs->pcs_cross = copy_npcs_from_complex_matrix(tmp, nPCs_max, nPCs, ncatrows, ncatcols);
    }

    if (!LALInferenceGetProcParamVal(commandLine,"--PCfile_pluspol") &&
            !LALInferenceGetProcParamVal(commandLine,"--PCfile_crosspol")){
        XLALPrintError("must supply at least one of --PCfile_pluspol and --PCfile_crosspol\n");
        exit(1);
    }

    REAL8 hrssMin=1.e-26, hrssMax=1.0e-16;
    REAL8 loghrssMin=log(hrssMin),loghrssMax=log(hrssMax);
    REAL8 psimin=0.0,psimax=LAL_PI,psi=0.0;
    REAL8 ramin=0.0,ramax=LAL_TWOPI,ra=0.0;
    REAL8 decmin=-LAL_PI/2.0,decmax=LAL_PI/2.0,dec=0.0;
    REAL8 endtime=0.0;
    REAL8 dt=0.1;
    REAL8 zero=0.0;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime"))) endtime=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--dt"))) dt=atof(ppt->value);

    REAL8 timeMin=endtime-0.5*dt;
    REAL8 timeMax=endtime+0.5*dt;



    if (nPCs>=1){

        REAL8 beta1_min=-250, beta1_max=250;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta1_min"))) beta1_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta1_max"))) beta1_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta1", zero, beta1_min, beta1_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=2){

        REAL8 beta2_min=-150, beta2_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta2_min"))) beta2_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta2_max"))) beta2_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta2", zero, beta2_min, beta2_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=3){

        REAL8 beta3_min=-150, beta3_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta3_min"))) beta3_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta3_max"))) beta3_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta3", zero, beta3_min, beta3_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=4){

        REAL8 beta4_min=-100, beta4_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta4_min"))) beta4_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta4_max"))) beta4_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta4", zero, beta4_min, beta4_max, LALINFERENCE_PARAM_LINEAR);
    }
    if (nPCs>=5){

        REAL8 beta5_min=-100, beta5_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta5_min"))) beta5_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta5_max"))) beta5_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta5", zero, beta5_min, beta5_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=6){

        REAL8 beta6_min=-100, beta6_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta6_min"))) beta6_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta6_max"))) beta6_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta6", zero, beta6_min, beta6_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=7){

        REAL8 beta7_min=-100, beta7_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta7_min"))) beta7_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta7_max"))) beta7_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta7", zero, beta7_min, beta7_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=8){

        REAL8 beta8_min=-100, beta8_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta8_min"))) beta8_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta8_max"))) beta8_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta8", zero, beta8_min, beta8_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=9){

        REAL8 beta9_min=-100, beta9_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta9_min"))) beta9_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta9_max"))) beta9_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta9", zero, beta9_min, beta9_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=10){

        REAL8 beta10_min=-100, beta10_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta10_min"))) beta10_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta10_max"))) beta10_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta10", zero, beta10_min, beta10_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=11){

        REAL8 beta11_min=-100, beta11_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta11_min"))) beta11_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta11_max"))) beta11_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta11", zero, beta11_min, beta11_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=12){

        REAL8 beta12_min=-100, beta12_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta12_min"))) beta12_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta12_max"))) beta12_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta12", zero, beta12_min, beta12_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=13){

        REAL8 beta13_min=-100, beta13_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta13_min"))) beta13_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta13_max"))) beta13_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta13", zero, beta13_min, beta13_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=14){

        REAL8 beta14_min=-100, beta14_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta14_min"))) beta14_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta14_max"))) beta14_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta14", zero, beta14_min, beta14_max, LALINFERENCE_PARAM_LINEAR);
    }
    if (nPCs>=15){

        REAL8 beta15_min=-50, beta15_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta15_min"))) beta15_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta15_max"))) beta15_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta15", zero, beta15_min, beta15_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=16){

        REAL8 beta16_min=-50, beta16_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta16_min"))) beta16_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta16_max"))) beta16_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta16", zero, beta16_min, beta16_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=17){

        REAL8 beta17_min=-50, beta17_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta17_min"))) beta17_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta17_max"))) beta17_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta17", zero, beta17_min, beta17_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=18){

        REAL8 beta18_min=-50, beta18_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta18_min"))) beta18_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta18_max"))) beta18_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta18", zero, beta18_min, beta18_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=19){

        REAL8 beta19_min=-50, beta19_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta19_min"))) beta19_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta19_max"))) beta19_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta19", zero, beta19_min, beta19_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs==20){

        REAL8 beta20_min=-50, beta20_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta20_min"))) beta20_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta20_max"))) beta20_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta20", zero, beta20_min, beta20_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (LALInferenceGetProcParamVal(commandLine,"--use-hrss")){
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--hrssMin"))) hrssMin=(atof(ppt->value));
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--hrssMax"))) hrssMax=(atof(ppt->value));
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "hrss",  zero,hrssMin, hrssMax,   LALINFERENCE_PARAM_LINEAR);
    }
    else
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssMin"))) loghrssMin=(atof(ppt->value));
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssMax"))) loghrssMax=(atof(ppt->value));
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "loghrss",  zero,loghrssMin, loghrssMax,   LALINFERENCE_PARAM_LINEAR);

    LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", zero, timeMin, timeMax, LALINFERENCE_PARAM_LINEAR);

    /* If we are marginalising over the time, remove that variable from the model (having set the prior above) */
    /* Also set the prior in model->params, since Likelihood can't access the state! (ugly hack) */
    if(LALInferenceGetProcParamVal(commandLine,"--margtime")){
        LALInferenceVariableItem *p=LALInferenceGetItem(state->priorArgs,"time_min");
        LALInferenceAddVariable(model->params,"time_min",p->value,p->type,p->vary);
        p=LALInferenceGetItem(state->priorArgs,"time_max");
        LALInferenceAddVariable(model->params,"time_max",p->value,p->type,p->vary);
        LALInferenceRemoveVariable(model->params,"time");
    }
    if (LALInferenceGetProcParamVal(commandLine, "--margtimephi") || LALInferenceGetProcParamVal(commandLine, "--margphi")) {
        fprintf(stderr,"ERROR: cannot use margphi or margtimephi with burst approximants. Please use margtime or no marginalization\n");
        exit(1);
    }

    if (LALInferenceGetProcParamVal(commandLine,"--fixsky")){
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--ra"))) ra=(atof(ppt->value));
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--dec"))) dec=(atof(ppt->value));
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--psi"))) psi=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", ra, ramin, ramax, LALINFERENCE_PARAM_FIXED);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", dec, decmin, decmax, LALINFERENCE_PARAM_FIXED);
       LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", psi, psimin, psimax, LALINFERENCE_PARAM_FIXED);
   } 
   else
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, ramin, ramax, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decmin, decmax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psimin, psimax, LALINFERENCE_PARAM_LINEAR);


    /* Set model sampling rates to be consistent with data */
    model->deltaT = state->data->timeData->deltaT;
    model->deltaF = state->data->freqData->deltaF;
    UINT4 nifo=0;
    LALInferenceIFOData *dataPtr = state->data;
    while (dataPtr != NULL)
    {
        dataPtr = dataPtr->next;
        nifo++;
    }

    model->freqhPlus = XLALCreateCOMPLEX16FrequencySeries("freqhPlus",
            &(state->data->freqData->epoch),
            0.0,
            model->deltaF,
            &lalDimensionlessUnit,
            state->data->freqData->data->length);

    model->freqhCross = XLALCreateCOMPLEX16FrequencySeries("freqhCross",
            &(state->data->freqData->epoch),
            0.0,
            model->deltaF,
            &lalDimensionlessUnit,
            state->data->freqData->data->length);

    /* Create arrays for holding single-IFO likelihoods, etc. */
    model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
    model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

    /* Choose proper template */
    model->templt = LALInferenceInitBurstTemplate(state);

    return model;


}

 /* Setup the variables to control Principal component template generation for Spectrograms */
 /* Includes specification of prior ranges */

LALInferenceModel * LALInferenceInitPrincipalCompSpecModel(LALInferenceRunState *state)
{
    printf("-----Using LALInferenceInitPrincipalCompSpecVariables!\n");
 
    LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
    model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
    memset(model->params, 0, sizeof(LALInferenceVariables));

    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;

    /* PC model configuration TODO: sanity check values */
    UINT4 nPCs = 0; // number of principle components to use 
    UINT4 nPCs_max = 20; // max number of PCs
    UINT4 ncatrows = 0; // number of rows in pc matrix (frequency samples) 
    UINT4 ncatcols = 0; // number of columns in pc matrix (waveforms) 
    char *PCfilePlus = NULL; // name of PC file
    char *PCfileCross = NULL; // name of PC file

    if((ppt=LALInferenceGetProcParamVal(commandLine,"--nPCs"))){
        nPCs=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --nPCs");
        exit(1);
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--ncatrows"))){
        ncatrows=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --ncatrows");
        exit(1);
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--ncatcols"))){
        ncatcols=atoi(ppt->value);
    }
    else{
        XLALPrintError("must supply --ncatcols");
        exit(1);
    }

    /* Check PC matrix is compatible with data */
    /*if (state->data->freqData->data->length != ncatrows){
        XLALPrintError("length of F-domain data (%d) does not match length of F-domain PCs (%d)\n",
                state->data->freqData->data->length, ncatrows);
        exit(1);
    }*/

    /* Retrieve PCs */
    model->pcs = XLALMalloc(sizeof(LALInferencePCsModel));
    model->pcs->nPCs = nPCs;

    fprintf(stdout, "PC memory allocated.\n");

    if((ppt=LALInferenceGetProcParamVal(commandLine,"--PCfile_pluspol"))){
        PCfilePlus=ppt->value;
        fprintf(stdout, "Loading PCs from + file\n");
        /* Read PC matrix */
        gsl_matrix_complex *tmp = get_complex_matrix_from_file(PCfilePlus, ncatrows, ncatcols);

        /* Now put first nPCs into a ncatrows x 10 matrix for easily
         * standardised template functions */
	model->pcs->pcs_plus = copy_npcs_from_complex_matrix(tmp, nPCs_max, nPCs, ncatrows, ncatcols);

	/*fprintf(stderr, "[0][0]real = %e , [0][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)));
        fprintf(stderr, "[0][1]real = %e , [0][1]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 1)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 1)));
        fprintf(stderr, "[0][2]real = %e , [0][2]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 2)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 2)));
        fprintf(stderr, "[0][3]real = %e , [0][3]imag = %e\n\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 3)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 3)));

        fprintf(stderr, "[0][0]real = %e , [0][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)));
        fprintf(stderr, "[1][0]real = %e , [1][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 0)));
        fprintf(stderr, "[2][0]real = %e , [2][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 0)));
        fprintf(stderr, "[3][0]real = %e , [3][0]imag = %e\n\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 0)));

	fprintf(stderr, "[0][2]real = %e , [0][2]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 2)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 2)));
        fprintf(stderr, "[1][2]real = %e , [1][2]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 2)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 2)));
        fprintf(stderr, "[2][2]real = %e , [2][2]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 2)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 2)));
        fprintf(stderr, "[3][2]real = %e , [3][2]imag = %e\n\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 2)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 2)));

	fprintf(stderr, "[0][3]real = %e , [0][3]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 3)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 3)));
        fprintf(stderr, "[1][3]real = %e , [1][3]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 3)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 3)));
        fprintf(stderr, "[2][3]real = %e , [2][3]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 3)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 3)));
        fprintf(stderr, "[3][3]real = %e , [3][3]imag = %e\n\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 3)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 3)));*/

//	fprintf(stderr, "Copied them\n");
    }
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--PCfile_crosspol"))){
        PCfileCross=ppt->value;
        fprintf(stdout, "Loading PCs from x file\n");
        /* Read PC matrix */
        gsl_matrix_complex *tmp = get_complex_matrix_from_file(PCfileCross, ncatrows, ncatcols);
        /* Now put first nPCs into a ncatrows x 10 matrix for easily
         * standardised template functions */
        model->pcs->pcs_cross = copy_npcs_from_complex_matrix(tmp, nPCs_max, nPCs, ncatrows, ncatcols);
    }

    fprintf(stderr, "PC files read in.\n");

    if (!LALInferenceGetProcParamVal(commandLine,"--PCfile_pluspol") &&
            !LALInferenceGetProcParamVal(commandLine,"--PCfile_crosspol")){
        XLALPrintError("must supply at least one of --PCfile_pluspol and --PCfile_crosspol\n");
        exit(1);
    }

    REAL8 hrssMin=1.0e-24, hrssMax=1.0e-16; //min was 1.e-26 , max 1.0e-16 
    REAL8 loghrssMin=log(hrssMin),loghrssMax=log(hrssMax);
    REAL8 psimin=0.0,psimax=LAL_PI;
    REAL8 ramin=0.0,ramax=LAL_TWOPI;
    REAL8 decmin=-LAL_PI/2.0,decmax=LAL_PI/2.0;
    REAL8 endtime=0.0;
    REAL8 dt=0.1;
    REAL8 zero=0.0;
    REAL8 ra = 4.6499;
    REAL8 dec = -0.5063;
    REAL8 psi = 0.3046;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime"))) endtime=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--dt"))) dt=atof(ppt->value);

    REAL8 timeMin=endtime-0.5*dt;
    //REAL8 timeMin = endtime - 0.5 * dt + 1.0;
    REAL8 timeMax=endtime+0.5*dt;
    //REAL8 timeMax = endtime + 0.5 * dt + 1.0;


    if (nPCs>=1){

      REAL8 beta1_min=-250, beta1_max=250;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta1_min"))) beta1_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta1_max"))) beta1_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta1", zero, beta1_min, beta1_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=2){

      REAL8 beta2_min=-250, beta2_max=250;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta2_min"))) beta2_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta2_max"))) beta2_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta2", zero, beta2_min, beta2_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=3){

      REAL8 beta3_min=-250, beta3_max=250;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta3_min"))) beta3_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta3_max"))) beta3_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta3", zero, beta3_min, beta3_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=4){

      REAL8 beta4_min=-150, beta4_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta4_min"))) beta4_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta4_max"))) beta4_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta4", zero, beta4_min, beta4_max, LALINFERENCE_PARAM_LINEAR);
    }
    if (nPCs>=5){

      REAL8 beta5_min=-150, beta5_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta5_min"))) beta5_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta5_max"))) beta5_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta5", zero, beta5_min, beta5_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=6){

      REAL8 beta6_min=-150, beta6_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta6_min"))) beta6_min=atof(ppt->value);
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta6_max"))) beta6_max=atof(ppt->value);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta6", zero, beta6_min, beta6_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=7){

      REAL8 beta7_min=-150, beta7_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta7_min"))) beta7_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta7_max"))) beta7_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta7", zero, beta7_min, beta7_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=8){

      REAL8 beta8_min=-150, beta8_max=150;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta8_min"))) beta8_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta8_max"))) beta8_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta8", zero, beta8_min, beta8_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=9){

      REAL8 beta9_min=-100, beta9_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta9_min"))) beta9_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta9_max"))) beta9_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta9", zero, beta9_min, beta9_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=10){

        REAL8 beta10_min=-100, beta10_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta10_min"))) beta10_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta10_max"))) beta10_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta10", zero, beta10_min, beta10_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=11){

        REAL8 beta11_min=-100, beta11_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta11_min"))) beta11_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta11_max"))) beta11_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta11", zero, beta11_min, beta11_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=12){

        REAL8 beta12_min=-100, beta12_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta12_min"))) beta12_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta12_max"))) beta12_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta12", zero, beta12_min, beta12_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=13){

        REAL8 beta13_min=-100, beta13_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta13_min"))) beta13_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta13_max"))) beta13_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta13", zero, beta13_min, beta13_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=14){

        REAL8 beta14_min=-100, beta14_max=100;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta14_min"))) beta14_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta14_max"))) beta14_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta14", zero, beta14_min, beta14_max, LALINFERENCE_PARAM_LINEAR);
    }
    if (nPCs>=15){

        REAL8 beta15_min=-50, beta15_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta15_min"))) beta15_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta15_max"))) beta15_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta15", zero, beta15_min, beta15_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=16){

        REAL8 beta16_min=-50, beta16_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta16_min"))) beta16_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta16_max"))) beta16_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta16", zero, beta16_min, beta16_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=17){

        REAL8 beta17_min=-50, beta17_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta17_min"))) beta17_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta17_max"))) beta17_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta17", zero, beta17_min, beta17_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=18){

        REAL8 beta18_min=-50, beta18_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta18_min"))) beta18_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta18_max"))) beta18_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta18", zero, beta18_min, beta18_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs>=19){

        REAL8 beta19_min=-50, beta19_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta19_min"))) beta19_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta19_max"))) beta19_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta19", zero, beta19_min, beta19_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (nPCs==20){

        REAL8 beta20_min=-50, beta20_max=50;

        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta20_min"))) beta20_min=(atof(ppt->value));
        if((ppt=LALInferenceGetProcParamVal(commandLine,"--beta20_max"))) beta20_max=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "beta20", zero, beta20_min, beta20_max, LALINFERENCE_PARAM_LINEAR);
    }

    if (LALInferenceGetProcParamVal(commandLine,"--use-hrss")){
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--hrssMin"))) hrssMin=(atof(ppt->value));
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--hrssMax"))) hrssMax=(atof(ppt->value));
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "hrss",  zero,hrssMin, hrssMax,   LALINFERENCE_PARAM_LINEAR);
    }
    else
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssMin"))) loghrssMin=(atof(ppt->value));
      if((ppt=LALInferenceGetProcParamVal(commandLine,"--loghrssMax"))) loghrssMax=(atof(ppt->value));
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "loghrss",  zero,loghrssMin, loghrssMax,   LALINFERENCE_PARAM_LINEAR);

      LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", zero, timeMin, timeMax, LALINFERENCE_PARAM_LINEAR);
      //LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", endtime, timeMin, timeMax, LALINFERENCE_PARAM_OUTPUT);

    if (LALInferenceGetProcParamVal(commandLine,"--fixsky")){
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--ra"))) ra=(atof(ppt->value));
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--dec"))) dec=(atof(ppt->value));
       if((ppt=LALInferenceGetProcParamVal(commandLine,"--psi"))) psi=(atof(ppt->value));
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", ra, ramin, ramax, LALINFERENCE_PARAM_FIXED);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", dec, decmin, decmax, LALINFERENCE_PARAM_FIXED);
       LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", psi, psimin, psimax, LALINFERENCE_PARAM_FIXED);
   } 
   else
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, ramin, ramax, LALINFERENCE_PARAM_CIRCULAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decmin, decmax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psimin, psimax, LALINFERENCE_PARAM_LINEAR);

    /* Set model sampling rates to be consistent with data */
    model->deltaT = state->data->timeData->deltaT;
    model->deltaF = state->data->freqData->deltaF;
    UINT4 nifo=0;
    LALInferenceIFOData *dataPtr = state->data;
    while (dataPtr != NULL)
    {
        dataPtr = dataPtr->next;
        nifo++;
    }

    struct stft_data *spec_hplus = malloc(sizeof(struct stft_data));
    double *freqs_hplus = (double *)malloc(sizeof(double)*65);
    double **mag_hplus = (double **)malloc(sizeof(double *)*65);
    for (int n = 0; n < 65; n++) {
      mag_hplus[n]=(double *)malloc(sizeof(double)*936);
    }
    spec_hplus->mag_data = mag_hplus;
    spec_hplus->freqs = freqs_hplus;
    spec_hplus->freq_N = 65;
    spec_hplus->times_N = 936;
    
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime"))) spec_hplus->trig_time=atof(ppt->value);
    
    model->spechPlus = spec_hplus;

    struct stft_data *spec_hcross = malloc(sizeof(struct stft_data));
    double *freqs_hcross = (double *)malloc(sizeof(double)*65);
    double **mag_hcross = (double **)malloc(sizeof(double *)*65);
    for (int n = 0; n < 65; n++) {
      //complex_result[i]=(complex double *)malloc(sizeof(complex double)*num_of_ffts);
      mag_hcross[n]=(double *)malloc(sizeof(double)*936);
    }
    spec_hcross->mag_data = mag_hcross;
    spec_hcross->freqs = freqs_hcross;
    spec_hcross->freq_N = 65;
    spec_hcross->times_N = 936;
    model->spechCross = spec_hcross;

    //fprintf(stderr, "Did the spec model stuff\n");

    /* Create arrays for holding single-IFO likelihoods, etc. */
    model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
    model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

    /* Choose proper template */
    model->templt = LALInferenceInitBurstTemplate(state);

    fprintf(stderr, "Finishing Spectrogram InitBurst.c\n");

    return model;


}





















