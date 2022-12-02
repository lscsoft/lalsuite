/*
 *  LALInferenceCBCInit.c:  Bayesian Followup initialisation routines.
 *
 *  Copyright (C) 2012 Vivien Raymond, John Veitch, Salvatore Vitale
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


#include <stdio.h>
#include <errno.h>
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
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceCalibrationErrors.h>
#include <lal/LALSimNeutronStar.h>

static int checkParamInList(const char *list, const char *param);
static int checkParamInList(const char *list, const char *param)
{
  /* Check for param in comma-seperated list */
  char *post=NULL,*pos=NULL;
  if (list==NULL) return 0;
  if (param==NULL) return 0;

  if(!(pos=strstr(list,param))) return 0;

  /* The string is a substring. Check that it is a token */
  /* Check the character before and after */
  if(pos!=list)
  if(*(pos-1)!=',')
  return 0;

  post=&(pos[strlen(param)]);
  if(*post!='\0')
  if(*post!=',')
  return 0;
  return 1;
}

static void print_flags_orders_warning(SimInspiralTable *injt, ProcessParamsTable *commline);
static void LALInferenceInitSpinVariables(LALInferenceRunState *state, LALInferenceModel *model);
static void LALInferenceInitMassVariables(LALInferenceRunState *state);
static void LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model);
static void LALInferenceCheckApproximantNeeds(LALInferenceRunState *state,Approximant approx);


/* Initialize a bare-bones run-state,
   calling the "ReadData()" function to gather data & PSD from files,
   sets up the random seed and rng, and initializes other variables accordingly.
*/
LALInferenceRunState *LALInferenceInitRunState(ProcessParamsTable *command_line) {
    ProcessParamsTable *ppt=NULL;
    INT4 randomseed;
    FILE *devrandom;
    struct timeval tv;
    LALInferenceRunState *run_state = XLALCalloc(1, sizeof(LALInferenceRunState));

    /* Check that command line is consistent first */
    LALInferenceCheckOptionsConsistency(command_line);
    run_state->commandLine = command_line;

    /* Initialize parameters structure */
    run_state->algorithmParams = XLALCalloc(1, sizeof(LALInferenceVariables));
    run_state->priorArgs = XLALCalloc(1, sizeof(LALInferenceVariables));
    run_state->proposalArgs = XLALCalloc(1, sizeof(LALInferenceVariables));

    /* Read data from files or generate fake data */
    run_state->data = LALInferenceReadData(command_line);
    if (run_state->data == NULL)
        return(NULL);

    /* Setup the random number generator */
    gsl_rng_env_setup();
    run_state->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);

    /* Use clocktime if seed isn't provided */
    ppt = LALInferenceGetProcParamVal(command_line, "--randomseed");
    if (ppt)
        randomseed = atoi(ppt->value);
    else {
        if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
            gettimeofday(&tv, 0);
            randomseed = tv.tv_sec + tv.tv_usec;
        } else {
            fread(&randomseed, sizeof(randomseed), 1, devrandom);
            fclose(devrandom);
        }
    }
    gsl_rng_set(run_state->GSLrandom, randomseed);

    /* Save the random seed */
    LALInferenceAddVariable(run_state->algorithmParams, "random_seed", &randomseed,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    return(run_state);
}

/* Draw initial parameters for each of the threads in run state */
void LALInferenceDrawThreads(LALInferenceRunState *run_state) {
    if (run_state == NULL)
        return;

    if (run_state->threads == NULL) {
        XLALPrintError("Error: LALInferenceDrawThreads expects initialized run_state->threads\n");
        XLAL_ERROR_VOID(XLAL_EINVAL);
    }

    LALInferenceThreadState *thread = &(run_state->threads[0]);
    INT4 t;

    /* If using a malmquist prior, force a strict prior window on distance for starting point, otherwise
     * the approximate prior draws are very unlikely to be within the malmquist prior */
    REAL8 dist_low, dist_high;
    REAL8 restricted_dist_low = log(10.0);
    REAL8 restricted_dist_high = log(100.0);
    INT4 changed_dist = 0;
    if (LALInferenceCheckVariable(run_state->priorArgs, "malmquist") &&
        LALInferenceCheckVariableNonFixed(thread->currentParams, "logdistance")) {
        changed_dist = 1;
        LALInferenceGetMinMaxPrior(run_state->priorArgs, "logdistance", &dist_low, &dist_high);
        LALInferenceRemoveMinMaxPrior(run_state->priorArgs, "logdistance");
        LALInferenceAddMinMaxPrior(run_state->priorArgs, "logdistance",
                                   &restricted_dist_low, &restricted_dist_high, LALINFERENCE_REAL8_t);
    }

    /* If the currentParams are not in the prior, overwrite and pick paramaters
     *   from the priors. OVERWRITE EVEN USER CHOICES.
     *   (necessary for complicated prior shapes where
     *   LALInferenceCyclicReflectiveBound() is not enough) */
    #pragma omp parallel for private(thread)
    for (t = 0; t < run_state->nthreads; t++) {
        LALInferenceVariables *priorDraw = XLALCalloc(1, sizeof(LALInferenceVariables));

        thread = &(run_state->threads[t]);

        /* Try not to clobber values given on the command line */
        LALInferenceCopyVariables(thread->currentParams, priorDraw);
        LALInferenceDrawApproxPrior(thread, priorDraw, priorDraw);
        LALInferenceCopyUnsetREAL8Variables(priorDraw, thread->currentParams,
                                            run_state->commandLine);

        while(isinf(run_state->prior(run_state,
                                     thread->currentParams,
                                     thread->model))) {
            LALInferenceDrawApproxPrior(thread,
                                        thread->currentParams,
                                        thread->currentParams);
        }

        /* Make sure that our initial value is within the
        *     prior-supported volume. */
        LALInferenceCyclicReflectiveBound(thread->currentParams, run_state->priorArgs);

        /* Initialize starting likelihood and prior */
        thread->currentPrior = run_state->prior(run_state,
                                                thread->currentParams,
                                                thread->model);

        thread->currentLikelihood = run_state->likelihood(thread->currentParams,
                                                          run_state->data,
                                                          thread->model);

        LALInferenceClearVariables(priorDraw);
        XLALFree(priorDraw);
    }

    /* Replace distance prior if changed for initial sample draw */
    if (changed_dist) {
        LALInferenceRemoveMinMaxPrior(run_state->priorArgs, "logdistance");
        LALInferenceAddMinMaxPrior(run_state->priorArgs, "logdistance",
                                   &dist_low, &dist_high, LALINFERENCE_REAL8_t);
    }
}

/*
 * Initialize threads in memory, using LALInferenceInitCBCModel() to init models.
 */
void LALInferenceInitCBCThreads(LALInferenceRunState *run_state, INT4 nthreads) {
  if (run_state == NULL){
    LALInferenceInitCBCModel(run_state);
    return;
  }

  ProcessParamsTable *commandLine=run_state->commandLine;

  LALInferenceThreadState *thread;
  INT4 t, nifo;
  INT4 randomseed;
  LALInferenceIFOData *data;
  run_state->nthreads = nthreads;
  run_state->threads = LALInferenceInitThreads(nthreads);

  for (t = 0; t < nthreads; t++) {
    thread = &(run_state->threads[t]);

    /* Link back to run-state */
    thread->parent = run_state;

    /* Set up CBC model and parameter array */
    thread->model = LALInferenceInitCBCModel(run_state);
    thread->model->roq_flag = 0;

    /* Allocate IFO likelihood holders */
    nifo = 0;
    data = run_state->data;
    while (data != NULL) {
        data = data->next;
        nifo++;
    }
    thread->currentIFOLikelihoods = XLALCalloc(nifo, sizeof(REAL8));

    /* Setup ROQ */
    if (LALInferenceGetProcParamVal(commandLine, "--roqtime_steps")){

        LALInferenceSetupROQmodel(thread->model, commandLine);
        fprintf(stderr, "done LALInferenceSetupROQmodel\n");

    }else{
      thread->model->roq_flag=0;
    }

    LALInferenceCopyVariables(thread->model->params, thread->currentParams);
    LALInferenceCopyVariables(run_state->proposalArgs, thread->proposalArgs);

    /* Link thread-state prior-args to the parent run-state's */
    thread->priorArgs = run_state->priorArgs;

    /* Use clocktime if seed isn't provided */
    thread->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
    randomseed = gsl_rng_get(run_state->GSLrandom);
    gsl_rng_set(thread->GSLrandom, randomseed);
  }

  return;
}


/* Setup the template generation */
/* Defaults to using LALSimulation */
LALInferenceTemplateFunction LALInferenceInitCBCTemplate(LALInferenceRunState *runState)
{
  char help[]="(--template [LAL,PhenSpin,LALGenerateInspiral,LALSim,multiband]\tSpecify template (default LAL)\n";
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  //Help is taken care of in LALInferenceInitCBCVariables
  //ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  //if(ppt)
  //{
  //	fprintf(stdout,"%s",help);
  //	return;
  //}
  /* This is the LAL template generator for inspiral signals */
  LALInferenceTemplateFunction templt = &LALInferenceTemplateXLALSimInspiralChooseWaveform;
  ppt=LALInferenceGetProcParamVal(commandLine,"--template");
  if(ppt) {
    if(!strcmp("LALSim",ppt->value))
      templt=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
    else if(!strcmp("null",ppt->value))
        templt=&LALInferenceTemplateNullFreqdomain;
	else if(!strcmp("multiband",ppt->value)){
        templt=&LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated;
        fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated\"\n");
    }
    else {
        XLALPrintError("Error: unknown template %s\n",ppt->value);
        XLALPrintError("%s", help);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
    fprintf(stderr,"Warning: --LALSimulation is deprecated, the LALSimulation package is now the default. To use LALInspiral specify:\n\
                    --template LALGenerateInspiral (for time-domain templates)\n\
                    --template LAL (for frequency-domain templates)\n");
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")){
  templt=&LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence;
        fprintf(stderr, "template is \"LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence\"\n");
  }
  else {
    fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveform\"\n");
  }
  return templt;
}

/* Setup the glitch model */
void LALInferenceInitGlitchVariables(LALInferenceRunState *runState, LALInferenceVariables *currentParams)
{
  ProcessParamsTable    *commandLine   = runState->commandLine;
  LALInferenceIFOData   *dataPtr       = runState->data;
  LALInferenceVariables *priorArgs     = runState->priorArgs;

  UINT4 i,nifo;
  UINT4 n = (UINT4)dataPtr->timeData->data->length;
  UINT4 gflag  = 1;
  REAL8 gmin   = 0.0;
  REAL8 gmax   = 20.0;

  //over-ride default gmax from command line
  if(LALInferenceGetProcParamVal(commandLine, "--glitchNmax"))
    gmax = (REAL8)atoi(LALInferenceGetProcParamVal(commandLine, "--glitchNmax")->value);

  //count interferometers in network before allocating memory
  //compute imin,imax for each IFO -- may be different
  nifo=0;
  dataPtr = runState->data;
  while (dataPtr != NULL)
  {
    dataPtr = dataPtr->next;
    nifo++;
  }
  dataPtr = runState->data;

  UINT4Vector *gsize  = XLALCreateUINT4Vector(nifo);
  //Meyer?? REAL8Vector *gprior = XLALCreateREAL8Vector((int)gmax+1);

  //Morlet??
  gsl_matrix *mAmp = gsl_matrix_alloc(nifo,(int)(gmax));
  gsl_matrix *mf0  = gsl_matrix_alloc(nifo,(int)(gmax));
  gsl_matrix *mQ   = gsl_matrix_alloc(nifo,(int)(gmax));
  gsl_matrix *mt0  = gsl_matrix_alloc(nifo,(int)(gmax));
  gsl_matrix *mphi = gsl_matrix_alloc(nifo,(int)(gmax));

  double Amin,Amax;
  double Qmin,Qmax;
  double f_min,f_max;
  double tmin,tmax;
  double pmin,pmax;
  double Anorm;

  REAL8 TwoDeltaToverN = 2.0 * dataPtr->timeData->deltaT / ((double) dataPtr->timeData->data->length);

  Anorm = sqrt(TwoDeltaToverN);
  Amin = 10.0/Anorm;
  Amax = 10000.0/Anorm;

  Qmin = 3.0;
  Qmax = 30.0;
  tmin = 0.0;
  tmax = dataPtr->timeData->data->length*dataPtr->timeData->deltaT;
  f_min = dataPtr->fLow;
  f_max = dataPtr->fHigh;
  pmin = 0.0;
  pmax = LAL_TWOPI;

  gsl_matrix_set_all(mAmp, Amin);
  gsl_matrix_set_all(mf0,  f_min);
  gsl_matrix_set_all(mQ,   Qmin);
  gsl_matrix_set_all(mt0,  tmin);
  gsl_matrix_set_all(mphi, pmin);

  gsl_matrix  *gFD       = gsl_matrix_alloc(nifo,(int)n); //store the Fourier-domain glitch signal

  for(i=0; i<nifo; i++) gsize->data[i]=0;

  //Morlet wavelet parameters
  LALInferenceAddVariable(currentParams, "morlet_FD",  &gFD,  LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "morlet_Amp", &mAmp, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "morlet_f0" , &mf0,  LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "morlet_Q"  , &mQ,   LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "morlet_t0" , &mt0,  LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "morlet_phi", &mphi, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);

  LALInferenceAddVariable(currentParams, "glitch_size", &gsize, LALINFERENCE_UINT4Vector_t, LALINFERENCE_PARAM_LINEAR);
  LALInferenceAddVariable(currentParams, "glitchFitFlag", &gflag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  LALInferenceAddMinMaxPrior(priorArgs, "morlet_Amp_prior", &Amin, &Amax, LALINFERENCE_REAL8_t);
  LALInferenceAddMinMaxPrior(priorArgs, "morlet_f0_prior" , &f_min, &f_max, LALINFERENCE_REAL8_t);
  LALInferenceAddMinMaxPrior(priorArgs, "morlet_Q_prior"  , &Qmin, &Qmax, LALINFERENCE_REAL8_t);
  LALInferenceAddMinMaxPrior(priorArgs, "morlet_t0_prior" , &tmin, &tmax, LALINFERENCE_REAL8_t);
  LALInferenceAddMinMaxPrior(priorArgs, "morlet_phi_prior", &pmin, &pmax, LALINFERENCE_REAL8_t);

  LALInferenceAddMinMaxPrior(priorArgs, "glitch_size", &gmin, &gmax, LALINFERENCE_REAL8_t);
  LALInferenceAddMinMaxPrior(priorArgs, "glitch_dim", &gmin, &gmax, LALINFERENCE_REAL8_t);

  LALInferenceAddVariable(priorArgs, "glitch_norm", &Anorm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
}

struct spcal_envelope
{
    gsl_spline  *amp_median,*amp_std,
                *phase_median,*phase_std;
    double      logf_min, logf_max;
};

/* Format string for the calibratino envelope file */
/* Frequency    Median Mag     Phase (Rad)    -1 Sigma Mag   -1 Sigma Phase +1 Sigma Mag   +1 Sigma Phase */

#define CAL_ENV_FORMAT "%lf %lf %lf %lf %lf %lf %lf\n"

static struct spcal_envelope *initCalibrationEnvelope(char *filename);

static struct spcal_envelope *initCalibrationEnvelope(char *filename)
{
    FILE *fp=fopen(filename,"r");
    char tmpline[1024];
    if(!fp) {fprintf(stderr,"Unable to open %s: Error %i %s\n",filename,errno,strerror(errno)); exit(1);}
    int Nlines=0;
    REAL8 freq, *logfreq=NULL, *mag_med=NULL, mag_low, mag_hi, *mag_std=NULL, *phase_med=NULL, phase_low, phase_hi, *phase_std=NULL;
    for(Nlines=0;fgets(tmpline,1024,fp); )
    {
        /* Skip header */
        if(tmpline[0]=='#') continue;
        /* Grow arrays */
        logfreq=realloc(logfreq,sizeof(*logfreq)*(Nlines+1));
        mag_med=realloc(mag_med,sizeof(*mag_med)*(Nlines+1));
        mag_std=realloc(mag_std,sizeof(*mag_std)*(Nlines+1));
        phase_med=realloc(phase_med,sizeof(*phase_med)*(Nlines+1));
        phase_std=realloc(phase_std,sizeof(*phase_std)*(Nlines+1));

        if((7!=sscanf(tmpline,CAL_ENV_FORMAT, &freq, &(mag_med[Nlines]), &(phase_med[Nlines]), &mag_low, &phase_low, &mag_hi, &phase_hi)))
        {
            fprintf(stderr,"Malformed input line in file %s: %s\n",filename,tmpline);
            exit(1);
        }
		mag_med[Nlines]-=1.0; /* Subtract off 1 to get delta */
        logfreq[Nlines]=log(freq);
        mag_std[Nlines]=(mag_hi - mag_low ) /2.0;
        phase_std[Nlines]=(phase_hi - phase_low) /2.0;
		Nlines++;
    }
    fprintf(stdout,"Read %i lines from calibration envelope %s\n",Nlines,filename);
    fclose(fp);

    struct spcal_envelope *env=XLALCalloc(1,sizeof(*env));
    env->amp_median = gsl_spline_alloc ( gsl_interp_cspline, Nlines);
    env->amp_std = gsl_spline_alloc ( gsl_interp_cspline, Nlines);
    env->phase_median = gsl_spline_alloc ( gsl_interp_cspline, Nlines);
    env->phase_std = gsl_spline_alloc ( gsl_interp_cspline, Nlines);

    gsl_spline_init(env->amp_median, logfreq, mag_med, Nlines);
    gsl_spline_init(env->amp_std, logfreq, mag_std, Nlines);
    gsl_spline_init(env->phase_median, logfreq, phase_med, Nlines);
    gsl_spline_init(env->phase_std, logfreq, phase_std, Nlines);
    env->logf_min = logfreq[0];
    env->logf_max = logfreq[Nlines-1];

    free(logfreq); free(mag_med); free(mag_std); free(phase_med); free(phase_std);

    return(env);
}

static int destroyCalibrationEnvelope(struct spcal_envelope *env);
static int destroyCalibrationEnvelope(struct spcal_envelope *env)
{
    if(!env) XLAL_ERROR(XLAL_EINVAL);
    if(env->amp_median) gsl_spline_free(env->amp_median);
    if(env->amp_std) gsl_spline_free(env->amp_std);
    if(env->phase_median) gsl_spline_free(env->phase_median);
    if(env->phase_std) gsl_spline_free(env->phase_std);
    XLALFree(env);
    return XLAL_SUCCESS;
}

void LALInferenceInitCalibrationVariables(LALInferenceRunState *runState, LALInferenceVariables *currentParams) {
  ProcessParamsTable *ppt = NULL;
  LALInferenceIFOData *ifo = NULL;
  LALInferenceIFOData *dataPtr=NULL;
  UINT4 calOn = 1;
  if ((ppt = LALInferenceGetProcParamVal(runState->commandLine, "--enable-spline-calibration"))){
    /* Use spline to marginalize*/
    UINT4 ncal = 5; /* Number of calibration nodes, log-distributed
		between fmin and fmax. */
    REAL8 ampUncertaintyPrior = 0.1; /* 10% amplitude */
    REAL8 phaseUncertaintyPrior = 5*M_PI/180.0; /* 5 degrees phase */
    if ((ppt = LALInferenceGetProcParamVal(runState->commandLine, "--spcal-nodes"))) {
      ncal = atoi(ppt->value);
      if (ncal < 3) { /* Cannot do spline with fewer than 3 points! */
	fprintf(stderr, "ERROR: given '--spcal-nodes %d', but cannot spline with fewer than 3\n", ncal);
	exit(1);
      }
    }

    LALInferenceAddVariable(currentParams, "spcal_active", &calOn, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(currentParams, "spcal_npts", &ncal, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    for(ifo=runState->data;ifo;ifo=ifo->next) {
      UINT4 i;

      char freqVarName[VARNAME_MAX];
      char ampVarName[VARNAME_MAX];
      char phaseVarName[VARNAME_MAX];

      REAL8 fMin = ifo->fLow;
      REAL8 fMax = ifo->fHigh;
      REAL8 logFMin = log(fMin);
      REAL8 logFMax = log(fMax);

      char amp_uncert_op[VARNAME_MAX];
      char pha_uncert_op[VARNAME_MAX];
      char env_uncert_op[VARNAME_MAX];
      struct spcal_envelope *env=NULL;

      if((VARNAME_MAX <= snprintf(amp_uncert_op, VARNAME_MAX, "--%s-spcal-amp-uncertainty", ifo->name))
          || (VARNAME_MAX <= snprintf(pha_uncert_op, VARNAME_MAX, "--%s-spcal-phase-uncertainty", ifo->name))
          || (VARNAME_MAX <= snprintf(env_uncert_op, VARNAME_MAX, "--%s-spcal-envelope",ifo->name)) )
      {
        fprintf(stderr,"variable name too long\n"); exit(1);
      }

      if( (ppt=LALInferenceGetProcParamVal(runState->commandLine, env_uncert_op)))
      {
          env = initCalibrationEnvelope(ppt->value);
          /* If the spline limits are within the data frequency bounds, update the bounds */
          logFMin = env->logf_min>logFMin ? env->logf_min : logFMin;
          logFMax = env->logf_max<logFMax ? env->logf_max : logFMax;
      }
      else
      {
        if ((ppt = LALInferenceGetProcParamVal(runState->commandLine, amp_uncert_op))) {
            ampUncertaintyPrior = atof(ppt->value);
        }
        else{
            fprintf(stderr,"Error, missing %s or %s\n",amp_uncert_op, env_uncert_op);
            exit(1);
        }

        if ((ppt = LALInferenceGetProcParamVal(runState->commandLine, pha_uncert_op))) {
            phaseUncertaintyPrior = M_PI/180.0*atof(ppt->value); /* CL arg in degrees, variable in radians */
        }
        else{
            fprintf(stderr,"Error, missing %s or %s\n",pha_uncert_op,env_uncert_op);
            exit(1);
        }
      }
      REAL8 dLogF = (logFMax - logFMin)/(ncal-1);
      /* Now add each spline node */
      for(i=0;i<ncal;i++)
	  {
			  if((VARNAME_MAX <= snprintf(freqVarName, VARNAME_MAX, "%s_spcal_logfreq_%i",ifo->name,i))
                  || (VARNAME_MAX <= snprintf(ampVarName, VARNAME_MAX, "%s_spcal_amp_%i", ifo->name,i))
                  || (VARNAME_MAX <= snprintf(phaseVarName, VARNAME_MAX, "%s_spcal_phase_%i", ifo->name,i)) )
              {
                  fprintf(stderr,"Variable name too long\n"); exit(1);
              }
              REAL8 amp_std=ampUncertaintyPrior,amp_mean=0.0;
			  REAL8 phase_std=phaseUncertaintyPrior,phase_mean=0.0;
			  REAL8 logFreq = logFMin + i*dLogF;
			  LALInferenceAddREAL8Variable(currentParams,freqVarName,logFreq,LALINFERENCE_PARAM_FIXED);
			  if(env)
			  {
					  amp_std = gsl_spline_eval(env->amp_std, logFreq, NULL);
					  amp_mean = gsl_spline_eval(env->amp_median, logFreq, NULL);
					  phase_std = gsl_spline_eval(env->phase_std, logFreq, NULL);
					  phase_mean = gsl_spline_eval(env->phase_median, logFreq, NULL);
			  }
			  LALInferenceRegisterGaussianVariableREAL8(runState, currentParams, ampVarName, 0, amp_mean, amp_std, LALINFERENCE_PARAM_LINEAR);
			  LALInferenceRegisterGaussianVariableREAL8(runState, currentParams, phaseVarName, 0, phase_mean, phase_std, LALINFERENCE_PARAM_LINEAR);
	  } /* End loop over spline nodes */

	  if(env) destroyCalibrationEnvelope(env);
	} /* End loop over IFOs */
  } /* End case of spline calibration error */
  else if(LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp") ||LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha")){
    /* Use constant (in frequency) approximation for the errors */
    if (LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp")){
      /*For the moment the prior ranges are the same for the three IFOs */
      REAL8 camp_max_A=0.25; /* plus minus 25% amplitude errors*/
      REAL8 camp_min_A=-0.25;
      REAL8 zero=0.0;
      dataPtr = runState->data;
      while (dataPtr != NULL){
        char CA_A[320];
        sprintf(CA_A,"%s_%s","calamp",dataPtr->name);
        LALInferenceRegisterUniformVariableREAL8(runState, currentParams, CA_A, zero, camp_min_A, camp_max_A, LALINFERENCE_PARAM_LINEAR);
        dataPtr = dataPtr->next;
      }

      LALInferenceAddVariable(currentParams, "constantcal_active", &calOn, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
      /*If user specifies a width for the error prior, a gaussian prior will be used, otherwise a flat prior will be used*/
      REAL8 ampUncertaintyPrior=-1.0;
      ppt = LALInferenceGetProcParamVal(runState->commandLine, "--constcal_ampsigma");
      if (ppt) {
        ampUncertaintyPrior = atof(ppt->value);
      }
      LALInferenceAddVariable(runState->priorArgs, "constcal_amp_uncertainty", &ampUncertaintyPrior, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    }
    if (LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha")){
      /* Add linear calibration phase errors to the measurement. For the moment the prior ranges are the same for the three IFOs */
      REAL8 cpha_max_A=0.349;  /* plus/minus 20 degs*/
      REAL8 cpha_min_A=-0.349;
      REAL8 zero=0.0;
      dataPtr = runState->data;
      while (dataPtr != NULL)
      {
        char CP_A[320];
        sprintf(CP_A,"%s_%s","calpha",dataPtr->name);
        LALInferenceRegisterUniformVariableREAL8(runState, currentParams, CP_A, zero, cpha_min_A, cpha_max_A, LALINFERENCE_PARAM_LINEAR);
        dataPtr = dataPtr->next;
      }
      LALInferenceAddVariable(currentParams, "constantcal_active", &calOn, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

     /*If user specifies a width for the error prior, a gaussian prior will be used, otherwise a flat prior will be used*/
      REAL8 phaseUncertaintyPrior=-1.0;
      ppt = LALInferenceGetProcParamVal(runState->commandLine, "--constcal_phasigma");
      if (ppt) {
        phaseUncertaintyPrior = M_PI/180.0*atof(ppt->value); /* CL arg in degrees, variable in radians */
      }
      LALInferenceAddVariable(runState->priorArgs, "constcal_phase_uncertainty", &phaseUncertaintyPrior, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    }
  }
  else{
    /* No calibration marginalization asked. Just exit */
    return;
  }
}

void LALInferenceRegisterGaussianVariableREAL8(LALInferenceRunState *state, LALInferenceVariables *var, const char *name, REAL8 startval, REAL8 mean, REAL8 stdev, LALInferenceParamVaryType varytype)
{
  char meanopt[VARNAME_MAX+8];
  char sigmaopt[VARNAME_MAX+9];
  char valopt[VARNAME_MAX+3];
  char fixopt[VARNAME_MAX+7];
  ProcessParamsTable *ppt=NULL;

  sprintf(meanopt,"--%s-mean",name);
  sprintf(sigmaopt,"--%s-sigma",name);
  sprintf(valopt,"--%s",name);
  sprintf(fixopt,"--fix-%s",name);

  if((ppt=LALInferenceGetProcParamVal(state->commandLine,meanopt))) mean=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,sigmaopt))) stdev=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,fixopt)))
  {
    varytype = LALINFERENCE_PARAM_FIXED;
    startval = atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,valopt))) startval=atof(ppt->value);

  XLAL_CHECK_ABORT(stdev>0);
  LALInferenceAddVariable(var,name,&startval,LALINFERENCE_REAL8_t,varytype);
  LALInferenceAddGaussianPrior(state->priorArgs, name, &mean, &stdev, LALINFERENCE_REAL8_t);

}

void LALInferenceRegisterUniformVariableREAL8(LALInferenceRunState *state, LALInferenceVariables *var, const char *name, REAL8 startval, REAL8 min, REAL8 max, LALInferenceParamVaryType varytype)
{
  char minopt[VARNAME_MAX+7];
  char maxopt[VARNAME_MAX+7];
  char valopt[VARNAME_MAX+3];
  char fixopt[VARNAME_MAX+7];
  ProcessParamsTable *ppt=NULL;

  sprintf(minopt,"--%s-min",name);
  sprintf(maxopt,"--%s-max",name);
  sprintf(valopt,"--%s",name);
  sprintf(fixopt,"--fix-%s",name);

  if((ppt=LALInferenceGetProcParamVal(state->commandLine,minopt))) min=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,maxopt))) max=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,fixopt))) {
		  varytype=LALINFERENCE_PARAM_FIXED;
		  startval=atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(state->commandLine,valopt))) startval=atof(ppt->value);
  else if(varytype!=LALINFERENCE_PARAM_FIXED) startval=min+(max-min)*gsl_rng_uniform(state->GSLrandom);

  /* Error checking */
  if(min>max) {
    fprintf(stderr,"ERROR: Prior for %s has min(%lf) > max(%lf)\n",name,min,max);
    exit(1);
  }
  if(startval<min || startval>max){
    fprintf(stderr,"ERROR: Initial value %lf for %s lies outwith prior (%lf,%lf)\n",startval,name,min,max);
    exit(1);
  }
  /* Mass parameters checks*/
  if (!strcmp(name,"eta"))
    if (max>0.25){
      fprintf(stderr,"ERROR: maximum of eta cannot be larger than 0.25. Check --eta-max\n");
      exit(1);
    }
  if (!strcmp(name,"q")){
    REAL8 qMin=min;
    REAL8 qMax=max;

    if (qMin <= 0.0 || qMin > 1.0)
    {
        fprintf(stderr,"ERROR: qMin must be between 0 and 1, got value qMin=%f\n",qMin);
		exit(1);
    }
    if (qMax > 1.0 || qMax <0.0 || qMax < qMin)
    {
      fprintf(stderr,"ERROR: qMax must be between 0 and 1, and qMax > qMin. Got value qMax=%f, qMin=%f\n",qMax,qMin);
	  exit(1);
    }
  }
  /*End of mass parameters check */

  LALInferenceAddVariable(var,name,&startval,LALINFERENCE_REAL8_t,varytype);
  LALInferenceAddMinMaxPrior(state->priorArgs, name, &min, &max, LALINFERENCE_REAL8_t);

}


/* Setup the variables to control template generation for the CBC model */
/* Includes specification of prior ranges. Returns address of new LALInferenceVariables */

LALInferenceModel *LALInferenceInitCBCModel(LALInferenceRunState *state) {
  char help[]="\
    ----------------------------------------------\n\
    --- Injection Arguments ----------------------\n\
    ----------------------------------------------\n\
    (--inj injections.xml) Injection XML file to use\n\
    (--event N)            Event number from Injection XML file to use\n\
\n\
    ----------------------------------------------\n\
    --- Template Arguments -----------------------\n\
    ----------------------------------------------\n\
    (--use-eta)            Jump in symmetric mass ratio eta, instead of q=m1/m2 (m1>m2)\n\
    (--approx)             Specify a template approximant and phase order to use\n\
                         (default TaylorF2threePointFivePN). Available approximants:\n\
                         default modeldomain=\"time\": GeneratePPN, TaylorT1, TaylorT2,\n\
                                                       TaylorT3, TaylorT4, EOB, EOBNR,\n\
                                                       EOBNRv2, EOBNRv2HM, SEOBNRv1,\n\
                                                       SpinTaylor, SpinQuadTaylor, \n\
                                                       SpinTaylorFrameless, SpinTaylorT4,\n\
                                                       PhenSpinTaylorRD, NumRel.\n\
                         default modeldomain=\"frequency\": TaylorF1, TaylorF2, TaylorF2RedSpin,\n\
                                                       TaylorF2RedSpinTidal, IMRPhenomA,\n\
                                                       IMRPhenomB, IMRPhenomP, IMRPhenomPv2.\n\
    (--amporder PNorder)            Specify a PN order in amplitude to use (defaults: LALSimulation: max available; LALInspiral: newtownian).\n\
    (--fref f_ref)                  Specify a reference frequency at which parameters are defined (default 100).\n\
    (--tidal)                   Enables tidal corrections, only with LALSimulation.\n\
    (--tidalT)                  Enables reparmeterized tidal corrections, only with LALSimulation.\n\
    (--4PolyEOS)                Enables 4-piece polytropic EOS parmeterization, only with LALSimulation.\n\
    (--4SpectralDecomp)         Enables 4-coeff. spectral decomposition EOS parmeterization, only with LALSimulation.\n\
    (--eos   EOS)               Fix the neutron star EOS. Use \"--eos help\" for allowed names\n\
    (--BinaryLove)              Enable the Binary Neutron Star common EoS tidal model, expressed through EOS-insensitive relations\n\
    (--spinOrder PNorder)           Specify twice the PN order (e.g. 5 <==> 2.5PN) of spin effects to use, only for LALSimulation (default: -1 <==> Use all spin effects).\n\
    (--tidalOrder PNorder)          Specify twice the PN order (e.g. 10 <==> 5PN) of tidal effects to use, only for LALSimulation (default: -1 <==> Use all tidal effects).\n\
    (--numreldata FileName)         Location of NR data file for NR waveforms (with NR_hdf5 approx).\n\
    (--modeldomain)                 domain the waveform template will be computed in (\"time\" or \"frequency\"). If not given will use LALSim to decide\n\
    (--spinAligned or --aligned-spin)  template will assume spins aligned with the orbital angular momentum.\n\
    (--singleSpin)                  template will assume only the spin of the most massive binary component exists.\n\
    (--noSpin, --disable-spin)      template will assume no spins (giving this will void spinOrder!=0) \n\
    (--no-detector-frame)              model will NOT use detector-centred coordinates and instead RA,dec\n\
    (--nonGR_alpha value) this is a LIV parameter which should only be passed when log10lambda_eff/lambda_eff is passed as a grtest-parameter for LIV test\n\
    (--LIV_A_sign) this is a LIV parameter determining if +A or -A is being tested; A occurs in the modified dispersion relation. LIV_A_sign has to be either +1 or -1 \n\
    (--liv) this flag is now needed for launching a LIV run\n\
    (--grtest-parameters dchi0,..,dxi1,..,dalpha1,..) template will assume deformations in the corresponding phase coefficients.\n\
    (--ppe-parameters aPPE1,....     template will assume the presence of an arbitrary number of PPE parameters. They must be paired correctly.\n\
    (--modeList lm,l-m...,lm,l-m)           List of modes to be used by the model. The chosen modes ('lm') should be passed as a ',' seperated list.\n\
    (--phenomXHMMband float)     Threshold parameter for the Multibanding of the non-precessing hlm modes in IMRPhenomXHM and IMRPhenomXPHM. If set to 0 then do not use multibanding.\n\
                                 Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html\n\
    (--phenomXPHMMband float)    Threshold parameter for the Multibanding of the Euler angles in IMRPhenomXPHM. If set to 0 then do not use multibanding.\n\
                                 Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html\n\
    (--phenomXPFinalSpinMod int) Change version for the final spin model used in IMRPhenomXP/IMRPhenomXPHM.\n\
                                 Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html\n\
    (--phenomXPrecVersion int)   Change version of the Euler angles for the twisting-up of IMRPhenomXP/IMRPhenomXPHM.\n\
                                 Options and default values can be found in https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom_x__c.html\n\
\n\
    ----------------------------------------------\n\
    --- Starting Parameters ----------------------\n\
    ----------------------------------------------\n\
    You can generally have MCMC chains to start from a given parameter value by using --parname VALUE. Names currently known to the code are:\n\
     time                         Waveform time (overrides random about trigtime).\n\
     chirpmass                    Chirpmass\n\
     eta                          Symmetric massratio (needs --use-eta)\n\
     q                            Asymmetric massratio (a.k.a. q=m2/m1 with m1>m2)\n\
     phase                        Coalescence phase.\n\
     costheta_jn                  Cosine of angle between J and line of sight [rads]\n\
     logdistance                  Log Distance (requires --use-logdistance)\n\
     rightascension               Rightascensions\n\
     declination                  Declination.\n\
     polarisation                 Polarisation angle.\n\
    * Spin Parameters:\n\
     a_spin1                      Spin1 magnitude\n\
     a_spin2                      Spin2 magnitude\n\
     tilt_spin1                   Angle between spin1 and orbital angular momentum\n\
     tilt_spin2                   Angle between spin2 and orbital angular momentum \n\
     phi_12                       Difference between spins' azimuthal angles \n\
     phi_jl                       Difference between total and orbital angular momentum azimuthal angles\n\
    * Equation of State parameters:\n\
     (requires --tidal)\n\
     lambda1                      lambda1.\n\
     lambda2                      lambda2.\n\
     (requires --tidalT)\n\
     lambdaT                      lambdaT.\n\
     dLambdaT                     dLambdaT.\n\
     (requires --4PolyEOS)\n\
     logp1                        logp1.\n\
     gamma1                       gamma1.\n\
     gamma2                       gamma2.\n\
     gamma3                       gamma3.\n\
     (requires --BinaryLove)\n\
     lambdaS                      Symmetric tidal deformability.\n\
     (requires --4SpectralDecomp):\n\
     SDgamma0                     SDgamma0.\n\
     SDgamma1                     SDgamma1.\n\
     SDgamma2                     SDgamma2.\n\
     SDgamma3                     SDgamma3.\n\
    * \n\
    * Spin-induced quadrupole moment test parameters:\n\
     (requires --dQuadMon12)\n\
     dQuadMon1                      dQuadMon1.\n\
     dQuadMon2                      dQuadMon2.\n\
     (requires --dQuadMonSA) \n\
     dQuadMonS                      dQuadMonS.\n\
     dQuadMonA                      dQuadMonA.\n\
    (dQuadMonS and dQuadMonA are the symmetric and antisymmetric combinations of dQuadMon1 and dQuadMon2).\n\
    ----------------------------------------------\n\
    --- Prior Ranges -----------------------------\n\
    ----------------------------------------------\n\
    You can generally use --paramname-min MIN --paramname-max MAX to set the prior range for the parameter paramname\n\
    The names known to the code are listed below.\n\
    Component masses, total mass and time have dedicated options listed here:\n\n\
    (--trigtime time)                       Center of the prior for the time variable.\n\
    (--comp-min min)                        Minimum component mass (1.0).\n\
    (--comp-max max)                        Maximum component mass (100.0).\n\
    (--mass1-min min, --mass1-max max)      Min and max for mass1 (default: same as comp-min,comp-max, will over-ride these.\n\
    (--mass2-min min, --mass2-max max)      Min and max for mass2 (default: same as comp-min,comp-max, will over-ride these.\n\
    (--mtotal-min min)                      Minimum total mass (2.0).\n\
    (--mtotal-max max)                      Maximum total mass (200.0).\n\
    (--dt time)                             Width of time prior, centred around trigger (0.2s).\n\
\n\
    (--ns-max-mass)                          Maximum observed NS mass used for EOS prior (none).\n\
    (--varyFlow, --flowMin, --flowMax)       Allow the lower frequency bound of integration to vary in given range.\n\
    (--pinparams)                            List of parameters to set to injected values [mchirp,asym_massratio,etc].\n\
    ----------------------------------------------\n\
    --- Fix Parameters ---------------------------\n\
    ----------------------------------------------\n\
    You can generally fix a parameter to be fixed to a given values by using --fix-paramname VALUE\n\
    where the known names have been listed above.\n\
\n";

  /* Print command line arguments if state was not allocated */
  if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  LALStatus status;
  memset(&status,0,sizeof(status));
  int errnum;
  SimInspiralTable *injTable=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;
  LALInferenceVariables *proposalArgs=state->proposalArgs;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *ppt_order=NULL;
  LALPNOrder PhaseOrder=-1;
  LALPNOrder AmpOrder=-1;
  Approximant approx=NumApproximants;
  REAL8 f_ref = 100.0;
  LALInferenceIFOData *dataPtr;
  UINT4 event=0;
  UINT4 i=0;
  /* Default priors */
  REAL8 Dmin=1.0;
  REAL8 Dmax=2000.0;
  REAL8 Dinitial = (Dmax + Dmin)/2.0;
  REAL8 mcMin=1.0;
  REAL8 mcMax=15.3;
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 qMin=1./30.; // The ratio between min and max component mass (see InitMassVariables)
  REAL8 qMax=1.0;
  REAL8 psiMin=0.0,psiMax=LAL_PI;
  REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0;
  REAL8 raMin=0.0,raMax=LAL_TWOPI;
  REAL8 phiMin=0.0,phiMax=LAL_TWOPI;
  REAL8 costhetaJNmin=-1.0 , costhetaJNmax=1.0;
  REAL8 dt=0.1;  /* Half the width of time prior */
  REAL8 lambda1Min=0.0;
  REAL8 lambda1Max=3000.0;
  REAL8 lambda2Min=0.0;
  REAL8 lambda2Max=3000.0;
  REAL8 lambdaTMin=0.0;
  REAL8 lambdaTMax=3000.0;
  REAL8 dLambdaTMin=-500.0;
  REAL8 dLambdaTMax=500.0;
  REAL8 dQuadMonMin=-200.0;
  REAL8 dQuadMonMax=200.0;
  REAL8 logp1Min=33.6;
  REAL8 logp1Max=35.4;
  REAL8 gamma1Min=2.0;
  REAL8 gamma1Max=4.5;
  REAL8 gamma2Min=1.1;
  REAL8 gamma2Max=4.5;
  REAL8 gamma3Min=1.1;
  REAL8 gamma3Max=4.5;
  REAL8 lambdaSMin=0.0;
  REAL8 lambdaSMax=5000.0;
  REAL8 SDgamma0Min=0.2;
  REAL8 SDgamma0Max=2.0;
  REAL8 SDgamma1Min=-1.6;
  REAL8 SDgamma1Max=1.7;
  REAL8 SDgamma2Min=-0.6;
  REAL8 SDgamma2Max=0.6;
  REAL8 SDgamma3Min=-0.02;
  REAL8 SDgamma3Max=0.02;
  gsl_rng *GSLrandom=state->GSLrandom;
  REAL8 endtime=0.0, timeParam=0.0;
  REAL8 timeMin=endtime-dt,timeMax=endtime+dt;
  REAL8 zero=0.0; /* just a number that will be overwritten anyway*/

  /* Over-ride prior bounds if analytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence(state));
  }
  else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence_bimod(state));
  }
  else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence(state)); /* CHECKME: Use the default prior for unimodal */
  }

  LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  memset(model->params, 0, sizeof(LALInferenceVariables));
  model->eos_fam = NULL;

  UINT4 signal_flag=1;
  ppt = LALInferenceGetProcParamVal(commandLine, "--noiseonly");
  if(ppt)signal_flag=0;
  LALInferenceAddVariable(model->params, "signalModelFlag", &signal_flag,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);

  /* Read injection XML file for parameters if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(ppt){
    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
    if(!injTable){
      fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
      exit(1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--event");
    if(ppt){
      event= atoi(ppt->value);
      fprintf(stderr,"Reading event %d from file\n",event);
      i=0;
      while(i<event) {i++; injTable=injTable->next;} /* select event */
    }
  }

  /* See if there are any parameters pinned to injection values */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    char *pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    char **strings=NULL;
    UINT4 N;
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    LALInferenceInjectionToVariables(injTable,&tempParams);
    LALInferenceVariableItem *node=NULL;
    while(N>0){
      N--;
      char *name=strings[N];
      fprintf(stdout,"Pinning parameter %s\n",name);
      node=LALInferenceGetItem(&tempParams,name);
      if(node) LALInferenceAddVariable(model->params,node->name,node->value,node->type,node->vary);
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",name);}
    }
  }

  /* Over-ride approximant if user specifies */
  ppt=LALInferenceGetProcParamVal(commandLine,"--approximant");
  if(ppt){
    approx = XLALGetApproximantFromString(ppt->value);
    ppt_order=LALInferenceGetProcParamVal(commandLine,"--order");
    if(ppt_order) PhaseOrder = XLALGetOrderFromString(ppt_order->value);
  }
  ppt=LALInferenceGetProcParamVal(commandLine,"--approx");
  if(ppt){
    approx = XLALGetApproximantFromString(ppt->value);
    XLAL_TRY(PhaseOrder = XLALGetOrderFromString(ppt->value),errnum);
    if( (int) PhaseOrder == XLAL_FAILURE || errnum) {
      PhaseOrder=-1;
    }
  }

  ppt=LALInferenceGetProcParamVal(commandLine,"--amporder");
  if(ppt) AmpOrder=atoi(ppt->value);

  if(approx==NumApproximants && injTable){ /* Read aproximant from injection file */
    approx=XLALGetApproximantFromString(injTable->waveform);
  }
  if(approx==NumApproximants){
       approx=TaylorF2; /* Defaults to TF2 */
       XLALPrintWarning("You did not provide an approximant for the templates. Using default %s, which might now be what you want!\n",XLALSimInspiralGetStringFromApproximant(approx));
  }

  /* Set the model domain appropriately */
  if (XLALSimInspiralImplementedFDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  } else if (XLALSimInspiralImplementedTDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_TIME;
  } else {
    fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
    exit(1);
  }

  ppt=LALInferenceGetProcParamVal(commandLine, "--fref");
  if (ppt) f_ref = atof(ppt->value);

  ppt=LALInferenceGetProcParamVal(commandLine,"--modeldomain");
  if(ppt){
    if ( ! strcmp( "time", ppt->value ) )
    {
      model->domain = LAL_SIM_DOMAIN_TIME;
    }
    else if ( ! strcmp( "frequency", ppt->value ) )
    {
      model->domain = LAL_SIM_DOMAIN_FREQUENCY;
    }
    else
    {
      fprintf( stderr, "invalid argument to --modeldomain:\n"
              "unknown domain specified: "
              "domain must be one of: time, frequency\n");
      exit( 1 );
    }
  }

  /* This sets the component masses and total mass priors, if given in command line.
   * The prior for other parameters are now read in in RegisterUniformVariable, if given by the user. */
  LALInferenceInitMassVariables(state);
  /* now we need to update the chirp mass and q limits accordingly */
  REAL8 m2_min=LALInferenceGetREAL8Variable(state->priorArgs,"mass2_min");
  REAL8 m1_max=LALInferenceGetREAL8Variable(state->priorArgs,"mass1_max");
  REAL8 mtot_min = *(REAL8 *)LALInferenceGetVariable(state->priorArgs,"MTotMin");
  REAL8 mtot_max = *(REAL8 *)LALInferenceGetVariable(state->priorArgs,"MTotMax");
  qMin = m2_min / m1_max;
  mcMin =mtot_min*pow(qMin/pow(1.+qMin,2.),3./5.);
  mcMax =mtot_max*pow(0.25,3./5.);

  /************ Initial Value Related Argument START *************/
  /* Read time parameter from injection file */
  if(injTable)
  {
    endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
    fprintf(stdout,"Using end time from injection file: %lf\n", endtime);
  }
  /* Over-ride end time if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if(ppt){
    endtime=atof(ppt->value);
    printf("Read end time %f\n",endtime);
  }
  /* Over-ride time prior window if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt)
    dt=atof(ppt->value);
  timeMin=endtime-dt; timeMax=endtime+dt;
  timeParam = timeMin + (timeMax-timeMin)*gsl_rng_uniform(GSLrandom);

  /* Initial Value Related END */
  LALInferenceAddVariable(model->params, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "LAL_AMPORDER",     &AmpOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "f_ref", &f_ref, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

  /* Get frequency bounds */
  REAL8 fLow = INFINITY; // lowest frequency being analyzed across the network
  REAL8 fHigh = -INFINITY; // highest frequency being analyzed across the network

  dataPtr = state->data;
  while (dataPtr != NULL)
  {
    if (dataPtr->fLow < fLow)
        fLow = dataPtr->fLow;
    if (dataPtr->fHigh > fHigh)
        fHigh = dataPtr->fHigh;

    dataPtr = dataPtr->next;
  }

  model->fLow = fLow;
  model->fHigh = fHigh;

  /* Check if flow is varying */
  ppt=LALInferenceGetProcParamVal(commandLine,"--vary-flow");
  if(ppt){
    REAL8 fLow_min = fLow;
    REAL8 fLow_max = 200.0;
    if(LALInferenceCheckVariable(model->params,"f_ref"))
      f_ref = *(REAL8*)LALInferenceGetVariable(model->params, "f_ref");
      if (f_ref > 0.0 && fLow_max > f_ref) {
        fprintf(stdout,"WARNING: flow can't go higher than the reference frequency.  Setting flow-max to %f\n",f_ref);
        fLow_max = f_ref;
      }
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "flow", fLow, fLow_min, fLow_max, LALINFERENCE_PARAM_LINEAR);
  } else {
    LALInferenceAddVariable(model->params, "flow", &fLow,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }
  /* Set up the variable parameters */

  /********************* TBL: Adding noise-fitting parameters  *********************/
  UINT4 nscale_block; //number of noise parameters per IFO (1 per frequency block)
  UINT4 nscale_bin;   //number of Fourier bins in each noise block
  REAL8 nscale_dflog; //logarithmic spacing for noise parameters
  REAL8 nscale_min;   //minimum value for psd scale parameter
  REAL8 nscale_max;   //maximum value for psd scale parameters
  UINT4 nscale_dim;   //total dimension of noise model (params X detectors)
  UINT4 nscale_flag;  //flag to tell likelihood function if psd fitting is in use

  REAL8Vector *nscale_prior = NULL; //std. dev. of prior distribution
  REAL8Vector *nscale_sigma = NULL; //std. dev. of prior distribution

  //assume no noise fitting
  nscale_flag=0;

  //set Nblock to default unless specified at command line
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdNblock");
  if(ppt) nscale_block = atoi(ppt->value);
  else nscale_block = 8;

  //First, figure out sizes of dataset to set up noise blocks
  UINT4 nifo; //number of data channels
  UINT4 imin; //minimum Fourier bin for integration in IFO
  UINT4 imax; //maximum Fourier bin for integration in IFO
  UINT4 f_min = 1; //minimum Fourier bin for integration over network
  UINT4 f_max = 1; //maximum Fourier bin for integration over network
  REAL8 df = 1.0; //frequency resolution

  /* Set model sampling rates to be consistent with data */
  model->deltaT = state->data->timeData->deltaT;
  model->deltaF = state->data->freqData->deltaF;

  /* Get number of interferometers */
  nifo=0;
  dataPtr = state->data;
  while (dataPtr != NULL)
  {
    df      = 1.0 / (((double)dataPtr->timeData->data->length) * model->deltaT);
    imin    = (UINT4)ceil( dataPtr->fLow  / df);
    imax    = (UINT4)floor(dataPtr->fHigh / df);

    if(nifo==0)
    {
      f_min=imin;
      f_max=imax;
    }
    else
    {
      if(imin<f_min)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_min=imin;
      }
      if(imax>f_max)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_max=imax;
      }
    }

    dataPtr = dataPtr->next;
    nifo++;
  }
  imin = f_min;
  imax = f_max;

  UINT4 j = 0;

  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFit");
  if (ppt) {
      printf("WARNING: --psdFit has been deprecated in favor of --psd-fit\n");
  } else {
      ppt = LALInferenceGetProcParamVal(commandLine, "--psd-fit");
  }
  if(ppt)//MARK: Here is where noise PSD parameters are being added to the model
  {

    printf("Setting up PSD fitting for %i ifos...\n",nifo);

    dataPtr = state->data;

    gsl_matrix *bands_min = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix *bands_max = gsl_matrix_alloc(nifo,nscale_block);

    i=0;
    while (dataPtr != NULL)
    {
      printf("ifo=%i  %s\n",i,dataPtr->name);fflush(stdout);

        nscale_bin   = (imax+1-imin)/nscale_block;
        nscale_dflog = log( (double)(imax+1)/(double)imin )/(double)nscale_block;

        int freq_min, freq_max;

        for (j = 0; j < nscale_block; j++)
        {

            freq_min = (int) exp(log((double)imin ) + nscale_dflog*j);
            freq_max = (int) exp(log((double)imin ) + nscale_dflog*(j+1));

            gsl_matrix_set(bands_min,i,j,freq_min);
            gsl_matrix_set(bands_max,i,j,freq_max);
        }


      dataPtr = dataPtr->next;
      i++;

    }

    printf("Running PSD fitting with bands (Hz)...\n");
    dataPtr = state->data;
    i=0;
    while (dataPtr != NULL)
    {
      printf("%s:",dataPtr->name);
      for (j = 0; j < nscale_block; j++)
      {
        printf(" %f-%f ",gsl_matrix_get(bands_min,i,j)*df,gsl_matrix_get(bands_max,i,j)*df);
      }
      printf("\n");
      dataPtr = dataPtr->next;
      i++;
    }

    nscale_bin   = (f_max+1-f_min)/nscale_block;
    nscale_dflog = log( (double)(f_max+1)/(double)f_min )/(double)nscale_block;

    nscale_min   = 1.0e-1;
    nscale_max   = 1.0e+1;
    nscale_dim   = nscale_block*nifo;
    nscale_flag  = 1;

    // Set noise parameter arrays.
    nscale_prior = XLALCreateREAL8Vector(nscale_block);
    nscale_sigma = XLALCreateREAL8Vector(nscale_block);
    for(i=0; i<nscale_block; i++)
    {
      nscale_prior->data[i] = 1.0/sqrt( gsl_matrix_get(bands_max,0,i)-gsl_matrix_get(bands_min,0,i) );
      nscale_sigma->data[i] = nscale_prior->data[i]/sqrt((double)(nifo*nscale_block));
    }

    gsl_matrix *nscale = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix_set_all(nscale, 1.0);

    LALInferenceAddVariable(model->params, "psdscale", &nscale, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(model->params, "logdeltaf", &nscale_dflog, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(model->params, "psdBandsMin", &bands_min, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(model->params, "psdBandsMax", &bands_max, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);

    //Set up noise priors
    LALInferenceAddVariable(priorArgs,      "psddim",   &nscale_dim,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddMinMaxPrior(priorArgs,   "psdscale", &nscale_min,  &nscale_max,   LALINFERENCE_REAL8_t);
    LALInferenceAddMinMaxPrior(priorArgs,   "psdrange", &nscale_min,  &nscale_max,   LALINFERENCE_REAL8_t);
    LALInferenceAddVariable(priorArgs,      "psdsigma", &nscale_prior, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

    //Store meta data for noise model in proposal
    LALInferenceAddVariable(proposalArgs, "psdblock", &nscale_block, LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(proposalArgs, "psdbin",   &nscale_bin,   LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(proposalArgs, "psdsigma", &nscale_sigma, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);


  }//End of noise model initialization
  LALInferenceAddVariable(model->params, "psdScaleFlag", &nscale_flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  UINT4 psdGaussianPrior=1;
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFlatPrior");
  if(ppt)psdGaussianPrior=0;
  LALInferenceAddVariable(priorArgs, "psdGaussianPrior", &psdGaussianPrior,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);

  ppt = LALInferenceGetProcParamVal(commandLine, "--glitchFit");
  if (ppt)
      printf("WARNING: --glitchFit has been deprecated in favor of --glitch-fit\n");
  else
      ppt = LALInferenceGetProcParamVal(commandLine, "--glitch-fit");
  if (ppt)
      LALInferenceInitGlitchVariables(state, model->params);

  /* Handle, if present, requests for calibration parameters. */
  LALInferenceInitCalibrationVariables(state, model->params);

  //Only add waveform parameters to model if needed
  if(signal_flag)
  {
    /* The idea here is the following:
     * We call RegisterUniformVariable with startval=0 and meanigful min and max values.
     * That function will then take care of setting startval to a random value between min and max, or read a value from command line (with --parname VALUE).
     * The user can fix the param to a given value with --fix-parname --parname VALUE
     * */
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "chirpmass", zero, mcMin, mcMax, LALINFERENCE_PARAM_LINEAR);
    /* Check if running with symmetric (eta) or asymmetric (q) mass ratio.*/
    ppt=LALInferenceGetProcParamVal(commandLine,"--use-eta");
    if(ppt)
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "eta", zero, etaMin, etaMax, LALINFERENCE_PARAM_LINEAR);
    else
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "q", zero, qMin, qMax, LALINFERENCE_PARAM_LINEAR);


    if(!LALInferenceGetProcParamVal(commandLine,"--margphi") && !LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phase", zero, phiMin, phiMax, LALINFERENCE_PARAM_CIRCULAR);
    }

  /* Check for distance prior for use if the user samples in logdistance */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--distance-max"))) Dmax=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--distance-min"))) Dmin=atof(ppt->value);
  LALInferenceParamVaryType distanceVary = LALINFERENCE_PARAM_LINEAR;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--fix-distance")))
  {
    Dinitial=atof(ppt->value);
    distanceVary = LALINFERENCE_PARAM_FIXED;
  }
  
  LALInferenceRegisterUniformVariableREAL8(state, model->params, "logdistance", log(Dinitial), log(Dmin), log(Dmax), distanceVary) ;
  if(LALInferenceGetProcParamVal(commandLine,"--margdist")||LALInferenceGetProcParamVal(commandLine,"--margdist-comoving"))
  {
      /* If using margdist, remove the distance parameters and add the ranges into the model params as a way of passing them in */
      REAL8 a = log(Dmin), b=log(Dmax);
      LALInferenceAddMinMaxPrior(model->params, "logdistance", &a, &b, LALINFERENCE_REAL8_t);
      UINT4 margdist=1;
      int cosmology=0;
      LALInferenceAddVariable(model->params, "MARGDIST", &margdist, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
	  LALInferenceRemoveVariable(model->params, "logdistance");
      if(LALInferenceGetProcParamVal(commandLine,"--margdist-comoving"))
          cosmology=1;
      LALInferenceAddINT4Variable(model->params,"MARGDIST_COSMOLOGY",cosmology, LALINFERENCE_PARAM_FIXED);
  }
  
  LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psiMin, psiMax, LALINFERENCE_PARAM_LINEAR);
  LALInferenceRegisterUniformVariableREAL8(state, model->params, "costheta_jn", zero, costhetaJNmin, costhetaJNmax,LALINFERENCE_PARAM_LINEAR);

  /* Option to use the detector-aligned frame */
  if(!LALInferenceGetProcParamVal(commandLine,"--no-detector-frame") && nifo >1)
  {
        printf("Using detector-based sky frame\n");
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"t0",timeParam,timeMin,timeMax,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"cosalpha",0,-1,1,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"azimuth",0.0,0.0,LAL_TWOPI,LALINFERENCE_PARAM_CIRCULAR);
        /* add the time parameter then remove it so that the prior is set up properly */
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", timeParam, timeMin, timeMax,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRemoveVariable(model->params,"time");
        INT4 one=1;
        LALInferenceAddVariable(model->params,"SKY_FRAME",&one,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  }
  else
  {
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, raMin, raMax,LALINFERENCE_PARAM_CIRCULAR);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decMin, decMax, LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", timeParam, timeMin, timeMax,LALINFERENCE_PARAM_LINEAR);
   }

  /* If we are marginalising over the time, remove that variable from the model (having set the prior above) */
  /* Also set the prior in model->params, since Likelihood can't access the state! (ugly hack) */
  if(LALInferenceGetProcParamVal(commandLine,"--margtime") || LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
	  LALInferenceVariableItem *p=LALInferenceGetItem(state->priorArgs,"time_min");
	  LALInferenceAddVariable(model->params,"time_min",p->value,p->type,p->vary);
	  p=LALInferenceGetItem(state->priorArgs,"time_max");
	  LALInferenceAddVariable(model->params,"time_max",p->value,p->type,p->vary);
	  if (LALInferenceCheckVariable(model->params,"time")) LALInferenceRemoveVariable(model->params,"time");
      if (LALInferenceCheckVariable(model->params,"t0")) LALInferenceRemoveVariable(model->params,"t0");
	  if (LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
		  UINT4 margphi = 1;
		  LALInferenceAddVariable(model->params, "margtimephi", &margphi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
	  }
  }

    /* If requested by the user populate the testing GR or PPE model parameters */
  if (LALInferenceGetProcParamVal(commandLine,"--grtest-parameters") || LALInferenceGetProcParamVal(commandLine,"--ppe-parameters"))
  {
    LALInferenceInitNonGRParams(state, model);
  }
  if (LALInferenceGetProcParamVal(commandLine,"--grtest-parameters"))
  {
    ppt=LALInferenceGetProcParamVal(commandLine,"--grtest-parameters");
    if ((checkParamInList(ppt->value,"log10lambda_eff")) || (checkParamInList(ppt->value,"lambda_eff")) ){
        if ( (approx != IMRPhenomPv2) && (approx != IMRPhenomD) && (approx != SEOBNRv4_ROM) && (approx != IMRPhenomPv2_NRTidal) && (approx != IMRPhenomPv2_NRTidalv2) && (approx != IMRPhenomPv3HM) && (approx != IMRPhenomHM) )
            {
            XLALPrintWarning("LIV parameters not compatible with approximant %s. Can be used only with IMRPhenomPv2, IMRPhenomD, SEOBNRv4_ROM, IMRPhenomPv2_NRTidalv2, IMRPhenomPv3HM, IMRPhenomHM, and IMRPhenomPv2_NRTidal.\n",XLALSimInspiralGetStringFromApproximant(approx));
            }
    }
  }
  /* PPE parameters */

  ppt=LALInferenceGetProcParamVal(commandLine, "--TaylorF2ppE");
  if(approx==TaylorF2 && ppt){

    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppealpha",zero, -1000.0 , 1000.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppebeta", zero, -1000.0 , 1000.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppeuppera", zero, -3.0, 3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppeupperb", zero, -3.0, 3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppelowera", zero, -3.0, 2.0/3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppelowerb", zero, -4.5, 1.0, LALINFERENCE_PARAM_LINEAR);

  }

  if((!!LALInferenceGetProcParamVal(commandLine,"--tidalT") + !!LALInferenceGetProcParamVal(commandLine,"--tidal")
    + !!LALInferenceGetProcParamVal(commandLine,"--4PolyEOS") + !!LALInferenceGetProcParamVal(commandLine,"--4SpectralDecomp")
    + !!LALInferenceGetProcParamVal(commandLine,"--eos") + !!LALInferenceGetProcParamVal(commandLine,"--BinaryLove")) > 1 )
  {
      XLALPrintError("Error: cannot use more than one of --tidalT, --tidal, --4PolyEOS, --4SpectralDecomp, --eos and --BinaryLove.\n");
      XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  if(LALInferenceGetProcParamVal(commandLine,"--tidalT")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambdaT", zero, lambdaTMin, lambdaTMax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "dLambdaT", zero, dLambdaTMin, dLambdaTMax, LALINFERENCE_PARAM_LINEAR);
  // Pull in tidal parameters (lambda1,lambda2)
  } else if(LALInferenceGetProcParamVal(commandLine,"--tidal")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambda1", zero, lambda1Min, lambda1Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambda2", zero, lambda2Min, lambda2Max, LALINFERENCE_PARAM_LINEAR);
  // Pull in 4-piece polytrope parameters (logp1,gamma1,gamma2,gamma3)
  } else if(LALInferenceGetProcParamVal(commandLine,"--4PolyEOS")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "logp1", zero, logp1Min, logp1Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "gamma1", zero, gamma1Min, gamma1Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "gamma2", zero, gamma2Min, gamma2Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "gamma3", zero, gamma3Min, gamma3Max, LALINFERENCE_PARAM_LINEAR);
  // Pull in spectral decomposition parameters (SDgamma0,SDgamma1,SDgamma2,SDgamma3)
  } else  if(LALInferenceGetProcParamVal(commandLine,"--4SpectralDecomp")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "SDgamma0", zero, SDgamma0Min, SDgamma0Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "SDgamma1", zero, SDgamma1Min, SDgamma1Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "SDgamma2", zero, SDgamma2Min, SDgamma2Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "SDgamma3", zero, SDgamma3Min, SDgamma3Max, LALINFERENCE_PARAM_LINEAR);
  } else if((ppt=LALInferenceGetProcParamVal(commandLine,"--eos"))){
    LALSimNeutronStarEOS *eos=NULL;
    errnum=XLAL_SUCCESS;
    XLAL_TRY(eos=XLALSimNeutronStarEOSByName(ppt->value), errnum);
    if(errnum!=XLAL_SUCCESS)
        XLAL_ERROR_NULL(errnum,"%s: %s",__func__,XLALErrorString(errnum));
    
    XLAL_TRY(model->eos_fam = XLALCreateSimNeutronStarFamily(eos),errnum);
    if(errnum!=XLAL_SUCCESS)
        XLAL_ERROR_NULL(errnum,"%s: %s",__func__,XLALErrorString(errnum));
    if(!model->eos_fam) XLAL_ERROR_NULL(XLAL_EINVAL, "Unable to initialise EOS family");

  // Pull in symmetric tidal deformability (lambdaS) and the uniform variable used to marginlise over the BinaryLove fit uncertainty
  } else if((ppt=LALInferenceGetProcParamVal(commandLine,"--BinaryLove"))){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambdaS", zero, lambdaSMin, lambdaSMax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "BLuni", 0.5, 0.0, 1.0, LALINFERENCE_PARAM_LINEAR);
  }    
  
  

    if(LALInferenceGetProcParamVal(commandLine,"--dQuadMon12")){
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "dQuadMon1", zero, dQuadMonMin, dQuadMonMax, LALINFERENCE_PARAM_LINEAR); 
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "dQuadMon2", zero, dQuadMonMin, dQuadMonMax, LALINFERENCE_PARAM_LINEAR); 
	}else if(LALInferenceGetProcParamVal(commandLine,"--dQuadMonSA")){
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "dQuadMonS", zero, dQuadMonMin, dQuadMonMax, LALINFERENCE_PARAM_LINEAR); 
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "dQuadMonA", zero, dQuadMonMin, dQuadMonMax, LALINFERENCE_PARAM_LINEAR); 
		}

  if((!!LALInferenceGetProcParamVal(commandLine,"--dQuadMon12") + !!LALInferenceGetProcParamVal(commandLine,"--dQuadMonSA")) > 1 )
  {    
	XLALPrintError("Error: cannot use more than one of --dQuadMon12 and --dQuadMonSA.\n");
	XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinOrder");
  if(ppt) {
    spinO = atoi(ppt->value);
    LALInferenceAddVariable(model->params, "spinO", &spinO,
        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }
  LALSimInspiralTidalOrder tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine, "--tidalOrder");
  if(ppt) {
    tideO = atoi(ppt->value);
    LALInferenceAddVariable(model->params, "tideO", &tideO,
        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }

  /* LALInference uses a proprietary spin parameterization, the conversion from which
   * assumes the LALSimulations default frame */
  LALSimInspiralFrameAxis frameAxis = LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;

  model->LALpars = XLALCreateDict();
  XLALSimInspiralWaveformParamsInsertPNSpinOrder(model->LALpars,  spinO);
  XLALSimInspiralWaveformParamsInsertPNTidalOrder(model->LALpars, tideO);
  XLALSimInspiralWaveformParamsInsertFrameAxis(model->LALpars,frameAxis);
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--numreldata"))) {
    XLALSimInspiralWaveformParamsInsertNumRelData(model->LALpars, ppt->value);
    fprintf(stdout,"Template will use %s.\n",ppt->value);
  }
  if (LALInferenceGetProcParamVal(commandLine, "--liv"))
    XLALSimInspiralWaveformParamsInsertEnableLIV(model->LALpars, 1);

  /* Pass custom mode array list to waveform generator */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--modeList"))) {
    char *end_str;
    char substring[] = "-";
    // strtok_r will modify argument in place, since MCMC calls this again later we need to copy  to tmp string to avoid chaing command line
    char* tempstr = calloc(strlen(ppt->value)+1, sizeof(char));
    strcpy(tempstr, ppt->value);
    char *modes = strtok_r(tempstr, ",", &end_str);
    int  ii=0, jj=0;
    // These two will be lists containing integers with the pairs of (l,m) to use
    // e.g. modes_l[0]=2 and modes_m[0]=-1 means the (l,m)=(2,-1) will be passed on to LALSim
    // Eventually the code below loops over the lenght of these array and pass all modes to LALSim
    int *modes_l=NULL;
    int *modes_m=NULL;
    // Counters used to realloc
    int m_size=1;
    int l_size=1;
    fprintf(stdout, "Template will use a custom mode array.\n");
    while (modes != NULL) {
        int k = 0;
	// reads in modes with negative m's first
        if (strstr(modes, substring) != NULL) {
            char *end_token;
            char *token = strtok_r(modes, "-", &end_token);
            while (token != NULL) {
                if ( k == 0 ) {
	            modes_l=realloc(modes_l,sizeof(*modes_l)*l_size);
		    l_size+=1;
		    modes_l[ii++]= atoi(token);
                } else {
		    modes_m=realloc(modes_m,sizeof(*modes_m)*m_size);
                    m_size+=1;
		    modes_m[jj++]= -1 * atoi(token);
                }
                k += 1;
                token = strtok_r(NULL, "-", &end_token);
            }
        } else {
	    // now read modes with positive m. These are read in as integer, then module 10 is used to extract the m and l
	    // FIXME this will break if l or m ever get >10.
            int value = atoi(modes);
            for (k=2; k>=0; k--) {
                if (k == 2) {
		    modes_m=realloc(modes_m,sizeof(*modes_m)*m_size);
                    m_size+=1;
		    modes_m[jj++] = value % 10;
                } else if (k == 1) {
		    modes_l=realloc(modes_l,sizeof(*modes_l)*l_size);
                    l_size+=1;
		    modes_l[ii++]= value % 10;
                }
                value /= 10;
            }
        }
        modes = strtok_r(NULL, ",", &end_str);
    }

    // Create empty modearray and fills it with the modes
    LALValue *ModeArray = XLALSimInspiralCreateModeArray();
    for (int idx_i=0; idx_i<ii; idx_i+=1){
       fprintf(stdout, "Template will use the l=%d m=%d mode \n", modes_l[idx_i], modes_m[idx_i]);
       XLALSimInspiralModeArrayActivateMode(ModeArray, modes_l[idx_i], modes_m[idx_i]);
    }
    XLALSimInspiralWaveformParamsInsertModeArray(model->LALpars, ModeArray);
    XLALDestroyValue(ModeArray);
  }

  /* Pass threshold for PhenomXHM Multibanding */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--phenomXHMMband"))) {
      REAL8 threshold = atof(ppt->value);
      XLALSimInspiralWaveformParamsInsertPhenomXHMThresholdMband(model->LALpars, threshold);
      fprintf(stdout,"Template will use phenomXHMMband %s.\n",ppt->value);
  }

  /* Pass threshold for PhenomXPHM Multibanding of the Euler angles. */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--phenomXPHMMband"))) {
      REAL8 threshold = atof(ppt->value);
      XLALSimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(model->LALpars, threshold);
      fprintf(stdout,"Template will use phenomXPHMMband %s.\n",ppt->value);
  }

  /* Version for the final spin model to use in IMRPhenomXP/IMRPhenomXPHM. */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--phenomXPFinalSpinMod"))) {
      REAL8 afversion = atoi(ppt->value);
      XLALSimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(model->LALpars, afversion);
      fprintf(stdout,"Template will use PhenomXPFinalSpinMod %s.\n",ppt->value);
  }

  /* Version of the Euler angles in the twisting-up of IMRPhenomXP/IMRPhenomXPHM. */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--phenomXPrecVersion"))) {
      REAL8 eulerangleversion = atoi(ppt->value);
      XLALSimInspiralWaveformParamsInsertPhenomXPrecVersion(model->LALpars, eulerangleversion);
      fprintf(stdout,"Template will use PhenomXPrecVersion %s.\n",ppt->value);
  }




  fprintf(stdout,"\n\n---\t\t ---\n");
  LALInferenceInitSpinVariables(state, model);
  LALInferenceCheckApproximantNeeds(state,approx);

  if (injTable)
     print_flags_orders_warning(injTable,commandLine);

     /* Print info about orders and waveflags used for templates */

     fprintf(stdout,"Templates will run using Approximant %i (%s), phase order %i, amp order %i, spin order %i tidal order %i in the %s domain.\n",approx,XLALSimInspiralGetStringFromApproximant(approx),PhaseOrder,AmpOrder,(int) spinO, (int) tideO, model->domain==LAL_SIM_DOMAIN_TIME?"time":"frequency");
     fprintf(stdout,"---\t\t ---\n\n");
  }//end of signal only flag
  else
  {
    /* Print info about orders and waveflags used for templates */
    fprintf(stdout,"\n\n------\n");
    fprintf(stdout,"Noise only run\n");
    fprintf(stdout,"------\n\n");
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

  model->freqhs = XLALCalloc(nifo, sizeof(COMPLEX16FrequencySeries *));
  for (i=0; i<nifo; i++)
      model->freqhs[i] = XLALCreateCOMPLEX16FrequencySeries("freqh",
                                                            &(state->data->freqData->epoch),
                                                            0.0,
                                                            model->deltaF,
                                                            &lalDimensionlessUnit,
                                                            state->data->freqData->data->length);

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  /* Choose proper template */
  model->templt = LALInferenceInitCBCTemplate(state);

  /* Use same window and FFT plans on model as data */
  model->window = state->data->window;
  model->timeToFreqFFTPlan = state->data->timeToFreqFFTPlan;
  model->freqToTimeFFTPlan = state->data->freqToTimeFFTPlan;

  /* Initialize waveform cache */
  model->waveformCache = XLALCreateSimInspiralWaveformCache();

  return(model);
}



/* Setup the variable for the evidence calculation test for review */
/* 5-sigma ranges for analytic likeliood function */
/* https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/LALInferenceReviewAnalyticGaussianLikelihood */
LALInferenceModel *LALInferenceInitModelReviewEvidence(LALInferenceRunState *state)
{
    LALInferenceIFOData *dataPtr;
    INT4 nifo=0;
    ProcessParamsTable *commandLine=state->commandLine;
    ProcessParamsTable *ppt=NULL;
    char **strings=NULL;
    char *pinned_params=NULL;
    UINT4 N=0,i;
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
            pinned_params=ppt->value;
            LALInferenceVariables tempParams;
            memset(&tempParams,0,sizeof(tempParams));
            LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    }

    LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
    model->params = XLALCalloc(1, sizeof(LALInferenceVariables));

    dataPtr = state->data;
    while (dataPtr != NULL) {
      nifo++;
      dataPtr = dataPtr->next;
    }

    /* Create arrays for holding single-IFO likelihoods, etc. */
    model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
    model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

	i=0;

  /* Parameter bounds at ±5 sigma */
  fprintf(stdout,"Setting up priors\n");
  LALInferenceParamVaryType type=LALINFERENCE_PARAM_LINEAR;
  for(i=0;i<15;i++)
  {
    REAL8 min = LALInferenceAnalyticMeansCBC[i] - 5.0/scaling[i]*sqrt(CM[i][i]);
    REAL8 max = LALInferenceAnalyticMeansCBC[i] + 5.0/scaling[i]*sqrt(CM[i][i]);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, LALInferenceAnalyticNamesCBC[i], LALInferenceAnalyticMeansCBC[i], min, max, type );
    fprintf(stdout,"%s: %e - %e\n",LALInferenceAnalyticNamesCBC[i],min,max);
  }

	return(model);
}


LALInferenceModel *LALInferenceInitModelReviewEvidence_bimod(LALInferenceRunState *state)
{
  LALInferenceIFOData *dataPtr;
  INT4 nifo=0;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  char **strings=NULL;
  char *pinned_params=NULL;
  UINT4 N=0,i;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
  }

  LALInferenceModel *model = XLALCalloc(1, sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));

  dataPtr = state->data;
  while (dataPtr != NULL) {
    nifo++;
    dataPtr = dataPtr->next;
  }

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  i=0;

  /* Parameter bounds ± 5 sigma from 2 modes separated by 8 sigma */
  LALInferenceParamVaryType type=LALINFERENCE_PARAM_LINEAR;
  fprintf(stdout,"Setting up priors\n");
  for(i=0;i<15;i++)
  {
    REAL8 min = LALInferenceAnalyticMeansCBC[i] - (4.0+5.0)/scaling[i]*sqrt(CM[i][i]);
    REAL8 max = LALInferenceAnalyticMeansCBC[i] + (4.0+5.0)/scaling[i]*sqrt(CM[i][i]);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, LALInferenceAnalyticNamesCBC[i], LALInferenceAnalyticMeansCBC[i], min, max, type );
    fprintf(stdout,"%s: %e - %e\n",LALInferenceAnalyticNamesCBC[i],min,max);
  }

  return(model);
}

LALInferenceModel *LALInferenceInitModelReviewEvidence_banana(LALInferenceRunState *state)
{
  LALInferenceIFOData *dataPtr;
  INT4 nifo = 0;
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

  dataPtr = state->data;
  while (dataPtr != NULL) {
      nifo++;
      dataPtr = dataPtr->next;
  }

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  i=0;

  struct varSettings {const char *name; REAL8 val, min, max;};

  struct varSettings setup[]=
  {
    {.name="time", .val=0.0, .min=-2., .max=2.},
    {.name="mass1", .val=16., .min=14., .max=18.},
    {.name="mass2", .val=7., .min=5., .max=9.},
    {.name="logdistance", .val=log(50.), .min=log(45.), .max=log(55.)},
    {.name="costheta_jn", .val=cos(LAL_PI/2.), .min=cos(3.570796327), .max=cos(-0.429203673)},
    {.name="phase", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="polarisation", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="rightascension", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="declination", .val=0., .min=-2., .max=2.},
    {.name="a_spin1", .val=0.5, .min=-1.5, .max=2.5},
    {.name="a_spin2", .val=0.5, .min=-1.5, .max=2.5},
    {.name="tilt_spin1", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="tilt_spin2", .val=LAL_PI/2., .min=-0.429203673, .max=3.570796327},
    {.name="phi12", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
    {.name="phi_jl", .val=LAL_PI, .min=1.141592654, .max=5.141592654},
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

static void print_flags_orders_warning(SimInspiralTable *injt, ProcessParamsTable *commline){

    /* If lalDebugLevel > 0, print information about:
     *
     * - Eventual injection/template mismatch on phase and amplitude orders, as well as on waveFlags
     * - Those fiels being set only for injection or template
     *
     **/
    XLALPrintWarning("\n");
    LALPNOrder PhaseOrder=-1;
    LALPNOrder AmpOrder=-1;
    LALSimInspiralSpinOrder default_spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;
    LALSimInspiralTidalOrder default_tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;
    Approximant approx=NumApproximants;
    ProcessParamsTable *ppt=NULL;
    ProcessParamsTable *ppt_order=NULL;
    int errnum;
    ppt=LALInferenceGetProcParamVal(commline,"--approximant");
    if(ppt){
        approx=XLALGetApproximantFromString(ppt->value);
        ppt=LALInferenceGetProcParamVal(commline,"--order");
        if(ppt) PhaseOrder = XLALGetOrderFromString(ppt->value);
    }
    ppt=LALInferenceGetProcParamVal(commline,"--approx");
    if(ppt){
       approx=XLALGetApproximantFromString(ppt->value);
       XLAL_TRY(PhaseOrder = XLALGetOrderFromString(ppt->value),errnum);
       if( (int) PhaseOrder == XLAL_FAILURE || errnum) {
          XLALPrintWarning("WARNING: No phase order given.  Using maximum available order for the template.\n");
          PhaseOrder=-1;
        }
     }
     /* check approximant is given */
    if (approx==NumApproximants){
        approx=XLALGetApproximantFromString(injt->waveform);
        XLALPrintWarning("WARNING: You did not provide an approximant for the templates. Using value in injtable (%s), which might not what you want!\n",XLALSimInspiralGetStringFromApproximant(approx));
     }

    /* check inj/rec amporder */
    ppt=LALInferenceGetProcParamVal(commline,"--amporder");
    if(ppt) AmpOrder=atoi(ppt->value);
    ppt=LALInferenceGetProcParamVal(commline,"--ampOrder");
    if(ppt) AmpOrder = XLALGetOrderFromString(ppt->value);
    if(AmpOrder!=(LALPNOrder)injt->amp_order)
       XLALPrintWarning("WARNING: Injection specified amplitude order %i. Template will use  %i\n",
               injt->amp_order,AmpOrder);

    /* check inj/rec phase order */
    if(PhaseOrder!=(LALPNOrder)XLALGetOrderFromString(injt->waveform))
        XLALPrintWarning("WARNING: Injection specified phase order %i. Template will use %i\n",\
             XLALGetOrderFromString(injt->waveform),PhaseOrder);

    /* check inj/rec spinflag */
    ppt=LALInferenceGetProcParamVal(commline, "--spinOrder");
    ppt_order=LALInferenceGetProcParamVal(commline, "--inj-spinOrder");
    if (ppt && ppt_order){
       if (!(atoi(ppt->value)== atoi(ppt_order->value)))
            XLALPrintWarning("WARNING: Set different spin orders for injection (%i ) and template (%i) \n",atoi(ppt_order->value),atoi(ppt->value));
    }
    else if (ppt || ppt_order){
        if (ppt)
            XLALPrintWarning("WARNING: You set the spin order only for the template (%i). Injection will use default value (%i). You can change that with --inj-spinOrder. \n",atoi(ppt->value),default_spinO);
        else
            XLALPrintWarning("WARNING: You set the spin order only for the injection (%i). Template will use default value (%i). You can change that with --spinOrder. \n",atoi(ppt_order->value),default_spinO);     }
    else
        XLALPrintWarning("WARNING: You did not set the spin order. Injection and template will use default values (%i). You change that using --inj-spinOrder (set injection value) and --spinOrder (set template value).\n",default_spinO);
    /* check inj/rec tidal flag */
    ppt=LALInferenceGetProcParamVal(commline, "--tidalOrder");
    ppt_order=LALInferenceGetProcParamVal(commline, "--inj-tidalOrder");
    if (ppt && ppt_order){
        if (!(atoi(ppt->value)==atoi( ppt_order->value)))
            XLALPrintWarning("WARNING: Set different tidal orders for injection (%i ) and template (%i) \n",atoi(ppt_order->value),atoi(ppt->value));
    }
    else if (ppt || ppt_order){
        if (ppt)
            XLALPrintWarning("WARNING: You set the tidal order only for the template (%d). Injection will use default value (%i). You can change that with --inj-tidalOrder. \n",atoi(ppt->value),default_tideO);
        else
            XLALPrintWarning("WARNING: You set the tidal order only for the injection (%i). Template will use default value (%i). You can  change that with --tidalOrder\n",atoi(ppt_order->value),default_tideO);
        }
    else
       XLALPrintWarning("WARNING: You did not set the tidal order. Injection and template will use default values (%i). You change that using --inj-tidalOrder (set injection value) and --tidalOrder (set template value).\n",default_tideO);
    return;
}

void LALInferenceCheckOptionsConsistency(ProcessParamsTable *commandLine)

{ /*
  Go through options and check for possible errors and inconsistencies (e.g. seglen < 0 )

  */

  ProcessParamsTable *ppt=NULL,*ppt2=NULL;
  REAL8 tmp=0.0;
  INT4 itmp=0;

  ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  if (ppt)
    return;

  // Check PSDlength > 0 if specified
  ppt=LALInferenceGetProcParamVal(commandLine,"--psdlength");
  if (ppt) {
      tmp=atof(ppt->value);
      if (tmp<0.0){
        fprintf(stderr,"ERROR: PSD length must be positive. Exiting...\n");
        exit(1);
      }
  }

  // Check seglen > 0
  REAL8 seglen=0.;
  ppt=LALInferenceGetProcParamVal(commandLine,"--seglen");
  if (!ppt){
    XLALPrintError("Must provide segment length with --seglen. Exiting...");
    exit(1);
  }
  else seglen=atof(ppt->value);

  tmp=atof(ppt->value);
  if (tmp<0.0){
    fprintf(stderr,"ERROR: seglen must be positive. Exiting...\n");
    exit(1);
  }
  REAL8 timeSkipStart=0.;
  REAL8 timeSkipEnd=0.;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--time-pad-start")))
     timeSkipStart=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--time-pad-end")))
     timeSkipEnd=atof(ppt->value);
  if(timeSkipStart+timeSkipEnd > seglen)
  {
    fprintf(stderr,"ERROR: --time-pad-start + --time-pad-end is greater than --seglen!");
    exit(1);
  }

  /* Flags consistency */
  ppt=LALInferenceGetProcParamVal(commandLine,"--disable-spin");
  ppt2=LALInferenceGetProcParamVal(commandLine,"--noSpin");
  if (ppt || ppt2){
    ppt2=LALInferenceGetProcParamVal(commandLine,"--spinO");
    if (ppt2){
      itmp=atoi(ppt2->value);
      if (itmp>0 || itmp==-1)
        XLALPrintWarning("--spinO > 0 or -1 will be ignored due to --disable-spin. If you want to include spin terms in the template, remove --disable-spin\n");
        exit(1);
      }
    if (!ppt2){
      XLALPrintWarning("--spinO defaulted to -1. This will be ignored due to --disable-spin. If you want to include spin terms in the template, remove --disable-spin\n");
      }
    ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
    ppt2=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
    if (ppt|| ppt2){
      fprintf(stderr,"--aligned-spin and --disable-spin are incompatible options. Exiting\n");
      exit(1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--singleSpin");
    if(ppt){
      fprintf(stderr,"--singleSpin and --disable-spin are incompatible options. Exiting\n");
      exit(1);
    }
  }

  /* lalinference_nest only checks */
  // Check live points
  ppt=LALInferenceGetProcParamVal(commandLine,"--nlive");
  if (ppt){
    itmp=atoi(ppt->value);
    if (itmp<0){
      fprintf(stderr,"ERROR: nlive must be positive. Exiting...\n");
      exit(1);
    }
    if (itmp<100){
      XLALPrintWarning("WARNING: Using %d live points. This is very little and may lead to unreliable results. Consider increasing.\n",itmp);
    }
    if (itmp>5000){
      XLALPrintWarning("WARNING: Using %d live points. This is a very large number and may lead to very long runs. Consider decreasing.\n",itmp);
    }
  }
  // Check nmcmc points
  ppt=LALInferenceGetProcParamVal(commandLine,"--nmcmc");
  if (ppt){
    itmp=atoi(ppt->value);
    if (itmp<0){
      fprintf(stderr,"ERROR: nmcmc must be positive (or omitted). Exiting...\n");
      exit(1);
    }
    if (itmp<100){
      XLALPrintWarning("WARNING: Using %d nmcmc. This is very little and may lead to unreliable results. Consider increasing.\n",itmp);
    }
  }

  /* Ensure that the user is not trying to marginalise the likelihood
     in an inconsistent way */
  if (LALInferenceGetProcParamVal(commandLine, "--margtime") && LALInferenceGetProcParamVal(commandLine, "--margphi")) {
    fprintf(stderr, "ERROR: trying to separately marginalise in time and phase.  Use '--margtimephi' instead");
    exit(1);
  }

  if (LALInferenceGetProcParamVal(commandLine, "--margtimephi") && LALInferenceGetProcParamVal(commandLine, "--margtime")) {
    fprintf(stderr, "ERROR: cannot marginalise in time and phi and separately time.  Pick either '--margtimephi' OR '--margtime'");
    exit(1);
  }

  if (LALInferenceGetProcParamVal(commandLine, "--margtimephi") && LALInferenceGetProcParamVal(commandLine, "--margphi")) {
    fprintf(stderr, "ERROR: cannot marginalise in time and phi and separately in phi.  Pick either '--margtimephi' OR '--margtime'");
    exit(1);
  }

  /* Check for small sample rates when margtime-ing. */
  if (LALInferenceGetProcParamVal(commandLine, "--margtime") || LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
    ppt = LALInferenceGetProcParamVal(commandLine, "--srate");
    if (ppt) {
      int srate = atoi(ppt->value);

      if (srate < 4096) {
	XLALPrintWarning("WARNING: you have chosen to marginalise in time with a sample rate of %d, but this typically gives incorrect results for CBCs; use at least 4096 Hz to be safe", srate);
      }
    }
  }

  return;
}

void LALInferenceInitSpinVariables(LALInferenceRunState *state, LALInferenceModel *model){


  LALStatus status;
  memset(&status,0,sizeof(status));

  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;

  Approximant approx= *(Approximant*) LALInferenceGetVariable(model->params, "LAL_APPROXIMANT");

  REAL8 a1min=0.0,a1max=1.0;
  REAL8 a2min=0.0,a2max=1.0;
  REAL8 tilt1min=0.0,tilt1max=LAL_PI;
  REAL8 tilt2min=0.0,tilt2max=LAL_PI;
  REAL8 phi12min=0.0,phi12max=LAL_TWOPI;
  REAL8 phiJLmin=0.0,phiJLmax=LAL_TWOPI;

  /* Default to precessing spins */
  UINT4 spinAligned=0;
  UINT4 singleSpin=0;
  UINT4 noSpin=0;

  /* Let's first check that the user asked, then we check the approximant can make it happen */
  /* Check for spin disabled */
  ppt=LALInferenceGetProcParamVal(commandLine, "--noSpin");
  if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--disable-spin");
  if (ppt) noSpin=1;

  /* Check for aligned spin */
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinAligned");
  if (!ppt) ppt=LALInferenceGetProcParamVal(commandLine,"--aligned-spin");
  if(ppt)
    spinAligned=1;

  /* Check for single spin */
  ppt=LALInferenceGetProcParamVal(commandLine,"--singleSpin");
  if(ppt){
    singleSpin=1;
  }

  SpinSupport spin_support=XLALSimInspiralGetSpinSupportFromApproximant(approx);

  /* Now check what the approx can do and eventually change user's choices to comply.
   * Also change the reference frame -- For the moment use default as the corresponding patch to LALSimulation has not been finished yet */
  if (spin_support==LAL_SIM_INSPIRAL_SPINLESS)
    noSpin=1;
  else if (spin_support==LAL_SIM_INSPIRAL_SINGLESPIN)
    singleSpin=1;
  else if (spin_support==LAL_SIM_INSPIRAL_ALIGNEDSPIN){
    spinAligned=1;
  }

  if (spinAligned){
  /* If spin aligned the magnitude is in the range [-1,1] */
    a1min=-1.0;
    a2min=-1.0;
  }

  /* IMRPhenomP only supports spins up to 0.9 and q > 1/10. Set prior consequently*/
  if (approx==IMRPhenomP && (a1max>0.9 || a2max>0.9)){
    a1max=0.9;
    a2max=0.9;
    if (spinAligned){
      a1min=-0.9;
      a2min=-0.9;
      }
    XLALPrintWarning("WARNING: IMRPhenomP only supports spin magnitude up to 0.9. Changing the a1max=a2max=0.9.\n");
  }

  /* IMRPhenomPv2 preliminary only supports spins up to 0.99 and q > 1/20. Set prior consequently*/
  if (approx==IMRPhenomPv2 && (a1max>0.99 || a2max>0.99)){
    a1max=0.99;
    a2max=0.99;
    /* If spin-aligned, IMRPhenomPv2 should reduce to IMRPhenomD with a lower spin magnitude of -1.0*/
    if (spinAligned){
      a1min=-1.0;
      a2min=-1.0;
      }
    XLALPrintWarning("WARNING: IMRPhenomPv2 preliminary only supports spin magnitude up to 0.99. Changing the a1max=a2max=0.99.\n");
  }

  /* Start with parameters that are free (or pinned if the user wants so). The if...else below may force some of them to be fixed or ignore some of them, depending on the spin configuration*/
  LALInferenceParamVaryType tilt1Vary = LALINFERENCE_PARAM_LINEAR;
  LALInferenceParamVaryType tilt2Vary = LALINFERENCE_PARAM_LINEAR;
  LALInferenceParamVaryType phi12Vary = LALINFERENCE_PARAM_CIRCULAR;
  LALInferenceParamVaryType spin1Vary = LALINFERENCE_PARAM_LINEAR;
  LALInferenceParamVaryType spin2Vary = LALINFERENCE_PARAM_LINEAR;
  LALInferenceParamVaryType phiJLVary = LALINFERENCE_PARAM_CIRCULAR;

  /* Add parameters depending on the values of noSpin, singleSpin and spinAligned
   * noSpin -> add nothing
   * spinAligned -> add a_spin1 and a_spin2 (if the approximant is spin aligned only use the names spin1,spin1 instead)
   * singleSpin -> add a_spin1, tilt_spin1, phi_JL
   * singleSpin+spinAligned -> add a_spin1
   * otherwise -> add everything */
   /* Note: To get spin aligned is sufficient to not add the spin angles because LALInferenceTemplate will default to aligned spin if it doesn't find spin angle params in model->params. */
  if (!noSpin){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "a_spin1", 0.0, a1min, a1max,spin1Vary);
    if (!singleSpin)
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "a_spin2", 0.0, a2min, a2max,spin2Vary);
    if (!spinAligned){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phi_jl", 0.0, phiJLmin,  phiJLmax, phiJLVary);
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "tilt_spin1", 0.0, tilt1min,tilt1max,tilt1Vary);
      if (!singleSpin){
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "tilt_spin2", 0.0, tilt2min,tilt2max,tilt2Vary);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "phi12", 0.0, phi12min,phi12max,phi12Vary);
      }
    }
  }

  /* Print to stdout what will be used */
  if (noSpin)
    fprintf(stdout,"Templates will run without spins\n");
  else{
    if (spinAligned && singleSpin)
      fprintf(stdout,"Templates will run with spin 1 aligned to L \n");
    if (spinAligned && !singleSpin)
      fprintf(stdout,"Templates will run with spins aligned to L \n");
    if (!spinAligned && singleSpin)
      fprintf(stdout,"Templates will run with precessing spin 1 \n");
    if (!spinAligned && !singleSpin)
      fprintf(stdout,"Templates will run with precessing spins \n");
  }
}

void LALInferenceInitMassVariables(LALInferenceRunState *state){

  LALStatus status;
  memset(&status,0,sizeof(status));

  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;

  REAL8 m1_min=1.0,m2_min=1.0;
  REAL8 m1_max=100.0,m2_max=100.0;
  REAL8 MTotMax=200.0;
  REAL8 MTotMin=2.0;

  /* Over-ride component masses */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--comp-min")))
  {
    m1_min=m2_min=atof(ppt->value);
  }

  if((ppt=LALInferenceGetProcParamVal(commandLine,"--comp-max")))
  {
    m1_max=m2_max=atof(ppt->value);
  }

  /* optional limits on individual masses */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mass1-min")))
  {
    m1_min=atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mass1-max")))
  {
    m1_max = atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mass2-min")))
  {
    m2_min=atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mass2-max")))
  {
    m2_max=atof(ppt->value);
  }

  /* Set the total mass bounds based on m1,m2 */
  MTotMin=m1_min + m2_min;
  MTotMax=m1_max + m2_max;

  /* Over-ride Mtotal bounds if requested */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mtotal-max")))
  {
    MTotMax=atof(ppt->value);
  }
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--mtotal-min")))
  {
    MTotMin=atof(ppt->value);
  }

  LALInferenceAddREAL8Variable(priorArgs,"mass1_min",m1_min,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddREAL8Variable(priorArgs,"mass1_max",m1_max,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddREAL8Variable(priorArgs,"mass2_min",m2_min,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddREAL8Variable(priorArgs,"mass2_max",m2_max,LALINFERENCE_PARAM_FIXED);

  LALInferenceAddVariable(priorArgs,"MTotMax",&MTotMax,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs,"MTotMin",&MTotMin,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);

  return;

}

void LALInferenceCheckApproximantNeeds(LALInferenceRunState *state,Approximant approx){

  REAL8 min,max;
  REAL8 a1min,a1max,a2min,a2max;
  UINT4 q=0;

  if (LALInferenceCheckVariable(state->priorArgs,"q_min")){
    LALInferenceGetMinMaxPrior(state->priorArgs, "q", &min, &max);
    q=1;
  }
  else if (LALInferenceCheckVariable(state->priorArgs,"eta_min"))
    LALInferenceGetMinMaxPrior(state->priorArgs, "eta", &min, &max);

  /* IMRPhenomP only supports q > 1/10. Set prior consequently  */
  if (q==1 && approx==IMRPhenomP && min<1./10.){
    min=1.0/10.;
    LALInferenceRemoveVariable(state->priorArgs,"q_min");
    LALInferenceAddVariable(state->priorArgs,"q_min",&min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"WARNING: IMRPhenomP only supports mass ratios up to 10 ( suggested max: 4). Changing the min prior for q to 1/10\n");
  }
  if (q==0 && approx==IMRPhenomP && min<0.08264462810){
     min=0.08264462810;  //(that is eta for a 1-10 system)
     LALInferenceRemoveVariable(state->priorArgs,"eta_min");
     LALInferenceAddVariable(state->priorArgs,"eta_min",&min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
     fprintf(stdout,"WARNING: IMRPhenomP only supports mass ratios up to 10 ( suggested max: 4). Changing the min prior for eta to 0.083\n");
  }

  /* In the IMRPhenomD review a region with possible issues  (q < 1/10 for maximal spins) was identified */
  /* This is not cause to restrict the waveform to exclude this region, but caution should be exercised */
  if (approx==IMRPhenomD){
    if (LALInferenceCheckMinMaxPrior(state->priorArgs,"a_spin1") && LALInferenceCheckMinMaxPrior(state->priorArgs,"a_spin2")){
       LALInferenceGetMinMaxPrior(state->priorArgs, "a_spin1", &a1min, &a1max);
       LALInferenceGetMinMaxPrior(state->priorArgs, "a_spin2", &a2min, &a2max);
       if (q==1 && min<1./10 && ((a1max==1.0 || a2max==1.0) || (a1min==-1.0 || a2min==-1.0))){
          fprintf(stdout,"WARNING: Based on the allowed prior volume (q < 1/10 for maximal spins), please consult\n");
          fprintf(stdout,"IMRPhenomD review wiki for further discussion on reliable parameter spaces\n");
          fprintf(stdout,"https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/PhenD_LargeNegativeSpins\n");
       }
       if (q==0 && min<0.082644628 && ((a1max==1.0 || a2max==1.0) || (a1min==-1.0 || a2min==-1.0))){
          fprintf(stdout,"WARNING: Based on the allowed prior volume (q < 1/10 for maximal spins), please consult\n");
          fprintf(stdout,"IMRPhenomD review wiki for further discussion on reliable parameter spaces\n");
          fprintf(stdout,"https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/PhenD_LargeNegativeSpins\n");
       }
    }
  }


  /* IMRPhenomPv2 only supports q > 1/20. Set prior consequently  */
  if (q==1 && approx==IMRPhenomPv2 && min<1./20.){
    min=1.0/20.;
    LALInferenceRemoveVariable(state->priorArgs,"q_min");
    LALInferenceAddVariable(state->priorArgs,"q_min",&min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"WARNING: IMRPhenomPv2 preliminary only supports mass ratios up to 20 ( suggested max: 18). Changing the min prior for q to 1/20\n");
  }
  if (q==0 && approx==IMRPhenomPv2 && min<0.04535147392){
     min=0.04535147392;  //(that is eta for a 1-20 system)
     LALInferenceRemoveVariable(state->priorArgs,"eta_min");
     LALInferenceAddVariable(state->priorArgs,"eta_min",&min,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_FIXED);
     fprintf(stdout,"WARNING: IMRPhenomPv2 preliminary only supports mass ratios up to 20 ( suggested max: 20). Changing the min prior for eta to 0.045\n");
  }

  (void) max;
  return;
}

/*******************************************************************
 * LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model)
 * Function to initialise either the TaylorF2Test of SpinTaylorT4Test waveform models
 * or the PPE waveform model
 *******************************************************************/
static void LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model)
{
    ProcessParamsTable *commandLine = state->commandLine;
    ProcessParamsTable *ppt=NULL;
    ProcessParamsTable *ppta=NULL;
    ProcessParamsTable *pptb=NULL;
    /* check that the user does not request both a TaylorF2Test and a PPE waveform model */
    if (LALInferenceGetProcParamVal(commandLine,"--grtest-parameters") && LALInferenceGetProcParamVal(commandLine,"--ppe-parameters"))
    {
        fprintf(stderr,"--grtest-parameters and --ppe-parameters are not simultaneously supported. Please choose one. Aborting\n");
        exit(-1);
    }


    ppt=LALInferenceGetProcParamVal(commandLine,"--grtest-parameters");
    if (ppt)
    {
        REAL8 dchi_max=1.;
        REAL8 dchi_min=-1.;
        REAL8 dxi_max=1.;
        REAL8 dxi_min=-1.;
        REAL8 dalpha_max=1.;
        REAL8 dalpha_min=-1.;
        REAL8 dbeta_max=1.;
        REAL8 dbeta_min=-1.;
        REAL8 dsigma_max=1.;
        REAL8 dsigma_min=-1.;
        REAL8 tmpVal=0.0;
        if ((pptb=LALInferenceGetProcParamVal(commandLine,"--LIV_A_sign"))) {
          REAL8 LIV_A_sign;
          LIV_A_sign = atof(pptb->value);
          if ((LIV_A_sign != -1) && (LIV_A_sign != 1)) {
            XLALPrintError("LIV_A_sign can only take either +1 or -1.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
          }
          else LALInferenceAddVariable(model->params,"LIV_A_sign", &LIV_A_sign, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
        if ((ppta=LALInferenceGetProcParamVal(commandLine,"--nonGR_alpha"))) {
          REAL8 nonGR_alpha;
          nonGR_alpha = atof(ppta->value);
          if (nonGR_alpha==2.0) {
            XLALPrintError("nonGR_alpha=2 is degenerate with p^2c^2 term in dispersion relation and cannot be tested. Please choose a different value.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
          }
          else LALInferenceAddVariable(model->params,"nonGR_alpha", &nonGR_alpha, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        }
	/* Relative shifts for inspiral phase PN coefficients (absolute value for dchi1) */
        if (checkParamInList(ppt->value,"dchi0")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi0", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi1", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi2", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi3", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi4", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi5", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi5l")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi5l", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi6")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi6", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi6l")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi6l", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi7")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi7", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for pre-merger phase coefficients (PhenomC/P) */
        if (checkParamInList(ppt->value,"dxi1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi1", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi2", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi3", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi4", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi5", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi6")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi6", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for merger-ringdown phase coefficients  (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dalpha1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha1", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha2", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha3", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha4", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha5", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for phenomenological inspiral phase coefficients (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dsigma1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma1", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma2", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma3", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma4", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for intermediate phase coefficients (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dbeta1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta1", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dbeta2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta2", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dbeta3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta3", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"lambda_eff")) {
          if (ppta==NULL) {
            XLALPrintError("A value for nonGR_alpha has to be passed with lambda_eff.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
           }
          else if (pptb==NULL) {
            XLALPrintError("A value of +1 or -1 for LIV_A_sign has to be passed with lambda_eff, respectively determing if +A or -A in the MDR is being tested.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
          }
          else {
            REAL8 lambda_eff_min = 1E-6;
            REAL8 lambda_eff_max = 1E13;
            LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambda_eff", 1E11, lambda_eff_min, lambda_eff_max, LALINFERENCE_PARAM_LINEAR);
          }
        }
        if (checkParamInList(ppt->value,"log10lambda_eff")) {
          if (ppta==NULL) {
            XLALPrintError("A value for nonGR_alpha has to be passed with log10lambda_eff. A value for LIV_A_sign also has to be passed!\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
           }
          else if (pptb==NULL) {
            XLALPrintError("A value of +1 or -1 for LIV_A_sign has to be passed with log10lambda_eff, respectively determing if +A or -A in the MDR is being tested.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
          }
          else {
            REAL8 log10lambda_eff_min = -6.0;
            REAL8 log10lambda_eff_max = 13.0;
            LALInferenceRegisterUniformVariableREAL8(state, model->params, "log10lambda_eff", 3.0, log10lambda_eff_min, log10lambda_eff_max, LALINFERENCE_PARAM_LINEAR);
          }
      }

  }

    ppt=LALInferenceGetProcParamVal(commandLine,"--ppe-parameters");
    if (ppt)
    {
        /* amplitude parameters */
        REAL8 appe_min = -5.0,appe_max=5.0;
        REAL8 alphappe_min = -1000.0,alphappe_max=1000.0;
        REAL8 bppe_min = -5.0,bppe_max=5.0;
        REAL8 betappe_min = -1000.0,betappe_max=1000.0;
        char aPPEparam[64]="";
        char alphaPPEparam[64]="";
        /* phase parameters */
        char bPPEparam[64]="";
        char betaPPEparam[64]="";
        int counters[4]={0};
        do
        {
            sprintf(aPPEparam, "%s%d","aPPE",++counters[0]);
            if (checkParamInList(ppt->value,aPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, aPPEparam, 0.0, appe_min, appe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(alphaPPEparam, "%s%d","alphaPPE",++counters[1]);
            if (checkParamInList(ppt->value,alphaPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, alphaPPEparam, 0.0, alphappe_min, alphappe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(bPPEparam, "%s%d","bPPE",++counters[2]);
            if (checkParamInList(ppt->value,bPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, bPPEparam, 0.0, bppe_min, bppe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(betaPPEparam, "%s%d","betaPPE",++counters[3]);
            if (checkParamInList(ppt->value,betaPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, betaPPEparam, 0.0, betappe_min, betappe_max, LALINFERENCE_PARAM_LINEAR);

        } while((checkParamInList(ppt->value,aPPEparam))||(checkParamInList(ppt->value,alphaPPEparam))||(checkParamInList(ppt->value,bPPEparam))||(checkParamInList(ppt->value,betaPPEparam)));
        if ((counters[0]!=counters[1])||(counters[2]!=counters[3])) {fprintf(stderr,"Unequal number of PPE parameters detected! Check your command line!\n"); exit(-1);}
    }

}
