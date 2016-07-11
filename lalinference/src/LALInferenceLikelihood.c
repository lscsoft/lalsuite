/*
 *  LALInferenceLikelihood.c:  Bayesian Followup likelihood functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys and John Veitch, Will M. Farr
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

#include <complex.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInference.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_complex_math.h>
#include <lal/LALInferenceTemplate.h>

#include "logaddexp.h"

typedef enum
{
  GAUSSIAN,
  STUDENTT,
  MARGPHI,
  MARGTIME,
  MARGTIMEPHI
} LALInferenceLikelihoodFlags;


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static REAL8 LALInferenceFusedFreqDomainLogLikelihood(LALInferenceVariables *currentParams,
                                               LALInferenceIFOData *data,
                                               LALInferenceModel *model,
                                               LALInferenceLikelihoodFlags marginalisationflags);

static double integrate_interpolated_log(double h, REAL8 *log_ys, size_t n, double *imean, size_t *imax);


void LALInferenceInitLikelihood(LALInferenceRunState *runState)
{
    char help[]="\
    ----------------------------------------------\n\
    --- Likelihood Arguments ---------------------\n\
    ----------------------------------------------\n\
    (--zeroLogLike)                  Use flat, null likelihood\n\
    (--studentTLikelihood)           Use the Student-T Likelihood that marginalizes over noise\n\
    (--correlatedGaussianLikelihood) Use analytic, correlated Gaussian for Likelihood Z=-21.3\n\
    (--bimodalGaussianLikelihood)    Use analytic, bimodal correlated Gaussian for Likelihood Z=-25.9\n\
    (--rosenbrockLikelihood)         Use analytic, Rosenbrock banana for Likelihood\n\
    (--noiseonly)                    Using noise-only likelihood\n\
    (--margphi)                      Using marginalised phase likelihood\n\
    (--margtime)                     Using marginalised time likelihood\n\
    (--margtimephi)                  Using marginalised in time and phase likelihood\n\
    \n";

    /* Print command line arguments if help requested */
    LALInferenceIFOData *ifo=NULL;
    if(runState == NULL || LALInferenceGetProcParamVal(runState->commandLine,"--help"))
    {
        fprintf(stdout,"%s",help);
        if (runState){
            ifo=runState->data;
            while(ifo) {
                fprintf(stdout,"(--dof-%s DoF)\tDegrees of freedom for %s\n",ifo->name,ifo->name);
                ifo=ifo->next;
            }
        }
        return;
    }

    ProcessParamsTable *commandLine=runState->commandLine;
    ifo=runState->data;

    LALInferenceThreadState *thread = runState->threads[0];

    REAL8 nullLikelihood = 0.0; // Populated if such a thing exists

   if (LALInferenceGetProcParamVal(commandLine, "--zeroLogLike")) {
    /* Use zero log(L) */
    runState->likelihood=&LALInferenceZeroLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood")) {
    runState->likelihood=&LALInferenceCorrelatedAnalyticLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood")) {
    runState->likelihood=&LALInferenceBimodalCorrelatedAnalyticLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood")) {
    runState->likelihood=&LALInferenceRosenbrockLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--studentTLikelihood")) {
    fprintf(stderr, "Using Student's T Likelihood.\n");
    runState->likelihood=&LALInferenceFreqDomainStudentTLogLikelihood;

    /* Set the noise model evidence to the student t model value */
    LALInferenceTemplateNullFreqdomain(thread->model);
    LALInferenceTemplateFunction temp = thread->model->templt;
    thread->model->templt = &LALInferenceTemplateNullFreqdomain;
    REAL8 noiseZ = LALInferenceFreqDomainStudentTLogLikelihood(thread->currentParams, runState->data, thread->model);
    thread->model->templt = temp;
    LALInferenceAddVariable(runState->algorithmParams, "logZnoise", &noiseZ, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    fprintf(stdout,"Student-t Noise evidence %lf\n", noiseZ);

   } else if (LALInferenceGetProcParamVal(commandLine, "--margphi")) {
    fprintf(stderr, "Using marginalised phase likelihood.\n");
    runState->likelihood=&LALInferenceMarginalisedPhaseLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--margtime")) {
    fprintf(stderr, "Using marginalised time likelihood.\n");
    runState->likelihood=&LALInferenceMarginalisedTimeLogLikelihood;
   } else if (LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
     fprintf(stderr, "Using marginalised in time and phase likelihood.\n");
     runState->likelihood=&LALInferenceMarginalisedTimePhaseLogLikelihood;
     //LALInferenceAddVariable(runState->currentParams, "margtimephi", &margphi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
   }
     else if (LALInferenceGetProcParamVal(commandLine, "--roqtime_steps")) {
     fprintf(stderr, "Using ROQ in likelihood.\n");
     runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
    }
     else if (LALInferenceGetProcParamVal(commandLine, "--fastSineGaussianLikelihood")){
      fprintf(stderr, "WARNING: Using Fast SineGaussian likelihood and WF for LIB.\n");
      runState->likelihood=&LALInferenceFastSineGaussianLogLikelihood;
   } else {
      runState->likelihood=&LALInferenceUndecomposedFreqDomainLogLikelihood;
   }

   /* Try to determine a model-less likelihood, if such a thing makes sense */
   if (runState->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood){

                nullLikelihood = LALInferenceNullLogLikelihood(runState->data);

    }
    else if (runState->likelihood==&LALInferenceFreqDomainStudentTLogLikelihood ||
       runState->likelihood==&LALInferenceMarginalisedTimeLogLikelihood ||
       runState->likelihood==&LALInferenceMarginalisedTimePhaseLogLikelihood ||
       runState->likelihood==&LALInferenceMarginalisedPhaseLogLikelihood) {

       void *oldtemplate = runState->threads[0]->model->templt;
       runState->threads[0]->model->templt = &LALInferenceTemplateNullFreqdomain;
       nullLikelihood = runState->likelihood(thread->currentParams, runState->data, thread->model);
       runState->threads[0]->model->templt = oldtemplate;
    
   }

    //null log likelihood logic doesn't work with noise parameters
    if (LALInferenceGetProcParamVal(runState->commandLine,"--psdFit") ||
       LALInferenceGetProcParamVal(runState->commandLine,"--psd-fit") ||
       LALInferenceGetProcParamVal(runState->commandLine,"--glitchFit") ||
       LALInferenceGetProcParamVal(runState->commandLine,"--glitch-fit")) {
           nullLikelihood = 0.0;
           ifo = runState->data;
           while (ifo != NULL) {
               ifo->nullloglikelihood = 0.0;
               ifo = ifo->next;
           }
    }

   INT4 t;
   for(t=0; t < runState->nthreads; t++)
       runState->threads[t]->nullLikelihood = nullLikelihood;

   LALInferenceAddVariable(runState->proposalArgs, "nullLikelihood", &nullLikelihood,
                           LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    return;
}


const char *non_intrinsic_params[] = {"rightascension", "declination", "polarisation", "time",
                                "deltaLogL", "logL", "deltaloglH1", "deltaloglL1", "deltaloglV1",
                                "logw", "logPrior","hrss","loghrss", NULL};

LALInferenceVariables LALInferenceGetInstrinsicParams(LALInferenceVariables *currentParams)
/***************************************************************/
/* Return a variables structure containing only intrinsic      */
/* parameters.                                                 */
/***************************************************************/
{
    // TODO: add pointer to template function here.
    // (otherwise same parameters but different template will lead to no re-computation!!)
    LALInferenceVariables intrinsicParams;
    const char **non_intrinsic_param = non_intrinsic_params;

    intrinsicParams.head      = NULL;
    intrinsicParams.dimension = 0;
    LALInferenceCopyVariables(currentParams, &intrinsicParams);

    while (*non_intrinsic_param) {
        if (LALInferenceCheckVariable(&intrinsicParams, *non_intrinsic_param))
            LALInferenceRemoveVariable(&intrinsicParams, *non_intrinsic_param);
        non_intrinsic_param++;
    }

    return intrinsicParams;
}

/* Check to see if item is in the NULL-terminated array.
 If so, return 1. Otherwise, add it to the array and return 0
 */
static int checkItemAndAdd(void *item, void **array);
static int checkItemAndAdd(void *item, void **array)
{
  UINT4 i=0;
  if(!array || !item) return 0;
  while(array[i])
  {
    if(array[i++]==item) return 1;
  }
  array[i]=item;
  return 0;
}

/* ============ Likelihood computations: ========== */

/**
 * For testing purposes (for instance sampling the prior), likelihood that returns 0.0 = log(1) every
 * time.  Activated with the --zeroLogLike command flag.
 */
REAL8 LALInferenceZeroLogLikelihood(LALInferenceVariables *currentParams,
                                    LALInferenceIFOData UNUSED *data,
                                    LALInferenceModel UNUSED *model) {

    INT4 SKY_FRAME=0;
    REAL8 ra,dec,GPSdouble;

    if(LALInferenceCheckVariable(currentParams,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(currentParams,"SKY_FRAME");

    if(SKY_FRAME==1)
    {
      REAL8 t0=LALInferenceGetREAL8Variable(currentParams,"t0");
      REAL8 alph=acos(LALInferenceGetREAL8Variable(currentParams,"cosalpha"));
      REAL8 theta=LALInferenceGetREAL8Variable(currentParams,"azimuth");
      LALInferenceDetFrameToEquatorial(data->detector,data->next->detector,
                                       t0,alph,theta,&GPSdouble,&ra,&dec);
      LALInferenceAddVariable(currentParams,"rightascension",&ra,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      LALInferenceAddVariable(currentParams,"declination",&dec,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      LALInferenceAddVariable(currentParams,"time",&GPSdouble,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    }

  REAL8 zero = 0.0;
  LALInferenceAddVariable(currentParams,"optimal_snr",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  LALInferenceAddVariable(currentParams,"matched_filter_snr",&zero,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  return 0.0;
}




REAL8 LALInferenceUndecomposedFreqDomainLogLikelihood(LALInferenceVariables *currentParams,
                                                      LALInferenceIFOData *data,
                                                      LALInferenceModel *model)
{
  return LALInferenceFusedFreqDomainLogLikelihood(currentParams,
                                                 data,
                                                 model,
                                                  GAUSSIAN);
}


static REAL8 LALInferenceFusedFreqDomainLogLikelihood(LALInferenceVariables *currentParams,
                                                        LALInferenceIFOData *data,
                                                        LALInferenceModel *model,
                                                        LALInferenceLikelihoodFlags marginalisationflags)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  double Fplus, Fcross;
  //double diffRe, diffIm;
  //double dataReal, dataImag;
  double glitchReal=0.0, glitchImag=0.0;
  //REAL8 plainTemplateReal, plainTemplateImag;
  //REAL8 templateReal=0.0, templateImag=0.0;
  int i, j, lower, upper, ifo;
  LALInferenceIFOData *dataPtr;
  double ra=0.0, dec=0.0, psi=0.0, gmst=0.0;
  double GPSdouble=0.0, t0=0.0;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift=0;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit=0.0, re, im, dre, dim, newRe, newIm;
  double timeTmp;
  double mc;
  /* Burst templates are generated at hrss=1, thus need to rescale amplitude */
  double amp_prefactor=1.0;

  COMPLEX16FrequencySeries *calFactor = NULL;
  COMPLEX16 calF = 0.0;

  char freqVarName[VARNAME_MAX];
  char ampVarName[VARNAME_MAX];
  char phaseVarName[VARNAME_MAX];

  REAL8Vector *logfreqs = NULL;
  REAL8Vector *amps = NULL;
  REAL8Vector *phases = NULL;

  UINT4 spcal_active = 0;
  REAL8 calamp=0.0;
  REAL8 calpha=0.0;
  REAL8 cos_calpha=cos(calpha);
  REAL8 sin_calpha=sin(calpha);
  UINT4 constantcal_active=0;

  /* ROQ likelihood stuff */
  REAL8 d_inner_h=0.0;


  if (LALInferenceCheckVariable(currentParams, "spcal_active") && (*(UINT4 *)LALInferenceGetVariable(currentParams, "spcal_active"))) {
    spcal_active = 1;
  }
  if (LALInferenceCheckVariable(currentParams, "constantcal_active") && (*(UINT4 *)LALInferenceGetVariable(currentParams, "constantcal_active"))) {
   constantcal_active = 1;
  }
  if (spcal_active && constantcal_active){
    fprintf(stderr,"ERROR: cannot use spline and constant calibration error marginalization together. Exiting...\n");
    exit(1);
  }
  if (model->roq_flag && constantcal_active){
    fprintf(stderr,"ERROR: cannot use ROQ likelihood and constant calibration error marginalization together. Exiting...\n");
    exit(1);
  }

  REAL8 degreesOfFreedom=2.0;
  REAL8 chisq=0.0;
  /* margphi params */
  //REAL8 Rre=0.0,Rim=0.0;
  REAL8 D=0.0,S=0.0;
  COMPLEX16 Rcplx=0.0;
  int margphi=0;
  int margtime=0;
  REAL8 desired_tc=0.0;
  if (marginalisationflags==MARGPHI || marginalisationflags==MARGTIMEPHI)
    margphi=1;
  if (marginalisationflags==MARGTIME || marginalisationflags==MARGTIMEPHI)
    margtime=1;

  LALStatus status;
  memset(&status,0,sizeof(status));

  if(data==NULL) {XLAL_ERROR_REAL8(XLAL_EINVAL,"ERROR: Encountered NULL data pointer in likelihood\n");}

  int Nifos=0;
  for(dataPtr=data;dataPtr;dataPtr=dataPtr->next) Nifos++;
  void *generatedFreqModels[1+Nifos];
  for(i=0;i<=Nifos;i++) generatedFreqModels[i]=NULL;

  //noise model meta parameters
  gsl_matrix *nparams = NULL;//pointer to matrix holding noise parameters

  gsl_matrix *psdBandsMin  = NULL;//pointer to matrix holding min frequencies for psd model
  gsl_matrix *psdBandsMax = NULL;//pointer to matrix holding max frequencies for psd model

  //different formats for storing glitch model for DWT, FFT, and integration
  gsl_matrix *glitchFD=NULL;

  int Nblock = 1;            //number of frequency blocks per IFO
  int psdFlag = 0;           //flag for including psd fitting
  int glitchFlag = 0;   //flag for including glitch model
  int signalFlag = 1;   //flag for including signal model

  //check if psd parameters are included in the model
  psdFlag = 0;
  if(LALInferenceCheckVariable(currentParams, "psdScaleFlag"))
    psdFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "psdScaleFlag"));
  if(psdFlag)
  {
    //if so, store current noise parameters in easily accessible matrix
    nparams = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "psdscale"));
    Nblock = (int)nparams->size2;

    psdBandsMin = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "psdBandsMin"));
    psdBandsMax = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "psdBandsMax"));

  }
  double alpha[Nblock];
  double lnalpha[Nblock];

  double psdBandsMin_array[Nblock];
  double psdBandsMax_array[Nblock];

  //check if glitch model is being used
  glitchFlag = 0;
  if(LALInferenceCheckVariable(currentParams,"glitchFitFlag"))
    glitchFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "glitchFitFlag"));
  if(glitchFlag)
    glitchFD = *((gsl_matrix **)LALInferenceGetVariable(currentParams, "morlet_FD"));

  //check if signal model is being used
  signalFlag=1;
  if(LALInferenceCheckVariable(currentParams, "signalModelFlag"))
    signalFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "signalModelFlag"));

  int freq_length=0,time_length=0;
  COMPLEX16Vector * dh_S_tilde=NULL;
  COMPLEX16Vector * dh_S_phase_tilde = NULL;
  REAL8Vector *dh_S=NULL;
  REAL8Vector *dh_S_phase = NULL;
  /* Setup times to integrate over */
  freq_length = data->freqData->data->length;
  time_length = 2*(freq_length-1);

  /* Desired tc == 2 seconds before buffer end.  Only used during
     margtime{phi} to try to place the waveform in a reasonable
     place before time-shifting */
  deltaT = data->timeData->deltaT;
  REAL8 epoch = XLALGPSGetREAL8(&(data->freqData->epoch));
  desired_tc = epoch + (time_length-1)*deltaT - 2.0;

  if(signalFlag)
  {
    if(LALInferenceCheckVariable(currentParams,"logmc")){
      mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
      LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    }
    if(LALInferenceCheckVariable(currentParams, "loghrss")){
      amp_prefactor = exp(*(REAL8*)LALInferenceGetVariable(currentParams,"loghrss"));
    }
    else if (LALInferenceCheckVariable(currentParams, "hrss")){
      amp_prefactor = (*(REAL8*)LALInferenceGetVariable(currentParams,"hrss"));
    }

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(currentParams,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(currentParams,"SKY_FRAME");

    if(SKY_FRAME==0){
      /* determine source's sky location & orientation parameters: */
      ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
      dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
    }
    else
    {
	    if(Nifos<2){
		    fprintf(stderr,"ERROR: Cannot use --detector-frame with less than 2 detectors!\n");
		    exit(1);
	    }
      if(!margtime)
      {
         t0=LALInferenceGetREAL8Variable(currentParams,"t0");
      }
      else /* Use the desired end time to compute the mapping to ra,dec */
      {
              t0=desired_tc;
      }
      REAL8 alph=acos(LALInferenceGetREAL8Variable(currentParams,"cosalpha"));
      REAL8 theta=LALInferenceGetREAL8Variable(currentParams,"azimuth");
      LALInferenceDetFrameToEquatorial(data->detector,data->next->detector,
                                       t0,alph,theta,&GPSdouble,&ra,&dec);
      LALInferenceAddVariable(currentParams,"rightascension",&ra,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      LALInferenceAddVariable(currentParams,"declination",&dec,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      if(!margtime) LALInferenceAddVariable(currentParams,"time",&GPSdouble,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    }
    psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
    if(!margtime)
	      GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
    else
	      GPSdouble = XLALGPSGetREAL8(&(data->freqData->epoch));


    // Add phase parameter set to 0 for calculation
    if(margphi ){
      REAL8 phi0=0.0;
      if(LALInferenceCheckVariable(currentParams,"phase")) LALInferenceRemoveVariable(currentParams,"phase");
      LALInferenceAddVariable(currentParams, "phase",&phi0,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    }
  }


  if(margtime)
  {
    GPSdouble = desired_tc;
    dh_S_tilde = XLALCreateCOMPLEX16Vector(freq_length);
    dh_S = XLALCreateREAL8Vector(time_length);

    if (dh_S_tilde ==NULL || dh_S == NULL)
      XLAL_ERROR_REAL8(XLAL_ENOMEM, "Out of memory in LALInferenceMarginalisedTimeLogLikelihood.");

    for (i = 0; i < freq_length; i++) {
      dh_S_tilde->data[i] = 0.0;
    }

    if (margphi) {
      dh_S_phase_tilde = XLALCreateCOMPLEX16Vector(freq_length);
      dh_S_phase = XLALCreateREAL8Vector(time_length);

      if (dh_S_phase_tilde == NULL || dh_S_phase == NULL) {
	XLAL_ERROR_REAL8(XLAL_ENOMEM, "Out of memory in time-phase marginalised likelihood.");
      }

      for (i = 0; i < freq_length; i++) {
	dh_S_phase_tilde->data[i] = 0.0;
      }
    }
  }

  /* figure out GMST: */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  chisquared = 0.0;
  REAL8 loglikelihood = 0.0;

  /* Reset SNR */
  model->SNR = 0.0;

  /* loop over data (different interferometers): */
  for(dataPtr=data,ifo=0; dataPtr; dataPtr=dataPtr->next,ifo++) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* model->params (set, e.g., from the trigger value).     */

    /* Reset log-likelihood */
    model->ifo_loglikelihoods[ifo] = 0.0;
    model->ifo_SNRs[ifo] = 0.0;
    COMPLEX16 this_ifo_d_inner_h = 0.0;
    REAL8 this_ifo_s = 0.0;
    // Check if student-t likelihood is being used
    if(marginalisationflags==STUDENTT)
    {
      /* extract the element from the "df" vector that carries the current Ifo's name: */
      CHAR df_variable_name[64];
      snprintf(df_variable_name,sizeof(df_variable_name),"df_%s",dataPtr->name);
      if(LALInferenceCheckVariable(currentParams,df_variable_name)){
        printf("Found variable %s\n",df_variable_name);
        degreesOfFreedom = *(REAL8*) LALInferenceGetVariable(currentParams,df_variable_name);
      }
      else {
        degreesOfFreedom = dataPtr->STDOF;
      }
      if (!(degreesOfFreedom>0)) {
        XLALPrintError(" ERROR in StudentTLogLikelihood(): degrees-of-freedom parameter must be positive.\n");
        XLAL_ERROR_REAL8(XLAL_EDOM);
      }
    }

    if(signalFlag){

      /* Check to see if this buffer has already been filled with the signal.
       Different dataPtrs can share the same signal buffer to avoid repeated
       calls to template */
      if(!checkItemAndAdd((void *)(model->freqhPlus), generatedFreqModels))
      {
        /* Compare parameter values with parameter values corresponding  */
        /* to currently stored template; ignore "time" variable:         */
        if (LALInferenceCheckVariable(model->params, "time")) {
          timeTmp = *(REAL8 *) LALInferenceGetVariable(model->params, "time");
          LALInferenceRemoveVariable(model->params, "time");
        }
        else timeTmp = GPSdouble;

        LALInferenceCopyVariables(currentParams, model->params);
        // Remove time variable so it can be over-written (if it was pinned)
        if(LALInferenceCheckVariable(model->params,"time")) LALInferenceRemoveVariable(model->params,"time");
        LALInferenceAddVariable(model->params, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);

        INT4 errnum=0;
        XLAL_TRY(model->templt(model),errnum);
        errnum&=~XLAL_EFUNC;
        if(errnum!=XLAL_SUCCESS)
        {
          switch(errnum)
          {
            case XLAL_EUSR0: /* Template generation failed in a known way, set -Inf likelihood */
              return (-DBL_MAX);
              break;
            default: /* Panic! */
              fprintf(stderr,"Unhandled error in template generation - exiting!\n");
              fprintf(stderr,"XLALError: %d, %s\n",errnum,XLALErrorString(errnum));
              exit(1);
              break;
          }

        }

        if (model->domain == LAL_SIM_DOMAIN_TIME) {
          /* TD --> FD. */
          LALInferenceExecuteFT(model);
        }
      }

        /* Template is now in model->timeFreqhPlus and hCross */

        /* Calibration stuff if necessary */
        /*spline*/
        if (spcal_active) {
          snprintf(freqVarName, VARNAME_MAX, "%s_spcal_logfreq", dataPtr->name);
          snprintf(ampVarName, VARNAME_MAX, "%s_spcal_amp", dataPtr->name);
          snprintf(phaseVarName, VARNAME_MAX, "%s_spcal_phase", dataPtr->name);

          logfreqs = *(REAL8Vector **)LALInferenceGetVariable(currentParams, freqVarName);
          amps = *(REAL8Vector **)LALInferenceGetVariable(currentParams, ampVarName);
          phases = *(REAL8Vector **)LALInferenceGetVariable(currentParams, phaseVarName);

	if (model->roq_flag) {

             LALInferenceSplineCalibrationFactorROQ(logfreqs, amps, phases,
						model->roq->frequencyNodesLinear,
						&(model->roq->calFactorLinear),
						model->roq->frequencyNodesQuadratic,
						&(model->roq->calFactorQuadratic));
          }

	else{
          if (calFactor == NULL) {
            calFactor = XLALCreateCOMPLEX16FrequencySeries("calibration factors",
                       &(dataPtr->freqData->epoch),
                       0, dataPtr->freqData->deltaF,
                       &lalDimensionlessUnit,
                       dataPtr->freqData->data->length);
          }
          LALInferenceSplineCalibrationFactor(logfreqs, amps, phases, calFactor);
	}

        }
        /*constant*/
        if (constantcal_active){
          char CA_A[10]="";
          sprintf(CA_A,"%s_%s","calamp",dataPtr->name);
          if (LALInferenceCheckVariable(currentParams, CA_A))
            calamp=(*(REAL8*) LALInferenceGetVariable(currentParams, CA_A));
          else
            calamp=0.0;
          char CP_A[10]="";
          sprintf(CP_A,"%s_%s","calpha",dataPtr->name);
          if (LALInferenceCheckVariable(currentParams, CP_A))
            calpha=(*(REAL8*) LALInferenceGetVariable(currentParams, CP_A));
          else
            calpha=0.0;
          cos_calpha=cos(calpha);
          sin_calpha=-sin(calpha);
        }
        /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
        XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, ra, dec, psi, gmst);

        /* signal arrival time (relative to geocenter); */
        timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location, ra, dec, &GPSlal);
        /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
        /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
        if (margtime)
          /* If we are marginalising over time, we want the
	      freq-domain signal to have tC = epoch, so we shift it
	      from the model's "time" parameter to epoch */
          timeshift =  (epoch - (*(REAL8 *) LALInferenceGetVariable(model->params, "time"))) + timedelay;
        else
          timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(model->params, "time"))) + timedelay;
        twopit    = LAL_TWOPI * timeshift;

        /* For burst, add the right hrss in the amplitude. */
        Fplus*=amp_prefactor;
        Fcross*=amp_prefactor;

        dataPtr->fPlus = Fplus;
        dataPtr->fCross = Fcross;
        dataPtr->timeshift = timeshift;
    }//end signalFlag condition

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);

    /* Employ a trick here for avoiding cos(...) and sin(...) in time
       shifting.  We need to multiply each template frequency bin by
       exp(-J*twopit*deltaF*i) = exp(-J*twopit*deltaF*(i-1)) +
       exp(-J*twopit*deltaF*(i-1))*(exp(-J*twopit*deltaF) - 1) .  This
       recurrance relation has the advantage that the error growth is
       O(sqrt(N)) for N repetitions. */

    /* See, for example,

       Press, Teukolsky, Vetteling & Flannery, 2007.  Numerical
       Recipes, Third Edition, Chapter 5.4.

       Singleton, 1967. On computing the fast Fourier
       transform. Comm. ACM, vol. 10, 647–654. */

    /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
    dim = -sin(twopit*deltaF);
    dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);

    //Set up noise PSD meta parameters
    for(i=0; i<Nblock; i++)
    {
      if(psdFlag)
      {
        alpha[i]   = gsl_matrix_get(nparams,ifo,i);
        lnalpha[i] = log(alpha[i]);

        psdBandsMin_array[i] = gsl_matrix_get(psdBandsMin,ifo,i);
        psdBandsMax_array[i] = gsl_matrix_get(psdBandsMax,ifo,i);
      }
      else
      {
        alpha[i]=1.0;
        lnalpha[i]=0.0;
      }
    }

    if (model->roq_flag) {

	double complex weight_iii;

	if (spcal_active){

	    for(unsigned int iii=0; iii < model->roq->frequencyNodesLinear->length; iii++){

			complex double template_EI = model->roq->calFactorLinear->data[iii] * (dataPtr->fPlus*model->roq->hptildeLinear->data->data[iii] + dataPtr->fCross*model->roq->hctildeLinear->data->data[iii] );

			weight_iii = gsl_spline_eval (dataPtr->roq->weights_linear[iii].spline_real_weight_linear, timeshift, dataPtr->roq->weights_linear[iii].acc_real_weight_linear) + I*gsl_spline_eval (dataPtr->roq->weights_linear[iii].spline_imag_weight_linear, timeshift, dataPtr->roq->weights_linear[iii].acc_imag_weight_linear);

			this_ifo_d_inner_h += ( weight_iii * ( conj( template_EI ) ) );
		}

		for(unsigned int jjj=0; jjj < model->roq->frequencyNodesQuadratic->length; jjj++){

			this_ifo_s += dataPtr->roq->weightsQuadratic[jjj] * creal( conj( model->roq->calFactorQuadratic->data[jjj] * (model->roq->hptildeQuadratic->data->data[jjj]*dataPtr->fPlus + model->roq->hctildeQuadratic->data->data[jjj]*dataPtr->fCross) ) * ( model->roq->calFactorQuadratic->data[jjj] * (model->roq->hptildeQuadratic->data->data[jjj]*dataPtr->fPlus + model->roq->hctildeQuadratic->data->data[jjj]*dataPtr->fCross) ) );
		}
	}

	else{

		for(unsigned int iii=0; iii < model->roq->frequencyNodesLinear->length; iii++){

			complex double template_EI = dataPtr->fPlus*model->roq->hptildeLinear->data->data[iii] + dataPtr->fCross*model->roq->hctildeLinear->data->data[iii];

			weight_iii = gsl_spline_eval (dataPtr->roq->weights_linear[iii].spline_real_weight_linear, timeshift, dataPtr->roq->weights_linear[iii].acc_real_weight_linear) + I*gsl_spline_eval (dataPtr->roq->weights_linear[iii].spline_imag_weight_linear, timeshift, dataPtr->roq->weights_linear[iii].acc_imag_weight_linear);

			this_ifo_d_inner_h += weight_iii*conj(template_EI) ;

		}

		for(unsigned int jjj=0; jjj < model->roq->frequencyNodesQuadratic->length; jjj++){
			complex double template_EI = model->roq->hptildeQuadratic->data->data[jjj]*Fplus + model->roq->hctildeQuadratic->data->data[jjj]*Fcross;

			this_ifo_s += dataPtr->roq->weightsQuadratic[jjj] * creal( conj(template_EI) * (template_EI) );
					}
	}

	d_inner_h += creal(this_ifo_d_inner_h);
	S += this_ifo_s;
	model->ifo_loglikelihoods[ifo] = creal(this_ifo_d_inner_h) - (0.5*this_ifo_s) + dataPtr->nullloglikelihood;

	loglikelihood += model->ifo_loglikelihoods[ifo];

	char varname[VARNAME_MAX];
    	sprintf(varname,"%s_optimal_snr",dataPtr->name);
    	REAL8 this_ifo_snr = sqrt(this_ifo_s);
	model->ifo_SNRs[ifo] = this_ifo_snr;
    	LALInferenceAddREAL8Variable(currentParams,varname,this_ifo_snr,LALINFERENCE_PARAM_OUTPUT);

	sprintf(varname,"%s_cplx_snr_amp",dataPtr->name);
	REAL8 cplx_snr_amp = cabs(this_ifo_d_inner_h)/this_ifo_snr;
    	LALInferenceAddREAL8Variable(currentParams,varname,cplx_snr_amp,LALINFERENCE_PARAM_OUTPUT);

    	sprintf(varname,"%s_cplx_snr_arg",dataPtr->name);
    	REAL8 cplx_snr_phase = carg(this_ifo_d_inner_h);
    	LALInferenceAddREAL8Variable(currentParams,varname,cplx_snr_phase,LALINFERENCE_PARAM_OUTPUT);


    }

    else{

    	REAL8 *psd=&(dataPtr->oneSidedNoisePowerSpectrum->data->data[lower]);
    COMPLEX16 *dtilde=&(dataPtr->freqData->data->data[lower]);
    COMPLEX16 *hptilde=&(model->freqhPlus->data->data[lower]);
    COMPLEX16 *hctilde=&(model->freqhCross->data->data[lower]);
    COMPLEX16 diff=0.0;
    COMPLEX16 template=0.0;
    REAL8 templatesq=0.0;
    REAL8 this_ifo_S=0.0;
    COMPLEX16 this_ifo_Rcplx=0.0;

    for (i=lower,chisq=0.0,re = cos(twopit*deltaF*i),im = -sin(twopit*deltaF*i);
         i<=upper;
         i++, psd++, hptilde++, hctilde++, dtilde++,
         newRe = re + re*dre - im*dim,
         newIm = im + re*dim + im*dre,
         re = newRe, im = newIm)
    {

      COMPLEX16 d=*dtilde;
      /* Normalise PSD to our funny standard (see twoDeltaTOverN
	 below). */
      REAL8 sigmasq=(*psd)*deltaT*deltaT;

      if (constantcal_active) {
        REAL8 dre_tmp= creal(d)*cos_calpha - cimag(d)*sin_calpha;
        REAL8 dim_tmp = creal(d)*sin_calpha + cimag(d)*cos_calpha;
        dre_tmp/=(1.0+calamp);
        dim_tmp/=(1.0+calamp);

        d=crect(dre_tmp,dim_tmp);
        sigmasq/=((1.0+calamp)*(1.0+calamp));
      }

      REAL8 singleFreqBinTerm;


      /* Add noise PSD parameters to the model */
      if(psdFlag)
      {
        for(j=0; j<Nblock; j++)
        {
          if (i >= psdBandsMin_array[j] && i <= psdBandsMax_array[j])
          {
            sigmasq  *= alpha[j];
            loglikelihood -= lnalpha[j];
          }
        }
      }

      //subtract GW model from residual
      diff = d;

      if(signalFlag){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      COMPLEX16 plainTemplate = Fplus*(*hptilde)+Fcross*(*hctilde);

      /* Do time shifting */
      template = plainTemplate * (re + I*im);

      if (spcal_active) {
          calF = calFactor->data->data[i];
          template = template*calF;
      }

      diff -= template;

      }//end signal subtraction

      //subtract glitch model from residual
      if(glitchFlag)
      {
        /* fourier amplitudes of glitches */
        glitchReal = gsl_matrix_get(glitchFD,ifo,2*i);
        glitchImag = gsl_matrix_get(glitchFD,ifo,2*i+1);
        COMPLEX16 glitch = glitchReal + I*glitchImag;
        diff -=glitch*deltaT;

      }//end glitch subtraction

      templatesq=creal(template)*creal(template) + cimag(template)*cimag(template);
      REAL8 datasq = creal(d)*creal(d)+cimag(d)*cimag(d);
      D+=TwoDeltaToverN*datasq/sigmasq;
      this_ifo_S+=TwoDeltaToverN*templatesq/sigmasq;
      COMPLEX16 dhstar = TwoDeltaToverN*d*conj(template)/sigmasq;
      this_ifo_Rcplx+=dhstar;
      Rcplx+=dhstar;

      switch(marginalisationflags)
      {
        case GAUSSIAN:
        {
          REAL8 diffsq = creal(diff)*creal(diff)+cimag(diff)*cimag(diff);
          chisq = TwoDeltaToverN*diffsq/sigmasq;
          singleFreqBinTerm = chisq;
          chisquared  += singleFreqBinTerm;
          model->ifo_loglikelihoods[ifo] -= singleFreqBinTerm;
          break;
        }
        case STUDENTT:
        {
          REAL8 diffsq = creal(diff)*creal(diff)+cimag(diff)*cimag(diff);
          chisq = TwoDeltaToverN*diffsq/sigmasq;
          singleFreqBinTerm = ((degreesOfFreedom+2.0)/2.0) * log(1.0 + chisq/degreesOfFreedom) ;
          chisquared  += singleFreqBinTerm;
          model->ifo_loglikelihoods[ifo] -= singleFreqBinTerm;
          break;
        }
        case MARGTIME:
        case MARGTIMEPHI:
        {
          loglikelihood+=-TwoDeltaToverN*(templatesq+datasq)/sigmasq;

          /* Note: No Factor of 2 here, since we are using the 2-sided
	     COMPLEX16FFT.  Also, we use d*conj(h) because we are
	     using a complex->real *inverse* FFT to compute the
	     time-series of likelihoods. */
          dh_S_tilde->data[i] += TwoDeltaToverN * d * conj(template) / sigmasq;

          if (margphi) {
            /* This is the other phase quadrature */
            dh_S_phase_tilde->data[i] += TwoDeltaToverN * d * conj(I*template) / sigmasq;
          }

          break;
        }
        case MARGPHI:
        {
          break;
        }
        default:
          break;
      }



    } /* End loop over freq bins */
    switch(marginalisationflags)
    {
    case GAUSSIAN:
    case STUDENTT:
      loglikelihood += model->ifo_loglikelihoods[ifo];
      break;
    case MARGTIME:
    case MARGPHI:
    case MARGTIMEPHI:
      /* These are non-separable likelihoods, so single IFO log(L)
	 doesn't make sense. */
      model->ifo_loglikelihoods[ifo] = 0.0;
      break;
    default:
      break;
    }
    S+=this_ifo_S;
    char varname[VARNAME_MAX];
    sprintf(varname,"%s_optimal_snr",dataPtr->name);
    LALInferenceAddREAL8Variable(currentParams,varname,sqrt(2.0*this_ifo_S),LALINFERENCE_PARAM_OUTPUT);

    sprintf(varname,"%s_cplx_snr_amp",dataPtr->name);
    REAL8 cplx_snr_amp=0.0;
    REAL8 cplx_snr_phase=carg(this_ifo_Rcplx);
    if(this_ifo_S > 0) cplx_snr_amp=2.0*cabs(this_ifo_Rcplx)/sqrt(2.0*this_ifo_S);

    LALInferenceAddREAL8Variable(currentParams,varname,cplx_snr_amp,LALINFERENCE_PARAM_OUTPUT);

    sprintf(varname,"%s_cplx_snr_arg",dataPtr->name);
    LALInferenceAddREAL8Variable(currentParams,varname,cplx_snr_phase,LALINFERENCE_PARAM_OUTPUT);

   /* Clean up calibration if necessary */
    if (!(calFactor == NULL)) {
      XLALDestroyCOMPLEX16FrequencySeries(calFactor);
      calFactor = NULL;
    }
  } /* end loop over detectors */

  }
  if (model->roq_flag){



	REAL8 OptimalSNR=sqrt(S);
        REAL8 MatchedFilterSNR = d_inner_h/OptimalSNR;
        /* fprintf(stderr, "%f\n", d_inner_h - 0.5*OptimalSNR); */
        LALInferenceAddVariable(currentParams,"optimal_snr",&OptimalSNR,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
        LALInferenceAddVariable(currentParams,"matched_filter_snr",&MatchedFilterSNR,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

	model->SNR = OptimalSNR;

	if ( model->roq->hptildeLinear ) XLALDestroyCOMPLEX16FrequencySeries(model->roq->hptildeLinear);
  	if ( model->roq->hctildeLinear ) XLALDestroyCOMPLEX16FrequencySeries(model->roq->hctildeLinear);
  	if ( model->roq->hptildeQuadratic ) XLALDestroyCOMPLEX16FrequencySeries(model->roq->hptildeQuadratic);
  	if ( model->roq->hctildeQuadratic ) XLALDestroyCOMPLEX16FrequencySeries(model->roq->hctildeQuadratic);

	return(loglikelihood); /* The ROQ isn't compatible with the stuff below, so we can just exit here */



  }

  // for models which are non-factorising
  switch(marginalisationflags)
  {
    case MARGPHI:
    {
      REAL8 R = 2.0*cabs(Rcplx);
      REAL8 phase_maxL = carg(Rcplx);
      if(phase_maxL<0.0) phase_maxL=LAL_TWOPI+phase_maxL;
      LALInferenceAddVariable(currentParams,"phase_maxl",&phase_maxL,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	  if(LALInferenceCheckVariable(currentParams,"phase")) LALInferenceRemoveVariable(currentParams,"phase");
      gsl_sf_result result;
      REAL8 I0x=0.0;
      if(GSL_SUCCESS==gsl_sf_bessel_I0_scaled_e(R, &result))
      {
        I0x=result.val;
      }
      else printf("ERROR: Cannot calculate I0(%lf)\n",R);
      /* This is marginalised over phase only for now */
      loglikelihood += -(S+D) + log(I0x) + R ;
      d_inner_h= 0.5*R;
      break;
    }
    case GAUSSIAN:
    {
      d_inner_h = creal(Rcplx);
      break;
    }
    case MARGTIMEPHI:
    case MARGTIME:
    {
      /* LALSuite only performs complex->real reverse-FFTs. */
      dh_S_tilde->data[0] = crect( creal(dh_S_tilde->data[0]), 0. );

      XLALREAL8ReverseFFT(dh_S, dh_S_tilde, data->margFFTPlan);

      if (margphi) {
          dh_S_phase_tilde->data[0] = crect( creal(dh_S_phase_tilde->data[0]), 0.0);
          XLALREAL8ReverseFFT(dh_S_phase, dh_S_phase_tilde, data->margFFTPlan);
      }

      REAL8 time_low,time_high;
      LALInferenceGetMinMaxPrior(currentParams,"time",&time_low,&time_high);
      t0 = XLALGPSGetREAL8(&(data->freqData->epoch));
      int istart = (UINT4)round((time_low - t0)/deltaT);
      int iend = (UINT4)round((time_high - t0)/deltaT);
      if(iend > (int) dh_S->length || istart < 0 ) {
              fprintf(stderr,"ERROR: integration over time extends past end of buffer! Is your time prior too wide?\n");
              exit(1);
      }
      UINT4 n = iend - istart;
      REAL8 xMax = -1.0;
      REAL8 angMax = 0.0;
      if (margphi) {
          /* We've got the real and imaginary parts of the FFT in the two
             arrays.  Now combine them into one Bessel function. */
          for (i = istart; i < iend; i++) {
              /* Note: No factor of 2 for x because the 2-sided FFT above introduces that for us */
              double x = sqrt(dh_S->data[i]*dh_S->data[i] + dh_S_phase->data[i]*dh_S_phase->data[i]);
              if (x > xMax) { /* Store the phase angle at max L */
                  angMax = atan2(dh_S_phase->data[i], dh_S->data[i]);
                  xMax=x;
              }
              double I0=log(gsl_sf_bessel_I0_scaled(x)) + fabs(x);
              dh_S->data[i] = I0;
          }
      }
      size_t imax;
      REAL8 imean;
      loglikelihood += integrate_interpolated_log(deltaT, dh_S->data + istart, n, &imean, &imax) - log(n*deltaT);

      REAL8 max_time=t0+((REAL8) imax + istart)*deltaT;
      REAL8 mean_time=t0+(imean+(double)istart)*deltaT;

      if(margphi){
        REAL8 phase_maxL=angMax;
        if(phase_maxL<0.0) phase_maxL=LAL_TWOPI+phase_maxL;
        LALInferenceAddVariable(currentParams,"phase_maxl",&phase_maxL,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	    if(LALInferenceCheckVariable(currentParams,"phase")) LALInferenceRemoveVariable(currentParams,"phase");
        d_inner_h= 0.5*xMax;
      }
      else
      {
        d_inner_h=0.5*dh_S->data[imax+istart];
      }
      LALInferenceAddVariable(currentParams,"time_maxl",&max_time,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      LALInferenceAddVariable(currentParams,"time_mean",&mean_time,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
      XLALDestroyCOMPLEX16Vector(dh_S_tilde);
      XLALDestroyREAL8Vector(dh_S);
      if (margphi) {
        XLALDestroyCOMPLEX16Vector(dh_S_phase_tilde);
        XLALDestroyREAL8Vector(dh_S_phase);
      }
      break;
    }
    default:
      break;

  }
  /* SNR variables */
  REAL8 OptimalSNR=sqrt(2.0*S);
  REAL8 MatchedFilterSNR = 0.;

  /* Avoid nan's, since noise-only model has OptimalSNR == 0. */
  if (OptimalSNR > 0.)
      MatchedFilterSNR = 2.0*d_inner_h/OptimalSNR;
  LALInferenceAddVariable(currentParams,"optimal_snr",&OptimalSNR,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  LALInferenceAddVariable(currentParams,"matched_filter_snr",&MatchedFilterSNR,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);

  //loglikelihood = -1.0 * chisquared; // note (again): the log-likelihood is unnormalised!

  return(loglikelihood);
}

/***************************************************************/
/* Student-t (log-) likelihood function                        */
/* as described in Roever/Meyer/Christensen (2011):            */
/*   "Modelling coloured residual noise                        */
/*   in gravitational-wave signal processing."                 */
/*   Classical and Quantum Gravity, 28(1):015010.              */
/*   http://dx.doi.org/10.1088/0264-9381/28/1/015010           */
/*   http://arxiv.org/abs/0804.3853                            */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "time"            (REAL8, GPS sec.)                     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This function is essentially the same as the                */
/* "UndecomposedFreqDomainLogLikelihood()" function.           */
/* The additional parameter to be supplied is the (REAL8)      */
/* degrees-of-freedom parameter (nu) for each Ifo.             */
/* The additional "df" argument gives the corresponding        */
/* d.f. parameter for each element of the "*data" list.        */
/* The names of "df" must match the "->name" slot of           */
/* the elements of "data".                                     */
/*                                                             */
/* (TODO: allow for d.f. parameter to vary with frequency,     */
/*        i.e., to be a set of vectors corresponding to        */
/*        frequencies)                                         */
/***************************************************************/

REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams,
                                                    LALInferenceIFOData *data,
                                                    LALInferenceModel *model)
{

  return LALInferenceFusedFreqDomainLogLikelihood(currentParams, data, model, STUDENTT);

}


REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * dataPtr,
                                                COMPLEX16Vector * freqData1,
                                                COMPLEX16Vector * freqData2)
{
  if (dataPtr==NULL || freqData1 ==NULL || freqData2==NULL){
  	XLAL_ERROR_REAL8(XLAL_EFAULT);
  	}

  int lower, upper, i;
  double deltaT, deltaF;

  double overlap=0.0;

  /* determine frequency range & loop over frequency bins: */
  deltaT = dataPtr->timeData->deltaT;
  deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
  lower = ceil(dataPtr->fLow / deltaF);
  upper = floor(dataPtr->fHigh / deltaF);

  for (i=lower; i<=upper; ++i){
    overlap  += ((4.0*deltaF*(creal(freqData1->data[i])*creal(freqData2->data[i])+cimag(freqData1->data[i])*cimag(freqData2->data[i])))
                 / dataPtr->oneSidedNoisePowerSpectrum->data->data[i]);
  }

  return overlap;
}

REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data)
/*Identical to FreqDomainNullLogLikelihood                        */
{
	REAL8 loglikelihood, totalChiSquared=0.0;
	LALInferenceIFOData *ifoPtr=data;

	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
          REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	loglikelihood = -0.5 * totalChiSquared; // note (again): the log-likelihood is unnormalised!
	return(loglikelihood);
}

REAL8 LALInferenceMarginalisedPhaseLogLikelihood(LALInferenceVariables *currentParams,
                                                    LALInferenceIFOData *data,
                                                    LALInferenceModel *model)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Analytically marginalised over phase and distance           */
/* See LIGO-T1300326 for details                               */
/* At a distance of 1 Mpc for phi_0=0                          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  return LALInferenceFusedFreqDomainLogLikelihood(currentParams, data, model, MARGPHI);
}

/** Integrate interpolated log, returns the mean index in *imax if it
 * is not a NULL pointer.  Stores the mean index in *imean (can be
 * fractional).
 *
 * The method used is the trapezoid method, which is quadratically
 * accurate.
 */
static double integrate_interpolated_log(double h, REAL8 *log_ys, size_t n, double *imean, size_t *imax) {
  size_t i;
  double log_integral = -INFINITY;
  double max=-INFINITY;
  size_t imax_l=0;
  double log_imean_l=-INFINITY;
  double log_h = log(h);

  for (i = 1; i < n-1; i++) {
    log_integral = logaddexp(log_integral, log_ys[i]);
    log_imean_l = logaddexp(log_imean_l, log(i) + log_ys[i]);

    if (log_ys[i] > max) {
      max = log_ys[i];
      imax_l = i;
    }
  }

  log_integral = logaddexp(log_integral, log(0.5) + log_ys[0]);
  log_integral = logaddexp(log_integral, log(0.5) + log_ys[n-1]);

  /* No contribution to mean index from i = 0 term! */
  log_imean_l = logaddexp(log_imean_l, log(0.5) + log(n-1) + log_ys[n-1]);

  log_integral += log_h;
  log_imean_l += log_h;

  if (creal(log_ys[0]) > max) {
    max = log_ys[0];
    imax_l = 0;
  }

  if (log_ys[n-1] > max) {
    max = log_ys[n-1];
    imax_l = n-1;
  }

  log_imean_l -= log_integral;

  if(imean) *imean=exp(log_imean_l-log_integral);
  if(imax) *imax=imax_l;

  return log_integral;
}

REAL8 LALInferenceMarginalisedTimeLogLikelihood(LALInferenceVariables *currentParams,
                                                LALInferenceIFOData *data,
                                                LALInferenceModel *model)
{

  return ( LALInferenceFusedFreqDomainLogLikelihood(currentParams,data,model,MARGTIME));


}

REAL8 LALInferenceMarginalisedTimePhaseLogLikelihood(LALInferenceVariables *currentParams,
                                                LALInferenceIFOData *data,
                                                LALInferenceModel *model)
{

  return ( LALInferenceFusedFreqDomainLogLikelihood(currentParams,data,model,MARGTIMEPHI));


}

REAL8 LALInferenceFastSineGaussianLogLikelihood(LALInferenceVariables *currentParams,
                                                        LALInferenceIFOData *data,
                                                        LALInferenceModel *model)
/***************************************************************/
/* This is basically a loglikelihood for LIB that does not do  */
/* any check for extra options (marginalization, calibration,  */
/* As such, it is slighly faster. Don't use if you don't know  */
/* what you are doing                                          */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  double Fplus, Fcross;
  int i, lower, upper, ifo;
  LALInferenceIFOData *dataPtr;
  double ra=0.0, dec=0.0, psi=0.0, gmst=0.0;
  double GPSdouble=0.0;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift=0;  /* time shift (not necessarily same as above)                   */
  double deltaT, TwoDeltaToverN, deltaF, twopit=0.0, re, im, dre, dim, newRe, newIm;
  double timeTmp;
  /* Burst templates are generated at hrss=1, thus need to rescale amplitude */
  double amp_prefactor=1.0;

  REAL8 chisq=0.0;
  LALStatus status;
  memset(&status,0,sizeof(status));

  if(data==NULL) {XLAL_ERROR_REAL8(XLAL_EINVAL,"ERROR: Encountered NULL data pointer in likelihood\n");}

  int Nifos=0;
  for(dataPtr=data;dataPtr;dataPtr=dataPtr->next) Nifos++;
  void *generatedFreqModels[1+Nifos];
  for(i=0;i<=Nifos;i++) generatedFreqModels[i]=NULL;

  if(LALInferenceCheckVariable(currentParams, "loghrss")){
    amp_prefactor = exp(*(REAL8*)LALInferenceGetVariable(currentParams,"loghrss"));
  }
  else if (LALInferenceCheckVariable(currentParams, "hrss")){
    amp_prefactor = (*(REAL8*)LALInferenceGetVariable(currentParams,"hrss"));
  }

  /* determine source's sky location & orientation parameters: */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */

  INT4 SKY_FRAME=0;
  if(LALInferenceCheckVariable(currentParams,"SKY_FRAME"))
    SKY_FRAME=*(INT4 *)LALInferenceGetVariable(currentParams,"SKY_FRAME");
  if(SKY_FRAME==0){
    /* determine source's sky location & orientation parameters: */
    ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
    dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  }
  else
  {
    if(Nifos<2){
      fprintf(stderr,"ERROR: Cannot use --detector-frame with less than 2 detectors!\n");
      exit(1);
    }
    REAL8 t0=LALInferenceGetREAL8Variable(currentParams,"t0");
    REAL8 alph=acos(LALInferenceGetREAL8Variable(currentParams,"cosalpha"));
    REAL8 theta=LALInferenceGetREAL8Variable(currentParams,"azimuth");
    LALInferenceDetFrameToEquatorial(data->detector,data->next->detector,
                                   t0,alph,theta,&GPSdouble,&ra,&dec);
    LALInferenceAddVariable(currentParams,"rightascension",&ra,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(currentParams,"declination",&dec,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(currentParams,"time",&GPSdouble,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  deltaT = data->timeData->deltaT;
  /* figure out GMST: */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  chisquared = 0.0;
  REAL8 loglikelihood = 0.0;

  /* Reset SNR */
  model->SNR = 0.0;

  /* loop over data (different interferometers): */
  for(dataPtr=data,ifo=0; dataPtr; dataPtr=dataPtr->next,ifo++) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* model->params (set, e.g., from the trigger value).     */

    /* Reset log-likelihood */
    model->ifo_loglikelihoods[ifo] = 0.0;
    model->ifo_SNRs[ifo] = 0.0;

    /* Check to see if this buffer has already been filled with the signal.
    Different dataPtrs can share the same signal buffer to avoid repeated
    calls to template */
      if(!checkItemAndAdd((void *)(model->freqhPlus), generatedFreqModels))
    {
      /* Compare parameter values with parameter values corresponding  */
      /* to currently stored template; ignore "time" variable:         */
      if (LALInferenceCheckVariable(model->params, "time")) {
          timeTmp = *(REAL8 *) LALInferenceGetVariable(model->params, "time");
          LALInferenceRemoveVariable(model->params, "time");
      }
      else timeTmp = GPSdouble;

      LALInferenceCopyVariables(currentParams, model->params);
      // Remove time variable so it can be over-written (if it was pinned)
      if(LALInferenceCheckVariable(model->params,"time")) LALInferenceRemoveVariable(model->params,"time");
      LALInferenceAddVariable(model->params, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);

      INT4 errnum=0;
      XLAL_TRY(model->templt(model),errnum);
      errnum&=~XLAL_EFUNC;
      if(errnum!=XLAL_SUCCESS)
      {
        switch(errnum)
        {
          case XLAL_EUSR0: /* Template generation failed in a known way, set -Inf likelihood */
            return (-DBL_MAX);
            break;
          default: /* Panic! */
            fprintf(stderr,"Unhandled error in template generation - exiting!\n");
            fprintf(stderr,"XLALError: %d, %s\n",errnum,XLALErrorString(errnum));
            exit(1);
            break;
        }

      }

    }
        /* Template is now in model->timeFreqhPlus and hCross */

      /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
      XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, ra, dec, psi, gmst);
      /* signal arrival time (relative to geocenter); */
      timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location, ra, dec, &GPSlal);
      /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */
      /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
      timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(model->params, "time"))) + timedelay;
      twopit    = LAL_TWOPI * timeshift;

      /* For burst the effect of windowing in amplitude is important. Add it here. */
      Fplus*=amp_prefactor;
      Fcross*=amp_prefactor;

      dataPtr->fPlus = Fplus;
      dataPtr->fCross = Fcross;
      dataPtr->timeshift = timeshift;


      /* determine frequency range & loop over frequency bins: */
      deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
      lower = (UINT4)ceil(dataPtr->fLow / deltaF);
      upper = (UINT4)floor(dataPtr->fHigh / deltaF);
      TwoDeltaToverN = 2.0 * deltaT / ((double) dataPtr->timeData->data->length);

    /* Employ a trick here for avoiding cos(...) and sin(...) in time
       shifting.  We need to multiply each template frequency bin by
       exp(-J*twopit*deltaF*i) = exp(-J*twopit*deltaF*(i-1)) +
       exp(-J*twopit*deltaF*(i-1))*(exp(-J*twopit*deltaF) - 1) .  This
       recurrance relation has the advantage that the error growth is
       O(sqrt(N)) for N repetitions. */

    /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2*/
    dim = -sin(twopit*deltaF);
    dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);

    REAL8 *psd=&(dataPtr->oneSidedNoisePowerSpectrum->data->data[lower]);
    COMPLEX16 *dtilde=&(dataPtr->freqData->data->data[lower]);
    COMPLEX16 *hptilde=&(model->freqhPlus->data->data[lower]);
    COMPLEX16 *hctilde=&(model->freqhCross->data->data[lower]);
    COMPLEX16 diff=0.0;
    COMPLEX16 template=0.0;
    INT4 upppone=upper+1;
    for (i=lower,chisq=0.0,re = cos(twopit*deltaF*i),im = -sin(twopit*deltaF*i);
         i<upppone;
         i++, psd++, hptilde++, hctilde++, dtilde++,
         newRe = re + re*dre - im*dim,
         newIm = im + re*dim + im*dre,
         re = newRe, im = newIm)
    {

      COMPLEX16 d=*dtilde;
      /* Normalise PSD to our funny standard (see twoDeltaTOverN below). */
      REAL8 sigmasq=(*psd)*deltaT*deltaT;
      //subtract GW model from residual
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      COMPLEX16 plainTemplate = Fplus*(*hptilde)+Fcross*(*hctilde);

      /* Do time shifting */
      template = plainTemplate * (re + I*im);
      diff = (d - template);

      REAL8 diffsq = creal(diff)*creal(diff)+cimag(diff)*cimag(diff);
      chisq = TwoDeltaToverN*diffsq/sigmasq;
      chisquared  += chisq;
      model->ifo_loglikelihoods[ifo] -= chisq;
    } /* End loop over freq bins */

    loglikelihood += model->ifo_loglikelihoods[ifo];

  } /* end loop over detectors */
 // printf("%10.10e\n",loglikelihood);
  return(loglikelihood);
}


void LALInferenceNetworkSNR(LALInferenceVariables *currentParams,
                            LALInferenceIFOData *data,
                            LALInferenceModel *model)
/***************************************************************/
/* Calculate the SNR across the network.                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
  double Fplus, Fcross;
  REAL8 plainTemplateReal, plainTemplateImag;
  int i, lower, upper, ifo;
  LALInferenceIFOData *dataPtr;
  double ra=0.0, dec=0.0, psi=0.0, gmst=0.0;
  double GPSdouble=0.0;
  LIGOTimeGPS GPSlal;
  double deltaT, TwoOverNDeltaT, deltaF;
  double timeTmp;
  double mc;
  LALStatus status;
  memset(&status,0,sizeof(status));
  INT4 remove_time = 0;

  int signalFlag = 1;   //flag for including signal model

  int Nifos=0;
  for(dataPtr=data;dataPtr;dataPtr=dataPtr->next) Nifos++;
  void *generatedFreqModels[1+Nifos];
  for(i=0;i<=Nifos;i++) generatedFreqModels[i]=NULL;

  // If time isn't in current params, remove it at the end
  if(!LALInferenceCheckVariable(currentParams, "time"))
      remove_time = 1;

  //check if signal model is being used
  signalFlag=1;
  if(LALInferenceCheckVariable(currentParams, "signalModelFlag"))
    signalFlag = *((INT4 *)LALInferenceGetVariable(currentParams, "signalModelFlag"));

  /* Reset SNRs in model struct */
  model->SNR = 0.0;

  dataPtr = data;
  ifo = 0;
  while (dataPtr != NULL) {
      model->ifo_SNRs[ifo] = 0.0;
      ifo++;
      dataPtr = dataPtr->next;
  }

  if (!signalFlag)
      return;

  if(LALInferenceCheckVariable(currentParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
    LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }


  INT4 SKY_FRAME=0;
  if(LALInferenceCheckVariable(currentParams,"SKY_FRAME"))
    SKY_FRAME=*(INT4 *)LALInferenceGetVariable(currentParams,"SKY_FRAME");
  if(SKY_FRAME==0){
    /* determine source's sky location & orientation parameters: */
    ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
    dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  }
  else
  {
    if(Nifos<2){
      fprintf(stderr,"ERROR: Cannot use --detector-frame with less than 2 detectors!\n");
      exit(1);
    }
    REAL8 t0=LALInferenceGetREAL8Variable(currentParams,"t0");
    REAL8 alph=acos(LALInferenceGetREAL8Variable(currentParams,"cosalpha"));
    REAL8 theta=LALInferenceGetREAL8Variable(currentParams,"azimuth");
    LALInferenceDetFrameToEquatorial(data->detector,data->next->detector,
                                   t0,alph,theta,&GPSdouble,&ra,&dec);
    LALInferenceAddVariable(currentParams,"rightascension",&ra,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(currentParams,"declination",&dec,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(currentParams,"time",&GPSdouble,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */

  if (LALInferenceCheckVariable(currentParams,"time"))
      GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
  else {
      UINT4 freq_length = data->freqData->data->length;
      UINT4 time_length = 2*(freq_length-1);
      REAL8 epoch = XLALGPSGetREAL8(&(data->freqData->epoch));
      GPSdouble = epoch + (time_length-1)*data->timeData->deltaT - 2.0;
  }

  /* figure out GMST: */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);

  ifo=0;
  dataPtr = data;
  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself   */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* model->params (set, e.g., from the trigger value).     */

    /* Check to see if this buffer has already been filled with the signal.
     Different dataPtrs can share the same signal buffer to avoid repeated
     calls to template */
    if(!checkItemAndAdd((void *)(model->freqhPlus), generatedFreqModels))
    {
      /* to currently stored template; ignore "time" variable:         */
      if (LALInferenceCheckVariable(model->params, "time")) {
        timeTmp = *(REAL8 *) LALInferenceGetVariable(model->params, "time");
        LALInferenceRemoveVariable(model->params, "time");
      }
      else timeTmp = GPSdouble;

      LALInferenceCopyVariables(currentParams, model->params);
      // Remove time variable so it can be over-written (if it was pinned)
      if(LALInferenceCheckVariable(model->params,"time")) LALInferenceRemoveVariable(model->params,"time");
      LALInferenceAddVariable(model->params, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      if (!LALInferenceCheckVariable(model->params, "phase")) {
        double pi2 = M_PI / 2.0;
        LALInferenceAddVariable(model->params, "phase", &pi2, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
      }
      INT4 errnum=0;
      XLAL_TRY(model->templt(model),errnum);
      errnum&=~XLAL_EFUNC;
      if(errnum!=XLAL_SUCCESS)
      {
        switch(errnum)
        {
          case XLAL_EUSR0: /* Template generation failed in a known way, set -Inf likelihood */
            return;
            break;
          default: /* Panic! */
            fprintf(stderr,"Unhandled error in template generation - exiting!\n");
            fprintf(stderr,"XLALError: %d, %s\n",errnum,XLALErrorString(errnum));
            exit(1);
            break;
        }

      }

      if (model->domain == LAL_SIM_DOMAIN_TIME) {
        /* TD --> FD. */
        LALInferenceExecuteFT(model);
      }
    }

    /* Template is now in dataPtr->timeFreqModelhPlus and hCross */

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, (const REAL4(*)[3])dataPtr->detector->response, ra, dec, psi, gmst);

    dataPtr->fPlus = Fplus;
    dataPtr->fCross = Fcross;

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    TwoOverNDeltaT = 2.0 / (deltaT * ((double) dataPtr->timeData->data->length));

    for (i=lower; i<=upper; ++i){
      //subtract GW model from residual
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = Fplus * creal(model->freqhPlus->data->data[i])
                          +  Fcross * creal(model->freqhCross->data->data[i]);
      plainTemplateImag = Fplus * cimag(model->freqhPlus->data->data[i])
                          +  Fcross * cimag(model->freqhCross->data->data[i]);

      /* un-do 1/deltaT scaling: */
      model->ifo_SNRs[ifo] += 2.0 * TwoOverNDeltaT * ( plainTemplateReal*plainTemplateReal + plainTemplateImag*plainTemplateImag ) / dataPtr->oneSidedNoisePowerSpectrum->data->data[i];
    }

    model->SNR += model->ifo_SNRs[ifo];
    model->ifo_SNRs[ifo] = sqrt(model->ifo_SNRs[ifo]);

    ifo++; //increment IFO counter for noise parameters
    dataPtr = dataPtr->next;
  }

  if (remove_time && LALInferenceCheckVariable(currentParams, "time"))
    LALInferenceRemoveVariable(currentParams, "time");

  model->SNR = sqrt(model->SNR);
}
