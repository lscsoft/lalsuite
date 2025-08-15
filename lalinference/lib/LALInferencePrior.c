/*
 *  LALInferencePrior.c:  Nested Sampling using LALInference
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <lal/LALSimBurst.h>
#include <lal/XLALGSL.h>

#include "logaddexp.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Private helper function prototypes */
static double qInnerIntegrand(double M2, void *viData);
static double etaInnerIntegrand(double M2, void *viData);
static double outerIntegrand(double M1, void *voData);

static REAL8 REAL8max(REAL8 a, REAL8 b);
static REAL8 REAL8max(REAL8 a, REAL8 b)
{
  return (a>b?a:b);
}

void LALInferenceInitCBCPrior(LALInferenceRunState *runState)
{
    char help[]="\
    ----------------------------------------------\n\
    --- Prior Arguments --------------------------\n\
    ----------------------------------------------\n\
    (--distance-prior-uniform)       Impose uniform prior on distance and not volume (False)\n\
    (--distance-prior-comoving-volume) Impose uniform prior on source-frame comoving volume (False)\n\
    (--malmquistprior)               Impose selection effects on the prior (False)\n\
    (--malmquist-loudest-snr)        Threshold SNR in the loudest detector (0.0)\n\
    (--malmquist-second-loudest-snr) Threshold SNR in the second loudest detector (5.0)\n\
    (--malmquist-network-snr)        Threshold network SNR (0.0)\n\
    (--analyticnullprior)            Use analytic null prior\n\
    (--nullprior)                    Use null prior in the sampled parameters\n\
    (--alignedspin-zprior)           Use prior on z component of spin that corresponds to fully precessing model\n\
    (--spin-volumetricprior)         Use prior on spin components that is uniform inside the sphere\n\
    \n";
    ProcessParamsTable *ppt = NULL;

    /* Print command line arguments if help requested */
    if(runState == NULL || LALInferenceGetProcParamVal(runState->commandLine, "--help")) {
        fprintf(stdout, "%s", help);
        return;
    }
    ProcessParamsTable *commandLine=runState->commandLine;

    /* Choose the proper prior */
    if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
        runState->prior = &LALInferenceAnalyticNullPrior;
        runState->CubeToPrior = &LALInferenceAnalyticCubeToPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
        runState->prior = &LALInferenceNullPrior;
        /* CubeToPrior missing for null prior */
    } else {
        runState->prior = &LALInferenceInspiralPrior;
        runState->CubeToPrior = &LALInferenceInspiralCubeToPrior;

    }

    /* Optional uniform prior on distance */
    INT4 uniform_distance = 0;
    INT4 src_comove_volume_distance = 0;
    if (LALInferenceGetProcParamVal(commandLine, "--distance-prior-uniform")) {
      uniform_distance = 1;
    }
    if (LALInferenceGetProcParamVal(commandLine, "--distance-prior-comoving-volume")) {
      src_comove_volume_distance = 1;
    }
    LALInferenceAddVariable(runState->priorArgs,
                                "uniform_distance", &uniform_distance,
                                LALINFERENCE_INT4_t,
                                LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->priorArgs,
                                "src_comove_volume_distance", &src_comove_volume_distance,
                                LALINFERENCE_INT4_t,
                                LALINFERENCE_PARAM_OUTPUT);

    /* Set up malmquist prior */
    INT4 malmquist = 0;
    if (LALInferenceGetProcParamVal(commandLine, "--malmquistprior")) {
        malmquist = 1;
        REAL8 malmquist_loudest = 0.0;
        REAL8 malmquist_second_loudest = 5.0;
        REAL8 malmquist_network = 0.0;

        ppt = LALInferenceGetProcParamVal(commandLine,
                                            "--malmquist-loudest-snr");
        if (ppt) malmquist_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(commandLine,
                                            "--malmquist-second-loudest-snr");
        if (ppt) malmquist_second_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(commandLine,
                                            "--malmquist-network-snr");
        if (ppt) malmquist_network = atof(ppt->value);

        LALInferenceAddVariable(runState->priorArgs,
                                "malmquist", &malmquist,
                                LALINFERENCE_INT4_t,
                                LALINFERENCE_PARAM_OUTPUT);

        LALInferenceAddVariable(runState->priorArgs,
                                "malmquist_loudest_snr",
                                &malmquist_loudest,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);

        LALInferenceAddVariable(runState->priorArgs,
                                "malmquist_second_loudest_snr",
                                &malmquist_second_loudest,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);

        LALInferenceAddVariable(runState->priorArgs,
                                "malmquist_network_snr",
                                &malmquist_network,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);
    }
    INT4 one=1;
    if(LALInferenceGetProcParamVal(commandLine,"--alignedspin-zprior"))
    {
      LALInferenceAddVariable(runState->priorArgs,"projected_aligned_spin",&one,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }
    if(LALInferenceGetProcParamVal(commandLine,"--spin-volumetricprior"))
    {
      LALInferenceAddVariable(runState->priorArgs,"volumetric_spin",&one,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
    }
    if(LALInferenceGetProcParamVal(commandLine,"--alignedspin-zprior")&&LALInferenceGetProcParamVal(commandLine,"--spin-volumetricprior"))
    {
        fprintf(stderr,"Error: You cannot use both --alignedspin-zprior and --spin-volumetricprior\n");
        abort();
    }
    if(LALInferenceGetProcParamVal(commandLine,"--distance-prior-comoving-volume")&&LALInferenceGetProcParamVal(commandLine,"--distance-prior-uniform"))
    {
        fprintf(stderr,"Error: You cannot use both --distance-prior-comoving-volume and --distance-prior-uniform\n");
        abort();
    }
}

void LALInferenceInitLIBPrior(LALInferenceRunState *runState)
{
    /*Call CBC prior first, in case CBC approx is used, then check for burst approx and eventually overwrite runState->prior */
    LALInferenceInitCBCPrior(runState);
    /*LIB specific call in case of burst approximant */
    ProcessParamsTable *commandLine=runState->commandLine;
    ProcessParamsTable *ppt = NULL;
    if ((ppt=LALInferenceGetProcParamVal(commandLine,"--approx"))){
      if ((XLALCheckBurstApproximantFromString(ppt->value)))
      {
         /* Choose the proper prior */
         if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") ||
         LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
         LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
            runState->prior = &LALInferenceAnalyticNullPrior;
         } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
            runState->prior = &LALInferenceNullPrior;
         } else {
           runState->prior = &LALInferenceSineGaussianPrior;
         }
      }
    }
}

static REAL8 LALInferenceConstantCalibrationPrior(LALInferenceRunState *runState, LALInferenceVariables *params) {

  LALInferenceIFOData *ifo = NULL;
  REAL8 ampWidth = -1.0;
  REAL8 phaseWidth = -1.0;
  REAL8 logPrior = 0.0;

  if (runState->commandLine == NULL || (!LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp") &&
      !LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha")))
  {
    return logPrior;
  }

  ifo = runState->data;
  do {

    if((LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp"))){
      ampWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "constcal_amp_uncertainty");
      if (ampWidth>0){
	      char ampVarName[VARNAME_MAX];
	      REAL8 amp = 0.0;
	      if((VARNAME_MAX <= snprintf(ampVarName, VARNAME_MAX, "calamp_%s", ifo->name)))
          {
                fprintf(stderr,"variable name too long\n"); abort();
          }
	      amp = *(REAL8*)LALInferenceGetVariable(params, ampVarName);
	      logPrior += -0.5*log(2.0*M_PI) - log(ampWidth) - 0.5*amp*amp/ampWidth/ampWidth;
      }
    }
    if((LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha"))){
      phaseWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "constcal_phase_uncertainty");
      if (phaseWidth>0){
	      char phaseVarName[VARNAME_MAX];
	      REAL8 phase = 0.0;
	      if((VARNAME_MAX <= snprintf(phaseVarName, VARNAME_MAX, "calpha_%s", ifo->name)))
          {
            fprintf(stderr,"variable name too long\n"); abort();
          }
	      phase = *(REAL8 *)LALInferenceGetVariable(params, phaseVarName);
	      logPrior += -0.5*log(2.0*M_PI) - log(phaseWidth) - 0.5*phase*phase/phaseWidth/phaseWidth;
      }
    }

    ifo = ifo->next;
  } while (ifo);

  return logPrior;
}

UINT4 LALInferenceCubeToConstantCalibrationPrior(LALInferenceRunState *runState, LALInferenceVariables *params, INT4 *idx, double *Cube, void UNUSED *context)
{
  LALInferenceIFOData *ifo = NULL;
  REAL8 ampWidth = -1.0;
  REAL8 phaseWidth = -1.0;

  if (runState->commandLine == NULL || (!LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp") &&
      !LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha")))
  {
    return 1;
  }
  ifo = runState->data;
  do {

    if((LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalAmp")))
    {
      ampWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "constcal_amp_uncertainty");
      if (ampWidth>0)
      {
        char ampVarName[VARNAME_MAX];
        REAL8 amp = 0.0;
        if((VARNAME_MAX <= snprintf(ampVarName, VARNAME_MAX, "calamp_%s", ifo->name)))
        {
            fprintf(stderr,"variable name too long\n"); abort();
        }
        amp = LALInferenceCubeToGaussianPrior(Cube[(*idx)++], 0.0, ampWidth);
        LALInferenceSetVariable(params, ampVarName, &amp);
      }
    }
    if((LALInferenceGetProcParamVal(runState->commandLine, "--MarginalizeConstantCalPha")))
    {
      phaseWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "constcal_phase_uncertainty");
      if (phaseWidth>0)
      {
        char phaseVarName[VARNAME_MAX];
        REAL8 phase = 0.0;
        if((VARNAME_MAX <= snprintf(phaseVarName, VARNAME_MAX, "calpha_%s", ifo->name)))
        {
            fprintf(stderr,"variable name too long\n"); abort();
        }
        phase = LALInferenceCubeToGaussianPrior(Cube[(*idx)++], 0.0, phaseWidth);
        LALInferenceSetVariable(params, phaseVarName, &phase);
      }
    }

    ifo = ifo->next;
  } while (ifo);

  return 1;

}



/* Return the log Prior for the glitch amplitude */
REAL8 logGlitchAmplitudeDensity(REAL8 A, REAL8 Q, REAL8 f)
{
  REAL8 SNR;
  /*
  REAL8 PIterm = 0.5*LAL_2_SQRTPI*LAL_SQRT1_2;
  */
  REAL8 PIterm = LAL_2_SQRTPI*LAL_SQRT1_2;
  REAL8 SNRPEAK = 5.0;

  SNR = A*sqrt( (PIterm*Q/f) );
  return log(SNR/(SNRPEAK*SNRPEAK))+(-SNR/SNRPEAK);
}

/* Return the log Prior for the glitch model */
static REAL8 LALInferenceGlitchPrior(LALInferenceRunState *runState, LALInferenceVariables *params) {
  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 logPrior=0.0;

  /* check if glitch model is being used */
  UINT4 glitchFlag = 0;
  if(LALInferenceCheckVariable(params,"glitchFitFlag"))
    glitchFlag = *((INT4 *)LALInferenceGetVariable(params, "glitchFitFlag"));

  if(glitchFlag)
  {
    UINT4 nifo,nglitch;
    // Get parameters for current glitch model
    UINT4Vector *gsize   = *(UINT4Vector **) LALInferenceGetVariable(params, "glitch_size");

    gsl_matrix *glitch_f = *(gsl_matrix **)LALInferenceGetVariable(params, "morlet_f0");
    gsl_matrix *glitch_Q = *(gsl_matrix **)LALInferenceGetVariable(params, "morlet_Q");
    gsl_matrix *glitch_A = *(gsl_matrix **)LALInferenceGetVariable(params, "morlet_Amp");
    gsl_matrix *gparams = NULL;

    REAL8 component_min=0.0;
    REAL8 component_max=0.0;
    REAL8 val=0.0;

    char priormin[512];
    char priormax[512];

    REAL8 Anorm = *(REAL8 *)LALInferenceGetVariable(priorParams,"glitch_norm");

    REAL8 A,f,Q;
    for(nifo=0; nifo<(UINT4)gsize->length; nifo++)
    {
      for(nglitch=0; nglitch<gsize->data[nifo]; nglitch++)
      {
        A = gsl_matrix_get(glitch_A,nifo,nglitch);
        Q = gsl_matrix_get(glitch_Q,nifo,nglitch);
        f = gsl_matrix_get(glitch_f,nifo,nglitch);

        logPrior += logGlitchAmplitudeDensity(A*Anorm,Q,f);
        //printf("logPrior=%g\n",logPrior);
      }
    }
    for(;item;item=item->next){
      if(!strcmp(item->name,"morlet_f0" ) || !strcmp(item->name,"morlet_Q"  ) || !strcmp(item->name,"morlet_t0" ) || !strcmp(item->name,"morlet_phi") )
      {
        gparams = *((gsl_matrix **)(item->value));

        sprintf(priormin,"%s_prior_min",item->name);
        sprintf(priormax,"%s_prior_max",item->name);
        component_min=*(REAL8 *)LALInferenceGetVariable(priorParams,priormin);
        component_max=*(REAL8 *)LALInferenceGetVariable(priorParams,priormax);

        for(UINT4 i=0; i<gsize->length; i++)
        {
          for(UINT4 j=0; j<gsize->data[i]; j++)
          {
            val = gsl_matrix_get(gparams,i,j);

            //rejection sample on prior
            if(val<component_min || val>component_max) return -INFINITY;
            else logPrior -= log(component_max-component_min);
          }
        }
      }//end morlet parameters prior
    }//end loop through params
  }//end check for glitch parameters

  return logPrior;
}

/* Return the log Prior for the psd model */
static REAL8 LALInferencePSDPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 logPrior=0.0;


  /* check if PSD model is being used */
  UINT4 psdFlag = 0;
  if(LALInferenceCheckVariable(params, "psdScaleFlag"))
    psdFlag = *((INT4 *)LALInferenceGetVariable(params, "psdScaleFlag"));

  /* PSD scale parameters */
  if(psdFlag)
  {
    UINT4 i;
    UINT4 j;

    REAL8 val   = 0.0;
    REAL8 var   = 1.0;
    REAL8 mean  = 1.0;

    REAL8Vector *sigma = *((REAL8Vector **)LALInferenceGetVariable(priorParams, "psdsigma"));
    gsl_matrix *nparams = *((gsl_matrix **)LALInferenceGetVariable(params, "psdscale"));

    REAL8 component_min=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_min");
    REAL8 component_max=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_max");

    UINT4 psdGaussianPrior = *(UINT4 *)LALInferenceGetVariable(priorParams,"psdGaussianPrior");

    //Loop over IFOs
    for(i=0; i<(UINT4)nparams->size1; i++)
    {
      //Loop over PSD windows
      for(j=0; j<(UINT4)nparams->size2; j++)
      {
        var = sigma->data[j]*sigma->data[j];
        val = gsl_matrix_get(nparams,i,j);

        //reject prior
        if(val < component_min || val > component_max) return -INFINITY;
        else if(psdGaussianPrior) logPrior += -0.5*( (mean-val)*(mean-val)/var + log(2.0*LAL_PI*var) );
      }//end loop over windows
    }//end loop over IFOs

  }//end psdFlag conditional

  return logPrior;
}


/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model)
{
  if (runState == NULL || runState->priorArgs == NULL || params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT, "Null arguments received.");

  REAL8 logPrior=0.0;

  LALInferenceVariableItem *item=NULL;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=-INFINITY, max=INFINITY;
  REAL8 mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;
  REAL8 c0=1.012306, c1=1.136740, c2=0.262462, c3=0.016732, c4=0.000387; /* fitting coefficients for Will's cosmological distance prior, see https://git.ligo.org/RatesAndPopulations/lalinfsamplereweighting/blob/master/ApproxPrior.ipynb */

  /* check if signal model is being used */
  UINT4 signalFlag=1;
  if(LALInferenceCheckVariable(params, "signalModelFlag"))
    signalFlag = *((INT4 *)LALInferenceGetVariable(params, "signalModelFlag"));

  if(signalFlag){

  /* Check boundaries for signal model parameters */
  for(item=params->head;item;item=item->next)
  {
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else if (LALInferenceCheckMinMaxPrior(priorParams, item->name))
    {
      if(item->type==LALINFERENCE_REAL8_t){
        LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
        if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -INFINITY;
      }
    }
    else if (LALInferenceCheckGaussianPrior(priorParams, item->name))
    {
      if(item->type==LALINFERENCE_REAL8_t){
	REAL8 mean,stdev,val;
	val = *(REAL8 *)item->value;
	LALInferenceGetGaussianPrior(priorParams, item->name, &mean, &stdev);
	logPrior+= -0.5*(mean-val)*(mean-val)/stdev/stdev - 0.5*log(LAL_TWOPI) - log(stdev);
      }
    }
  }
  if(LALInferenceCheckVariable(params, "flow") &&
          LALInferenceCheckVariableNonFixed(params, "flow")) {
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"flow"));
  }


  if(LALInferenceCheckVariable(params,"logdistance"))
  {
    REAL8 log_dist = *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
    if ((LALInferenceCheckVariable(priorParams,"uniform_distance") && LALInferenceGetINT4Variable(priorParams,"uniform_distance"))) {
      logPrior+=log_dist;
    }
    else if ((LALInferenceCheckVariable(priorParams,"src_comove_volume_distance") && LALInferenceGetINT4Variable(priorParams,"src_comove_volume_distance"))) {
	  REAL8 dist_Gpc = exp(log_dist)/1000.0;
      REAL8 dist_Gpc2= dist_Gpc*dist_Gpc;
      REAL8 dist_Gpc3= dist_Gpc2*dist_Gpc;
      REAL8 dist_Gpc4= dist_Gpc3*dist_Gpc;
      REAL8 denominator = c0+c1*dist_Gpc+c2*dist_Gpc2+c3*dist_Gpc3+c4*dist_Gpc4;
	  logPrior+=3.0* log_dist-log(denominator);
	}
	else {
      logPrior+=3.0* log_dist;
    }
  }
  else if(LALInferenceCheckVariable(params,"distance"))
  {
    if (!(LALInferenceCheckVariable(priorParams,"uniform_distance")&&LALInferenceGetINT4Variable(priorParams,"uniform_distance"))) {
      REAL8 dist = *(REAL8 *)LALInferenceGetVariable(params,"distance");
      if ((LALInferenceCheckVariable(priorParams,"src_comove_volume_distance") && LALInferenceGetINT4Variable(priorParams,"src_comove_volume_distance"))) {
        REAL8 dist_Gpc = dist/1000.0;
        REAL8 dist_Gpc2= dist_Gpc*dist_Gpc;
        REAL8 dist_Gpc3= dist_Gpc2*dist_Gpc;
        REAL8 dist_Gpc4= dist_Gpc3*dist_Gpc;
        REAL8 denominator = c0+c1*dist_Gpc+c2*dist_Gpc2+c3*dist_Gpc3+c4*dist_Gpc4;
        logPrior+=2.0*log(dist)-log(denominator);
      }
	  else {
        logPrior+=2.0*log(dist);
      }
    }
  }
  if(LALInferenceCheckVariable(params,"declination"))
  {
    /* Check that this is not an output variable */
    if(LALInferenceGetVariableVaryType(params,"declination")==LALINFERENCE_PARAM_LINEAR)
      logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
  }


  if(LALInferenceCheckVariable(params,"logmc")) {
    mc=exp(*(REAL8 *)LALInferenceGetVariable(params,"logmc"));
  } else if(LALInferenceCheckVariable(params,"chirpmass")) {
    mc=(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass"));
  }

  if(LALInferenceCheckVariable(params,"q")) {
    q=*(REAL8 *)LALInferenceGetVariable(params,"q");
    LALInferenceMcQ2Masses(mc,q,&m1,&m2);
  } else if(LALInferenceCheckVariable(params,"eta")) {
    eta=*(REAL8 *)LALInferenceGetVariable(params,"eta");
    LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
  }

  if(LALInferenceCheckVariable(params,"logmc")) {
    if(LALInferenceCheckVariable(params,"q"))
      logPrior+=log(m1*m1);
    else
      logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
  } else if(LALInferenceCheckVariable(params,"chirpmass")) {
    if(LALInferenceCheckVariable(params,"q"))
      logPrior+=log(m1*m1/mc);
    else
      logPrior+=log(((m1+m2)*(m1+m2))/((m1-m2)*pow(eta,3.0/5.0)));
  }

  /* Check for individual mass priors */
  if(LALInferenceCheckVariable(priorParams,"mass1_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass1_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass2_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass2_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
				  return -INFINITY;


  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      return -INFINITY;

  if(LALInferenceCheckVariable(priorParams,"MTotMin"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
      return -INFINITY;

  if(model != NULL &&
        LALInferenceCheckVariable(priorParams,"malmquist") &&
        *(UINT4 *)LALInferenceGetVariable(priorParams,"malmquist") &&
        !within_malmquist(runState, params, model))
      return -INFINITY;

  UINT4 volumetric_spins = LALInferenceCheckVariable(runState->priorArgs,"volumetric_spin") && LALInferenceGetVariable(runState->priorArgs,"volumetric_spin");
  /* Apply spin priors for precessing case */
  if(LALInferenceCheckVariable(params,"tilt_spin1"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"tilt_spin1");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      if(volumetric_spins)
      {
              /* homogenous inside spin bound */
              /* V = (4/3)*pi*(a_max^3 - a_min^3) */
              REAL8 a = LALInferenceGetREAL8Variable(params,"a_spin1");
              REAL8 a_max,a_min;
              LALInferenceGetMinMaxPrior(runState->priorArgs,"a_spin1",&a_min,&a_max);
              REAL8 V = (4./3.)*LAL_PI * (a_max*a_max*a_max - a_min*a_min*a_min);
              logPrior+=log(fabs(a*a))-log(fabs(V));
      }
      /* Usual case has uniform in a, but both cases have sin(tilt) from volume element */
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin1"))));
    }
  }
  else
  {
     if(volumetric_spins)
     {
            /* Volumetric prior marginalised onto z component */
              REAL8 a = LALInferenceGetREAL8Variable(params,"a_spin1");
              REAL8 a_max,a_min;
              LALInferenceGetMinMaxPrior(runState->priorArgs,"a_spin1",&a_min,&a_max);
              REAL8 V = (4./3.)*LAL_PI * (a_max*a_max*a_max);
              logPrior+=log(fabs((3./4.)*(a_max*a_max - a*a)))-log(fabs(V));
     }
  }
  if(LALInferenceCheckVariable(params,"tilt_spin2"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"tilt_spin2");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      if(volumetric_spins)
      {
              REAL8 a = LALInferenceGetREAL8Variable(params,"a_spin2");
              REAL8 a_max,a_min;
              LALInferenceGetMinMaxPrior(runState->priorArgs,"a_spin2",&a_min,&a_max);
              REAL8 V = (4./3.)*LAL_PI * (a_max*a_max*a_max - a_min*a_min*a_min);
              logPrior+=log(fabs(a*a))-log(fabs(V));
      }
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin2"))));
    }
  }
  else
  {
     if(volumetric_spins)
     {
            /* Volumetric prior marginalised onto z component */
              REAL8 a = LALInferenceGetREAL8Variable(params,"a_spin2");
              REAL8 a_max,a_min;
              LALInferenceGetMinMaxPrior(runState->priorArgs,"a_spin2",&a_min,&a_max);
              REAL8 V = (4./3.)*LAL_PI * (a_max*a_max*a_max);
              logPrior+=log(fabs((3./4.)*(a_max*a_max - a*a)))-log(fabs(V));
     }
  }

  /* Optional prior on aligned spin component that corresponds to the effective prior on
   that component when using a precessing spin model. p(z) = (1/2)(1/R)log(|z|/R)
   Where R is the maximum magnitude of the spin vector max(|a_spin1_max|,|a_spin1_min|).
   */
  if (LALInferenceCheckVariable(priorParams,"projected_aligned_spin") && LALInferenceGetVariable(priorParams,"projected_aligned_spin"))
  {
    REAL8 z=0.0;
    /* Double-check for tilts to prevent accidental double-prior */
    if(LALInferenceCheckVariable(params,"a_spin1") && !LALInferenceCheckVariable(params,"tilt_spin1"))
    {
      REAL8 R = REAL8max(fabs(LALInferenceGetREAL8Variable(priorParams,"a_spin1_max")),fabs(LALInferenceGetREAL8Variable(priorParams,"a_spin1_min")));
      z=LALInferenceGetREAL8Variable(params,"a_spin1");
      logPrior += -log(2.0) - log(R) + log(-log(fabs(z) / R));
    }
    if(LALInferenceCheckVariable(params,"a_spin2")&& !LALInferenceCheckVariable(params,"tilt_spin2"))
    {
      REAL8 R = REAL8max(fabs(LALInferenceGetREAL8Variable(priorParams,"a_spin2_max")),fabs(LALInferenceGetREAL8Variable(priorParams,"a_spin2_min")));
      z=LALInferenceGetREAL8Variable(params,"a_spin2");
      logPrior += -log(2.0) - log(R) + log(-log(fabs(z) / R));
    }

  }

  // Apply additional prior if using eos parameters
  if((LALInferenceCheckVariable(params,"logp1")&&LALInferenceCheckVariable(params,"gamma1")&&LALInferenceCheckVariable(params,"gamma2")&&LALInferenceCheckVariable(params,"gamma3")))
  {
    /*If EOS params and masses are aphysical, return -INFINITY to ensure point is rejected*/
    if(LALInferenceEOSPhysicalCheck(params,runState->commandLine)==XLAL_FAILURE){
       return -INFINITY;
    }
  }
  else if((LALInferenceCheckVariable(params,"SDgamma0")&&LALInferenceCheckVariable(params,"SDgamma1")&&LALInferenceCheckVariable(params,"SDgamma2")&&LALInferenceCheckVariable(params,"SDgamma3")))
  {
    /*If EOS params and masses are aphysical, return -INFINITY to ensure point is rejected*/
    if(LALInferenceEOSPhysicalCheck(params,runState->commandLine)==XLAL_FAILURE){
       return -INFINITY;
    }
  }

  }/* end prior for signal model parameters */


  /* Calibration priors. */
  /* Disabled as this is now handled automatically */
  //logPrior += LALInferenceSplineCalibrationPrior(runState, params);
  logPrior += LALInferenceConstantCalibrationPrior(runState, params);
  /* Evaluate PSD prior (returns 0 if no PSD model) */
  logPrior += LALInferencePSDPrior(runState, params);

  /* Evaluate glitch prior (returns 0 if no glitch model) */
  logPrior += LALInferenceGlitchPrior(runState, params);

  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context)
{
    REAL8 min=-INFINITY, max=INFINITY;
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *timeID = &info[2][0];
    int i = 0;

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(params,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(params,"SKY_FRAME");
    double azimuth, cosalpha, lat, longitude, t0, tc;

    if(SKY_FRAME==1)
    {
      // detector frame azimuth
      if (LALInferenceCheckVariable(params, "azimuth"))
      {
        item = LALInferenceGetItem(params, "azimuth");
        if (item->vary != LALINFERENCE_PARAM_FIXED)
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, "azimuth", (void *)&min, (void *)&max);
          azimuth = LALInferenceCubeToFlatPrior(Cube[i], min, max);
          LALInferenceSetVariable(params, "azimuth", &azimuth);
          i++;
        }
      }

      // detector frame cosalpha
      if (LALInferenceCheckVariable(params, "cosalpha"))
      {
        item = LALInferenceGetItem(params, "cosalpha");
        if (item->vary != LALINFERENCE_PARAM_FIXED)
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, "cosalpha", (void *)&min, (void *)&max);
          cosalpha = LALInferenceCubeToFlatPrior(Cube[i], min, max);
          LALInferenceSetVariable(params, "cosalpha", &cosalpha);
          i++;
        }
      }
    }
    else
    {
      // latitude
      if (LALInferenceCheckVariable(params, "declination"))
      {
        item = LALInferenceGetItem(params, "declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            lat = asin(2.0 * Cube[i] - 1.0);
            LALInferenceSetVariable(params, "declination", &lat);
            i++;
        }
      }

      // longitude
      if (LALInferenceCheckVariable(params, "rightascension"))
      {
        item = LALInferenceGetItem(params, "rightascension");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "rightascension", (void *)&min, (void *)&max);
            longitude = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &longitude);
            i++;
        }
      }
    }

    // phi
    if (LALInferenceCheckVariable(params, "phase"))
    {
        item = LALInferenceGetItem(params, "phase");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phase", (void *)&min, (void *)&max);
            double phi = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phase", &phi);
            i++;
        }
    }

    // psi
    if (LALInferenceCheckVariable(params, "polarisation"))
    {
        item = LALInferenceGetItem(params, "polarisation");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "polarisation", (void *)&min, (void *)&max);
            double psi = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "polarisation", &psi);
            i++;
        }
    }

    if(SKY_FRAME==1)
    {
      // detector frame t0
      if (LALInferenceCheckVariable(params, "t0"))
      {
        item = LALInferenceGetItem(params, "t0");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "t0", (void *)&min, (void *)&max);
            t0 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "t0", &t0);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }
    else
    {
      // time
      if (LALInferenceCheckVariable(params, "time"))
      {
        item = LALInferenceGetItem(params, "time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }


    // mass variables
    double mc = 0.0, eta = 0.0, q = 0.0, m1 = 0.0, m2 = 0.0, m = 0.0;

    // check if mchirp is fixed
    if( LALInferenceCheckVariable(params,"logmc") )
    {
        item = LALInferenceGetItem(params, "logmc");
        if(item->vary == LALINFERENCE_PARAM_FIXED) mc = exp(*(REAL8 *)LALInferenceGetVariable(params, "logmc"));
    }
    else if( LALInferenceCheckVariable(params,"chirpmass") )
    {
        item = LALInferenceGetItem(params, "chirpmass");
        if(item->vary == LALINFERENCE_PARAM_FIXED) mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
    }

    // check if eta is fixed
    if( LALInferenceCheckVariable(params,"eta") )
    {
        item = LALInferenceGetItem(params, "eta");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "eta");
            if( mc != 0.0 ) LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        }
    }
    else if( LALInferenceCheckVariable(params,"q") )
    {
        item = LALInferenceGetItem(params, "q");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "q");
            if( mc != 0.0 ) LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        }
    }

    //m1 & m2
    if( m1 == 0.0 && m2 == 0.0 )
    {
        REAL8 m1_min = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass1_min");
        REAL8 m2_min = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass2_min");
        REAL8 m1_max = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass1_max");
        REAL8 m2_max = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass2_max");

        if( m1_min == m1_max && m2_min==m2_max)
        {
          m1 = m1_min;
          m2 = m2_min;
          m = m1 + m2;
          eta = m1 * m2 / (m*m);
          mc = pow(eta,0.6) * m;
          q = m2 / m1; // asymmetric mass ratio, m1 >= m2
        }
        else
        {
            m1 = LALInferenceCubeToFlatPrior(Cube[i], m1_min, m1_max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i+1], m2_min, m2_max);
            if(m1<m2)
            {
                double temp = m2;
                m2 = m1;
                m1 = temp;
            }

            m = m1 + m2;
            eta = m1 * m2 / (m*m);
            mc = pow(eta,0.6) * m;
            q = m2 / m1; // asymmetric mass ratio, m1 >= m2
            i++;
            i++;
        }

        // chirp mass and eta/q
        if(LALInferenceCheckVariable(params,"eta")||LALInferenceCheckVariable(params,"q"))
        {
            if(LALInferenceCheckVariable(params,"logmc"))
            {
                double logmc = log(mc);
                LALInferenceSetVariable(params, "logmc", &logmc);
            }
            else if(LALInferenceCheckVariable(params,"chirpmass"))
            {
                LALInferenceSetVariable(params, "chirpmass", &mc);
            }

                  if(LALInferenceCheckVariable(params,"q"))
                LALInferenceSetVariable(params, "q", &q);
            else if(LALInferenceCheckVariable(params,"eta"))
                    LALInferenceSetVariable(params, "eta", &eta);
        }
    }

    // distance
    double dist;
    if( LALInferenceCheckVariable(params,"logdistance") )
    {
        item = LALInferenceGetItem(params, "logdistance");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "logdistance", (void *)&min, (void *)&max);
            min = exp(min); max = exp(max);
            dist = LALInferenceCubeToPowerPrior(2.0, Cube[i], min, max);
            double logdist = log(dist);
            LALInferenceSetVariable(params, "logdistance", &logdist);
            i++;
        }
    }
    else if( LALInferenceCheckVariable(params,"distance") )
    {
        item = LALInferenceGetItem(params, "distance");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "distance", (void *)&min, (void *)&max);
            dist = LALInferenceCubeToPowerPrior(2.0, Cube[i], min, max);
            LALInferenceSetVariable(params, "distance", &dist);
            i++;
        }
    }

    // a_spin1
    if(LALInferenceCheckVariable(params,"a_spin1"))
    {
        item = LALInferenceGetItem(params, "a_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin1", &a_spin1);
            i++;
        }
    }
    // a_spin2
    if(LALInferenceCheckVariable(params,"a_spin2"))
    {
        item = LALInferenceGetItem(params, "a_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin2", &a_spin2);
            i++;
        }
    }

    // cos theta_JN for system-frame parameters
    if(LALInferenceCheckVariable(params,"costheta_jn"))
    {
        item = LALInferenceGetItem(params, "costheta_jn");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "costheta_jn", (void *)&min, (void *)&max);
            double costheta_JN = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "costheta_jn", &costheta_JN);
            i++;
        }
    }

    // phi_JL for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi_jl"))
    {
        item = LALInferenceGetItem(params, "phi_jl");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_jl", (void *)&min, (void *)&max);
            double phi_JL = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_jl", &phi_JL);
            i++;
        }
    }

    // tilt of spin 1 for system-frame parameters
    if(LALInferenceCheckVariable(params,"tilt_spin1"))
    {
        item = LALInferenceGetItem(params, "tilt_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "tilt_spin1", (void *)&min, (void *)&max);
            double tilt_spin1 = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "tilt_spin1", &tilt_spin1);
            i++;
        }
    }

    // tilt of spin 2 for system-frame parameters
    if(LALInferenceCheckVariable(params,"tilt_spin2"))
    {
        item = LALInferenceGetItem(params, "tilt_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "tilt_spin2", (void *)&min, (void *)&max);
            double tilt_spin2 = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "tilt_spin2", &tilt_spin2);
            i++;
        }
    }

    // phi12 for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi12"))
    {
        item = LALInferenceGetItem(params, "phi12");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi12", (void *)&min, (void *)&max);
            double phi12 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi12", &phi12);
            i++;
        }
    }

    UINT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);
    UINT4 ConstCalib = LALInferenceCubeToConstantCalibrationPrior(runState, params, &i, Cube, context);
    //UINT4 SplineCalib = LALInferenceCubeToSplineCalibrationPrior(runState, params, &i, Cube, context);

    /* Check boundaries */
    if (ScaleTest==0 || ConstCalib==0 /*|| SplineCalib==0*/) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else if (item->type != LALINFERENCE_gslMatrix_t && item->type != LALINFERENCE_REAL8Vector_t)
        {
            LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
            if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return 0 ;
        }
    }

    /* Check for real mass values */
    if(LALInferenceCheckVariable(params,"logmc"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"logmc")))
            return 0;
    if(LALInferenceCheckVariable(params,"chirpmass"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass")))
            return 0;
    if(LALInferenceCheckVariable(params,"eta"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"eta"))
           ||*(REAL8 *)LALInferenceGetVariable(params,"eta") < 0.0
           || *(REAL8 *)LALInferenceGetVariable(params,"eta") > 0.25)
            return 0;
    if(LALInferenceCheckVariable(params,"q"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"q"))
           ||*(REAL8 *)LALInferenceGetVariable(params,"q") < 0.0
           || *(REAL8 *)LALInferenceGetVariable(params,"q") > 1.0)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMin"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMax"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
            return 0;

	/* Check for individual mass priors */
	if(LALInferenceCheckVariable(priorParams,"mass1_min"))
			if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
					return 0;
	if(LALInferenceCheckVariable(priorParams,"mass1_max"))
			if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
					return 0;
	if(LALInferenceCheckVariable(priorParams,"mass2_min"))
			if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
					return 0;
	if(LALInferenceCheckVariable(priorParams,"mass2_max"))
			if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
					return 0;

    if(model != NULL &&
        LALInferenceCheckVariable(priorParams,"malmquist") &&
        *(UINT4 *)LALInferenceGetVariable(priorParams,"malmquist") &&
        !within_malmquist(runState, params, model))
      return 0;

    return 1;
}

void LALInferenceCyclicReflectiveBound(LALInferenceVariables *parameter,
                                       LALInferenceVariables *priorArgs){
  REAL8 val, offset, delta, min, max;
  UINT8 n;

  if (parameter == NULL || priorArgs == NULL)
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");

  /* Apply cyclic and reflective boundaries to parameter to bring it back
     within the prior */
  LALInferenceVariableItem *paraHead=NULL;

  for (paraHead=parameter->head; paraHead; paraHead=paraHead->next) {
    if( paraHead->vary==LALINFERENCE_PARAM_FIXED ||
        paraHead->vary==LALINFERENCE_PARAM_OUTPUT ||
        !LALInferenceCheckMinMaxPrior(priorArgs, paraHead->name) ) continue;

    LALInferenceGetMinMaxPrior(priorArgs,paraHead->name, &min, &max);

    /* Check that the minimum and maximum make sense. */
    if (min >= max) {
      XLAL_ERROR_VOID(XLAL_EINVAL, "Minimum %f for variable '%s' is not less than maximum %f.", min, paraHead->name, max);
    }
    delta = max - min;

    val = *(REAL8 *)paraHead->value;

    if (val == INFINITY)
      return;

    // Nothing to do if between bounds
    if ((val >= min) && (val <= max))
      continue;

    if(paraHead->vary==LALINFERENCE_PARAM_CIRCULAR) {
      /* For cyclic boundaries, mod out by range. */
      if (val > max) {
	offset = val - min;
	*(REAL8 *)paraHead->value = min + fmod(offset, delta);
      } else {
	offset = max - val;
	*(REAL8 *)paraHead->value = max - fmod(offset, delta);
      }
    } else if (paraHead->vary==LALINFERENCE_PARAM_LINEAR && paraHead->type==LALINFERENCE_REAL8_t) {
      /* For linear boundaries, reflect about endpoints of range until
	 within range. SKIP NOISE PARAMETERS (ONLY CHECK REAL8) */
      /*
	This reflection is accomplished by also modding out the 'excess' (that is,
	the difference between the value and the upper or lower bound, as appropriate)
	as done above for cyclic parameters.  However, unlike in the cyclic case, to actually
	agree with what the reflection would give, we must also note whether the range divides
	into the excess an even or an odd number of times, and therefore whether we add the
	remainder to the lower bound or subtract it from the upper bound. Which we do also
	depends on whether the value was below the minimum or above the maximum.
       */
      if (val > max) {
	offset = val - max;
	// How many times do we fold?
	n = (UINT8) (offset/delta);
	// Now make offset the remainder
	offset = fmod(offset, delta);
	if (n % 2){
	  // If we fold an odd number of times,
	  // then add offset to min
	  *(REAL8 *)paraHead->value = min + offset;
	} else {
	  // Otherwise, subtract it from max
	  *(REAL8 *)paraHead->value = max - offset;
	}
      } else {
	// We only get here if val <  min
	offset = min - val;
	// How many times do we fold?
	n = (UINT8) (offset/delta);
	// Now make offset the remainder
	offset = fmod(offset, delta);
	if (n % 2){
	  // If we fold an odd number of times,
	  // then subtract offset from max
	  *(REAL8 *)paraHead->value = max - offset;
	} else {
	  // Otherwise, add it to min
	  *(REAL8 *)paraHead->value = min + offset;
	}
      }
    }
  }
  return;
}


void LALInferenceRotateInitialPhase( LALInferenceVariables *parameter){

  if (parameter == NULL)
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");

  LALInferenceVariableItem *paraHead = NULL;
  LALInferenceVariableItem *paraPhi0 = parameter->head;
  REAL8 rotphi0 = 0.;
  UINT4 idx1 = 0, idx2 = 0;

  for(paraHead=parameter->head;paraHead;paraHead=paraHead->next){
    if (paraHead->vary == LALINFERENCE_PARAM_CIRCULAR &&
        !strcmp(paraHead->name, "psi") ){
      /* if psi is outside the -pi/4 -> pi/4 boundary the set to rotate phi0
         by pi (psi will have been rescaled to be between 0 to 2pi as a
        circular parameter).*/
      if (*(REAL8 *)paraHead->value > LAL_TWOPI ||
        *(REAL8 *)paraHead->value < 0. ) rotphi0 = LAL_PI;

      /* Psi has been found; set flag. */
      idx1++;
    }

    /* If phi0 is found, set some flags to indicate so. */
    if (paraHead->vary == LALINFERENCE_PARAM_CIRCULAR &&
        !strcmp(paraHead->name, "phi0") ){
      idx1++;
      idx2++;
    }

    /* Advance paraPhi0 if phi0 hasn't yet been found. */
    if ( idx2 == 0 ) paraPhi0 = paraPhi0->next;

    /* If both variables have been found, stop here. */
    if ( idx1 == 2 ) break;
  }

  if( rotphi0 != 0. ) *(REAL8 *)paraPhi0->value += rotphi0;

  return;
}

/* Return the log Prior of the variables specified for the sky localisation project, ref: https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/SkyLocComparison#priors, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralSkyLocPrior(LALInferenceRunState *runState, LALInferenceVariables *params,  UNUSED LALInferenceModel *model)
{
  REAL8 logPrior=0.0;
  REAL8 val=0.0;
  static int SkyLocPriorWarning = 0;
  (void)runState;
  LALInferenceVariableItem *item=NULL;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=-INFINITY, max=INFINITY;
  REAL8 logmc=0.0,mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;

  if (!SkyLocPriorWarning ) {
    SkyLocPriorWarning  = 1;
    fprintf(stderr, "SkyLocalization priors are being used. (in %s, line %d)\n", __FILE__, __LINE__);
  }
  /* Check boundaries for signal model parameters */
  for(item=params->head;item;item=item->next)
  {
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else if (LALInferenceCheckMinMaxPrior(priorParams, item->name))
    {
      if(item->type==LALINFERENCE_REAL8_t){
        LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
        if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -INFINITY;
      }
    }
    else if (LALInferenceCheckGaussianPrior(priorParams, item->name))
    {
      if(item->type==LALINFERENCE_REAL8_t){
	REAL8 mean,stdev;
	val = *(REAL8 *)item->value;
	LALInferenceGetGaussianPrior(priorParams, item->name, &mean, &stdev);
	logPrior+= -0.5*(mean-val)*(mean-val)/stdev/stdev - 0.5*log(LAL_TWOPI) - log(stdev);
      }
    }
  }

  /*Use a uniform in log D distribution*/

  if(LALInferenceCheckVariable(params, "flow") &&
          LALInferenceCheckVariableNonFixed(params, "flow")) {
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"flow"));
  }

  if(LALInferenceCheckVariable(params,"distance"))
    logPrior-=log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));
  if(LALInferenceCheckVariable(params,"theta_jn"))
    logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_jn"))));
  if(LALInferenceCheckVariable(params,"declination"))
  {
    /* Check that this is not an output variable */
    if(LALInferenceGetVariableVaryType(params,"declination")==LALINFERENCE_PARAM_LINEAR)
      logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
  }
  if(LALInferenceCheckVariable(params,"tilt_spin1"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"tilt_spin1");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin1"))));
    }
  }
  if(LALInferenceCheckVariable(params,"tilt_spin2"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"tilt_spin2");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin2"))));
    }
  }
  /*priors uniform in the individual masses. Not taking into account if mtot_max < m1_max+m2_max */
  if(LALInferenceCheckVariable(params,"eta")||LALInferenceCheckVariable(params,"q")) {
    if(LALInferenceCheckVariable(params,"logmc")) {
      logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
      if(LALInferenceCheckVariable(params,"q")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"q");
        LALInferenceMcQ2Masses(exp(logmc),q,&m1,&m2);
        logPrior+=log(m1*m1);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"eta");
        LALInferenceMcEta2Masses(exp(logmc),eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
      }
      /*careful using LALInferenceMcEta2Masses, it returns m1>=m2*/
    } else if(LALInferenceCheckVariable(params,"chirpmass")) {
      mc=*(REAL8 *)LALInferenceGetVariable(params,"chirpmass");
      if(LALInferenceCheckVariable(params,"q")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"q");
        LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        logPrior+=log(m1*m1/mc);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"eta");
        LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2))/((m1-m2)*pow(eta,3.0/5.0)));
      }
    }
  }

  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      return -INFINITY;

  if(LALInferenceCheckVariable(priorParams,"MTotMin"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
      return -INFINITY;
  /* Check for individual mass priors */
  if(LALInferenceCheckVariable(priorParams,"mass1_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass1_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass2_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
				  return -INFINITY;
  if(LALInferenceCheckVariable(priorParams,"mass2_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
				  return -INFINITY;

  //PSD priors are Gaussian
  if(LALInferenceCheckVariable(params, "psdscale"))
  {
    UINT4 i;
    UINT4 j;

    REAL8 var;
    REAL8 mean = 1.0;
    REAL8 prior= 0.0;
    UINT4 psdGaussianPrior;

    REAL8Vector *sigma = *((REAL8Vector **)LALInferenceGetVariable(priorParams, "psdsigma"));
    gsl_matrix *nparams = *((gsl_matrix **)LALInferenceGetVariable(params,"psdscale"));

    min=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_min");
    max=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_max");

    psdGaussianPrior = *(UINT4 *)LALInferenceGetVariable(priorParams,"psdGaussianPrior");

    for(i=0; i<(UINT4)nparams->size1; i++)
    {
      for(j=0; j<(UINT4)nparams->size2; j++)
      {
        var = sigma->data[j]*sigma->data[j];
        val = gsl_matrix_get(nparams,i,j);
        //reject prior
        if(val < min || val > max) return -INFINITY;
        else if(psdGaussianPrior)prior += -0.5*( (mean-val)*(mean-val)/var + log(2.0*LAL_PI*var) );
      }
    }
    logPrior+=prior;
  }

  /* Calibration parameters */
  //logPrior += LALInferenceSplineCalibrationPrior(runState, params);
  logPrior += LALInferenceConstantCalibrationPrior(runState, params);
  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralSkyLocCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, UNUSED LALInferenceModel *model, double *Cube, void *context)
{
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *timeID = &info[2][0];
    REAL8 min=-INFINITY, max=INFINITY;
    int i = 0;

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(params,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(params,"SKY_FRAME");
    double azimuth, cosalpha, lat, longitude, t0, tc;

    if(SKY_FRAME==1)
    {
      // detector frame azimuth
      if (LALInferenceCheckVariable(params, "azimuth"))
      {
        item = LALInferenceGetItem(params, "azimuth");
        if (item->vary != LALINFERENCE_PARAM_FIXED)
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, "azimuth", (void *)&min, (void *)&max);
          azimuth = LALInferenceCubeToFlatPrior(Cube[i], min, max);
          LALInferenceSetVariable(params, "azimuth", &azimuth);
          i++;
        }
      }

      // detector frame cosalpha
      if (LALInferenceCheckVariable(params, "cosalpha"))
      {
        item = LALInferenceGetItem(params, "cosalpha");
        if (item->vary != LALINFERENCE_PARAM_FIXED)
        {
          LALInferenceGetMinMaxPrior(runState->priorArgs, "cosalpha", (void *)&min, (void *)&max);
          cosalpha = LALInferenceCubeToFlatPrior(Cube[i], min, max);
          LALInferenceSetVariable(params, "cosalpha", &cosalpha);
          i++;
        }
      }
    }
    else
    {
      // latitude
      if (LALInferenceCheckVariable(params, "declination"))
      {
        item = LALInferenceGetItem(params, "declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            lat = asin(2.0 * Cube[i] - 1.0);
            LALInferenceSetVariable(params, "declination", &lat);
            i++;
        }
      }

      // longitude
      if (LALInferenceCheckVariable(params, "rightascension"))
      {
        item = LALInferenceGetItem(params, "rightascension");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "rightascension", (void *)&min, (void *)&max);
            longitude = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &longitude);
            i++;
        }
      }
    }

    // phi
    if (LALInferenceCheckVariable(params, "phase"))
    {
        item = LALInferenceGetItem(params, "phase");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phase", (void *)&min, (void *)&max);
            double phi = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phase", &phi);
            i++;
        }
    }

    // psi
    if (LALInferenceCheckVariable(params, "polarisation"))
    {
        item = LALInferenceGetItem(params, "polarisation");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "polarisation", (void *)&min, (void *)&max);
            double psi = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "polarisation", &psi);
            i++;
        }
    }

    if(SKY_FRAME==1)
    {
      // detector frame t0
      if (LALInferenceCheckVariable(params, "t0"))
      {
        item = LALInferenceGetItem(params, "t0");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "t0", (void *)&min, (void *)&max);
            t0 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "t0", &t0);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }
    else
    {
      // time
      if (LALInferenceCheckVariable(params, "time"))
      {
        item = LALInferenceGetItem(params, "time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }

    // check if mchirp is fixed
    double mc = 0.0, eta = 0.0, q = 0.0, m1 = 0.0, m2 = 0.0, m = 0.0;
    if( LALInferenceCheckVariable(params,"logmc") )
    {
        item = LALInferenceGetItem(params, "logmc");
        if(item->vary == LALINFERENCE_PARAM_FIXED) mc = exp(*(REAL8 *)LALInferenceGetVariable(params, "logmc"));
    }
    else if( LALInferenceCheckVariable(params,"chirpmass") )
    {
        item = LALInferenceGetItem(params, "chirpmass");
        if(item->vary == LALINFERENCE_PARAM_FIXED) mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
    }

    // check if eta is fixed
    if( LALInferenceCheckVariable(params,"eta") )
    {
        item = LALInferenceGetItem(params, "eta");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "eta");
            if( mc != 0.0 ) LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        }
    }
    else if( LALInferenceCheckVariable(params,"q") )
    {
        item = LALInferenceGetItem(params, "q");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "q");
            if( mc != 0.0 ) LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        }
    }

    //m1 & m2
    if( m1 == 0.0 && m2 == 0.0 )
    {
      REAL8 m1_min = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass1_min");
      REAL8 m1_max = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass1_max");
      REAL8 m2_min = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass2_min");
      REAL8 m2_max = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass2_max");
      if( m1_min == m1_max && m2_min == m2_max)
        {
            m1 = m1_min;
            m2 = m2_min;
            m = m1 + m2;
            eta = m1 * m2 / (m*m);
            mc = pow(eta,0.6) * m;
            q = m2 / m1; // asymmetric mass ratio, m1 >= m2
        }
        else
        {
            m1 = LALInferenceCubeToFlatPrior(Cube[i], m1_min, m1_max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i+1], m2_min, m2_max);
            if(m1<m2)
            {
                double temp = m2;
                m2 = m1;
                m1 = temp;
            }

            m = m1 + m2;
            eta = m1 * m2 / (m*m);
            mc = pow(eta,0.6) * m;
            q = m2 / m1; // asymmetric mass ratio, m1 >= m2
            i++;
            i++;
        }

        // chirp mass and eta/q
        if(LALInferenceCheckVariable(params,"eta")||LALInferenceCheckVariable(params,"q"))
        {
            if(LALInferenceCheckVariable(params,"logmc"))
            {
                double logmc = log(mc);
                LALInferenceSetVariable(params, "logmc", &logmc);
            }
            else if(LALInferenceCheckVariable(params,"chirpmass"))
            {
                LALInferenceSetVariable(params, "chirpmass", &mc);
            }

                  if(LALInferenceCheckVariable(params,"q"))
                LALInferenceSetVariable(params, "q", &q);
            else if(LALInferenceCheckVariable(params,"eta"))
                    LALInferenceSetVariable(params, "eta", &eta);
        }
    }

    // distance
    double dist;
    if( LALInferenceCheckVariable(params,"logdistance") )
    {
        item = LALInferenceGetItem(params, "logdistance");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "logdistance", (void *)&min, (void *)&max);
            min = exp(min); max = exp(max);
            dist = LALInferenceCubeToLogFlatPrior(Cube[i], min, max);
            double logdist = log(dist);
            LALInferenceSetVariable(params, "logdistance", &logdist);
            i++;
        }
    }
    else if( LALInferenceCheckVariable(params,"distance") )
    {
        item = LALInferenceGetItem(params, "distance");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "distance", (void *)&min, (void *)&max);
            dist = LALInferenceCubeToLogFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "distance", &dist);
            i++;
        }
    }

    // a_spin1
    if(LALInferenceCheckVariable(params,"a_spin1"))
    {
        item = LALInferenceGetItem(params, "a_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin1", &a_spin1);
            i++;
        }
    }

    // a_spin2
    if(LALInferenceCheckVariable(params,"a_spin2"))
    {
        item = LALInferenceGetItem(params, "a_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin2", &a_spin2);
            i++;
        }
    }

    // theta_JN for system-frame parameters
    if(LALInferenceCheckVariable(params,"costheta_jn"))
    {
        item = LALInferenceGetItem(params, "costheta_jn");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "costheta_jn", (void *)&min, (void *)&max);
            double costheta_JN = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "costheta_jn", &costheta_JN);
            i++;
        }
    }

    // phi_JL for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi_jl"))
    {
        item = LALInferenceGetItem(params, "phi_jl");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_jl", (void *)&min, (void *)&max);
            double phi_JL = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_jl", &phi_JL);
            i++;
        }
    }

    // tilt of spin 1 for system-frame parameters
    if(LALInferenceCheckVariable(params,"tilt_spin1"))
    {
        item = LALInferenceGetItem(params, "tilt_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "tilt_spin1", (void *)&min, (void *)&max);
            double tilt_spin1 = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "tilt_spin1", &tilt_spin1);
            i++;
        }
    }

    // tilt of spin 2 for system-frame parameters
    if(LALInferenceCheckVariable(params,"tilt_spin2"))
    {
        item = LALInferenceGetItem(params, "tilt_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "tilt_spin2", (void *)&min, (void *)&max);
            double tilt_spin2 = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "tilt_spin2", &tilt_spin2);
            i++;
        }
    }

    // phi12 for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi12"))
    {
        item = LALInferenceGetItem(params, "phi12");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi12", (void *)&min, (void *)&max);
            double phi12 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi12", &phi12);
            i++;
        }
    }

    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);
    UINT4 ConstCalib = LALInferenceCubeToConstantCalibrationPrior(runState, params, &i, Cube, context);
    //UINT4 SplineCalib = LALInferenceCubeToSplineCalibrationPrior(runState, params, &i, Cube, context);

    /* Check boundaries */
    if (ScaleTest==0 || ConstCalib==0 /* || SplineCalib==0 */) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else
        {
            LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
            if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return 0 ;
        }
    }

    /* Check for real mass values */
    if(LALInferenceCheckVariable(params,"logmc"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"logmc")))
            return 0;
    if(LALInferenceCheckVariable(params,"chirpmass"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass")))
            return 0;
    if(LALInferenceCheckVariable(params,"eta"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"eta"))
           ||*(REAL8 *)LALInferenceGetVariable(params,"eta") < 0.0
           || *(REAL8 *)LALInferenceGetVariable(params,"eta") > 0.25)
            return 0;
    if(LALInferenceCheckVariable(params,"q"))
        if(isnan(*(REAL8 *)LALInferenceGetVariable(params,"q"))
           ||*(REAL8 *)LALInferenceGetVariable(params,"q") < 0.0
           || *(REAL8 *)LALInferenceGetVariable(params,"q") > 1.0)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMin"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMax"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
            return 0;
  /* Check for individual mass priors */
  if(LALInferenceCheckVariable(priorParams,"mass1_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
				  return 0;
  if(LALInferenceCheckVariable(priorParams,"mass1_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
				  return 0;
  if(LALInferenceCheckVariable(priorParams,"mass2_min"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
				  return 0;
  if(LALInferenceCheckVariable(priorParams,"mass2_max"))
		  if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
				  return 0;

    return 1;
}


typedef struct {
  double M1;
  double McMin;
  double McMax;
  double massRatioMin;
  double massRatioMax;
} innerData;

static double qInnerIntegrand(double M2, void *viData) {
  innerData *iData = (innerData *)viData;
  double Mc = pow(M2*iData->M1, 3.0/5.0)/pow(M2+iData->M1, 1.0/5.0);
  double q = M2/iData->M1;
  if (Mc < iData->McMin || Mc > iData->McMax || q < iData->massRatioMin || q > iData->massRatioMax) {
    return 0.0;
  } else {
    return pow(Mc, -11.0/6.0);
  }
}

#define LALINFERENCE_PRIOR_SQR(x) ((x)*(x))

static double etaInnerIntegrand(double M2, void *viData) {
  innerData *iData = (innerData *)viData;
  double Mc = pow(M2*iData->M1, 3.0/5.0)/pow(M2+iData->M1, 1.0/5.0);
  double eta = M2*iData->M1/LALINFERENCE_PRIOR_SQR(M2+iData->M1);
  if (Mc < iData->McMin || Mc > iData->McMax || eta < iData->massRatioMin || eta > iData->massRatioMax) {
    return 0.0;
  } else {
    return pow(Mc, -11.0/6.0);
  }
}

#undef LALINFERENCE_PRIOR_SQR

typedef struct {
    gsl_integration_workspace *wsInner;
    size_t wsInnerSize;
    double McMin;
    double McMax;
    double massRatioMin;
    double massRatioMax;
    double MTotMax;
    double MMin;
    double epsabs;
    double epsrel;
    gsl_function innerIntegrand;
} outerData;

#define LALINFERENCE_PRIOR_MIN(x, y) ((x) < (y) ? (x) : (y))

static double outerIntegrand(double M1, void *voData) {

    outerData *oData = (outerData *)voData;
    gsl_function f;
    innerData iData;
    double result, err;

    iData.M1 = M1;
    iData.McMin = oData->McMin;
    iData.McMax = oData->McMax;
    iData.massRatioMin = oData->massRatioMin;
    iData.massRatioMax = oData->massRatioMax;

    f.function=(oData->innerIntegrand.function);

    f.params = &iData;

    gsl_integration_qag(&f, oData->MMin, LALINFERENCE_PRIOR_MIN(M1, oData->MTotMax-M1), oData->epsabs, oData->epsrel,
                        oData->wsInnerSize, GSL_INTEG_GAUSS61, oData->wsInner, &result, &err);

    return result;
}

#undef LALINFERENCE_PRIOR_MIN

REAL8 LALInferenceComputePriorMassNorm(const double MMin, const double MMax, const double MTotMax,
                    const double McMin, const double McMax,
                    const double massRatioMin, const double massRatioMax, const char *massRatioName) {
    const double epsabs = 1e-8;
    const double epsrel = 1e-8;
    const size_t wsSize = 10000;
    double result, err;
    outerData oData;
    gsl_function f;

    if(!massRatioName)
        XLAL_ERROR_REAL8(XLAL_EFAULT, "Null arguments received.");
    else if(!strcmp(massRatioName,"q"))
        oData.innerIntegrand.function = &qInnerIntegrand;
    else if(!strcmp(massRatioName,"eta"))
        oData.innerIntegrand.function = &etaInnerIntegrand;
    else
        XLAL_ERROR_REAL8(XLAL_ENAME, "Invalid mass ratio name specified");

    // Disable GSL error reporting in favour of XLAL (the integration routines are liable to fail).
    gsl_error_handler_t *oldHandler = gsl_set_error_handler_off();

    gsl_integration_workspace *wsOuter = gsl_integration_workspace_alloc(wsSize);
    gsl_integration_workspace *wsInner = gsl_integration_workspace_alloc(wsSize);

    oData.wsInnerSize = wsSize;
    oData.wsInner = wsInner;
    oData.McMin = McMin;
    oData.McMax = McMax;
    oData.massRatioMin  = massRatioMin;
    oData.massRatioMax  = massRatioMax;
    oData.MTotMax = MTotMax;
    oData.epsabs = epsabs;
    oData.epsrel = epsrel;
    oData.MMin = MMin;

    f.function = &outerIntegrand;
    f.params = &oData;

    int status = gsl_integration_qag(&f, MMin, MMax, epsabs, epsrel, wsSize, GSL_INTEG_GAUSS61, wsOuter,
                        &result, &err);

    if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

    gsl_set_error_handler(oldHandler);

    gsl_integration_workspace_free(wsOuter);
    gsl_integration_workspace_free(wsInner);

    return result;
}

/* Function to add the min and max values for the prior onto the priorArgs */
void LALInferenceAddMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 *min, REAL8 *max, LALInferenceVariableType type){
  if (*min >= *max)
    XLAL_ERROR_VOID(XLAL_EINVAL, "Minimum must be less than maximum, but %f >= %f.", *min, *max);

  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];

  sprintf(minName,"%s_min",name);
  sprintf(maxName,"%s_max",name);

  LALInferenceAddVariable(priorArgs,minName,min,type,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs,maxName,max,type,LALINFERENCE_PARAM_FIXED);
  return;
}

/* Function to remove the min and max values for the prior onto the priorArgs */
void LALInferenceRemoveMinMaxPrior(LALInferenceVariables *priorArgs, const char *name){
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];

  sprintf(minName,"%s_min",name);
  sprintf(maxName,"%s_max",name);

  LALInferenceRemoveVariable(priorArgs, minName);
  LALInferenceRemoveVariable(priorArgs, maxName);
  return;
}

/* Check for a min/max prior set */
int LALInferenceCheckMinMaxPrior(LALInferenceVariables *priorArgs, const char *name)
{
  char minName[VARNAME_MAX+4];
  char maxName[VARNAME_MAX+4];
  sprintf(minName,"%s_min",name);
  sprintf(maxName,"%s_max",name);

  return (LALInferenceCheckVariable(priorArgs,minName) && LALInferenceCheckVariable(priorArgs,maxName));
}

/* Get the min and max values of the prior from the priorArgs list, given a name */
void LALInferenceGetMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 *min, REAL8 *max)
{
    char minName[VARNAME_MAX+4];
    char maxName[VARNAME_MAX+4];
    void *ptr=NULL;
    sprintf(minName,"%s_min",name);
    sprintf(maxName,"%s_max",name);

    ptr=LALInferenceGetVariable(priorArgs,minName);
    if(ptr) *min=*(REAL8*)ptr;
    else XLAL_ERROR_VOID(XLAL_EFAILED);
    ptr=LALInferenceGetVariable(priorArgs,maxName);
    if(ptr) *max=*(REAL8*)ptr;
    else XLAL_ERROR_VOID(XLAL_EFAILED);
    return;
}


/* Check for a Gaussian Prior of the standard form */
int LALInferenceCheckGaussianPrior(LALInferenceVariables *priorArgs, const char *name)
{
  char meanName[VARNAME_MAX+14];
  char sigmaName[VARNAME_MAX+15];
  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);
  return (LALInferenceCheckVariable(priorArgs,meanName) && LALInferenceCheckVariable(priorArgs,sigmaName));
}

/* Function to add the min and max values for the prior onto the priorArgs */
void LALInferenceAddGaussianPrior( LALInferenceVariables *priorArgs,
                                   const char *name, REAL8 *mu, REAL8 *sigma,
                                   LALInferenceVariableType type ){
  char meanName[VARNAME_MAX+14];
  char sigmaName[VARNAME_MAX+15];

  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);

  LALInferenceAddVariable(priorArgs,meanName,mu,type,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs,sigmaName,sigma,type,LALINFERENCE_PARAM_FIXED);
  return;
}

/* Function to remove the min and max values for the prior onto the priorArgs */
void LALInferenceRemoveGaussianPrior(LALInferenceVariables *priorArgs, const char *name){
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];

  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);

  LALInferenceRemoveVariable(priorArgs, meanName);
  LALInferenceRemoveVariable(priorArgs, sigmaName);
  return;
}

/* Get the min and max values of the prior from the priorArgs list, given a name
*/
void LALInferenceGetGaussianPrior(LALInferenceVariables *priorArgs,
                                  const char *name, REAL8 *mu, REAL8 *sigma)
{
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  void *ptr=NULL;

  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);

  ptr = LALInferenceGetVariable(priorArgs, meanName);
  if ( ptr ) *mu = *(REAL8*)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, sigmaName);
  if ( ptr ) *sigma = *(REAL8*)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  return;
}

/* Function to add a correlation matrix to the prior onto the priorArgs */
void LALInferenceAddCorrelatedPrior(LALInferenceVariables *priorArgs,
                                    const char *name, gsl_matrix **cor,
                                    REAL8 *mu, REAL8 *sigma, UINT4 *idx){
  char corName[VARNAME_MAX];
  char invName[VARNAME_MAX];
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"correlation_matrix");
  sprintf(invName,"inverse_correlation_matrix");

  /* add correlation matrix if not already added */
  if ( !LALInferenceCheckVariable( priorArgs, corName ) ){
    LALInferenceAddVariable(priorArgs, corName, cor, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);

    /* get the matrix we've just added */
    gsl_matrix *thiscor = NULL, *usecor = NULL;
    thiscor = *(gsl_matrix **)LALInferenceGetVariable( priorArgs, corName );
    usecor = gsl_matrix_alloc( thiscor->size1, thiscor->size2 );
    XLAL_CALLGSL( gsl_matrix_memcpy( usecor, thiscor ) );

    gsl_matrix *invcor = gsl_matrix_alloc( usecor->size1, usecor->size2 );
    gsl_permutation *p = gsl_permutation_alloc( usecor->size1 );
    INT4 s;

    /* check correlation matrix is positive definite */
    if( !LALInferenceCheckPositiveDefinite( usecor, usecor->size1 ) ){
      XLAL_ERROR_VOID( XLAL_EFUNC | XLAL_EINVAL, "Error... matrix is not positive definite!"  );
    }

    /* invert correlation matrix */
    XLAL_CALLGSL( gsl_linalg_LU_decomp( usecor, p, &s ) );
    XLAL_CALLGSL( gsl_linalg_LU_invert( usecor, p, invcor ) );
    XLAL_CALLGSL( gsl_permutation_free( p ) );
    XLAL_CALLGSL( gsl_matrix_free( usecor) );

    LALInferenceAddVariable(priorArgs, invName, &invcor, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
  }

  sprintf(meanName, "%s_correlation_mean", name);
  sprintf(sigmaName, "%s_correlation_sigma", name);
  sprintf(idxName,"%s_index",name);

  LALInferenceAddVariable(priorArgs, meanName, mu, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs, sigmaName, sigma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs, idxName, idx, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  return;
}

/* Get the correlation matrix and parameter index */
void LALInferenceGetCorrelatedPrior(LALInferenceVariables *priorArgs,
                                    const char *name, gsl_matrix **cor, gsl_matrix **invcor,
                                    REAL8 *mu, REAL8 *sigma, UINT4 *idx){
  char corName[VARNAME_MAX];
  char invName[VARNAME_MAX];
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  char idxName[VARNAME_MAX];
  void *ptr = NULL;

  sprintf(corName, "correlation_matrix");
  sprintf(invName,"inverse_correlation_matrix");
  sprintf(meanName, "%s_correlation_mean", name);
  sprintf(sigmaName, "%s_correlation_sigma", name);
  sprintf(idxName,"%s_index",name);

  ptr = LALInferenceGetVariable(priorArgs, corName);
  if ( ptr ) *cor = *(gsl_matrix **)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, invName);
  if ( ptr ) *invcor = *(gsl_matrix **)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, meanName);
  if ( ptr ) *mu = *(REAL8 *)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, sigmaName);
  if ( ptr ) *sigma = *(REAL8 *)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, idxName);
  if ( ptr ) *idx = *(UINT4 *)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  return;
}

/* Remove the correlated prior (for all variables) */
void LALInferenceRemoveCorrelatedPrior(LALInferenceVariables *priorArgs){
  char corName[VARNAME_MAX];
  char invName[VARNAME_MAX];
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"correlation_matrix");
  sprintf(invName,"inverse_correlation_matrix");
  sprintf(meanName, "_correlation_mean");
  sprintf(sigmaName, "_correlation_sigma");
  sprintf(idxName,"_index");

  LALInferenceRemoveVariable(priorArgs, corName);
  LALInferenceRemoveVariable(priorArgs, invName);

  /* remove correlated prior parameters */
  LALInferenceVariableItem *item = priorArgs->head;
  for( ; item; item = item->next ){
    if ( strstr(item->name, meanName) != NULL && strstr(item->name, sigmaName) != NULL && strstr(item->name, idxName) != NULL ){
      LALInferenceRemoveVariable(priorArgs, meanName);
      LALInferenceRemoveVariable(priorArgs, sigmaName);
      LALInferenceRemoveVariable(priorArgs, idxName);
    }
  }

  return;
}

/* Check for a correlated prior of the standard form */
int LALInferenceCheckCorrelatedPrior(LALInferenceVariables *priorArgs,
                                     const char *name){
  char corName[VARNAME_MAX];
  char invName[VARNAME_MAX];
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"correlation_matrix");
  sprintf(invName,"inverse_correlation_matrix");
  sprintf(meanName, "%s_correlation_mean", name);
  sprintf(sigmaName, "%s_correlation_sigma", name);
  sprintf(idxName,"%s_index",name);

  return (LALInferenceCheckVariable(priorArgs,corName) &&
          LALInferenceCheckVariable(priorArgs,invName) &&
          LALInferenceCheckVariable(priorArgs,idxName) &&
          LALInferenceCheckVariable(priorArgs,meanName) &&
          LALInferenceCheckVariable(priorArgs,sigmaName));
}


/* Function to add a Gaussian Mixture Model prior (without any hyperparameters)
 * \c name is either a single parameter name (e.g. "H0" for a 1D GMM) or a set of parameter
 * names separated by colons, e.g. "H0:COSIOTA" for a multivariate GMM. The means, covariance
 * matrices, and weights of each mode must be supplied, along with vectors containg the
 * minium and maximum allowed prior extent for each parameter. */
void LALInferenceAddGMMPrior( LALInferenceVariables *priorArgs, const char *name,
                              REAL8Vector ***mus, gsl_matrix ***covs,
                              REAL8Vector **weights, REAL8Vector **minrange, REAL8Vector **maxrange ){
  char gmmParsName[VARNAME_MAX] = "gmm_parameter_lists"; // contains list of ':'-separated lists of GMM parameters
  LALStringVector *parLists = NULL;

  /* check if any other GMM prior already exist */
  if ( !LALInferenceCheckVariable( priorArgs, gmmParsName ) ){
    /* if not then create a new list */
    parLists = XLALAppendString2Vector( parLists, name );
  }
  else{
    /* otherwise add to existing list */
    parLists = *(LALStringVector **)LALInferenceGetVariable( priorArgs, gmmParsName );
    parLists = XLALAppendString2Vector( parLists, name );
  }
  LALInferenceAddVariable( priorArgs, gmmParsName, &parLists, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );

  char musName[VARNAME_MAX];     /* array of means for each parameter for each GMM component */
  char sigmasName[VARNAME_MAX];  /* array of standard deviations from each parameter for each GMM component */
  char weightsName[VARNAME_MAX]; /* weights for each GMM component */
  char corName[VARNAME_MAX];     /* correlation matrices of GMM components */
  char invcorName[VARNAME_MAX];  /* inverse correlation matrices of GMM components */
  char detName[VARNAME_MAX];     /* determinants of the GMM components */
  char minName[VARNAME_MAX];     /* minimum extent of prior for each parameter */
  char maxName[VARNAME_MAX];     /* maximum extent of prior for each parameter */

  sprintf(musName, "%s_gmm_sigmas", name);
  sprintf(sigmasName, "%s_gmm_mus", name);
  sprintf(weightsName, "%s_gmm_weights", name);
  sprintf(corName, "%s_gmm_cors", name);
  sprintf(invcorName, "%s_gmm_invcors", name);
  sprintf(detName, "%s_gmm_dets", name);
  sprintf(minName, "%s_gmm_min", name);
  sprintf(maxName, "%s_gmm_max", name);

  /* get number of parameters */
  TokenList *toks = NULL;
  if ( XLALCreateTokenList( &toks, name, ":" ) != XLAL_SUCCESS ){
    XLAL_ERROR_VOID( XLAL_EINVAL, "Could not create ':' separated token list for GMM prior" );
  }
  UINT4 npars = toks->nTokens;

  if ( !mus[0][0] || !covs[0][0] || !weights[0] ){
    XLAL_ERROR_VOID( XLAL_EINVAL, "GMM means, covariance matrices and weights must all be specified" );
  }

  /* get number of modes from weights vector */
  UINT4 nmodes = weights[0]->length;

  /* convert covariance matrix to correlation matrix (this will avoid numerical dynamic range issues when
   * drawing from multivariate Gaussians) */
  REAL8Vector **sigmas = XLALCalloc(nmodes, sizeof(REAL8Vector *)); /* get standard deviations of parameters for each mode */
  UINT4 j = 0, i = 0;
  gsl_matrix **cormat = XLALCalloc(nmodes, sizeof(gsl_matrix *)); // correlation matrices
  gsl_matrix **invcor = XLALCalloc(nmodes, sizeof(gsl_matrix *)); // inverse correlation matrices
  REAL8Vector *dets = XLALCreateREAL8Vector(nmodes);              // determinants of covariance matrices
  for ( j = 0; j < nmodes; j++ ){
    sigmas[j] = XLALCreateREAL8Vector( npars );
    cormat[j] = gsl_matrix_calloc(npars, npars); /* initialise to all zeros */
    gsl_matrix *invSig = gsl_matrix_calloc(npars, npars); /* diagonal matrix containing inverse of standard deviations */
    gsl_matrix *tmpMat = gsl_matrix_calloc(npars, npars);

    if ( covs[0][j]->size1 != covs[0][j]->size2 || covs[0][j]->size1 != npars ){
      XLAL_ERROR_VOID( XLAL_EINVAL, "GMM covariance matrices must be square and have the correct number of parameters" );
    }

    for ( i = 0; i < npars; i++ ){
      sigmas[j]->data[i] = sqrt(gsl_matrix_get(covs[0][j], i, i));
      gsl_matrix_set(invSig, i, i, 1./sigmas[j]->data[i]);
    }
    /* get correlation coefficient matrix */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invSig, covs[0][j], 0., tmpMat);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmpMat, invSig, 0., cormat[j]);

    gsl_matrix_free(invSig);

    /* get inverse of correlation matrix */
    XLAL_CALLGSL( gsl_matrix_memcpy( tmpMat, cormat[j] ) ); // copy of correlation matrix

    invcor[j] = gsl_matrix_alloc( npars, npars );
    gsl_permutation *p = gsl_permutation_alloc( npars );
    INT4 s;

    /* check correlation matrix is positive definite */
    if( !LALInferenceCheckPositiveDefinite( tmpMat, tmpMat->size1 ) ){
      XLAL_ERROR_VOID( XLAL_EFUNC | XLAL_EINVAL, "Error... matrix is not positive definite!"  );
    }

    /* invert correlation matrix */
    XLAL_CALLGSL( gsl_linalg_LU_decomp( tmpMat, p, &s ) );
    XLAL_CALLGSL( gsl_linalg_LU_invert( tmpMat, p, invcor[j] ) );

    /* get determinant of covariance matrix */
    dets->data[j] = gsl_linalg_LU_det( tmpMat, s );

    XLAL_CALLGSL( gsl_permutation_free( p ) );
    XLAL_CALLGSL( gsl_matrix_free( tmpMat ) );
  }

  /* add values */
  LALInferenceAddVariable( priorArgs, musName, mus, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( priorArgs, sigmasName, &sigmas, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( priorArgs, corName, &cormat, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( priorArgs, invcorName, &invcor, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( priorArgs, detName, &dets, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

  /* make sure weights are normalised to 1 */
  REAL8 weightsum = 0.;
  for ( i = 0; i < weights[0]->length; i++ ){ weightsum += weights[0]->data[i]; }
  for ( i = 0; i < weights[0]->length; i++ ){ weights[0]->data[i] /= weightsum; }

  LALInferenceAddVariable( priorArgs, weightsName, weights, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );

  /* check if ranges are set and if not set to +/- infinity */
  if ( !minrange[0] ){
    minrange[0] = XLALCreateREAL8Vector( npars );
    for ( i = 0; i < npars; i++ ){ minrange[0]->data[i] = -INFINITY; } /* set lower bound to -infinity */
  }
  if ( !maxrange[0] ){
    maxrange[0] = XLALCreateREAL8Vector( npars );
    for ( i = 0; i < npars; i++ ){ maxrange[0]->data[i] = INFINITY; } /* set upper bound to infinity */
  }

  for ( i = 0; i < npars; i++ ){
    if ( minrange[0]->data[i] > maxrange[0]->data[i] ){
      XLAL_ERROR_VOID( XLAL_EINVAL, "Bounds for GMM are wrong" );
    }
  }

  LALInferenceAddVariable( priorArgs, minName, minrange, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
  LALInferenceAddVariable( priorArgs, maxName, maxrange, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED );
  return;
}

/* get the standard deviations, correlation matrices, inverse correlation matrices, means, weights and bounds
 * of the GMM prior, the index of the parameter given in name, and the full (i.e. including all parameters)
 * name of the GMM prior */
void LALInferenceGetGMMPrior( LALInferenceVariables *priorArgs, const char *name,
                              REAL8Vector ***mus, REAL8Vector ***sigmas, gsl_matrix ***cors, gsl_matrix ***invcors,
                              REAL8Vector **weights, REAL8Vector **minrange, REAL8Vector **maxrange,
                              REAL8Vector **dets, UINT4 *idx, CHAR **fullname ){
  /* find list of GMM parameters that this parameter lives in */
  char gmmParsName[VARNAME_MAX] = "gmm_parameter_lists"; // contains list of ':'-separated lists of GMM parameters

  /* check if the given variable is in any of the lists */
  LALStringVector *parLists = *(LALStringVector **)LALInferenceGetVariable( priorArgs, gmmParsName );
  UINT4 numlists = parLists->length;

  UINT4 found = 0, i = 0, j = 0;
  for ( i = 0; i < numlists; i++ ){
    TokenList *toks = NULL;
    if ( XLALCreateTokenList( &toks, parLists->data[i], ":" ) != XLAL_SUCCESS ){ XLAL_ERROR_VOID(XLAL_EFAILED); }
    for ( j = 0; j < toks->nTokens; j++ ){
      if ( !strcmp(toks->tokens[j], name) ){
        found = 1;
        *idx = j; // set index of parameter
        break;
      }
    }
    XLALDestroyTokenList( toks );
    if ( found ){ break; }
  }
  if ( !found ){ XLAL_ERROR_VOID(XLAL_EFAILED); } // parameter could not be found

  *fullname = parLists->data[i];

  char musName[VARNAME_MAX];
  char sigmasName[VARNAME_MAX];
  char weightsName[VARNAME_MAX];
  char corName[VARNAME_MAX];
  char invcorName[VARNAME_MAX];
  char detName[VARNAME_MAX];
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];
  void *ptr = NULL;

  sprintf(musName, "%s_gmm_sigmas", parLists->data[i]);
  sprintf(sigmasName, "%s_gmm_mus", parLists->data[i]);
  sprintf(weightsName, "%s_gmm_weights", parLists->data[i]);
  sprintf(corName, "%s_gmm_cors", parLists->data[i]);
  sprintf(invcorName, "%s_gmm_invcors", parLists->data[i]);
  sprintf(detName, "%s_gmm_dets", parLists->data[i]);
  sprintf(minName, "%s_gmm_min", parLists->data[i]);
  sprintf(maxName, "%s_gmm_max",parLists->data[i]);

  ptr = LALInferenceGetVariable(priorArgs, musName);
  if ( ptr ){ *mus = *(REAL8Vector ***)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, sigmasName);
  if ( ptr ){ *sigmas = *(REAL8Vector ***)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, weightsName);
  if ( ptr ){ *weights = *(REAL8Vector **)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, corName);
  if ( ptr ){ *cors = *(gsl_matrix ***)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, invcorName);
  if ( ptr ){ *invcors = *(gsl_matrix ***)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, detName);
  if ( ptr ){ *dets = *(REAL8Vector **)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, minName);
  if ( ptr ){ *minrange = *(REAL8Vector **)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }

  ptr = LALInferenceGetVariable(priorArgs, maxName);
  if ( ptr ){ *maxrange = *(REAL8Vector **)ptr; }
  else{ XLAL_ERROR_VOID(XLAL_EFAILED); }
  return;
}

/* check for GMM prior */
int LALInferenceCheckGMMPrior(LALInferenceVariables *priorArgs, const char *name){
  char gmmParsName[VARNAME_MAX] = "gmm_parameter_lists"; // contains list of ':'-separated lists of GMM parameters

  /* check if there is a list of GMM parameters */
  if ( !LALInferenceCheckVariable(priorArgs, gmmParsName) ){ return 0; }

  /* check if the given variable is in any of the lists */
  LALStringVector *parLists = *(LALStringVector **)LALInferenceGetVariable( priorArgs, gmmParsName );
  UINT4 numlists = parLists->length;

  if ( numlists == 0 ){ return 0; }

  UINT4 found = 0, i = 0, j = 0;
  for ( i = 0; i < numlists; i++ ){
    TokenList *toks = NULL;
    if ( XLALCreateTokenList( &toks, parLists->data[i], ":" ) != XLAL_SUCCESS ){ return 0; }
    for ( j = 0; j < toks->nTokens; j++ ){
      if ( !strcmp(toks->tokens[j], name) ){
        found = 1;
        break;
      }
    }
    XLALDestroyTokenList( toks );
    if ( found ){ break; }
  }

  if ( !found ){ return 0; } // parameter could not be found

  /* check for required items for the GMM containing the found parameter */
  char musName[VARNAME_MAX];
  char sigmasName[VARNAME_MAX];
  char weightsName[VARNAME_MAX];
  char corName[VARNAME_MAX];
  char invcorName[VARNAME_MAX];
  char detName[VARNAME_MAX];
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];

  sprintf(musName, "%s_gmm_sigmas", parLists->data[i]);
  sprintf(sigmasName, "%s_gmm_mus", parLists->data[i]);
  sprintf(weightsName, "%s_gmm_weights", parLists->data[i]);
  sprintf(corName, "%s_gmm_cors", parLists->data[i]);
  sprintf(invcorName, "%s_gmm_invcors", parLists->data[i]);
  sprintf(detName, "%s_gmm_dets", parLists->data[i]);
  sprintf(minName, "%s_gmm_min", parLists->data[i]);
  sprintf(maxName, "%s_gmm_max", parLists->data[i]);

  return (LALInferenceCheckVariable(priorArgs, musName) &&
          LALInferenceCheckVariable(priorArgs, sigmasName) &&
          LALInferenceCheckVariable(priorArgs, weightsName) &&
          LALInferenceCheckVariable(priorArgs, corName) &&
          LALInferenceCheckVariable(priorArgs, invcorName) &&
          LALInferenceCheckVariable(priorArgs, minName) &&
          LALInferenceCheckVariable(priorArgs, maxName) &&
          LALInferenceCheckVariable(priorArgs, detName));
}


/* remove GMM Prior (either using the full name [i.e. comma separated list of GMM parameters], or
 * a single parameter names, from the GMM prior. Note that the corresponding entry in the "gmm_parameter_lists"
 * will be set to a string termination character, or if it just contains one entry it will be removed. */
void LALInferenceRemoveGMMPrior( LALInferenceVariables *priorArgs, const char *name ){
  char gmmParsName[VARNAME_MAX] = "gmm_parameter_lists"; // contains list of ':'-separated lists of GMM parameters

  /* check if the given variable is in any of the lists */
  LALStringVector *parLists = *(LALStringVector **)LALInferenceGetVariable( priorArgs, gmmParsName );
  UINT4 numlists = parLists->length;

  UINT4 found = 0, i = 0, j = 0;
  for ( i = 0; i < numlists; i++ ){
    if ( !strcmp(name, parLists->data[i]) ){
      found = 1;
      break;
    }
    TokenList *toks = NULL;
    if ( XLALCreateTokenList( &toks, parLists->data[i], ":" ) != XLAL_SUCCESS ){ XLAL_ERROR_VOID(XLAL_EFAILED); }
    for ( j = 0; j < toks->nTokens; j++ ){
      if ( !strcmp(toks->tokens[j], name) ){
        found = 1;
        break;
      }
    }
    XLALDestroyTokenList( toks );
    if ( found ){ break; }
  }

  if ( !found ){ return; } // parameter doesn't seem to be present anyway!

  char musName[VARNAME_MAX];
  char sigmasName[VARNAME_MAX];
  char weightsName[VARNAME_MAX];
  char corName[VARNAME_MAX];
  char invcorName[VARNAME_MAX];
  char detName[VARNAME_MAX];
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];

  sprintf(musName, "%s_gmm_sigmas", parLists->data[i]);
  sprintf(sigmasName, "%s_gmm_mus", parLists->data[i]);
  sprintf(weightsName, "%s_gmm_weights", parLists->data[i]);
  sprintf(corName, "%s_gmm_cors", parLists->data[i]);
  sprintf(invcorName, "%s_gmm_invcors", parLists->data[i]);
  sprintf(detName, "%s_gmm_dets", parLists->data[i]);
  sprintf(minName, "%s_gmm_min", parLists->data[i]);
  sprintf(maxName, "%s_gmm_max", parLists->data[i]);

  LALInferenceRemoveVariable(priorArgs, musName);
  LALInferenceRemoveVariable(priorArgs, sigmasName);
  LALInferenceRemoveVariable(priorArgs, weightsName);
  LALInferenceRemoveVariable(priorArgs, corName);
  LALInferenceRemoveVariable(priorArgs, invcorName);
  LALInferenceRemoveVariable(priorArgs, detName);
  LALInferenceRemoveVariable(priorArgs, minName);
  LALInferenceRemoveVariable(priorArgs, maxName);

  if ( numlists == 1 ){ /* just remove list */
    LALInferenceRemoveVariable(priorArgs, gmmParsName);
  }
  else{ /* set found entry in "gmm_parameter_lists" to string termination character */
    parLists->data[i][0] = '\0';
  }

  return;
}


/* Check for a Fermi-Dirac Prior */
int LALInferenceCheckFermiDiracPrior(LALInferenceVariables *priorArgs, const char *name)
{
  char rName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  sprintf(rName,"%s_fermi_r",name);
  sprintf(sigmaName,"%s_fermi_sigma",name);
  return (LALInferenceCheckVariable(priorArgs,rName) && LALInferenceCheckVariable(priorArgs,sigmaName));
}

/* Function to add the r and sigma values for the prior onto the priorArgs */
void LALInferenceAddFermiDiracPrior( LALInferenceVariables *priorArgs,
                                   const char *name, REAL8 *sigma, REAL8 *r,
                                   LALInferenceVariableType type ){
  char rName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];

  sprintf(rName,"%s_fermi_r",name);
  sprintf(sigmaName,"%s_fermi_sigma",name);

  LALInferenceAddVariable(priorArgs,rName,r,type,LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs,sigmaName,sigma,type,LALINFERENCE_PARAM_FIXED);
  return;
}

/* Function to remove the r and sigma values for the prior onto the priorArgs */
void LALInferenceRemoveFermiDiracPrior(LALInferenceVariables *priorArgs, const char *name){
  char rName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];

  sprintf(rName,"%s_fermi_r",name);
  sprintf(sigmaName,"%s_fermi_sigma",name);

  LALInferenceRemoveVariable(priorArgs, rName);
  LALInferenceRemoveVariable(priorArgs, sigmaName);
  return;
}

/* Get the r and sigma values of the prior from the priorArgs list, given a name */
void LALInferenceGetFermiDiracPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *sigma, REAL8 *r)
{
  char rName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  void *ptr=NULL;

  sprintf(rName,"%s_fermi_r",name);
  sprintf(sigmaName,"%s_fermi_sigma",name);

  ptr = LALInferenceGetVariable(priorArgs, rName);
  if ( ptr ) *r = *(REAL8*)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, sigmaName);
  if ( ptr ) *sigma = *(REAL8*)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  return;
}


/* Check for a prior uniform in the log */
int LALInferenceCheckLogUniformPrior(LALInferenceVariables *priorArgs,
                                     const char *name)
{
  char xminName[VARNAME_MAX];
  char xmaxName[VARNAME_MAX];

  sprintf(xminName, "%s_loguniform_xmin", name);
  sprintf(xmaxName, "%s_loguniform_xmax", name);

  return (LALInferenceCheckVariable(priorArgs, xminName) &&
          LALInferenceCheckVariable(priorArgs, xmaxName));
}

/* Add domain boundary values for the prior onto the priorArgs */
void LALInferenceAddLogUniformPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *xmin, REAL8 *xmax,
                                    LALInferenceVariableType type )
{
  if (*xmin >= *xmax || *xmin < 0. ){
    XLAL_ERROR_VOID(XLAL_EINVAL, "Minimum must be less than maximum (and minumum must be greater than zero), but %f >= %f.", *xmin, *xmax);
  }

  char xminName[VARNAME_MAX];
  char xmaxName[VARNAME_MAX];

  sprintf(xminName, "%s_loguniform_xmin", name);
  sprintf(xmaxName, "%s_loguniform_xmax", name);

  LALInferenceAddVariable(priorArgs, xminName, xmin, type,
                          LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs, xmaxName, xmax, type,
                          LALINFERENCE_PARAM_FIXED);
  return;
}

/* Remove the domain boundary values for the prior onto the priorArgs */
void LALInferenceRemoveLogUniformPrior(LALInferenceVariables *priorArgs,
                                       const char *name)
{
  char xminName[VARNAME_MAX];
  char xmaxName[VARNAME_MAX];

  sprintf(xminName, "%s_loguniform_xmin", name);
  sprintf(xmaxName, "%s_loguniform_xmax", name);

  LALInferenceRemoveVariable(priorArgs, xminName);
  LALInferenceRemoveVariable(priorArgs, xmaxName);
  return;
}

/* Get domain boundary values of the prior from priorArgs list, given a name */
void LALInferenceGetLogUniformPrior(LALInferenceVariables *priorArgs,
                                    const char *name, REAL8 *xmin, REAL8 *xmax)
{
  char xminName[VARNAME_MAX];
  char xmaxName[VARNAME_MAX];
  void *ptr=NULL;

  sprintf(xminName, "%s_loguniform_xmin", name);
  sprintf(xmaxName, "%s_loguniform_xmax", name);

  ptr = LALInferenceGetVariable(priorArgs, xminName);
  if ( ptr ) {
    *xmin = *(REAL8*)ptr;
  } else {
    XLAL_ERROR_VOID(XLAL_EFAILED);
  }

  ptr = LALInferenceGetVariable(priorArgs, xmaxName);
  if ( ptr ) {
    *xmax = *(REAL8*)ptr;
  } else {
    XLAL_ERROR_VOID(XLAL_EFAILED);
  }

  return;
}

void LALInferenceDrawFromPrior( LALInferenceVariables *output,
                                LALInferenceVariables *priorArgs,
                                gsl_rng *rdm) {
  if (output == NULL || priorArgs == NULL || rdm == NULL)
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");

  LALInferenceVariableItem *item = output->head;

  /* check if using a k-D tree as the prior */
  if( LALInferenceCheckVariable( priorArgs, "kDTreePrior" ) ){
    LALInferenceKDTree *tree =
      *(LALInferenceKDTree **)LALInferenceGetVariable(priorArgs, "kDTreePrior");

    /* get parameter template */
    LALInferenceVariables *templt =
      *(LALInferenceVariables **)LALInferenceGetVariable(priorArgs,
                                                         "kDTreePriorTemplate");

    UINT4 Ncell = 8; /* number of points in a prior cell - i.e. controls
                        how fine or coarse the prior looks (default to 8) */

    if( LALInferenceCheckVariable( priorArgs, "kDTreePriorNcell" ) )
      Ncell = *(UINT4 *)LALInferenceGetVariable( priorArgs,"kDTreePriorNcell");

    /* draw all points from the prior distribution */
    REAL8 *proposedPt = XLALCalloc(tree->dim, sizeof(REAL8));

    /* A randomly-chosen point from those in the tree. */
    LALInferenceKDDrawEigenFrame(rdm, tree, proposedPt, Ncell);
    LALInferenceKDREAL8ToVariables(output, proposedPt, templt);
  }
  else{
    for(;item;item=item->next){
      if(item->vary==LALINFERENCE_PARAM_CIRCULAR ||
         item->vary==LALINFERENCE_PARAM_LINEAR)
        LALInferenceDrawNameFromPrior( output, priorArgs, item->name,
                                       item->type, rdm );
    }

    /* remove multivariate deviates value if set */
    if ( LALInferenceCheckVariable( priorArgs, "multivariate_deviates" ) )
      LALInferenceRemoveVariable( priorArgs, "multivariate_deviates" );

  }
}

void LALInferenceDrawNameFromPrior( LALInferenceVariables *output,
                                    LALInferenceVariables *priorArgs,
                                    char *name, LALInferenceVariableType type,
                                    gsl_rng *rdm) {
  if (output == NULL || priorArgs == NULL || name == NULL || rdm == NULL)
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");

  REAL8 tmp = 0.;

  /* test for a Gaussian prior */
  if( LALInferenceCheckGaussianPrior( priorArgs, name ) ){
    REAL8 mu = 0., sigma = 0.;

    LALInferenceGetGaussianPrior( priorArgs, name, (void *)&mu, (void *)&sigma );
    tmp = mu + gsl_ran_gaussian(rdm, (double)sigma);
  }
  /* test for uniform prior */
  else if( LALInferenceCheckMinMaxPrior( priorArgs, name ) ){
    REAL8 min = 0., max = 0.;

    LALInferenceGetMinMaxPrior(priorArgs, name, &min, &max);
    tmp = min + (max-min)*gsl_rng_uniform( rdm );
  }
  /* test for a Fermi-Dirac prior */
  else if( LALInferenceCheckFermiDiracPrior( priorArgs, name ) ){
    REAL8 r = 0., sigma = 0., cp;

    LALInferenceGetFermiDiracPrior(priorArgs, name, &sigma, &r);

    /* use the inverse sampling transform to draw a new sample */
    do { /* numerical issues mean that the analytic solution to this equation can go negative, so make sure that is not the case */
      cp = gsl_rng_uniform( rdm ); /* draw a point uniformly between 0 and 1 */
      tmp = log(-exp(-r) + pow(1. + exp(r), -cp) + exp(-r)*pow(1. + exp(r), -cp));
      tmp *= -sigma;
    } while ( tmp < 0. );
  }
  /* test for a prior uniform in the log */
  else if( LALInferenceCheckLogUniformPrior( priorArgs, name ) ){
    REAL8 xmin = 0., xmax = 0., cp;

    LALInferenceGetLogUniformPrior(priorArgs, name, &xmin, &xmax);

    if( xmin <= 0 ) {
      XLAL_ERROR_VOID(XLAL_EDOM, "Log-uniform min value not positive.");
    } else if( xmax < xmin ) {
      XLAL_ERROR_VOID(XLAL_EDOM, "Log-uniform min value greater than max.");
    }

    /* use the inverse sampling transform to draw a new sample */
    cp = gsl_rng_uniform( rdm ); /* draw a point uniformly between 0 and 1 */

    /* the percentage-point function (PPF, or inverse CDF) for PDF~1/x is:
    \f$x = (\frac{x_{\rm max}}{x_{\rm min}})^{\rm cdf} x_{\rm min}\f$ */
    tmp = xmin * pow(xmax/xmin, cp);
  }
  /* test for a prior drawn from correlated values */
  else if( LALInferenceCheckCorrelatedPrior( priorArgs, name ) ){
    gsl_matrix *cor = NULL, *invcor = NULL;
    REAL8 mu = 0, sigma = 0.;
    UINT4 idx = 0, dims = 0;
    REAL4Vector *tmps = NULL;

    LALInferenceGetCorrelatedPrior( priorArgs, name, &cor, &invcor, &mu, &sigma, &idx );
    dims = cor->size1;

    /* to avoid unnecessary repetition the multivariate deviates are be
       added as a new variable that can be extracted during multiple calls.
       This will be later removed by the LALInferenceDrawFromPrior function. */
    if ( LALInferenceCheckVariable( priorArgs, "multivariate_deviates" ) ){
      tmps = *(REAL4Vector **)LALInferenceGetVariable(priorArgs,
                                                      "multivariate_deviates");
    }
    else{
      RandomParams *randParam;
      UINT4 randomseed = gsl_rng_get(rdm);

      /* check matrix for positive definiteness */
      if( !LALInferenceCheckPositiveDefinite( cor, dims ) ){
        XLAL_ERROR_VOID(XLAL_EFUNC | XLAL_EINVAL, "Matrix is not positive-definite!");
      }

      /* draw values from the multivariate Gaussian described by the correlation matrix */
      tmps = XLALCreateREAL4Vector( dims );
      randParam = XLALCreateRandomParams( randomseed );
      XLALMultiNormalDeviates( tmps, cor, dims, randParam );

      LALInferenceAddVariable( priorArgs, "multivariate_deviates", &tmps,
                               LALINFERENCE_REAL8Vector_t,
                               LALINFERENCE_PARAM_FIXED );
      XLALDestroyRandomParams( randParam );
    }

    /* set random number for given parameter index (converted to a draw from the covariance matrix) */
    tmp = mu + sigma*tmps->data[idx];

    /* free tmps */
    if ( !LALInferenceCheckVariable( priorArgs, "multivariate_deviates" ) )
      XLALDestroyREAL4Vector( tmps );
  }
  /* test if a prior drawn from a Gaussian Mixture Model */
  else if ( LALInferenceCheckGMMPrior( priorArgs, name ) ){
    gsl_matrix **cor, UNUSED **invcor;
    REAL8Vector **mus, **sigmas, *weights = NULL, *minrange = NULL, *maxrange = NULL, UNUSED *dets = NULL;
    REAL8 cp = 0., cumweights = 0.;
    REAL4Vector *tmps = NULL;
    UINT4 idx = 0, dims = 0;
    CHAR *fullname = NULL;

    LALInferenceGetGMMPrior( priorArgs, name, &mus, &sigmas, &cor, &invcor, &weights, &minrange, &maxrange, &dets, &idx, &fullname );

    dims = cor[0]->size1;

    /* set mode to use and multivariate deviates for it  */
    CHAR gmmmvd[VARNAME_MAX], gmmmode[VARNAME_MAX];
    UINT4 thismode = 0;
    sprintf(gmmmvd, "%s_gmm_mvd", fullname);
    sprintf(gmmmode, "%s_gmm_mode", fullname);
    if ( LALInferenceCheckVariable( priorArgs, gmmmvd ) ){ /* get values if already set */
      tmps = *(REAL4Vector **)LALInferenceGetVariable(priorArgs, gmmmvd);
      thismode = LALInferenceGetUINT4Variable( priorArgs, gmmmode);
    }
    else{ /* otherwise create array */
      tmps = XLALCreateREAL4Vector( dims );

      /* get GMM mode to use */
      cp = gsl_rng_uniform( rdm );
      cumweights = weights->data[0];
      /* get index of mode */
      while( cp > cumweights ){
        thismode++;
        cumweights += weights->data[thismode];
      }

      RandomParams *randParam;
      UINT4 randomseed = gsl_rng_get(rdm);
      randParam = XLALCreateRandomParams( randomseed );
      XLALMultiNormalDeviates( tmps, cor[thismode], dims, randParam );
      XLALDestroyRandomParams( randParam );

      /* add to priorArgs */
      LALInferenceAddVariable( priorArgs, gmmmvd, &tmps, LALINFERENCE_void_ptr_t, LALINFERENCE_PARAM_FIXED );

      /* add the mode to the prior args */
      LALInferenceAddVariable( priorArgs, gmmmode, &thismode, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );
    }

    /* set random number for given parameter index (converted to a draw from the covariance matrix) */
    tmp = mus[thismode]->data[idx] + sigmas[thismode]->data[idx]*tmps->data[idx];

    /* get number of parameter from GMM already passed */
    CHAR nparsdonename[VARNAME_MAX];
    UINT4 npars = 1;
    sprintf(nparsdonename, "%s_gmm_npars", fullname);
    if ( LALInferenceCheckVariable( priorArgs, nparsdonename ) ){
      UINT4 nparsdone = LALInferenceGetUINT4Variable( priorArgs, nparsdonename );
      npars = nparsdone;
      npars++;
      LALInferenceRemoveVariable( priorArgs, nparsdonename );
    }

    LALInferenceAddVariable( priorArgs, nparsdonename, &npars, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED );

    /* remove GMM multivariate deviates if all parameter in the GMM have been set */
    if ( npars == dims ){
      LALInferenceRemoveVariable( priorArgs, gmmmvd );
      LALInferenceRemoveVariable( priorArgs, gmmmode );
      LALInferenceRemoveVariable( priorArgs, nparsdonename );
    }
  }
  /* not a recognised prior type */
  else{
    return;
  }

  switch ( type ){
    case LALINFERENCE_REAL4_t:
    {
      REAL4 val = (REAL4)tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }
    case LALINFERENCE_REAL8_t:
    {
      REAL8 val = tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }
    case LALINFERENCE_INT4_t:
    {
      INT4 val = (INT4)tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }
    case LALINFERENCE_INT8_t:
    {
      INT8 val = (INT8)tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }
    //LALInferenceDrawFromPrior() does not handle gsl matrices properly at the moment
    /*case LALINFERENCE_gslMatrix_t:
    {
      REAL8 val = tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }*/
    default:
      XLALPrintWarning("Trying to randomise a non-numeric parameter %s!\n",name);
      break;
  }
}


/* Switch reads true if parameters lie within Malmquist prior */
UINT4 within_malmquist(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model) {
    UINT4 i=0, nifo=0;

    LALInferenceIFOData *ifo = runState->data;
    while (ifo != NULL) {
        nifo++;
        ifo = ifo->next;
    }

    LALInferenceNetworkSNR(params, runState->data, model);
    REAL8 loudest_snr=0.0, second_loudest_snr=0.0;
    for (i=0; i<nifo; i++) {
        if (model->ifo_SNRs[i] > second_loudest_snr) {
            if (model->ifo_SNRs[i] > loudest_snr) {
                second_loudest_snr = loudest_snr;
                loudest_snr = model->ifo_SNRs[i];
            } else {
                second_loudest_snr = model->ifo_SNRs[i];
            }
        }
    }

    REAL8 malmquist_loudest = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_loudest_snr"));
    REAL8 malmquist_second_loudest = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_second_loudest_snr"));
    REAL8 malmquist_network = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_network_snr"));

    if (loudest_snr < malmquist_loudest
          || second_loudest_snr < malmquist_second_loudest
          || model->SNR < malmquist_network)
        return(0);
    else
        return(1);
}


REAL8 LALInferenceAnalyticNullPrior(LALInferenceRunState UNUSED *runState, LALInferenceVariables *params, LALInferenceModel UNUSED *model) {
  REAL8 logPrior=0.0;
  REAL8 logmc=0.0,mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;

  if(LALInferenceCheckVariable(params,"eta")||LALInferenceCheckVariable(params,"q")) {
    if(LALInferenceCheckVariable(params,"logmc")) {
      logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
      if(LALInferenceCheckVariable(params,"q")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"q");
        LALInferenceMcQ2Masses(exp(logmc),q,&m1,&m2);
        logPrior+=log(m1*m1);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"eta");
        LALInferenceMcEta2Masses(exp(logmc),eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
      }
      /*careful using LALInferenceMcEta2Masses, it returns m1>=m2*/
    } else if(LALInferenceCheckVariable(params,"chirpmass")) {
      mc=*(REAL8 *)LALInferenceGetVariable(params,"chirpmass");
      if(LALInferenceCheckVariable(params,"q")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"q");
        LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        logPrior+=log(m1*m1/mc);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"eta");
        LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2))/((m1-m2)*pow(eta,3.0/5.0)));
      }
    }
  }
  logPrior+=LALInferenceFlatBoundedPrior(runState, params);
  return(logPrior);
}

UINT4 LALInferenceAnalyticCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, UNUSED LALInferenceModel *model, double *Cube, void *context) {
    int i = 0;
    REAL8 min=-INFINITY,max=INFINITY;
    REAL8 m1=1.,m2=1.;
    LALInferenceVariableItem *item;

    char **info = (char **)context;
    char *timeID = &info[2][0];

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(params,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(params,"SKY_FRAME");

    if(LALInferenceCheckVariable(params,"mass1"))
    {
        if(LALInferenceGetItem(params,"mass1")->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "mass1", (void *)&min, (void *)&max);
            m1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "mass1", &m1);
        }
        else
        {
            m1=(*(REAL8 *)LALInferenceGetVariable(params,"mass1"));
        }
    }
    if(LALInferenceCheckVariable(params,"mass2"))
    {
        if(LALInferenceGetItem(params,"mass2")->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "mass2", (void *)&min, (void *)&max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "mass2", &m2);
        }
        else
        {
            m2=(*(REAL8 *)LALInferenceGetVariable(params,"mass2"));
        }
    }

    if(LALInferenceCheckVariable(params,"phase"))
    {
        item = LALInferenceGetItem(params,"phase");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phase", (void *)&min, (void *)&max);
            double phase = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phase", &phase);
            i++;
        }
    }

    if(LALInferenceCheckVariable(params,"polarisation"))
    {
        item = LALInferenceGetItem(params,"polarisation");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "polarisation", (void *)&min, (void *)&max);
            double polarisation = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "polarisation", &polarisation);
            i++;
        }
    }

    if(SKY_FRAME==1)
    {
      if(LALInferenceCheckVariable(params,"azimuth"))
      {
        item = LALInferenceGetItem(params,"azimuth");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "azimuth", (void *)&min, (void *)&max);
            double azimuth = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "azimuth", &azimuth);
            i++;
        }
      }
    }
    else
    {
      if(LALInferenceCheckVariable(params,"rightascension"))
      {
        item = LALInferenceGetItem(params,"rightascension");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "rightascension", (void *)&min, (void *)&max);
            double rightascension = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &rightascension);
            i++;
        }
      }
    }

    if(SKY_FRAME==1)
    {
      if(LALInferenceCheckVariable(params,"cosalpha"))
      {
        item = LALInferenceGetItem(params,"cosalpha");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "cosalpha", (void *)&min, (void *)&max);
            double cosalpha = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "cosalpha", &cosalpha);
            i++;
        }
      }
    }
    else
    {
      if(LALInferenceCheckVariable(params,"declination"))
      {
        item = LALInferenceGetItem(params,"declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "declination", (void *)&min, (void *)&max);
            double declination = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "declination", &declination);
            i++;
        }
      }
    }

    if(LALInferenceCheckVariable(params,"distance"))
    {
        item = LALInferenceGetItem(params,"distance");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "distance", (void *)&min, (void *)&max);
            double distance = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "distance", &distance);
            i++;
        }
    }

    if(SKY_FRAME==1)
    {
      if(LALInferenceCheckVariable(params,"t0"))
      {
        item = LALInferenceGetItem(params,"t0");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "t0", (void *)&min, (void *)&max);
            double t0 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "t0", &t0);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }
    else
    {
      if(LALInferenceCheckVariable(params,"time"))
      {
        item = LALInferenceGetItem(params,"time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            double tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            sprintf(timeID,"%d",i);
            i++;
        }
      }
    }

    if(LALInferenceCheckVariable(params,"a_spin1"))
    {
        item = LALInferenceGetItem(params,"a_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin1", &a_spin1);
            i++;
        }
    }

    if(LALInferenceCheckVariable(params,"a_spin2"))
    {
        item = LALInferenceGetItem(params,"a_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "a_spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "a_spin2", &a_spin2);
            i++;
        }
    }

    LALInferenceVariables *priorParams=runState->priorArgs;
    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);
    UINT4 ConstCalib = LALInferenceCubeToConstantCalibrationPrior(runState, params, &i, Cube, context);
    //UINT4 SplineCalib = LALInferenceCubeToSplineCalibrationPrior(runState, params, &i, Cube, context);

    if (ScaleTest==0 || ConstCalib==0 /* || SplineCalib==0 */) return 0;
    else return 1;
}

REAL8 LALInferenceFlatBoundedPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
  LALInferenceVariableItem *cur;
  REAL8 min,max;
  if(!params||!runState) XLAL_ERROR(XLAL_EFAULT, "Encountered NULL pointer in prior");
  if(!runState->priorArgs) return 0.0; /* no prior ranges specified */
  for(cur=params->head;cur;cur=cur->next)
  {
    if(cur->type!=LALINFERENCE_REAL8_t) continue;
    if(LALInferenceCheckMinMaxPrior(runState->priorArgs, cur->name))
    {
      LALInferenceGetMinMaxPrior(runState->priorArgs, cur->name, &min, &max);
      if (min>*(REAL8 *)cur->value || max<*(REAL8 *)cur->value) return -INFINITY;
    }
  }
  return 0.0;
}

REAL8 LALInferenceNullPrior(LALInferenceRunState UNUSED *runState, LALInferenceVariables UNUSED *params, LALInferenceModel UNUSED *model) {
  return 0.0;
}

UINT4 LALInferenceCubeToPSDScaleParams(LALInferenceVariables *priorParams, LALInferenceVariables *params, INT4 *idx, double *Cube, void UNUSED *context)
{
  //char **info = (char **)context;
  //char *header = &info[1][0];
  //char tempstr[50];

  //PSD priors are Gaussian
  if(LALInferenceCheckVariable(params, "psdscale"))
  {
    UINT4 i;
    UINT4 j;

    REAL8 val,min,max;
    REAL8 mean = 1.0;
    UINT4 psdGaussianPrior;

    REAL8Vector *sigma = *((REAL8Vector **)LALInferenceGetVariable(priorParams, "psdsigma"));
    gsl_matrix *nparams = *((gsl_matrix **)LALInferenceGetVariable(params,"psdscale"));

    min=0.1;//*(REAL8 *)LALInferenceGetVariable(priorParams,"psdscale_min");
    max=10.0;//*(REAL8 *)LALInferenceGetVariable(priorParams,"psdscale_max");

    psdGaussianPrior = *(UINT4 *)LALInferenceGetVariable(priorParams,"psdGaussianPrior");

    for(i=0; i<(UINT4)nparams->size1; i++)
    {
      for(j=0; j<(UINT4)nparams->size2; j++)
      {
        // calculate value
        if (psdGaussianPrior)
            val = LALInferenceCubeToGaussianPrior(Cube[(*idx)],mean,sigma->data[j]);
        else
            val = LALInferenceCubeToFlatPrior(Cube[(*idx)],min,max);

        // set value
        gsl_matrix_set(nparams,i,j,val);
        (*idx)++;
      }
    }

    for(i=0; i<(UINT4)nparams->size1; i++)
    {
      for(j=0; j<(UINT4)nparams->size2; j++)
      {
        //reject prior
        val = gsl_matrix_get(nparams,i,j);
        if(val < min || val > max) return 0;
      }
    }
  }

  return 1;
}

/* A simple SineGaussianPrior -- will also work for other simple burst templates (gaussians) */
REAL8 LALInferenceSineGaussianPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model)
{
  REAL8 logPrior=0.0;
  (void)runState;
  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min, max;
  (void) model;
  /* Check boundaries */
  for(;item;item=item->next)
  {
    // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else
    {
      LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
      if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -INFINITY;
    }
  }
  /*Use a distribution uniform in space volume */
  if(LALInferenceCheckVariable(params,"loghrss"))
    logPrior+=-3.0* *(REAL8 *)LALInferenceGetVariable(params,"loghrss");
  else if(LALInferenceCheckVariable(params,"hrss"))
    logPrior+=-4.0* log(*(REAL8 *)LALInferenceGetVariable(params,"hrss"));
  if(LALInferenceCheckVariable(params,"declination"))
    logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
  logPrior += LALInferenceConstantCalibrationPrior(runState, params);
  return(logPrior);
}

/**
 * Prior that converts from a Cube parameter in [0,1] to the flat prior bounded by x1
 * and x2.
 */
REAL8 LALInferenceCubeToFlatPrior(double r, double x1, double x2)
{
    return x1 + r * ( x2 - x1 );
}

/**
 * Prior that converts from a Cube parameter in [0,1] to the flat in log prior bounded
 * by x1 and x2.
 */
REAL8 LALInferenceCubeToLogFlatPrior(double r, double x1, double x2)
{
    double lx1, lx2;
    lx1 = log( x1 );
    lx2 = log( x2 );
    return exp( lx1 + r * ( lx2 - lx1 ) );
}

/**
 * Prior that converts from a Cube parameter in [0,1] to the power prior bounded by x1
 * and x2 with power p.
 */
REAL8 LALInferenceCubeToPowerPrior(double p, double r, double x1, double x2)
{
    double pp = p + 1.0;
    return pow(r * pow(x2, pp) + (1.0 - r) * pow(x1, pp), 1.0 / pp);
}

/**
 * Prior that converts from a Cube parameter in [0,1] to the Gaussian prior with given
 * mean and standard deviation.
 */
REAL8 LALInferenceCubeToGaussianPrior(double r, double mean, double sigma)
{
    return gsl_cdf_gaussian_Pinv(r,sigma) + mean;
}

/**
 * Prior that converts from a Cube parameter in [0,1] to the sine prior with given
 * min (x1) and max (x2) values
 */
REAL8 LALInferenceCubeToSinPrior(double r, double x1, double x2)
{
    return acos((1.0-r)*cos(x1)+r*cos(x2));
}

/* Functions for Fermi-Dirac prior distribution */

/** \brief Return the Fermi-Dirac distribution log prior
 *
 * The function returns the log of the prior for a Fermi-Dirac distribution
 * \f[p(h|\sigma, r, I) = \frac{1}{\sigma\log{\left(1+e^{r} \right)}}\left(e^{((h/\sigma) - r)} + 1\right)^{-1},\f]
 * where \f$r = \mu/\sigma\f$ to give a more familiar form of the function. Given how it is used the function
 * does not actually compute the normalisation factor in the prior.
 */
REAL8 LALInferenceFermiDiracPrior( LALInferenceVariables *priorArgs, const char *name, REAL8 value ){
  if ( !LALInferenceCheckFermiDiracPrior( priorArgs, name ) ){
    XLAL_ERROR_REAL8( XLAL_EINVAL, "No Fermi-Dirac prior given for parameter '%s'", name);
  }

  REAL8 r = 0., sigma = 0.;
  LALInferenceGetFermiDiracPrior(priorArgs, name, &sigma, &r);

  if ( value < 0. ){ return -INFINITY; } /* value must be positive */
  else{ return -logaddexp((value/sigma)-r, 0.); } /* log of Fermi-Dirac distribution (normalisation not required) */
}

/* Return the log Prior for a Gaussian Mixture Model given a value - the function returns zero
 * unless all required parameter values for a given GMM prior have already been passed to the
 * function */
REAL8 LALInferenceGMMPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 value){
  REAL8 logPrior = 0.0;

  if( !LALInferenceCheckGMMPrior( priorArgs, name ) ){
    XLAL_ERROR_REAL8( XLAL_EINVAL, "No Gaussian Mixture Model prior given for parameter '%s'", name);
  }

  char gmmParsName[VARNAME_MAX] = "gmm_parameter_lists"; // contains list of ':'-separated lists of GMM parameters

  /* check if the given variable is in any of the lists */
  LALStringVector *parLists = *(LALStringVector **)LALInferenceGetVariable( priorArgs, gmmParsName );
  UINT4 numlists = parLists->length, npars = 0;

  UINT4 found = 0, i = 0, j = 0;
  TokenList *toks = NULL;
  for ( i = 0; i < numlists; i++ ){
    if ( XLALCreateTokenList( &toks, parLists->data[i], ":" ) != XLAL_SUCCESS ){ XLAL_ERROR_REAL8(XLAL_EFAILED); }
    for ( j = 0; j < toks->nTokens; j++ ){
      if ( !strcmp(toks->tokens[j], name) ){
        found = 1;
        break;
      }
    }
    if ( found ){ break; }
  }

  /* set value of (normalised) current parameter in priorArgs */
  CHAR parname[VARNAME_MAX];
  sprintf(parname, "%s_gmm_value", name);
  LALInferenceAddVariable( priorArgs, parname, &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED );

  /* check if other required parameter values have already been set */
  REAL8Vector *allvalues = NULL;
  npars = toks->nTokens;
  allvalues = XLALCreateREAL8Vector( npars );
  for ( j = 0; j < npars; j++ ){
    sprintf(parname, "%s_gmm_value", toks->tokens[j]);
    if ( LALInferenceCheckVariable( priorArgs, parname ) ){
      allvalues->data[j] = LALInferenceGetREAL8Variable( priorArgs, parname );
    }
    else{ /* not all values are available yet, so return zero */
      XLALDestroyREAL8Vector( allvalues );
      return 0.;
    }
  }

  REAL8Vector **gmmsigmas = NULL, **gmmmus = NULL, *gmmweights = NULL, *gmmdets = NULL;
  gsl_matrix **cor, **invcor;
  REAL8Vector *gmmlow = NULL, *gmmhigh = NULL;
  UINT4 idx = 0;
  CHAR *fullname = NULL;

  /* get GMM parameters */
  LALInferenceGetGMMPrior( priorArgs, name, &gmmmus, &gmmsigmas, &cor, &invcor, &gmmweights, &gmmlow, &gmmhigh, &gmmdets, &idx, &fullname );

  /* check values are within limits */
  for ( j = 0; j < npars; j++ ){
    if ( allvalues->data[j] < gmmlow->data[j] || allvalues->data[j] > gmmhigh->data[j] ){
      logPrior = -INFINITY;
      break;
    }
  }

  if ( isfinite( logPrior ) ){
    REAL8 thisGauss = 0.;
    logPrior = -INFINITY;
    for ( i = 0; i < gmmweights->length; i++ ){
      thisGauss = 0.;
      gsl_vector *vmu = gsl_vector_calloc( npars ); /* vector to contain value-mu */
      gsl_matrix *invcov = gsl_matrix_calloc( npars, npars ); /* inverse covariance matrix */
      gsl_matrix *tmpMat = gsl_matrix_calloc( npars, npars );
      for ( j = 0; j < npars; j++ ){
        // normalise values to be from zero mean unit variance Gaussian
        gsl_vector_set(vmu, j, (allvalues->data[j]-gmmmus[i]->data[j])/gmmsigmas[i]->data[j]);
      }

      /* calculate log probability */
      gsl_vector *tmpVec = gsl_vector_calloc( npars );
      gsl_blas_dgemv (CblasNoTrans, 1.0, invcor[i], vmu, 0.0, tmpVec);
      gsl_blas_ddot( vmu, tmpVec, &thisGauss );
      thisGauss *= -0.5;

      thisGauss += log(gmmweights->data[i]);
      thisGauss -= (0.5*(REAL8)npars*(LAL_LNPI + LAL_LN2) + log(gmmdets->data[i])); /* normalisation */
      logPrior = logaddexp(logPrior, thisGauss); /* sum Gaussianians */

      gsl_matrix_free( invcov );
      gsl_matrix_free( tmpMat );
      gsl_vector_free( vmu );
      gsl_vector_free( tmpVec );
    }
  }

  /* remove all values that have been set */
  for ( j = 0; j < npars; j++ ){
    sprintf(parname, "%s_gmm_value", toks->tokens[j]);
    LALInferenceRemoveVariable( priorArgs, parname );
  }
  XLALDestroyTokenList( toks );
  XLALDestroyREAL8Vector( allvalues );

  return logPrior;
}

/* Return the log Prior for a parameter that has a prior that is uniform in log space */
REAL8 LALInferenceLogUniformPrior( LALInferenceVariables *priorArgs, const char *name, REAL8 value ){
  if ( !LALInferenceCheckLogUniformPrior( priorArgs, name ) ){
    XLAL_ERROR_REAL8( XLAL_EINVAL, "No log uniform prior given for parameter '%s'", name);
  }

  REAL8 min = 0., max = 0., lrat = 0.;
  LALInferenceGetLogUniformPrior( priorArgs, name, &min, &max );

  if ( value < 0. || value < min || value > max ){ return -INFINITY; }
  lrat = log(max/min);

  return -log(value*lrat);
}
