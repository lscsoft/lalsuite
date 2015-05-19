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
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Private helper function prototypes */
static double qInnerIntegrand(double M2, void *viData);
static double etaInnerIntegrand(double M2, void *viData);
static double outerIntegrand(double M1, void *voData);

void LALInferenceInitCBCPrior(LALInferenceRunState *runState)
{
    char help[]="\
                 ------------------------------------------------------------------------------------------------------------------\n\
                 --- Prior Arguments     ------------------------------------------------------------------------------------------\n\
                 ------------------------------------------------------------------------------------------------------------------\n\
                 (--malmquistprior)               Impose selection effects on the prior (False).\n\
                 (--malmquist-loudest-snr)        Threshold SNR in the loudest detector (0.0).\n\
                 (--malmquist-second-loudest-snr) Threshold SNR in the second loudest detector (5.0).\n\
                 (--malmquist-network-snr)        Threshold network SNR (0.0).\n\
                 (--analyticnullprior)            Use analytic null prior.\n\
                 (--nullprior)                    Use null prior in the sampled parameters.\n";
    ProcessParamsTable *commandLine=runState->commandLine;
    ProcessParamsTable *ppt = NULL;

    /* Print command line arguments if help requested */
    if(LALInferenceGetProcParamVal(commandLine, "--help"))
    {
        fprintf(stdout,"%s",help);
        return;
    }

    /* Choose the proper prior */
    if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
        runState->prior = &LALInferenceAnalyticNullPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
        runState->prior = &LALInferenceNullPrior;
    } else {
        runState->prior = &LALInferenceInspiralPrior;
    }

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
}


static REAL8 LALInferenceSplineCalibrationPrior(LALInferenceRunState *runState, LALInferenceVariables *params) {
  LALInferenceIFOData *ifo = NULL;
  REAL8 ampWidth = -1.0;
  REAL8 phaseWidth = -1.0;
  REAL8 logPrior = 0.0;

  if (runState->commandLine == NULL || !(LALInferenceGetProcParamVal(runState->commandLine, "--enable-spline-calibration"))) {
    return logPrior;
  }

  ampWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "spcal_amp_uncertainty");
  phaseWidth = *(REAL8 *)LALInferenceGetVariable(runState->priorArgs, "spcal_phase_uncertainty");

  ifo = runState->data;
  do {
    size_t i;

    char ampVarName[VARNAME_MAX];
    char phaseVarName[VARNAME_MAX];

    REAL8Vector *amps = NULL;
    REAL8Vector *phase = NULL;

    snprintf(ampVarName, VARNAME_MAX, "%s_spcal_amp", ifo->name);
    snprintf(phaseVarName, VARNAME_MAX, "%s_spcal_phase", ifo->name);

    amps = *(REAL8Vector **)LALInferenceGetVariable(params, ampVarName);
    phase = *(REAL8Vector **)LALInferenceGetVariable(params, phaseVarName);

    for (i = 0; i < amps->length; i++) {
      logPrior += -log(2.0*M_PI) - log(ampWidth) - 0.5*amps->data[i]*amps->data[i]/ampWidth/ampWidth;
      logPrior += -log(2.0*M_PI) - log(phaseWidth) - 0.5*phase->data[i]*phase->data[i]/phaseWidth/phaseWidth;
    }

    ifo = ifo->next;
  } while (ifo);

  return logPrior;
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

    char priormin[100];
    char priormax[100];

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
            if(val<component_min || val>component_max) return -DBL_MAX;
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
        if(val < component_min || val > component_max) return -DBL_MAX;
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

  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=-INFINITY, max=INFINITY;
  REAL8 mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;

  /* check if signal model is being used */
  UINT4 signalFlag=1;
  if(LALInferenceCheckVariable(params, "signalModelFlag"))
    signalFlag = *((INT4 *)LALInferenceGetVariable(params, "signalModelFlag"));

  if(signalFlag){

  /* Check boundaries for signal model parameters */
  for(;item;item=item->next)
  {
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else if (LALInferenceCheckMinMaxPrior(priorParams, item->name))
    {
      if(item->type==LALINFERENCE_REAL8_t){
        LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
        if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
      }
    }
  }
  if(LALInferenceCheckVariable(params,"flow"))
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"flow"));

  if(LALInferenceCheckVariable(params,"logdistance"))
    logPrior+=3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
  else if(LALInferenceCheckVariable(params,"distance"))
    logPrior+=2.0*log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));
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

  /* Check for component masses in range, if specified */
  if(LALInferenceCheckVariable(priorParams,"component_min"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m1
       || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m2)
      return -DBL_MAX;

  if(LALInferenceCheckVariable(priorParams,"component_max"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m1
       || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m2)
      return -DBL_MAX;

  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      return -DBL_MAX;

  if(model != NULL &&
        LALInferenceCheckVariable(priorParams,"malmquist") &&
        *(UINT4 *)LALInferenceGetVariable(priorParams,"malmquist") &&
        !within_malmquist(runState, params, model))
      return -DBL_MAX;

  }/* end prior for signal model parameters */

  /* Calibration priors. */
  logPrior += LALInferenceSplineCalibrationPrior(runState, params);

  /* Evaluate PSD prior (returns 0 if no PSD model) */
  logPrior += LALInferencePSDPrior(runState, params);

  /* Evaluate glitch prior (returns 0 if no glitch model) */
  logPrior += LALInferenceGlitchPrior(runState, params);

  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context)
{
    REAL8 min=-INFINITY, max=INFINITY, logPrior=0.;
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
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
          Cube[i] = azimuth;
          i++;
          strcat(header, "azimuth ");
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
          Cube[i] = cosalpha;
          i++;
          strcat(header, "cosalpha ");
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
            Cube[i] = lat;
            logPrior += log(fabs(cos(lat)));
            i++;
            strcat(header,"dec ");
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
            Cube[i] = longitude;
            i++;
            strcat(header,"RA ");
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
            Cube[i] = phi;
            i++;
            strcat(header,"phi_orb ");
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
            Cube[i] = psi;
            i++;
            strcat(header,"psi ");
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
            Cube[i] = t0;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"t0 ");
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
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"time ");
        }
      }
    }

    /***** OLD INCORRECT MASS PRIOR
    // mchirp
    double mc=0.0;
    if( LALInferenceCheckVariable(params,"logmc") )
    {
        item = LALInferenceGetItem(params, "logmc");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "logmc", (void *)&min, (void *)&max);
            min = exp(min); max = exp(max);
            mc = LALInferenceCubeToPowerPrior(-5.0/6.0, Cube[i], min, max);
            double logmc = log(mc);
            LALInferenceSetVariable(params, "logmc", &logmc);
            Cube[i] = mc;
            logPrior -= 5.0/6.0*logmc;
            i++;
            strcat(header,"mc ");
        }
        else
        {
            mc = exp(*(REAL8 *)LALInferenceGetVariable(params, "logmc"));
        }
    }
    else if( LALInferenceCheckVariable(params,"chirpmass") )
    {
        item = LALInferenceGetItem(params, "chirpmass");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "chirpmass", (void *)&min, (void *)&max);
            mc = LALInferenceCubeToPowerPrior(-5.0/6.0, Cube[i], min, max);
            LALInferenceSetVariable(params, "chirpmass", &mc);
            Cube[i] = mc;
            logPrior -= 5.0/6.0*log(mc);
            i++;
            strcat(header,"mc ");
        }
        else
        {
            mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
        }
    }

    // eta
    double eta=0.25, q=1.0, m1=0.0, m2=0.0;
    if( LALInferenceCheckVariable(params,"eta") )
    {
        item = LALInferenceGetItem(params, "eta");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "eta", (void *)&min, (void *)&max);
            eta = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            Cube[i]=eta;
            i++;
            strcat(header,"eta ");
            LALInferenceSetVariable(params, "eta", &eta);
        }
        else
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "eta");
        }
        LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
    }
    else if( LALInferenceCheckVariable(params,"q") )
    {
        item = LALInferenceGetItem(params, "q");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            min = 0.0, max = 0.25;
            eta = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
            q=m2/m1;
            Cube[i]=q;
            i++;
            strcat(header,"q ");
            LALInferenceSetVariable(params, "q", &q);
        }
        else
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "q");
            LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        }
    }
    ****/

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
        min = *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min");
        max = *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max");
        if( min == max )
        {
            m1 = m2 = min;
            m = m1 + m2;
            eta = m1 * m2 / (m*m);
            mc = pow(eta,0.6) * m;
            q = m2 / m1; // asymmetric mass ratio, m1 >= m2
        }
        else
        {
            m1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i+1], min, max);
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
            Cube[i] = mc; i++;
            strcat(header,"mc ");
            Cube[i] = eta; i++;
            strcat(header,"eta ");
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
            Cube[i] = dist;
            logPrior += 2.0*logdist;
            i++;
            strcat(header,"dist ");
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
            Cube[i] = dist;
            logPrior += 2.0*log(dist);
            i++;
            strcat(header,"dist ");
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
            Cube[i] = a_spin1;
            i++;
            strcat(header,"a1 ");
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
            Cube[i] = a_spin2;
            i++;
            strcat(header,"a2 ");
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
            Cube[i] = costheta_JN;
            i++;
            strcat(header,"costheta_jn ");
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
            Cube[i] = phi_JL;
            i++;
            strcat(header,"phi_jl ");
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
            Cube[i] = tilt_spin1;
            i++;
            strcat(header,"tilt1 ");
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
            Cube[i] = tilt_spin2;
            i++;
            strcat(header,"tilt2 ");
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
            Cube[i] = phi12;
            i++;
            strcat(header,"phi12 ");
        }
    }

    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);

    Cube[i] = m1; i++; strcat(header,"m1 ");
    Cube[i] = m2; i++; strcat(header,"m2 ");

    if(SKY_FRAME==1)
    {
        LALInferenceDetFrameToEquatorial(runState->data->detector, runState->data->next->detector,
                                         t0, acos(cosalpha), azimuth, &tc, &longitude, &lat);
        Cube[i] = tc; i++; strcat(header, "time ");
        Cube[i] = longitude; i++; strcat(header, "RA ");
        Cube[i] = lat; i++; strcat(header, "dec ");
    }

    Cube[i] = LALInferenceInspiralPrior(runState,params,model);
    i++; strcat(header,"logprior ");

    // f_ref for system-frame parameters
    if(LALInferenceCheckVariable(params,"f_ref"))
    {
        Cube[i] = *(REAL8 *)LALInferenceGetVariable(params,"f_ref");
        i++;
        strcat(header,"f_ref ");
    }

    strcat(header,"logl");

    /* Check boundaries */
    if (ScaleTest==0) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else if (item->type != LALINFERENCE_gslMatrix_t)
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

    /* Check for component masses in range, if specified */
    if(LALInferenceCheckVariable(priorParams,"component_min"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m1
           || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"component_max"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m1
           || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMin"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMax"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
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
  if (parameter == NULL || priorArgs == NULL)
    XLAL_ERROR_VOID(XLAL_EFAULT, "Null arguments received.");

  /* Apply cyclic and reflective boundaries to parameter to bring it back
     within the prior */
  LALInferenceVariableItem *paraHead=NULL;
  REAL8 min,max;
  for (paraHead=parameter->head;paraHead;paraHead=paraHead->next) {
    if( paraHead->vary==LALINFERENCE_PARAM_FIXED ||
        paraHead->vary==LALINFERENCE_PARAM_OUTPUT ||
        !LALInferenceCheckMinMaxPrior(priorArgs, paraHead->name) ) continue;

    LALInferenceGetMinMaxPrior(priorArgs,paraHead->name, &min, &max);

  /* Check that the minimum and maximum make sense. */
  if (min >= max) {
    XLAL_ERROR_VOID(XLAL_EINVAL, "Minimum %f for variable '%s' is not less than maximum %f.", min, paraHead->name, max);
  }

    if(paraHead->vary==LALINFERENCE_PARAM_CIRCULAR) {
      /* For cyclic boundaries, mod out by range. */

      REAL8 val = *(REAL8 *)paraHead->value;

      if (val > max) {
        REAL8 offset = val - min;
        REAL8 delta = max-min;

        *(REAL8 *)paraHead->value = min + fmod(offset, delta);
      } else {
        REAL8 offset = max - val;
        REAL8 delta = max - min;

        *(REAL8 *)paraHead->value = max - fmod(offset, delta);
      }
    } else if (paraHead->vary==LALINFERENCE_PARAM_LINEAR && paraHead->type==LALINFERENCE_REAL8_t) {
      /* For linear boundaries, reflect about endpoints of range until
         withoun range.
         SKIP NOISE PARAMETERS (ONLY CHECK REAL8) */
      while(1) {
        /* Loop until broken. */
        REAL8 val = *(REAL8 *)paraHead->value;
        if (val > max) {
          /* val <-- max - (val - max) */
          *(REAL8 *)paraHead->value = 2.0*max - val;
        } else if (val < min) {
          /* val <-- min + (min - val) */
          *(REAL8 *)paraHead->value = 2.0*min - val;
        } else {
          /* In range, so break. */
          break;
        }
      }
    }
  }
  return;
}


/**
 * \brief Rotate initial phase if polarisation angle is cyclic around ranges
 * If the polarisation angle parameter \f$\psi\f$ is cyclic about its upper and
 * lower ranges of \f$-\pi/4\f$ to \f$\pi/4\f$ then the transformation for
 * crossing a boundary requires the initial phase parameter \f$\phi_0\f$ to be
 * rotated through \f$\pi\f$ radians. The function assumes the value of
 * \f$\psi\f$ has been rescaled to be between 0 and \f$2\pi\f$ - this is a
 * requirement of the covariance matrix routine \c LALInferenceNScalcCVM
 * function.
 *
 * This is particularly relevant for pulsar analyses.
 *
 * \param parameter [in] Pointer to an array of parameters
 */
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
  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=-INFINITY, max=INFINITY;
  REAL8 logmc=0.0,mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;

  if (!SkyLocPriorWarning ) {
    SkyLocPriorWarning  = 1;
    fprintf(stderr, "SkyLocalization priors are being used. (in %s, line %d)\n", __FILE__, __LINE__);
  }
  /* Check boundaries */
  for(;item;item=item->next)
  {
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else
    {
            val = 0.0;
            min =-DBL_MAX;
            max = DBL_MAX;
            if(strcmp(item->name,"psdscale"))
            {
                LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
                val = *(REAL8 *)item->value;
            }

            if(val<min || val>max) return -DBL_MAX;
    }
  }

  /*Use a uniform in log D distribution*/

  if(LALInferenceCheckVariable(params,"flow"))
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"flow"));

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

  /* Check for component masses in range, if specified */
  if(LALInferenceCheckVariable(priorParams,"component_min"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m1
       || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m2)
      return -DBL_MAX;

  if(LALInferenceCheckVariable(priorParams,"component_max"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m1
       || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m2)
      return -DBL_MAX;

  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      return -DBL_MAX;

  if(LALInferenceCheckVariable(priorParams,"MTotMin"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
      return -DBL_MAX;

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
        if(val < min || val > max) return -DBL_MAX;
        else if(psdGaussianPrior)prior += -0.5*( (mean-val)*(mean-val)/var + log(2.0*LAL_PI*var) );
      }
    }
    logPrior+=prior;
  }

  /* Calibration parameters */
  logPrior += LALInferenceSplineCalibrationPrior(runState, params);

  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralSkyLocCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context)
{
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];
    REAL8 min=-INFINITY, max=INFINITY, logPrior=0.;
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
          Cube[i] = azimuth;
          i++;
          strcat(header, "azimuth ");
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
          Cube[i] = cosalpha;
          i++;
          strcat(header, "cosalpha ");
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
            Cube[i] = lat;
            logPrior += log(fabs(cos(lat)));
            i++;
            strcat(header,"dec ");
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
            Cube[i] = longitude;
            i++;
            strcat(header,"RA ");
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
            Cube[i] = phi;
            i++;
            strcat(header,"phi_orb ");
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
            Cube[i] = psi;
            i++;
            strcat(header,"psi ");
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
            Cube[i] = t0;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"t0 ");
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
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"time ");
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
        min = *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min");
        max = *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max");
        if( min == max )
        {
            m1 = m2 = min;
            m = m1 + m2;
            eta = m1 * m2 / (m*m);
            mc = pow(eta,0.6) * m;
            q = m2 / m1; // asymmetric mass ratio, m1 >= m2
        }
        else
        {
            m1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i+1], min, max);
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
            Cube[i] = mc; i++;
            strcat(header,"mc ");
            Cube[i] = eta; i++;
            strcat(header,"eta ");
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
            Cube[i] = dist;
            logPrior -= logdist;
            i++;
            strcat(header,"dist ");
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
            Cube[i] = dist;
            logPrior -= log(dist);
            i++;
            strcat(header,"dist ");
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
            Cube[i] = a_spin1;
            i++;
            strcat(header,"a1 ");
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
            Cube[i] = a_spin2;
            i++;
            strcat(header,"a2 ");
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
            Cube[i] = costheta_JN;
            i++;
            strcat(header,"costheta_jn ");
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
            Cube[i] = phi_JL;
            i++;
            strcat(header,"phi_jl ");
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
            Cube[i] = tilt_spin1;
            i++;
            strcat(header,"tilt1 ");
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
            Cube[i] = tilt_spin2;
            i++;
            strcat(header,"tilt2 ");
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
            Cube[i] = phi12;
            i++;
            strcat(header,"phi12 ");
        }
    }

    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);

    Cube[i] = m1; i++; strcat(header,"m1 ");
    Cube[i] = m2; i++; strcat(header,"m2 ");

    if (SKY_FRAME==1)
    {
        LALInferenceDetFrameToEquatorial(runState->data->detector, runState->data->next->detector,
                                         t0, acos(cosalpha), azimuth, &tc, &longitude, &lat);
        Cube[i] = tc; i++; strcat(header, "time ");
        Cube[i] = longitude; i++; strcat(header, "RA ");
        Cube[i] = lat; i++; strcat(header, "dec ");
    }

    Cube[i] = LALInferenceInspiralSkyLocPrior(runState,params,model);
    i++; strcat(header,"logprior ");

    // f_ref for system-frame parameters
    if(LALInferenceCheckVariable(params,"f_ref"))
    {
        Cube[i] = *(REAL8 *)LALInferenceGetVariable(params,"f_ref");
        i++;
        strcat(header,"f_ref ");
    }

    strcat(header,"logl");

    /* Check boundaries */
    if (ScaleTest==0) return 0;
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

    /* Check for component masses in range, if specified */
    if(LALInferenceCheckVariable(priorParams,"component_min"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m1
           || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_min") > m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"component_max"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m1
           || *(REAL8 *)LALInferenceGetVariable(priorParams,"component_max") < m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMin"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
            return 0;

    if(LALInferenceCheckVariable(priorParams,"MTotMax"))
        if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
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
  char minName[VARNAME_MAX];
  char maxName[VARNAME_MAX];
  sprintf(minName,"%s_min",name);
  sprintf(maxName,"%s_max",name);

  return (LALInferenceCheckVariable(priorArgs,minName) && LALInferenceCheckVariable(priorArgs,maxName));
}

/* Get the min and max values of the prior from the priorArgs list, given a name */
void LALInferenceGetMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, REAL8 *min, REAL8 *max)
{
    char minName[VARNAME_MAX];
    char maxName[VARNAME_MAX];
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
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);
  return (LALInferenceCheckVariable(priorArgs,meanName) && LALInferenceCheckVariable(priorArgs,sigmaName));
}

/* Function to add the min and max values for the prior onto the priorArgs */
void LALInferenceAddGaussianPrior( LALInferenceVariables *priorArgs,
                                   const char *name, REAL8 *mu, REAL8 *sigma,
                                   LALInferenceVariableType type ){
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];

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
                                    UINT4 *idx){
  char corName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"%s_correlation_matrix",name);
  sprintf(idxName,"%s_index",name);

  LALInferenceAddVariable(priorArgs, corName, cor,
                          LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(priorArgs, idxName, idx, LALINFERENCE_UINT4_t,
                          LALINFERENCE_PARAM_FIXED) ;
  return;
}

/* Get the correlation matrix and parameter index */
void LALInferenceGetCorrelatedPrior(LALInferenceVariables *priorArgs,
                                    const char *name, gsl_matrix **cor,
                                    UINT4 *idx){
  char corName[VARNAME_MAX];
  char idxName[VARNAME_MAX];
  void *ptr = NULL;

  sprintf(corName,"%s_correlation_matrix",name);
  sprintf(idxName,"%s_index",name);

  ptr = LALInferenceGetVariable(priorArgs, corName);
  if ( ptr ) *cor = *(gsl_matrix **)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  ptr = LALInferenceGetVariable(priorArgs, idxName);
  if ( ptr ) *idx = *(UINT4 *)ptr;
  else XLAL_ERROR_VOID(XLAL_EFAILED);

  return;
}

/* Remove the correlated prior */
void LALInferenceRemoveCorrelatedPrior(LALInferenceVariables *priorArgs,
                                       const char *name){
  char corName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"%s_correlation_matrix",name);
  sprintf(idxName,"%s_index",name);

  LALInferenceRemoveVariable(priorArgs, corName);
  LALInferenceRemoveVariable(priorArgs, idxName);
  return;
}

/* Check for a correlated prior of the standard form */
int LALInferenceCheckCorrelatedPrior(LALInferenceVariables *priorArgs,
                                     const char *name){
  char corName[VARNAME_MAX];
  char idxName[VARNAME_MAX];

  sprintf(corName,"%s_correlation_matrix",name);
  sprintf(idxName,"%s_index",name);

  return (LALInferenceCheckVariable(priorArgs,corName) &&
          LALInferenceCheckVariable(priorArgs,idxName));
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
  /* test for a prior drawn from correlated values */
  else if( LALInferenceCheckCorrelatedPrior( priorArgs, name ) ){
    gsl_matrix *cor = NULL;
    UINT4 idx = 0, dims = 0;
    REAL4Vector *tmps = NULL;

    LALInferenceGetCorrelatedPrior( priorArgs, name, &cor, &idx );
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

      /* draw values from the multivariate Gaussian described by the correlation
         matrix */
      tmps = XLALCreateREAL4Vector( dims );
      randParam = XLALCreateRandomParams( randomseed );
      XLALMultiNormalDeviates( tmps, cor, dims, randParam );

      LALInferenceAddVariable( priorArgs, "multivariate_deviates", &tmps,
                               LALINFERENCE_REAL8Vector_t,
                               LALINFERENCE_PARAM_FIXED );
    }

    /* set random number for given parameter index */
    tmp = tmps->data[idx];

    /* free tmps */
    if ( !LALInferenceCheckVariable( priorArgs, "multivariate_deviates" ) )
      XLALDestroyREAL4Vector( tmps );
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
    case LALINFERENCE_gslMatrix_t:
    {
      REAL8 val = tmp;
      LALInferenceSetVariable(output, name, &val);
      break;
    }
    default:
      XLAL_ERROR_VOID(XLAL_EINVAL, "Trying to randomise a non-numeric parameter!");
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

UINT4 LALInferenceAnalyticCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, LALInferenceModel *model, double *Cube, void *context) {
    int i = 0;
    REAL8 logPrior=0.,min=-INFINITY,max=INFINITY;
    REAL8 m1=1.,m2=1.,mc,eta,m;
    LALInferenceVariableItem *item;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];

    INT4 SKY_FRAME=0;
    if(LALInferenceCheckVariable(params,"SKY_FRAME"))
      SKY_FRAME=*(INT4 *)LALInferenceGetVariable(params,"SKY_FRAME");

    if(LALInferenceCheckVariable(params,"m1"))
    {
        if(LALInferenceGetItem(params,"m1")->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "m1", (void *)&min, (void *)&max);
            m1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "m1", &m1);
            Cube[i] = m1; i++;
            logPrior -= log(max-min);
            strcat(header,"m1 ");
        }
        else
        {
            m1=(*(REAL8 *)LALInferenceGetVariable(params,"m1"));
        }
    }
    if(LALInferenceCheckVariable(params,"m2"))
    {
        if(LALInferenceGetItem(params,"m2")->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "m2", (void *)&min, (void *)&max);
            m2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "m2", &m2);
            Cube[i] = m2; i++;
            logPrior -= log(max-min);
            strcat(header,"m2 ");
        }
        else
        {
            m2=(*(REAL8 *)LALInferenceGetVariable(params,"m2"));
        }
    }
    m = m1 + m2;
    eta = m1 * m2 / (m*m);
    mc = pow(eta,0.6) * m;

    if(LALInferenceCheckVariable(params,"phase"))
    {
        item = LALInferenceGetItem(params,"phase");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phase", (void *)&min, (void *)&max);
            double phase = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phase", &phase);
            Cube[i] = phase;
            i++;
            logPrior -= log(max-min);
            strcat(header,"phi_orb ");
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
            Cube[i] = polarisation;
            i++;
            logPrior -= log(max-min);
            strcat(header,"psi ");
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
            Cube[i] = azimuth;
            i++;
            logPrior -= log(max-min);
            strcat(header,"azimuth ");
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
            Cube[i] = rightascension;
            i++;
            logPrior -= log(max-min);
            strcat(header,"RA ");
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
            Cube[i] = cosalpha;
            i++;
            logPrior -= log(max-min);
            strcat(header,"cosalpha ");
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
            Cube[i] = declination;
            i++;
            logPrior -= log(max-min);
            strcat(header,"dec ");
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
            Cube[i] = distance;
            i++;
            logPrior -= log(max-min);
            strcat(header,"dist ");
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
            Cube[i] = t0;
            sprintf(timeID,"%d",i);
            i++;
            logPrior -= log(max-min);
            strcat(header,"t0 ");
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
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            logPrior -= log(max-min);
            strcat(header,"time ");
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
            Cube[i] = a_spin1;
            i++;
            logPrior -= log(max-min);
            strcat(header,"a1 ");
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
            Cube[i] = a_spin2;
            i++;
            logPrior -= log(max-min);
            strcat(header,"a2 ");
        }
    }

    LALInferenceVariables *priorParams=runState->priorArgs;
    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);

    Cube[i] = mc; i++; strcat(header,"mc ");
    Cube[i] = eta; i++; strcat(header,"eta ");
    Cube[i] = LALInferenceAnalyticNullPrior(runState,params,model);
    i++; strcat(header,"logprior ");

    strcat(header,"logl");

    if (ScaleTest==0) return 0;
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
      if (min>*(REAL8 *)cur->value || max<*(REAL8 *)cur->value) return -DBL_MAX;
    }
  }
  return 0.0;
}

REAL8 LALInferenceNullPrior(LALInferenceRunState UNUSED *runState, LALInferenceVariables UNUSED *params, LALInferenceModel UNUSED *model) {
  return 0.0;
}

UINT4 LALInferenceCubeToPSDScaleParams(LALInferenceVariables *priorParams, LALInferenceVariables *params, INT4 *idx, double *Cube, void *context)
{
  char **info = (char **)context;
  char *header = &info[1][0];
  char tempstr[50];

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
        Cube[(*idx)] = val;
        (*idx)++;
        sprintf(tempstr,"psdscale_%d_%d ",i,j);
        strcat(header,tempstr);
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
