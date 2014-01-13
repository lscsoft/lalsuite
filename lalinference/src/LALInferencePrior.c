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

/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
  if (params == NULL || runState == NULL || runState->priorArgs == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT, "Null arguments received.");

  REAL8 logPrior=0.0;

  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=-INFINITY, max=INFINITY;
  REAL8 mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;
  //REAL8 tmp=0.;
  /* Check boundaries */
  for(;item;item=item->next)
  {
    // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
    if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
    else
    {
      LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
      if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
    }
  }
  if(LALInferenceCheckVariable(params,"fLow"))
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"fLow"));

  if(LALInferenceCheckVariable(params,"logdistance"))
    logPrior+=3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
  else if(LALInferenceCheckVariable(params,"distance"))
    logPrior+=2.0*log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));

  if(LALInferenceCheckVariable(params,"inclination"))
    logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))));
  if(LALInferenceCheckVariable(params,"theta_JN"))
    logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_JN"))));
  if(LALInferenceCheckVariable(params,"declination"))
    logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
  if(LALInferenceCheckVariable(params,"theta_spin1"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"theta_spin1");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
     logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))));
    }
  }
  if(LALInferenceCheckVariable(params,"theta_spin2"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"theta_spin2");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))));
    }
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

  if(LALInferenceCheckVariable(params,"asym_massratio")) {
    q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
    LALInferenceMcQ2Masses(mc,q,&m1,&m2);
  } else if(LALInferenceCheckVariable(params,"massratio")) {
    eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
    LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
  }

  if(LALInferenceCheckVariable(params,"logmc")) {
    if(LALInferenceCheckVariable(params,"asym_massratio"))
      logPrior+=log(m1*m1);
    else
      logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
  } else if(LALInferenceCheckVariable(params,"chirpmass")) {
    if(LALInferenceCheckVariable(params,"asym_massratio"))
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

  if(LALInferenceCheckVariable(priorParams,"malmquist") &&
        *(UINT4 *)LALInferenceGetVariable(priorParams,"malmquist") &&
        !within_malmquist(runState, params))
      return -DBL_MAX;

  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, double *Cube, void *context)
{
    REAL8 min=-INFINITY, max=INFINITY, logPrior=0.;
    LALInferenceVariableItem *item;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];
    int i = 0;

    // latitude
    if (LALInferenceCheckVariable(params, "declination"))
    {
        item = LALInferenceGetItem(params, "declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double lat = asin(2.0 * Cube[i] - 1.0);
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
            double longitude = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &longitude);
            Cube[i] = longitude;
            i++;
            strcat(header,"RA ");
        }
    }

    // inclination angle
    if (LALInferenceCheckVariable(params, "inclination"))
    {
        item = LALInferenceGetItem(params, "inclination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double iota;
            iota = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "inclination", &iota);
            Cube[i] = iota;
            logPrior += log(fabs(sin(iota)));
            i++;
            strcat(header,"iota ");
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
            strcat(header,"phi0 ");
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

    // time
    if (LALInferenceCheckVariable(params, "time"))
    {
        item = LALInferenceGetItem(params, "time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            double tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"time ");
        }
    }

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
            strcat(header,"mchirp ");
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
            strcat(header,"mchirp ");
        }
        else
        {
            mc = *(REAL8 *)LALInferenceGetVariable(params, "chirpmass");
        }
    }

    // eta
    double eta=0.25, q=1.0, m1=0.0, m2=0.0;
    if( LALInferenceCheckVariable(params,"massratio") )
    {
        item = LALInferenceGetItem(params, "massratio");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "massratio", (void *)&min, (void *)&max);
            eta = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            Cube[i]=eta;
            i++;
            strcat(header,"eta ");
            LALInferenceSetVariable(params, "massratio", &eta);
        }
        else
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "massratio");
        }
        LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
    }
    else if( LALInferenceCheckVariable(params,"asym_massratio") )
    {
        item = LALInferenceGetItem(params, "asym_massratio");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            min = 0.0, max = 0.25;
            eta = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
            q=m2/m1;
            Cube[i]=q;
            i++;
            strcat(header,"q ");
            LALInferenceSetVariable(params, "asym_massratio", &q);
        }
        else
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "asym_massratio");
            LALInferenceMcQ2Masses(mc,q,&m1,&m2);
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
    else if(LALInferenceCheckVariable(params,"spin1"))
    {
        item = LALInferenceGetItem(params, "spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin1", &a_spin1);
            Cube[i] = a_spin1;
            i++;
            strcat(header,"a1 ");
        }
    }

    // theta_spin1
    if(LALInferenceCheckVariable(params,"theta_spin1"))
    {
        item = LALInferenceGetItem(params, "theta_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin1;
            theta_spin1 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin1", &theta_spin1);
            Cube[i] = theta_spin1;
            logPrior += log(fabs(sin(theta_spin1)));
            i++;
            strcat(header,"theta1 ");
        }
    }

    // phi_spin1
    if(LALInferenceCheckVariable(params,"phi_spin1"))
    {
        item = LALInferenceGetItem(params, "phi_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin1", (void *)&min, (void *)&max);
            double phi_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin1", &phi_spin1);
            Cube[i] = phi_spin1;
            i++;
            strcat(header,"phi1 ");
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
    else if(LALInferenceCheckVariable(params,"spin2"))
    {
        item = LALInferenceGetItem(params, "spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin2", &a_spin2);
            Cube[i] = a_spin2;
            i++;
            strcat(header,"a2 ");
        }
    }

    // theta_spin2
    if(LALInferenceCheckVariable(params,"theta_spin2"))
    {
        item = LALInferenceGetItem(params, "theta_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin2;
            theta_spin2 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin2", &theta_spin2);
            Cube[i] = theta_spin2;
            logPrior += log(fabs(sin(theta_spin2)));
            i++;
            strcat(header,"theta2 ");
        }
    }

    // phi_spin2
    if(LALInferenceCheckVariable(params,"phi_spin2"))
    {
        item = LALInferenceGetItem(params, "phi_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin2", (void *)&min, (void *)&max);
            double phi_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin2", &phi_spin2);
            Cube[i] = phi_spin2;
            i++;
            strcat(header,"phi2 ");
        }

    }

    // theta_JN for system-frame parameters
    if(LALInferenceCheckVariable(params,"theta_JN"))
    {
        item = LALInferenceGetItem(params, "theta_JN");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "theta_JN", (void *)&min, (void *)&max);
            double theta_JN = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "theta_JN", &theta_JN);
            Cube[i] = theta_JN;
            i++;
            strcat(header,"theta_JN ");
        }
    }

    // phi_JL for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi_JL"))
    {
        item = LALInferenceGetItem(params, "phi_JL");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_JL", (void *)&min, (void *)&max);
            double phi_JL = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_JL", &phi_JL);
            Cube[i] = phi_JL;
            i++;
            strcat(header,"phi_JL ");
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
            strcat(header,"tilt_spin1 ");
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
            strcat(header,"tilt_spin2 ");
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

    LALInferenceVariables *priorParams=runState->priorArgs;
    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);

    Cube[i] = m1; i++; strcat(header,"m1 ");
    Cube[i] = m2; i++; strcat(header,"m2 ");
    Cube[i] = LALInferenceInspiralPrior(runState,params);
    i++; strcat(header,"logprior ");

    // fRef for system-frame parameters
    if(LALInferenceCheckVariable(params,"fRef"))
    {
        Cube[i] = *(REAL8 *)LALInferenceGetVariable(params,"fRef");
        i++;
        strcat(header,"f_ref ");
    }

    strcat(header,"logl");

    /* Check boundaries */
    if (ScaleTest==0) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else if (item->type != LALINFERENCE_gslMatrix_t)
        {
            LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
            if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return 0 ;
        }
    }

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

    if(LALInferenceCheckVariable(priorParams,"malmquist") &&
        *(UINT4 *)LALInferenceGetVariable(priorParams,"malmquist") &&
        !within_malmquist(runState, params))
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
  /* REAL8 mu, sigma; */
  for (paraHead=parameter->head;paraHead;paraHead=paraHead->next) {
    if( paraHead->vary==LALINFERENCE_PARAM_FIXED ||
        paraHead->vary==LALINFERENCE_PARAM_OUTPUT ||
        !LALInferenceCheckMinMaxPrior(priorArgs, paraHead->name) ) continue;

    LALInferenceGetMinMaxPrior(priorArgs,paraHead->name, &min, &max);

  /* Check that the minimum and maximum make sense. */
  if (min >= max) {
    XLAL_ERROR_VOID(XLAL_EINVAL, "Minimum %f for variable '%s' is not less than maximum %f.", paraHead->name, min, max);
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
    } else if (paraHead->vary==LALINFERENCE_PARAM_LINEAR) {
      /* For linear boundaries, reflect about endpoints of range until
         withoun range. */
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
REAL8 LALInferenceInspiralSkyLocPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
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
    // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
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
  //if(LALInferenceCheckVariable(params,"logdistance"))
  //  logPrior+=3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");

  if(LALInferenceCheckVariable(params,"fLow"))
    logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"fLow"));

  if(LALInferenceCheckVariable(params,"distance"))
    logPrior-=log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));
  if(LALInferenceCheckVariable(params,"inclination"))
    logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))));
  if(LALInferenceCheckVariable(params,"theta_JN"))
    logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_JN"))));
  if(LALInferenceCheckVariable(params,"declination"))
    logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
  if(LALInferenceCheckVariable(params,"theta_spin1"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"theta_spin1");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))));
    }
  }
  if(LALInferenceCheckVariable(params,"theta_spin2"))
  {
    LALInferenceParamVaryType vtype=LALInferenceGetVariableVaryType(params,"theta_spin2");
    if(vtype!=LALINFERENCE_PARAM_FIXED && vtype!=LALINFERENCE_PARAM_OUTPUT)
    {
      logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))));
    }
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
  if(LALInferenceCheckVariable(params,"massratio")||LALInferenceCheckVariable(params,"asym_massratio")) {
    if(LALInferenceCheckVariable(params,"logmc")) {
      logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        LALInferenceMcQ2Masses(exp(logmc),q,&m1,&m2);
        logPrior+=log(m1*m1);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
        LALInferenceMcEta2Masses(exp(logmc),eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
      }
      /*careful using LALInferenceMcEta2Masses, it returns m1>=m2*/
    } else if(LALInferenceCheckVariable(params,"chirpmass")) {
      mc=*(REAL8 *)LALInferenceGetVariable(params,"chirpmass");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        logPrior+=log(m1*m1/mc);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
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

    //REAL8 val;
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

  return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralSkyLocCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, double *Cube, void *context)
{
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];
    REAL8 min=-INFINITY, max=INFINITY, logPrior=0.;
    int i = 0;

    // latitude
    if (LALInferenceCheckVariable(params, "declination"))
    {
        item = LALInferenceGetItem(params, "declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double lat = asin(2.0 * Cube[i] - 1.0);
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
            double longitude = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &longitude);
            Cube[i] = longitude;
            i++;
            strcat(header,"RA ");
        }
    }

    // inclination angle
    if (LALInferenceCheckVariable(params, "inclination"))
    {
        item = LALInferenceGetItem(params, "inclination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double iota;
            iota = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "inclination", &iota);
            Cube[i] = iota;
            logPrior += log(fabs(sin(iota)));
            i++;
            strcat(header,"iota ");
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
            strcat(header,"phi0 ");
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

    // time
    if (LALInferenceCheckVariable(params, "time"))
    {
        item = LALInferenceGetItem(params, "time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            double tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"time ");
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
    if( LALInferenceCheckVariable(params,"massratio") )
    {
        item = LALInferenceGetItem(params, "massratio");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "massratio");
            if( mc != 0.0 ) LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        }
    }
    else if( LALInferenceCheckVariable(params,"asym_massratio") )
    {
        item = LALInferenceGetItem(params, "asym_massratio");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "asym_massratio");
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
            strcat(header,"mchirp ");
            Cube[i] = eta; i++;
            strcat(header,"eta ");
        }

        // chirp mass and eta/q
        if(LALInferenceCheckVariable(params,"massratio")||LALInferenceCheckVariable(params,"asym_massratio"))
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

                  if(LALInferenceCheckVariable(params,"asym_massratio"))
                LALInferenceSetVariable(params, "asym_massratio", &q);
            else if(LALInferenceCheckVariable(params,"massratio"))
                    LALInferenceSetVariable(params, "massratio", &eta);
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
    else if(LALInferenceCheckVariable(params,"spin1"))
    {
        item = LALInferenceGetItem(params, "spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin1", &a_spin1);
            Cube[i] = a_spin1;
            i++;
            strcat(header,"a1 ");
        }
    }

    // theta_spin1
    if(LALInferenceCheckVariable(params,"theta_spin1"))
    {
        item = LALInferenceGetItem(params, "theta_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin1;
            theta_spin1 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin1", &theta_spin1);
            Cube[i] = theta_spin1;
            logPrior += log(fabs(sin(theta_spin1)));
            i++;
            strcat(header,"theta1 ");
        }
    }

    // phi_spin1
    if(LALInferenceCheckVariable(params,"phi_spin1"))
    {
        item = LALInferenceGetItem(params, "phi_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin1", (void *)&min, (void *)&max);
            double phi_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin1", &phi_spin1);
            Cube[i] = phi_spin1;
            i++;
            strcat(header,"phi1 ");
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
    else if(LALInferenceCheckVariable(params,"spin2"))
    {
        item = LALInferenceGetItem(params, "spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin2", &a_spin2);
            Cube[i] = a_spin2;
            i++;
            strcat(header,"a2 ");
        }
    }

    // theta_spin2
    if(LALInferenceCheckVariable(params,"theta_spin2"))
    {
        item = LALInferenceGetItem(params, "theta_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin2;
            theta_spin2 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin2", &theta_spin2);
            Cube[i] = theta_spin2;
            logPrior += log(fabs(sin(theta_spin2)));
            i++;
            strcat(header,"theta2 ");
        }
    }

    // phi_spin2
    if(LALInferenceCheckVariable(params,"phi_spin2"))
    {
        item = LALInferenceGetItem(params, "phi_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin2", (void *)&min, (void *)&max);
            double phi_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin2", &phi_spin2);
            Cube[i] = phi_spin2;
            i++;
            strcat(header,"phi2 ");
        }

    }

    // theta_JN for system-frame parameters
    if(LALInferenceCheckVariable(params,"theta_JN"))
    {
        item = LALInferenceGetItem(params, "theta_JN");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "theta_JN", (void *)&min, (void *)&max);
            double theta_JN = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "theta_JN", &theta_JN);
            Cube[i] = theta_JN;
            i++;
            strcat(header,"theta_JN ");
        }
    }

    // phi_JL for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi_JL"))
    {
        item = LALInferenceGetItem(params, "phi_JL");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_JL", (void *)&min, (void *)&max);
            double phi_JL = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_JL", &phi_JL);
            Cube[i] = phi_JL;
            i++;
            strcat(header,"phi_JL ");
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
            strcat(header,"tilt_spin1 ");
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
            strcat(header,"tilt_spin2 ");
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
    Cube[i] = LALInferenceInspiralSkyLocPrior(runState,params);
    i++; strcat(header,"logprior ");

    // fRef for system-frame parameters
    if(LALInferenceCheckVariable(params,"fRef"))
    {
        Cube[i] = *(REAL8 *)LALInferenceGetVariable(params,"fRef");
        i++;
        strcat(header,"f_ref ");
    }

    strcat(header,"logl");

    /* Check boundaries */
    if (ScaleTest==0) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else
        {
            LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
            if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return 0 ;
        }
    }

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


/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPriorNormalised(LALInferenceRunState *runState, LALInferenceVariables *params)
{
  static int S6PEpriorWarning = 0;
  REAL8 logPrior=0.0;
  REAL8 val;
(void)runState;
LALInferenceVariableItem *item=params->head;
    LALInferenceVariables *priorParams=runState->priorArgs;
    REAL8 min=-INFINITY, max=INFINITY;
    REAL8 mc=0.0;
  REAL8 eta=0.0;
    REAL8 m1,m2;
    REAL8 massRatioMin=0.0, massRatioMax=0.0; // min,max for q or eta
    REAL8 MTotMax=0.0;
  REAL8 component_max, component_min;
    char normName[VARNAME_MAX];
    char massRatioName[VARNAME_MAX];
    REAL8 norm=0.0;

  UINT4 signalFlag=1;
  if(LALInferenceCheckVariable(params, "signalModelFlag"))
    signalFlag = *((UINT4 *)LALInferenceGetVariable(params, "signalModelFlag"));

  if (!S6PEpriorWarning) {
    S6PEpriorWarning = 1;
    fprintf(stderr, "S6PEpaper priors are being used. (in %s, line %d)\n", __FILE__, __LINE__);
  }

  if(signalFlag){
    if(LALInferenceCheckVariable(params,"asym_massratio")){
      LALInferenceGetMinMaxPrior(priorParams, "asym_massratio", (void *)&massRatioMin, (void *)&massRatioMax);
      strcpy(massRatioName,"asym_massratio");
    }else{
      LALInferenceGetMinMaxPrior(priorParams, "massratio", (void *)&massRatioMin, (void *)&massRatioMax);
      strcpy(massRatioName,"massratio");
    }
  }
    /* Check boundaries */
    for(;item;item=item->next)
    {
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT) continue;
        if(item->type!=LALINFERENCE_REAL8_t && item->type!=LALINFERENCE_gslMatrix_t) continue;
        else
        {

            val = 0;
            min =-DBL_MAX;
            max = DBL_MAX;
            /* generic prior check won't work as-is for nonREAL8 noise parameters */
            if(item->type==LALINFERENCE_REAL8_t){
            if(LALInferenceCheckMinMaxPrior(priorParams,item->name))
            {
                LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
                val = *(REAL8 *)item->value;
            }
            }
            if(val<min || val>max) return -DBL_MAX;
            else
            {

                    if(!strcmp(item->name, "chirpmass") || !strcmp(item->name, "logmc")){
                      if( LALInferenceCheckVariable(priorParams,"component_max") && LALInferenceCheckVariable(priorParams,"component_min")
                        && LALInferenceCheckVariable(priorParams,"MTotMax")
                        && (LALInferenceCheckVariable(params,"asym_massratio") || LALInferenceCheckVariable(params,"massratio")) ){

                            MTotMax=*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax");
                            component_min=*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min");
                            component_max=*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max");

                            if(LALInferenceCheckVariable(priorParams,"mass_norm")) {
                                norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass_norm");
                            }else{
                                  if( MTotMax > component_max && MTotMax < 2.0*component_max - component_min ) {
                                     norm = -log( (pow(MTotMax-2.0*component_min,2)/4.0) - (pow(MTotMax-component_max-component_min,2)/2.0) );
                        }else if(MTotMax >= component_max - component_min){
                          norm = -log( pow(MTotMax-2.0*component_min,2)/4.0 );
                        }else if(2.0*MTotMax < component_max){
                          norm = -log( pow(component_max-component_min,2)/2.0 );
                        }else{
                          norm = 0.0; //no prior area for the masses !!
                        }
                                 LALInferenceAddVariable(priorParams, "mass_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                            }
                        logPrior+=norm;
                        if(!strcmp(item->name, "chirpmass")){
                              mc=(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass"));
                            }
                        else if(!strcmp(item->name, "logmc")){
                              mc=exp(*(REAL8 *)LALInferenceGetVariable(params,"logmc"));
                           }

                              if(LALInferenceCheckVariable(params,"asym_massratio"))
                                LALInferenceMcQ2Masses(mc,*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio"),&m1,&m2);
                              else if(LALInferenceCheckVariable(params,"massratio")){
                                eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
                                LALInferenceMcEta2Masses(mc,*(REAL8 *)LALInferenceGetVariable(params,"massratio"),&m1,&m2);
                        }
                           if(component_min > m1 || component_min > m2)
                              return -DBL_MAX;
                            if(component_max < m1 || component_max < m2)
                              return -DBL_MAX;
                           if(MTotMax < m1+m2)
                              return -DBL_MAX;

                           if(LALInferenceCheckVariable(params,"logmc")) {
                                if(LALInferenceCheckVariable(params,"asym_massratio")) {
                                      logPrior+=log(m1*m1);
                                } else {
                        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
                                    }
                           } else if(LALInferenceCheckVariable(params,"chirpmass")) {
                                 if(LALInferenceCheckVariable(params,"asym_massratio")) {
                                  logPrior+=log(m1*m1/mc);
                                    } else {
                                  logPrior+=log(((m1+m2)*(m1+m2))/((m1-m2)*pow(eta,3.0/5.0)));
                        }
                    }
                     }else{
                        fprintf(stderr,"ERROR; component_max, component_min and MTotMax required for this prior.\n");
                        exit(1);
                     }
                }
            else if(!strcmp(item->name, "massratio") || !strcmp(item->name, "asym_massratio")) continue;
                else if(!strcmp(item->name, "mass1")){
                      if( LALInferenceCheckVariable(priorParams,"component_max") && LALInferenceCheckVariable(priorParams,"component_min")
                         && LALInferenceCheckVariable(priorParams,"MTotMax")
                         && LALInferenceCheckVariable(params,"mass2") ){
                        m1=(*(REAL8 *)LALInferenceGetVariable(params,"mass1"));
                        m2=(*(REAL8 *)LALInferenceGetVariable(params,"mass2"));
                        MTotMax=*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax");
                        component_min=*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min");
                        component_max=*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max");
                        if(component_min > m1 || component_min > m2)
                          return -DBL_MAX;
                        if(component_max < m1 || component_max < m2)
                          return -DBL_MAX;
                         if(MTotMax < m1+m2)
                          return -DBL_MAX;

                        if(LALInferenceCheckVariable(priorParams,"mass_norm")) {
                          norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass_norm");
                        }else{
                          if( MTotMax < component_max || MTotMax > 2.0*component_max - component_min ) {
                                fprintf(stderr,"ERROR; MTotMax < component_max || MTotMax > 2.0*component_max - component_min\n");
                                fprintf(stderr,"MTotMax = %lf, component_max=%lf, component_min=%lf\n",MTotMax,component_min,component_max);
                                exit(1);
                          }
                          norm = -log( (pow(MTotMax-component_min,2.0)/4.0) - (pow(MTotMax-component_max,2.0)/2.0) );
                      //printf("norm@%s=%f\n",item->name,norm);
                          LALInferenceAddVariable(priorParams, "mass_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                            }
                            logPrior+=norm;
                       }else{
                        fprintf(stderr,"ERROR; mass2, component_max, component_min and MTotMax required for this prior.\n");
                        exit(1);
                    }
                 }
            else if(!strcmp(item->name, "distance")){
                    if(LALInferenceCheckVariable(priorParams,"distance_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"distance_norm");
                    }
                    else
                    {
                        norm = +1.09861228867-log(max*max*max-min*min*min);
                        LALInferenceAddVariable(priorParams, "distance_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += 2.0*log(*(REAL8 *)LALInferenceGetVariable(params,"distance"))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "logdistance")){
                    if(LALInferenceCheckVariable(priorParams,"logdistance_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"logdistance_norm");
                    }
                    else
                    {
                        norm = 1.38629436112-log(max*max*max*max-min*min*min*min);
                        LALInferenceAddVariable(priorParams, "logdistance_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += 3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance")+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "fLow")){
                    if(LALInferenceCheckVariable(priorParams,"fLow_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"fLow_norm");
                    }
                    else
                    {
                        norm = 0.69314718056-log(pow(max,2.)-pow(min,2.));
                        LALInferenceAddVariable(priorParams, "fLow_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
          logPrior+=log(*(REAL8 *)LALInferenceGetVariable(params,"fLow"));
        }


                else if(!strcmp(item->name, "inclination")){
                    if(LALInferenceCheckVariable(priorParams,"inclination_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"inclination_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "inclination_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "theta_JN")){
                    if(LALInferenceCheckVariable(priorParams,"theta_JN_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"theta_JN_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "theta_JN_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_JN"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "declination")){
                    if(LALInferenceCheckVariable(priorParams,"declination_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"declination_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( -sin(LAL_PI*fractpart_min)+sin(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "declination_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "theta_spin1")){
                    if(LALInferenceCheckVariable(priorParams,"theta_spin1_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"theta_spin1_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "theta_spin1_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "tilt_spin1")){
                    if(LALInferenceCheckVariable(priorParams,"tilt_spin1_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"tilt_spin1_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "tilt_spin1_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin1"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "theta_spin2")){
                    if(LALInferenceCheckVariable(priorParams,"theta_spin2_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"theta_spin2_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "theta_spin2_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
                else if(!strcmp(item->name, "tilt_spin2")){
                    if(LALInferenceCheckVariable(priorParams,"tilt_spin2_norm")) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"tilt_spin2_norm");
                    }
                    else
                    {
                        REAL8 intpart_min=0.0;
                        REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
                        REAL8 intpart_max=0.0;
                        REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
                        norm = -log( cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min) );
                        LALInferenceAddVariable(priorParams, "tilt_spin2_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"tilt_spin2"))))+norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }
        //PSD priors are Gaussian
                else if(!strcmp(item->name, "psdscale"))
        {
          UINT4 i;
          UINT4 j;

          //REAL8 val;
          REAL8 var;
          REAL8 mean = 1.0;
          REAL8 prior= 0.0;
          UINT4 psdGaussianPrior;

          REAL8Vector *sigma = *((REAL8Vector **)LALInferenceGetVariable(priorParams, "psdsigma"));
          gsl_matrix *nparams = *((gsl_matrix **)(item->value));

          component_min=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_min");
          component_max=*(REAL8 *)LALInferenceGetVariable(priorParams,"psdrange_max");

          psdGaussianPrior = *(UINT4 *)LALInferenceGetVariable(priorParams,"psdGaussianPrior");

          for(i=0; i<(UINT4)nparams->size1; i++)
          {
            for(j=0; j<(UINT4)nparams->size2; j++)
            {
              var = sigma->data[j]*sigma->data[j];
              val = gsl_matrix_get(nparams,i,j);

              //reject prior
              if(val < component_min || val > component_max) return -DBL_MAX;
              else if(psdGaussianPrior)prior += -0.5*( (mean-val)*(mean-val)/var + log(2.0*LAL_PI*var) );
            }
          }
          logPrior+=prior;
        }
        else if(!strcmp(item->name,"morlet_Amp") || !strcmp(item->name,"morlet_f0" ) || !strcmp(item->name,"morlet_Q"  ) || !strcmp(item->name,"morlet_t0" ) || !strcmp(item->name,"morlet_phi") )
        {
          REAL8 prior = 0.0;
          gsl_matrix *gparams = *((gsl_matrix **)(item->value));
          UINT4Vector *gsize  = *((UINT4Vector **)LALInferenceGetVariable(params,"glitch_size"));

          char priormin[100];
          char priormax[100];
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
              else prior -= log(component_max-component_min);
            }
          }
          logPrior += prior;//-log(prior);
        }//end morlet parameters prior
        else if(item->type==LALINFERENCE_REAL8_t){
                    sprintf(normName,"%s_norm",item->name);
                    if(LALInferenceCheckVariable(priorParams,normName)) {
                        norm = *(REAL8 *)LALInferenceGetVariable(priorParams,normName);
                    }
                    else
                    {
                        norm = -log(max-min);
                        LALInferenceAddVariable(priorParams, normName, &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
                    }
                    logPrior += norm;
                    //printf("logPrior@%s=%f\n",item->name,logPrior);
                }//end final else of prior

            }//if(val<min || val>max) return -DBL_MAX;

        }//if(item->type!=LALINFERENCE_REAL8_t....
    }//for(;item;item=item->next)

    return(logPrior);
}

/* Convert the hypercube parameter to physical parameters, for the non-spinning inspiral signal case */
UINT4 LALInferenceInspiralPriorNormalisedCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, double *Cube, void *context)
{
    LALInferenceVariableItem *item;
    LALInferenceVariables *priorParams=runState->priorArgs;

    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];
    REAL8 min=-INFINITY, max=INFINITY, logPrior=0.;
    int i = 0;

    // latitude
    if (LALInferenceCheckVariable(params, "declination"))
    {
        item = LALInferenceGetItem(params, "declination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double lat = asin(2.0 * Cube[i] - 1.0);
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
            double longitude = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "rightascension", &longitude);
            Cube[i] = longitude;
            i++;
            strcat(header,"RA ");
        }
    }

    // inclination angle
    if (LALInferenceCheckVariable(params, "inclination"))
    {
        item = LALInferenceGetItem(params, "inclination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double iota;
            iota = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "inclination", &iota);
            Cube[i] = iota;
            logPrior += log(fabs(sin(iota)));
            i++;
            strcat(header,"iota ");
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
            strcat(header,"phi0 ");
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

    // time
    if (LALInferenceCheckVariable(params, "time"))
    {
        item = LALInferenceGetItem(params, "time");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "time", (void *)&min, (void *)&max);
            double tc = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "time", &tc);
            Cube[i] = tc;
            sprintf(timeID,"%d",i);
            i++;
            strcat(header,"time ");
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
    if( LALInferenceCheckVariable(params,"massratio") )
    {
        item = LALInferenceGetItem(params, "massratio");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            eta = *(REAL8 *)LALInferenceGetVariable(params, "massratio");
            if( mc != 0.0 ) LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        }
    }
    else if( LALInferenceCheckVariable(params,"asym_massratio") )
    {
        item = LALInferenceGetItem(params, "asym_massratio");
        if(item->vary == LALINFERENCE_PARAM_FIXED)
        {
            q = *(REAL8 *)LALInferenceGetVariable(params, "asym_massratio");
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
            q = m2/m1; // asymmetric mass ratio, m1 >= m2
            Cube[i] = mc; i++;
            strcat(header,"mchirp ");
            Cube[i] = eta; i++;
            strcat(header,"eta ");
        }

        // chirp mass and eta/q
        if(LALInferenceCheckVariable(params,"massratio")||LALInferenceCheckVariable(params,"asym_massratio"))
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

                  if(LALInferenceCheckVariable(params,"asym_massratio"))
                LALInferenceSetVariable(params, "asym_massratio", &q);
            else if(LALInferenceCheckVariable(params,"massratio"))
                    LALInferenceSetVariable(params, "massratio", &eta);
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
    else if(LALInferenceCheckVariable(params,"spin1"))
    {
        item = LALInferenceGetItem(params, "spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin1", &a_spin1);
            Cube[i] = a_spin1;
            i++;
            strcat(header,"a1 ");
        }
    }

    // theta_spin1
    if(LALInferenceCheckVariable(params,"theta_spin1"))
    {
        item = LALInferenceGetItem(params, "theta_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin1;
            theta_spin1 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin1", &theta_spin1);
            Cube[i] = theta_spin1;
            logPrior += log(fabs(sin(theta_spin1)));
            i++;
            strcat(header,"theta1 ");
        }
    }

    // phi_spin1
    if(LALInferenceCheckVariable(params,"phi_spin1"))
    {
        item = LALInferenceGetItem(params, "phi_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin1", (void *)&min, (void *)&max);
            double phi_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin1", &phi_spin1);
            Cube[i] = phi_spin1;
            i++;
            strcat(header,"phi1 ");
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
    else if(LALInferenceCheckVariable(params,"spin2"))
    {
        item = LALInferenceGetItem(params, "spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin2", &a_spin2);
            Cube[i] = a_spin2;
            i++;
            strcat(header,"a2 ");
        }
    }

    // theta_spin2
    if(LALInferenceCheckVariable(params,"theta_spin2"))
    {
        item = LALInferenceGetItem(params, "theta_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            double theta_spin2;
            theta_spin2 = acos(1.0 - 2.0 * Cube[i]);
            LALInferenceSetVariable(params, "theta_spin2", &theta_spin2);
            Cube[i] = theta_spin2;
            logPrior += log(fabs(sin(theta_spin2)));
            i++;
            strcat(header,"theta2 ");
        }
    }

    // phi_spin2
    if(LALInferenceCheckVariable(params,"phi_spin2"))
    {
        item = LALInferenceGetItem(params, "phi_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin2", (void *)&min, (void *)&max);
            double phi_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin2", &phi_spin2);
            Cube[i] = phi_spin2;
            i++;
            strcat(header,"phi2 ");
        }
    }

    // theta_JN for system-frame parameters
    if(LALInferenceCheckVariable(params,"theta_JN"))
    {
        item = LALInferenceGetItem(params, "theta_JN");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "theta_JN", (void *)&min, (void *)&max);
            double theta_JN = LALInferenceCubeToSinPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "theta_JN", &theta_JN);
            Cube[i] = theta_JN;
            i++;
            strcat(header,"theta_JN ");
        }
    }

    // phi_JL for system-frame parameters
    if(LALInferenceCheckVariable(params,"phi_JL"))
    {
        item = LALInferenceGetItem(params, "phi_JL");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_JL", (void *)&min, (void *)&max);
            double phi_JL = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_JL", &phi_JL);
            Cube[i] = phi_JL;
            i++;
            strcat(header,"phi_JL ");
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
            strcat(header,"tilt_spin1 ");
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
            strcat(header,"tilt_spin2 ");
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
    Cube[i] = LALInferenceInspiralPriorNormalised(runState,params);
    i++; strcat(header,"logprior ");

    // fRef for system-frame parameters
    if(LALInferenceCheckVariable(params,"fRef"))
    {
        Cube[i] = *(REAL8 *)LALInferenceGetVariable(params,"fRef");
        i++;
        strcat(header,"f_ref ");
    }

    strcat(header,"logl");

    /* Check boundaries */
    if (ScaleTest==0) return 0;
    item=params->head;
    for(;item;item=item->next)
    {
        // if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
        if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
        else
        {
            LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
            if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return 0 ;
        }
    }

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

#define SQR(x) ((x)*(x))

static double etaInnerIntegrand(double M2, void *viData) {
  innerData *iData = (innerData *)viData;
  double Mc = pow(M2*iData->M1, 3.0/5.0)/pow(M2+iData->M1, 1.0/5.0);
  double eta = M2*iData->M1/SQR(M2+iData->M1);
  if (Mc < iData->McMin || Mc > iData->McMax || eta < iData->massRatioMin || eta > iData->massRatioMax) {
    return 0.0;
  } else {
    return pow(Mc, -11.0/6.0);
  }
}

#undef SQR

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

#define MIN(x, y) ((x) < (y) ? (x) : (y))

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

    gsl_integration_qag(&f, oData->MMin, MIN(M1, oData->MTotMax-M1), oData->epsabs, oData->epsrel,
                        oData->wsInnerSize, GSL_INTEG_GAUSS61, oData->wsInner, &result, &err);

    return result;
}

#undef MIN

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
    else if(!strcmp(massRatioName,"asym_massratio"))
        oData.innerIntegrand.function = &qInnerIntegrand;
    else if(!strcmp(massRatioName,"massratio"))
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
    //LALInferenceKDDrawFromBox(rdm, tree, proposedPt, Ncell);
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
UINT4 within_malmquist(LALInferenceRunState *runState, LALInferenceVariables *params) {
    UINT4 i=0, nifo=0;

    LALInferenceIFOData *ifo = runState->data;
    while (ifo != NULL) {
        nifo++;
        ifo = ifo->next;
    }

    REAL8 *SNRs = XLALMalloc(nifo * sizeof(REAL8));
    LALInferenceNetworkSNR(params, runState->data, runState->templt, SNRs);
    REAL8 loudest_snr=0.0, second_loudest_snr=0.0, network_snr=0.0;
    for (i=0; i<nifo; i++) {
        if (SNRs[i] > second_loudest_snr) {
            if (SNRs[i] > loudest_snr) {
                second_loudest_snr = loudest_snr;
                loudest_snr = SNRs[i];
            } else {
                second_loudest_snr = SNRs[i];
            }
        }
        network_snr += SNRs[i]*SNRs[i];
    }
    XLALFree(SNRs);
    network_snr = sqrt(network_snr);

    REAL8 malmquist_loudest = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_loudest_snr"));
    REAL8 malmquist_second_loudest = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_second_loudest_snr"));
    REAL8 malmquist_network = (*(REAL8 *)LALInferenceGetVariable(runState->priorArgs,"malmquist_network_snr"));

    if (loudest_snr < malmquist_loudest
          || second_loudest_snr < malmquist_second_loudest
          || network_snr < malmquist_network)
        return(0);
    else
        return(1);
}


REAL8 LALInferenceAnalyticNullPrior(LALInferenceRunState UNUSED *runState, LALInferenceVariables *params) {
  REAL8 logPrior=0.0;
  REAL8 logmc=0.0,mc=0.0;
  REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;

  if(LALInferenceCheckVariable(params,"massratio")||LALInferenceCheckVariable(params,"asym_massratio")) {
    if(LALInferenceCheckVariable(params,"logmc")) {
      logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        LALInferenceMcQ2Masses(exp(logmc),q,&m1,&m2);
        logPrior+=log(m1*m1);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
        LALInferenceMcEta2Masses(exp(logmc),eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
      }
      /*careful using LALInferenceMcEta2Masses, it returns m1>=m2*/
    } else if(LALInferenceCheckVariable(params,"chirpmass")) {
      mc=*(REAL8 *)LALInferenceGetVariable(params,"chirpmass");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        LALInferenceMcQ2Masses(mc,q,&m1,&m2);
        logPrior+=log(m1*m1/mc);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
        LALInferenceMcEta2Masses(mc,eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2))/((m1-m2)*pow(eta,3.0/5.0)));
      }
    }
  }
  logPrior+=LALInferenceFlatBoundedPrior(runState, params);
  return(logPrior);
}

UINT4 LALInferenceAnalyticCubeToPrior(LALInferenceRunState *runState, LALInferenceVariables *params, double *Cube, void *context) {
    int i = 0;
    REAL8 logPrior=0.,min=-INFINITY,max=INFINITY;
    REAL8 m1=1.,m2=1.,mc,eta,m;
    LALInferenceVariableItem *item;

    //char *header=(char *)LALInferenceGetVariable(runState->algorithmParams,"HeaderLine");
    char **info = (char **)context;
    char *header = &info[1][0];
    strcpy(header,"");
    char *timeID = &info[2][0];

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

    if(LALInferenceCheckVariable(params,"inclination"))
    {
        item = LALInferenceGetItem(params,"inclination");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "inclination", (void *)&min, (void *)&max);
            double inclination = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "inclination", &inclination);
            Cube[i] = inclination;
            i++;
            logPrior -= log(max-min);
            strcat(header,"iota ");
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
            Cube[i] = phase;
            i++;
            logPrior -= log(max-min);
            strcat(header,"phi0 ");
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
    else if(LALInferenceCheckVariable(params,"spin1"))
    {
        item = LALInferenceGetItem(params, "spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin1", (void *)&min, (void *)&max);
            double a_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin1", &a_spin1);
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
    else if(LALInferenceCheckVariable(params,"spin2"))
    {
        item = LALInferenceGetItem(params, "spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "spin2", (void *)&min, (void *)&max);
            double a_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "spin2", &a_spin2);
            Cube[i] = a_spin2;
            i++;
            logPrior -= log(max-min);
            strcat(header,"a2 ");
        }
    }

    if(LALInferenceCheckVariable(params,"theta_spin1"))
    {
        item = LALInferenceGetItem(params,"theta_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "theta_spin1", (void *)&min, (void *)&max);
            double theta_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "theta_spin1", &theta_spin1);
            Cube[i] = theta_spin1;
            i++;
            logPrior -= log(max-min);
            strcat(header,"theta1 ");
        }
    }

    if(LALInferenceCheckVariable(params,"theta_spin2"))
    {
        item = LALInferenceGetItem(params,"theta_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "theta_spin2", (void *)&min, (void *)&max);
            double theta_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "theta_spin2", &theta_spin2);
            Cube[i] = theta_spin2;
              i++;
              logPrior -= log(max-min);
            strcat(header,"theta2 ");
        }
    }

    if(LALInferenceCheckVariable(params,"phi_spin1"))
    {
        item = LALInferenceGetItem(params,"phi_spin1");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin1", (void *)&min, (void *)&max);
            double phi_spin1 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin1", &phi_spin1);
            Cube[i] = phi_spin1;
            i++;
            logPrior -= log(max-min);
            strcat(header,"phi1 ");
        }
    }

    if(LALInferenceCheckVariable(params,"phi_spin2"))
    {
        item = LALInferenceGetItem(params,"phi_spin2");
        if(item->vary != LALINFERENCE_PARAM_FIXED)
        {
            LALInferenceGetMinMaxPrior(runState->priorArgs, "phi_spin2", (void *)&min, (void *)&max);
            double phi_spin2 = LALInferenceCubeToFlatPrior(Cube[i], min, max);
            LALInferenceSetVariable(params, "phi_spin2", &phi_spin2);
            Cube[i] = phi_spin2;
            i++;
            logPrior -= log(max-min);
            strcat(header,"phi2 ");
        }
    }

    LALInferenceVariables *priorParams=runState->priorArgs;
    INT4 ScaleTest = LALInferenceCubeToPSDScaleParams(priorParams, params, &i, Cube, context);

    Cube[i] = mc; i++; strcat(header,"mchirp ");
    Cube[i] = eta; i++; strcat(header,"eta ");
    Cube[i] = LALInferenceAnalyticNullPrior(runState,params);
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

REAL8 LALInferenceNullPrior(LALInferenceRunState UNUSED *runState, LALInferenceVariables UNUSED *params) {
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
