
/*
 *  LALInferencePriorVolumes.c
 *
 *  Copyright (C) 2016 Salvatore Vitale
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
#include <math.h>
#include <gsl/gsl_integration.h>
#include <lal/LALInference.h>
#include <lal/LALInferencePriorVolumes.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALCosmologyCalculator.h>

typedef struct {
  double mc;
  LALInferenceVariables *prior;
} innerParams;

typedef struct {
  LALInferenceVariables *prior;
} outerParams;


static double chirp_to_comp_jacobian(double mc,double q);
static double mass_indicator(double mc, double q,LALInferenceVariables *priorParams);
static double integrand(double q,void *params);
static double inner_integral(double mc, void *params);
static double mass_outer_integral(LALInferenceVariables *priorArgs);
static double loudness_volume(LALInferenceRunState *state);

static double chirp_to_comp_jacobian(double mc,double q){

  double factor = mc * pow(1.0 +q, 0.2);
  double  m1 = factor * pow(q, -0.6);
  return m1*m1/mc;

}

static double mass_indicator(double mc, double q,LALInferenceVariables *priorParams){
  /*
   * 
   * This function gets mchirp and q
   * Calculates component and total mass and check if all parameters
   * are within their priors.
   * 
   * Return 1 or 0
   * 
   * TODO: make it work with eta
   * 
   * */
  double factor = mc * pow(1.0 +q, 0.2);
  double m1 = factor * pow(q, -0.6);
  double m2 = factor * pow(q, 0.4);

  /* Check for individual mass priors */
  if(LALInferenceCheckVariable(priorParams,"mass1_min"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass1_min") > m1)
      return 0.0;
  if(LALInferenceCheckVariable(priorParams,"mass1_max"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass1_max") < m1)
      return 0.0;
  if(LALInferenceCheckVariable(priorParams,"mass2_min"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass2_min") > m2)
      return 0.0;
  if(LALInferenceCheckVariable(priorParams,"mass2_max"))
    if(LALInferenceGetREAL8Variable(priorParams,"mass2_max") < m2)
      return 0.0;
  if(LALInferenceCheckVariable(priorParams,"MTotMax"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax") < m1+m2)
      return 0.0;
  if(LALInferenceCheckVariable(priorParams,"MTotMin"))
    if(*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMin") > m1+m2)
      return 0.0;
  if (LALInferenceCheckVariable(priorParams,"q_min")){
    double q_min,q_max;
    LALInferenceGetMinMaxPrior(priorParams, "q", &q_min, &q_max);
    if (q<q_min || q>q_max) return 0.0;
  }
  if (LALInferenceCheckVariable(priorParams,"chirpmass_min")){
    double mc_min,mc_max;
    LALInferenceGetMinMaxPrior(priorParams, "chirpmass", &mc_min, &mc_max);
    if (mc<mc_min || mc>mc_max) return 0.0;
  }
  return 1.0;
  /*
  if (LALInferenceCheckVariable(priorParams,"eta_min")){
    double eta_min,eta_max;
    LALInferenceGetMinMaxPrior(state->priorArgs, "eta", &eta_min, &eta_max);
    if (eta<eta_min || eta>eta_max) return 0.0;
  }*/
}

static double integrand(double q,void *params){
  
  /* Integrand for the dobule integral over Mc and q
   * 
   * This is the jacobian within the support of the prior, zero otherwise
   * 
   * */
  innerParams iData = *(innerParams *)params;
  double mc= iData.mc;
  
  return chirp_to_comp_jacobian(mc,q)*mass_indicator(mc,q, iData.prior);
  
}

static double inner_integral(double mc, void *params){
  
  // q comes through params
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;
  gsl_function F;
  F.function = &integrand;
  innerParams iParams;
  F.params = &iParams;
  outerParams oParams=*(outerParams *) params;
  iParams.mc=mc;
  iParams.prior= (LALInferenceVariables *) oParams.prior;
  
  double q_min,q_max;
  if (LALInferenceCheckVariable(iParams.prior,"q_min")){
    
    LALInferenceGetMinMaxPrior(iParams.prior, "q", &q_min, &q_max);
  }
  else{
    fprintf(stderr,"ERROR: q doesn't seem to be a valid param. Exiting\n");
    exit(1);
    }
  // TODO: make it work with eta
  int status = gsl_integration_qags (&F, q_min, q_max,epsabs, epsrel, wsSize,
                        w, &result, &error); 
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

  gsl_integration_workspace_free(w);
  return result;  
  }

  
static double mass_outer_integral(LALInferenceVariables *priorArgs){
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;

  gsl_function F;
  F.function = &inner_integral;
  outerParams oParams;
  F.params = &oParams;
  oParams.prior=priorArgs;
  
  double mc_min,mc_max;
  if (LALInferenceCheckVariable(priorArgs,"chirpmass_min")){
    
    LALInferenceGetMinMaxPrior(priorArgs, "chirpmass", &mc_min, &mc_max);
  }
  else{
    fprintf(stderr,"ERROR: chirpmass doesn't seem to be a valid param. Exiting\n");
    exit(1);
    }
  /* this integrates on chirpmass */
  int status = gsl_integration_qags (&F, mc_min, mc_max, epsabs, epsrel, wsSize, w, &result, &error); 
  
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

  gsl_integration_workspace_free(w);
  return result;
  }

typedef REAL8 (*LALInferenceLoudnessPriorFunction) (double x,LALCosmologicalParameters *omega);
  
typedef struct {
  LALInferenceVariables *priorargs;
  LALInferenceLoudnessPriorFunction priorfunc;
  LALCosmologicalParameters *omega;
} loudnessParams;

static double distance_prior(double d,LALCosmologicalParameters *omega){
  (void) omega;
  return d*d;
}
static double logdistance_prior(double ld,LALCosmologicalParameters *omega){
  (void) omega;
  return ld*ld*ld;
}
static double redshift_prior(double z,LALCosmologicalParameters *omega){
  return XLALUniformComovingVolumeDensity(z,omega);
}
  
static double loudness_integrand(double x,void *params){
  
  loudnessParams lParams=*(loudnessParams *) params;
  LALInferenceLoudnessPriorFunction priorf= lParams.priorfunc;
  LALCosmologicalParameters *omega=lParams.omega;
  return priorf(x,omega);
}
  
static double loudness_volume(LALInferenceRunState *state){


  LALInferenceVariables *priorArgs=state->priorArgs;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
  double result, error;
  const double epsabs = 1e-4;
  const double epsrel = 1e-4;
  const size_t wsSize = 10000;

  gsl_function F;
  F.function = &loudness_integrand;
  loudnessParams lParams;
  F.params = &lParams;
  double intmin,intmax;
  if (LALInferenceCheckVariable(priorArgs,"redshift_min")){
    lParams.priorfunc=&redshift_prior;
    // WILL need to enable this after we push the cosmology patch, for the moment set to null
    //lParams.omega=state->omega;
    lParams.omega=NULL;

    LALInferenceGetMinMaxPrior(priorArgs, "redshift", &intmin, &intmax);
  }
  else if (LALInferenceCheckVariable(priorArgs,"distance_min")){
    lParams.priorfunc=&distance_prior;
    lParams.omega=NULL;
    LALInferenceGetMinMaxPrior(priorArgs, "distance", &intmin, &intmax);
  }
  else if (LALInferenceCheckVariable(priorArgs,"logdistance_min")){
    lParams.priorfunc=&logdistance_prior;
    lParams.omega=NULL;
    LALInferenceGetMinMaxPrior(priorArgs, "logdistance", &intmin, &intmax);
  }
  else XLAL_ERROR_REAL8(XLAL_EINVAL,"No known distance parameter found, unable to proceed\n");

  lParams.priorargs=priorArgs;
  int status = gsl_integration_qags (&F, intmin, intmax, epsabs, epsrel, wsSize, w, &result, &error); 
  
  if (status)
        XLAL_ERROR_REAL8(XLAL_EFUNC | XLAL_EDATA, "Bad data; GSL integration failed.");

  gsl_integration_workspace_free(w);
  return result;
  
}
  
double LALInferenceMassPriorVolume(LALInferenceRunState *state){

  LALInferenceVariables *priorArgs=state->priorArgs;
  return mass_outer_integral(priorArgs);
}

double LALInferenceMassDistancePriorVolume(LALInferenceRunState *state){
  
  return LALInferenceMassPriorVolume(state)*loudness_volume(state);
}
