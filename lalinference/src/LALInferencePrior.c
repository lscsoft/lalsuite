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
#include <math.h>
#include <gsl/gsl_integration.h>

/* Private helper function prototypes */
static double qInnerIntegrand(double M2, void *viData);
static double etaInnerIntegrand(double M2, void *viData);
static double outerIntegrand(double M1, void *voData);
static double computePriorMassNorm(const double MMin, const double MMax, const double MTotMax, const double McMin, const double McMax, const double massRatioMin, const double massRatioMax, const char *massRatioName);

static void mc2masses(double mc, double eta, double *m1, double *m2);
static void q2masses(double mc, double q, double *m1, double *m2);

static void mc2masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m1 >= m2).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m2 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m1 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
  return;
}

static void q2masses(double mc, double q, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & asymmetric mass   */
/*  ratio (q).  note: q = m2/m1, where m1 >= m2    */
{
  double factor = pow( pow(mc,5.0)/(1+q), 1.0/7.0 );
  *m1 = factor * pow(q, -3.0/7.0);
  *m2 = factor * pow(q, +4.0/7.0);
  return;
}

/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
	REAL8 logPrior=0.0;
	
	(void)runState;
	LALInferenceVariableItem *item=params->head;
	LALInferenceVariables *priorParams=runState->priorArgs;
	REAL8 min, max;
	REAL8 logmc=0.0;
	REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;
	REAL8 tmp=0.;
	/* Check boundaries */
	for(;item;item=item->next)
	{
		// if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
		if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
                        continue;
		else
		{
			LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
		}
	}
	if(LALInferenceCheckVariable(params,"logdistance"))
		logPrior+=3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
	else if(LALInferenceCheckVariable(params,"distance"))
		logPrior+=2.0*log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));
	
	if(LALInferenceCheckVariable(params,"inclination"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))));
	if(LALInferenceCheckVariable(params,"declination"))
		logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
	if(LALInferenceCheckVariable(params,"theta_spin1"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))));
	if(LALInferenceCheckVariable(params,"theta_spin2"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))));
	if(LALInferenceCheckVariable(params,"a_spin1") && LALInferenceCheckVariable(params,"a_spin2")){
		
		if(*(REAL8 *)LALInferenceGetVariable(params,"a_spin2") > *(REAL8 *)LALInferenceGetVariable(params,"a_spin1")){
		 	tmp = *(REAL8 *)LALInferenceGetVariable(params,"a_spin1") ;
			*(REAL8 *)LALInferenceGetVariable(params,"a_spin1") = *(REAL8 *)LALInferenceGetVariable(params,"a_spin2");
			*(REAL8 *)LALInferenceGetVariable(params,"a_spin2") = tmp;
		}
	}
	if(LALInferenceCheckVariable(params,"logmc")) {
          logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
          /* Assume jumping in log(Mc), so use prior that works out to p(Mc) ~ Mc^-11/6 */
          logPrior+=-(5./6.)*logmc;
        } else if(LALInferenceCheckVariable(params,"chirpmass")) {
          logmc=log(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass"));
          /* Assume jumping in Mc, so can implement the Mc^-11/6 directly. */
          logPrior+=-(11./6.)*logmc;
        }
		
	if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        q2masses(exp(logmc),q,&m1,&m2);
    }
	else if(LALInferenceCheckVariable(params,"massratio")) {
		eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
		mc2masses(exp(logmc),eta,&m1,&m2);
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

	return(logPrior);
}

void LALInferenceCyclicReflectiveBound(LALInferenceVariables *parameter,
                                       LALInferenceVariables *priorArgs){
  /* Apply cyclic and reflective boundaries to parameter to bring it back
     within the prior */
  LALInferenceVariableItem *paraHead=NULL;
  REAL8 min,max;
  /* REAL8 mu, sigma; */
  for (paraHead=parameter->head;paraHead;paraHead=paraHead->next) {
    if( paraHead->vary==LALINFERENCE_PARAM_FIXED ||
        paraHead->vary==LALINFERENCE_PARAM_OUTPUT ||
        !LALInferenceCheckMinMaxPrior(priorArgs, paraHead->name) ) continue;
    
    LALInferenceGetMinMaxPrior(priorArgs,paraHead->name, (void *)&min, (void *)&max);
    
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


/** \brief Rotate initial phase if polarisation angle is cyclic around ranges
 * 
 * If the polarisation angle parameter \f$\psi\f$ is cyclic about its upper and
 * lower ranges of \f$-\pi/4\f$ to \f$\psi/4\f$ then the transformation for
 * crossing a boundary requires the initial phase parameter \f$\phi_0\f$ to be
 * rotated through \f$\pi\f$ radians. The function assumes the value of
 * \f$\psi\f$ has been rescaled to be between 0 and \f$2\pi\f$ - this is a
 * requirement of the covariance matrix routine \c LALInferenceNScalcCVM
 * function.  
 * 
 * This is particularly relevant for pulsar analyses.
 * 
 * \param parameter [in] Pointer to an array of parameters
 * \param priorArgs [in] Pointer to an array of prior ranges
 */
void LALInferenceRotateInitialPhase( LALInferenceVariables *parameter){
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
    
      idx1++;
    }
    
    if (paraHead->vary == LALINFERENCE_PARAM_CIRCULAR &&
        !strcmp(paraHead->name, "phi0") ){
      idx1++; 
      idx2++;
    }
    
    if ( idx2 == 0 ) paraPhi0 = paraPhi0->next;
    if ( idx1 == 2 ) break;
  }

  if( rotphi0 != 0. ) *(REAL8 *)paraPhi0->value += rotphi0;
 
  return;
}


/* Return the log Prior of the variables specified for the sky localisation project, ref: https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/SkyLocComparison#priors, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralSkyLocPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
	REAL8 logPrior=0.0;
	
	(void)runState;
	LALInferenceVariableItem *item=params->head;
	LALInferenceVariables *priorParams=runState->priorArgs;
	REAL8 min, max;
	REAL8 logmc=0.0,mc=0.0;
	REAL8 m1=0.0,m2=0.0,q=0.0,eta=0.0;
	/* Check boundaries */
	for(;item;item=item->next)
	{
		// if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR)
		if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT)
      continue;
		else
		{
			LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
		}
	}
  /*Use a uniform in log D distribution*/
	//if(LALInferenceCheckVariable(params,"logdistance"))
	//	logPrior+=3.0* *(REAL8 *)LALInferenceGetVariable(params,"logdistance");
    if(LALInferenceCheckVariable(params,"distance"))
      logPrior-=log(*(REAL8 *)LALInferenceGetVariable(params,"distance"));
	
	if(LALInferenceCheckVariable(params,"inclination"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))));
	if(LALInferenceCheckVariable(params,"declination"))
		logPrior+=log(fabs(cos(*(REAL8 *)LALInferenceGetVariable(params,"declination"))));
	if(LALInferenceCheckVariable(params,"theta_spin1"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))));
	if(LALInferenceCheckVariable(params,"theta_spin2"))
		logPrior+=log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))));

	/*priors uniform in the individual masses. Not taking into account if mtot_max < m1_max+m2_max */
  if(LALInferenceCheckVariable(params,"massratio")||LALInferenceCheckVariable(params,"asym_massratio")) {
    if(LALInferenceCheckVariable(params,"logmc")) {
      logmc=*(REAL8 *)LALInferenceGetVariable(params,"logmc");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        q2masses(exp(logmc),q,&m1,&m2);
        logPrior+=log(m1*m1);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
        mc2masses(exp(logmc),eta,&m1,&m2);
        logPrior+=log(((m1+m2)*(m1+m2)*(m1+m2))/(m1-m2));
      }
      /*careful using mc2masses, it returns m1>=m2*/
    } else if(LALInferenceCheckVariable(params,"chirpmass")) {
      mc=*(REAL8 *)LALInferenceGetVariable(params,"chirpmass");
      if(LALInferenceCheckVariable(params,"asym_massratio")) {
        q=*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio");
        q2masses(mc,q,&m1,&m2);
        logPrior+=log(m1*m1/mc);
      } else {
        eta=*(REAL8 *)LALInferenceGetVariable(params,"massratio");
        mc2masses(mc,eta,&m1,&m2);
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
  
	return(logPrior);
}


/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPriorNormalised(LALInferenceRunState *runState, LALInferenceVariables *params)
{
	REAL8 logPrior=0.0;
	
	(void)runState;
	LALInferenceVariableItem *item=params->head;
	LALInferenceVariables *priorParams=runState->priorArgs;
	REAL8 min, max;
	REAL8 logmc=0.0;
	REAL8 m1,m2; 
	REAL8 massRatioMin=0.0, massRatioMax=0.0; // min,max for q or eta
	REAL8 MTotMax=0.0;
	char normName[VARNAME_MAX];
	char massRatioName[VARNAME_MAX];
	REAL8 norm=0.0;

    if(LALInferenceCheckVariable(params,"asym_massratio")){
        LALInferenceGetMinMaxPrior(priorParams, "asym_massratio", (void *)&massRatioMin, (void *)&massRatioMax);
        strcpy(massRatioName,"asym_massratio");
    }
	else
    {
		LALInferenceGetMinMaxPrior(priorParams, "massratio", (void *)&massRatioMin, (void *)&massRatioMax);
        strcpy(massRatioName,"massratio");
    }
	/* Check boundaries */
	for(;item;item=item->next)
	{
		//if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR) continue;
		if(item->vary==LALINFERENCE_PARAM_FIXED || item->vary==LALINFERENCE_PARAM_OUTPUT) continue;
		else
		{
			LALInferenceGetMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
			else
			{
				if(!strcmp(item->name, "chirpmass") || !strcmp(item->name, "logmc")){
					if(LALInferenceCheckVariable(priorParams,"mass_norm")) {
						norm = *(REAL8 *)LALInferenceGetVariable(priorParams,"mass_norm");
					}
					else
					{
						if( LALInferenceCheckVariable(priorParams,"component_max") && LALInferenceCheckVariable(priorParams,"component_min") 
						   && (LALInferenceCheckVariable(params,"asym_massratio") || LALInferenceCheckVariable(params,"massratio")) ){
							if(LALInferenceCheckVariable(priorParams,"MTotMax")) 
								MTotMax=*(REAL8 *)LALInferenceGetVariable(priorParams,"MTotMax");
							else 
								MTotMax=2.0*(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max"));

							//LALInferenceGetMinMaxPrior(priorParams, "massratio", (void *)&etaMin, (void *)&etaMax); - already done before for loop
							norm = -log(computePriorMassNorm(*(REAL8 *)LALInferenceGetVariable(priorParams,"component_min"),
														*(REAL8 *)LALInferenceGetVariable(priorParams,"component_max"),
														MTotMax, min, max, massRatioMin, massRatioMax, massRatioName));
							//printf("norm@%s=%f\n",item->name,norm);
						}
						else {
							norm = -1.79175946923-log(pow(max,0.166666666667)-pow(min,0.166666666667));
						}
						LALInferenceAddVariable(priorParams, "mass_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					}
					if(!strcmp(item->name, "chirpmass")){
						logmc=log(*(REAL8 *)LALInferenceGetVariable(params,"chirpmass"));
					}
					else if(!strcmp(item->name, "logmc")){
						logmc=(*(REAL8 *)LALInferenceGetVariable(params,"logmc"));
					}
					
                    if(LALInferenceCheckVariable(params,"asym_massratio") || LALInferenceCheckVariable(params,"massratio")){
                        if(LALInferenceCheckVariable(params,"asym_massratio"))
                            q2masses(exp(logmc),*(REAL8 *)LALInferenceGetVariable(params,"asym_massratio"),&m1,&m2);
                        else if(LALInferenceCheckVariable(params,"massratio"))
                            mc2masses(exp(logmc),*(REAL8 *)LALInferenceGetVariable(params,"massratio"),&m1,&m2);
                        
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
					}
					
					logPrior += -(11./6.)*logmc+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "massratio") || !strcmp(item->name, "asym_massratio")) continue;

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
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						LALInferenceAddVariable(priorParams, "inclination_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"inclination"))))+norm;
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
						norm = -sin(LAL_PI*fractpart_min)+sin(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
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
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						LALInferenceAddVariable(priorParams, "theta_spin1_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin1"))))+norm;
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
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						LALInferenceAddVariable(priorParams, "theta_spin2_norm", &norm, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)LALInferenceGetVariable(params,"theta_spin2"))))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				
				else{
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
				}						
										
			}
			
		}
	}
	
	return(logPrior);
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
    char   massRatioName[VARNAME_MAX];
	double MTotMax;
	double MMin;
	double epsabs;
	double epsrel;
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

    if(!strcmp(oData->massRatioName,"asym_massratio"))
        f.function = &qInnerIntegrand;
    else if(!strcmp(oData->massRatioName,"massratio"))
        f.function = &etaInnerIntegrand;
	f.params = &iData;
	
	gsl_integration_qag(&f, oData->MMin, MIN(M1, oData->MTotMax-M1), oData->epsabs, oData->epsrel, 
						oData->wsInnerSize, GSL_INTEG_GAUSS61, oData->wsInner, &result, &err);
	
	return result;
}

#undef MIN

static double computePriorMassNorm(const double MMin, const double MMax, const double MTotMax, 
                    const double McMin, const double McMax,
                    const double massRatioMin, const double massRatioMax, const char *massRatioName) {
	const double epsabs = 1e-8;
	const double epsrel = 1e-8;
	const size_t wsSize = 10000;
	double result, err;
	outerData oData;
	gsl_function f;
	
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
    strcpy(oData.massRatioName,massRatioName);	

	f.function = &outerIntegrand;
	f.params = &oData;
	
	gsl_integration_qag(&f, MMin, MMax, epsabs, epsrel, wsSize, GSL_INTEG_GAUSS61, wsOuter, 
						&result, &err);
	
	gsl_integration_workspace_free(wsOuter);
	gsl_integration_workspace_free(wsInner);
	
	return result;
}


/* Function to add the min and max values for the prior onto the priorArgs */
void LALInferenceAddMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, void *min, void *max, LALInferenceVariableType type){
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
void LALInferenceGetMinMaxPrior(LALInferenceVariables *priorArgs, const char *name, void *min, void *max)
{
		char minName[VARNAME_MAX];
		char maxName[VARNAME_MAX];
		
		sprintf(minName,"%s_min",name);
		sprintf(maxName,"%s_max",name);
    
		*(REAL8 *)min=*(REAL8 *)LALInferenceGetVariable(priorArgs,minName);
		*(REAL8 *)max=*(REAL8 *)LALInferenceGetVariable(priorArgs,maxName);
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
void LALInferenceAddGaussianPrior(LALInferenceVariables *priorArgs, const char *name, void *mu,
  void *sigma, LALInferenceVariableType type){
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
void LALInferenceGetGaussianPrior(LALInferenceVariables *priorArgs, const char *name, void *mu,
  void *sigma)
{
  char meanName[VARNAME_MAX];
  char sigmaName[VARNAME_MAX];
                
  sprintf(meanName,"%s_gaussian_mean",name);
  sprintf(sigmaName,"%s_gaussian_sigma",name);
    
  *(REAL8 *)mu=*(REAL8 *)LALInferenceGetVariable(priorArgs,meanName);
  *(REAL8 *)sigma=*(REAL8 *)LALInferenceGetVariable(priorArgs,sigmaName);
  return;
                
}

void LALInferenceDrawFromPrior( LALInferenceVariables *output, 
                                LALInferenceVariables *priorArgs, 
                                gsl_rng *rdm) {  
  LALInferenceVariableItem *item = output->head;
  
  for(;item;item=item->next){
    if(item->vary==LALINFERENCE_PARAM_CIRCULAR || item->vary==LALINFERENCE_PARAM_LINEAR)
      LALInferenceDrawNameFromPrior( output, priorArgs, item->name, item->type, rdm );
  }
}

void LALInferenceDrawNameFromPrior( LALInferenceVariables *output, 
                                    LALInferenceVariables *priorArgs, 
                                    char *name, LALInferenceVariableType type, 
                                    gsl_rng *rdm) {  
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
    
    LALInferenceGetMinMaxPrior(priorArgs, name, (void *)&min, (void *)&max);
    tmp = min + (max-min)*gsl_rng_uniform( rdm );
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
    default:
      XLALPrintError ("%s: Trying to randomise a non-numeric \
parameter!\n", __func__ );
      XLAL_ERROR_VOID ( XLAL_EFUNC );
      break;
  }
}
