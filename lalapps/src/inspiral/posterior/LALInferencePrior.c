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

#include "LALInferencePrior.h"
#include <math.h>
#include<gsl/gsl_integration.h>

/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPrior(LALInferenceRunState *runState, LALVariables *params)
{
	REAL8 logPrior=0.0;
	
	(void)runState;
	LALVariableItem *item=params->head;
	LALVariables *priorParams=runState->priorArgs;
	REAL8 min, max;
	REAL8 logmc=0.0;
	REAL8 m1,m2,eta=0.0;
	/* Check boundaries */
	for(;item;item=item->next)
	{
		//if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR) continue;
		if(item->vary==PARAM_FIXED || item->vary==PARAM_OUTPUT) continue;
		else
		{
			getMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
		}
	}
	if(checkVariable(params,"logdistance"))
		logPrior+=3.0* *(REAL8 *)getVariable(params,"logdistance");
	else if(checkVariable(params,"distance"))
		logPrior+=2.0*log(*(REAL8 *)getVariable(params,"distance"));
	
	if(checkVariable(params,"inclination"))
		logPrior+=log(fabs(sin(*(REAL8 *)getVariable(params,"inclination"))));
	if(checkVariable(params,"declination"))
		logPrior+=log(fabs(cos(*(REAL8 *)getVariable(params,"declination"))));
	if(checkVariable(params,"theta_spin1"))
		logPrior+=log(fabs(sin(*(REAL8 *)getVariable(params,"theta_spin1"))));
	if(checkVariable(params,"theta_spin2"))
		logPrior+=log(fabs(sin(*(REAL8 *)getVariable(params,"theta_spin2"))));
	
	if(checkVariable(params,"logmc")) {
          logmc=*(REAL8 *)getVariable(params,"logmc");
          /* Assume jumping in log(Mc), so use prior that works out to p(Mc) ~ Mc^-11/6 */
          logPrior+=-(5./6.)*logmc;
        } else if(checkVariable(params,"chirpmass")) {
          logmc=log(*(REAL8 *)getVariable(params,"chirpmass"));
          /* Assume jumping in Mc, so can implement the Mc^-11/6 directly. */
          logPrior+=-(11./6.)*logmc;
        }
		
	if(checkVariable(params,"massratio"))
	{
		eta=*(REAL8 *)getVariable(params,"massratio");
		mc2masses(exp(logmc),eta,&m1,&m2);
	}
	
	/* Check for component masses in range, if specified */
	if(checkVariable(priorParams,"component_min"))
		if(*(REAL8 *)getVariable(priorParams,"component_min") > m1
		   || *(REAL8 *)getVariable(priorParams,"component_min") > m2)
			return -DBL_MAX;
	
	if(checkVariable(priorParams,"component_max"))
		if(*(REAL8 *)getVariable(priorParams,"component_max") < m1
		   || *(REAL8 *)getVariable(priorParams,"component_max") < m2)
			return -DBL_MAX;

	return(logPrior);
}

void LALInferenceCyclicReflectiveBound(LALVariables *parameter, LALVariables *priorArgs){
	/* Apply cyclic and reflective boundaries to parameter to bring it back within
	 the prior */
	LALVariableItem *paraHead=NULL;
	REAL8 delta;
	REAL8 min,max;
	for (paraHead=parameter->head;paraHead;paraHead=paraHead->next)
	{
		if(paraHead->vary==PARAM_FIXED || paraHead->vary==PARAM_OUTPUT) continue;
		getMinMaxPrior(priorArgs,paraHead->name, (void *)&min, (void *)&max);
		if(paraHead->vary==PARAM_CIRCULAR) /* For cyclic boundaries */
		{
			delta = max-min;
			while ( *(REAL8 *)paraHead->value > max) 
				*(REAL8 *)paraHead->value -= delta;
			while ( *(REAL8 *)paraHead->value < min) 
				*(REAL8 *)paraHead->value += delta;
		}
		else if(paraHead->vary==PARAM_LINEAR) /* Use reflective boundaries */
		{
			while(max<*(REAL8 *)paraHead->value || min>*(REAL8 *)paraHead->value){
				/*	printf("%s: max=%lf, min=%lf, val=%lf\n",paraHead->name,max,min,*(REAL8 *)paraHead->value); */
				if(max < *(REAL8 *)paraHead->value) *(REAL8 *)paraHead->value-=2.0*(*(REAL8 *)paraHead->value - max);
				if(min > *(REAL8 *)paraHead->value) *(REAL8 *)paraHead->value+=2.0*(min - *(REAL8 *)paraHead->value);
			}
		}
	}	
	return;
}


/* Return the log Prior of the variables specified, for the non-spinning/spinning inspiral signal case */
REAL8 LALInferenceInspiralPriorNormalised(LALInferenceRunState *runState, LALVariables *params)
{
	REAL8 logPrior=0.0;
	
	(void)runState;
	LALVariableItem *item=params->head;
	LALVariables *priorParams=runState->priorArgs;
	REAL8 min, max;
	REAL8 logmc=0.0;
	REAL8 m1,m2,eta=0.0;
	REAL8 etaMin=0.0, etaMax=0.0;
	REAL8 MTotMax=0.0;
	char normName[VARNAME_MAX];
	REAL8 norm=0.0;
	
	if(checkVariable(params,"massratio")){
		eta=*(REAL8 *)getVariable(params,"massratio");
		getMinMaxPrior(priorParams, "massratio", (void *)&etaMin, (void *)&etaMax);
	}
	
	/* Check boundaries */
	for(;item;item=item->next)
	{
		//if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR) continue;
		if(item->vary==PARAM_FIXED || item->vary==PARAM_OUTPUT) continue;
		else
		{
			getMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
			else
			{
				if(!strcmp(item->name, "chirpmass") || !strcmp(item->name, "logmc")){
					if(checkVariable(priorParams,"mass_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"mass_norm");
					}
					else
					{
						if(checkVariable(priorParams,"component_max") && checkVariable(priorParams,"component_min") 
						   && checkVariable(params,"massratio")){
							if(checkVariable(priorParams,"MTotMax")) 
								MTotMax=*(REAL8 *)getVariable(priorParams,"MTotMax");
							else 
								MTotMax=2.0*(*(REAL8 *)getVariable(priorParams,"component_max"));

							getMinMaxPrior(priorParams, "massratio", (void *)&etaMin, (void *)&etaMax);
							norm = -log(computePriorMassNorm(*(REAL8 *)getVariable(priorParams,"component_min"),
														*(REAL8 *)getVariable(priorParams,"component_max"),
														MTotMax, min, max, etaMin, etaMax));
							//printf("norm@%s=%f\n",item->name,norm);
						}
						else {
							norm = -1.79175946923-log(pow(max,0.166666666667)-pow(min,0.166666666667));
						}
						addVariable(priorParams, "mass_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					if(!strcmp(item->name, "chirpmass")){
						logmc=log(*(REAL8 *)getVariable(params,"chirpmass"));
					}
					else if(!strcmp(item->name, "logmc")){
						logmc=(*(REAL8 *)getVariable(params,"logmc"));
					}
					
					if(checkVariable(params,"massratio")){
						mc2masses(exp(logmc),*(REAL8 *)getVariable(params,"massratio"),&m1,&m2);
					
						if(checkVariable(priorParams,"component_min"))
							if(*(REAL8 *)getVariable(priorParams,"component_min") > m1
							   || *(REAL8 *)getVariable(priorParams,"component_min") > m2)
								return -DBL_MAX;
					
						if(checkVariable(priorParams,"component_max"))
							if(*(REAL8 *)getVariable(priorParams,"component_max") < m1
							   || *(REAL8 *)getVariable(priorParams,"component_max") < m2)
								return -DBL_MAX;

						if(checkVariable(priorParams,"MTotMax"))
							if(*(REAL8 *)getVariable(priorParams,"MTotMax") < m1+m2)
								return -DBL_MAX;
					}
					
					logPrior += -(11./6.)*logmc+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "massratio")) continue;

				else if(!strcmp(item->name, "distance")){
					if(checkVariable(priorParams,"distance_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"distance_norm");
					}
					else
					{
						norm = +1.09861228867-log(max*max*max-min*min*min);
						addVariable(priorParams, "distance_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += 2.0*log(*(REAL8 *)getVariable(params,"distance"))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "logdistance")){
					if(checkVariable(priorParams,"logdistance_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"logdistance_norm");
					}
					else
					{
						norm = 1.38629436112-log(max*max*max*max-min*min*min*min);
						addVariable(priorParams, "logdistance_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += 3.0* *(REAL8 *)getVariable(params,"logdistance")+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "inclination")){
					if(checkVariable(priorParams,"inclination_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"inclination_norm");
					}
					else
					{
						REAL8 intpart_min=0.0;
						REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
						REAL8 intpart_max=0.0;
						REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						addVariable(priorParams, "inclination_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)getVariable(params,"inclination"))))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "declination")){
					if(checkVariable(priorParams,"declination_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"declination_norm");
					}
					else
					{
						REAL8 intpart_min=0.0;
						REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
						REAL8 intpart_max=0.0;
						REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
						norm = -sin(LAL_PI*fractpart_min)+sin(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						addVariable(priorParams, "declination_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += log(fabs(cos(*(REAL8 *)getVariable(params,"declination"))))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "theta_spin1")){
					if(checkVariable(priorParams,"theta_spin1_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"theta_spin1_norm");
					}
					else
					{
						REAL8 intpart_min=0.0;
						REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
						REAL8 intpart_max=0.0;
						REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						addVariable(priorParams, "theta_spin1_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)getVariable(params,"theta_spin1"))))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				else if(!strcmp(item->name, "theta_spin2")){
					if(checkVariable(priorParams,"theta_spin2_norm")) {
						norm = *(REAL8 *)getVariable(priorParams,"theta_spin2_norm");
					}
					else
					{
						REAL8 intpart_min=0.0;
						REAL8 fractpart_min = modf(min/LAL_PI , &intpart_min);
						REAL8 intpart_max=0.0;
						REAL8 fractpart_max = modf(max/LAL_PI , &intpart_max);
						norm = cos(LAL_PI*fractpart_min)-cos(LAL_PI*fractpart_max)+2.0*(intpart_max-intpart_min);
						addVariable(priorParams, "theta_spin2_norm", &norm, REAL8_t, PARAM_FIXED);
					}
					logPrior += log(fabs(sin(*(REAL8 *)getVariable(params,"theta_spin2"))))+norm;
					//printf("logPrior@%s=%f\n",item->name,logPrior);
				}
				
				else{
					sprintf(normName,"%s_norm",item->name);
					if(checkVariable(priorParams,normName)) {
						norm = *(REAL8 *)getVariable(priorParams,normName);
					}
					else
					{
						norm = -log(max-min);
						addVariable(priorParams, normName, &norm, REAL8_t, PARAM_FIXED);
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
	double etaMin;
	double etaMax;
} innerData;

#define SQR(x) ((x)*(x))

double innerIntegrand(double M2, void *viData) {
	innerData *iData = (innerData *)viData;
	double Mc = pow(M2*iData->M1, 3.0/5.0)/pow(M2+iData->M1, 1.0/5.0);
	double eta = M2*iData->M1/SQR(M2+iData->M1);
	if (Mc < iData->McMin || Mc > iData->McMax || eta < iData->etaMin || eta > iData->etaMax) {
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
	double etaMin;
	double etaMax;
	double MTotMax;
	double MMin;
	double epsabs;
	double epsrel;
} outerData;

#define MIN(x, y) ((x) < (y) ? (x) : (y))

double outerIntegrand(double M1, void *voData) {
	outerData *oData = (outerData *)voData;
	gsl_function f;
	innerData iData;
	double result, err;
	
	iData.M1 = M1;
	iData.McMin = oData->McMin;
	iData.McMax = oData->McMax;
	iData.etaMin = oData->etaMin;
	iData.etaMax = oData->etaMax;
	
	f.function = &innerIntegrand;
	f.params = &iData;
	
	gsl_integration_qag(&f, oData->MMin, MIN(M1, oData->MTotMax-M1), oData->epsabs, oData->epsrel, 
						oData->wsInnerSize, GSL_INTEG_GAUSS61, oData->wsInner, &result, &err);
	
	return result;
}

#undef MIN

double computePriorMassNorm(const double MMin, const double MMax, const double MTotMax, 
                    const double McMin, const double McMax,
                    const double etaMin, const double etaMax) {
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
	oData.etaMin = etaMin;
	oData.etaMax = etaMax;
	oData.MTotMax = MTotMax;
	oData.epsabs = epsabs;
	oData.epsrel = epsrel;
	oData.MMin = MMin;
	
	f.function = &outerIntegrand;
	f.params = &oData;
	
	gsl_integration_qag(&f, MMin, MMax, epsabs, epsrel, wsSize, GSL_INTEG_GAUSS61, wsOuter, 
						&result, &err);
	
	gsl_integration_workspace_free(wsOuter);
	gsl_integration_workspace_free(wsInner);
	
	return result;
}


