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

/* Return the log Prior of the variables specified, for the non-spinning inspiral signal case */
REAL8 LALInferenceInspiralPriorNonSpinning(LALInferenceRunState *runState, LALVariables *params)
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
	
	if(checkVariable(params,"logmc"))
		logmc=*(REAL8 *)getVariable(params,"logmc");
	else if(checkVariable(params,"chirpmass"))
		logmc=log(*(REAL8 *)getVariable(params,"chirpmass"));
	
	logPrior+=-(5./6.)*logmc;
	
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
