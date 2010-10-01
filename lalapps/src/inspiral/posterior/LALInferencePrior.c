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
	/* Check boundaries */
	for(;item;item=item->next)
	{
		if(item->vary!=PARAM_LINEAR || item->vary!=PARAM_CIRCULAR) continue;
		else
		{
			getMinMaxPrior(priorParams, item->name, (void *)&min, (void *)&max);
			if(*(REAL8 *) item->value < min || *(REAL8 *)item->value > max) return -DBL_MAX;
		}
	}
	
	if(checkVariable(params,"distance"))
		logPrior+=2.0*log(*(REAL8 *)getVariable(params,"distance"));
	else if(checkVariable(params,"logdistance"))
		logPrior+=3.0* *(REAL8 *)getVariable(params,"logdistance");
	
	if(checkVariable(params,"inclination"))
		logPrior+=log(fabs(*(REAL8 *)getVariable(params,"inclination")));
	if(checkVariable(params,"declination"))
		logPrior+=log(fabs(*(REAL8 *)getVariable(params,"declination")));
	
	
	return(logPrior);
}
