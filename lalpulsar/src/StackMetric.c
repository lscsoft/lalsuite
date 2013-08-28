/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/AVFactories.h>
#include<lal/StackMetric.h>

/**
 * \brief Computes the parameter space metric for a coherent pulsar search.
 * \author Creighton, T. D.
 * \date 2000, 2001
 * \ingroup StackMetric_h
 *
 * \heading{Description}
 *
 * This function computes the metric \f$g_{\alpha\beta}(\mathbf{\lambda})\f$, as
 * discussed in \ref StackMetric_h, under the assumption
 * that the detected power is constructed from the incoherent sum of \f$N\f$
 * separate power spectrum, each derived from separate time intervals of
 * length \f$\Delta t\f$.  The indecies \f$\alpha\f$ and \f$\beta\f$ are assumed to
 * run from 0 to \f$n\f$, where \f$n\f$ is the total number of ``shape''
 * parameters.
 *
 * This routine has exactly the same calling structure and data storage
 * as the LALCoherentMetric() function.  Thus, the argument
 * \a *metric is a vector of length \f$(n+1)(n+2)/2\f$ storing all
 * non-redundant coefficients of \f$g_{\alpha\beta}\f$, or twice this length
 * if \a params->errors is nonzero.  See LALCoherentMetric() for
 * the indexing scheme.  The argument \a lambda is another vector,
 * of length \f$n+1\f$, storing the components of
 * \f$\mathbf{\lambda}=(\lambda^0,\ldots,\lambda^n)\f$ for the parameter space
 * point at which the metric is being evaluated.  The argument
 * \a *params stores the remaining parameters for computing the
 * metric, as given in the Structures section of \ref StackMetric_h.
 *
 * \heading{Algorithm}
 *
 * Most of what this routine does is set up arguments to be passed to the
 * function LALCoherentMetric().  Each metric component in the stack
 * metric is given simply by:
 * \f[
 * g_{\alpha\beta}(\mathbf{\lambda}) = \frac{1}{N} \sum_{k=1}^N
 * g^{(k)}_{\alpha\beta}(\mathbf{\lambda}) \; ,
 * \f]
 * where \f$g^{(k)}_{\alpha\beta}\f$ is just the coherent metric computed on
 * the time interval \f$[t_\mathrm{start}+(k-1)\Delta t,
 * t_\mathrm{start}+k\Delta t]\f$.  The estimated uncertainty
 * \f$s_{\alpha\beta}\f$ in this component is taken to be:
 * \f[
 * s_{\alpha\beta} = \frac{1}{\sqrt{N}} \max_{k\in\{1,\ldots,N\}}
 * s^{(k)}_{\alpha\beta} \; .
 * \f]
 * There are no clever tricks involved in any of these computations.
 */
void
LALStackMetric( LALStatus        *stat,
		REAL8Vector      *metric,
		REAL8Vector      *lambda,
		MetricParamStruc *params )
{
  INT4 n;  /* Number of metric coefficients. */
  INT4 i;  /* An index. */
  INT4 j;  /* Another index. */
  REAL8 t; /* Time at start of each coherent stretch. */
  REAL8Vector *subMetric=NULL; /* Coherent metric on each stack. */
  MetricParamStruc subParams;  /* Parameters passed to coherent metric
				  computation. */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT(metric,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(metric->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);

  /* Make sure that metric length is positive. */
  n=metric->length;
  ASSERT(n>0,stat,STACKMETRICH_EBAD,STACKMETRICH_MSGEBAD);

  /* Set up parameters for coherent metric computation. */
  memset(metric->data,0,n*sizeof(REAL8));
  memcpy(&subParams,params,sizeof(MetricParamStruc));
  subParams.n=1;
  TRY(LALDCreateVector(stat->statusPtr,&subMetric,n),stat);
  t=params->start;

  /* Compute coherent metrics and accumulate them. */
  i=params->n;
  while(i--){
    subParams.start=t;
    LALCoherentMetric(stat->statusPtr,subMetric,lambda,&subParams);
    BEGINFAIL(stat)
      TRY(LALDDestroyVector(stat->statusPtr,&subMetric),stat);
    ENDFAIL(stat);
    if(params->errors)
      for(j=0;j<n;j+=2){
	metric->data[j]+=subMetric->data[j];
	if(metric->data[j+1]<subMetric->data[j+1])
	  metric->data[j+1]=subMetric->data[j+1];
      }
    else
      for(j=0;j<n;j++)
	metric->data[j]+=subMetric->data[j];
    t+=params->deltaT;
  }

  /* Compute the average. */
  if(params->errors)
    for(j=0;j<n;j+=2){
      metric->data[j]/=params->n;
      metric->data[j+1]/=sqrt(params->n);
    }
  else
    for(j=0;j<n;j++)
      metric->data[j]/=params->n;

  /* Cleanup and exit. */
  TRY(LALDDestroyVector(stat->statusPtr,&subMetric),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
