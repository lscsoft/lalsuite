/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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
#include<lal/LALConstants.h>
#include<lal/AVFactories.h>
#include<lal/SeqFactories.h>
#include<lal/StackMetric.h>

#define COHERENTMETRICC_NPTS 100000	/** Number of points to average per time interval. */


/* #define APPROXIMATE Whether to use the approximation in the last
		       form of the averaging integral */

/** Local function to perform the averaging. */
static REAL8 Average(REAL8Vector *integrand, REAL8 *uncertainty);


/**
 * \brief Computes the parameter space metric for a coherent pulsar search.
 * \author Creighton, T. D., Jolien Creighton
 * \date 2000 - 2003
 * \ingroup StackMetric_h
 *
 * This function computes the metric \f$g_{\alpha\beta}(\mathbf{\lambda})\f$, as
 * discussed in \ref StackMetric_h, under the assumption
 * that the search consists of scanning the Fourier power spectrum
 * constructed from a single time interval \f$\Delta t\f$.  The indecies
 * \f$\alpha\f$ and \f$\beta\f$ are assumed to run from 0 to \f$n\f$, where \f$n\f$ is
 * the total number of shape parameters.
 * The argument \c metric is normally a vector of length
 * \f$(n+1)(n+2)/2\f$ storing all non-redundant coefficients of
 * \f$g_{\alpha\beta}\f$.  The indexing scheme is as follows: Let us assume
 * that \f$\alpha\geq\beta\f$.  Then
 * \f$g_{\alpha\beta} = g_{\beta\alpha} = \f$metric->data[\f$\beta + \alpha(\alpha+1)/2]\f$.
 * If <tt>params->errors</tt> is nonzero, then \a *metric must be double
 * this length, and LALCoherentMetric() will store metric
 * components and their estimated uncertainty in alternate slots; i.e.
 * the metric component
 * \f$g_{\alpha\beta}=\f$metric->data[\f$2\beta+\alpha(\alpha+1)]\f$,
 * and the uncertainty
 * \f$s_{\alpha\beta}=\f$metric->data[\f$1+2\beta+\alpha(\alpha+1)]\f$.
 *
 * The argument \a lambda is another vector, of length \f$n+1\f$,
 * storing the components of \f$\mathbf{\lambda}=(\lambda^0,\ldots,\lambda^n)\f$
 * for the parameter space point at which the metric is being evaluated.
 * The argument \a *params stores the remaining parameters for
 * computing the metric, as given in the Structures section of StackMetric.h.
 *
 * \heading{Algorithm}
 *
 * This routine simply computes the function given in Eq.\eqref{eq_gab_phi}
 * of \ref StackMetric_h.  Most of the work is done by the
 * function \a params->dtCanon(), which computes the value and
 * parameter derivatives of the canonical time coordinate
 * \f$\tau[t;\vec{\lambda}]\f$.  Since the phase function is simply
 * \f$\phi[t;\mathbf{\lambda}]=2\pi\lambda^0 \tau[t;\vec\lambda]\f$, where
 * \f$\lambda^0=f_0\f$, the metric components are simply:
 * \f{eqnarray}
 * g_{00}(\mathbf{\lambda}) & = &
 * 4\pi^2\bigg\langle\!(\tau[t;\vec\lambda])^2\bigg\rangle -
 * 4\pi^2\bigg\langle\!\tau[t;\vec\lambda]\bigg\rangle^2 \; ,
 * \nonumber\\
 * g_{i0}(\mathbf{\lambda}) \;\; = \;\;
 * g_{0i}(\mathbf{\lambda}) & = &
 * 4\pi^2 f_0\bigg\langle\!
 * \tau[t;\vec\lambda]
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
 * \bigg\rangle
 * -
 * 4\pi^2 f_0\bigg\langle\!
 * \tau[t;\vec\lambda]
 * \bigg\rangle
 * \bigg\langle
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
 * \bigg\rangle \; , \nonumber\\
 * g_{ij}(\mathbf{\lambda}) & = &
 * 4\pi^2 f_0^2\bigg\langle
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^j}
 * \bigg\rangle
 * -
 * 4\pi^2 f_0^2\bigg\langle
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
 * \bigg\rangle
 * \bigg\langle
 * \frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^j}
 * \bigg\rangle \; , \nonumber
 * \f}
 * where the indecies \f$i\f$ and \f$j\f$ run from 1 to \f$n\f$.
 *
 * In the rigorous definition of the metric, the angle brackets denote an
 * average over the \em canonical time, not the detector time, so we have:
 * \f{eqnarray}
 * \bigg\langle\ldots\bigg\rangle
 * & = & \frac{1}{\tau(t_\mathrm{start}+\Delta t)-\tau(t_\mathrm{start})}
 * \int_{\tau(t_\mathrm{start})}^{\tau(t_\mathrm{start}+\Delta t)}
 * \ldots\,d\tau \nonumber\\
 * & = & \frac{1}{\tau(t_\mathrm{start}+\Delta t)-\tau(t_\mathrm{start})}
 * \int_{t_\mathrm{start}}^{t_\mathrm{start}+\Delta t}
 * \ldots\frac{\partial\tau}{\partial t}\,dt \nonumber\\
 * & \approx & \frac{1}{\Delta t}
 * \int_{t_\mathrm{start}}^{t_\mathrm{start}+\Delta t}
 * \ldots\,dt \nonumber \; ,
 * \f}
 * where the approximation is good to order \f$\epsilon\f$ equal to the
 * maximum difference between 1 and \f$\partial\tau/\partial t\f$.  For an
 * Earth-motion barycentred time coordinate, for instance, \f$\epsilon\f$ is
 * of order \f$10^{-4}\f$ and the approximation is quite good.  However, the
 * current implementation of the metric algorithm uses the second, exact
 * formula.  If speed considerations become important, this is one
 * obvious area for simplification.
 *
 * At present, the time averaging is performed by evaluating \f$\tau\f$ and
 * its derivatives at equally-spaced points over \f$\Delta t\f$ and applying
 * a trapezoidal method; i.e.
 * \f[
 * \frac{1}{T}\int_0^T F(t)\,dt \approx \frac{1}{N}\left(\frac{1}{2}F_0
 * + F_1 + F_2 + \ldots + F_{N-1} + \frac{1}{2}F_N \right) \; ,
 * \f]
 * where \f$F_k=F(kT/N)\f$ are the \f$N+1\f$ sample points of the integrand.  The
 * number of points is a compiled-in constant.  The error is on the order
 * of \f$F''T^2/N^2\f$, where \f$F''\f$ is the second derivative of \f$F(t)\f$ at
 * some point in the interval; we liberally estimate the error to be:
 * \f[
 * \sigma = \max_{k\in\{1,\ldots,N-1\}}
 * \left\{F_{k-1} - 2F_k + F_{k+1}\right\} \; .
 * \f]
 * Other more sophisticated integration techniques may be considered in
 * future.  To save on recomputing costs, the derivatives of \f$\tau\f$ at
 * all times are stored in a large vector sequence.
 */
void
LALCoherentMetric( LALStatus        *stat,
		   REAL8Vector      *metric,
		   REAL8Vector      *lambda,
		   MetricParamStruc *params )
{
  UINT4 s = 0;      /* The number of parameters. */
  UINT4 i, j, k, l; /* Indecies. */
  REAL8 dt;         /* The sampling interval, in s. */
  REAL8 f0;         /* The frequency scale, lambda->data[0]. */
  REAL8 *data;      /* Multipurpose pointer to vector data. */
  REAL8Vector *a=NULL;  /* Data for the first time average. */
  REAL8Vector *b=NULL;  /* Data for the second time average. */
  REAL8Vector *c=NULL;  /* Data for the third time average. */
  REAL8Vector d;        /* Temporary variable. */
  REAL8Vector *variables=NULL;    /* Input to time function. */
  REAL8VectorSequence *dSeq=NULL; /* Phase derivatives. */
  CreateVectorSequenceIn in; /* Input structure. */
  PulsarTimesParamStruc *constants; /* Timing constants. */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR(stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT(metric,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(metric->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(lambda,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(lambda->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(params,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(params->dtCanon,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);

  /* Make sure parameters are valid. */
  s=lambda->length;
  ASSERT(s,stat,STACKMETRICH_EBAD,STACKMETRICH_MSGEBAD);
  if(params->errors)
  {
    ASSERT(metric->length==s*(s+1),stat,STACKMETRICH_EBAD,
	   STACKMETRICH_MSGEBAD);
  }
  else
  {
    ASSERT(metric->length==(s*(s+1))>>1,stat,STACKMETRICH_EBAD,
	   STACKMETRICH_MSGEBAD);
  }
  ASSERT(params->n==1,stat,STACKMETRICH_EBAD,
	 STACKMETRICH_MSGEBAD);

  /* Set vector sequence creation structure, and compute sampling
     rate. */
  in.vectorLength=s+1;
  in.length=COHERENTMETRICC_NPTS;
  dt=params->deltaT/(COHERENTMETRICC_NPTS-1.0);

  /* Create temporary storage. */
  TRY(LALDCreateVector(stat->statusPtr,&a,COHERENTMETRICC_NPTS),stat);
  LALDCreateVector(stat->statusPtr,&b,COHERENTMETRICC_NPTS);
  BEGINFAIL(stat)
    TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
  ENDFAIL(stat);
  LALDCreateVector(stat->statusPtr,&c,COHERENTMETRICC_NPTS);
  BEGINFAIL(stat) {
    TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
    TRY(LALDDestroyVector(stat->statusPtr,&b),stat);
  } ENDFAIL(stat);
  LALDCreateVectorSequence(stat->statusPtr,&dSeq,&in);
  BEGINFAIL(stat) {
    TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
    TRY(LALDDestroyVector(stat->statusPtr,&b),stat);
    TRY(LALDDestroyVector(stat->statusPtr,&c),stat);
  } ENDFAIL(stat);
  LALDCreateVector(stat->statusPtr,&variables,s);
  BEGINFAIL(stat) {
    TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
    TRY(LALDDestroyVector(stat->statusPtr,&b),stat);
    TRY(LALDDestroyVector(stat->statusPtr,&c),stat);
    TRY(LALDDestroyVectorSequence(stat->statusPtr,&dSeq),stat);
  } ENDFAIL(stat);

  /* Set up passing structure for the canonical time function. */
  memcpy(variables->data,lambda->data,s*sizeof(REAL8));
  *(variables->data)=params->start;
  f0=lambda->data[0];
  constants=params->constants;
  d.length=s+1;
  d.data=dSeq->data;
  l=COHERENTMETRICC_NPTS;

  /* Compute canonical time and its derivatives. */
  while(l--){
    (params->dtCanon)(stat->statusPtr,&d,variables,constants);
    BEGINFAIL(stat) {
      TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
      TRY(LALDDestroyVector(stat->statusPtr,&b),stat);
      TRY(LALDDestroyVector(stat->statusPtr,&c),stat);
      TRY(LALDDestroyVectorSequence(stat->statusPtr,&dSeq),stat);
      TRY(LALDDestroyVector(stat->statusPtr,&variables),stat);
    } ENDFAIL(stat);
    *(variables->data)+=dt;
    d.data+=s+1;
  }

  /* Compute the overall correction factor to convert the canonical
     time interval into the detector time interval. */
#ifndef APPROXIMATE
  dt=params->deltaT/(dSeq->data[(s+1)*(COHERENTMETRICC_NPTS-1)]
		     - dSeq->data[0]);
#else
  dt=1.0;
#endif

  /* Compute metric components.  First generate the g00 component. */
  /* Fill integrand arrays. */
  data=dSeq->data+1;
  l=COHERENTMETRICC_NPTS;
  k=0;
  while(l--){
#ifndef APPROXIMATE
    a->data[l]=*data*data[-1]*data[-1];
    b->data[l]=c->data[l]=*data*data[-1];
#else
    a->data[l]=data[-1]*data[-1];
    b->data[l]=c->data[l]=data[-1];
#endif
    data+=s+1;
  }
  /* Integrate to get the metric component. */
  if(params->errors){
    REAL8 aErr=0,bErr=0,cErr=0;
    REAL8 aAvg=Average(a,&aErr);
    REAL8 bAvg=Average(b,&bErr);
    REAL8 cAvg=Average(c,&cErr);
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-dt*bAvg*cAvg);
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
      (aErr + dt*bErr*fabs(cAvg) + dt*cErr*fabs(bAvg));
  }else
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
      (Average(a,0)-dt*Average(b,0)*Average(c,0));

  for(i=1;i<s;i++){
    /* Now generate the gi0 components. */
    /* Fill integrand arrays. */
    data=dSeq->data+1;
    l=COHERENTMETRICC_NPTS;
    while(l--){
#ifndef APPROXIMATE
      a->data[l]=*data*data[-1]*f0*data[i];
      b->data[l]=*data*data[-1];
      c->data[l]=*data*f0*data[i];
#else
      a->data[l]=data[-1]*f0*data[i];
      b->data[l]=data[-1];
      c->data[l]=f0*data[i];
#endif
      data+=s+1;
    }
    /* Integrate to get the metric component. */
    if(params->errors){
      REAL8 aErr=0,bErr=0,cErr=0;
      REAL8 aAvg=Average(a,&aErr);
      REAL8 bAvg=Average(b,&bErr);
      REAL8 cAvg=Average(c,&cErr);
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-dt*bAvg*cAvg);
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	(aErr + dt*bErr*fabs(cAvg) + dt*cErr*fabs(bAvg));
    }else
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	(Average(a,0)-dt*Average(b,0)*Average(c,0));

    for(j=1;j<=i;j++){
      /* Now generate the gij components. */
      /* Fill integrand arrays. */
      data=dSeq->data+1;
      l=COHERENTMETRICC_NPTS;
      while(l--){
#ifndef APPROXIMATE
	a->data[l]=*data*f0*f0*data[i]*data[j];
	b->data[l]=*data*f0*data[i];
	c->data[l]=*data*f0*data[j];
#else
	a->data[l]=f0*f0*data[i]*data[j];
	b->data[l]=f0*data[i];
	c->data[l]=f0*data[j];
#endif
	data+=s+1;
      }
      /* Integrate to get the metric component. */
      if(params->errors){
	REAL8 aErr=0,bErr=0,cErr=0;
	REAL8 aAvg=Average(a,&aErr);
	REAL8 bAvg=Average(b,&bErr);
	REAL8 cAvg=Average(c,&cErr);
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-dt*bAvg*cAvg);
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	  (aErr + dt*bErr*fabs(cAvg) + dt*cErr*fabs(bAvg));
      }else
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	  (Average(a,0)-dt*Average(b,0)*Average(c,0));
    }
  }

  /* Destroy temporary storage, and exit. */
  TRY(LALDDestroyVector(stat->statusPtr,&a),stat);
  TRY(LALDDestroyVector(stat->statusPtr,&b),stat);
  TRY(LALDDestroyVector(stat->statusPtr,&c),stat);
  TRY(LALDDestroyVectorSequence(stat->statusPtr,&dSeq),stat);
  TRY(LALDDestroyVector(stat->statusPtr,&variables),stat);
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


static REAL8
Average( REAL8Vector *integrand, REAL8 *uncertainty )
{
  INT4 i=integrand->length;
  REAL8 *data=integrand->data;
  REAL8 integral=0, error=0;

  /* Make sure i is at least 2 and data exist. */
  if((--i<=0) || (!data))
    return 0;

  /* Return the average, weighting ends by half. */
  integral=*(data++)*0.5;
  while(--i)
    integral+=*(data++);
  integral+=*data*0.5;

  /* If *uncertainty is non-null, and if there are at least three
     points, estimate the error in the trapezoidal integral. */
  if(uncertainty){
    if((i=integrand->length-2)<=0)
      *uncertainty=0.0;
    else{
      data=integrand->data+1;
      while(i--){
	REAL8 concavity=fabs(data[-1]-2.0*data[0]+data[1]);
	if(error<concavity)
	  error=concavity;
	data++;
      }
      *uncertainty=error;
    }
  }

  return integral/(integrand->length-1);
}

#undef COHERENTMETRICC_NPTS
