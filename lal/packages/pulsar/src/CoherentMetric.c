/******************************* <lalVerbatim file="CoherentMetricCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{CoherentMetric.c}}
\label{ss:CoherentMetric.c}

Computes the parameter space metric for a coherent pulsar search.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CoherentMetricCP}
\index{\texttt{LALCoherentMetric()}}

\subsubsection*{Description}

This function computes the metric $g_{\alpha\beta}(\bm{\lambda})$, as
discussed in the header \verb@StackMetric.h@, under the assumption
that the search consists of scanning the Fourier power spectrum
constructed from a single time interval $\Delta t$.  The indecies
$\alpha$ and $\beta$ are assumed to run from 0 to $n$, where $n$ is
the total number of ``shape'' parameters.

The argument \verb@*metric@ is normally a vector of length
$(n+1)(n+2)/2$ storing all non-redundant coefficients of
$g_{\alpha\beta}$.  The indexing scheme is as follows: Let us assume
that $\alpha\geq\beta$.  Then $g_{\alpha\beta} = g_{\beta\alpha} =
$\verb@metric->data[@$\beta + \alpha(\alpha+1)/2$\verb@]@.  If
\verb@params->errors@ is nonzero, then \verb@*metric@ must be double
this length, and \verb@LALCoherentMetric()@ will store metric
components and their estimated uncertainty in alternate slots; i.e.\
the metric component
$g_{\alpha\beta}=$\verb@metric->data[@$2\beta+\alpha(\alpha+1)$@, and
the uncertainty
$s_{\alpha\beta}=$\verb@metric->data[@$1+2\beta+\alpha(\alpha+1)$@.

The argument \verb@*lambda@ is another vector, of length $n+1$,
storing the components of $\bm{\lambda}=(\lambda^0,\ldots,\lambda^n)$
for the parameter space point at which the metric is being evaluated.
The argument \verb@*params@ stores the remaining parameters for
computing the metric, as given in the Structures section of
\verb@StackMetric.h@.

\subsubsection*{Algorithm}

This routne simply computes the function given in Eq.~\ref{eq:gab-phi}
of header \verb@StackMetric.h@.  Most of the work is done by the
function \verb@params->dtCanon()@, which computes the value and
parameter derivatives of the canonical time coordinate
$\tau[t;\vec\lambda]$.  Since the phase function is simply
$\phi[t;\bm{\lambda}]=\lambda^0 \tau[t;\vec\lambda]$, the metric
components are simply:
\begin{eqnarray}
g_{00}(\bm\lambda) & = &
	4\pi^2\bigg\langle\!(\tau[t;\vec\lambda])^2\bigg\rangle -
	4\pi^2\bigg\langle\!\tau[t;\vec\lambda]\bigg\rangle^2 \; ,
	\nonumber\\
g_{i0}(\bm\lambda) \;\; = \;\;
g_{0i}(\bm\lambda) & = &
	4\pi^2 f_0\bigg\langle\!
	\tau[t;\vec\lambda]
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
	\bigg\rangle
	-
	4\pi^2 f_0\bigg\langle\!
	\tau[t;\vec\lambda]
	\bigg\rangle
	\bigg\langle
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
	\bigg\rangle \; , \nonumber\\
g_{ij}(\bm\lambda) & = &
	4\pi^2 f_0^2\bigg\langle
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^j}
	\bigg\rangle
	-
	4\pi^2 f_0^2\bigg\langle
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^i}
	\bigg\rangle
	\bigg\langle
	\frac{\partial\tau[t;\vec\lambda]}{\partial\lambda^j}
	\bigg\rangle \; , \nonumber
\end{eqnarray}
where the indecies $i$ and $j$ run from 1 to $n$.

In the rigorous definition of the metric, the angle brackets denote an
average over the \emph{canonical} time, not the detector time, so we
have:
\begin{eqnarray}
\bigg\langle\ldots\bigg\rangle
& = & \frac{1}{\tau(t_\mathrm{start}+\Delta t)-\tau(t_\mathrm{start})}
	\int_{\tau(t_\mathrm{start})}^{\tau(t_\mathrm{start}+\Delta t)}
	\ldots\,d\tau \nonumber\\
& = & \frac{1}{\tau(t_\mathrm{start}+\Delta t)-\tau(t_\mathrm{start})}
	\int_{t_\mathrm{start}}^{t_\mathrm{start}+\Delta t}
	\ldots\frac{\partial\tau}{\partial t}\,dt \nonumber\\
& \approx & \frac{1}{\Delta t}
	\int_{t_\mathrm{start}}^{t_\mathrm{start}+\Delta t}
	\ldots\,dt \nonumber \; ,
\end{eqnarray}
where the approximation is good to order $\epsilon$ equal to the
maximum difference between 1 and $\partial\tau/\partial t$.  For an
Earth-motion barycentred time coordinate, for instance, $\epsilon$ is
of order $10^{-4}$ and the approximation is quite good.  However, the
current implementation of the metric algorithm uses the second, exact
formula.  If speed considerations become important, this is one
obvious area for simplification.

At present, the time averaging is performed by evaluating $\tau$ and
its derivatives at equally-spaced points over $\Delta t$ and applying
a trapezoidal method; i.e.
$$
\frac{1}{T}\int_0^T F(t)\,dt \approx \frac{1}{N}\left(\frac{1}{2}F_0
	+ F_1 + F_2 + \ldots + F_{N-1} + \frac{1}{2}F_N \right) \; ,
$$
where $F_k=F(kT/N)$ are the $N+1$ sample points of the integrand.  The
number of points is a compiled-in constant.  The error is on the order
of $F''T^2/N^2$, where $F''$ is the second derivative of $F(t)$ at
some point in the interval; we liberally estimate the error to be:
$$
\sigma = \max_{k\in\{1,\ldots,N-1\}}
	\left\{F_{k-1} - 2F_k + F_{k+1}\right\} \; .
$$
Other more sophisticated integration techniques may be considered in
future.  To save on recomputing costs, the derivatives of $\tau$ at
all times are stored in a large vector sequence.

\subsubsection*{Uses}
\begin{verbatim}
LALDCreateVector()              LALDDestroyVector()
LALDCreateVectorSequence()      LALDDestroyVectorSequence()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CoherentMetricCV}}

******************************************************* </lalLaTeX> */

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/LALConstants.h>
#include<lal/AVFactories.h>
#include<lal/SeqFactories.h>
#include<lal/StackMetric.h>

NRCSID(COHERENTMETRICC,"$Id$");

#define COHERENTMETRICC_NPTS 100000
/* Number of points to average per time interval. */

static REAL8
Average(REAL8Vector *integrand, REAL8 *uncertainty);
/* Local function to perform the averaging. */

/* <lalVerbatim file="CoherentMetricCP"> */
void
LALCoherentMetric( LALStatus        *stat,
		   REAL8Vector      *metric,
		   REAL8Vector      *lambda,
		   MetricParamStruc *params )
{ /* </lalVerbatim> */
  UINT4 s;          /* The number of parameters. */
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

  INITSTATUS(stat,"LALCoherentMetric",COHERENTMETRICC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT(metric,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(metric->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(lambda,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(lambda->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(params,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(params->dtCanon,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);

  /* Make sure parameters are valid. */
  ASSERT(s=lambda->length,stat,STACKMETRICH_EBAD,STACKMETRICH_MSGEBAD);
  /* (yes, that is an assignment, not a comparison) */
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
    TRY((params->dtCanon)(stat->statusPtr,&d,variables,constants),
	stat);
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
  dt=params->deltaT/(dSeq->data[(s+1)*(COHERENTMETRICC_NPTS-1)]
		     - dSeq->data[0]);

  /* Compute metric components.  First generate the g00 component. */
  /* Fill integrand arrays. */
  data=dSeq->data+1;
  l=COHERENTMETRICC_NPTS;
  k=0;
  while(l--){
    a->data[l]=*data*data[-1]*data[-1];
    b->data[l]=c->data[l]=*data*data[-1];
    data+=s+1;
  }
  /* Integrate to get the metric component. */
  if(params->errors){
    REAL8 aErr,bErr,cErr;
    REAL8 aAvg=Average(a,&aErr);
    REAL8 bAvg=Average(b,&bErr);
    REAL8 cAvg=Average(c,&cErr);
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-bAvg*cAvg);
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
      (aErr + bErr*fabs(cAvg) + cErr*fabs(bAvg));
  }else
    metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
      (Average(a,0)-Average(b,0)*Average(c,0));

  for(i=1;i<s;i++){
    /* Now generate the gi0 components. */
    /* Fill integrand arrays. */
    data=dSeq->data+1;
    l=COHERENTMETRICC_NPTS;
    while(l--){
      a->data[l]=*data*data[-1]*f0*data[i];
      b->data[l]=*data*data[-1];
      c->data[l]=*data*f0*data[i];
      data+=s+1;
    }
    /* Integrate to get the metric component. */
    if(params->errors){
      REAL8 aErr,bErr,cErr;
      REAL8 aAvg=Average(a,&aErr);
      REAL8 bAvg=Average(b,&bErr);
      REAL8 cAvg=Average(c,&cErr);
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-bAvg*cAvg);
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	(aErr + bErr*fabs(cAvg) + cErr*fabs(bAvg));
    }else
      metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	(Average(a,0)-Average(b,0)*Average(c,0));

    for(j=1;j<=i;j++){
      /* Now generate the gij components. */
      /* Fill integrand arrays. */
      data=dSeq->data+1;
      l=COHERENTMETRICC_NPTS;
      while(l--){
	a->data[l]=*data*f0*f0*data[i]*data[j];
	b->data[l]=*data*f0*data[i];
	c->data[l]=*data*f0*data[j];
	data+=s+1;
      }
      /* Integrate to get the metric component. */
      if(params->errors){
	REAL8 aErr,bErr,cErr;
	REAL8 aAvg=Average(a,&aErr);
	REAL8 bAvg=Average(b,&bErr);
	REAL8 cAvg=Average(c,&cErr);
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*(aAvg-bAvg*cAvg);
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	  (aErr + bErr*fabs(cAvg) + cErr*fabs(bAvg));
      }else
	metric->data[k++]=LAL_TWOPI*LAL_TWOPI*dt*
	  (Average(a,0)-Average(b,0)*Average(c,0));
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
