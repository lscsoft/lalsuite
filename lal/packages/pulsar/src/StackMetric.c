/********************************** <lalVerbatim file="StackMetricCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{StackMetric.c}}
\label{ss:StackMetric.c}

Computes the parameter space metric for a coherent pulsar search.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{StackMetricCP}
\idx{LALStackMetric()}

\subsubsection*{Description}

This function computes the metric $g_{\alpha\beta}(\bm{\lambda})$, as
discussed in the header \verb@StackMetric.h@, under the assumption
that the detected power is constructed from the incoherent sum of $N$
separate power spectrum, each derived from separate time intervals of
length $\Delta t$.  The indecies $\alpha$ and $\beta$ are assumed to
run from 0 to $n$, where $n$ is the total number of ``shape''
parameters.

This routine has exactly the same calling structure and data storage
as the \verb@CoherentMetric()@ function.  Thus, the argument
\verb@*metric@ is a vector of length $(n+1)(n+2)/2$ storing all
non-redundant coefficients of $g_{\alpha\beta}$, or twice this length
if \verb@params->errors@ is nonzero.  See \verb@CoherentMetric.c@ for
the indexing scheme.  The argument \verb@*lambda@ is another vector,
of length $n+1$, storing the components of
$\bm{\lambda}=(\lambda^0,\ldots,\lambda^n)$ for the parameter space
point at which the metric is being evaluated.  The argument
\verb@*params@ stores the remaining parameters for computing the
metric, as given in the Structures section of \verb@StackMetric.h@.

\subsubsection*{Algorithm}

Most of what this routine does is set up arguments to be passed to the
function \verb@CoherentMetric()@.  Each metric component in the stack
metric is given simply by:
$$
g_{\alpha\beta}(\bm{\lambda}) = \frac{1}{N} \sum_{k=1}^N
	g^{(k)}_{\alpha\beta}(\bm{\lambda}) \; ,
$$
where $g^{(k)}_{\alpha\beta}$ is just the coherent metric computed on
the time interval $[t_\mathrm{start}+(k-1)\Delta t,
t_\mathrm{start}+k\Delta t]$.  The estimated uncertainty
$s_{\alpha\beta}$ in this component is taken to be:
$$
s_{\alpha\beta} = \frac{1}{\sqrt{N}} \max_{k\in\{1,\ldots,N\}}
	s^{(k)}_{\alpha\beta} \; .
$$
There are no clever tricks involved in any of these computations.

\subsubsection*{Uses}
\begin{verbatim}
LALDCreateVector()
LALDDestroyVector()
LALCoherentMetric()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StackMetricCV}}

******************************************************* </lalLaTeX> */

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/AVFactories.h>
#include<lal/StackMetric.h>

NRCSID(STACKMETRICC,"$Id$");

/* <lalVerbatim file="StackMetricCP"> */
void
LALStackMetric( LALStatus        *stat,
		REAL8Vector      *metric,
		REAL8Vector      *lambda,
		MetricParamStruc *params )
{ /* </lalVerbatim> */
  INT4 n;  /* Number of metric coefficients. */
  INT4 i;  /* An index. */
  INT4 j;  /* Another index. */
  REAL8 t; /* Time at start of each coherent stretch. */
  REAL8Vector *subMetric=NULL; /* Coherent metric on each stack. */
  MetricParamStruc subParams;  /* Parameters passed to coherent metric
				  computation. */

  INITSTATUS(stat,"StackMetric",STACKMETRICC);
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
