/******************************** <lalVerbatim file="ProjectMetricCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{ProjectMetric.c}}
\label{ss:ProjectMetric.c}

Projects out the zeroth dimension of a metric.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ProjectMetricCP}
\idx{LALProjectMetric()}

\subsubsection*{Description}

This function takes a metric $g_{\alpha\beta}$, where
$\alpha,\beta=0,1,\ldots,n$, and computes the projected metric
$\gamma_{ij}$ on the subspace $i,j=1,\ldots,n$, as described in the
header \verb@StackMetric.h@.

The argument \verb@*metric@ stores the metric components in the manner
used by the functions \verb@CoherentMetric()@ and
\verb@StackMetric()@, and \verb@errors@ indicates whether error
estimates are included in \verb@*metric@.  Thus \verb@*metric@ is a
vector of length $(n+1)(n+2)/2$ if \verb@errors@ is zero, or of length
$(n+1)(n+2)$ if \verb@errors@ is nonzero; see \verb@CoherentMetric.c@
for the indexing scheme.

Upon return, \verb@*metric@ stores the components of $\gamma_{ij}$ in
the same manner as above, with the physically meaningless components
$\gamma_{\alpha0} = \gamma_{0\alpha}$ (and their uncertainties) set
identically to zero.

\subsubsection*{Algorithm}

The function simply implements Eq.~\ref{eq:gij-gab} in header
\verb@StackMetric.h@.  The formula used to convert uncertainties
$s_{\alpha\beta}$ in the metric components $g_{\alpha\beta}$ into
uncertainties $\sigma_{ij}$ in $\gamma_{ij}$ is:
$$
\sigma_{ij} = s_{ij}
	+ s_{0i}\left|\frac{g_{0j}}{g_{00}}\right|
	+ s_{0j}\left|\frac{g_{0i}}{g_{00}}\right|
	+ s_{00}\left|\frac{g_{0i}g_{0j}}{(g_{00})^2}\right| \; .
$$
Note that if the metric is highly degenerate, one may find that one or
more projected components are comparable in magnitude to their
estimated numerical uncertainties.  This can occur when the
observation times are very short or very long compared to the
timescales over which the timing derivatives are varying.  In the
former case, one is advised to use analytic approximations or a
different parameter basis.  In the latter case, the degenerate
components are often not relevant for data analysis, and can be
effectively set to zero.

Technically, starting from a full metric
$g_{\alpha\beta}(\bm{\lambda})$, the projection
$\gamma_{ij}(\vec\lambda)$ is the metric of a subspace
$\{\vec\lambda\}$ passing through the point \bm{\lambda} on a plane
orthogonal to the $\lambda^0$ axis.  In order for $\gamma_{ij}$ to
measure the \emph{maximum} distance between points $\vec\lambda$, it
is important to evaluate $g_{\alpha\beta}$ at the value of $\lambda^0$
that gives the largest possible separations.  For the pulsar search
formalism discussed in the header \verb@StackMetric.h@, this is always
achieved by choosing the largest value of $\lambda^0=f_\mathrm{max}$
that is to be covered in the search.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ProjectMetricCV}}

******************************************************* </lalLaTeX> */

#include<math.h>
#include<lal/LALStdlib.h>
#include<lal/StackMetric.h>

NRCSID(PROJECTMETRICC,"$Id$");

/* <lalVerbatim file="ProjectMetricCP"> */
void
LALProjectMetric( LALStatus *stat, REAL8Vector *metric, BOOLEAN errors )
{ /* </lalVerbatim> */
  UINT4 s;     /* The number of parameters before projection. */
  UINT4 i, j;  /* Indecies. */
  REAL8 *data; /* The metric data array. */

  INITSTATUS(stat,"ProjectMetric",PROJECTMETRICC);

  /* Check that data exist. */
  ASSERT(metric,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  ASSERT(metric->data,stat,STACKMETRICH_ENUL,STACKMETRICH_MSGENUL);
  data=metric->data;

  /* Make sure the metric's length is compatible with some
     dimensionality s. */
  for(s=0,i=1;i<metric->length;s++){
    i=s*(s+1);
    if(!errors)
      i=i>>1;
  }
  s--; /* counteract the final s++ */
  ASSERT(i==metric->length,stat,STACKMETRICH_EBAD,
	 STACKMETRICH_MSGEBAD);

  /* Now project out the zeroth component of the metric. */
  for(i=1;i<s;i++)
    for(j=1;j<=i;j++){
      INT4 i0 = (i*(i+1))>>1;
      INT4 myj0 = (j*(j+1))>>1;
      INT4 ij = i0+j;
      if(errors){
	data[2*ij]-=data[2*i0]*data[2*myj0]/data[0];
	data[2*ij+1]+=data[2*i0+1]*fabs(data[2*myj0]/data[0])
	  + data[2*myj0+1]*fabs(data[2*i0]/data[0])
	  + data[1]*fabs(data[2*i0]/data[0])*fabs(data[2*myj0]/data[0]);
      } else
	data[ij]-=data[i0]*data[myj0]/data[0];
    }

  /* Set all irrelevant metric coefficients to zero. */
  for(i=0;i<s;i++){
    INT4 i0 = i*(i+1)>>1;
    if(errors)
      data[2*i0]=data[2*i0+1]=0.0;
    else
      data[i0]=0.0;
  }

  /* And that's it! */
  RETURN(stat);
}
