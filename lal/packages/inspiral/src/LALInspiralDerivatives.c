/*  <lalVerbatim file="LALInspiralDerivativesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralDerivatives.c}}

Module to calculate the first order differential equations needed in the phasing formula.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralDerivativesCP}
\index{\verb&LALInspiralDerivatives()&}

\subsubsection*{Description}

The code \texttt{LALInspiralDerivatives.c} calculates the two coupled first--order differential equations
which we use to solve the gravitational wave phasing equation, as described in the documentation for the
function \texttt{TimeDomain2}.

The equations are

\begin{equation}
\frac{dv}{dt} = - \frac{\mathcal{F}(v)}{m E^{\prime}(v)} \,\,.
\label{ode1}
\end{equation}

and

\begin{equation}
\frac{d \phi(t)}{dt} = \frac{2v^{3}}{m}
\label{ode2}
\end{equation}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

This function has been made non--LAL compliant in the sense that it has no status structure.
If we had included the status structure, the code \texttt{RungeKutta4} which calls \texttt{LALInspiralDerivatives} would have needed an \texttt{ATTATCHSTATUSPTR} command, which
would have been called for every data point on the waveform. This was proving inhibitively slow, so we
removed them.

\vfill{\footnotesize\input{LALInspiralDerivativesCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALDERIVATIVESC, "$Id$");

/*  <lalVerbatim file="LALInspiralDerivativesCP"> */
void LALInspiralDerivatives (REAL8Vector *values,
			     REAL8Vector *dvalues,
			     void *params)
 { /* </lalVerbatim> */

  InspiralDerivativesIn *ak;

  ak = (InspiralDerivativesIn *) params;


   *(dvalues->data) = -ak->flux(*(values->data), ak->coeffs)/
                      (ak->totalmass*ak->dEnergy(*(values->data), ak->coeffs));

   *(dvalues->data+1) = 2.* *(values->data)* *(values->data)* *(values->data)/ak->totalmass;

  return;


}
