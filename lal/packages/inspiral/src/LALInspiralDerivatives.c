/*  <lalVerbatim file="LALInspiralDerivativesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralDerivatives.c}}

Module to calculate the RHS of the differential equations 
in Eq.~(\ref{eq:ode2}).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralDerivativesCP}
\idx{LALInspiralDerivatives()}
\begin{itemize}
\item {\tt values:} Input containing the values of the variables $v$ and $\phi$ at the
current time.
\item {\tt dvalues:} Output containing the derivatives $dv/dt$ and $d\phi/dt$ at the
current time.
\item {\tt params:} Input  of type {\tt InspiralDerivativesIn} that must be
cast to a {\tt void.}\\

\end{itemize}

\subsubsection*{Description}

This module calculates the right-hand sides of
the follwoing two coupled first-order differential equations which are
solved to obtain the gravitational wave phasing equation, 
as described in the documentation for the function \texttt{LALInspiralWave1}:
The equations are
\begin{equation}
\frac{dv}{dt} = - \frac{\mathcal{F}(v)}{m E^{\prime}(v)},\ \ \ \ 
\frac{d \phi(t)}{dt} = \frac{2v^{3}}{m}.
\label{ode2}
\end{equation}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\begin{itemize}
\item This function has been intentionally made non-LAL compliant in the sense that it 
has no status structure.  This is because this code
outputs the RHS of the differential equations
and is called repeatedly by a function that integrates the two differential
equations and should therefore not suffer from undue overheads. 
\item The input {\tt params} is of type {\tt InspiralDerivativesIn} and must
be cast to a void before calling this function. For example,\\[5pt]
\texttt {
   InspiralDerivativesIn in3;\\
   void *funcParams;\\[5pt]
   in3.totalmass = totalmass;\\
   $\ldots$\\
   funcParams = (void *) \&in3;
   }

\end{itemize}

\vfill{\footnotesize\input{LALInspiralDerivativesCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALDERIVATIVESC, "$Id$");

/*  <lalVerbatim file="LALInspiralDerivativesCP"> */
void 
LALInspiralDerivatives (
   REAL8Vector *values,
   REAL8Vector *dvalues,
   void        *params
   )
 { /* </lalVerbatim> */

  InspiralDerivativesIn *ak;
  REAL8 v;

  ak = (InspiralDerivativesIn *) params;
  v = *(values->data);

  dvalues->data[0] = -ak->flux(v, ak->coeffs)/ (ak->totalmass*ak->dEnergy(v, ak->coeffs));
  dvalues->data[1] = 2.* pow(v,3.)/ak->totalmass;

  return;
}
