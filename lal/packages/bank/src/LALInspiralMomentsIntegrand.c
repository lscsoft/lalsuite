/*  <lalVerbatim file="LALInspiralMomentsIntegrandCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralMomentsIntegrand.c}}


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralMomentsIntegrandCP}
\index{\texttt{LALInspiralMomentsIntegrand()}}

\subsubsection*{Description}

The moments of the noise curve are defined as 
\begin{equation}
I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{s_{h}(x f_{0})} \, dx \,.
\end{equation}
This function calculates the integrand of this integral, i.e.\ for a given $x$ it calculates
\begin{equation}
\frac{x^{-q/3}}{s_{h}(xf_{0})} \,\,.
\end{equation}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralMomentsIntegrandCV}}

</lalLaTeX>  */




#include <lal/LALInspiralBank.h>

/*  <lalVerbatim file="LALInspiralMomentsIntegrandCP"> */
void LALInspiralMomentsIntegrand (LALStatus *status,
                                  REAL8 *integrand,
                                  REAL8 x,
                                  void *params)
{ /* </lalVerbatim> */

   InspiralMomentsIntegrandIn *pars;
   status = NULL;
   pars = (InspiralMomentsIntegrandIn *) params;
   *integrand = pow(x,-pars->ndx) / pars->NoisePsd(x);
}
