/*  <lalVerbatim file="LALInspiralMomentsCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralMoments.c}}

Module to calculate the moment of the noise power spectral density. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralMomentsCP}
\index{\verb&LALInspiralMoments()&}

\subsubsection*{Description}

The moments of the noise curve are defined as
\begin{equation}
I(q)  \equiv S_{h}(f_{0}) \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{s_{h}(x f_{0})} \, dx \,.
\end{equation}
Because in practice we will always divide one of these moments by another, we do not need to include
the $S_{h}(f_{0})$ term, which always cancels.
This function calculates the integral
\begin{equation}
I = \int^{f_{c}/f_{0}}_{f_{s}/f_{0}} \frac{x^{-q/3}}{s_{h}(x f_{0})} \, dx \,.
\end{equation}
It then divides this quantity by a normalisation constant which has been passed to the function. In
the case of calculating the components of the metric for the signal manifold for the purpose of
generating a template bank, this constant is given by $I(7)$, because of the definition of the quantity
\begin{equation}
J(q) \equiv \frac{I(q)}{I(7)} \,.
\end{equation}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
LALInspiralMomentsIntegrand
LALDRombergIntegrate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralMomentsCV}}

</lalLaTeX>  */




/* ---------------------------------------------------------------------
  Author: B.S.Sathyaprakash, U.W.C., Cardiff.
  Revision History: Last Update 15.7.97.
  Purpose: Given the exponent {\sf ndx} and limits of integration 
     {\sf xs} and {\sf xc} this function returns the moment of the 
     power spectral density specified by {\sf choice}.
  Inputs: 
  Outputs:
     ${\sf moments} = \int_{xs}^{xc} x^{\sf -ndx}~dx/S_h(x).$
*----------------------------------------------------------------------- */

#include <lal/LALInspiralBank.h>
#include <lal/Integrate.h>

NRCSID(LALINSPIRALMOMENTSC, "$Id$");

/*  <lalVerbatim file="LALInspiralMomentsCP"> */

void LALInspiralMoments(LALStatus         *status,
                        REAL8             *moment,
                        InspiralMomentsIn pars)
{ /* </lalVerbatim> */

   DIntegrateIn In;
   void *funcparams;
   InspiralMomentsIntegrandIn intIn;

   INITSTATUS (status, "LALInspiralMoments", LALINSPIRALMOMENTSC);
   ATTATCHSTATUSPTR(status);
   intIn.ndx = pars.ndx;

   ASSERT (pars.NoisePsd, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   intIn.NoisePsd = pars.NoisePsd;
   funcparams = (void *) &(intIn);
   In.function = LALInspiralMomentsIntegrand;
   In.xmin = pars.xmin;
   In.xmax = pars.xmax;
   In.type = ClosedInterval;

   LALDRombergIntegrate (status->statusPtr, moment, &In, funcparams);
   CHECKSTATUSPTR(status);

   *moment /= pars.norm;
   DETATCHSTATUSPTR(status);
   RETURN (status);
}
