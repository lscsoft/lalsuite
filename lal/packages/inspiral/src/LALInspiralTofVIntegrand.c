/*  <lalVerbatim file="LALInspiralTofVIntegrandCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTofVIntegrand.c}}

The function \texttt{LALInspiralTofVIntegrand} calculates the quantity $E^{\prime}(v)/\mathcal{F}(v)$.
These are the energy and flux functions which are used in the gravitational wave phasing formula.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralTofVIntegrandCP}
\index{\verb&LALInspiralTofVIntegrand()&}

\subsubsection*{Description}

The function \texttt{LALInspiralTofVIntegrand} calculates the quantity $E^{\prime}(v)/\mathcal{F}(v)$. 
These are the energy and flux functions which are used in the gravitational wave phasing formula, which is
defined as 

\begin{eqnarray}
t(v) & = & t_{\rm ref} + m \int_v^{v_{\rm ref}} \,
\frac{E'(v)}{{\cal F}(v)} \, dv, \nonumber \\
\phi (v) & = & \phi_{\rm ref} + 2 \int_v^{v_{\rm ref}}  v^3 \,
\frac{E'(v)}{{\cal F}(v)} \, dv,
\label{phasing formula}
\end{eqnarray}

where $v=(\pi m F)^{1/3}$ is an invariantly defined velocity, $F$ is the instantaneous GW frequency, and
$m$ is the total mass of the binary.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

This function calls the function which represents $E^{\prime}(v)$ and $\mathcal{F}(v)$. The pointer to each
of these functions is set by a call to the function \texttt{CHooseModel}.
 
\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralTofVIntegrandCV}}

</lalLaTeX>  */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALTOFVINTEGRANDC, "$Id$"); 

/*  <lalVerbatim file="LALInspiralTofVIntegrandCP"> */
void LALInspiralTofVIntegrand (LALStatus *status,
		       REAL8 *integrand,
		       REAL8 v,
		       void *params)
{ /* </lalVerbatim>  */

   TofVIntegrandIn *ak;

   INITSTATUS (status, "LALInspiralTofVIntegrand", LALINSPIRALTOFVINTEGRANDC);
   ASSERT (integrand, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(v > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(v < 1., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   ak = (TofVIntegrandIn *) params;
   *integrand = ak->dEnergy(v, ak->coeffs)/ak->flux(v, ak->coeffs);

   RETURN (status);
}
