/*  <lalVerbatim file="LALInspiralPhiofVIntegrandCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralPhiofVIntegrand.c}}

The function \texttt{LALInspiralPhiofVIntegrandIn} calculates the quantity $v^{3} E^{\prime}(v)/\mathcal{F}(v)$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralPhiofVIntegrandCP}
\index{\verb&LALInspiralPhiofVIntegrand()&}

\subsubsection*{Description}

The function \texttt{LALInspiralPhiofVIntegrandIn} calculates the quantity $v^{3} E^{\prime}(v)/\mathcal{F}(v)$.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

This function calls the function which represents $E^{\prime}(v)$ and $\mathcal{F}(v)$. The pointer to each
of these functions is set by a call to the function \texttt{ChooseModel}.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralPhiofVIntegrandCV}}

</lalLaTeX>  */












#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHIOFVINTEGRANDC, "$Id$"); 

/*  <lalVerbatim file="LALInspiralPhiofVIntegrandCP"> */
void LALInspiralPhiofVIntegrand (LALStatus *status,
		         REAL8 *integrand,
		         REAL8 v,
		         void *params)
{ /* </lalVerbatim>  */

  PhiofVIntegrandIn *in;

  INITSTATUS (status, "LALInspiralPhiofVIntegrand", LALINSPIRALPHIOFVINTEGRANDC);
  ATTATCHSTATUSPTR(status);

  ASSERT (integrand, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (v>0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (v<1, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  in = (PhiofVIntegrandIn *) params;

  *integrand = pow (v, 3.) * in->dEnergy(v,in->coeffs)/in->flux(v,in->coeffs);

  DETATCHSTATUSPTR(status);
  RETURN(status);


}

