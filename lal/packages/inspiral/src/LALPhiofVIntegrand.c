/*  <lalVerbatim file="LALPhiofVIntegrandCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALPhiofVIntegrand.c}}

The function \texttt{LALPhiofVIntegrandIn} calculates the quantity $v^{3} E^{\prime}(v)/\mathcal{F}(v)$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALPhiofVIntegrandCP}
\index{\texttt{LALPhiofVIntegrand()}}

\subsubsection*{Description}

The function \texttt{LALPhiofVIntegrandIn} calculates the quantity $v^{3} E^{\prime}(v)/\mathcal{F}(v)$.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

This function calls the function which represents $E^{\prime}(v)$ and $\mathcal{F}(v)$. The pointer to each
of these functions is set by a call to the function \texttt{ChooseModel}.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALPhiofVIntegrandCV}}

</lalLaTeX>  */












#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALPHIOFVINTEGRANDC, "$Id$"); 

/*  <lalVerbatim file="LALPhiofVIntegrandCP"> */
void LALPhiofVIntegrand (LALStatus *status,
		         REAL8 *integrand,
		         REAL8 v,
		         void *params)
{ /* </lalVerbatim>  */

  PhiofVIntegrandIn *in;

  INITSTATUS (status, "LALPhiofVIntegrand", LALPHIOFVINTEGRANDC);

  ASSERT (integrand, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  in = (PhiofVIntegrandIn *) params;

  *integrand = pow (v, 3.) * in->dEnergy(v,in->coeffs)/in->flux(v,in->coeffs);

  RETURN(status);


}

