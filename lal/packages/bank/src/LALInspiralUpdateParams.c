/*  <lalVerbatim file="LALInspiralUpdateParamsCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralUpdateParams.c}}

Module to update the parameters used in creating a coarse bank. 
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralUpdateParamsCP}
\index{\verb&LALInspiralUpdateParams()&}

\subsubsection*{Description}
While scanning the $\tau_0$-direction after reaching the 
boundary of the parameter space, we have to return to the 
starting point of the same line and use the metric there
to increment one step upwards in the direction of $\tau_{2(3)}.$ 
to a {\it template list}.


\subsubsection*{Algorithm}

Copy the parameters in the temporary parameter structure
to the current parameter structure.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralUpdateParamsCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

NRCSID (LALINSPIRALUPDATEPARAMSC, "Id: $");

/*  <lalVerbatim file="LALInspiralUpdateParamsCP"> */
void LALInspiralUpdateParams(LALStatus          *status,
                             InspiralBankParams *bankParams,
                             InspiralMetric     metric,
                             REAL8              minimalmatch)
{  /*  </lalVerbatim>  */
   REAL8 dx0, dx1, myphi, theta, fac;

   INITSTATUS (status, "LALInspiralUpdateParams", LALINSPIRALUPDATEPARAMSC);
   ATTATCHSTATUSPTR(status);
   ASSERT (bankParams,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (metric.g00 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (metric.g11 > 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (minimalmatch < 1., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (minimalmatch > 0., status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (metric.theta != LAL_PI_2, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (metric.theta != -LAL_PI_2, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   dx0 = sqrt(2. * (1. - minimalmatch)/metric.g00 );
   dx1 = sqrt(2. * (1. - minimalmatch)/metric.g11 );
   myphi = atan2(dx1, dx0);
   theta = fabs(metric.theta);
   if (theta <= myphi) {
      fac = cos(theta);
      bankParams->dx0 = dx0 / fac;
      bankParams->dx1 = dx1 * fac;
   } 
   else {
      fac = sin(theta);
      bankParams->dx0 = dx1 / fac;
      bankParams->dx1 = dx0 * fac;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);

}
