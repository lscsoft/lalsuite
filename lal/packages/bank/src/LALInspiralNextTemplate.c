/*  <lalVerbatim file="LALInspiralNextTemplateCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralNextTemplate.c}}

Module to calculate the coordinates of the next template in the bank.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralNextTemplateCP}
\index{\verb&LALInspiralNextTemplate()&}

\subsubsection*{Description}

Module to calculate the coordinates of the next template in the bank.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralNextTemplateCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

NRCSID (LALINSPIRALNEXTTEMPLATEC, "Id: $");

/*  <lalVerbatim file="LALInspiralNextTemplateCP"> */
void LALInspiralNextTemplate(LALStatus *status, 
                             InspiralBankParams *bankPars,
                             InspiralMetric metric)
{ /* </lalVerbatim> */

   REAL8 x0tmp, myphi, theta;
   INT4 k;

   INITSTATUS (status, "LALInspiralNextTemplate", LALINSPIRALNEXTTEMPLATEC);
   myphi = fabs(atan(bankPars->dx1/bankPars->dx0));
   theta = metric.theta;
   if (theta > myphi) {
      x0tmp = bankPars->x0 + bankPars->dx0 * cos(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = x0tmp - k * bankPars->dx0;
   } else if (theta > 0 && theta < myphi) {
      x0tmp = bankPars->x0 - bankPars->dx1 * sin(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = (k+1) * bankPars->dx0 - x0tmp;
   } else if (-theta < myphi) {
      x0tmp = bankPars->x0 - bankPars->dx1 * sin(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = x0tmp - k * bankPars->dx0;
   } else {
      x0tmp = bankPars->x0 - bankPars->dx1 * cos(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = (k+1) * bankPars->dx0 - x0tmp;
   }
   RETURN(status);
}
