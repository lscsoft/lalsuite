/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*  <lalVerbatim file="LALInspiralNextTemplateCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralNextTemplate.c}}

Routine to compute the parameters of the template next to the
current template, but in the positive $\tau_{2(3)}$ axis.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralNextTemplateCP}
\idx{LALInspiralNextTemplate()}
\begin{itemize}
   \item \texttt{bankPars,} Output, the parameters of the bank at the next grid point; the
   point may, indeed, lay outside the parameter space.
   \item \texttt{matrix} Input, the value of the metric which would allow computation of the
   next lattice point (in the $\tau_{2(3)}$ direction).
\end{itemize}

\subsubsection*{Description}

The coarse grid algorithm works by starting at one corner of the
parameter space, incrementing along positive $\tau_0$ direction,
with increments determined by the local value of the metric,
till the boundary of the parameter space is reached. It then gets back to
the starting point and increments along positive $\tau_{2(3)}$
direction, with an increment defined by the metric defined locally;
it starts at the first point inside the parameter space but
{\it consistent} with a square lattice. This routine is called each
time a translation along the $\tau_{2(3)}$ direction is required.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralNextTemplateCV}}

</lalLaTeX>  */



#include <lal/LALInspiralBank.h>

NRCSID (LALINSPIRALNEXTTEMPLATEC, "Id: $");

/*  <lalVerbatim file="LALInspiralNextTemplateCP"> */

void LALInspiralNextTemplate(LALStatus          *status,
                             InspiralBankParams *bankPars,
                             InspiralMetric     metric)
{ /* </lalVerbatim> */

   REAL8 x0tmp, myphi, theta;
   INT4 k;

   INITSTATUS (status, "LALInspiralNextTemplate", LALINSPIRALNEXTTEMPLATEC);
   ATTATCHSTATUSPTR(status);
   ASSERT (bankPars,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (bankPars->dx0 != 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

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
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
