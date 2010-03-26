/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralUpdateParamsCV">
Author: Sathyaprakash, B. S., T. Cokelaer
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralUpdateParams.c}}

Module to update the parameters used in creating a coarse bank
based on a square lattice.
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralUpdateParamsCP}
\idx{LALInspiralUpdateParams()}
\begin{itemize}
   \item \texttt{bankParams} Output, refreshed to get the next location
   \item \texttt{metric} Input, metric at the current location
   \item \texttt{minimalmatch} Input, the minimal match
\end{itemize}

\subsubsection*{Description}
While scanning the $\tau_0$-direction after reaching the
boundary of the parameter space, we have to return to the
starting point of the same line and use the metric there
to increment one step upwards in the direction of $\tau_{2(3)}.$
to a {\it template list}.


The $dx_i$ returned by this function gives the spacing for a
square lattice (e.g., $dx_i$ as given in Owen PRD 53).

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
   ASSERT (metric.theta < LAL_PI_2, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);
   ASSERT (metric.theta > -LAL_PI_2, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   /* This dx0, dx1 are linked to a square placement only !! */
   dx0 = sqrt(2.L * (1.L - minimalmatch)/metric.g00 );
   dx1 = sqrt(2.L * (1.L - minimalmatch)/metric.g11 );

   if (metric.theta==0.L)
   {
	   bankParams->dx0 = dx0;
	   bankParams->dx1 = dx1;
   }
   else
   {
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
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);

}
