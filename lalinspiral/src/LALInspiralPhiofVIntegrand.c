/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash
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

This function calls {\tt dEnergy} and {\tt flux} functions that are defined in the
{\tt expnFunc} structure  and represent $E^{\prime}(v)$ and $\mathcal{F}(v)$, respectively,
and pointed to the appropriate PN functions with a call to \texttt{LALInspiralChooseModel.}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralPhiofVIntegrandCV}}

</lalLaTeX>  */












#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHIOFVINTEGRANDC, "$Id$");

/*  <lalVerbatim file="LALInspiralPhiofVIntegrandCP"> */
void
LALInspiralPhiofVIntegrand (
   LALStatus  *status,
   REAL8      *integrand,
   REAL8       v,
   void       *params
   )
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

