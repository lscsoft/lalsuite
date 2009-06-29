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

/*  <lalVerbatim file="LALInspiralPhasing2CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralPhasing2.c}}

The code \texttt{LALInspiralPhasing2.c} calculates the phase of an inspiral
waveform as a function of the
instantaneous frequency of the wave, up to $2^{nd}$ post--Newtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralPhasing2CP}
\index{\verb&LALInspiralPhasing2()&}
\begin{itemize}
\item {\tt phase:} Output, the phase of the wave at the current epoch.
\item {\tt v:} Input, the PN expansion parameter at the current epoch.
\item {\tt ak:} Input containing PN expansion coefficients.
\end{itemize}

\subsubsection*{Description}

The phase of the inspiral wave corresponding to the {\tt Approximant} {\tt TaylorT2}
as in Equation~{eq:InspiralPhasing2}.

\subsubsection*{Algorithm}
None.

\subsubsection*{Uses}
None.

\subsubsection*{Notes}
None.

\vfill{\footnotesize\input{LALInspiralPhasing2CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHASING2C, "$Id$");

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_0PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v5;

  INITSTATUS (status, "LALInspiralPhasing2_0PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v5 = pow(v,5.);
  *phase = ak->phiC
         + ak->pvaN / v5;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_2PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v5;

  INITSTATUS (status, "LALInspiralPhasing2_2PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v5 = v2*v2*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_3PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v5;

  INITSTATUS (status, "LALInspiralPhasing2_3PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v5 = v3*v2;

  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_4PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5;

  INITSTATUS (status, "LALInspiralPhasing2_4PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_5PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5;

  INITSTATUS (status, "LALInspiralPhasing2_5PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_6PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5,v6;

  INITSTATUS (status, "LALInspiralPhasing2_6PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5
         + (ak->pva6 + ak->pvl6*log(4*v)) * v6);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void
LALInspiralPhasing2_7PN (
   LALStatus  *status,
   REAL8      *phase,
   REAL8       v,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5,v6,v7;

  INITSTATUS (status, "LALInspiralPhasing2_7PN", LALINSPIRALPHASING2C);
  ATTATCHSTATUSPTR(status);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. +
         + ak->pva2 * v2
         + ak->pva3 * v3
         + ak->pva4 * v4
         + ak->pva5 * log(v/ak->vlso) * v5
         + (ak->pva6 + ak->pvl6*log(4*v)) * v6
         + ak->pva7 * v7);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
