/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash, Duncan Brown
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

/*  <lalVerbatim file="LALInspiralFrequency3CV">

Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralFrequency3.c}}

The code \texttt{LALInspiralFrequency3.c} calculates the frequency the
waveform from an inspiralling binary system as a function of time up to 3.5
post-Nowtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralFrequency3CP}
\index{\verb&LALInspiralFrequency3()&}
\begin{itemize}
\item {\tt frequency:} Output containing the inspiral waveform.
\item {\tt td:} Input containing PN expansion coefficients $F_k$ (cf. Table \ref{table:flux})
of frequency as a function of time.
\item {\tt ak:} Input containing all PN expansion coefficients.
\end{itemize}

\subsubsection*{Description}

This module computes the instantaneous frequency of an inspiral wave using
\begin{equation}
F(t) = F_N(\theta) \sum F_k \theta^k,
\end{equation}
where the expansion coefficients $F_k,$ Newtonian value $F_N$ and the
time-variable $\theta$ are defined in Table \ref{table:flux}.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}
The frequency evolution defined by post-Newtonian expansion is not monotonic.
Indeed, the equations become highly inaccurate close to the last stable orbit (lso)
and breakdown at or slightly after lso, and the frequency begins to decrease at later times.
It turns out that the evolution is monotonic at least up to lso.

\vfill{\footnotesize\input{LALInspiralFrequency3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALFREQUENCY3C, "$Id$");

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_0PN (
   LALStatus  *status,
   REAL8      *frequency,
   REAL8      td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta3;

  INITSTATUS (status, "LALInspiralFrequency3_0PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta3 = theta*theta*theta;

  *frequency = theta3*ak->ftaN;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_2PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3;

  INITSTATUS (status, "LALInspiralFrequency3_2PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2);


  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_3PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3;

  INITSTATUS (status, "LALInspiralFrequency3_3PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_4PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4;

  INITSTATUS (status, "LALInspiralFrequency3_4PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_5PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5;

  INITSTATUS (status, "LALInspiralFrequency3_5PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_6PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6;

  INITSTATUS (status, "LALInspiralFrequency3_6PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void
LALInspiralFrequency3_7PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak
   )
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;

  INITSTATUS (status, "LALInspiralFrequency3_7PN", LALINSPIRALFREQUENCY3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4
             + ak->fta5*theta5
             + (ak->fta6 + ak->ftl6*log(td))*theta6
             + ak->fta7*theta7);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
