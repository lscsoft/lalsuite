/*  <lalVerbatim file="LALInspiralPhasing3CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralPhasing3.c}}

The code \texttt{LALInspiralPhasing3.c} calculates the phase the waveform 
from an inspiralling binary system as a function of time up to second post-Nowtonian 
order.

\subsubsection*{Prototypes}

\vspace{0.1in}

\input{LALInspiralPhasing3CP}

\index{\verb&LALInspiralPhasing3()&}

\subsubsection*{Description}

The code \texttt{LALInspiralPhasing3.c} calculates the phase the waveform 
from an inspiralling binary system as a function of time up to second 
post-Nowtonian order. The method used is as follows.

From Blanchet, Iyer, Will and Wiseman, CQG \textbf{13}, 575, 1996, the 
instantaneous orbital phase $\phi$ is given in terms of the dimensionless 
time variable $\Theta$  by

\begin{eqnarray}
\phi(t) = \phi_{c} - &  \frac{1}{\eta} \left\{ \Theta^{5/8} + 
          \left( \frac{3715}{8064} + \frac{55}{96} \eta \right) \Theta^{3/8} - 
          \frac{3 \pi}{4} \Theta^{1/4} \right.  \\
        &  + \left. \left( \frac{9275495}{14450688} + \frac{284875}{258048} \eta +
          \frac{1855}{2048} \eta^{2} \right) \Theta^{1/8} \right\}
\label{phioft}
\end{eqnarray}
where $\phi_{c}$ is a constant which represents the value of the phase at 
instant $t_{c}$, which is the instant of coalescence of the two point-masses 
which constitute the binary.  The dimensionless time variable $\Theta$ is given by
\begin{equation}
\Theta = \frac{c^{3}_{0} \eta}{5Gm} (t_{c} - t) \,\,.
\end{equation}

All of the equations presented so far have included explicitly their 
dependence upon $G$ and $c$. The code uses units where $G=c=1$.


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralPhasing3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHASING3C, "$Id$");


/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_0PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta5;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta5 = pow(td,-0.625);
  *phase = (ak->ptaN/theta5);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_2PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta5;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta5 = theta2*theta2*theta;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_3PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta5;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta5 = theta2*theta3;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2 
         + ak->pta3*theta3);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_4PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2 
         + ak->pta3*theta3
         + ak->pta4*theta4);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_5PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2 
         + ak->pta3*theta3
         + ak->pta4*theta4 
         + ak->pta5 * log(td/ak->tn) * theta5);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_6PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2 
         + ak->pta3*theta3
         + ak->pta4*theta4 
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing3CP"> */

void LALInspiralPhasing3_7PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;

  INITSTATUS (status, "LALInspiralPhasing3", LALINSPIRALPHASING3C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;
  theta5 = theta4*theta;
  theta6 = theta5*theta;
  theta7 = theta6*theta;

  *phase = (ak->ptaN/theta5) * (1. 
         + ak->pta2*theta2 
         + ak->pta3*theta3
         + ak->pta4*theta4 
         + ak->pta5*log(td/ak->tn)*theta5
         +(ak->ptl6*log(td/256.) + ak->pta6)*theta6
         + ak->pta7*theta7);
  RETURN(status);
}
