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

\subsubsection*{Description}

The equation which the code evaluates is as follows:

\begin{eqnarray}
\phi(f) & = & \frac{16 \pi f_{a} \tau_{0}}{5} 
\left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] + 4
\pi f_{a}\tau_{2} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1} \right]
\nonumber  \\
        & - & 5 \pi f_{a} \tau_{3} \left[ 1 - 
              \left( \frac{f}{f_{a}} \right)^{-2/3} \right] + 8 \pi
              f_{a} \tau_{4} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1/3} \right] + \Phi \,.
\label{phioff}
\end{eqnarray}

The terms in this equation are defined as follows:
$\tau_{0}$ is the Newtonian chirp time, $\tau_{2}$ is the 
first post--Newtonian chirp time, and so on for $\tau_{3}$ and $\tau_{4}$. 
The parameter $f$ is the instantaneous frequency of the
gravitational wave, and $f_{a}$ is the value of this frequency when 
the wave enters the lower end of the detectors' bandwidth.
Here we will neglect the $\Phi$ term, which is an initial phase 
which can be defined alsewhere. If we define
\begin{equation}
\phi_{0} = \frac{16 \pi f_{a} \tau_{0}}{5} \,\,,
\end{equation}
\begin{equation}
\phi_{2} = 4 \pi f_{a} \tau_{2} \,\,,
\end{equation}
\begin{equation}
\phi_{3} = 5 \pi f_{a} \tau_{3} \,\,,
\end{equation}
\begin{equation}
\phi_{4} = 8 \pi f_{a} \tau_{4}
\end{equation}
and
\begin{equation}
\phi_{c} = \phi_{0} + \phi_{2} - \phi_{3} + \phi_{4}
\end{equation}
then Eq.(\ref{phioff}) becomes
\begin{eqnarray}
\phi(f) & = & \phi_{0} - \phi_{0}
\left( \frac{f}{f_{a}} \right)^{-5/3} + \phi_{2} - \phi_{2} 
\left( \frac{f}{f_{a}} \right)^{-1} \nonumber \\
      & - & \phi_{3} + \phi_{3} \left( \frac{f}{f_{a}} \right)^{-2/3} + 
            \phi_{4} - \phi_{4} \left( \frac{f}{f_{a}} \right)^{-1/3} + \Phi
\end{eqnarray}
i.e.\
\begin{eqnarray}
\phi(f) & = & - \phi_{0}\left( \frac{f}{f_{a}} \right)^{-5/3}- 
                \phi_{2} \left( \frac{f}{f_{a}} \right)^{-1} \nonumber \\
        & + & \phi_{3} \left( \frac{f}{f_{a}} \right)^{-2/3} - \phi_{4} 
              \left( \frac{f}{f_{a}} \right)^{-1/3} + \phi_{c} + \Phi
\label{phioff2}
\end{eqnarray}

As $f \rightarrow \infty$, then $\phi(f) \rightarrow \phi_{c} + \Phi$. 
If $f=f_{a}$ then $\phi(f=f_{a}) = \Phi$. If $f=f_{lso}$ then $\phi(f=f_{lso})$ 
is given by Eq(\ref{phioff2}). evaluated with $f=f_{lso}$.
 
\subsubsection*{Algorithm}

\subsubsection*{Uses}
For a fuller description of how this function is used in the generation of 
an inspiral waveform, see the documentation for the function \texttt{LALInspiralWave2()}.
The nomenclature adopted is the same as that used in Sathyaprakash, 
PRD, 50, R7111, 1994, which may be consulted for further details.

\subsubsection*{Notes}

See Sathyaprakash, PRD, 50, R7111, 1994, Eq.(5) or document for
module \texttt{LALInspiralWave2} for further details.

\vfill{\footnotesize\input{LALInspiralPhasing2CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHASING2C, "$Id$");

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_0PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v5;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v5 = pow(v,5.);
  *phase = ak->phiC
         + ak->pvaN / v5;

  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_2PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v5;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v5 = v2*v2*v;
  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. + 
         + ak->pva2 * v2);

  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_3PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v5;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);


  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  v2 = v*v;
  v3 = v2*v;
  v5 = v3*v2;

  *phase = ak->phiC
         + ak->pvaN / v5 * ( 1. + 
         + ak->pva2 * v2
         + ak->pva3 * v3);

  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_4PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
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
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_5PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
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

  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_6PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5,v6;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
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

  RETURN(status);

}

/*  <lalVerbatim file="LALInspiralPhasing2CP"> */

void LALInspiralPhasing2_7PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak)
{ /* </lalVerbatim>  */

  REAL8 v2,v3,v4,v5,v6,v7;

  INITSTATUS (status, "LALInspiralPhasing2", LALINSPIRALPHASING2C);
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

  RETURN(status);
}
