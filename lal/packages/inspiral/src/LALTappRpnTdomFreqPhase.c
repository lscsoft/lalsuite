/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomFreqPhase.c}}

The code \texttt{LALTappRpnTdomFreqPhase.c} calculates the phase of an inspiral waveform as a function of the
instantaneous frequency of the wave, up to $2^{nd}$ post--Newtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomFreqPhaseCP}
\index{\verb&LALTappRpnTdomFreqPhase()&}

\subsubsection*{Description}

The equation which the code evaluates is as follows:

\begin{eqnarray}
\phi(f) & = & \frac{16 \pi f_{a} \tau_{N}}{5} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] + 4
\pi f_{a}\tau_{P^{1}N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1} \right]
\nonumber  \\
        & - & 5 \pi f_{a} \tau_{P^{1.5}N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-2/3} \right] + 8 \pi
f_{a} \tau_{P^{2}N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1/3} \right] + \Phi \,.
\label{phioff}
\end{eqnarray}

The terms in this equation are defined as follows:
$\tau_{N}$ is the Newtonian chirp time, $\tau_{P^{1}N}$ is the first post--Newtonian chirp time, and so on
for $\tau_{P^{1.5}N}$ and $\tau_{P^{2}N}$. The parameter $f$ is the instantaneous frequency of the
gravitational wave, and $f_{a}$ is the value of this frequency when the wave enters the lower end of the
detectors' bandwidth.
Here we will neglect the $\Phi$ term, which is an initial phase which can be defined alsewhere. If we
define
\begin{equation}
\phi_{N} = \frac{16 \pi f_{a} \tau_{N}}{5} \,\,,
\end{equation}
\begin{equation}
\phi_{P^{1}N} = 4 \pi f_{a} \tau_{P^{1}N} \,\,,
\end{equation}
\begin{equation}
\phi_{P^{1.5}N} = 5 \pi f_{a} \tau_{P^{1.5}N} \,\,,
\end{equation}
\begin{equation}
\phi_{P^{2}N} = 8 \pi f_{a} \tau_{P^{2}N}
\end{equation}
and
\begin{equation}
\phi_{c} = \phi_{N} + \phi_{P^{1}N} - \phi_{P^{1.5}N} + \phi_{P^{2}N}
\end{equation}
then Eq.(\ref{phioff}) becomes
\begin{eqnarray}
\phi(f) & = & \phi_{N} - \phi_{N}\left( \frac{f}{f_{a}} \right)^{-5/3} + \phi_{P^{1}N} - \phi_{P^{1}N} \left(
\frac{f}{f_{a}} \right)^{-1} \nonumber \\
      & - & \phi_{P^{1.5}N} + \phi_{P^{1.5}N} \left( \frac{f}{f_{a}} \right)^{-2/3} + \phi_{P^{2}N} -
\phi_{P^{2}N} \left( \frac{f}{f_{a}} \right)^{-1/3} + \Phi
\end{eqnarray}
i.e.\
\begin{eqnarray}
\phi(f) & = & - \phi_{N}\left( \frac{f}{f_{a}} \right)^{-5/3}- \phi_{P^{1}N} \left( \frac{f}{f_{a}}
\right)^{-1} \nonumber \\
   & + & \phi_{P^{1.5}N} \left( \frac{f}{f_{a}} \right)^{-2/3} - \phi_{P^{2}N} \left( \frac{f}{f_{a}}
\right)^{-1/3} +
\phi_{c} + \Phi
\label{phioff2}
\end{eqnarray}

As $f \rightarrow \infty$, then $\phi(f) \rightarrow \phi_{c} + \Phi$. If $f=f_{a}$ then $\phi(f=f_{a}) =
\Phi$. If $f=f_{lso}$ then $\phi(f=f_{lso})$ is given by Eq(\ref{phioff2}). evaluated with $f=f_{lso}$.
 
\subsubsection*{Algorithm}

\subsubsection*{Uses}

For a fuller description of how this function is used in the generation of an inspiral waveform, see the
documentation for the function \texttt{TappRPNTdomFreq()}.
The nomenclature adopted is the same as that used in Sathyaprakash, PRD, 50, R7111, 1994, which may be
consulted for further details.

\subsubsection*{Notes}

See Sathyaprakash, PRD, 50, R7111, 1994, Eq.(5) or document for
module \texttt{TappRpnTdomFreq} for further details.

\vfill{\footnotesize\input{LALTappRpnTdomFreqPhaseCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALTAPPRPNTDOMFREQPHASEC, "$Id$");

/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCP"> */
void LALTappRpnTdomFreqPhase0PN (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params) 
{ /* </lalVerbatim>  */

  REAL8 f;

  INITSTATUS (status, "LALTappRpnTdomFreqPhase", LALTAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 * pow(f, fiveby3)
           + params->pc;

  RETURN(status);

}


/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCP"> */
void LALTappRpnTdomFreqPhase1PN (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params) 
{ /* </lalVerbatim>  */

  REAL8 f;

  INITSTATUS (status, "LALTappRpnTdomFreqPhase", LALTAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 / pow(f, fiveby3)
           + params->pc;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCP"> */
void LALTappRpnTdomFreqPhase2PN (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params) 
{ /* </lalVerbatim>  */

  REAL8 f;

  INITSTATUS (status, "LALTappRpnTdomFreqPhase", LALTAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = params->f;
  
  *phase = - params->p0 / pow(f, fiveby3)
           - params->p2 / f
           + params->pc;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCP"> */
void LALTappRpnTdomFreqPhase3PN (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params) 
{ /* </lalVerbatim>  */

  REAL8 f,f2,f4;

  INITSTATUS (status, "LALTappRpnTdomFreqPhase", LALTAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = pow(params->f,oneby3);
  f2 = f*f;
  f4 = f2*f2;

  *phase = - params->p0 / f4*f
           - params->p2 / f2*f
           + params->p3 / f2
           + params->pc;


  RETURN(status);

}


/*  <lalVerbatim file="LALTappRpnTdomFreqPhaseCP"> */
void LALTappRpnTdomFreqPhase4PN (LALStatus *status,
                                 REAL8 *phase, 
                                 InspiralPhasesInput *params) 
{ /* </lalVerbatim>  */

  REAL8 f,f2,f4;

  INITSTATUS (status, "LALTappRpnTdomFreqPhase", LALTAPPRPNTDOMFREQPHASEC);

  ASSERT(phase, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = pow(params->f,oneby3) ;
  f2 = f*f;
  f4 = f2*f2;

  *phase = - params->p0 / (f4*f)
           - params->p2 / (f2*f)
           + params->p3 / f2
           - params->p4 / f
           + params->pc;


  RETURN(status);

}
