/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomTimePhase.c}}

The code \texttt{LALTappRpnTdomTimePhase.c} calculates the phase the waveform from an inspiralling binary
system as a function of time up to second post-Nowtonian order.

\subsubsection*{Prototypes}

\vspace{0.1in}

\input{LALTappRpnTdomTimePhaseCP}

\index{\verb&LALTappRpnTdomTimePhase()&}

\subsubsection*{Description}

The code \texttt{LALTappRpnTdomTimePhase.c} calculates the phase the waveform from an inspiralling binary
system as a function of time up to second post-Nowtonian order. The method used is as follows.

From Blanchet, Iyer, Will and Wiseman, CQG \textbf{13}, 575, 1996, the instantaneous orbital phase $\phi$
is given in terms of the dimensionless time variable $\Theta$  by

\begin{eqnarray}
\phi(t) = \phi_{c} - &  \frac{1}{\eta} \left\{ \Theta^{5/8} + \left( \frac{3715}{8064} + \frac{55}{96} \eta
\right) \Theta^{3/8} - \frac{3 \pi}{4} \Theta^{1/4} \right.  \\
                     &  + \left. \left( \frac{9275495}{14450688} + \frac{284875}{258048} \eta +
\frac{1855}{2048} \eta^{2} \right) \Theta^{1/8} \right\}
\label{phioft}
\end{eqnarray}
where $\phi_{c}$ is a constant which represents the value of the phase at instant $t_{c}$, which is the
instant of coalescence of the two point--masses which constitute the binary.
The dimensionless time variable $\Theta$ is given by

\begin{equation}
\Theta = \frac{c^{3}_{0} \eta}{5Gm} (t_{c} - t) \,\,.
\end{equation}

All of the equations presented so far have included explicitly their dependence upon $G$ and $c$. The code
uses units where $G=c=1$.


\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALTappRpnTdomTimePhaseCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALTAPPRPNTDOMTIMEPHASEC, "$Id$");

/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCP"> */
void LALTappRpnTdomTimePhase0PN (LALStatus *status,
                                 InspiralwavePhaseOutput *output,
                                 InspiralwavePhaseInput *params) 
{ /* </lalVerbatim>  */

  INITSTATUS (status, "LALTappRpnTdomTimePhase", LALTAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  output->phase = -(pow(params->td,0.625))/params->etaby2;

  RETURN(status);
}

/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCP"> */
void LALTappRpnTdomTimePhase1PN (LALStatus *status,
                                 InspiralwavePhaseOutput *output,
			         InspiralwavePhaseInput *params) 

{ /* </lalVerbatim>  */

  INITSTATUS (status, "LALTappRpnTdomTimePhase", LALTAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  output->phase = -(pow(params->td,0.625))/params->etaby2;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCP"> */
void LALTappRpnTdomTimePhase2PN (LALStatus *status,
                                 InspiralwavePhaseOutput *output,
			         InspiralwavePhaseInput *params) 

{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimePhase", LALTAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  output->phase = -((theta4*theta)+(params->a2*theta2*theta))/params->etaby2;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCP"> */
void LALTappRpnTdomTimePhase3PN (LALStatus *status,
                                 InspiralwavePhaseOutput *output,
			         InspiralwavePhaseInput *params) 

{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimePhase", LALTAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  output->phase = -((theta4*theta)+(params->a2*theta2*theta)-(params->a3*theta2))
                   /params->etaby2;


  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomTimePhaseCP"> */
void LALTappRpnTdomTimePhase4PN (LALStatus *status,
                                 InspiralwavePhaseOutput *output,
			         InspiralwavePhaseInput *params) 

{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimePhase", LALTAPPRPNTDOMTIMEPHASEC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  output->phase = -((theta4*theta)+(params->a2*theta2*theta)-(params->a3*theta2)+(params->a4*theta))
                   /params->etaby2;


  RETURN(status);
}

