/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomTimeFrequency.c}}

The code \texttt{LALTappRpnTdomTimeFrequency.c} calculates the frequency the 
waveform from an inspiralling binary system as a function of time up to second 
post-Nowtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomTimeFrequencyCP}
\index{\verb&LALTappRpnTdomTimeFrequency()&}

\subsubsection*{Description}

The code \texttt{LALTappRpnTdomTimeFrequency.c} calculates the frequency 
the waveform from an inspiralling binary system as a function of time up to 
second post-Nowtonian order. The method used is as follows.


The frequency of a gravitational wave is related to the parameter $x$, which 
is defined as
\begin{equation}
x(t) \equiv \left( \frac{Gm \omega(t)}{c^{3}_{0}} \right)^{2/3}
\end{equation}
in the following way
\begin{equation}
f_{GW}(t) = \frac{c^{3}_{0}}{G m \pi} \, x(t)^{3/2} \,\,.
\label{fofx}
\end{equation}
Now, $x(t)$ is given by
\begin{eqnarray}
x(t) & =  &  \frac{\Theta^{-1/4}}{4} \left\{  1 + \left(\frac{743}{4032} + 
          \frac{11}{48} \eta \right) \Theta^{-1/4} - 
          \frac{\pi}{5} \Theta^{-3/8} \right. \nonumber \\
     &  + & \left. \left( \frac{19583}{254016} + \frac{24401}{193536} \eta + 
            \frac{31}{288} \eta^{2} \right) \Theta^{-1/2} \right\}
\label{xoft}
\end{eqnarray}

All of these equations have included explicitly their dependence upon $G$ 
and $c$.  The code uses units where $G=c=1$.  Alternatively, one can work 
in terms of $f_{GW}(t)$ itself rather than first calculating $x(t)$ and then
finding $f_{GW}(t)$ from that. An equivalent expression derived in this way is
\begin{eqnarray}
\omega (t) & =  & \frac{c_{0}^{3}}{8 G m} \left\{ \Theta^{-3/8} + 
                  \left( \frac{743}{2688} + \frac{11}{32} \eta \right) \Theta^{-5/8} - 
                  \frac{3 \pi}{10} \Theta^{-3/4} \right. \nonumber \\
            & \times  & \left. \left( \frac{1\,855\,099}{14\,450\,688} + 
              \frac{56\,975}{258\,048} \eta + \frac{371}{2048} \eta^{2} \right) \Theta^{-7/8} \right\}
\end{eqnarray}
where $\omega$ here is the \emph{orbital} frequency $\omega_{orb}$. 
Then one simply uses $f_{GW}=2f_{orb}$ and $f_{orb}=\omega_{orb}/2 \pi$, 
and so $f_{GW}=\omega_{orb}/ \pi$.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}
See Blanchett, Iyer, Will and Wiseman, CQG, 13, 575, 1996 for further details.

\vfill{\footnotesize\input{LALTappRpnTdomTimeFrequencyCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALTAPPRPNTDOMTIMEFREQUENCYC, "$Id$");

/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCP"> */
void LALTappRpnTdomTimeFrequency0PN (LALStatus *status,
                                     InspiralwaveFrequencyOutput *output,
                                     InspiralwaveFrequencyInput *params) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimeFrequency", LALTAPPRPNTDOMTIMEFREQUENCYC);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  output->frequency = pow(params->td,-threeby8)/params->eightPiM;

  RETURN(status);
}

/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCP"> */
void LALTappRpnTdomTimeFrequency1PN (LALStatus *status,
                                     InspiralwaveFrequencyOutput *output,
			             InspiralwaveFrequencyInput *params) 

{ /* </lalVerbatim>  */


  INITSTATUS (status, "LALTappRpnTdomTimeFrequency", LALTAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  output->frequency = pow(params->td,-threeby8)/params->eightPiM;

  RETURN(status);
}

/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCP"> */
void LALTappRpnTdomTimeFrequency2PN (LALStatus *status,
                                     InspiralwaveFrequencyOutput *output,
			             InspiralwaveFrequencyInput *params) 

{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimeFrequency", LALTAPPRPNTDOMTIMEFREQUENCYC);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  output->frequency = ( 1./(theta2*theta) 
                      + params->a2/(theta4*theta) )
                      / params->eightPiM;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCP"> */
void LALTappRpnTdomTimeFrequency3PN (LALStatus *status,
                                     InspiralwaveFrequencyOutput *output,
			             InspiralwaveFrequencyInput *params) 

{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimeFrequency", LALTAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  output->frequency = ( 1./(theta2*theta) 
                      + params->a2/(theta4*theta)
                      - params->a3/(theta4*theta2) )
                      / params->eightPiM;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomTimeFrequencyCP"> */
void LALTappRpnTdomTimeFrequency4PN (LALStatus *status,
                                     InspiralwaveFrequencyOutput *output,
			             InspiralwaveFrequencyInput *params) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta4;

  INITSTATUS (status, "LALTappRpnTdomTimeFrequency", LALTAPPRPNTDOMTIMEFREQUENCYC);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params->td > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  theta = pow(params->td,0.125);
  theta2 = theta*theta;
  theta4 = theta2*theta2;

  output->frequency = ( 1./(theta2*theta) 
                      + params->a2/(theta4*theta)
                      - params->a3/(theta4*theta2)
                      + params->a4/(theta4*theta2*theta) )
                      / params->eightPiM;

  RETURN(status);
}
