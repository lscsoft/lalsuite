/*  <lalVerbatim file="LALInspiralFrequency3CV">

Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralFrequency3.c}}

The code \texttt{LALInspiralFrequency3.c} calculates the frequency the 
waveform from an inspiralling binary system as a function of time up to second 
post-Nowtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralFrequency3CP}
\index{\verb&LALInspiralFrequency3()&}

\subsubsection*{Description}

The code \texttt{LALInspiralFrequency3.c} calculates the frequency 
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

\vfill{\footnotesize\input{LALInspiralFrequency3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALFREQUENCY3C, "$Id$");

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_0PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta3;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta3 = theta*theta*theta;

  *frequency = theta3*ak->ftaN;

  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_2PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2);
                      

  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_3PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_4PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

  ASSERT(ak, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  theta = pow(td,-0.125);
  theta2 = theta*theta;
  theta3 = theta2*theta;
  theta4 = theta3*theta;

  *frequency = theta3*ak->ftaN * (1.
             + ak->fta2*theta2
             + ak->fta3*theta3
             + ak->fta4*theta4);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_5PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

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
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_6PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

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
             + ak->fta6*theta6);
  RETURN(status);
}

/*  <lalVerbatim file="LALInspiralFrequency3CP"> */

void LALInspiralFrequency3_7PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak) 
{ /* </lalVerbatim>  */

  REAL8 theta,theta2,theta3,theta4,theta5,theta6,theta7;

  INITSTATUS (status, "LALInspiralFrequency3", LALINSPIRALFREQUENCY3C);

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
             + ak->fta6*theta6
             + ak->fta7*theta7);
  RETURN(status);
}
