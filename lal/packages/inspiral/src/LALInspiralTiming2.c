/*  <lalVerbatim file="LALInspiralTiming2CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTiming2.c}}

Module to calculate the LHS of Eq.(4) in the documentation for
\texttt{LALInspiralTiming2}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralTiming2CP}
\index{\verb&LALInspiralTiming2()&}

\subsubsection*{Description}

The module \texttt{LALInspiralTiming2} calculates the quantity which we 
may call toff, which is given by the following equation:

\begin{eqnarray}
\mathrm{toff} = t - t_{a}  &  - \tau_{0} 
                \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-8/3} \right] -
                \tau_{2} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-2} \right] \\
                & + \tau_{3} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] - 
               \tau_{4} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-4/3} \right]  \,.
\label{toff}
\end{eqnarray}

The terms in this equation are defined as follows: The parameter $t$ 
represents time, $\tau_{0}$ is the Newtonian chirp time, $\tau_{2}$ is the first
post--Newtonian chirp time, and so on for $\tau_{3}$ and $\tau_{4}$. 
The parameter $f$ is the instantaneous frequency of the gravitational wave, 
and $f_{a}$ is the value of this frequency when the wave enters the lower end of
the detectors' bandwidth. Here we will neglect the $t_{a}$ term, which is the 
instant at which the wave enters the lower end of the detectors bandwidth, 
and which can be defined alsewhere. If we define
\begin{equation}
\tau_{c} = \tau_{0} + \tau_{2} - \tau_{3} + \tau_{4}
\end{equation}
then Eq.(\ref{toff}) becomes
\begin{eqnarray}
\mathrm{toff} &  = & t - t_{a}  - \tau_{0} + \tau_{0}
\left( \frac{f}{f_{a}} \right)^{-8/3} - \tau_{2} +
\tau_{2} \left( \frac{f}{f_{a}} \right)^{-2} \nonumber \\
      & + & \tau_{3} - \tau_{3} \left( \frac{f}{f_{a}} \right)^{-5/3} - \tau_{4} +
\tau_{4} \left( \frac{f}{f_{a}} \right)^{-4/3}
\end{eqnarray}
i.e.\
\begin{eqnarray}
\mathrm{toff} &  = &  t - t_{a}  + \tau_{0} \left( \frac{f}{f_{a}} \right)^{-8/3} + 
\tau_{2} \left( \frac{f}{f_{a}} \right)^{-2} \nonumber \\
   & - & \tau_{3} \left( \frac{f}{f_{a}} \right)^{-5/3} + \tau_{4} \left( \frac{f}{f_{a}}
\right)^{-4/3} - \tau_{c}
\label{toff2}
\end{eqnarray}


\subsubsection*{Algorithm}


\subsubsection*{Uses}

\subsubsection*{Notes}

All frequencies are expressed in units of fs, the frequency of the
waveform when it first enters the detectable part of the detector's
bandwidth.  For a fuller description of how this function is used in 
the generation of an inspiral waveform, see the
documentation for the function \texttt{LALInspiralWave2}.
The nomenclature adopted is the same as that used in 
Sathyaprakash, PRD, 50, R7111, 1994, which may be
consulted for further details.

\vfill{\footnotesize\input{LALInspiralTiming2CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALTIMING2C, "$Id$");

/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_0PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v8 = pow(v,8.);

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8;

  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_2PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v8 = v2*v2*v2*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2);


  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_3PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v8 = v3*v3*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3);

  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_4PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4);

  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_5PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5);

  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_6PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v8 = v6*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(v)) * v6);

  RETURN(status);
}
/*  <lalVerbatim file="LALInspiralTiming2CP"> */   
void LALInspiralTiming2_7PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v7, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM*f, oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  v8 = v7*v;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1. 
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(v)) * v6
        + toffIn->t7 * v7);

  RETURN(status);
}
