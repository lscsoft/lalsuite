/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomFreqTofF.c}}

Module to calculate the LHS of Eq.(4) in the documentation for
\texttt{LALTappRpnTdomFreqTofF}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomFreqTofFCP}
\index{\texttt{LALTappRpnTdomFreqTofF()}}

\subsubsection*{Description}

The module \texttt{LALTappRpnTdomFreqToff} calculates the quantity which we may call toff, which is given by
the following equation:

\begin{eqnarray}
\mathrm{toff} = t - t_{a}  &  - \tau_{N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-8/3} \right] -
\tau_{P^{1}N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-2} \right] \\
    &  \\
          & + \tau_{P^{1.5}N} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] - \tau_{P^{2}N}
\left[ 1 - \left( \frac{f}{f_{a}} \right)^{-4/3} \right]  \,.
\label{toff}
\end{eqnarray}

The terms in this equation are defined as follows:
The parameter $t$ represents time, $\tau_{N}$ is the Newtonian chirp time, $\tau_{P^{1}N}$ is the first
post--Newtonian chirp time, and so on for $\tau_{P^{1.5}N}$ and $\tau_{P^{2}N}$. The parameter $f$ is the
instantaneous frequency
of the gravitational wave, and $f_{a}$ is the value of this frequency when the wave enters the lower end of
the detectors' bandwidth. Here we will neglect the $t_{a}$ term, which is the instant at which the wave
enters the lower end of the detectors bandwidth, and which can be defined alsewhere. If we define
\begin{equation}
\tau_{c} = \tau_{N} + \tau_{P^{1}N} - \tau_{P^{1.5}N} + \tau_{P^{2}N}
\end{equation}
then Eq.(\ref{toff}) becomes
\begin{eqnarray}
\mathrm{toff} &  = & t - t_{a}  - \tau_{N} + \tau_{N}\left( \frac{f}{f_{a}} \right)^{-8/3} - \tau_{P^{1}N} +
\tau_{P^{1}N} \left( \frac{f}{f_{a}} \right)^{-2} \nonumber \\
      & + & \tau_{P^{1.5}N} - \tau_{P^{1.5}N} \left( \frac{f}{f_{a}} \right)^{-5/3} - \tau_{P^{2}N} +
\tau_{P^{2}N} \left( \frac{f}{f_{a}} \right)^{-4/3}
\end{eqnarray}
i.e.\
\begin{eqnarray}
\mathrm{toff} &  = &  t - t_{a}  + \tau_{N} \left( \frac{f}{f_{a}} \right)^{-8/3} + \tau_{P^{1}N} \left(
\frac{f}{f_{a}} \right)^{-2} \nonumber \\
   & - & \tau_{P^{1.5}N} \left( \frac{f}{f_{a}} \right)^{-5/3} + \tau_{P^{2}N} \left( \frac{f}{f_{a}}
\right)^{-4/3} -
\tau_{c}
\label{toff2}
\end{eqnarray}


\subsubsection*{Algorithm}


\subsubsection*{Uses}

\subsubsection*{Notes}

All frequencies are expressed in units of fs, the frequency of the
waveform when it first enters the detectable part of the detector's
bandwidth.
For a fuller description of how this function is used in the generation of an inspiral waveform, see the
documentation for the function \texttt{LALTappRpnTdomFreq}.
The nomenclature adopted is the same as that used in Sathyaprakash, PRD, 50, R7111, 1994, which may be
consulted for
further details.

\vfill{\footnotesize\input{LALTappRpnTdomFreqTofFCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALTAPPRPNTDOMFREQTOFFC, "$Id$");

/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCP"> */   
void LALTappRpnTdomFreqTofF0PN(LALStatus *status,
                               REAL8 *toff,
                               REAL8 f,
                               void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;

  INITSTATUS (status, "LALTappRpnTdomFreqTofF", LALTAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  toffIn = (InspiralToffInput *) params;
  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        - toffIn->tc;

  RETURN(status);
}

/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCP"> */   
void LALTappRpnTdomFreqTofF1PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;

  INITSTATUS (status, "LALTappRpnTdomFreqToff", LALTAPPRPNTDOMFREQTOFFC);
  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  *toff = toffIn->t 
        + toffIn->t0 / pow(f, (eightby3)) 
        - toffIn->tc;

  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCP"> */   
void LALTappRpnTdomFreqTofF2PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 f2,f4,f8;

  INITSTATUS (status, "LALTappRpnTdomFreqToff", LALTAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = pow(f,oneby3);
  f2 = f*f;
  f4 = f2*f2;
  f8 = f4*f4;

  *toff = toffIn->t 
        + toffIn->t0 / f8 
        + toffIn->t2 / f4*f2
        - toffIn->tc;


  RETURN(status);
}


/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCP"> */   
void LALTappRpnTdomFreqTofF3PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 f2,f4,f8;

  INITSTATUS (status, "LALTappRpnTdomFreqToff", LALTAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  f = pow(f,oneby3);
  f2 = f*f;
  f4 = f2*f2;
  f8 = f4*f4;

  *toff = toffIn->t 
        + toffIn->t0 / f8 
        + toffIn->t2 / f4*f2
        - toffIn->t3 / f4*f 
        - toffIn->tc;


  RETURN(status);
}

   
/*  <lalVerbatim file="LALTappRpnTdomFreqTofFCP"> */   
void LALTappRpnTdomFreqTofF4PN (LALStatus *status,
                                REAL8 *toff,
			        REAL8 f,
                                void *params) 
{ /* </lalVerbatim>  */

  InspiralToffInput *toffIn;
  REAL8 f2,f4,f8;

  INITSTATUS (status, "LALTappRpnTdomFreqToff", LALTAPPRPNTDOMFREQTOFFC);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  f = pow(f,oneby3);
  f2 = f*f;
  f4 = f2*f2;
  f8 = f4*f4;

  *toff = toffIn->t 
        + toffIn->t0 / f8 
        + toffIn->t2 / (f4*f2)
        - toffIn->t3 / (f4*f) 
        + toffIn->t4 / f4
        - toffIn->tc;


  RETURN(status);
}
