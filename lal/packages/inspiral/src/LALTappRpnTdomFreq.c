/*  <lalVerbatim file="LALTappRpnTdomFreqCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTappRpnTdomFreq.c}}

The module \texttt{LALTappRpnTdomFreq} generates a chirp waveform for a binary system consisting of two
non--spinning point--mass stars in quasi--circular orbits, up to second post--Newtonian order.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTappRpnTdomFreqCP}
\index{\verb&LALTappRpnTdomFreq()&}

\subsubsection*{Description}

The code generates the waveform by the following method:
The waveform may be written in the form
\begin{equation}
h(t) = A \left[ \pi f(t) \right]^{2/3} \cos \left[ \phi (t) \right] \,,
\label{waveform1}
\end{equation}
where $f(t)$ is the instantaneous gravitational wave frequency, $\phi(t)$ is the  phase, and the detailed
form of the amplitude term $A$ is a constant that depends on the detectors' antenna pattern, the masses of
the two stars and the distance to the binary.

We begin with $t$ as a function of $f$, which may be derived by from the knowledge that we may write the
time parameter $t$ as

\begin{equation}
t = t_{a} + \int^{t}_{t_{a}} dt^{\prime}
\end{equation}
where $t_{a}$ is the instant when the gravitational wave has a frequency $f_{a}$, which is equal to the
lowest frequency below which the detectors' noise rises steeply, making it difficult to extract signal
power below $f_{a}$. This equation may be re--written as

\begin{equation}
t = t_{a} + \int^{f}_{f_{a}} \frac{dt^{\prime}}{dE} \frac{dE}{df} df
\end{equation}
where the gravitational wave luminosity $F(f)$ is given by
\begin{equation}
F(f) = -\frac{dE}{dt^{\prime}} \,.
\end{equation}
Writing
\begin{equation}
\frac{dE}{df} = E^{\prime}(f)
\end{equation}
we then obtain
\begin{equation}
t - t_{a} = - \int^{f}_{f_{a}} \frac{E^{\prime}(f)}{F(f)} df \,.
\label{toff_int}
\end{equation}
In order to evaluate the integral we need to know the detailed form of $E^{\prime}(f)$ and $F(f)$. These
are given by

\begin{eqnarray}
E^{\prime} (f) & = & \eta m \left[ - \frac{1}{3} (\pi m)^{2/3} f^{-1/3} + \frac{1}{18} (9+\eta) (\pi m)^{4/3}
f^{1/3} \right. \nonumber \\
             &  + &   \left. \frac{1}{8} (27-19 \eta + \eta^{2}/3) (\pi m)^{2} f \right]
\label{Eoff}
\end{eqnarray}
and
\begin{eqnarray}
F(f) &  = & - \frac{32}{5} \eta^{2} (\pi m f)^{10/3} 
\left[ 1 + \left( \frac{- 1247}{336} - \frac{35}{12} \eta \right)
(\pi m f)^{2/3} + 4 \pi (\pi m f) \right. \nonumber \\
     &  + & \left. \left( -\frac{44711}{9072} + \frac{9271}{504} \eta + 
\frac{65}{18} \eta^{2} \right) (\pi m f)^{4/3}
\right] \,.
\label{Foff}
\end{eqnarray}
This leads us to
\begin{eqnarray}
t - t_{a} & = & 
\tau_{0} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-8/3} \right] + 
\tau_{2} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-2} \right] \nonumber \\
          & - & \tau_{3} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] + \tau_{4}
\left[ 1 - \left( \frac{f}{f_{a}} \right)^{-4/3} \right]
\label{toff1}
\end{eqnarray}
where $\tau_{0}$ is usually referred to as the Newtonian chirp time, $\tau_{2}$ is the first
post--Newtonian chirp time, and so on for $\tau_{3}$ and $\tau_{4}$. Because the $f/f_{a}$ term
appears throughout the above equations, the code works with the instantaneous frequency expressed in units
of $f_{a}$, $(f/f_{a})$ rather than in Hz. The chirp times are related to the masses of the stars and
$f_{a}$ in the following way:
\begin{equation}
\tau_{0} = \frac{5}{256} \eta^{-1} m^{-5/3} (\pi f_{a})^{-8/3} \,,
\end{equation}
 
\begin{equation}
\tau_{2} = \frac{3715+4620 \eta}{64512 \eta m (\pi f_{a})^{2}} \,,
\end{equation}
 
\begin{equation}
\tau_{3} = \frac{\pi}{8 \eta m^{2/3} (\pi f_{a})^{5/3}}
\end{equation}

end

\begin{equation}
\tau_{4} = \frac{5}{128 \eta m^{1/3} (\pi f_{a})^{4/3}} \left[ \frac{3058673}{1016064} +
\frac{5429}{1008} \eta
+ \frac{617}{144} \eta^{2} \right] \,.
\end{equation}

Eq.(\ref{waveform1}), however, requires $f$ as a function of $t$. In order to obtain this we re--arrange
Eq.(\ref{toff1})  into the form
 
\begin{eqnarray}
t - t_{a} & - & \tau_{0} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-8/3} \right] - \tau_{2} \left[ 1 -
\left( \frac{f}{f_{a}} \right)^{-2} \right] \nonumber \\
          & + & \tau_{3} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] - \tau_{4}
\left[ 1 - \left( \frac{f}{f_{a}} \right)^{-4/3} \right] = 0 \,.
\label{TappRpnTdomFreqtoff2}
\end{eqnarray}
The left--hand--side of Eq.(\ref{TappRpnTdomFreqtoff2}) is fed into a root--finding function,  which finds the value of
$(f/f_{a})$
which solves the equation for a given $t$, $t_{a}$, $\tau_{0}$, $\tau_{2}$, $\tau_{3}$ and
$\tau_{4}$.
Once we have obtained $(f/f_{a})(t)$ in this way we use it to calculate the phase $\phi(t)$.
 
In the same way that we derived the expression for the time formula, we may write the phase of the GW at
time $t$ as:\begin{equation}
\phi(t) = \Phi + 2 \pi \int^{t}_{t_{a}} f(t^{\prime}) dt^{\prime}
\end{equation}
where $\Phi$ is the phase of the wave at $t=t_{a}$. This equation may be re--written in the form
\begin{equation}
\phi(t) = \Phi + 2 \pi \int^{f}_{f_{a}} f \frac{dt}{dE} \frac{dE}{df} df
\end{equation}
or,
\begin{equation}
\phi(t) = \Phi - 2 \pi \int^{f}_{f_{a}} f \frac{E^{\prime}(f)}{F(f)} df \,.
\end{equation}
Once more, knowledge of the functions $E^{\prime}(f)$ (Eq.(\ref{Eoff})) and $F(f)$ (Eq.(\ref{Foff}))
enables us to obtain the time dependence of the phase $\phi(t)$,
 
\begin{eqnarray}
\phi(t) & = & \frac{16 \pi f_{a} \tau_{0}}{5} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-5/3} \right] + 4
\pi f_{a}\tau_{2} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1} \right]  
\nonumber \\
        & - & 5 \pi f_{a} \tau_{3} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-2/3} \right] + 8 \pi
f_{a} \tau_{4} \left[ 1 - \left( \frac{f}{f_{a}} \right)^{-1/3} \right] + \Phi \,.
\label{TappRpnTdomFreqphioff}
\end{eqnarray}
 
This part of the calculation is performed by the function \\ \texttt{TappRpnTdomFreqPhase}. The notation
$\phi(t)$ on the LHS of this equation is used because although the formula actually tells us
$\phi(f/f_{a})$, we have calculated
$(f/f_{a})(t)$ and used this value in the formula.
 
To summarise, we first of all find the value of $(f/f_{a})$ which solves Eq.(\ref{TappRpnTdomFreqtoff2}). This is
$(f/f_{a})(t)$. This value of $(f/f_{a})(t)$ is then substituted into Eq.(\ref{TappRpnTdomFreqphioff}) to give us the
phase $\phi(t) = \phi((f/f_{a})(t))$. Then both $(f/f_{a})(t)$ and $\phi(t)$ may be substituted into
Eq.(\ref{waveform1}) to calculate the waveform $h(t)$.
 
It must be noted that it is possible to expand the energy and flux functions in terms of the GW frequency
as we have
done here, or as a function of relative velocity $v$ of the stars. The connection between the two
formulations is
\begin{equation}
v = (\pi m f)^{1/3} \,.
\end{equation}

\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc}\\
\texttt{LALDBisectionFindRoot}\\
\texttt{LALTappRpnTdomFreqPhase}\\

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALTappRpnTdomFreqCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

NRCSID (LALTAPPRPNTDOMFREQC, "$Id$");

/*  <lalVerbatim file="LALTappRpnTdomFreqCP"> */
void LALTappRpnTdomFreq(LALStatus *status, 
                        REAL4Vector *output, 
                        InspiralTemplate *params)
{ /* </lalVerbatim>  */
  void (*LALTappRpnTdomFreqPhase) (LALStatus *, REAL8 *, InspiralPhasesInput *);
  void (*LALTappRpnTdomFreqTofF) (LALStatus *, REAL8 *, REAL8, void *);

  REAL8 dt, fs, fu, fsPi, fHigh, phase0;
  REAL8 phase, v, totalMass, fLso, freq;
  INT4 i, startShift, endShift, count;
  InspiralPhasesInput phaseIn;
  DFindRootIn rootIn;
  InspiralToffInput toffIn;
  void *funcParams;


  INITSTATUS (status, "LALTappRpnTdomFreq", LALTAPPRPNTDOMFREQC);
  ATTATCHSTATUSPTR(status);

  ASSERT(output,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(output->data,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  dt = 1.0/(params->tSampling);
  fs = params->fLower;
  fu = params->fCutoff;
  fsPi = fs*LAL_PI;
  startShift = params->nStartPad;
  endShift = params->nEndPad;
  phase0 = params->startPhase;


  switch (params->order) {
     case newtonian:
     case oneHalfPN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase0PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF0PN;
          rootIn.function = LALTappRpnTdomFreqTofF0PN;
          break;
     case onePN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase2PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF2PN;
          rootIn.function = LALTappRpnTdomFreqTofF2PN;
          break;
     case onePointFivePN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase3PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF3PN;
          rootIn.function = LALTappRpnTdomFreqTofF3PN;
          break;
     case twoPN:
          LALTappRpnTdomFreqPhase = &LALTappRpnTdomFreqPhase4PN;
          LALTappRpnTdomFreqTofF = &LALTappRpnTdomFreqTofF4PN;
          rootIn.function = LALTappRpnTdomFreqTofF4PN;
          break;
     default:
          fprintf(stderr, "LALTappRpnTdomFreq: No order selected ... exiting\n");
          exit(0);
     }

/* Calculate the three unknown paramaters from (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);


  ASSERT(params->totalMass > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->totalMass < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta <=0.25, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->mu >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);



  toffIn.t0 = params->t0;
  toffIn.t2 = params->t2;
  toffIn.t3 = params->t3;
  toffIn.t4 = params->t4;
  toffIn.tc = params->tC;

  phaseIn.p0 = 3.2 * fsPi * params->t0;
  phaseIn.p2 = 4.0 * fsPi * params->t2;
  phaseIn.p3 = 5.0 * fsPi * params->t3;
  phaseIn.p4 = 8.0 * fsPi * params->t4;
  phaseIn.pc = phaseIn.p0 + phaseIn.p2 - phaseIn.p3 + phaseIn.p4;


  totalMass = params->totalMass*LAL_MTSUN_SI;

  fLso = 1.0/(LAL_PI*totalMass*pow(6.0, 1.5));

  if (fu) 
  fHigh = (fu < fLso) ? fu/fs : fLso/fs; 
  else 
  fHigh = fLso/fs;

  ASSERT(fHigh*fs < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh*fs > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  toffIn.t = 0.0;

  rootIn.xmax = 1.1*fu/fs;
  rootIn.xacc = 1.0e-8;
  rootIn.xmin = 0.999999;


  i=0;
  while (i<startShift) output->data[i++] = 0.0;


/* Now cast the input structure to argument 4 of BisectionFindRoot so that it 
  of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;

  LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
  CHECKSTATUSPTR(status);

  
  phaseIn.f=freq;
  count=1;

  while (freq < fHigh) {
    v = pow(freq*fs*LAL_PI*totalMass, oneby3);
    LALTappRpnTdomFreqPhase(status->statusPtr, &phase, &phaseIn);
    CHECKSTATUSPTR(status);
    output->data[i++] = (REAL4) params->signalAmplitude * v*v * cos(phase+phase0);
    toffIn.t=count*dt;
    ++count;
    funcParams = (void *) &toffIn;
    LALDBisectionFindRoot(status->statusPtr, &freq, &rootIn, funcParams);
    CHECKSTATUSPTR(status);
    phaseIn.f=freq;
  }

  while (i<(int)output->length) output->data[i++]=0.0;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
