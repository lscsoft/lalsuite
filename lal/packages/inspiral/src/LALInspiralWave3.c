/*  <lalVerbatim file="LALInspiralWave3CV">

Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave3.c}}

The code \texttt{LALInspiralWave3.c} generates a chirp waveform for a binary system consisting of two
non--spinning point--mass stars in quasi--circular orbits, up to second post--Newtonian order. The method
used is as follows.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave3CP}
\index{\verb&LALInspiralWave3()&}

\subsubsection*{Description}

The code \texttt{LALInspiralWave3.c} generates a chirp waveform for a binary system consisting of two
non--spinning point--mass stars in quasi--circular orbits, up to second post--Newtonian order. The method
used is as follows.

A gravitational wave detector is sensitive to the linear combination of the two polarizations $h_{+}$ and
$h_{\times}$ of the gravitational wave, in the following way,
\begin{equation}
h(t) = F_{+} h_{+} (t) + F_{\times} h_{\times} (t) \,\,,
\label{hoft}
\end{equation}
where $F_{+}$ and $F_{\times}$ are the \emph{beam--pattern functions} of the detector.
The two polarizations $h_{+}$ and $h_{\times}$ of the gravitational wave are given by the following
equation, which is taken from Blanchet, Iyer, Will and Wiseman, CQG \textbf{13}, 575, 1996, which we will
abbreviate as BIWW from now
on.

\begin{eqnarray}
h_{+,\times} &= & \frac{2Gm \eta}{c^{2}_{0} r} \left(\frac{Gm \omega}{c^{3}_{0}} \right)^{2/3} \left\{
\right. H^{(0)}_{+,\times} + x^{1/2} H^{(1/2)}_{+,\times}  \nonumber \\
             & +& x H^{(1)}_{+,\times} + x^{3/2} H^{(3/2)}_{+,\times} + x^{2} H^{(2)}_{+,\times} \left.
\right\}.
\label{H+cross}
\end{eqnarray}

In this Eq. $x$ is defined as $x \equiv (Gm \omega/c^{3}_{0})^{2/3}$ where $\omega$ is the 2PN--accurate
orbital frequency of the circular orbit, and $\omega=2 \pi/P$ where $P$ is the orbital period. $m \equiv
m_{1} + m_{2}$ is the total mass of the binary, and $c_{0}$ is the speed of light. $\eta$ is defined as
$\eta \equiv m_{1} m_{2}/m^{2}$.

In this code, we only include the $H^{(0)}_{+,\times}$ terms, which are given by (BIWW)

\begin{equation}
H^{(0)}_{+} = -(1+c^{2}) \cos 2\psi
\label{H0+}
\end{equation}
and
\begin{equation}
H^{(0)}_{\times} = -2c \sin 2 \psi
\label{H0cross}
\end{equation}

The notation used in these equations is as follows. The vector along the line of sight from the binary to
the detector defines the inclination angle $i$ with respect to the normal to the orbital plane of the
binary. In Eq.(\ref{H0+})
and Eq.(\ref{H0cross}) the variable c is defined as the cosine of the inclination angle $i$, $c=\cos i$.
The phase variable $\psi$ is defined by
\begin{equation}
\psi = \phi - \frac{2Gm \omega}{c^{3}_{0}} \ln \left( \frac{\omega}{\omega_{0}} \right)
\end{equation}
where $\phi$ is the actual orbital phase of the binary and $\omega_{0}$ is a constant which can be chosen
arbitrarily. In this code, we neglect the second term in this expression, leaving us with $\psi=\phi$.

We may also write the orbital angular velocity $\omega$ in terms of the orbital frequency $f_{orb}$ as
$\omega=2 \pi
f_{orb}$. The gravitational wave frequency $f_{GW}$ has a value which is twice that of the orbital
frequency, and so
$\omega= \pi f_{GW}$. If we substitute this into Eq.(\ref{H+cross}) then we obtain
\begin{equation}
h_{+} = A f_{GW}(t)^{2/3} (1+c^{2}) \cos 2\phi
\label{hplus}
\end{equation}
and
\begin{equation}
h_{\times} =  2A f_{GW}(t)^{2/3} 2c \sin 2 \phi
\label{htimes}
\end{equation}
where the constant $A$ is given by
\begin{equation}
A = \frac{- 2G^{5/3}m^{5/3} \eta \pi^{2/3}}{c^{4}_{0} r} \,\,.
\end{equation}
 
The two polarizations (\ref{hplus}) and (\ref{htimes}) are combined as given by Eq.(\ref{hoft}) to yield
\begin{equation}
h(t) = A f_{GW}(t)^{2/3} \left\{ F_{+} (1+c^{2}) \cos 2 \phi + 2F_{\times} c \sin 2 \phi \right\}
\end{equation}
This equation may be re--written in the form
\begin{equation}
h(t) = A \left[F^{2}_{+} (1+c^{2})^{2} + 4 F^{2}_{\times} c^{2} \right]^{1/2} \, f_{GW}(t)^{2/3} \, \cos
\left[ 2 \phi - \phi^{\prime} \right]
\end{equation}
where
\begin{equation}
\phi^{\prime} = \tan^{-1} \left( \frac{2 F_{\times} c}{F_{+} (1+c^{2})} \right)
\end{equation}
Therefore we arrive at an equation of the form
\begin{equation}
h(t) = A^{\prime} \, f_{GW}(t)^{2/3} \, \cos \left[ 2 \phi - \phi^{\prime} \right]
\label{hoft2}
\end{equation}
where
\begin{equation}
A^{\prime} = A \left[F^{2}_{+} (1+c^{2})^{2} + 4 F^{2}_{\times} c^{2} \right]^{1/2}
\end{equation}


In order to compute the time variation of $h(t)$ we need expressions for the time dependence of
$\psi(t)=\phi(t)$ and $f_{GW}(t)$. These are conveniently given in terms of the dimensionless time variable
 
\begin{equation}
\Theta = \frac{c^{3}_{0} \eta}{5Gm} (t_{c} - t)
\end{equation}
where $t_{c}$ is a constant which represents the instant of coalescence of the two point--masses which
constitute the binary. From BIWW we have the instantaneous orbital phase $\phi$ in terms of $\Theta$ as
given by
 
\begin{eqnarray}
\phi(t) & = & \phi_{c} - \frac{1}{\eta} \left\{ \Theta^{5/8} + \left( \frac{3715}{8064} + \frac{55}{96} \eta
\right) \Theta^{3/8} - \frac{3 \pi}{4} \Theta^{1/4} \right.  \nonumber \\
   &  + &  \left. \left( \frac{9275495}{14450688} + \frac{284875}{258048} \eta +
\frac{1855}{2048} \eta^{2} \right) \Theta^{1/8} \right\}
\label{TappRpnTdomTimephioft}
\end{eqnarray}
where $\phi_{c}$ is a constant representing the value of the phase at instant $t_{c}$.
 
 
We had the variable $x(t) \equiv (Gm \omega(t)/c^{3}_{0})^{2/3}$. This may be re--arranged to give
 
\begin{equation}
f_{GW}(t) = \frac{c^{3}_{0}}{G m \pi} \, x(t)^{3/2}
\label{TappRpnTdomTimefofx}
\end{equation}
where $x(t)$ is given by
\begin{eqnarray}
x(t) & =  &  \frac{\Theta^{-1/4}}{4} \left\{  1 + \left(\frac{743}{4032} + \frac{11}{48} \eta \right)
\Theta^{-1/4} - \frac{\pi}{5} \Theta^{-3/8} \right. \nonumber \\
     &  + & \left. \left( \frac{19583}{254016} + \frac{24401}{193536} \eta + \frac{31}{288} \eta^{2} \right)
\Theta^{-1/2} \right\}
\label{TappRpnTdomTimexoft}
\end{eqnarray}
 
All of the equations presented so far have included explicitly their dependence upon $G$ and $c$. The code
uses units where $G=c=1$.
 
To summarise, equations (\ref{TappRpnTdomTimexoft}) and (\ref{TappRpnTdomTimefofx}) are used to determine $f_{GW}(t)$ and
Eq.(\ref{TappRpnTdomTimephioft}) is used to determine $\phi(t)$. These quantities are then substituted into
Eq.(\ref{hoft2}) which defines $h(t)$.
 
\subsubsection*{Algorithm}

\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc} \\
\texttt{LALInspiralChooseModel} \\
\texttt{LALInspiralSetup} \\
\texttt{LALInspiralPhasing3} (via expnFunc)\\
\texttt{LALInspiralFrequency3}. (via expnFunc)\\


\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALWAVE3C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3CP"> */

void LALInspiralWave3 (
   LALStatus        *status,
   REAL4Vector      *output, 
   InspiralTemplate *params)

{ /* </lalVerbatim>  */

  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  expnFunc func;
  expnCoeffs ak;


  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  LALInspiralSetup (status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass >0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  tc=params->tC;       /* Instant of coalescence of the compact objects */
  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */


/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
   then the waveforn is truncated at f=flso.  
*/

  if (fu) 
     fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
     fHigh = ak.flso;

/* 
   Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum 
*/

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0; while (i<startShift) 
  {
      output->data[i] = 0.0;
      i++;
  }

  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = v*v; 
    output->data[i++] = (REAL4) (params->signalAmplitude * amp * cos(phase+phi0));
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    func.phasing3(status->statusPtr, &phase, td, &ak);
    CHECKSTATUSPTR(status);
    func.frequency3(status->statusPtr, &f, td, &ak);
    CHECKSTATUSPTR(status); 
  }
  params->fFinal = fOld;
  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
  while (i < (int)output->length) 
  {
      output->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
