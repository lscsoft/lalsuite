/************************************ <lalVerbatim file="SimulatePulsarSignalCV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{SimulatePulsarSignal.c}}
\label{ss:SimulatePulsarSignal.c}

Routines to simulate pulsar-signals "exactly".
UNDER CONSTRUCTION. NOT USEABLE YET!

\subsubsection*{Prototypes}

\input{SimulatePulsarSignalCP}

\subsubsection*{Description}

The motivation for this module is to provide functions to
simulate pulsar signals \emph{with the best possible accuracy}, 
i.e. using no approximations, contrary to \verb+GeneratePulsarSignal+. 

Obviously this is not meant as a fast code to be used in a Monte-Carlo
simulation, but rather as a \emph{reference} to compare other (faster)
functions agains, in order to be able to gauge the quality of a given
signal-generation routine.

We want to calculate $h(t)$, given by
\begin{equation}
\label{eq:1}
	h(t) = F_+(t)\, h_+(t) + F_\times(t) \,h_\times(t)\,,
\end{equation}
where $F_+$ and $F_x$ are called the \emph{beam-pattern} functions, 
which depend of the wave polarization $\psi$,
the source position $\alpha$, $\delta$ and the detector position and
orientation ($\gamma$, $\lambda$, $L$ and $\xi$). The expressions for 
the beam-pattern functions are given in \cite{JKS98}, which we write as
\begin{eqnarray}
F_+(t) = \sin \zeta \cos 2\psi \, a(t)  + \sin \zeta \sin 2\psi \, b(t)\,,\\
F_\times(t) = \sin\zeta  \cos 2\psi \,b(t) - \sin\zeta \sin 2\psi \, a(t) \,,
\end{eqnarray}
where $\zeta$ is the angle between the interferometer arms, and 
\begin{eqnarray}
a(t) &=& a_1 \cos[ 2 (\alpha - T)) ] + a_2 \sin[ 2(\alpha - T)]
+ a_3 \cos[ \alpha - T ] + a_4 \sin [ \alpha - T ] + a_5\,,\\
b(t) &=& b_1 \cos[ 2(\alpha - T)] + b_2 \sin[ 2(\alpha - T) ]
+ b_3 \cos[ \alpha - T ] + b_4 \sin[ \alpha - T]\,,
\end{eqnarray}
where $T$ is the local (mean) sidereal time of the detector, and the 
time-independent coefficients $a_i$ and $b_i$ are given by
\begin{eqnarray}
a_1 &=& {1\over 16} \sin 2\gamma \,(3- \cos 2\lambda)\,(3 - \cos 2\delta)\,,\\
a_2 &=& -{1\over 4}\cos 2\gamma \,\sin \lambda \,(3 - \cos 2\delta) \,,\\
a_3 &=& {1\over 4} \sin 2\gamma \,\sin 2\lambda \,\sin 2\delta  \,\\
a_4 &=& -{1\over2} \cos 2\gamma \,\cos \lambda \,\sin 2 \delta\,,\\
a_5 &=& {3\over4} \sin 2\gamma \, \cos^2 \lambda \,\cos^2 \delta\,,
\end{eqnarray}
and 
\begin{eqnarray}
b_1 &=& \cos 2\gamma \,\sin \lambda \,\sin \delta\,,\\
b_2 &=& {1\over 4} \sin 2\gamma \,(3-\cos 2\lambda)\, \sin \delta\,,\\
b_3 &=& \cos 2\gamma \,\cos \lambda \,\cos\delta \,, \\
b_4 &=& {1\over 2} \sin2\gamma \,\sin 2\lambda \,\cos\delta\,,
\end{eqnarray}

The source model considered is a plane-wave
\begin{eqnarray}
h_+(t) &=& A_+\, \cos \Psi(t)\,,\\
h_\times(t) &=& A_\times \, \sin \Psi(t)\,,
\end{eqnarray}
where the wave-phase is $\Psi(t) = \Phi_0 + \Phi(t)$, and for an
isolated pulsar we have
\begin{equation}
\Phi(t) = 2\pi \left[\sum_{s=0} {f^{(s)}(\tau_\mathrm{ref}) \over
(k+1)!} \left( \tau(t) - \tau_\mathrm{ref} \right)^{k+1} \right]\,,
\end{equation}
where $\tau_\mathrm{ref}$ is the "reference time" for the definition
of the pulsar-parameters $f^{(s)}$ in the solar-system barycenter
(SSB), and $\tau(t)$ is the SSB-time of the phase arriving at the
detector at UTC-time $t$, which depends on the source-position
($\alpha$, $\delta$) and the detector-position, namely
\begin{equation}
  \tau (t) = t + { \vec{r}(t)\cdot\vec{n} \over c}\,,
\end{equation}
where $\vec{r}(t)$ is the vector from SSB to the detector, and $\vec{n}$ 
is the unit-vector pointing \emph{to} the source.

\subsubsection*{Algorithm}

This is a standalone "clean-room" implementation using no other 
outside-functions \emph{except} for \verb+LALGPStoLMST1()+ to calculate 
the local (mean) sidereal time at the detector for given GPS-time, 
(which I double-checked with an independent Mathematica script),
and and \verb+LALBarycenter()+ to calculate $\tau(t)$.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SimulatePulsarSignalCV}}

******************************************************* </lalLaTeX> */

#include "AVFactories.h"
#include <lal/GeneratePulsarSignal.h>

NRCSID( SIMULATEPULSARSIGNALC, "$Id$");

extern INT4 lalDebugLevel;

#define TRUE  (1==1)
#define FALSE (1==0)

#define oneBillion 1000000000L



/*--------------------------------------------------------------------------------
 * Simulate a pulsar signal to best accuracy possible
 *--------------------------------------------------------------------------------*/
/* <lalVerbatim file="SimulatePulsarSignalCP"> */
void
LALSimulatePulsarSignal (LALStatus *stat, 
			 REAL8TimeSeries **timeSeries, 
			 const PulsarSignalParams *params)
{ /* </lalVerbatim> */
  REAL8 LMST;		/* local mean sidereal time */
  LALPlaceAndGPS place_and_gps;
  LALMSTUnitsAndAcc units_and_acc;
  LALDate date;
  CHARVector *dateString = NULL;
  REAL8 a1, a2, a3, a4, a5;
  REAL8 b1, b2, b3, b4;
  REAL8 Zeta;
  REAL8 SinZetaCos2Psi, SinZetaSin2Psi;
  REAL8 gamma, alphaBi, alphaEast, xAzi, yAzi;
  LALFrDetector *det = &(params->site->frDetector);
  REAL8 lambda, delta, alpha;
  UINT4 i, Nsteps;	/* time-counter, number of steps */
  INT4 step_ns;	/* stepsize delta t in nanoseconds */
  LIGOTimeGPS t_i; /* time-step t_i */

  INITSTATUS( stat, "LALSimulatePulsarSignal", SIMULATEPULSARSIGNALC );
  ATTATCHSTATUSPTR(stat);

  /* orientation of detector arms */
  xAzi = det->xArmAzimuthRadians;
  yAzi = det->yArmAzimuthRadians;

  /* first calculate all quantities that don't depend on time */

  /* prefactor: angle between detector arms and polarization */
  Zeta =  xAzi - yAzi;
  if (Zeta < 0) Zeta = -Zeta;
  if(params->site->type == LALDETECTORTYPE_CYLBAR) Zeta = LAL_PI_2;

  SinZetaCos2Psi = sin(Zeta)*cos(2.0 * params->pulsar.psi);
  SinZetaSin2Psi = sin(Zeta)*sin(2.0 * params->pulsar.psi);

  /* get detector orientation gamma */
  gamma = atan2 ( cos(xAzi) + cos(yAzi), sin(xAzi)+sin(yAzi) );
  if (gamma < 0) gamma += LAL_TWOPI;	/* make sure it's positive (?do we need that?) */

  /*  printf ("\nDEBUG: gamma = %f deg\n", gamma * (REAL8)LAL_180_PI ); */

  /* the factors a_i and b_i */
  a1 = (1.0/16.0) * sin(2.0 * gamma) * (3.0 - cos(2.0*lambda)) * (3.0 - cos(2.0*delta));
  a2 = -(1.0/4.0) * cos(2.0 * gamma) * sin(lambda) * (3.0 - cos(2.0*delta));
  a3 = (1.0/4.0)  * sin(2.0 * gamma) * sin(2.0*lambda) * sin(2.0*delta);
  a4 = -(1.0/2.0) * cos(2.0 * gamma) * cos(lambda) * sin(2.0*delta);
  a5 = (3.0/4.0)  * sin(2.0 * gamma) * cos(lambda)*cos(lambda) * cos(delta)*cos(delta);


  b1 =            cos(2.0*gamma) * sin(lambda) * sin(delta);
  b2 = (1.0/4.0)* sin(2.0*gamma) * (3.0- cos(2.0*lambda)) * sin(delta);
  b3 =            cos(2.0*gamma) * cos(lambda) * cos(delta);
  b4 = (1.0/2.0)* sin(2.0*gamma) * sin(2.0*lambda) * cos(delta);


  /* get stepsize delta t */
  step_ns = (INT4)( (1.0/params->samplingRate) * oneBillion + 0.5 );	/* round to ns */

  Nsteps = (UINT4)( 1.0 * params->duration * params->samplingRate);
  t_i = params->startTimeGPS;

  place_and_gps.p_detector = params->site;
  units_and_acc.accuracy = LALLEAPSEC_STRICT;
  units_and_acc.units =   MST_RAD;	/* return LMST in radians */


  /* main loop: generate time-series */
  for (i=0; i < Nsteps; i++)
    {
      REAL8 AlphaMinusT;
      REAL8 cos2amT, sin2amT, cosamT, sinamT;

      place_and_gps.p_gps = &t_i;

      TRY (LALGPStoLMST1(stat->statusPtr, &LMST, &place_and_gps, &units_and_acc), stat);

      AlphaMinusT = params->pulsar.position.longitude - LMST;
      cos2amT = cos (2.0*AlphaMinusT);
      sin2amT = sin (2.0*AlphaMinusT);
      cosamT =  cos (AlphaMinusT);
      sinamT =  sin (AlphaMinusT);


    } /* for i < Nsteps */



  

  DETATCHSTATUSPTR(stat);
  RETURN(stat);

} /* LALSimulatePulsarSignal() */

