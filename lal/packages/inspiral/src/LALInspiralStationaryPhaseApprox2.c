/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox2CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox2.c}}
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform.
 * %% A one-line description of the function(s) defined in this module.
 *
 * \subsubsection*{Prototypes}
 * \input{LALInspiralStationaryPhaseApprox2CP}
 * \idx{LALInspiralStationaryPhaseApprox2()}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
Consider a gravitational wave (GW) signal of the form,
\begin {equation}
h(t)=2a(t)\cos\phi(t)= a(t) \left [ e^{-i \phi(t)} + e^{i \phi(t)} \right ],
\end {equation}
where $\phi(t)$ is the phasing formula, either specified as an explicit
function of time or given implicitly by a set of differential equations
\cite{dis3}.  The quantity $2\pi F(t) = {d\phi(t)}/{dt}$ defines the instantaneous
GW frequency $F(t)$, and is assumed to be
continuously increasing. (We assume $F(t)>0$.)
 Now the Fourier transform $\tilde h(f)$ of $h(t)$ is defined as
\begin {equation}
\tilde{h}(f) \equiv \int_{-\infty}^{\infty} dt e^{2\pi ift} h(t)
            = \int_{-\infty}^{\infty}\,dt\, a(t) 
              \left[ e^{2\pi i f t - \phi(t)}  +  e^{2\pi ift +\phi(t)}\right ].
\end {equation}
The above transform can be computed in the stationary
phase approximation (SPA). For positive frequencies only the first term
on the right  contributes and yields the following {\it usual} SPA:
\begin {equation}
\tilde{h}^{\rm uspa}(f)= \frac {a(t_f)} {\sqrt {\dot{F}(t_f)}}
e^{ i\left[ \psi_f(t_f) -\pi/4\right]},\ \ 
\psi_f(t) \equiv  2 \pi f t -\phi(t), 
\label{eq:inspiralspa1}
\end {equation}
and $t_f$ is the saddle point defined by solving for $t$, $ d \psi_f(t)/d t = 0$,
i.e. the time $t_f$ when the GW frequency $F(t)$ becomes equal to the 
Fourier variable $f$. In the adiabatic approximation where 
the value of $t_f$ is given by the following integral:
\begin{equation}
t_f = t_{\rm ref} + m \int_{v_f}^{v_{\rm ref}} \frac{E'(v)}{{\cal F}(v)} dv,
\phi (v) = \phi_{\rm ref} + 2 \int_v^{v_{\rm ref}} dv v^3 \, \frac{E'(v)}{{\cal F}(v)},
\label{eq:InspiralTimeAndPhaseFuncs}
\end{equation}
where $v_{\rm ref}$ is a fiducial reference point that sets the origin of
time, $v_f \equiv (\pi m f)^{1/3},$ $E'(v)\equiv dE/dv$ is the derivative of
the binding energy of the system and ${\cal F}(v)$ is the gravitational wave
flux. 
Using $t_f$ and $\phi(t_f)$ in the above equation and 
using it in the expression for $\psi_f(t)$ we find
\begin{equation}
 \psi_f(t_f) = 2 \pi f t_{\rm ref} - \phi_{\rm ref} + 2\int_{v_f}^{v_{\rm ref}} 
(v_f^3 - v^3)
\frac{E'(v)}{{\cal {\cal F}}(v)} dv .
\label{eq:InspiralFourierPhase}
\end{equation}
This is the general form of the stationary phase approximation which
can be applied to {\it all} time-domain signals, including the P-approximant
and effective one-body waveforms. In some cases the Fourier domain phasing
can be worked out explicitly, which we now give:

Using PN expansions of energy and flux but
re-expanding the ratio $E'(v)/{\cal F}(v)$ in Eq.~(\ref{eq:InspiralFourierPhase}) one
can solve the integral explicitly. This leads to the following
explicit, Taylor-like, Fourier domain phasing formula:
\begin{equation}
 \psi_f(t_f) = 2 \pi f t_{\rm ref} - \phi_{\rm ref} + 
 \psi_N \sum_{k=0}^5 {\psi}_k (\pi m f)^{(k-5)/3} 
\label{eq:InspiralFourierPhase:f2}
\end{equation}
where the coefficients ${\psi}_k$ up to 2.5 post-Newtonian approximation are given by:
$$\psi_N =  \frac{3}{128\eta},\ \ \ \psi_0 = 1,\ \ \ \psi_1 = 0,\ \ \   
\psi_2 =  \frac{5}{9} \left ( \frac{743}{84} + 11\eta\right ),\ \ \ 
\psi_3 =  -16\pi,$$
$$\psi_4 =  \frac{9275495}{7225344}+\frac{284875\eta}{129024 } + \frac{1855\eta^2}{1024},$$
$$\psi_5 =  \frac{5}{3} \left ( \frac{7729}{252} + \eta \right ) \pi +
   \frac{8}{3} \left ( \frac{38645}{672} + \frac{15}{8} \eta \right ) 
	\ln \left ( \frac{v}{v_{\rm ref}} \right )\pi.$$
Eq.~(\ref{eq:InspiralFourierPhase:f2}) is (one of) the  standardly used frequency-domain phasing formulas.
This is what is implimented in {\tt LALInspiralStationaryPhaseApproximation2()}

Alternatively, substituting (without doing any re-expansion or re-summation) 
for the energy and flux functions their PN expansions
or the P-approximants of energy and flux functions 
and solving the integral in Eq.~(\ref{eq:InspiralFourierPhase}) numerically
one obtains the T-approximant SPA or P-approximant SPA, respectively.
However, just as in the time-domain, the frequency-domain phasing is 
most efficiently computed by a pair of coupled, non-linear, ODE's:
\begin{equation}
\frac{d\psi}{df} - 2\pi t = 0, \ \ \ \
\frac{dt}{df} + \frac{\pi m^2}{3v^2} \frac{E'(f)}{{\cal F}(f)} = 0,
\label {eq:frequencyDomainODE}
\end{equation}
rather  than by numerically computing the integral in  
Eqs.~(\ref{eq:InspiralFourierPhase}).
This is what is implimented in {\tt LALInspiralStationaryPhaseApproximation1()}

 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * The standard SPA is given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 * We define a variable function pointer {\tt LALInspiralTaylorF2Phasing} and point
 * it to one of the {\texttt static} functions defined within this function
 * that explicitly calculate the Fourier phase at the PN order chosen by the user.
 * The function returns the Frequency domain waveform in the convention of fftw.
 * Morever the reference points are chosen so that on inverse Fourier transforming
 * the time-domain will 
 * \begin{itemize}
 * \item be padded with zeroes in the first {\tt nStartPad} bins,
 * \item begin with a phase shift of {\tt nStartPhase} radians,
 * \item have an amplitude of $Nv_0^2$ at time $t_0$ (=0 when nStartPad=0)
	 where $v_0$ is the PN expansion parameter and $N$ is the number of 
	 data points.
 * \end{itemize}
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
   The derivative of the frequency that occurs in the
   amplitude of the Fourier transform, namely
   $1/\sqrt{\dot{F}(t)},$ in Eq.~(\ref{eq:inspiralspa1}), 
   is computed using
   \begin{eqnarray}
   \dot{F(t)} & = & \frac{dF}{dt} \\
              & = & \frac{dF}{dv}\frac{dv}{dE}\frac{dE}{dt}\\
              & = & \frac{3v^2}{\pi m}\frac{{-\cal F}(v)}{E'(v)},
   \end{eqnarray}
   where we have used the fact that the gravitational wave flux
   is related to the binding energy $E$ via energy balance equation
   ${\cal F} = -dE/dt$ and $F=v^3/m.$
   At the Newtonian order $E=-\eta m v^2/2,$ and ${\cal F} = 32\eta^2 v^{10}/5,$
   giving $\dot{F}(t(v)) = \frac{96\eta}{5\pi m^2} v^{11}.$ Taking
   $2a(t(v)) = v^2$ (i.e., $h(t) = v^2 \cos (\phi(t)),$ this gives, the 
   total amplitude of the Fourier transform to be 
   $$\frac{a(t(v))}{\sqrt{\dot{F}(t(v))}} =  \sqrt{\frac{5\pi m^2}{384\eta}} v_f^{-7/2}.$$
   This is the amplitude used in most of literature. However, including the
   full PN expansion of $1/\sqrt{\dot{F}(t)},$ gives a better agreement 
   (Damour, Iyer, Sathyaprakash, PRD 62, 084036 (2000)) between
   the time-domain and Fourier domains signals and this code therefore uses the full
   PN expansion.
     
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{LALInspiralStationaryPhaseApprox2CV}}
 *
 **** </lalLaTeX> */

#include "LALInspiral.h"

static void LALInspiralTaylorF2Phasing0PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing2PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing3PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing4PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing5PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing6PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);
static void LALInspiralTaylorF2Phasing7PN (REAL8 v, REAL8 *phase, expnCoeffs *ak);

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX2C, "$Id$");

/*  <lalVerbatim file="LALInspiralStationaryPhaseApprox2CP"> */
void 
LALInspiralStationaryPhaseApprox2 
(
 LALStatus *status,
 REAL4Vector *signal,
 InspiralTemplate *params) 
{ /* </lalVerbatim>  */
   REAL8 Oneby3, h1, h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;
   void (*LALInspiralTaylorF2Phasing)(REAL8 v, REAL8 *phase, expnCoeffs *ak);

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox2", LALINSPIRALSTATIONARYPHASEAPPROX2C);
   ATTATCHSTATUSPTR(status);
   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);

   switch (params->order) 
   {
	   case newtonian:
	   case oneHalfPN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing0PN;
		   break;
	   case onePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing2PN;
		   break;
	   case onePointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing3PN;
		   break;
	   case twoPN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing4PN;
		   break;
	   case twoPointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing5PN;
		   break;
	   case threePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing6PN;
		   break;
	   case threePointFivePN:
		   LALInspiralTaylorF2Phasing = LALInspiralTaylorF2Phasing7PN;
		   break;
   }
   n = signal->length;
   nby2 = n/2;
   memset( &ak, 0, sizeof( ak ) );
   LALInspiralSetup(status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   Oneby3 = 1.L/3.L;
   df = params->tSampling/signal->length;
   pimmc = LAL_PI * params->totalMass * LAL_MTSUN_SI;
   f0 = params->fLower;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn; 
   v = (pimmc*f0, Oneby3); 

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.
    */
   shft = 2.L*LAL_PI * (ak.tn + params->nStartPad/params->tSampling);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * df;
/* 
   Compute the standard stationary phase approximation. 
*/
   h1 = signal->data[0] = 0.L;
   h2 = signal->data[nby2] = 0.L;
   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn)
      {
	      /* 
	       * All frequency components below f0 and above fn are set to zero  
	       */
	      signal->data[i] = 0.;
	      signal->data[n-i] = 0.;
      }
      else
      {
	      v = pow(pimmc * f, Oneby3);
	      LALInspiralTaylorF2Phasing(v, &psif, &ak);
	      psi = shft * f + phi + psif;
	      /* 
		 dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
		 dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(pi m^2)
		 Note that our energy is defined as E=-eta v^2/2 and NOT -eta m v^2/2.
		 This is the reason why there is an extra m in the last equation above
		 amp = amp0 * pow(-dEnergybyFlux(v)/v^2, 0.5) * v^2;
		     = amp0 * pow(-dEnergybyFlux(v), 0.5) * v;
		     
	       */
	      amp = amp0 * pow(-func.dEnergy(v,&ak)/func.flux(v,&ak),0.5L) * v;
	      signal->data[i] = (REAL4) (amp * cos(psi));
	      signal->data[n-i] = (REAL4) (-amp * sin(psi));

      }

      /*
	 printf ("%e %e \n", v, psif);
	 printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
	 printf ("&\n");
       */
   
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing0PNCP"> */
static void LALInspiralTaylorF2Phasing0PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */

   REAL8 x;
   x = v*v;
   *phase = ak->pfaN/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing2PNCP"> */
static void LALInspiralTaylorF2Phasing2PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing3PNCP"> */
static void LALInspiralTaylorF2Phasing3PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing4PNCP"> */
static void LALInspiralTaylorF2Phasing4PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x;
   x = v*v;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * x*x)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing5PNCP"> */
static void LALInspiralTaylorF2Phasing5PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y;
   x = v*v;
   y = x*x;
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing6PNCP"> */
static void LALInspiralTaylorF2Phasing6PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y, z;
   x = v*v;
   y = x*x;
   z = y*x;
   /*
   fprintf(stderr, "This order is currently not implemented\n");
   exit(0);
   */
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v + ak->pfa6 * z)/pow(v,5.);
}

/*  <lalVerbatim file="LALInspiralTaylorF2Phasing7PNCP"> */
static void LALInspiralTaylorF2Phasing7PN (REAL8 v, REAL8 *phase, expnCoeffs *ak) {
/* </lalVerbatim>  */
   REAL8 x, y, z;
   x = v*v;
   y = x*x;
   z = y*x;
   /*
   fprintf(stderr, "This order is currently not implemented\n");
   exit(0);
   */
   *phase = ak->pfaN * (1. + ak->pfa2 * x + ak->pfa3 * v*x + ak->pfa4 * y
         + (ak->pfa5 + ak->pfl5 * log(v/ak->v0)) * y*v + ak->pfa6 * z + ak->pfa7*z*v)/pow(v,5.);
}
