/*  <lalVerbatim file="LALInspiralSetupCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralSetup.c}}

Module to generate all the Taylor and Pade coefficients needed in 
waveform generation. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralSetupCP}
\index{\verb&LALInspiralSetup()&}

\subsubsection*{Description}

Module to generate all the coefficiants needed in the Taylor and Pade expansions of the energy and
flux functions $E^{\prime}(v)$ and $\mathcal{F}(v)$. 
These are used to solve the gravitational wave phasing formula,

\begin{eqnarray}
t(v) & = & t_{0} - m \int_{v_{0}}^{v} \,
\frac{E'(v)}{{\cal F}(v)} \, dv, \nonumber \\
\phi (v) & = & \phi_{0} - 2 \int_{v_{0}}^{v}  v^3 \,
\frac{E'(v)}{{\cal F}(v)} \, dv,
\end{eqnarray}

where $v=(\pi m F)^{1/3}$ is an invariantly defined velocity, $F$ is the instantaneous GW
frequency, and
$m$ is the total mass of the binary.

The expressions for $E^{\prime}(v)$ and $\mathcal{F}(v)$ may be written either in the form of a
Taylor
expansion or a P--Approximant expansion.
The Taylor expansion is only known up to second post--Newtonian order,

\begin{equation}
E^{\prime}(v) = - \eta v \left[1 - \frac{1}{6} (9+\eta) v^{2} -\frac{3}{8}(27-19 \eta +
\frac{\eta^{2}}{3}) v^{4} \right]
\end{equation}

The $n$--PN expression includes terms up to and including $v^{2n}$.
 
The P--Approximant expansion for $E^{\prime}(v)$ is given by

\begin{equation}
E^{\prime}_{P_{n}}(v) =\frac{\eta v}{\left[ 1 + E_{P_{n}} (x) \right] \sqrt{1 + e_{P_{n}} (x)}}
\frac{d e_{P_{n}} (x)}{dx}
\end{equation}
where
\begin{equation}
E_{P_{n}} (x) = \left[ 1 + 2 \eta ( \sqrt{1 + e_{P_{n}} (x)} -1) \right]^{1/2} -1
\end{equation}
and
\begin{equation}
e_{P_{2n}} (x) = -x P^{m}_{m + \epsilon} \left[ \sum^{n}_{k=0} a_{k} x^{k} \right]
\end{equation}
where $x=v^{2}$, and $ P^{m}_{m + \epsilon} \left[ \sum^{n}_{k=0} a_{k} x^{k} \right]$ is the Pade
Approximant of the Taylor series $\sum^{n}_{k=0} a_{k} x^{k}$.
An example of this is
\begin{equation}
e_{P_{4}}(x) = \frac{-x c_{0}}{1 + \frac{c_{1}x}{1 + c_{2}x}} = \frac{-c_{0} x (1 + c_{2}x)}
               {1 + (c_{1} + c_{2}) x}
\end{equation}


The Taylor series expansion of the flux function in the $\eta \neq 0$ case is known up to order
$v^{5}$,
\begin{equation}
\mathcal{F}(v) = \frac{32}{5} \eta^{2} v^{10} \left[ \sum^{5}_{k=0} A_{k} v^{k} \right]
\end{equation}
where
\begin{equation}
A_{0} = 1
\end{equation}
\begin{equation}
A_{1} = 0
\end{equation}
\begin{equation}
A_{2} = - = \frac{1247}{336} - \frac{35}{12} \eta
\end{equation}
\begin{equation}
A_{3} = 4 \pi
\end{equation}
\begin{equation}
A_{4} = - \frac{44711}{9072} + \frac{9271}{504} \eta + \frac{65}{18} \eta^{2}
\end{equation}
\begin{equation}
A_{5} = - \left( \frac{8191}{672} + \frac{535}{24} \eta \right) \pi
\end{equation}
If we now introduce the following "factored" flux function
\begin{equation}
f(v) = \left( 1 - \frac{v}{v_{pole} (\eta) } \right) \,\,  \mathcal{F}(v)
\end{equation}
which we write using Pade Approximants,
\begin{equation}
f_{P_{n}} (v) = \frac{32}{5} \eta^{2} v^{10} P^{m}_{m+\epsilon} \left[ \sum^{5}_{k=0} f_{k} v^{k}
\right]
\end{equation}
where the coefficients $f_{k}$ are defined in terms if the original coefficients $A_{k}$ of the
Taylor
series, $f_{k}=A_{k} - F_{k-1}/v_{pole}$.
Then we re--arrange to define a new $\mathcal{F}(v)$ given by
\begin{equation}
\mathcal{F}_{P_{n}}(v) = \frac{1}{\left( 1 - \frac{v}{v^{P_{n}}_{pole}(\eta)} \right)} \,\,
f_{P_{n}}(v)
\end{equation}
where $v^{P_{n}}_{pole}$ denotes the pole velocity defined by the $v^{n}$ approximant of $e(x)$.
 
This module sets up all the coefficients so that they may be used by the function
\texttt{ChooseModel} to define the energy and flux functions.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralSetupCV}}

</lalLaTeX>  */


   
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALSETUPC, "$Id$");

/*  <lalVerbatim file="LALInspiralSetupCP"> */
void LALInspiralSetup (LALStatus *status,
                       expnCoeffs *ak,
                       InspiralTemplate *params) 
{  /* </lalVerbatim>  */
   
   int ieta;
   REAL8 lso, eta, vpole;
   REAL8 a1, a2, a3, a4, a5;
   REAL8 c1, c2, c3, c4, c5;
   REAL8 oneby6=1.0/6.0;

   INITSTATUS (status, "LALInspiralSetup", LALINSPIRALSETUPC);

   ASSERT (ak,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->mass1 > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->mass2 > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fCutoff > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fCutoff > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 2*params->fCutoff, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   ak->ieta = params->ieta;
   ak->t0 = params->startTime;
   ak->m1 = params->mass1;
   ak->m2 = params->mass2;
   ak->f0 = params->fLower;
   ak->fn = params->fCutoff;
   ak->samplingrate = params->tSampling;
   ak->samplinginterval = 1./ak->samplingrate;

   lso = sqrt(oneby6);
/* 
   ieta determines the nature of the waveforms: 
   ieta=0 testmass waveforms 
   ieta=1 comparable mass waveforms.
*/
   ieta = ak->ieta;

/* Compute the total mass and eta from m1 and m2 */
   ak->totalmass = (ak->m1 + ak->m2);
   ak->eta = (ak->m1*ak->m2) / (ak->totalmass*ak->totalmass);
   eta = ak->eta;
   ak->totalmass = ak->totalmass * LAL_MTSUN_SI;

/* Set initial velocity according to initial frequency */

   ak->v0 = pow (LAL_PI * ak->totalmass * ak->f0, oneby3);

/* Taylor coefficients of E(x) */
   ak->ETaN = -eta/2.;
   ak->ETa1 = -(9.+ieta*eta)/12.;
   ak->ETa2 = -(27.-19*ieta*eta+ieta*eta*eta/3.)/8.;

/* Taylor coefficients of e(x) */
   ak->eTaN = -1.;
   ak->eTa1 = -1.-ieta*eta/3.;
   ak->eTa2 = -3.+35.*ieta*eta/12.;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
   ak->dETaN = -eta;
   ak->dETa1 = -(9. + ieta*eta)/6.;
   ak->dETa2 = -(3./8.) * (27. - 19.*ieta*eta + ieta*eta*eta/3.);

/* Pade coefficients of e(x)-function. */
   ak->ePaN = -1.;
   ak->ePa1 = 1.+ieta*eta/3.;
   ak->ePa2 = -(144. - 81.*ieta*eta + 4.*ieta*eta*eta) / (36.+12*ieta*eta);

/* Location of the 2PN T- and P-approximant last stable orbit and pole: */
   ak->vlsoT0 = lso;
   ak->vlsoP0 = lso;
   ak->vlsoP2 = lso;
/* 
   vlsoT2 =  6./(9.+ieta*eta); 
   This correct value makes vlso too large for vlsoT2 hence use 1/sqrt(6)
*/
   ak->vlsoT2 = lso;   
   ak->vlsoT4 = pow(-ak->ETa1 + pow(ak->ETa1*ak->ETa1 - 3*ak->ETa2,0.5)/(3*ak->ETa2), 0.5);
   ak->vlsoP4 = pow((-1.+pow(-ak->ePa1/ak->ePa2,0.5))/(ak->ePa1 + ak->ePa2), 0.5);
   ak->vpoleP4 = vpole = pow(4.*(3.+ieta*eta)/(36.-35.*ieta*eta), 0.5);

/* Taylor coefficients of flux. */
   ak->fTaN = ak->fPaN = ak->FTaN = 32.*eta*eta/5.;
   ak->FTa1 = 0.;
   ak->FTa2 = -1247./336.-35.*ieta*eta/12.;
   ak->FTa3 = 4.*LAL_PI;
   ak->FTa4 = -44711./9072.+9271.*ieta*eta/504.+65.*ieta*eta*eta/18.;
   ak->FTa5 = -(8191./672.+535./24.*ieta*eta)*LAL_PI;

/* Taylor coefficients of f(v)=(1-v/vpole)F(v) */

   ak->fTa1 = ak->FTa1 - 1./vpole;
   ak->fTa2 = ak->FTa2 - ak->FTa1/vpole;
   ak->fTa3 = ak->FTa3 - ak->FTa2/vpole;
   ak->fTa4 = ak->FTa4 - ak->FTa3/vpole;
   ak->fTa5 = ak->FTa5 - ak->FTa4/vpole;

/* Pade coefficients of f(v);  assumes that a0=1 => c0=1 */

   a1 = ak->fTa1;
   a2 = ak->fTa2;
   a3 = ak->fTa3;
   a4 = ak->fTa4;
   a5 = ak->fTa5;

   c1 = -a1;
   c2 = -(c1*c1 - a2)/c1;
   c3 = -(c1*pow(c2+c1,2.) + a3)/(c1*c2);
   c4 = -(c1*pow(c2+c1,3.) + c1*c2*c3*(c3+2*c2+2*c1) - a4)/(c1*c2*c3);
   c5 = -(c1*pow(pow(c1+c2,2.)+c2*c3,2.) + c1*c2*c3*pow(c1+c2+c3+c4,2.) + a5)/(c1*c2*c3*c4);

   ak->fPa1 = c1;
   ak->fPa2 = c2;
   ak->fPa3 = c3;
   ak->fPa4 = c4;
   ak->fPa5 = c5;

   RETURN (status);

}
