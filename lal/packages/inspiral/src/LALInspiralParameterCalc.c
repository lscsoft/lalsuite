/*  <lalVerbatim file="LALInspiralParameterCalcCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralParameterCalc.c}}


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralParameterCalcCP}
\index{\verb&LALInspiralParameterCalc()&}

\subsubsection*{Description}

The code \texttt{LALInspiralParameterCalc.c} takes as its input one pair of parameters which may be chosen
from the following set of five: $(m_{1}, m_{2}, m, \eta, \mu)$, where $m_{1}$ and $m_{2}$ are the masses of
the two compact objects, $m=m_{1}+m_{2}$ is their combined mass, $\eta=m_{1}m_{2}/(m_{1}+m_{2})^{2}$ is the
symmetric mass ratio, and $\mu=m_{1}m_{2}/(m_{1}+m_{2})$.

Whichever pair of parameters is given to the function as an input, the function calculates the other three.
It also calculates the Newtonian chirp time $\tau_{0}$, the first post--Newtonian chirp time
$\tau_{2}$, the 1.5 order post--Newtonian chirp time $\tau_{3}$, the second order
post--Newtonian chirptime $\tau_{4}$ and the chirp
mass $\mathcal{M}$ which is defined as $\mathcal{M}=(\mu^{3} m^{2})^{1/5}$.

The chirp times are related to the masses of the stars and $f_{a}$ in the following way:
\begin{equation}
\tau_{0} = \frac{5}{256} \eta^{-1} m^{-5/3} (\pi f_{a})^{-8/3} \,,
\end{equation}

\begin{equation}
\tau_{2} = \frac{3715+4620 \eta}{64512 \eta m (\pi f_{a})^{2}} \,,
\end{equation}

\begin{equation}
\tau_{3} = \frac{\pi}{8 \eta m^{2/3} (\pi f_{a})^{5/3}}
\end{equation}

and
\begin{equation}
\tau_{4} = \frac{5}{128 \eta m^{1/3} (\pi f_{a})^{4/3}} \left[ \frac{3058673}{1016064} +
\frac{5429}{1008} \eta
+ \frac{617}{144} \eta^{2} \right] \,.
\end{equation}

These formulas show that an additional parameter $f_{a}$ is needed. This is the frequency at which the
detectors' noise curve rises steeply (the seismic limit).

\subsubsection*{Algorithm}


\subsubsection*{Uses}
None.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralParameterCalcCV}}

</lalLaTeX>  */



#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

NRCSID (LALINSPIRALPARAMETERCALCC, "$Id$");

/*  <lalVerbatim file="LALInspiralParameterCalcCP"> */
void 
LALInspiralParameterCalc (
   LALStatus *status, 
   InspiralTemplate *params)
{ /* </lalVerbatim> */

   REAL8 m1, m2, totalMass, eta, mu, piFl, etamin, tiny;
   REAL8 x1, x2, A0, A2, A3, A4, B2, B4, C4, dumm_const1, dumm_const2, dumm_const3;
   static REAL8 oneby4;
   void *pars;
   DFindRootIn rootIn;
   EtaTau02In Tau2In;
   EtaTau04In Tau4In;

   INITSTATUS (status, "LALInspiralParameterCalc", LALINSPIRALPARAMETERCALCC );
   ATTATCHSTATUSPTR(status);
 
   ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params->massChoice >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->massChoice <= 6, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   oneby4 = 1./4.;
   etamin = 1.e-10;
   tiny = 1.e-10;
   piFl = LAL_PI * params->fLower;
   switch(params->massChoice) {

      case m1Andm2:

         ASSERT(params->mass1 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mass2 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         m1 = params->mass1;
         m2 = params->mass2;
         params->totalMass = totalMass = m1+m2;
         params->eta = eta = m1*m2/pow(totalMass,2);
         if (params->eta > oneby4) params->eta -= tiny;
         params->mu = mu = m1*m2/totalMass;
         params->chirpMass = pow(mu,0.6)*pow(totalMass,0.4);

      break;

      case totalMassAndEta:

         ASSERT(params->totalMass > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->eta > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         if (params->eta > oneby4) params->eta -= tiny;
         ASSERT(params->eta <= oneby4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         totalMass = params->totalMass;
         eta = params->eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass + (totalMass*sqrt(oneby4-eta));
            params->mass2 = 0.5*totalMass - (totalMass*sqrt(oneby4-eta));
         }
         params->mu = eta*totalMass;
         params->chirpMass = pow(eta,0.6)*totalMass;

      break;

      case totalMassAndMu:

         ASSERT(params->totalMass > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mu > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mu < params->totalMass, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         totalMass = params->totalMass;
         mu = params->mu;
         eta =  (params->mu)/totalMass;
         if (eta > oneby4) eta -= tiny;
            params->eta = eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass + (totalMass*sqrt(oneby4-eta));
            params->mass2 = 0.5*totalMass - (totalMass*sqrt(oneby4-eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t02:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t2 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         A0 = 5./ pow(piFl, eightby3)/256.;
         A2 = 3715.0/(64512.0*pow(piFl,2.0));
         B2 = 4620.0/3715;
         Tau2In.t2 = params->t2;
         Tau2In.A2 = A2 * pow(params->t0/A0, 0.6);
         Tau2In.B2 = B2;
         pars = (void *) &Tau2In;
         rootIn.function = &LALEtaTau02;
         rootIn.xmax = oneby4+tiny;
         rootIn.xmin = etamin;
         rootIn.xacc = 1.e-8;
         LALEtaTau02(status->statusPtr, &x1, rootIn.xmax, pars);
         CHECKSTATUSPTR(status);
         LALEtaTau02(status->statusPtr, &x2, rootIn.xmin, pars);
         CHECKSTATUSPTR(status);
         if (x1*x2 > 0) {
            params->eta = 0.;
            DETATCHSTATUSPTR(status);
            RETURN(status);
         } else {
            LALDBisectionFindRoot(status->statusPtr, &eta, &rootIn, pars);
            CHECKSTATUSPTR(status);
         }
         if (eta > oneby4) eta-=tiny;
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass + (totalMass*sqrt(oneby4-eta));
            params->mass2 = 0.5*totalMass - (totalMass*sqrt(oneby4-eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t03:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t3 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         A0 = 5./ pow(piFl, eightby3)/256.;
         A3 = LAL_PI / pow(piFl, fiveby3)/8.;
         totalMass = A0 * params->t3/(A3 * params->t0);
         eta = A0/(params->t0 * pow(totalMass, fiveby3));
         if (eta > oneby4) eta-=tiny;
         params->eta = eta;
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass + (totalMass*sqrt(oneby4-eta));
            params->mass2 = 0.5*totalMass - (totalMass*sqrt(oneby4-eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;
 
      case t04:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t4 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         A0 = 5./(256. * pow(piFl, eightby3));
         A4 = 5./(128.0 * pow(piFl,fourby3)) * 3058673./1016064.;
         B4 = 5429./1008 * 1016064./3058673.;
         C4 = 617./144. * 1016064./3058673.;
         Tau4In.t4 = params->t4;
         Tau4In.A4 = A4 * pow(params->t0/A0, 0.2);
         Tau4In.B4 = B4;
         Tau4In.C4 = C4;
         pars = (void *) &Tau4In;
         rootIn.function = &LALEtaTau04;
         rootIn.xmax = oneby4+tiny;
         rootIn.xmin = etamin;
         rootIn.xacc = 1.e-8;
         LALEtaTau04(status->statusPtr, &x1, rootIn.xmax, pars);
         CHECKSTATUSPTR(status);
         LALEtaTau04(status->statusPtr, &x2, rootIn.xmin, pars);
         CHECKSTATUSPTR(status);
         if (x1*x2 > 0) {
            params->eta = 0.;
            DETATCHSTATUSPTR(status);
            RETURN(status);
         } else {
            LALDBisectionFindRoot(status->statusPtr, &eta, &rootIn, pars);
            CHECKSTATUSPTR(status);
         }
         if (eta > oneby4) eta-=tiny;
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass + (totalMass*sqrt(oneby4-eta));
            params->mass2 = 0.5*totalMass - (totalMass*sqrt(oneby4-eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;
   }
   if (params->eta > oneby4) params->eta-=tiny;
   totalMass = totalMass*LAL_MTSUN_SI;
   params->t0 = 5.0/(256.0*eta*pow(totalMass,fiveby3)*pow(piFl,eightby3));
   params->t2 = (3715.0 + (4620.0*eta))/(64512.0*eta*totalMass*pow(piFl,2.0));
   params->t3 = LAL_PI/(8.0*eta*pow(totalMass,twoby3)*pow(piFl,fiveby3));
   dumm_const1 = 3058673.0/1016064.0;
   dumm_const2 = (5429.0*eta)/1008.0;
   dumm_const3 = (617.0*pow(eta,2.0))/144.0;
   params->t4 = (5.0/(128.0*eta*pow(totalMass,oneby3)*pow(piFl,fourby3)))
                *(dumm_const1 + dumm_const2 + dumm_const3);

   switch (params->order) {
                        
      case newtonian:
      case oneHalfPN:
         params->t2=0.0;
         params->t3=0.0;
         params->t4=0.0;
         params->tC = params->t0;
      break;

      case onePN:
         params->t3=0.0;
         params->t4=0.0;
         params->tC = params->t0 + params->t2;
      break;

      case onePointFivePN:
         params->t4=0.0;
         params->tC = params->t0 + params->t2 - params->t3;
      break;

      case twoPN:
      default:
         params->tC = params->t0 + params->t2 - params->t3 + params->t4;
      break;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

