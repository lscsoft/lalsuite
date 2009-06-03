/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*  <lalVerbatim file="LALInspiralParameterCalcCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralParameterCalc.c}}
Given a pair of masses (or other equivalent parameters) compute
related chirp parameters.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralParameterCalcCP}
\idx{LALInspiralParameterCalc()}
\begin{itemize}
\item\texttt{params:} Input/Output, given a pair of binary parameters and a lower
frequency cutoff, other equivalent parameters are computed by this function.
\end{itemize}

\subsubsection*{Description}

The code takes as its input {\tt params->fLower} in Hz and
a pair of masses (in units of $M_\odot$) or chirptimes (in seconds measured from {\tt params->fLower})
and computes all the other {\em mass} parameters in the {\tt params} structure.
Users choice of input pair of {\em masses} should be specified by appropriately setting
the variable {\tt params->massChoice} as described in the Table below:
\begin{table}[h]
\begin{center}
\caption{For a given {\tt params->massChoice} in column 1 the user should specify the
parameters as in column 2, in units as in column 3. Column 4 gives the conventional meaning
of the parameters. Chirp times are measured from a lower frequency cutoff given
in {\tt params->fLower.}}
\begin{tabular}{cccc}
\hline
{\tt params->massChoice} & User should set & in units & which means \\
\hline
{\tt m1Andm2}         & ({\tt mass1, mass2})   & $(M_\odot, M_\odot)$          & $(m_1,m_2)$ \\
{\tt totalMassAndEta} & ({\tt totalmass, eta}) & $(M_\odot, 0 < \eta \le 1/4)$ & $(m, \eta)$\\
{\tt totalMassAndMu}  & ({\tt totalmass, mu})  & $(M_\odot, M_\odot)$          & $(m, \mu)$ \\
{\tt t02}             & ({\tt t0, t2})         & (sec, sec) & $(\tau_0, \tau_2)$ \\
{\tt t03}             & ({\tt t0, t3})         & (sec, sec) & $(\tau_0, \tau_3)$ \\
{\tt t04}             & ({\tt t0, t4})         & (sec, sec) & $(\tau_0, \tau_4)$ \\
\hline
\end{tabular}
\end{center}
\end{table}

If \texttt{massChoice} is not set properly an error condition will occur and
the function is aborted with a return value 999.
In the above list $m_{1}$ and $m_{2}$ are the masses of
the two compact objects, $m=m_{1}+m_{2}$ is the total
mass, $\eta=m_{1}m_{2}/(m_{1}+m_{2})^{2}$ is the
symmetric mass ratio, $\mu=m_{1}m_{2}/(m_{1}+m_{2})$ is
the reduced mass and $\tau$'s are the chirptimes
defined in terms of $f_{a}$={\tt fLower} by:
\begin{eqnarray}
\tau_{0} = \frac{5}{256 \eta m^{5/3} (\pi f_{a})^{8/3}}, \ \ \
\tau_{2} = \frac{(3715 + 4620 \eta)}{64512 \eta m (\pi f_{a})^{2}}, \ \ \
\tau_{3} = \frac{\pi}{8 \eta m^{2/3} (\pi f_{a})^{5/3}}\nonumber \\
\tau_{4} = \frac{5}{128 \eta m^{1/3} (\pi f_{a})^{4/3}} \left[ \frac{3058673}{1016064} +
\frac{5429}{1008} \eta + \frac{617}{144} \eta^{2} \right],\ \ \
\tau_5 = \frac {5}{256\eta f_a}  \left (\frac {7729}{252} + \eta \right ).
\end{eqnarray}
%% Beyond 2.5 PN order, chirp times do not have an
%% explicit expression in terms of the masses and $f_a.$
Whichever pair of parameters is given to the function as an input, the function
calculates the rest.  Apart from the various masses and chirptimes the function
also calculates the chirp mass $\mathcal{M}=(\mu^{3} m^{2})^{1/5}$ and
the total chirp time $\tau_C$ consistent with the approximation chosen:
\begin{table}[h]
\begin{center}
\caption{$t_C$ will be set according to the PN order chosen in {\tt params->approximant.}}
\begin{tabular}{cccccc}
\hline
& {\tt Newtonian} & {\tt onePN} & {\tt onePointFivePN} & {\tt twoPN} & {\tt twoPointFivePN}\\
\hline
  $\tau_C$
& $\tau_0$
& $\tau_0 + \tau_2$
& $\tau_0 + \tau_2-\tau_3$
& $\tau_0 + \tau_2-\tau_3 + \tau_4$
& $\tau_0 + \tau_2-\tau_3 + \tau_4 - \tau_5$ \\
\hline
\end{tabular}
\end{center}
\end{table}

\subsubsection*{Algorithm}
Root finding by bisection method is used to solve for mass ratio $\eta$ when
chirptimes $(\tau_0,\, \tau_2)$ or $(\tau_0,\, \tau_4)$ is input.

\subsubsection*{Uses}
When appropriate this function calls:\\
\texttt{
LALDBisectionFindRoot\\
LALEtaTau02\\
LALEtaTau04\\
}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralParameterCalcCV}}

</lalLaTeX>  */



#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

NRCSID (LALINSPIRALPARAMETERCALCC, "$Id$");

/*  <lalVerbatim file="LALInspiralParameterCalcCP"> */
void
LALInspiralParameterCalc (
   LALStatus        *status,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   REAL8 m1, m2, totalMass, eta, mu, piFl, etamin, tiny, ieta;
   REAL8 x1, x2, A0, A2, A3, A4, B2, B4, C4,v,tN;
   REAL8 theta = -11831.L/9240.L;
   REAL8 lambda = -1987.L/3080.L;
   static REAL8 oneby4;
   void *pars;
   DFindRootIn rootIn;
   EtaTau02In Tau2In;
   EtaTau04In Tau4In;

   INITSTATUS (status, "LALInspiralParameterCalc", LALINSPIRALPARAMETERCALCC );
   ATTATCHSTATUSPTR(status);

   ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT((INT4)params->massChoice >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->massChoice <= 15, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   totalMass 	= 0.0;
   ieta 	= params->ieta;
   ieta 	= 1.;
   oneby4 	= 1./4.;
   etamin 	= 1.e-10;
   tiny 	= 1.e-10;
   piFl 	= LAL_PI * params->fLower;

   switch(params->massChoice)
   {
      case massesAndSpin:
      /*case spinOnly:*/
      case minmaxTotalMass:
      case m1Andm2:
      case fixedMasses:

         ASSERT(params->mass1 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->mass2 > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

         m1 = params->mass1;
         m2 = params->mass2;
         params->totalMass = totalMass = m1+m2;
         params->eta = eta = m1*m2/pow(totalMass,2);
         if (params->eta > oneby4) {
      		 params->eta -= tiny;
         }
         params->mu = mu = m1*m2/totalMass;
         params->chirpMass = pow(mu,0.6)*pow(totalMass,0.4);
         params->psi0 = 3./128./params->eta
	                * 1. * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-5./3.) ;
         params->psi3 = -3./128./params->eta
	                * (16 * LAL_PI) * pow((LAL_PI * params->totalMass * LAL_MTSUN_SI),-2./3.);

      break;

      case totalMassAndEta:
      case totalMassUAndEta:

         ASSERT(params->totalMass > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->eta > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

         if (params->eta > oneby4) {
		params->eta -= tiny;
   	}
         ASSERT(params->eta <= oneby4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

         totalMass = params->totalMass;
         eta = params->eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
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
         if (eta > oneby4) {
		 eta -= tiny;
	 }
            params->eta = eta;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t02:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t2 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

         A0 = 5./ pow(piFl, eightby3)/256.;
         A2 = 3715.0/(64512.0*pow(piFl,2.0));
         B2 = 4620.0/3715 * ieta;
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
         if (eta > oneby4) {
		 eta-=tiny;
   	}
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
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

	 if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;

      case t04:

         ASSERT(params->t0 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
         ASSERT(params->t4 > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

	 A0 = 5./(256. * pow(piFl, eightby3));
         A4 = 5./(128.0 * pow(piFl,fourby3)) * 3058673./1016064.;
         B4 = 5429./1008 * 1016064./3058673. * ieta;
         C4 = 617./144. * 1016064./3058673. * ieta;
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
         if (eta > oneby4) {
		 eta-=tiny;
	 }
         params->eta = eta;
         totalMass = pow(A0/(eta*params->t0), 0.6);
         totalMass = params->totalMass = totalMass/LAL_MTSUN_SI;
         if (eta <= oneby4) {
            params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
            params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
         }
         params->chirpMass = pow(eta,0.6)*totalMass;
         params->mu = eta*totalMass;

      break;


      case psi0Andpsi3:
      if (params->psi0 > 0 && params->psi3 < 0)
      {
	      params->totalMass = totalMass = -params->psi3/(16.L * LAL_PI * LAL_PI * params->psi0)/LAL_MTSUN_SI;
	      params->eta = eta = 3.L/(128.L * params->psi0 * pow (LAL_PI * totalMass*LAL_MTSUN_SI, fiveby3));

	      /* if eta < 1/4 amd M > 0 then physical values*/
	      if (eta <= oneby4)
	      {
		      params->mass1 = 0.5*totalMass * ( 1.L + sqrt(1.L - 4.L*eta));
		      params->mass2 = 0.5*totalMass * ( 1.L - sqrt(1.L - 4.L*eta));
		      params->mu = eta*totalMass;
		      params->chirpMass = pow(eta,0.6)*totalMass;
	      }
      }
      else
      {
	      params->eta = 0.;
	      DETATCHSTATUSPTR(status);
	      RETURN(status);
      }
      break;

     default:
      ABORT (status, 999, "Improper choice for massChoice in LALInspiralParameterCalc\n");
      break;
   }

   if (params->eta > oneby4) {
	   params->eta-=tiny;
	}
   totalMass 	= totalMass*LAL_MTSUN_SI;

   /* Should use the coefficients from LALInspiraSetup.c to avoid errors.
    * */
   v = pow(piFl * totalMass, 1.L/3.L);
   tN = 5.L/256.L / eta * totalMass / pow(v,8.L);

   params->t0 	= 5.0L/(256.0L*eta*pow(totalMass,fiveby3)*pow(piFl,eightby3));
   params->t2 	= (3715.0L + (4620.0L*ieta*eta))/(64512.0*eta*totalMass*pow(piFl,2.0));
   params->t3 	= LAL_PI/(8.0*eta*pow(totalMass,twoby3)*pow(piFl,fiveby3));
   params->t4 	= (5.0/(128.0*eta*pow(totalMass,oneby3)*pow(piFl,fourby3)))
              	* (3058673./1016064. + 5429.*ieta*eta/1008. +617.*ieta*eta*eta/144.);
   params->t5 	= -5.*(7729./252. - 13./3.*ieta*eta)/(256.*eta*params->fLower);
   /* This is a ddraft. t6 and t7 need to be checked propely*/
   params->t6 =  -10052469856691./23471078400. + 128./3.*LAL_PI*LAL_PI
     +(15335597827.L/15240960.L-451.L/12.L*LAL_PI*LAL_PI+352./3.*theta-2464.L/9.L*lambda)*ieta*eta
     +6848.L/105.L* LAL_GAMMA
     -15211.L/1728.L*ieta*eta*eta+25565.L/1296.L*eta*eta*eta*ieta;
   params->t6 = tN * (params->t6  + 6848.L/105.L*log(4.*v)) * pow(v,6);
   params->t7 = (-15419335.L/127008.L-75703.L/756.L*ieta*eta+14809.L/378.L*ieta*eta*eta) * LAL_PI * tN * pow(v,7);

   params->psi0 = 3.L/(128.L * eta * pow(LAL_PI * totalMass, fiveby3));
   params->psi3 = -3.L * LAL_PI/(8.L * eta * pow(LAL_PI * totalMass, twoby3));

   switch (params->order) {

      case newtonian:
      case oneHalfPN:
         params->t2=0.0;
/*       params->t3=0.0;*/
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0;
      break;

      case onePN:
         params->t3=0.0;
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2;
      break;

      case onePointFivePN:
         params->t4=0.0;
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3;
      break;

      case twoPN:
         params->t5=0.0;
         params->t6=0.0;
         params->t7=0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4;
      break;

      case twoPointFivePN:
         params->t6 = 0.0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5;

      case threePN:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6;

      case threePointFivePN:
      default:
         /*check the initialisation and then comment the next line. For now we
          * set t6=0 and t7=0*/
         params->t6 = 0;
         params->t7 = 0.0;
         params->tC = params->t0 + params->t2 - params->t3 + params->t4 - params->t5 + params->t6 - params->t7;
      break;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


