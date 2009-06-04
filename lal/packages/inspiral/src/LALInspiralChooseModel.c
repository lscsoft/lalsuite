/*
*  Copyright (C) 2007 David Churches, Jolien Creighton, David McKechan, B.S. Sathyaprakash, Thomas Cokelaer, Duncan Brown
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

/*  <lalVerbatim file="LALInspiralChooseModelCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralChooseModel.c}}

Module to set the pointers to the required energy and flux functions.
Normally, a user is not required to call this function to generate a waveform.
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralChooseModelCP}
\index{\verb&LALInspiralChooseModel()&}
\begin{itemize}
\item {\tt f:} Output containing the pointers to the appropriate
energy, flux, frequency, timing and phasing functions.
\item {\tt ak:} Output containing the PN expnasion coefficients.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}
This module gives the post-Newtonian expansions and/or P-approximants
to the energy, its derivative and gravitational-wave flux functions. More
specifically, the {\tt static REAL8} functions below give Taylor expansions
of $dE/dv,$ and ${\cal F}(v),$ P-approximants of $e(v),$ $dE/dv$
(derived from $e(v)$) and ${\cal F}(v).$

{\tt LALInspiralChooseModel}
is used to set pointers to the required energy and flux functions
$E^{\prime}_T(v),$ $\mathcal{F}_T(v),$ $E^{\prime}_P(v)$ and $\mathcal{F}_P(v),$
in {\tt expnFunc,} as also the GW phasing and frequency fucntions used in
the various approximants to generate the waveform.
More specifically pointers are set to the following functions in the structure
{\tt expnFunc}:
\begin{itemize}
  \item {\tt EnergyFunction *dEnergy}
  \item {\tt FluxFunction *flux}
  \item {\tt InspiralTiming2 *timing2}
  \item {\tt InspiralPhasing2 *phasing2}
  \item {\tt InspiralPhasing3 *phasing3}
  \item {\tt InspiralFrequency3 *frequency3}
\end{itemize}
{\tt LALInspiralChooseModel} also outputs in {\tt ak} the
last stable orbit (LSO) velocity $v_{\rm LSO}$ (as {\tt ak->vn})
defined by the equation $E'(v_{\rm LSO})=0,$
the values of the GW frequency $f_{\rm LSO}=v_{\rm LSO}^3/(\pi m)$
(as {\tt ak->fn}) and time (as {\tt ak->tn}) elapsed from {\tt params->fLower}
to smaller of {\tt fCutOff} and {\tt ak->fn} by evaluating the integral
\begin{equation}
t_n = t_{0} - m \int^{v_n}_{v_0} \frac{E^{\prime}(v)}{\mathcal{F}(v)} \, dv\,,
\end{equation}
where $t_{0}$ (usually equal to zero) is the user specified starting
time for the waveform when the wave frequency reaches {\tt params->fLower}
and $v_{0}= (\pi m f)^{1/3}$ (with $f={\tt params->fLower}$) is the  velocity
at time $t_{0}.$  Note that $E'(v)$ and ${\cal F}(v)$ are defined in
{\tt f->dEnergy} and {\tt f->flux.}

\subsubsection*{Algorithm}
Numerical integration is used to compute {\tt ak->tn.}
\subsubsection*{Uses}
LALInspiralTofV

\subsubsection*{Notes}
\begin{itemize}
\item See Damour, Iyer and Sathyaprakash, PRD 57, 885, 1998 for further details.
Damour, Iyer and Sathyaprakash, PRD 63, 044023, 2001 is a resource paper that
summarizes how to generate waveforms in different approximations to the dynamics
of a compact binary under radiation reaction.
\item The Pade Approximant for the 1PN expansion is undefined as also
EOB at orders less than 2PN. BCV is independent of the PN order.
Spinning waveforms are only defined at the highest PN order.
\end{itemize}


\vfill{\footnotesize\input{LALInspiralChooseModelCV}}

</lalLaTeX>  */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALCHOOSEMODELC, "$Id$");

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEt0(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 dEnergy;
   dEnergy = ak->dETaN * v;
   return (dEnergy);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEt2(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x);
   return (dEnergy);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEt4(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x);
   return (dEnergy);
}


/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEt6(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 dEnergy, x;
   x = v*v;
   dEnergy = ak->dETaN * v * (1. + ak->dETa1*x + ak->dETa2*x*x + ak->dETa3*x*x*x);
   return (dEnergy);
}





/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft0(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10;
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft2(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2);
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft3(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v);
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft4(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1. + ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4);
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft5(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v);
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft6(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6);
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Ft7(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->FTaN * v10 * (1.+ ak->FTa2*v2 + ak->FTa3*v2*v + ak->FTa4*v4
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(v))*v6 + ak->FTa7*v6*v);
   return (flux);
}


#if 0 /* NOT USED */
/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 ep0(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 x, energy;
   ak = NULL;
   x = v*v;
   energy = -x;
   return (energy);
}
#endif

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 ep2(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1 * x);
   return (energy);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 ep4(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x));
   return (energy);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 ep6(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 x, energy;
   x = v*v;
   energy = -x / (1. + ak->ePa1*x /(1. + ak->ePa2*x /(1. + ak->ePa3*x)));
   return (energy);
}

#if 0 /* NOT USED */
/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEp0(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 energy, denergy, Energy, dEnergy, y;
   energy = ep0(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1;
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}
#endif

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEp2(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep2(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = -1. / ((1. + ak->ePa1*x)*(1. + ak->ePa1*x));
   dEnergy = v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEp4(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep4(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = (1. + 2.*ak->ePa2*x + ((ak->ePa1 + ak->ePa2) * ak->ePa2 * x*x))/pow(1. + (ak->ePa1 + ak->ePa2) * x ,2.);
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}


/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 dEp6(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 energy, denergy, Energy, dEnergy, x, y;
   x = v*v;
   energy = ep6(v, ak);
   y = sqrt(1.+energy);
   Energy = sqrt(1. + 2.* ak->eta * (y - 1.)) - 1.;
   denergy = (1. + 2.*(ak->ePa2+ak->ePa3)*x + (ak->ePa1*ak->ePa2
           + ak->ePa2*ak->ePa2 + 2.* ak->ePa2*ak->ePa3
           + ak->ePa3*ak->ePa3) * x*x)
           /pow(1. + (ak->ePa1 + ak->ePa2 + ak->ePa3) * x
           + ak->ePa1*ak->ePa3*x*x,2.);
   dEnergy = - v * ak->eta * denergy /((1.+Energy) * y);
   return(dEnergy);
}



#if 0 /* NOT USED */
/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp0(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10;
   return (flux);
}
#endif

#if 0 /* NOT USED */
/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp1(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v) * (1.-v/ak->vpoleP4));
   return (flux);
}
#endif

#if 0 /* NOT USED */
/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp2(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v / (1.+ak->fPa2*v)) * (1.-v/ak->vpoleP4));
   return (flux);
}
#endif

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp3(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v)))
	* (1.-v/ak->vpoleP4));
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp4(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v)))) * (1.-v/ak->vpoleP4));
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp5(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
	/ (1.+ak->fPa4*v / (1.+ak->fPa5*v))))) * (1.-v/ak->vpoleP4));
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp6(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v))))))
        * (1.-v/ak->vpoleP6));
   /* */
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp7(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v6,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v)))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * ak->FTl6*v6) ;
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp8(REAL8 v, expnCoeffs *ak)
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v6,v8,v10, l6, l8;
   v2 = v*v;
   v4 = v2*v2;
   v6 = v4*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   l6 = ak->FTl6;
   l8 = ak->FTl8 - ak->FTa2*ak->FTl6;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v
	/ (1.+ak->fPa8*v))))))))
        * (1.-v/ak->vpoleP6));
   flux *= (1.+  log(v/ak->vlsoP4) * (l6*v6 + l8*v8) ) ;
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
void
LALInspiralChooseModel(
   LALStatus        *status,
   expnFunc         *f,
   expnCoeffs       *ak,
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */

   REAL8 vn, vlso;
   TofVIn in1;
   REAL8 tofv;
   void *in2;

   INITSTATUS (status, "LALInspiralChooseModel", LALINSPIRALCHOOSEMODELC);
   ATTATCHSTATUSPTR(status);

   ASSERT (f,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (ak,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->order != oneHalfPN,status,LALINSPIRALH_ENULL,LALINSPIRALH_MSGENULL);
   ASSERT((INT4)params->order >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT((INT4)params->order <= 8, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   vlso = 0;

   switch (params->order)
   {
      case newtonian:
      switch (params->approximant)
      {
         case AmpCorPPN:
         case Eccentricity:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
         case IMRPhenomA:
            ak->vn = ak->vlso = vlso = ak->vlsoT0;
            f->dEnergy = dEt0;
            f->flux = Ft0;
            f->phasing2 = &LALInspiralPhasing2_0PN;
            f->timing2 = &LALInspiralTiming2_0PN;
            f->phasing3 = &LALInspiralPhasing3_0PN;
            f->frequency3 = &LALInspiralFrequency3_0PN;
            break;
         case PadeT1:
         case PadeF1:
         case EOB:
         case EOBNR:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case oneHalfPN:
        ABORT(status, LALINSPIRALH_ECHOICE, "OneHalfPN is not valid");
        break;
      case onePN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
         case IMRPhenomA:

            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft2;
            f->phasing2 = &LALInspiralPhasing2_2PN;
            f->timing2 = &LALInspiralTiming2_2PN;
            f->phasing3 = &LALInspiralPhasing3_2PN;
            f->frequency3 = &LALInspiralFrequency3_2PN;
            break;
         case PadeT1:
         case PadeF1:
         case EOB:
         case EOBNR:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case onePointFivePN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
         case IMRPhenomA:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft3;
            f->phasing3 = &LALInspiralPhasing3_3PN;
            f->frequency3 = &LALInspiralFrequency3_3PN;
            f->phasing2 = &LALInspiralPhasing2_3PN;
            f->timing2 = &LALInspiralTiming2_3PN;
            break;
         case PadeT1:
            ak->vn = ak->vlso = vlso = ak->vlsoP0;
            f->dEnergy = dEp2;
            f->flux = Fp3;
            break;
         case PadeF1:
         case EOB:
         case EOBNR:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case twoPN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
/*
   The value vlsoT4 is too large and doesn't work sometimes;
   so we use vlsoT2.
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt4;
            f->flux = Ft4;
            f->phasing2 = &LALInspiralPhasing2_4PN;
            f->timing2 = &LALInspiralTiming2_4PN;
            f->phasing3 = &LALInspiralPhasing3_4PN;
            f->frequency3 = &LALInspiralFrequency3_4PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case IMRPhenomA:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp4;
            break;
         case PadeF1:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case twoPointFivePN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
/*
   The value vlsoT4 is too large and doesn't work with 2.5 PN
   Taylor approximant; so we use vlsoT2.
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt4;
            f->flux = Ft5;
            f->phasing2 = &LALInspiralPhasing2_5PN;
            f->timing2 = &LALInspiralTiming2_5PN;
            f->phasing3 = &LALInspiralPhasing3_5PN;
            f->frequency3 = &LALInspiralFrequency3_5PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case IMRPhenomA:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp5;
            break;
         case PadeF1:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            /* FIXME: TODO: DO SOMETHING HERE!!!! */
            break;
      }
      break;
      case threePN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
      case SpinTaylor:
/*
   vlsoT6 is as yet undetermined and vlsoT4 is too large in
   certain cases (TaylorT2 crashes for (1.4,10)); using vlsoT2;
*/
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt6;
            f->flux = Ft6;
            f->phasing2 = &LALInspiralPhasing2_6PN;
            f->timing2 = &LALInspiralTiming2_6PN;
            f->phasing3 = &LALInspiralPhasing3_6PN;
            f->frequency3 = &LALInspiralFrequency3_6PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case IMRPhenomA:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp6;
            break;
         case PadeF1:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            /* FIXME: TODO: DO SOMETHING HERE!!!! */
            break;
      }
      break;
      case threePointFivePN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt6;
            f->flux = Ft7;
            f->phasing2 = &LALInspiralPhasing2_7PN;
            f->timing2 = &LALInspiralTiming2_7PN;
            f->phasing3 = &LALInspiralPhasing3_7PN;
            f->frequency3 = &LALInspiralFrequency3_7PN;
            break;
         case PadeT1:
         case EOB:
         case EOBNR:
         case IMRPhenomA:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp7;
            break;
         case PadeF1:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case pseudoFourPN:
      switch (params->approximant)
      {
         case Eccentricity:
            ABORT(status, LALINSPIRALH_EORDERMISSING, LALINSPIRALH_MSGEORDERMISSING);
            break;
         case EOB:
         case EOBNR:
         case IMRPhenomA:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp7;
            break;
         case AmpCorPPN:
         case TaylorT1:
         case TaylorT2:
         case TaylorT3:
         case TaylorF1:
         case TaylorF2:
         case SpinTaylorT3:
         case SpinTaylor:
         case PadeT1:
         case PadeF1:
         case TaylorEt:
         case TaylorT4:
         case TaylorN:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      default:
         break;
   }


   switch (params->approximant){
   case AmpCorPPN:
   case TaylorT1:
   case TaylorT2:
   case TaylorT3:
   case TaylorF1:
   case EOB:
   case EOBNR:
   case PadeT1:
   case PadeF1:
   case TaylorF2:
   case SpinTaylorT3:
   case SpinTaylor:
   case TaylorEt:
   case TaylorT4:
   case TaylorN:
     ak->flso = pow(ak->vlso,3.)/(LAL_PI * ak->totalmass);

     if (ak->fn) {
       vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
       ak->vn = (vn < vlso) ? vn :  vlso;
     }

     in1.t=0.0;
     in1.v0=ak->v0;
     in1.t0=ak->t0;
     in1.vlso=ak->vlso;
     in1.totalmass = ak->totalmass;
     in1.dEnergy = f->dEnergy;
     in1.flux = f->flux;
     in1.coeffs = ak;

     in2 = (void *) &in1;

     LALInspiralTofV(status->statusPtr, &tofv, ak->vn, in2);
     CHECKSTATUSPTR(status);

     ak->tn = -tofv - ak->samplinginterval;
     params->fCutoff = ak->fn = pow(ak->vn, 3.)/(LAL_PI * ak->totalmass);
     /*
       for (v=0; v<ak->vn; v+=0.001)
       {
       FtN = Ft0(v,ak);
       printf("%e %e %e %e %e %e %e\n", v,
       Ft2(v,ak)/FtN, Ft3(v,ak)/FtN, Ft4(v,ak)/FtN, Ft5(v,ak)/FtN,
       Ft6(v,ak)/FtN, Ft7(v,ak)/FtN);
       }
       exit(0);
     */
     break;
 case BCV:
 case BCVSpin:
 case IMRPhenomA:
   ak->tn = 100.;
   break;
 case Eccentricity:
   /* The eccentric waveforms contain harmonic, so similarly to amplitude corrected waveforms
    * the duration are longer than non eccentric waveform and starts at 2fl/3*/
   ak->tn = 5.*ak->totalmass/256./ak->eta/pow(LAL_PI*ak->totalmass*params->fLower/3.*2.,8./3.);
   ak->flso = pow(ak->vlso,3.)/(LAL_PI * ak->totalmass);
   break;
 default:
   ABORT( status, 9999, "Unknown case in switch." );
}

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

