/*  <lalVerbatim file="LALInspiralChooseModelCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralChooseModel.c}}

Module to set the pointers to the required energy and flux functions.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralChooseModelCP}
\index{\verb&LALInspiralChooseModel()&}

\subsubsection*{Description}

Module to set the pointers to the required energy and flux functions $E^{\prime}(v)$ and $\mathcal{F}(v)$.
These are used to solve the gravitational wave phasing formula,

\begin{eqnarray}
t(v) & = & t_{0} - m \int_{v_{0}}^{v} \,
\frac{E'(v)}{{\cal F}(v)} \, dv, \nonumber \\
\phi (v) & = & \phi_{0} - 2 \int_{v_{0}}^{v}  v^3 \,
\frac{E'(v)}{{\cal F}(v)} \, dv,
\end{eqnarray}

where $v=(\pi m F)^{1/3}$ is an invariantly defined velocity, $F$ is the instantaneous GW frequency, and
$m$ is the total mass of the binary.

The expressions for $E^{\prime}(v)$ and $\mathcal{F}(v)$ may be written either in the form of a Taylor
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

The Taylor series expansion of the flux function in the $\eta \neq 0$ case is known up to order $v^{5}$,
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
where the coefficients $f_{k}$ are defined in terms if the original coefficients $A_{k}$ of the Taylor 
series, $f_{k}=A_{k} - F_{k-1}/v_{pole}$.
Then we re--arrange to define a new $\mathcal{F}(v)$ given by
\begin{equation}
\mathcal{F}_{P_{n}}(v) = \frac{1}{\left( 1 - \frac{v}{v^{P_{n}}_{pole}(\eta)} \right)} \,\, f_{P_{n}}(v)
\end{equation}
where $v^{P_{n}}_{pole}$ denotes the pole velocity defined by the $v^{n}$ approximant of $e(x)$.

This module defines all the energy and flux functions in terms of the coefficients which need to have been
previously defined by the function \texttt{InspiralSetup}.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
None.

\subsubsection*{Notes}
See Damour, Iyer and Sathyaprakash, PRD 57, 885, 1998 for further details.
The Pade Approximant for the 1PN expansion is undefined, and so cannot be used.


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
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(4.*v))*v6);
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
        + ak->FTa5*v4*v + (ak->FTa6 + ak->FTl6*log(4.*v))*v6 + ak->FTa7*v6*v);
   return (flux);
}


/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 ep0(REAL8 v, expnCoeffs *ak) 
{ /* </lalVerbatim>  */
   REAL8 x, energy;
   ak = NULL;
   x = v*v;
   energy = -x;
   return (energy);
}

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

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp3(REAL8 v, expnCoeffs *ak) 
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v))) * (1.-v/ak->vpoleP4));
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
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v / (1.+ak->fPa4*v)))) * (1.-v/ak->vpoleP4));
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
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v / (1.+ak->fPa4*v / (1.+ak->fPa5*v))))) * (1.-v/ak->vpoleP4));
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
static REAL8 Fp6(REAL8 v, expnCoeffs *ak) 
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v 
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v)))))) 
        * (1.-v/ak->vpoleP6));
   return (flux);
}

static REAL8 Fp7(REAL8 v, expnCoeffs *ak) 
{ /* </lalVerbatim>  */
   REAL8 flux,v2,v4,v8,v10;
   v2 = v*v;
   v4 = v2*v2;
   v8 = v4*v4;
   v10 = v8*v2;
   flux = ak->fPaN * v10/ ((1.+ak->fPa1*v/(1.+ak->fPa2*v/ (1.+ak->fPa3*v 
        / (1.+ak->fPa4*v / (1.+ak->fPa5*v / (1.+ak->fPa6*v / (1.+ak->fPa7*v))))))) 
        * (1.-v/ak->vpoleP6));
   return (flux);
}

/*  <lalVerbatim file="LALInspiralChooseModelCP"> */
void 
LALInspiralChooseModel(
   LALStatus *status,
   expnFunc *f,
   expnCoeffs *ak,
   InspiralTemplate *params)
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
   ASSERT((INT4)params->order <= 7, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   vlso = 0;

   switch (params->order) 
   {
      case newtonian:
      case oneHalfPN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
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
            ak->vn = ak->vlso = vlso = ak->vlsoP0;
            f->dEnergy = dEp0;
            f->flux = Fp0;
            break;
         case EOB:
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case onePN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
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
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case onePointFivePN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorT3:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft3;
            f->phasing3 = &LALInspiralPhasing3_3PN;
            f->frequency3 = &LALInspiralFrequency3_3PN;
            break;
         case TaylorF2:
         case TaylorT2:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt2;
            f->flux = Ft3;
            f->phasing2 = &LALInspiralPhasing2_3PN;
            f->timing2 = &LALInspiralTiming2_3PN;
            break;
         case PadeT1:
         case PadeF1:
            ak->vn = ak->vlso = vlso = ak->vlsoP0;
            f->dEnergy = dEp2;
            f->flux = Fp3;
            break;
         case EOB:
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case twoPN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
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
         case PadeF1:
         case EOB:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp4;
            break;
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      break;
      case twoPointFivePN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
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
         case PadeF1:
         case EOB:
            ak->vn = ak->vlso = vlso = ak->vlsoP4;
            f->dEnergy = dEp4;
            f->flux = Fp5;
            break;
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
      }
      break;
      case threePN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
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
         case PadeF1:
         case EOB:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp6;
            break;
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
      }
      break;
      case threePointFivePN:
      switch (params->approximant) 
      {
         case TaylorT1:
         case TaylorF1:
         case TaylorF2:
         case TaylorT2:
         case TaylorT3:
            ak->vn = ak->vlso = vlso = ak->vlsoT2;
            f->dEnergy = dEt6;
            f->flux = Ft7;
            f->phasing2 = &LALInspiralPhasing2_7PN;
            f->timing2 = &LALInspiralTiming2_7PN;
            f->phasing3 = &LALInspiralPhasing3_7PN;
            f->frequency3 = &LALInspiralFrequency3_7PN;
            break;
         case PadeT1:
         case PadeF1:
         case EOB:
            ak->vn = ak->vlso = vlso = ak->vlsoP6;
            f->dEnergy = dEp6;
            f->flux = Fp7;
            break;
         case DJS:
         case INSPA:
         case IRSPA:
            ABORT(status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
            break;
         default:
            break;
      }
      default:
         break;
   }

   ak->flso = pow(ak->vlso,3.)/(LAL_PI * ak->totalmass);

   if (ak->fn) {
      vn = pow(LAL_PI * ak->totalmass * ak->fn, oneby3);
      ak->vn = (vn < vlso) ? vn :  vlso;
   } 

   in1.t=0.0;
   in1.v0=ak->v0;
   in1.t0=ak->t0;
   in1.totalmass = ak->totalmass;
   in1.dEnergy = f->dEnergy;
   in1.flux = f->flux;
   in1.coeffs = ak;
   
   in2 = (void *) &in1;      
   
   LALInspiralTofV(status->statusPtr, &tofv, ak->vn, in2);
   CHECKSTATUSPTR(status);
   
   ak->tn = -tofv - ak->samplinginterval;
   ak->fn = pow(ak->vn, 3.)/(LAL_PI * ak->totalmass);
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
   DETATCHSTATUSPTR(status);
   RETURN (status); 
}

