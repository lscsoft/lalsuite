/*  <lalVerbatim file="LALInspiralMomentsCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralMoments.c}}

Module to calculate the moment of the noise power spectral density. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralMomentsCP}
\index{\verb&LALInspiralMoments()&}

\subsubsection*{Description}

Purpose: Given the exponent {\sf ndx} and limits of integration 
{\sf xs} and {\sf xc} this function returns the moment of the 
power spectral density specified by {\sf choice}.

\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
LALGEOPsd
LALLIGOIPsd
TAMAPsd
VIRGOPsd
LALInspiralMomentsIntegrand
LALDRombergIntegrate
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralMomentsCV}}

</lalLaTeX>  */




/* ---------------------------------------------------------------------
  Author: B.S.Sathyaprakash, U.W.C., Cardiff.
  Revision History: Last Update 15.7.97.
  Purpose: Given the exponent {\sf ndx} and limits of integration 
     {\sf xs} and {\sf xc} this function returns the moment of the 
     power spectral density specified by {\sf choice}.
  Inputs: 
  Outputs:
     ${\sf moments} = \int_{xs}^{xc} x^{\sf -ndx}~dx/S_h(x).$
*----------------------------------------------------------------------- */

#include <lal/LALInspiralBank.h>
#include <lal/Integrate.h>
void qromb (LALStatus *status, REAL8 *moment, DIntegrateIn *In, void *funcparams);

NRCSID(LALINSPIRALMOMENTSC, "$Id$");

/*  <lalVerbatim file="LALInspiralMomentsCP"> */
void LALInspiralMoments(LALStatus *status,
                        REAL8 *moment,
                        InspiralMomentsIn pars)
{ /* </lalVerbatim> */

   DIntegrateIn In;
   void *funcparams;
   InspiralMomentsIntegrandIn intIn;

   INITSTATUS (status, "LALInspiralMoments", LALINSPIRALMOMENTSC);
   ATTATCHSTATUSPTR(status);
   intIn.ndx = pars.ndx;
   switch (pars.detector) {
      case geo:
         intIn.NoisePsd = &LALGEOPsd;
         break;
      case ligo:
         intIn.NoisePsd = &LALLIGOIPsd;
         break;
      case tama:
         intIn.NoisePsd = &LALTAMAPsd;
         break;
      case virgo:
         intIn.NoisePsd = &LALVIRGOPsd;
         break;
   }
   funcparams = (void *) &(intIn);
   In.function = LALInspiralMomentsIntegrand;
   In.xmin = pars.xmin;
   In.xmax = pars.xmax;
   In.type = ClosedInterval;

/*   qromb (status->statusPtr, moment, &In, funcparams);
*/
   LALDRombergIntegrate (status->statusPtr, moment, &In, funcparams);

   *moment /= pars.norm;
/*
   DRombergIntegrate (status->statusPtr, moment, &In, funcparams);
   fprintf(stderr, "%e %e %e %e\n", pars.xmin, pars.xmax, 3*pars.q, *moment);
*/
   DETATCHSTATUSPTR(status);
   RETURN (status);
}
