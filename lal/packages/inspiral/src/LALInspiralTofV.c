/*  <lalVerbatim file="LALInspiralTofVCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTofV.c}}

Module which calculates the quantity which we denote tofv, which is given by the following
equation
\begin{equation}
\mathrm{tofv} = t(v) - t_{0} + m \int_{v_{0}}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv \,\,.
\end{equation}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralTofVCP}
\index{\verb&LALInspiralTofV()&}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALDRombergIntegrate}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralTofVCV}}

</lalLaTeX>  */



#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

NRCSID (LALINSPIRALTOFVC, "$Id$");

/*  <lalVerbatim file="LALInspiralTofVCP"> */
void LALInspiralTofV (LALStatus *status,
	      REAL8 *tofv,
	      REAL8 v,
	      void *params)
{ /* </lalVerbatim>  */

   void *funcParams;
   DIntegrateIn intinp;
   TofVIntegrandIn in2;
   TofVIn *in1;
   REAL8 answer;
   REAL8 sign;


   INITSTATUS (status, "LALInspiralTofV", LALINSPIRALTOFVC);
   ATTATCHSTATUSPTR(status);

   ASSERT (tofv, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(v > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(v < 1., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   sign = 1.0;


   in1 = (TofVIn *) params;

   intinp.function = LALInspiralTofVIntegrand;
   intinp.xmin = in1->v0;
   intinp.xmax = v;
   intinp.type = ClosedInterval;


   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;

   if (v==in1->v0)
   {
     *tofv = in1->t - in1->t0;
     DETATCHSTATUSPTR(status);
     RETURN (status);
   }

   if(in1->v0 > v)
   {
      intinp.xmin = v;
      intinp.xmax = in1->v0;
      sign = -1.0;
   }
	
   LALDRombergIntegrate (status->statusPtr, &answer, &intinp, funcParams);
   CHECKSTATUSPTR(status);

   *tofv = in1->t - in1->t0 + in1->totalmass*answer*sign;

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

