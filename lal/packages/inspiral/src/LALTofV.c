/*  <lalVerbatim file="LALTofVCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTofV.c}}

Module which calculates the quantity which we denote tofv, which is given by the following
equation
\begin{equation}
\mathrm{tofv} = t(v) - t_{0} + m \int_{v_{0}}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv \,\,.
\end{equation}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTofVCP}
\index{\verb&LALTofV()&}

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALDRombergIntegrate}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALTofVCV}}

</lalLaTeX>  */



#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

NRCSID (LALTOFVC, "$Id$");

/*  <lalVerbatim file="LALTofVCP"> */
void LALTofV (LALStatus *status,
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


   INITSTATUS (status, "LALTofV", LALTOFVC);
   ATTATCHSTATUSPTR(status);

   ASSERT (tofv, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   sign = 1.0;


   in1 = (TofVIn *) params;

   intinp.function = LALTofVIntegrand;
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

