/*  <lalVerbatim file="LALTimeDomain2TemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTimeDomain2Templates.c}}

Module to generate two inspiral waveforms simultaneously, which differ in phase by $\pi/2$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALTimeDomain2TemplatesCP}
\index{\texttt{LALTimeDomain2Templates()}}

\subsubsection*{Description}

Module to generate two inspiral waveforms simultaneously, which differ in phase by $\pi/2$.
This is so that when we need to maximise over the phase at time of arrival by seperately correlating a
signal with a zero and $\pi/2$ phase waveform, we do not need to call the waveform generation function
twice.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{InspiralSetup}\\
\texttt{ChooseModel}\\
\texttt{TappRpnTdomTime}\\
\texttt{InspiralVelocity}\\
\texttt{InspiralPhase}\\
\texttt{InspiralDerivatives}\\
\texttt{RungeKutta4}.

\subsubsection*{Notes}

See the documentation for the module \texttt{LALTimeDomain2} for further details.

\vfill{\footnotesize\input{LALTimeDomain2TemplatesCV}}

</lalLaTeX>  */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALTIMEDOMAIN2TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALTimeDomain2TemplatesCP"> */
void LALTimeDomain2Templates(LALStatus *status,
		             REAL4Vector *signal1,
		             REAL4Vector *signal2,
		             InspiralTemplate *params)
 { /* </lalVerbatim> */

   INT4 n=2, count;
   double m, dt, t, v, p, h1, h2, f;
   REAL8Vector values;
   REAL8Vector dvalues;
   REAL8Vector valuesNew;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   REAL8Vector yt;
   REAL8Vector dym;
   REAL8Vector dyt;
   expnCoeffs ak;
   expnFunc func;


   INITSTATUS(status, "LALTimeDomain2Templates", LALTIMEDOMAIN2TEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   values.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
   dvalues.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
   valuesNew.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);

   m = ak.totalmass;
   dt = 1./params->tSampling;

   t = 0.0;
   in1.t = t;
   in1.t0=ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;

   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;

   funcParams = (void *) &in3;
   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   f = (v*v*v)/(LAL_PI*m);

   LALInspiralPhase(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);


   *(values.data) = v; 
   *(values.data+1) = p;

   dym.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
   dyt.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
   yt.data = (REAL8 *)LALMalloc(sizeof(REAL8)*n);
   
   in4.function = LALInspiralDerivatives;
   in4.x = t;
   in4.y = &values;
   in4.h = dt;
   in4.n = n;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   count = 0;
   while (count < params->nStartPad) 
	 *(signal1->data + count) = *(signal2->data + count++) = 0.;

   t = 0.0;
   do {
      h1 = params->signalAmplitude * v*v * cos(p);
      h2 = params->signalAmplitude * v*v * cos(p+LAL_PI_2);
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      *(signal1->data+count) = (REAL4) h1;
      *(signal2->data+count) = (REAL4) h2;
      t = (++count-params->nStartPad) * dt;
   } while (t < ak.tn);

   while (count < (int)signal1->length) 
	 *(signal1->data + count) = *(signal2->data + count++) = 0.;


   LALFree(yt.data);
   LALFree(dyt.data);
   LALFree(dym.data);
   LALFree(values.data);
   LALFree(dvalues.data);
   LALFree(valuesNew.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}
