/*  <lalVerbatim file="LALInspiralWave1TemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1Templates.c}}

This is exactly the same as \texttt{LALInspiralWave1,} except that
it generates two templates one for which the starting phase is 
\texttt{params.startPhase} and the other for which the phase is
\texttt{params.startPhase+$\pi/2$}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave1TemplatesCP}
\index{\verb&LALInspiralWave1Templates()&}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\texttt{LALInspiralSetup}\\
\texttt{LALInspiralChooseModel}\\
\texttt{LALInspiralVelocity}\\
\texttt{LALInspiralPhasing1}\\
\texttt{LALInspiralDerivatives}\\
\texttt{LALRungeKutta4}.

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave1TemplatesCV}}

</lalLaTeX>  */

/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVE1TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1TemplatesCP"> */
void LALInspiralWave1Templates(LALStatus *status,
		    REAL4Vector *signal1,
		    REAL4Vector *signal2,
		    InspiralTemplate *params)
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 amp, m, dt, t, v, p, h1, h2, f;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS(status, "LALInspiralWave1Templates", LALINSPIRALWAVE1TEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = n;
   dummy.length = n * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[n];
   valuesNew.data = &dummy.data[2*n];
   yt.data = &dummy.data[3*n];
   dym.data = &dummy.data[4*n];
   dyt.data = &dummy.data[5*n];

   m = ak.totalmass;
   dt = 1./params->tSampling;

   ASSERT(ak.totalmass > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

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

   LALInspiralPhasing1(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);


   *(values.data) = v; 
   *(values.data+1) = p;

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
      amp = params->signalAmplitude * v*v;
      h1 = amp * cos(p);
      h2 = amp * cos(p+LAL_PI_2);
      *(signal1->data + count) = (REAL4) h1;
      *(signal2->data + count) = (REAL4) h2;
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      t = (++count-params->nStartPad) * dt;
   } while (t < ak.tn);

   while (count < (int)signal1->length) 
      *(signal1->data + count) = *(signal2->data + count++) = 0.;

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}
