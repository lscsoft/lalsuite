/*  <lalVerbatim file="LALInspiralWave1CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1.c} and \texttt{LALInspiralWave1Templates.c}}

The code \texttt{LALInspiralWave1} generates an time-domain inspiral waveform corresponding to the 
\texttt{approximant} \texttt{TaylorT1} and \texttt{PadeT1} as outlined in the
documentation for the function \texttt{LALInspiralWave}. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave1CP}
\index{\verb&LALInspiralWave1()&}
\begin{itemize}
\item {\tt signal:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALInspiralWave1TemplatesCP}
\index{\verb&LALInspiralWave1Templates()&}
\begin{itemize}
\item {\tt signal1:} Output containing the 0-phase inspiral waveform.
\item {\tt signal2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralWave1} is called if the user has specified the 
\texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.
{\tt LALInspiralWave1Templates} is exactly the same as \texttt{LALInspiralWave1,} except that
it generates two templates one for which the starting phase is 
\texttt{params.startPhase} and the other for which the phase is
\texttt{params.startPhase + $\pi/2$}.


\subsubsection*{Algorithm}
This code uses a fourth-order Runge-Kutta algorithm to solve the ODEs 
in Equation (\ref{eq:ode2}).

\subsubsection*{Uses}

\texttt{LALInspiralSetup}\\
\texttt{LALInspiralChooseModel}\\
\texttt{LALInspiralVelocity}\\
\texttt{LALInspiralPhasing1}\\
\texttt{LALInspiralDerivatives}\\
\texttt{LALRungeKutta4}.
 

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave1CV}}

</lalLaTeX>  */

/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVE1C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1CP"> */
void 
LALInspiralWave1(
   LALStatus        *status,
   REAL4Vector      *signal,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 m, dt, t, v, p, h, f, fu, fHigh, piM;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS(status, "LALInspiralWave1", LALINSPIRALWAVE1C);
   ATTATCHSTATUSPTR(status);

   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
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
   fu = params->fCutoff;
   if (fu) 
      fHigh = (fu < ak.flso) ? fu : ak.flso; 
   else 
      fHigh = ak.flso;
/* 
    Check that the highest frequency is less than half 
    the sampling frequency - the Nyquist theorem 
*/

   ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


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
   {
       *(signal->data + count) = 0.;
       count++;
   }

   t = 0.0;
   piM = LAL_PI * m;
   do {
      h = params->signalAmplitude * v*v * cos(p);
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      *(signal->data+count) = (REAL4) h;
      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;
   } while (t < ak.tn && f<fHigh);
   params->fFinal = f;

   while (count < (int)signal->length) 
   {
       *(signal->data + count) = 0.;
       count++;
   }

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}



/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

NRCSID (LALINSPIRALWAVE1TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1TemplatesCP"> */
void 
LALInspiralWave1Templates(
   LALStatus        *status,
   REAL4Vector      *signal1,
   REAL4Vector      *signal2,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 amp, m, dt, t, v, p, h1, h2, f, fu, fHigh, piM;
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

   piM = LAL_PI * m;
   f = (v*v*v)/piM;

   fu = params->fCutoff;
   if (fu) 
      fHigh = (fu < ak.flso) ? fu : ak.flso; 
   else 
      fHigh = ak.flso;
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
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

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
      f = v*v*v/piM;
   } while (t < ak.tn && f<fHigh);
   params->fFinal = f;

   while (count < (int)signal1->length) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}

/* 
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

NRCSID (LALINSPIRALWAVE1FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1ForInjectionCP"> */
void 
LALInspiralWave1ForInjection(
   LALStatus        *status,
   REAL4Vector      *inject_hc,
   REAL4Vector      *inject_hp,
   REAL4Vector      *inject_phase,
   REAL4Vector      *inject_freq,			   
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 n=2, count;
   REAL8 amp, m, dt, t, v, p, h1, h2, f, fu, fHigh, piM;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

   REAL8 unitHz; 
   unitHz = (params->mass1 +params->mass2) *LAL_MTSUN_SI*(REAL8)LAL_PI;


   INITSTATUS(status, "LALInspiralWave1ForInjection", LALINSPIRALWAVE1TEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(inject_hc,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_hp,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_hc->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_hp->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   
   ASSERT(inject_phase,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_freq,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_freq->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(inject_phase->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  
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

   piM = LAL_PI * m;
   f = (v*v*v)/piM;

   fu = params->fCutoff;
   if (fu) 
      fHigh = (fu < ak.flso) ? fu : ak.flso; 
   else 
      fHigh = ak.flso;
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
   {

     *(inject_hc->data + count) = 
       *(inject_hp->data + count) = 
       *(inject_phase ->data +count) =
       *(inject_freq->data +count) = 0.;
      count++;
   }

   t = 0.0;
   do {
      amp = params->signalAmplitude * v*v;
      h1 = amp * cos(p);
      h2 = amp * cos(p+LAL_PI_2);


      *(inject_hc->data + count) = (REAL4) amp;
      *(inject_hp->data + count) = (REAL4) amp;
      *(inject_phase->data + count) = (REAL8) (p);
      *(inject_freq->data + count) = (REAL4) v*v*v / unitHz;
         
      LALInspiralDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      in4.dydx = &dvalues;
      in4.x=t;
      LALRungeKutta4(status->statusPtr, &valuesNew, &in4, funcParams);
      CHECKSTATUSPTR(status);
      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);
      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;
   } while (t < ak.tn && f<fHigh);
   params->fFinal = f;

   while (count < (INT4)inject_hc->length) 
   {
     *(inject_hc->data + count) = 
       *(inject_hp->data + count) = 
       *(inject_phase ->data +count) =
       *(inject_freq->data +count) =0.;
     count++;
   }
   
   LALFree(dummy.data);
   
   DETATCHSTATUSPTR(status);
   RETURN (status);

}
