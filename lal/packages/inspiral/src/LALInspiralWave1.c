/*  <lalVerbatim file="LALInspiralWave1CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1.c}}

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

\subsubsection*{Description}

This function is called if the user has specified the \texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.

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
