/*  <lalVerbatim file="LALInspiralWave1CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1.c}}

The code \texttt{LALInspiralWave1.c} generates an inspiral waveform using method \texttt{one} as outlined in the
documentation for the function \texttt{InspiralWave}. This means that the integrals involving the energy
and flux function are solved numerically.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave1CP}
\index{\verb&LALInspiralWave1()&}

\subsubsection*{Description}

This function is called if the user has specified method \texttt{one} in the \texttt{enum} which is of type
\texttt{enum Method}. This will be set in the input structure, which is of type \texttt{InspiralTemplate}.

This means that the the gravitational wave phasing
formula is solved as follows. The formula is
 
\begin{eqnarray}
t(v) & = & t_{0} - m \int_{v_{0}}^{v} \,
\frac{E'(v)}{{\cal F}(v)} \, dv, \nonumber \\
\phi (v) & = & \phi_{0} - 2 \int_{v_{0}}^{v}  v^3 \,
\frac{E'(v)}{{\cal F}(v)} \, dv,
\label{TimeDomain2PhasingFormula}
\end{eqnarray}
 
where $v=(\pi m F)^{1/3}$ is an invariantly defined velocity, $F$ is the instantaneous GW frequency, and
$m$ is the total mass of the binary.
 
The variable $t$ is in our hands, and the phasing formula is solved by finding the value of $v$ (for a
given $t$) which satisfies the first of the equations in Eq. (\ref{TimeDomain2PhasingFormula}). This value of $v$ is
then substituted into the second equation to yield $\phi(t)$.
 
Method \texttt{one} leaves  $E^{\prime}(v)/\mathcal{F}(v)$ as a series in terms of ascending powers of $x$
and integrate the equations numerically. This numerical integration is done by replacing the integrals by a
pair of first order coupled differential equations, as follows.
 
Instead of writing the first of the equations in Eq. (\ref{TimeDomain2PhasingFormula}), we may write
 
\begin{equation}
\frac{dt(v)}{dv} = -m \frac{E^{\prime}(v)}{\mathcal{F}(v)} \,\,,
\end{equation}
 
therefore
 
\begin{equation}
\frac{dv}{dt} = - \frac{\mathcal{F}(v)}{m E^{\prime}(v)} \,\,.
\label{TimeDomain2Ode1}
\end{equation}           
 
Instead of the second equation in Eq. (\ref{TimeDomain2PhasingFormula}) we may write
 
\begin{equation}
\frac{d \phi (t)}{dt} = 2 \pi F(t)
\end{equation}
 
where
 
\begin{equation}
F = \frac{v^{3}}{\pi m} \,\,,
\end{equation}
 
therefore
 
\begin{equation}
\frac{d \phi(t)}{dt} = \frac{2v^{3}}{m}
\label{TimeDomain2Ode2}
\end{equation}
 
So Eq. (\ref{TimeDomain2Ode1}) and Eq. (\ref{TimeDomain2Ode2}) are our coupled first order differential equations to solve. This
we do using a fourth--order Runge--Kutta algorithm.
 
 
 
Inside the input structure of type \texttt{InspiralTemplate}, the parameter \texttt{order}, which is of type \texttt{enum Order}, lets the user choose to which order of
post--Newtonian expansion they would like to go. The parameter \texttt{domain}, which is of type
\texttt{enum Domain}, lets the user to choose to generate the waveform in the time or frequency domain. The
user has have the additional choice of specifying the form of the expansion which defines
$E^{\prime}(v)/\mathcal{F}(v)$. This ratio can be expressed in the form of a \emph{Taylor series}
(T--Approximants) or by using \emph{P--Approximants}. This choice is handled by the parameter
\texttt{approximant}, which is of type \texttt{enum Approximant}. 

\subsubsection*{Algorithm}
This code uses a fourth--order Runge--Kutta algorithm to solve the integrals as a pair of coupled
first--order differential equations.

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
void LALInspiralWave1(LALStatus *status,
		    REAL4Vector *signal,
		    InspiralTemplate *params)
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

   while (count < (int)signal->length) 
   {
       *(signal->data + count) = 0.;
       count++;
   }

   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}
