/*  <lalVerbatim file="LALInspiralWave3TemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave3Templates.c}}

Module to generate two inspiral waveforms simultaneously, 
which differ in phase by $\pi/2$.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave3TemplatesCP}
\index{\verb&LALInspiralWave3Templates()&}
\begin{itemize}
\item {\tt output1:} Output containing the 0-phase inspiral waveform.
\item {\tt output2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

Module to generate two inspiral waveforms simultaneously, which differ in 
phase by $\pi/2$.  This is so that when we need to maximise over the phase 
at time of arrival by seperately correlating a signal with a zero and $\pi/2$ 
phase waveform, we do not need to call the waveform generation function
twice.

\subsubsection*{Algorithm}


\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc} \\
\texttt{LALInspiralPhasing3} \\
\texttt{LALInspiralFrequency3}. \\

\subsubsection*{Notes}

See the documentation for the module \texttt{LALInspiralWave3} 
for further details.

\vfill{\footnotesize\input{LALInspiralWave3TemplatesCV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALWAVE3TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3TemplatesCP"> */

void 
LALInspiralWave3Templates (
   LALStatus        *status,
   REAL4Vector      *output1, 
   REAL4Vector      *output2, 
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi1, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  expnFunc func;
  expnCoeffs ak;


  INITSTATUS (status, "LALInspiralWave3Templates", LALINSPIRALWAVE3TEMPLATESC);
  ATTATCHSTATUSPTR(status);

  ASSERT(output1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL); 
  ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  LALInspiralSetup (status->statusPtr, &ak, params);
  CHECKSTATUSPTR(status);
  LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
  CHECKSTATUSPTR(status);

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */

/* Calculate the three unknown paramaters in (m1,m2,M,eta,mu) from the two
   which are given.  */

  LALInspiralParameterCalc (status->statusPtr, params);
  CHECKSTATUSPTR(status);

  ASSERT(params->totalMass >0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(params->eta >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  tc=params->tC;       /* Instant of coalescence of the compact objects */
  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */


/* 
   If flso is less than the user inputted upper frequency cutoff fu, 
   then the waveforn is truncated at f=flso.  
*/

  if (fu) 
     fHigh = (fu < ak.flso) ? fu : ak.flso; 
  else 
     fHigh = ak.flso;

/* 
   Check that the highest frequency is less than half the sampling frequency - 
   the Nyquist theorum 
*/

  ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);
  i=0; while (i<startShift) 
  {
      output1->data[i] = output2->data[i] = 0.0;
      i++;
  }

  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;
  phi1=phi0+LAL_PI_2;

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = params->signalAmplitude * v*v; 
    output1->data[i] = (REAL4) amp * cos(phase+phi0);
    output2->data[i] = (REAL4) amp * cos(phase+phi1);
    ++i;
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    func.phasing3(status->statusPtr, &phase, td, &ak);
    CHECKSTATUSPTR(status);
    func.frequency3(status->statusPtr, &f, td, &ak);
    CHECKSTATUSPTR(status); 
  }
  params->fFinal = fOld;
  
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/
  while (i < (int)output1->length) 
  {
      output1->data[i]=output2->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
