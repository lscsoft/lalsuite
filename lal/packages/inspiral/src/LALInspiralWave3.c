/*  <lalVerbatim file="LALInspiralWave3CV">

Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave3.c}}

This module generates a time-domain chirp waveform of type {\tt TaylorT3}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave3CP}
\index{\verb&LALInspiralWave3()&}
\begin{itemize}
\item {\tt output:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}
The {\tt TaylorT3} approximant corresponds to the case wherein 
the phase of the waveform is given as an explicit function of time
as in Equation (\ref{eq:InspiralWavePhase3}). 

 
\subsubsection*{Algorithm}

\subsubsection*{Uses}

\texttt{LALInspiralParameterCalc} \\
\texttt{LALInspiralChooseModel} \\
\texttt{LALInspiralSetup} \\
\texttt{LALInspiralPhasing3} (via expnFunc)\\
\texttt{LALInspiralFrequency3}. (via expnFunc)\\


\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave3CV}}

</lalLaTeX>  */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

typedef struct
{
	void (*func)(LALStatus *status, REAL8 *f, REAL8 tC, expnCoeffs *ak);
	expnCoeffs ak;
} 
ChirptimeFromFreqIn;

static void LALInspiralFrequency3Wrapper(LALStatus *status, REAL8 *f, REAL8 tC, void *pars);

NRCSID (LALINSPIRALWAVE3C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave3CP"> */

void 
LALInspiralWave3 (
   LALStatus        *status,
   REAL4Vector      *output, 
   InspiralTemplate *params
   )

{ /* </lalVerbatim>  */

  INT4 i, startShift, count;
  REAL8 dt, fu, ieta, eta, tc, totalMass, t, td, c1, phi0, phi;
  REAL8 v, f, fHigh, amp, tmax, fOld, phase;
  DFindRootIn rootIn;
  expnFunc func;
  expnCoeffs ak;
  ChirptimeFromFreqIn timeIn;
  void *pars;


  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  ASSERT(output, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(output->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
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


  eta = params->eta;   /* Symmetric mass ratio  */
  ieta = params->ieta;
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */
  /* constant that appears in the definition of the time parameter Theta */
  c1 = eta/(5.*totalMass);
 
  /*
   * In Jan 2003 we realized that the tC determined as a sum of chirp times is
   * not quite the tC that should enter the definition of Theta in the expression
   * for the frequency as a function of time (see DIS3, 2000). This is because
   * chirp times are obtained by inverting t(f). Rather tC should be obtained by
   * solving the equation f0 - f(tC) = 0. This is what is implemented below.
   */

  timeIn.func = func.frequency3;
  timeIn.ak = ak;
  rootIn.function = &LALInspiralFrequency3Wrapper;
  rootIn.xmin = c1*params->tC/2.;
  rootIn.xmax = c1*params->tC;
  rootIn.xacc = 1.e-6;
  pars = (void*) &timeIn;
  /* tc is the instant of coalescence */
  LALDBisectionFindRoot (status->statusPtr, &tc, &rootIn, pars);

  tc /= c1;

  tc += params->startTime;       /* Add user given startTime to instant of 
				     coalescence of the compact objects */

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

  i=0; while (i<startShift) 
  {
      output->data[i] = 0.0;
      i++;
  }

  t=0.0;
  td = c1*(tc-t);
  func.phasing3(status->statusPtr, &phase, td, &ak);
  CHECKSTATUSPTR(status);
  func.frequency3(status->statusPtr, &f, td, &ak);
  CHECKSTATUSPTR(status);
  phi0=-phase+phi;

  /*
  fprintf(stderr, "Starting frequency=%e\n", f);
   */

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f<fHigh && t<tmax && f>fOld) 
  {
    fOld = f; 
    v = pow(f*LAL_PI*totalMass, oneby3);
    amp = v*v; 
    output->data[i++] = (REAL4) (params->signalAmplitude * amp * cos(phase+phi0));
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
  while (i < (int)output->length) 
  {
      output->data[i]=0.0;
      i++;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

static void LALInspiralFrequency3Wrapper(LALStatus *status, REAL8 *f, REAL8 tC, void *pars)
{

  ChirptimeFromFreqIn *in;
  REAL8 freq;
  INITSTATUS (status, "LALInspiralWave3", LALINSPIRALWAVE3C);
  ATTATCHSTATUSPTR(status);

  in = (ChirptimeFromFreqIn *) pars;
  in->func(status->statusPtr, &freq, tC, &(in->ak));
  *f = freq - in->ak.f0;

  /*
  fprintf(stderr, "Here freq=%e f=%e tc=%e f0=%e\n", freq, *f, tC, in->ak.f0);
   */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
