/******************************** <lalVerbatim file="LALRandomSignalCV">
Author: Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{RandomSignal.c}}
\label{ss:RandomSignal.c}

Test code for the inspiral wave generation and noisemodels modules.

\subsubsection*{Usage}
\begin{verbatim}
RandomSignal
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how one might generate an inspiral
waveform and compute its overlap with simulated detector noise
(with or witnout the signal present). The parameter
{\texttt randIn.type=0} generates only signal
{\texttt randIn.type=1} generates only noise
{\texttt randIn.type=2} generates {\texttt randIn.SignalAmp * signal(t)
+ randIn.NoiseAmp * noise(t)}
Note that one must calculate the length of the waveform and allocate memory for it
\emph{before} calling \texttt{InspiralWave}. The length of the waveform can 
be calculated by calling the function \texttt{InspiralWaveLength} beforehand, as shown.

There are only two functions which one can call to generate waveforms. These are
\texttt{InspiralWave},
which will return a \emph{single} waveform, and \texttt{InspiralWaveTemplates}, which
returns a \emph{pair}
of waveforms which have phases that differ by $\pi/2$.

\subsubsection*{Exit codes}
\input{RandomSignalCE}

\subsubsection*{Uses}
This code directly uses the following functions and macros (see those functions
to find out what they call in turn):
\begin{verbatim}
lalDebugLevel
LALEstimateFwdRealFFTPlan
LALEstimateInvRealFFTPlan
LALInspiralParameterCalc
LALInspiralWaveLength
LALInspiralWaveOverlap
LALNoiseSpectralDensity 
LALRandomInspiralSignal
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{RandomSignalCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="RandomSignalCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
 
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=1;
     
int 
main ( void ) 
{
   static LALStatus status;
   INT4 i, j;
   REAL8 dt, df, t0;
   REAL4Vector signal, correlation;
   RandomInspiralSignalIn randIn;
   OverlapIn overlapin;
   OverlapOut overlapout, overlapoutmax;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *RandomSignal;
/*---------------------------------------------------------------------------*/

   randIn.type = 2;
   randIn.SignalAmp = 10.;
   randIn.NoiseAmp = 1.0;
   randIn.useed = 610903;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=1000;
   randIn.mMin = 4.0;
   randIn.MMax = 10.0;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.ieta = 1; 
   randIn.param.fLower = 40.;
   randIn.param.fCutoff = 1000.;
   randIn.param.tSampling = 4000.;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 0;
   randIn.param.method = one;
   randIn.param.order = twoPN;
   randIn.param.domain = TimeDomain;
   randIn.param.approximant = pade;
   randIn.param.massChoice = m1Andm2;

   dt = 1./randIn.param.tSampling;
   t0 = 0.0;
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   correlation.length = signal.length;
   randIn.psd.length = signal.length/2 + 1;

   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &randIn.psd, &LALLIGOIPsd, df);

   overlapin.psd = randIn.psd;
/*------------------------ estimate the fft plans --------------------*/
   LALEstimateFwdRealFFTPlan(&status, &fwdp, signal.length);
   LALEstimateInvRealFFTPlan(&status, &revp, signal.length);

   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;

   i=1;
   while (i--) {
      LALRandomInspiralSignal(&status, &signal, &randIn);
      fprintf(stderr, "%d %e %e\n", i, randIn.param.mass1, randIn.param.mass2);
      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
      printf_timeseries (correlation.length, correlation.data, dt, t0) ;
/*
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
      printf_timeseries (correlation.length, correlation.data, dt, t0) ;
*/
     }

/*
      fprintf(stderr, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e\n",  
        randIn.param.t0, 
        randIn.param.t2, 
        randIn.param.t3, 
        randIn.param.t4, 
        randIn.param.mass1, 
        randIn.param.mass2, 
        randIn.param.totalMass, 
        overlapin.param.t0,
        overlapin.param.t2,
        overlapin.param.t3,
        overlapin.param.t4,
        overlapin.param.totalMass,
        overlapin.param.mass1,
        overlapin.param.mass2,
        overlapoutmax.phase,
        overlapoutmax.bin,
        overlapoutmax.max
      );

*/
/* destroy the plans, correlation and signal */
   LALFree(signal.data);
   LALFree(randIn.psd.data);
   LALFree(correlation.data);
   LALDestroyRealFFTPlan (&status, &fwdp);
   LALDestroyRealFFTPlan (&status, &revp);
   LALCheckMemoryLeaks();
   return(0);
}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
{
  int i=0;
  FILE *file;

  file = fopen("RandomSignal.out", "w");

  do 
     fprintf (file, "%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 

  fprintf(file, "&\n");
  fclose(file);
}

