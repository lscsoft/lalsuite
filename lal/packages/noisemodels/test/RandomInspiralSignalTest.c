/******************************** <lalVerbatim file="RandomInspiralSignalTestCV">
Author: Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{RandomInspiralSignalTest.c}}
\label{ss:RandomInspiralSignalTest.c}

Test code for the inspiral wave generation and noisemodels modules.

\subsubsection*{Usage}
\begin{verbatim}
RandomInspiralSignalTest
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
\input{RandomInspiralSignalTestCE}

\subsubsection*{Uses}
This code directly uses the following functions and macros (see those functions
to find out what they call in turn):
\begin{verbatim}
lalDebugLevel
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALInspiralParameterCalc
LALInspiralWaveLength
LALInspiralWaveOverlap
LALNoiseSpectralDensity 
LALRandomInspiralSignal
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{RandomInspiralSignalTestCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="RandomInspiralSignalTestCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
 
void printf_timeseries(INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0, FILE *file);

INT4 lalDebugLevel=0;
     
int 
main ( void ) 
{
   static LALStatus status;
   INT4 i;
   REAL8 dt;
   REAL8 df;
   REAL8 t0;
   REAL4Vector signal;
   REAL4Vector correlation;
   RandomInspiralSignalIn randIn;
   InspiralWaveOverlapIn overlapin;
   InspiralWaveOverlapOut overlapout;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *file;
/*---------------------------------------------------------------------------*/

   fprintf(stderr, "This test code does the following things:\n");
   fprintf(stderr, "(a) generates a signal and computes its correlation\n");
   fprintf(stderr, "  with two orthogonal templates with same parameter \n");
   fprintf(stderr, "  values as the signal\n");
   fprintf(stderr, "(b) generates noise expected in an interferometer and\n");
   fprintf(stderr, "  computes its correlation with two arbitrary \n");
   fprintf(stderr, "  orthogonal templates \n");
   fprintf(stderr, "(c) generates signal + noise expected in an \n");
   fprintf(stderr, "  interferometer and computes its correlation \n");
   fprintf(stderr, "  with two orthogonal templates with same\n");
   fprintf(stderr, "  parametter values as the signal\n");
   fprintf(stderr, "All correlations are written in RandomInspiralSignalTest.out in a\n");
   fprintf(stderr, "format suitable for display with xmgr/xgrace\n");

   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
   randIn.SignalAmp = 10.;
   randIn.NoiseAmp = 1.0;
   randIn.useed = 92048;
   randIn.mMin = 4.0;
   /* randIn.mMax = 5.0; */
   randIn.MMax = 10.0;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=1000;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.ieta = 1; 
   randIn.param.fLower = 40.;
   randIn.param.fCutoff = 1000.;
   randIn.param.tSampling = 4000.;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 0;
   randIn.param.OmegaS = 0.0;
   randIn.param.Theta = 0.0;
   randIn.param.order = twoPN;
   randIn.param.approximant = TaylorT2;
   randIn.param.massChoice = m1Andm2;
   randIn.etaMin = randIn.mMin * ( randIn.MMax - randIn.mMin) /
      pow(randIn.MMax,2.);

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
/*------------------------ create the fft plans --------------------*/
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);

   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;
   file = fopen("RandomInspiralSignalTest.out", "w");

   i=1;
   while (i--) {
      randIn.type = 0;
      LALRandomInspiralSignal(&status, &signal, &randIn);
/*
      fprintf(stderr, "%d %e %e\n", i, randIn.param.mass1, randIn.param.mass2);
      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
*/
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;

      randIn.type = 1;
      LALRandomInspiralSignal(&status, &signal, &randIn);
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;

      randIn.type = 2;
      LALRandomInspiralSignal(&status, &signal, &randIn);
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
   }

   fprintf(stderr, "Parameters of the signal and max correlation obtained in case (c)\n\n");
   fprintf(stderr, "m1=%e m2=%e t0=%e t2=%e \n phase_max=%e bin_max=%d overlap_max=%e\n",  
      randIn.param.mass1, 
      randIn.param.mass2, 
      randIn.param.t0, 
      randIn.param.t2, 
      overlapout.phase,
      overlapout.bin,
      overlapout.max
   );

/* destroy the plans, correlation and signal */

   LALFree(signal.data);
   LALFree(randIn.psd.data);
   LALFree(correlation.data);
   LALDestroyRealFFTPlan (&status, &fwdp);
   LALDestroyRealFFTPlan (&status, &revp);
   LALCheckMemoryLeaks();
   fclose(file);
   return(0);
}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0, FILE *file) 
{
   int i=0;
 
 
   do 
      fprintf (file, "%e %e\n", i*delta+t0, *(signal+i));
   while (n-++i); 
 
   fprintf(file, "&\n");
}

