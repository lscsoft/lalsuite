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
(with or witnout the signal present).  Recently, the code
was modified so that users can output either the time-domain signal,
noise, or signal+noise and not the correlated output. This is
done via the variable \texttt{TimeDomain} which controls the type of output
If \texttt{TimeDomain=0} then the code fitlers the random signal
with a template of the same parameters and outputs the
results of the correlation.
If \texttt{TimeDomain=1} then the code outputs the time-domain
signal/noise/signal+noise 

The parameter
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
 
INT4 lalDebugLevel=1;
     
int 
main ( void ) 
{
   UINT4 TimeDomain=0;
/*
 * The int TimeDomain is used to control the type of output
 * If TimeDomain=0 then the code fitlers the random signal
 * with a template of the same parameters and outputs the
 * results of the correlation.
 * If TimeDomain=1 then the code outputs the time-domain
 * signal/noise/signal+noise
 */

   static LALStatus status;
   INT4 i, j ;
   REAL8 dt;
   REAL8 df;
   REAL8 t0;
   REAL4Vector signal;
   REAL4Vector correlation;
   REAL8Vector psd;
   REAL8 mMin, mMax, MMax;
   UINT4 useed;

   InspiralTemplate param;
   InspiralWaveOverlapIn overlapin;
   InspiralWaveOverlapOut overlapout;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   REAL8 chirpMass, spin1Frac, spin2Frac, mass1Sq, mass2Sq, spin1Theta, spin2Theta, spin1Phi, spin2Phi;
   FILE *file;
/*---------------------------------------------------------------------------*/

   fprintf(stderr, "This test code does the following things:\n");
   fprintf(stderr, "(a) generates a signal from a spinning BH binary and \n");
   fprintf(stderr, "    computes its correlationwith two orthogonal templates \n");
   fprintf(stderr, "    with same parameter values as the signal\n");

   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
   useed = 92048;
   mMin = 5.0;
   mMax = 10.0;
   MMax = mMax*2.;
   param.startTime=0.0; 
   param.startPhase=0.88189; 
   param.nStartPad=1000;
   param.mass1 = mMin;
   param.mass2 = mMin;
   param.ieta = 1; 
   param.fLower = 40.;
   param.fCutoff = 1000.;
   param.tSampling = 4000.;
   param.signalAmplitude = 1.0;
   param.nEndPad = 0;
   param.OmegaS = 0.0;
   param.Theta = 0.0;
   param.order = twoPN;
   param.approximant = SpinTaylorT3;
   param.massChoice = m1Andm2;
   param.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;

   srandom(useed);
   spin1Frac = 0.9;
   spin2Frac = 0.9;
   dt = 1./param.tSampling;
   t0 = 0.0;
   param.approximant = TaylorT3;
   LALInspiralWaveLength (&status, &signal.length, param);
   LALInspiralParameterCalc(&status, &param);
   fprintf(stderr, "#signal length=%d\n", signal.length);
   param.approximant = SpinTaylorT3;
   correlation.length = signal.length;
   psd.length = signal.length/2 + 1;

   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*psd.length);
   df = param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &psd, &LALLIGOIPsd, df);

   overlapin.psd = psd;
/*------------------------ create the fft plans --------------------*/
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);

   overlapin.fwdp = fwdp;
   overlapin.revp = revp;
   file = fopen("SpinOverlaps.out", "a");
      
   param.sourceTheta =  LAL_PI/4.L;
   param.sourcePhi = LAL_PI/5.L;

   param.orbitTheta0 = LAL_PI/3.L;
   param.orbitPhi0 = LAL_PI/2.L;

   spin1Theta = LAL_PI/4.;
   spin1Phi = LAL_PI/8.;

   spin2Theta = -LAL_PI/4.;
   spin2Phi = LAL_PI/4.;
		   
      
   i=4;
   while (i--) {
      REAL8 chirpMass, e1, e2, norm;
   
      e1 = random()/(float)RAND_MAX;
      e2 = random()/(float)RAND_MAX;
		   
      param.mass1 = mMin + (mMax - mMin) * e1;
      param.mass2 = mMin + (mMax - mMin) * e2;
   
      mass1Sq = pow(param.mass1*LAL_MTSUN_SI,2.L);
      mass2Sq = pow(param.mass2*LAL_MTSUN_SI,2.L);
      param.spin1[0] =  mass1Sq * spin1Frac * sin(spin1Theta) * cos(spin1Phi);
      param.spin1[1] =  mass1Sq * spin1Frac * sin(spin1Theta) * sin(spin1Phi);
      param.spin1[2] =  mass1Sq * spin1Frac * cos(spin1Theta);
      param.spin2[0] =  mass2Sq * spin2Frac * sin(spin2Theta) * cos(spin2Phi);
      param.spin2[1] =  mass2Sq * spin2Frac * sin(spin2Theta) * sin(spin2Phi);
      param.spin2[2] =  mass2Sq * spin2Frac * cos(spin2Theta);
      
      param.approximant = SpinTaylorT3;
      LALInspiralParameterCalc(&status, &param);
      LALInspiralWave(&status, &correlation, &param);
      /*
      for (j=0; j<correlation.length; j++) printf("%e\n", correlation.data[j]);
      printf("&\n");
      */
      LALREAL4VectorFFT(&status, &signal, &correlation, fwdp);
      LALInspiralWaveNormalise(&status, &signal, &norm, psd);
      chirpMass = pow(param.eta,0.6L) * (param.mass1+param.mass2);

      param.approximant = TaylorT3;
      overlapin.signal = signal;
      overlapin.param = param;
      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
      fprintf(stderr, "%e %e %e %d\n", chirpMass, overlapout.max*signal.length, overlapout.phase, overlapout.bin);
      fprintf(file, "%e %e %e %d\n", chirpMass, overlapout.max*signal.length, overlapout.phase, overlapout.bin);

   }
   fprintf(file, "&\n");


/* destroy the plans, correlation and signal */

   LALFree(signal.data);
   LALFree(psd.data);
   LALFree(correlation.data);
   LALDestroyRealFFTPlan (&status, &fwdp);
   LALDestroyRealFFTPlan (&status, &revp);
   LALCheckMemoryLeaks();
   fclose(file);
   return(0);
}
