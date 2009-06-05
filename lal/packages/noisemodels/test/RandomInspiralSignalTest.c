/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
LALInspiralFindEventsCluster
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

NRCSID (RANDOMINSPIRALSIGNALTESTC,"$Id$");

void printf_timeseries(INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0, FILE *file);

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
   INT4 i, j;
   REAL8 dt;
   REAL8 df;
   REAL8 t0;
   REAL4Vector signal;
   REAL4Vector correlation;
   static RandomInspiralSignalIn randIn;
   InspiralWaveOverlapIn overlapin;
   InspiralWaveOverlapOut overlapout;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   InspiralEventsList *eventlist=NULL;

   InspiralFindEventsIn findeventsin;
   INT4 nEvents=0;

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
   randIn.SignalAmp = 10.L;
   randIn.NoiseAmp = 1.;
   randIn.useed = 83275488;
   randIn.mMin = 3.0;
   randIn.mMax = 20.0;
   randIn.MMax = 2.*randIn.mMax;
   randIn.param.startTime=0.0;
   randIn.param.startPhase=0.88189;
   randIn.param.nStartPad = 3000;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.ieta = 1;
   randIn.param.fLower = 40.;
   randIn.param.fCutoff = 2000.;
   randIn.param.tSampling = 4096.;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 0;
   randIn.param.OmegaS = 0.0;
   randIn.param.Theta = 0.0;
   randIn.param.approximant = EOB;
   randIn.param.order = twoPN;
   randIn.param.massChoice = totalMassAndEta;
   randIn.etaMin = randIn.mMin * ( randIn.MMax - randIn.mMin) /
      pow(randIn.MMax,2.);

   randIn.param.psi0 = 72639.;
   randIn.param.psi3 = -768.78;
   randIn.param.alpha = 0.766;
   randIn.param.fFinal = 331.4;

   randIn.psi0Min = 10000.L;
   randIn.psi0Max = 100000.L;
   randIn.psi3Min = -1000.L;
   randIn.psi3Max = 1000.L;

   findeventsin.currentGPSTime = 0;
   findeventsin.Threshold = 4.0;
   findeventsin.ClusterThreshold = 2.;
   findeventsin.displayPSD = 0;
   findeventsin.displayData = 0;
   findeventsin.displayTemplates = 0;
   findeventsin.displayCorrelation = 1;
   findeventsin.displayCorrelationStats = 0;
   findeventsin.nBegin = 0;
   findeventsin.nEnd = 0;

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

   findeventsin.fwdp = fwdp;
   findeventsin.revp = revp;
   findeventsin.psd = randIn.psd;

   file = fopen("RandomInspiralSignalTest.out", "w");

   randIn.param.approximant = PadeT1;
   i=4;
   while (i--) {
      randIn.type = 0;
      LALRandomInspiralSignal(&status, &signal, &randIn);
      fprintf(stderr, "%d %e %e\n", i, randIn.param.mass1, randIn.param.mass2);
      if (TimeDomain)
      {
	      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
      }
      else
      {
	      randIn.param.approximant = PadeT1;
	      LALInspiralParameterCalc(&status, &randIn.param);
	      overlapin.signal = signal;
	      randIn.param.fCutoff = 1./(pow(6.,1.5) * LAL_PI * randIn.param.totalMass * LAL_MTSUN_SI);
	      overlapin.param = randIn.param;
              overlapin.ifExtOutput = 0;
	      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
	      fprintf(stderr, "phase_max=%e bin_max=%d overlap_max=%e\n",
			      overlapout.phase, overlapout.bin, overlapout.max);
      }

      randIn.type = 1;
      LALRandomInspiralSignal(&status, &signal, &randIn);
      if (TimeDomain)
      {
	      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
      }
      else
      {
	      LALInspiralParameterCalc(&status, &randIn.param);
	      overlapin.signal = signal;
	      overlapin.param = randIn.param;
	      randIn.param.fCutoff = 1./(pow(6.,1.5) * LAL_PI * randIn.param.totalMass * LAL_MTSUN_SI);
              overlapin.ifExtOutput = 0;
	      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
	      fprintf(stderr, "phase_max=%e bin_max=%d overlap_max=%e\n",
			      overlapout.phase, overlapout.bin, overlapout.max);
      }

      randIn.type = 2;
      randIn.param.approximant = TaylorT1;
      randIn.param.massChoice = m1Andm2;
      randIn.param.nStartPad = 3000;
      LALRandomInspiralSignal(&status, &signal, &randIn);
      randIn.param.nStartPad = 0;
      if (TimeDomain)
      {
	      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
      }
      else
      {
	      if (randIn.param.approximant != BCV) LALInspiralParameterCalc(&status, &randIn.param);
	      overlapin.signal = signal;
	      overlapin.param = randIn.param;
	      if (randIn.param.approximant !=BCV)
		      randIn.param.fCutoff = 1./(pow(6.,1.5) * LAL_PI * randIn.param.totalMass * LAL_MTSUN_SI);
	      else
		      randIn.param.fCutoff = randIn.param.fFinal;
              overlapin.ifExtOutput = 0;
	      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
	      printf_timeseries (correlation.length, correlation.data, dt, t0, file) ;
	      fprintf(stderr, "m1=%e m2=%e t0=%e t2=%e psi0=%e psi3=%e phase_max=%e bin_max=%d overlap_max=%e\n",
			      randIn.param.mass1, randIn.param.mass2, randIn.param.t0, randIn.param.t2,
			      randIn.param.psi0, randIn.param.psi3,
			      overlapout.phase, overlapout.bin, overlapout.max);
	      findeventsin.signal = signal;
	      findeventsin.param = randIn.param;
	      LALInspiralFindEventsCluster (&status, &nEvents, &eventlist, &findeventsin);
	      fprintf(stderr, "Number of Events=%d\n", nEvents);
	      if (nEvents)
	      {

		      for (j=0; j<nEvents; j++)
		      {

			      fprintf(stderr, "%d %d %e %e\n",
					      eventlist[j].impulseTime, eventlist[j].endTime, eventlist[j].snr, eventlist[j].chisq);
		      }
		      LALFree(eventlist);
		      eventlist = NULL;
	      }
      }
   }


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

