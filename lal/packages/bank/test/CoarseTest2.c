/******************************** <lalVerbatim file="LALCoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALCoarseTest.c}}
\label{ss:LALCoarseTest.c}

Test code for the inspiral modules.

\subsubsection*{Usage}
\begin{verbatim}
LALCoarseTest
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how one might generate inspiral
waveforms and use them to compute the overlap of a random signal
(with or witnout simulated noise) of a certain strength. The parameter
{\texttt randIn.type=0} generates only signal
{\texttt randIn.type=1} generates only noise
{\texttt randIn.type=2} generates {\texttt randIn.SignalAmp * signal(t)
+ randIn.NoiseAmp * noise(t)}
Note that one must calculate the length of the waveform and allocate memory for it
\emph{before} calling
\texttt{InspiralWave}. The length of the waveform can be calculated by calling the function
\texttt{InspiralWaveLength} beforehand, as shown.

There are only two functions which one can call to generate waveforms. These are
\texttt{InspiralWave},
which will return a \emph{single} waveform, and \texttt{InspiralWaveTemplates}, which
returns a \emph{pair}
of waveforms which have phases which differ by $\pi/2$.

\subsubsection*{Exit codes}
\input{LALCoarseTestCE}

\subsubsection*{Uses}
\begin{verbatim}
laldebuglevel
InspiralWaveLength
InspiralWave
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALCoarseTestCV}}
******************************************************* </lalLaTeX> */


/***************************** <lalErrTable file="LALCoarseTestCE"> */
/***************************** </lalErrTable> */





#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
 
void printf_timeseries (INT4, REAL4 *, REAL8, REAL8);

INT4 lalDebugLevel=1;
     
int 
main ( void )
{
/* 
   It is assumed that the templateList is no longer than 50000; change it
   if you expect your lattice to be larger than that
*/
#if FIXME
   static InspiralTemplateList list[5000];
   static LALStatus status;
   InspiralCoarseBankIn coarseIn;
   INT4 i, j, nlist, pad;
   REAL8 dt, df, t0;
   REAL4Vector signal, correlation;
   RandomInspiralSignalIn randIn;
   Detector choice;
   OverlapIn overlapin;
   OverlapOut overlapout;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *CoarseOut;

   CoarseOut = fopen("CoarseTest2.out", "w");

/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   coarseIn.mMin = 4.0;
   coarseIn.MMax = 40.0;
   coarseIn.mmCoarse = 0.90;
   coarseIn.mmFine = 0.97;
   coarseIn.fLower = 40.;
   coarseIn.fUpper = 1000;
   coarseIn.iflso = 0;
   coarseIn.tSampling = 4000.;
   coarseIn.detector = virgo;
   coarseIn.method = one;
   coarseIn.order = twoPN;
   coarseIn.approximant = taylor;
   coarseIn.domain = TimeDomain;
   coarseIn.space = Tau0Tau2;

   randIn.type = 2;
   randIn.SignalAmp = 10.0;
   randIn.NoiseAmp = 1.0;
   randIn.useed = 99184045;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.0; 
   pad = randIn.param.nStartPad=1000;

/*---------------------------------------------------------------------------*/
/*CHANGE NOTHING BELOW THIS LINE IF YOU ARE UNSURE OF WHAT YOU ARE DOING*/
/*---------------------------------------------------------------------------*/
   LALInspiralCreateCoarseBank(&status, list, &nlist, coarseIn);
   fprintf(CoarseOut, "Number of Coarse Bank Templates=%d\n",nlist);
   for (i=0; i<nlist; i++) {
      fprintf(CoarseOut, "%e %e %e %e %e %e %e\n", 
         list[i].params.t0, 
         list[i].params.t2, 
         list[i].params.t3, 
         list[i].params.totalMass,
         list[i].params.eta, 
         list[i].params.mass1, 
         list[i].params.mass2);
   }
   fprintf(CoarseOut, "&\n");

   randIn.mMin = coarseIn.mMin;
   randIn.MMax = coarseIn.MMax;
   randIn.param.mass1=randIn.mMin;
   randIn.param.mass2=randIn.mMin;
   randIn.param.ieta=1; 
   randIn.param.fLower=coarseIn.fLower;
   randIn.param.fCutoff=coarseIn.fUpper;
   randIn.param.tSampling=coarseIn.tSampling;
   randIn.param.signalAmplitude=1.0;
   randIn.param.nEndPad=0;
   randIn.param.method=coarseIn.method;
   randIn.param.order=coarseIn.order;
   randIn.param.domain=coarseIn.domain;
   randIn.param.approximant=coarseIn.approximant;
   randIn.param.massChoice=m1Andm2;

   dt = 1./randIn.param.tSampling;
   t0 = 0.0;
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   correlation.length = signal.length;
   randIn.psd.length = signal.length/2 + 1;

   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   choice = virgo;
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &randIn.psd, choice, df);
/*
   printf_timeseries (randIn.psd.length, randIn.psd.data, df, t0);
   exit(0);
*/
   overlapin.psd = randIn.psd;
   /* Measure the plans */
   LALMeasureFwdRealFFTPlan(&status, &fwdp, signal.length);
   LALMeasureInvRealFFTPlan(&status, &revp, signal.length);
   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;
   i=5;
   while (i--) {
      LALRandomInspiralSignal(&status, &signal, &randIn);
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      for (j=0; j<nlist; j++) {
         REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
/*    
      The "distance" between the random signal and the template
      on the lattice can be computed using dt0 and dt1 and the metric
*/
         dt0 = -list[j].params.t0 + randIn.param.t0;
         dt1 = -list[j].params.t2 + randIn.param.t2;
/*   
      Compute the approximate match of the template and the
      random signal using the metric; if that is larger than about 
      50 % then compute the exact match by generating the template. 
      (Sometime in the future we may want the following six lines 
      to be replaced by a routine like this.
      ApproximateMatch (status, &match, list[j], dt0, dt1);)
*/ 
         g00=list[j].metric.g00; 
         g11=list[j].metric.g11; 
         ang=list[j].metric.theta;
         h00 = g00 * cos(ang)*cos(ang) + g11 * sin(ang)*sin(ang);
         h11 = g11 * cos(ang)*cos(ang) + g00 * sin(ang)*sin(ang);
         h01 = (g00-g11) * cos(ang)*sin(ang);
         match = 1. - h00*dt0*dt0 - h11*dt1*dt1 - 2.*h01*dt0*dt1;
         if (match > 0.8) {
            overlapin.param = list[j].params;
            LALInspiralWaveOverlap(&status, &correlation, &overlapout, &overlapin);
            fprintf(CoarseOut, "%e %e %d %e %e %e %e\n",  
                    randIn.param.t0, 
                    randIn.param.t2, 
                    i, 
                    overlapout.max, 
                    match, 
                    list[j].params.t0, 
                    list[j].params.t2
            );
/*
            printf_timeseries (correlation.length, correlation.data, dt, t0);
*/
         }
      }
   }
   fclose(CoarseOut);
   /* destroy the plans, correlation and signal */
   LALFree(signal.data);
   LALFree(correlation.data);
   LALFree(fwdp);
   LALFree(revp);
/*
   LALDestroyFFTPlan(&status,&fwdp);   
   LALDestroyFFTPlan(&status,&revp);
*/
   return(1);
#else
   return 77;
#endif
}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t_0) 
{
  int i=0;

  do 
     printf ("%e %e\n", i*delta+t_0, *(signal+i));

  while (n-++i); 
  printf("&\n");
}
