/******************************** <lalVerbatim file="FilterTestCV">
Author: Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{FilterTest.c}}
\label{ss:FilterTest.c}

Test code for the inspiral bank modules.

\subsubsection*{Usage}
\begin{verbatim}
FilterTest
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
\input{FilterTestCE}

\subsubsection*{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\begin{verbatim}
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALRandomInspiralSignal
LALInspiralParameterCalc
LALNoiseSpectralDensity 
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALForwardRealFFT
LALReverseRealFFT
LALDestroyRealFFTPlan
LALInspiralWaveOverlap
LALInspiralParameterCalc
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FilterTestCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="FilterTestCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
 
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=1;
     
int 
main ( void ) 
{
   InspiralTemplateList *list=NULL;
   static LALStatus status;
   InspiralCoarseBankIn coarseIn;
   InspiralBankParams bankPars;
   INT4 ntrials, i, j, nlist, jmax;
   REAL8 dt, df, t0, omax;
   REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
   REAL4Vector signal, correlation;
   UINT4 numPSDpts = 1048576;
   void *noisemodel = LALLIGOIPsd;
        
   static RandomInspiralSignalIn randIn;
   InspiralWaveOverlapIn overlapin;
   InspiralWaveOverlapOut overlapout, overlapoutmax;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *FilterTest;
   FilterTest = fopen("FilterTest.out", "w");

   fprintf(stderr, "This test code does three things:\n");
   fprintf(stderr, "(a) Creates a filter bank at a certain minimal match\n");
   fprintf(stderr, "(b) Generates signals with random parmeters\n");
   fprintf(stderr, "(c) Filters each of these signals with the a \n");
   fprintf(stderr, "subset of templates close to the random signal \n");
   fprintf(stderr, "and reports the best SNR achived\n");
   fprintf(stderr, "Results of the run are written in FilterTest.out\n");

   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   coarseIn.mMin = 3.0;
   coarseIn.mMax = 10.0;
   coarseIn.MMax = 2.*coarseIn.mMax;
   coarseIn.massRange = MinMaxComponentMass;
   coarseIn.mmCoarse = 0.90;
   coarseIn.mmFine = 0.95;
   coarseIn.fLower = 40.;
   coarseIn.fUpper = 1000;
   coarseIn.iflso = 0;
   coarseIn.tSampling = 2048.;
   coarseIn.order = 4;
   coarseIn.approximant = TaylorT1;
   coarseIn.space = Tau0Tau3;
/* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) /
      pow(coarseIn.MMax,2.);
   /* fill the psd */
   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
   coarseIn.shf.f0 = 0;
   LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
   coarseIn.shf.deltaF = coarseIn.tSampling / (REAL8) coarseIn.shf.data->length;
   LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );

   fprintf(FilterTest, "#mMin mMax mmCoarse mmFine fLower fUpper tSampling order approximant space\n");
   fprintf(FilterTest, "#%e %e %e %e %e %e %e %d %d %d\n", 
      coarseIn.mMin,
      coarseIn.MMax,
      coarseIn.mmCoarse,
      coarseIn.mmFine,
      coarseIn.fLower,
      coarseIn.fUpper,
      coarseIn.tSampling,
      coarseIn.order,
      coarseIn.approximant,
      coarseIn.space
   );
   fflush(stdout);
   randIn.type = 2;
   randIn.SignalAmp = 10.;
   randIn.NoiseAmp = 1.e0;
   randIn.useed = 218092;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=2000;

   fprintf(FilterTest, "#type SignalAmp NoiseAmp startTime startPhase startpad\n");
   fprintf(FilterTest, "#%d %e %e %e %e %d\n",
      randIn.type,
      randIn.SignalAmp,
      randIn.NoiseAmp,
      randIn.param.startTime,
      randIn.param.startPhase,
      randIn.param.nStartPad
   );

/*---------------------------------------------------------------------------*/
/*CHANGE NOTHING BELOW THIS LINE IF YOU ARE UNSURE OF WHAT YOU ARE DOING     */
/*---------------------------------------------------------------------------*/
   LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);
/* REPORTSTATUS(&status); */
   if (nlist==0) exit(0);

   fprintf(stderr, "#Number of Coarse Bank Templates=%d\n",nlist);
   fprintf(FilterTest, "Number of Coarse Bank Templates=%d\n",nlist);

   for (i=0; i<nlist; i++) {
      fprintf(FilterTest, "%e %e %e %e %e %e %e %e\n", 
         list[i].params.t0, 
         list[i].params.t2, 
         list[i].params.t3, 
         list[i].params.t4, 
         list[i].params.totalMass,
         list[i].params.eta, 
         list[i].params.mass1, 
         list[i].params.mass2);
   }
   fprintf(FilterTest, "&\n");
   LALInspiralSetSearchLimits(&status, &bankPars, coarseIn);
   randIn.mMin = coarseIn.mMin;
   randIn.mMax = coarseIn.mMax; 
   randIn.MMax = coarseIn.MMax;
   randIn.t0Min = bankPars.x0Min;
   randIn.t0Max = bankPars.x0Max;
   randIn.tnMin = bankPars.x1Min;
   randIn.tnMax = bankPars.x1Max;
   randIn.etaMin = coarseIn.etamin;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.ieta = 1; 
   randIn.param.fLower = coarseIn.fLower;
   randIn.param.fCutoff = coarseIn.fUpper;
   randIn.param.tSampling = coarseIn.tSampling;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 2000;
   randIn.param.order = 6;
   randIn.param.approximant = PadeT1;
   randIn.param.massChoice = t03;
   if (randIn.param.massChoice != m1Andm2) 
   {
      if (coarseIn.space==Tau0Tau2) 
         randIn.param.massChoice = t02;
      else
         randIn.param.massChoice = t03;
   }
   dt = 1./randIn.param.tSampling;
   t0 = 0.0;
   signal.length = 0.;
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   fprintf(FilterTest, "signal length = %d\n", signal.length);
/* REPORTSTATUS(&status); */
   correlation.length = signal.length;
   randIn.psd.length = signal.length/2 + 1;

   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &randIn.psd, noisemodel, df);
/* REPORTSTATUS(&status); */

   overlapin.psd = randIn.psd;
/*--------------------------
   Estimate the plans 
--------------------------*/
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
/* REPORTSTATUS(&status); */
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);
/* REPORTSTATUS(&status); */
   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;
   ntrials=5;
   fprintf(FilterTest, "#mass1 mass2 totalMass t0 t2 t3 t4 mass1 mass2 totalMass t0 t2 t3 t4 omax \n");
   fprintf(FilterTest, "#Signal Length=%d Number of sims=%d\n", signal.length, i);
   fprintf(stderr,"----------------------------------------------\n");
   fprintf(stderr, "   t0             t2        Overlap/SNR\n");
   fprintf(stderr,"----------------------------------------------\n");
   fprintf(FilterTest, "   t0            t2        Overlap/SNR\n");
   printf("t0Min=%e, t0Max=%e, tnMin=%e, tnMax=%e\n",randIn.t0Min,randIn.t0Max,randIn.tnMin,randIn.tnMax);
   fflush(stdout);
   while (ntrials--) {
      LALRandomInspiralSignal(&status, &signal, &randIn);
/*
      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
      for (j=0; j<signal.length; j++) printf("%e\n", correlation.data[j]);
      break;
*/
      LALInspiralParameterCalc(&status, &randIn.param);
      overlapin.signal = signal;
      jmax = 0;
      omax = 0.0;
      for (j=0; j<nlist; j++) {
/*-----------------------------------------------------------------------   
   The "distance" between the random signal and the template
   on the lattice can be computed using dt0 and dt1 and the metric
------------------------------------------------------------------------*/
         dt0 = -list[j].params.t0 + randIn.param.t0;
         dt1 = -list[j].params.t2 + randIn.param.t2;
/*-----------------------------------------------------------------------   
   Compute the approximate match of the template and the
   random signal using the metric; if that is larger than about 
   50 % then compute the exact match by generating the template. 
   (Sometime in the future we may want the following six lines 
   to be replaced by a routine like this.
   ApproximateMatch (status, &match, list[j], dt0, dt1);)
------------------------------------------------------------------------*/
         g00=list[j].metric.g00; 
         g11=list[j].metric.g11; 
         ang=abs(list[j].metric.theta);
         h00 = g00 * cos(ang)*cos(ang) + g11 * sin(ang)*sin(ang);
         h11 = g11 * cos(ang)*cos(ang) + g00 * sin(ang)*sin(ang);
         h01 = (g00-g11) * cos(ang)*sin(ang);
         match = 1. - h00*dt0*dt0 - h11*dt1*dt1 - 2.*h01*dt0*dt1;
         if (match > coarseIn.mmCoarse/2.) {
           overlapin.param = list[j].params;
           LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
           if (omax < overlapout.max) {
                omax = overlapout.max;
                overlapoutmax = overlapout;
                jmax = j;
            }
         }
     }
     printf("%e %e %e %e %e %e %e %e\n", randIn.param.t0, randIn.param.t2, randIn.param.t3, randIn.param.mass1, randIn.param.mass2, randIn.param.totalMass, randIn.param.eta, omax*signal.length);
     fprintf(FilterTest,"%e %e %e %e %e %e %e %e\n", randIn.param.t0, randIn.param.t2, randIn.param.t3, randIn.param.mass1, randIn.param.mass2, randIn.param.totalMass, randIn.param.eta, omax*signal.length);
   fflush(stdout);
   }
   fclose(FilterTest);

   /* destroy the plans, correlation and signal */

   if (signal.data != NULL) LALFree(signal.data);
   if (correlation.data != NULL) LALFree(correlation.data);
   if (randIn.psd.data != NULL) LALFree(randIn.psd.data);
   if (list!= NULL) LALFree(list);
   LALDestroyRealFFTPlan(&status,&fwdp);   
   LALDestroyRealFFTPlan(&status,&revp);
   LALDDestroyVector( &status, &(coarseIn.shf.data) );

   LALCheckMemoryLeaks();
   return(0);
}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
{
  int i=0;

  do 
     printf ("%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 
  printf("&\n");
}
