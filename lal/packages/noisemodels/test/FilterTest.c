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
#include <lal/SeqFactories.h>
 
void  LALInspiralCreateFlatBank(LALStatus *status, REAL4VectorSequence *list, InspiralBankParams *bankParams);
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=1;
     
int 
main ( void ) 
{
   REAL4VectorSequence *list=NULL;
   static LALStatus status;
   static InspiralBankParams bankParams;
   INT4 ntrials, i, j, nlist, jmax;
   REAL8 df, omax;
   REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
   REAL4Vector signal, correlation;
   UINT4 numPSDpts = 1048576;
   void *noisemodel = LALLIGOIPsd;
   static RandomInspiralSignalIn randIn;
   static InspiralWaveOverlapIn overlapin;
   static InspiralWaveOverlapOut overlapout, overlapoutmax;
   static CreateVectorSequenceIn in; 
  REAL8FrequencySeries shf;
   static InspiralMetric metric;
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

   bankParams.x0Min = 0.0;
   bankParams.x0Max = 2.5e5;
   bankParams.x1Min = -2.2e3;
   bankParams.x1Max = 8.e2;
   bankParams.minimalMatch = 0.95;

   in.length = 1;
   in.vectorLength = 2;
   LALSCreateVectorSequence(&status, &list, &in);
   list->vectorLength = 2;
/*---------------------------------------------------------------------------*/
   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   randIn.type = 0;
   randIn.SignalAmp = 10.0;
   randIn.NoiseAmp = 1.0;
   randIn.useed = 128092;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=2000;
   randIn.psi0Min = bankParams.x0Min;
   randIn.psi0Max = bankParams.x0Max;
   randIn.psi3Max = bankParams.x1Max;
   randIn.psi3Min = bankParams.x1Min;
   randIn.param.fLower = 40.;
   randIn.param.fCutoff = 1000.;
   randIn.param.tSampling = 2048;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 2000;
   randIn.param.order = 4;
   randIn.param.approximant = PadeT1;
   randIn.param.massChoice = m1Andm2;
   randIn.mMin = 5.0;
   randIn.mMax = 20.0;
   randIn.MMax = 2.*randIn.mMax;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.fendBCV = 300.;
   randIn.param.alpha = 0.;

   signal.length = 0.;
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   fprintf(FilterTest, "signal length = %d\n", signal.length);
/* REPORTSTATUS(&status); */
   correlation.length = signal.length;
   randIn.psd.length = signal.length/2 + 1;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
   shf.f0 = 0;
   LALDCreateVector( &status, &(shf.data), randIn.psd.length );
   shf.deltaF = randIn.param.tSampling / signal.length;
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, shf.data, noisemodel, df);
   LALInspiralComputeMetricBCV(&status, &metric, &shf, &randIn.param);
   bankParams.metric = &metric;
   LALInspiralCreateFlatBank(&status, list, &bankParams);
   nlist = list->length;
/* REPORTSTATUS(&status); */
   if (nlist==0) exit(0);

   fprintf(stderr, "#Number of Coarse Bank Templates=%d\n",nlist);
   fprintf(FilterTest, "Number of Coarse Bank Templates=%d\n",nlist);

   for (i=0; i<nlist; i++)
   {
	   int i2;
	   i2=i*2;
	   fprintf(FilterTest, "%e %e\n", list->data[i2], list->data[i2+1]);
   }

/* REPORTSTATUS(&status); */

   randIn.psd = *shf.data;
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
   ntrials=10;
   fprintf(FilterTest, "#Signal Length=%d Number of sims=%d\n", signal.length, i);
   fprintf(stderr,"----------------------------------------------\n");
   fprintf(stderr, "   psi0             psi3        Overlap/SNR\n");
   fprintf(stderr,"----------------------------------------------\n");
   fprintf(FilterTest, "   psi0            psi3        Overlap/SNR\n");
   printf("psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n",
		   randIn.psi0Min,randIn.psi0Max,randIn.psi3Min,randIn.psi3Max);
   while (ntrials--) {
	      
      int numTemplates=0;
      randIn.type = 0;
      /*
      randIn.param.massChoice = m1Andm2;
      */
      randIn.param.massChoice = psi0psi3;
      randIn.param.approximant = BCV;
      LALRandomInspiralSignal(&status, &signal, &randIn);
/*
      LALREAL4VectorFFT(&status, &correlation, &signal, revp);
      for (j=0; j<signal.length; j++) printf("%e\n", correlation.data[j]);
      break;
      LALInspiralParameterCalc(&status, &randIn.param);
*/
      randIn.param.massChoice = m1Andm2;
      LALInspiralParameterCalc(&status, &randIn.param);
      /* We shall choose the cutoff at r=3M for EOB waveforms */
      switch (randIn.param.approximant)
      {
	   case EOB:
	      randIn.param.fCutoff = 1./(pow(3.,1.5) * LAL_PI * randIn.param.totalMass * LAL_MTSUN_SI);
	      break;
	   case BCV:
	      randIn.param.fCutoff = randIn.param.fendBCV;
	      break;
	   default:
	      randIn.param.fendBCV = randIn.param.fCutoff;
	      break;
      }
      randIn.param.fendBCV = randIn.param.fCutoff;

      randIn.param.approximant = BCV;
      overlapin.signal = signal;
      overlapin.param = randIn.param;
      jmax = 0;
      omax = 0.0;
      for (j=0; j<nlist; j++) {
	      int j2;
	      j2 = j*2;
/*-----------------------------------------------------------------------   
   The "distance" between the random signal and the template
   on the lattice can be computed using dt0 and dt1 and the metric
------------------------------------------------------------------------*/
         dt0 = -list->data[j2] + randIn.param.psi0;
         dt1 = -list->data[j2+1] + randIn.param.psi3;
/*-----------------------------------------------------------------------   
   Compute the approximate match of the template and the
   random signal using the metric; if that is larger than about 
   50 % then compute the exact match by generating the template. 
   (Sometime in the future we may want the following six lines 
   to be replaced by a routine like this.
   ApproximateMatch (status, &match, list[j], dt0, dt1);)
------------------------------------------------------------------------*/
         match = 1. - metric.G00*dt0*dt0 - metric.G11*dt1*dt1 - 2.*metric.G01*dt0*dt1;
         /* if (match > bankParams.minimalMatch/2.) */
	 {
		 numTemplates++;
		 overlapin.param.psi0 = list->data[j2];
		 overlapin.param.psi3 = list->data[j2+1];
		 LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
		 if (omax < overlapout.max) {
			 omax = overlapout.max;
			 overlapoutmax = overlapout;
			 jmax = j;
		 }
         }
     }
     j = 2*jmax;
     printf("%e %e %e %e %e %d\n", list->data[j], list->data[j+1], randIn.param.mass1, randIn.param.mass2, omax, numTemplates);
     fprintf(FilterTest,"%e %e %e %e %e %d\n", list->data[j], list->data[j+1], randIn.param.mass1, randIn.param.mass2, omax, numTemplates);
     /*
     printf("%e %e %e %e %e %d\n", list->data[j], list->data[j+1], randIn.param.psi0, randIn.param.psi3, omax, numTemplates);
     fprintf(FilterTest,"%e %e %e %e %e %d\n", list->data[j], list->data[j+1], randIn.param.psi0, randIn.param.psi3, omax, numTemplates);
     */
   fflush(stdout);
   }
   fclose(FilterTest);

   /* destroy the plans, correlation and signal */

   /*
   LALDDestroyVector( &status, &(shf.data) );
   if (signal.data != NULL) LALFree(signal.data);
   if (correlation.data != NULL) LALFree(correlation.data);
   if (randIn.psd.data != NULL) LALFree(randIn.psd.data);
   if (list!= NULL) LALFree(list);
   LALDestroyRealFFTPlan(&status,&fwdp);   
   LALDestroyRealFFTPlan(&status,&revp);

   LALCheckMemoryLeaks();
   */
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
