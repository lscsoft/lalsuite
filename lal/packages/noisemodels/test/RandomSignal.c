/******************************** <lalVerbatim file="LALRandomSignalCV">
Author: Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALRandomSignal.c}}
\label{ss:LALRandomSignal.c}

Test code for the inspiral bank modules.

\subsubsection*{Usage}
\begin{verbatim}
LALRandomSignal
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
\input{LALRandomSignalCE}

\subsubsection*{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\begin{verbatim}
laldebuglevel
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALRandomInspiralSignal
LALInspiralParameterCalc
LALNoiseSpectralDensity 
LALMeasureFwdRealFFTPlan
LALMeasureInvRealFFTPlan
LALInspiralWaveOverlap
LALInspiralParameterCalc
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRandomSignalCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="LALRandomSignalCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
 
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=1;
     
int 
main ( void ) 
{
/* 
   It is assumed that InspiralTemplateList is no larger than 5000; 
   change it if you expect your lattice to be larger than that
*/
   static InspiralTemplateList list[5000];
   static LALStatus status;
   InspiralCoarseBankIn coarseIn;
   INT4 i, j, nlist, jmax;
   REAL8 dt, df, t0, omax;
   REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
   REAL4Vector signal, correlation;
   RandomInspiralSignalIn randIn;
   Detector choice;
   OverlapIn overlapin;
   OverlapOut overlapout, overlapoutmax;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *RandomSignal;
   RandomSignal = fopen("RandomSignal.out", "w");

/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   coarseIn.mMin = 1.0;
   coarseIn.MMax = 40.0;
   coarseIn.mmCoarse = 0.90;
   coarseIn.mmFine = 0.97;
   coarseIn.fLower = 40.;
   coarseIn.fUpper = 1000;
   coarseIn.iflso = 0;
   coarseIn.tSampling = 4000.;
   coarseIn.detector = ligo;
   coarseIn.method = one;
   coarseIn.order = twoPN;
   coarseIn.approximant = pade;
   coarseIn.domain = TimeDomain;
   coarseIn.space = Tau0Tau2;
/* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) /
      pow(coarseIn.MMax,2.);

   fprintf(RandomSignal, "#mMin mMax mmCoarse mmFine fLower fUpper tSampling detector method order approximant domain space\n");
   fprintf(RandomSignal, "#%e %e %e %e %e %e %e %d %d %d %d %d %d\n", 
      coarseIn.mMin,
      coarseIn.MMax,
      coarseIn.mmCoarse,
      coarseIn.mmFine,
      coarseIn.fLower,
      coarseIn.fUpper,
      coarseIn.tSampling,
      coarseIn.detector,
      coarseIn.method,
      coarseIn.order,
      coarseIn.approximant,
      coarseIn.domain,
      coarseIn.space
   );

   randIn.type = 0;
   randIn.SignalAmp = 10.;
   randIn.NoiseAmp = 1.0;
   randIn.useed = 610903;
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=10000;

   fprintf(RandomSignal, "#type SignalAmp NoiseAmp startTime startPhase startpad\n");
   fprintf(RandomSignal, "%d %e %e %e %e %d\n",
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
   LALInspiralCreateCoarseBank(&status, list, &nlist, coarseIn);
   fprintf(RandomSignal, "#Number of Coarse Bank Templates=%d\n",nlist);
   fprintf(stderr, "Number of Coarse Bank Templates=%d\n",nlist);
   for (i=0; i<nlist; i++) {
      fprintf(RandomSignal, "%e %e %e %e %e %e %e %e\n", 
         list[i].params.t0, 
         list[i].params.t2, 
         list[i].params.t3, 
         list[i].params.t4, 
         list[i].params.totalMass,
         list[i].params.eta, 
         list[i].params.mass1, 
         list[i].params.mass2);
   }
   fprintf(RandomSignal, "&\n");
   randIn.mMin = coarseIn.mMin;
   randIn.MMax = coarseIn.MMax;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;
   randIn.param.ieta = 1; 
   randIn.param.fLower = coarseIn.fLower;
   randIn.param.fCutoff = coarseIn.fUpper;
   randIn.param.tSampling = coarseIn.tSampling;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 0;
   randIn.param.method = coarseIn.method;
   randIn.param.order = coarseIn.order;
   randIn.param.domain = coarseIn.domain;
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
   choice = coarseIn.detector;
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, &randIn.psd, choice, df);

   overlapin.psd = randIn.psd;
/*--------------------------
   Measure the plans 
--------------------------*/
   LALMeasureFwdRealFFTPlan(&status, &fwdp, signal.length);
   LALMeasureInvRealFFTPlan(&status, &revp, signal.length);
   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;
   i=10;
   fprintf(RandomSignal, "#mass1 mass2 totalMass t0 t2 t3 t4 mass1 mass2 totalMass t0 t2 t3 t4 omax \n");
   fprintf(RandomSignal, "#Signal Length=%d Number of sims=%d\n", signal.length, i);
   while (i--) {
      LALRandomInspiralSignal(&status, &signal, &randIn);
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
         ang=list[j].metric.theta;
         h00 = g00 * cos(ang)*cos(ang) + g11 * sin(ang)*sin(ang);
         h11 = g11 * cos(ang)*cos(ang) + g00 * sin(ang)*sin(ang);
         h01 = (g00-g11) * cos(ang)*sin(ang);
         match = 1. - h00*dt0*dt0 - h11*dt1*dt1 - 2.*h01*dt0*dt1;
         if (match > 0.9) {
           overlapin.param = list[j].params;
           LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
            if (omax < overlapout.max) {
                omax = overlapout.max;
                overlapoutmax = overlapout;
                jmax = j;
            }
         }
     }

/*
     printf_timeseries (correlation.length, correlation.data, dt, t0);
     exit(0);
*/

      fprintf(RandomSignal, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e\n",  
        randIn.param.t0, 
        randIn.param.t2, 
        randIn.param.t3, 
        randIn.param.t4, 
        randIn.param.mass1, 
        randIn.param.mass2, 
        randIn.param.totalMass, 
        list[jmax].params.mass1, 
        list[jmax].params.mass2,
        list[jmax].params.totalMass,
        list[jmax].params.t0, 
        list[jmax].params.t2,
        list[jmax].params.t3,
        list[jmax].params.t4,
        overlapoutmax.phase,
        overlapoutmax.bin,
        overlapoutmax.max
      );
      fprintf(stderr, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e\n",  
        randIn.param.t0, 
        randIn.param.t2, 
        randIn.param.t3, 
        randIn.param.t4, 
        randIn.param.mass1, 
        randIn.param.mass2, 
        randIn.param.totalMass, 
        list[jmax].params.mass1, 
        list[jmax].params.mass2,
        list[jmax].params.totalMass,
        list[jmax].params.t0, 
        list[jmax].params.t2,
        list[jmax].params.t3,
        list[jmax].params.t4,
        overlapoutmax.phase,
        overlapoutmax.bin,
        overlapoutmax.max
      );
   }
   fclose(RandomSignal);
   /* destroy the plans, correlation and signal */
   LALFree(signal.data);
   LALFree(correlation.data);
   LALFree(fwdp);
   LALFree(revp);
/*
   LALDestroyFFTPlan(&status,&fwdp);   
   LALDestroyFFTPlan(&status,&revp);
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
