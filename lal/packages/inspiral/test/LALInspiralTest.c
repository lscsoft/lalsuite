/*  <lalVerbatim file="LALInspiralTestCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTest.c}}
Test routine for wave generation codes. First set the 
\texttt{InspiralTemplate} structure (example given below). Call the function
\begin{verbatim}
LALInspiralWaveLength (&status, &n, params)
\end{verbatim}
to measure the length of the array required which will be returned in 
\texttt{n}. Then call the function 
\begin{verbatim}
LALInspiralWave (&status, &signal, &params)
\end{verbatim}
to generate the wave, which will be returned in \texttt{signal}.
Example values of the parameters that can be set (with options in 
brackets):

\begin{verbatim}
   params.ieta=1;          (1 for comparable masses model, 0 for test mass model)
   params.mass1=1.4;       (masses of the component stars in solar masses) 
   params.mass2=1.4; 
   params.startTime=0.0;   (defined so that the instantaneous GW frequency 
                            is params.fLower at params.startTime)
   params.startPhase=0.0;  (0 to $\pi/2$)
   params.fLower=40.0;     (in Hz)
   params.fCutoff=1000.0;  (in Hz)
   params.tSampling=4000.; (in Hz; should be larger than $2\times \rm fCutoff$
                            or $2\times f_{\rm lso}$, whichever is smaller)
   params.signalAmplitude=1.0; 
   params.nStartPad=0;     (number of leading zero bins)
   params.nEndPad=0;       (number of trailing zero bins)
   params.method=one;      (the method in the convention of DIS 2001; 
                            one=ODE solver, 
                            two=implicit phasing formula solved in quadrature, 
                            three=explicit time-domain phasing)
                            eob=effective-one-body approach
   params.order=twoPN;     (also Newtonian, onePN, oneAndHalfPN, twoAndHalfPN)
   params.domain=TimeDomain;
   params.approximant=pade;(also taylor)
   params.massChoice=m1Andm2; (also t0t2, t0t3, t0t4, totalMassAndEta,totalMassAndMu) 
\end{verbatim}

\vfill{\footnotesize\input{LALInspiralTestCV}}

</lalLaTeX> */

#include <stdio.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel=0;

void printf_timeseries (int n, float *signal, double delta, double t0) ;
void printf_timeseries (int n, float *signal, double delta, double t0) 
{
  int i=0;
  FILE *outfile1;

  outfile1=fopen("wave1.dat","w");

  do 
     fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i));

  while (n-++i); 
  fprintf(outfile1,"&\n");
}


int main () {
   REAL4Vector signal;
   InspiralTemplate params;
   double dt;
   static LALStatus status;
   INT4 n;

   params.ieta=1; 
   params.mass1=1.40; 
   params.mass2=1.40; 
   params.startTime=0.0; 
   params.startPhase=0.0; 
   params.fLower=40.0; 
   params.fCutoff=1000.0;
   params.tSampling=4000.0;
   params.signalAmplitude=1.0;
   params.nStartPad=1000;
   params.nEndPad=1000;
   params.method=two;
   params.order=twoPN;
   params.domain=TimeDomain;
   params.approximant=pade;
   params.massChoice=m1Andm2;

   LALInspiralWaveLength (&status, &n, params);
   LALInspiralParameterCalc (&status, &params);
   fprintf(stderr, "signal length=%d\n", n);
   signal.length = n;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*n);
   LALInspiralWave (&status, &signal, &params);
   dt = 1./params.tSampling;
   printf_timeseries(signal.length, signal.data, dt, params.startTime);
   LALFree(signal.data);
   signal.data = NULL;
   LALCheckMemoryLeaks();

   return 0;
}
