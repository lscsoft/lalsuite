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
   params.approximant=TaylorT1; (approximant in the convention of DIS 2001; 
                            TaylorT1, PadeT1=ODE solver, 
                            TaylorT2=implicit phasing formula solved in quadrature, 
                            TaylorT3=explicit time-domain phasing)
                            EOB=effective-one-body approach
   params.order=twoPN;     (also newtonian, onePN, oneAndHalfPN, twoPN, 
                            twoAndHalfPN, threePN, threeAndHalfPN)
   params.massChoice=m1Andm2; (also t0t2, t0t3, t0t4, totalMassAndEta,totalMassAndMu) 
\end{verbatim}

\vfill{\footnotesize\input{LALInspiralTestCV}}

</lalLaTeX> */

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
INT4 lalDebugLevel=1;

void printf_timeseries (int n, float *signal, double delta, double t0) ;
void printf_timeseries (int n, float *signal, double delta, double t0) 
{
  int i=0;
  FILE *outfile1;

  outfile1=fopen("wave1.dat","a");

  do 
     fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 

  fprintf(outfile1,"&\n");
  fclose(outfile1);
}


int main () {
   static REAL4Vector *signal1, *signal2;
   static LALStatus status;
   InspiralTemplate params;
   REAL8 dt;
   UINT4 n;

   params.ieta=1; 
   params.mass1=10.; 
   params.mass2=1.4; 
   params.startTime=0.0; 
   params.startPhase=0.83038459; 
   params.fLower=40.0; 
   params.fCutoff=1000.00;
   params.tSampling=2048.0;
   params.signalAmplitude=1.0;
   params.nStartPad=8000;
   params.nEndPad=1000;
   params.order=6;
   params.approximant=TaylorT1;
   params.massChoice=m1Andm2;
   params.OmegaS=0.;
   params.Theta=0.;
   dt = 1./params.tSampling;

   LALInspiralWaveLength(&status, &n, params);
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "signal length=%d\n", n);
   LALCreateVector(&status, &signal1, n);
   LALCreateVector(&status, &signal2, n);
   LALInspiralWave(&status, signal1, &params);
   LALInspiralWaveTemplates (&status, signal1, signal2, &params);

   if (params.approximant==TaylorF2)
   {
      static RealFFTPlan *fwdp,*revp;
/*
      LALCreateForwardRealFFTPlan(&status, &fwdp, n, 0);
*/
      LALCreateReverseRealFFTPlan(&status, &revp, n, 0);
      LALREAL4VectorFFT(&status, signal2, signal1, revp);
/*
      LALDestroyRealFFTPlan (&status, &fwdp);
*/
      LALDestroyRealFFTPlan (&status, &revp);
   }
   printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
   printf_timeseries(signal2->length, signal2->data, dt, params.startTime);

   LALDestroyVector(&status, &signal1);
   LALDestroyVector(&status, &signal2);
   LALCheckMemoryLeaks();

   return 0;
}
