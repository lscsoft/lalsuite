/*  <lalVerbatim file="LALInspiralTestCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTest.c}}
Test routine for codes that generate inspiral waveform from non-spinning black 
hole binaries.  Time domain signals are returned when the {\tt approximant} is
one of {\tt TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB, SpinTaylorT3}
and frequency domain signals are returned when the {\tt approximant} is
one of {\tt TaylorF1, TaylorF2, BCV.} This code checks every available approximant
at every order and reports whether or not there was any problem with the
generation codes.

To generate a waveform first set the \texttt{InspiralTemplate} structure (see
below for an example).  Next, to measure the length of the array required 
call the function\\
\texttt{
	LALInspiralWaveLength (\&status, \&n, params)
	}\\
The length will be returned in \texttt{n}. Finally, call the function \\
\texttt{
	LALInspiralWave(\&status, signal1, params);
	}\\
to generate the wave, which will be returned in \texttt{signal}.

Example values of the parameters that can be set (with options in brackets) is:
\begin{verbatim}
   params.OmegaS = 0.;     (Unknown 3PN parameter in energy; shown to be 0 by DJS)
   params.Theta = 0.;      (Unknown 3PN parameter in flux; arbitrarily set to 0)
   params.ieta=1;          (1 for comparable masses model, 0 for test mass model)
   params.mass1=1.4;       (masses of the component stars in solar masses) 
   params.mass2=1.4; 
   params.startTime=0.0;   (defined so that the instantaneous GW frequency 
                            is params.fLower at params.startTime)
   params.startPhase=0.0;  (0 to LAL_PI_2)
   params.fLower=40.0;     (in Hz)
   params.fCutoff=1000.0;  (in Hz)
   params.tSampling=4000.; (in Hz; should be larger than 2 fCutoff or 2 flso, 
                            whichever is smaller)
   params.signalAmplitude=1.0; 
   params.nStartPad=0;     (number of leading zero bins)
   params.nEndPad=0;       (number of trailing zero bins)
   params.approximant=TaylorT1; (TaylorT1, PadeT1=ODE solver; 
                            TaylorT2=implicit phasing formula solved in quadrature; 
                            TaylorT3=explicit time-domain phasing;
                            TaylorF1=stationary phase approx. using ODEs;
                            TaylorF2=usual stationary phase approx.;
                            EOB=effective-one-body approach)
   params.order=twoPN;     (also newtonian, onePN, oneAndHalfPN, twoPN, 
                            twoAndHalfPN, threePN, threeAndHalfPN)
   params.massChoice=m1Andm2; (also t0t2, t0t3, t0t4, totalMassAndEta,totalMassAndMu) 
   params.psi0 = 132250.;   (parameter required to generate BCV detection templates)
   params.psi3 = -1014.2;   (parameter required to generate BCV detection templates)
   params.alpha = 0.528;    (amplitude correction used in BCV templates)
   params.fFinal = 868.7;  (frequency at which the BCV template is terminated)

\end{verbatim}

\vfill{\footnotesize\input{LALInspiralTestCV}}

</lalLaTeX> */

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
INT4 lalDebugLevel=0;

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


int main (void) {
   static REAL4Vector *signal1;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt;
   UINT4 n;

   params.OmegaS = 0.;
   params.Theta = 0.;
   params.ieta=1; 
   params.mass1=5.; 
   params.mass2=5.; 
   params.startTime=0.0; 
   params.startPhase=0.0;
   params.fLower=85.0; 
   params.fCutoff=2047.00;
   params.tSampling=4096.0;
   params.order=6;
   params.approximant=EOB;
/* SpinTaylorT3 */
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=0;
   params.massChoice=m1Andm2;
   dt = 1./params.tSampling;


   LALInspiralWaveLength(&status, &n, params);
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "Testing Inspiral Signal Generation Codes:\n");
   fprintf(stderr, "Signal length=%d, m1=%e, m2=%e, fLower=%e, fUpper=%e\n", 
		   n, params.mass1, params.mass2, params.fLower, params.fCutoff);
   LALCreateVector(&status, &signal1, n);

   LALInspiralWave(&status, signal1, &params);
   fprintf(stderr,"fFinal = %e\n", params.fFinal);
	      printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
REPORTSTATUS(&status);

   printf("duration is %f \n", params.tC);
   printf("final frequency is %f \n", params.fFinal);

   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();
   return 0;
}
