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
   static REAL4Vector *signal1, *signal2;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt;
   UINT4 n, i;

   params.approximant=TaylorF2;
   params.OmegaS = 0.;
   params.Theta = 0.;
   params.ieta=1; 
   params.mass1=10.; 
   params.mass2=10.; 
   params.startTime=0.0; 
   params.startPhase=0.;
   params.fLower=40.0; 
   params.fCutoff=1000.0;
   params.tSampling=4096.0;
   params.order=4;
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=0;
   params.massChoice=m1Andm2;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   dt = 1./params.tSampling;

   LALInspiralWaveLength(&status, &n, params);
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "Testing Inspiral Signal Generation Codes:\n");
   fprintf(stderr, "Signal length=%d, t0=%e, t2=%e, m1=%e, m2=%e, fLower=%e, fUpper=%e\n", n, params.t0, params.t2, params.mass1, params.mass2, params.fLower, params.fCutoff);
   LALCreateVector(&status, &signal1, n);
   LALCreateVector(&status, &signal2, n);

   /*
   params.psi0 = 72639.;
   params.psi3 = -768.78;
   params.alpha = 0.766;
   params.fendBCV = 331.4;
   params.approximant = BCV;
   */
			   
   if (params.approximant==TaylorF1 || params.approximant==TaylorF2 || params.approximant==BCV) 
   {
	RealFFTPlan *revp = NULL;
	COMPLEX8Vector *Signal1 = NULL;
	LALInspiralWave(&status, signal1, &params);
	/*
	   REPORTSTATUS(&status);
	 */
	LALCreateReverseRealFFTPlan(&status, &revp, n, 0);

	LALCCreateVector(&status, &Signal1, n/2+1);
	for (i=1; i<n/2; i++) 
	{
		Signal1->data[i].re = signal1->data[i];
		Signal1->data[i].im = signal1->data[n-i];
	}
	Signal1->data[0].re = 0.;
	Signal1->data[0].re = 0.;
	Signal1->data[n/2].re = 0.;
	Signal1->data[n/2].re = 0.;

	LALReverseRealFFT(&status, signal2, Signal1, revp);
	LALCDestroyVector (&status, &Signal1);
	printf_timeseries(signal2->length, signal2->data, dt, params.startTime);

	LALREAL4VectorFFT(&status, signal2, signal1, revp);
	LALDestroyRealFFTPlan (&status, &revp);
	printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
   }
   else
   {
	LALInspiralWave(&status, signal2, &params);
	/*
	   REPORTSTATUS(&status);
	 */
   }
   fprintf(stderr, "approximant=%d order=%d,", params.approximant, params.order);

   if (status.statusCode) 
	   fprintf(stderr, " not available\n");
   else 
	   fprintf(stderr, " successful\n");
   /* 
    * static RealFFTPlan *fwdp;
    * LALCreateForwardRealFFTPlan(&status, &fwdp, n, 0);
    */
   /* 
    * LALInspiralWaveTemplates (&status, signal1, signal2, &params); 
    * LALDestroyRealFFTPlan (&status, &fwdp); 
    *
    */
   /* 
    * printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
    */
   LALDestroyVector(&status, &signal2);
   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();
   return 0;
}
