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
   static REAL4Vector *signal1;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt, mass1Sq, mass2Sq, spin1Frac, spin2Frac, spin1Theta, spin1Phi, spin2Theta, spin2Phi;
   static REAL8 dTheta, dPhi;
   UINT4 n, i, count=0;

   dTheta = LAL_PI_2/2.0;
   dPhi = LAL_TWOPI/2.0;
   params.ieta=1; 
   params.mass1=10.0; 
   params.mass2=1.40; 
   params.startTime=0.0; 
   params.startPhase=0.0;
   params.fLower=40.0; 
   params.fCutoff=2000.00;
   params.tSampling=4096.0;
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=1000;
   params.order=4;
   params.approximant=TaylorT3;
   params.massChoice=m1Andm2;
   params.OmegaS = 0.;
   params.Theta = 0.;
   params.sourceTheta = LAL_PI/6.L;
   params.sourcePhi = LAL_PI/6.L;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   dt = 1./params.tSampling;
   spin1Frac = 0.9;
   spin2Frac = 0.9;

   mass1Sq = pow(params.mass1*LAL_MTSUN_SI,2.L);
   mass2Sq = pow(params.mass2*LAL_MTSUN_SI,2.L);

   LALInspiralWaveLength(&status, &n, params);
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "#signal length=%d\n", n);
   LALCreateVector(&status, &signal1, n);
   params.orbitTheta0 = LAL_PI_2/3.;
   params.orbitPhi0 = LAL_PI/6.;

   /*
   for (params.orbitTheta0=0; params.orbitTheta0<LAL_PI_2; params.orbitTheta0+=dTheta)
   for (params.orbitPhi0=0; params.orbitPhi0<LAL_TWOPI; params.orbitPhi0+=dPhi)
   */
   for (spin1Theta=LAL_PI/6.; spin1Theta<LAL_PI_2; spin1Theta+=dTheta)
   for (spin1Phi=LAL_PI/5.; spin1Phi<LAL_TWOPI; spin1Phi+=dPhi)
   for (spin2Theta=LAL_PI/4.; spin2Theta<LAL_PI_2; spin2Theta+=dTheta)
   for (spin2Phi=LAL_PI/3.; spin2Phi<LAL_TWOPI; spin2Phi+=dPhi)
   {
	   fprintf(stderr, "%d %e %e %e %e %e %e\n",
			   count++, params.orbitTheta0, params.orbitPhi0, spin1Theta, spin1Phi, spin2Theta, spin2Phi);
	   params.spin1[0] =  mass1Sq * spin1Frac * sin(spin1Theta) * cos(spin1Phi);
	   params.spin1[1] =  mass1Sq * spin1Frac * sin(spin1Theta) * sin(spin1Phi);
	   params.spin1[2] =  mass1Sq * spin1Frac * cos(spin1Theta);
	   params.spin2[0] =  mass2Sq * spin2Frac * sin(spin2Theta) * cos(spin2Phi);
	   params.spin2[1] =  mass2Sq * spin2Frac * sin(spin2Theta) * sin(spin2Phi);
	   params.spin2[2] =  mass2Sq * spin2Frac * cos(spin2Theta);
   
	   LALInspiralSpinModulatedWave(&status, signal1, &params);
   
	   /*
	   printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
	    */
   }

   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();

   return 0;
}
