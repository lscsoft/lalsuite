/*  <lalVerbatim file="LALInspiralSpinningBHBinariesTestCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralSpinningBHBinariesTest.c}}
Test routine for spin-modulted inspiral waveform generation code. First set the 
\texttt{InspiralTemplate} structure (example given below). Then call the function\\
\texttt{
	LALInspiralWaveLength (\&status, \&n, params)
	}\\
to measure the length {\tt n} of the array required. Then call the function \\
\texttt{
	LALInspiralWave(\&status, signal1, params);
	}\\
to generate the wave, which will be returned in \texttt{signal}.
Example values of the parameters that can be set (with options in 
brackets):

\begin{verbatim}
   params.OmegaS = 0.;     (Unknown 3PN parameter in energy shown to be zero by DJS)
   params.ieta=1;          (1 for comparable masses model, 0 for test mass model)
   params.mass1=1.4;       (masses of the component stars in solar masses) 
   params.mass2=1.4; 
   params.startTime=0.0;   (in sec. Defined so that the instantaneous GW frequency 
                            is params.fLower at params.startTime)
   params.startPhase=0.0;  (0 to LAL_PI_2)
   params.fLower=40.0;     (in Hz)
   params.fCutoff=1000.0;  (in Hz)
   params.tSampling=4000.; (in Hz; should be larger than 2*fCutoff or 2*flso, 
                            whichever is smaller)
   params.signalAmplitude=1.0; 
   params.nStartPad=0;     (number of leading zero bins)
   params.nEndPad=0;       (number of trailing zero bins)
   params.approximant=SpinTaylorT3; 
   params.order=twoPN;     (also newtonian, onePN, oneAndHalfPN, twoPN, 
                            twoAndHalfPN, threePN, threeAndHalfPN)
   params.massChoice=m1Andm2; (also t0t2, t0t3, t0t4, totalMassAndEta,totalMassAndMu) 
   params.Theta = 0.;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI; (distance in sec)
   params.sourceTheta = LAL_PI/6.L;   (Source co-latitute)
   params.sourcePhi = LAL_PI/6.L;   (Source azimuth)

   mass1Sq = pow(params.mass1*LAL_MTSUN_SI,2.L); (mass of the 1st body in sec.)
   mass2Sq = pow(params.mass2*LAL_MTSUN_SI,2.L); (mass of the 2nd body in sec.)
   spin1Frac = 0.9;        (spin of body 1 in units of its mass)
   spin2Frac = 0.3;        (spin of body 2 in units of its mass)

   params.orbitTheta0 = LAL_PI_2/3.; (initial co-latitute orientation of the orbit)
   params.orbitPhi0 = LAL_PI/6.; (initial azimuth orientation of the orbit)

   (spin of body 1)
   params.spin1[0] =  mass1Sq * spin1Frac * sin(spin1Theta) * cos(spin1Phi);  
   params.spin1[1] =  mass1Sq * spin1Frac * sin(spin1Theta) * sin(spin1Phi);
   params.spin1[2] =  mass1Sq * spin1Frac * cos(spin1Theta);

   (spin of body 2)
   params.spin2[0] =  mass2Sq * spin2Frac * sin(spin2Theta) * cos(spin2Phi);  
   params.spin2[1] =  mass2Sq * spin2Frac * sin(spin2Theta) * sin(spin2Phi);
   params.spin2[2] =  mass2Sq * spin2Frac * cos(spin2Theta);
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
   UINT4 n,  count=0;

   dTheta = LAL_PI_2/1.9;
   dPhi = LAL_TWOPI/1.9;
   params.OmegaS = 0.;
   params.Theta = 0.;
   params.ieta=1; 
   params.mass1=10.0; 
   params.mass2=1.0; 
   params.startTime=0.0; 
   params.startPhase=0.0;
   params.fLower=40.0; 
   params.fCutoff=2000.00;
   params.tSampling=4096.0;
   params.order=4;
   params.approximant=TaylorT3;
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=1000;
   params.massChoice=m1Andm2;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   params.sourceTheta = LAL_PI/6.L;
   params.sourcePhi = LAL_PI/6.L;
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
   params.orbitTheta0 = 0.;
   params.orbitPhi0 = 0.;

   /*
   for (params.orbitTheta0=0; params.orbitTheta0<LAL_PI_2; params.orbitTheta0+=dTheta)
   for (params.orbitPhi0=0; params.orbitPhi0<LAL_TWOPI; params.orbitPhi0+=dPhi)
   */
   for (spin1Theta=0.; spin1Theta<LAL_PI_2; spin1Theta+=dTheta)
   for (spin1Phi=0.; spin1Phi<LAL_TWOPI; spin1Phi+=dPhi)
   for (spin2Theta=0.; spin2Theta<LAL_PI_2; spin2Theta+=dTheta)
   for (spin2Phi=0.; spin2Phi<LAL_TWOPI; spin2Phi+=dPhi)
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
	   printf_timeseries(signal1->length, signal1->data, dt, params.startTime);
   
	   /*
	    */
   }

   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();

   return 0;
}
