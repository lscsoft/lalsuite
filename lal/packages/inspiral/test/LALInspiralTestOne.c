/*  <lalVerbatim file="LALInspiralTestOneCV">
Author: Sathyaprakash, B. S.
$Id$
$Name$
$Author$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralTestOne.c}}
Test routine for wave generation codes. Exactly as in
{\tt LALInspiralTest} except that this code tests only 
one function as chosen by the user in {\tt Approximant}
and {\tt Order}.

\vfill{\footnotesize\input{LALInspiralTestOneCV}}

</lalLaTeX> */

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
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


int main (int argc , char **argv) {
   static REAL4Vector *signal1, *signal2;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt;
   UINT4 n, i;

   if (argc==1)
     {
       params.psi0 = 100000;
       params.psi3 = -1000;
       params.mass1= 10;
       params.mass2= 10;          
     }
   else if (argc==3)
     {
       params.psi0 = atof(argv[1]);
       params.psi3 = atof(argv[2]);
       params.mass1= atof(argv[1]);
       params.mass2= atof(argv[2]);
     }
   
   params.approximant=EOB;
   params.OmegaS = 0.;
   params.Zeta2  = 0.; /*use by EOB @ 3PN*/
   params.Theta  = 0.;
   params.ieta   = 1; 
   params.mass1  = 3; 
   params.mass2  = 3; 
   params.startTime=0.0; 
   params.startPhase=0.;
   params.fLower  = 20.0; 
   params.fCutoff = 1000.0;
   params.tSampling = 2048.0;
   params.order = 4;
   params.signalAmplitude=1.0;
   params.nStartPad = 0;
   params.nEndPad = 1200;
   params.massChoice=m1Andm2;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   dt = 1./params.tSampling;

   LALInspiralWaveLength(&status, &n, params);
  
   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "Testing Inspiral Signal Generation Codes:\n");
   fprintf(stderr, "Signal length=%d, t0=%e, t2=%e, m1=%e, m2=%e, fLower=%e, fUpper=%e\n", n, params.t0, params.t2, params.mass1, params.mass2, params.fLower, params.fCutoff);
   LALCreateVector(&status, &signal1, n);
   LALCreateVector(&status, &signal2, n);

   
   
   params.alpha = 0.;
   params.beta  = .5;
   params.alpha1 = 0.3;
   params.alpha2= 0.3;
   
   params.fFinal = 1000.;
   params.approximant = BCV;
			   
   if (params.approximant==TaylorF1 || params.approximant==TaylorF2 || params.approximant==BCV  || params.approximant==BCVSpin) 
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
	Signal1->data[0].im = 0.;
	Signal1->data[n/2].re = 0.;
	Signal1->data[n/2].im = 0.;

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
	printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
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
