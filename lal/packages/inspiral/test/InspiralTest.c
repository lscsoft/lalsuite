/* 
   Testing the interface routine inspiralwave.c 
   to generate t- or f-domain wave using either
   integrals or ODEs.
   April 5 , 00.
*/
#include <stdio.h>
#include <lal/LALInspiral.h>

INT4 lalDebugLevel=1;

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
   params.mass1=1.4; 
   params.mass2=10.; 
   params.startTime=0.0; 
   params.startPhase=0.0; 
   params.fLower=40.0; 
   params.fCutoff=1000000.0;
   params.tSampling=4000.0;
   params.signalAmplitude=1.0;
   params.nStartPad=0;
   params.nEndPad=0;
   params.method=one;
   params.order=twoPN;
   params.domain=TimeDomain;
   params.approximant=pade;
   params.massChoice=m1Andm2;

   LALInspiralWaveLength (&status, &n, params);
   fprintf(stderr, "signal length=%d\n", n);
   signal.length = n;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*n);
   LALInspiralWave (&status, &signal, &params);
   dt = 1./params.tSampling;
   printf_timeseries(signal.length, signal.data, dt, params.startTime);

   return 0;
}
