/* 
   Testing the interface routine inspiralwave.c 
   to generate t- or f-domain wave using either
   integrals or ODEs.
   April 5 , 00.
*/
#include <stdio.h>
#include "Inspiral.h"
#include "LALStdlib.h"

INT4 LALDebugLevel=1;

void printf_timeseries (int n, double *signal, double delta, double t0) 
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
   REAL8Vector signal;
   InspiralTemplate params;
   /*
   double dt;
   */
   static LALStatus status;


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
   params.order=twoPointFivePN;
   params.domain=time;
   params.approximant=pade;
   params.massChoice=m1Andm2;

   LALInspiralWave (&status, &signal, &params);
   /*
   dt = 1./params.tSampling;
   printf_timeseries(signal.length, signal.data, dt, 0.0);
   */

   return 0;
}
