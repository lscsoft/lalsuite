/*  <lalVerbatim file="LALInspiralTestOneCV">
Author: Cokelaer T., Jones G..
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{InspiralTools}}
Test routine for wave generation codes. Exactly as in
{\tt LALInspiralTest} except that this code could test every
waveforms as chosen by the user with different kind of input 
parameters.

\vfill{\footnotesize\input{LALInspiralTestOneCV}}

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


int main (int argc, char **argv ) {
   static REAL4Vector *signal1, *signal2;
   static LALStatus status;
   static InspiralTemplate params;
   static REAL8 dt;
   UINT4 n, i, tool;

   params.approximant=EOB;
   params.OmegaS = 0.;
   params.Zeta2 = 0.; /*use by EOB @ 3PN*/
   params.Theta = 0.;
   params.ieta=1; 
   params.mass1=10; 
   params.mass2=1.4; 
   params.startTime=0.0; 
   params.startPhase=0.;
   params.fLower=40.0; 
   params.fCutoff=2000.0;
   params.tSampling=4096.0;
   params.order=4;
   params.signalAmplitude=1.0;
   params.nStartPad=1000;
   params.nEndPad=1200;
   params.massChoice=m1Andm2;
   params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   dt = 1./params.tSampling;
   params.alpha = 0;
   params.alpha1 = 0;
   params.alpha2 = 0;
   params.beta = 0;
   params.psi0 = 72639.;
   params.psi3 = -768.78;
   params.fendBCV = 331.4;
   tool = 0; /*0 -->waveform generation*/
   i=1;

/*TODO 
  1- in the switch approximant, change the number with the correct string of characters
  2 -add option for m1 and m2
  3 - add the help documenation
  4 - genenral documenation
*/
   while(i <argc)
   {
     if (strcmp(argv[i],"-tool")==0)
       tool = atoi(argv[++i]);
     else if (strcmp(argv[i],"-fl")==0)
       params.fLower = atof(argv[++i]);
     else if (strcmp(argv[i],"-alpha")==0)
       params.alpha = atof(argv[++i]); 
     else if (strcmp(argv[i],"-alpha1")==0)
       params.alpha1 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-alpha2")==0)
       params.alpha2 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-beta")==0)
       params.beta = atof(argv[++i]); 	  
     else if (strcmp(argv[i],"-order")==0)
       params.order = atoi(argv[++i]); 	   
     else if (strcmp(argv[i],"-psi0")==0)
       params.psi0 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-psi3")==0)
       params.psi3 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-fbcv")==0)
       params.fendBCV = atof(argv[++i]); 


	   else if (strcmp(argv[i], "-approximant")==0)
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   params.approximant = TaylorT1;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   params.approximant = 1;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   params.approximant = 2;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   params.approximant = 3;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   params.approximant = 4;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   params.approximant = 5;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   params.approximant = 6;
		   else if (strcmp(argv[i],"EOB")==0)
			   params.approximant = 7;
		   else if (strcmp(argv[i],"BCV")==0)
			   params.approximant = 8;
		   else if (strcmp(argv[i],"BCVSpin")==0)
			   params.approximant = 9;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   params.approximant = 10;

	   }	
	   else 
	   {
		   fprintf(stderr,"\nUSAGE: %s [options]\n", argv[0]);
		   fprintf(stderr,"The options are (with default values in brackets)\n");
		   fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n",params.alpha);
		   fprintf(stderr,"-approximant : Post-Newtonian model such as TaylorT1(2,3), TaylorF1(2), PadeT1, PadeF1, EOB, BCV, BCVSpin and SpinTaylorT3  (PadeT1)\n");
		   fprintf(stderr,"       -psi0 : First parameter for BCV template (%7.2f)\n",params.psi0);
		   fprintf(stderr,"       -psi3 : Second parameter for BCV template \n");


		   return 1;	

	   }
	   i++;       
   }







   /*should replace all the following code by some functions */

   if (tool==0)
     {

       LALInspiralWaveLength(&status, &n, params);
       if (params.approximant ==BCV || params.approximant==BCVSpin)
	 {
       params.massChoice = psi0Andpsi3;
       n = 32768;       

     }

   LALInspiralParameterCalc(&status, &params);
   fprintf(stderr, "Testing Inspiral Signal Generation Codes:\n");
   fprintf(stderr, "Signal length=%d, t0=%e, t2=%e, m1=%e, m2=%e, fLower=%e, fUpper=%e\n", n, params.t0, params.t2, params.mass1, params.mass2, params.fLower, params.fCutoff);
   LALCreateVector(&status, &signal1, n);
   LALCreateVector(&status, &signal2, n);


   /*   params.approximant = BCV;*/
			   
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


     }
   else {printf("no such tool implemneted\n");}

   LALDestroyVector(&status, &signal2);
   LALDestroyVector(&status, &signal1);
   LALCheckMemoryLeaks();
   return 0;
}



/*void GenerateInspiralWaveform(){}*/
