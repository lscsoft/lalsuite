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
INT4 lalDebugLevel=1;


typedef enum{
  WAVEFORM  = 0, Waveform = 0, waveform = 0,
    OVERLAP = 1, Overlap = 1, overlap = 1 ,
    LENGTH  = 2, Length =2, 
    CONVERSION = 3, Conversion = 3, conversion = 3,
    LASTFREQ = 4
} Tools; 


void printf_timeseries (int n, float *signal, double delta, double t0) ;
void SetDefault(InspiralTemplate *params);
void OtherDefault(UINT4 *tool);
void ParseParameters(int , char **, InspiralTemplate *, UINT4 *);
void PrintParams(InspiralTemplate params);
void PrintParamsMass(InspiralTemplate params);
void GenerateWaveform(LALStatus *status, InspiralTemplate params);
void LastFreq(LALStatus *status, InspiralTemplate params);

/*/NRCSID(INSPIRALTOOLSC, "InspiralTools, main");*/

int main (int argc, char **argv ) {
  static InspiralTemplate params;
  UINT4   tool;
  static LALStatus status;
  
  
  /*/INITSTATUS(&status, "Main InspiralTools", INSPIRALTOOLSC );8/*/
  /*/ATTATCHSTATUSPTR(&status); */

  SetDefault(&params);
  OtherDefault(&tool);
  ParseParameters(argc, argv, &params, &tool);
   
   
   /*TODO 
     3 - add the help documenation
     4 - genenral documenation
   */
   

   switch(tool)
     {
     case Waveform:
       GenerateWaveform(&status, params);       
/*       CHECKSTATUSPTR(&status);*/
       break;
       /*TODO */ 
     case Overlap:
       break;
     case LASTFREQ:
       LastFreq(&status, params);
      break;
     case Length:
       break;
     case Conversion:
       LALInspiralParameterCalc(&status, &params);
       PrintParamsMass(params);
     }

/*   DETATCHSTATUSPTR(status);
   return (status);*/
  return 0;
}



/*void GenerateInspiralWaveform(){}*/



void 
SetDefault(InspiralTemplate *params)
{
  params->approximant=EOB;
  params->OmegaS = 0.;
  params->Zeta2 = 0.; /*use by EOB @ 3PN*/
  params->Theta = 0.;
  params->ieta=1; 
  params->mass1=10; 
  params->mass2=1.4; 
  params->startTime=0.0; 
  params->startPhase=0.;
  params->fLower=40.0; 
  params->fCutoff=2000.0;
  params->tSampling=4096.0;
  params->order=4;
  params->signalAmplitude=1.0;
  params->nStartPad=1000;
  params->nEndPad=1200;
  params->massChoice=m1Andm2;
  params->distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
  params->alpha = 0;
  params->alpha1 = 0;
  params->alpha2 = 0;
  params->beta = 0;
  params->psi0 = 80000.;
  params->psi3 = -100;
  params->fendBCV = 331.4;
}

void
OtherDefault(UINT4 *tool)
{
  *tool = WAVEFORM;
}



void 
ParseParameters(int argc, 
		char **argv,
		InspiralTemplate *params,
		UINT4 *tool)
{
  int i=1;

   while(i <argc)
   {
     if (strcmp(argv[i],"-tool")==0)
       {
	 if (strcmp(argv[++i],"WAVEFORM")==0
	     || 
	     strcmp(argv[i],"Waveform")==0
	     ||
	     strcmp(argv[i],"waveform")==0
	     )
	   *tool = WAVEFORM;
	 else if (strcmp(argv[i],"CONVERSION")==0
		  || 
		  strcmp(argv[i],"Conversion")==0
		  ||
		  strcmp(argv[i],"conversion")==0
		  )
	   *tool = CONVERSION;
	 else if (strcmp(argv[i],"LASTFREQ")==0) *tool = LASTFREQ;
       }
     else if (strcmp(argv[i],"-fl")==0)
       params->fLower = atof(argv[++i]);
     else if (strcmp(argv[i],"-nstartpad")==0)
       params->nStartPad = atoi(argv[++i]); 
     else if (strcmp(argv[i],"-alpha")==0)
       params->alpha = atof(argv[++i]); 
     else if (strcmp(argv[i],"-alpha1")==0)
       params->alpha1 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-alpha2")==0)
       params->alpha2 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-beta")==0)
       params->beta = atof(argv[++i]); 	  
     else if (strcmp(argv[i],"-order")==0)
       params->order = atoi(argv[++i]); 	   
     else if (strcmp(argv[i],"-psi0")==0)
       params->psi0 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-psi3")==0)
       params->psi3 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-sampling")==0)
       params->tSampling = atof(argv[++i]); 
     else if (strcmp(argv[i],"-m1")==0)
       params->mass1 = atof(argv[++i]); 
     else if (strcmp(argv[i],"-m2")==0)
       params->mass2 = atof(argv[++i]); 


     else if (strcmp(argv[i],"-fbcv")==0)
       params->fendBCV = atof(argv[++i]); 


	   else if (strcmp(argv[i], "-approximant")==0)
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   params->approximant = TaylorT1;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   params->approximant = TaylorT2;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   params->approximant = TaylorT3;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   params->approximant = TaylorF1;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   params->approximant = TaylorF2;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   params->approximant = PadeT1;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   params->approximant = PadeF1;
		   else if (strcmp(argv[i],"EOB")==0)
			   params->approximant = EOB;
		   else if (strcmp(argv[i],"BCV")==0)
			   params->approximant = BCV;
		   else if (strcmp(argv[i],"BCVSpin")==0)
			   params->approximant = BCVSpin;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   params->approximant = SpinTaylorT3;

	   }	
	   else 
	   {
		   fprintf(stderr,"\nUSAGE: %s [options]\n", argv[0]);
		   fprintf(stderr,"The options are (with default values in brackets)\n");
		   fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n",params->alpha);
		   fprintf(stderr,"-approximant : Post-Newtonian model such as TaylorT1(2,3), TaylorF1(2), PadeT1, PadeF1, EOB, BCV, BCVSpin and SpinTaylorT3  (PadeT1)\n");
		   fprintf(stderr,"       -psi0 : First parameter for BCV template (%7.2f)\n",params->psi0);
		   fprintf(stderr,"       -psi3 : Second parameter for BCV template \n");
		   exit(0);

	   }
	   i++;       
   }
}

void
PrintParamsMass(InspiralTemplate params)
{
  printf("m1= %15.12e\n", params.mass1);
  printf("m2= %15.12e\n",params.mass2);
  printf("psi0= %15.12e\n",params.psi0);
  printf("psi3= %15.12e\n",params.psi3);
  printf("M= %15.12e\n",params.totalMass);
  printf("eta= %15.12e\n",params.eta);
}  


void
PrintParams(InspiralTemplate params)
{
  printf("\nieta     = %15.12d\n", params.ieta);
  printf("level      = %15.12d\n",params.level);
  printf("nStartPad  = %15.12d\n",params.nStartPad);
  printf("nEndPad    = %15.12d\n",params.nEndPad);
  printf("minMatch   = %15.12f\n",params.minMatch);
  printf("mass1      = %15.12f\n",params.mass1); 
  printf("mass2      = %15.12f\n",params.mass2);
  printf("totalMass  = %15.12f\n", params.totalMass); 
  printf("chirpmass  = %15.12f\n", params.chirpMass); 
  printf("psi0       = %15.12f\n",params.psi0);
  printf("psi3       = %15.12f\n ",params.psi3);
  printf("fendBCV    = %15.12f\n",params.fendBCV);
  printf("alpha      = %15.12f\n",params.alpha);
  printf("alpha1     = %15.12f\n",params.alpha1);
  printf("alpha2     = %15.12f\n",params.alpha2);
  printf("beta       = %15.12f\n",params.beta);
  printf("tc         = %15.12f\n",params.tC); 
  printf("eta        = %15.12f\n",params.eta);
  printf("fLower     = %15.12f\n",params.fLower);
  printf("fcutoff    = %15.12f\n",params.fCutoff);
  printf("tsampling  = %15.12f\n",params.tSampling);
  printf("phase0     = %15.12f\n",params.startPhase);
  printf("t0         = %15.12f\n ",params.startTime);
  printf("fFinal     = %15.12f\n",params.fFinal);
  printf("zeta2      = %15.12f\n",params.Zeta2);
  printf("omegaS      = %15.12f\n",params.OmegaS);
  printf("order   %d\n",params.order);
  printf("approximant %15.12d\n",   params.approximant);

}




NRCSID(GENERATEWAVEFORMC, "GenerateWaveform in InspoiralTools");
void GenerateWaveform(LALStatus *status, InspiralTemplate params)
{
   static REAL4Vector *signal1, *signal2;
   static REAL8 dt;
   int n, i;
   dt = 1./params.tSampling;
   
   INITSTATUS(status,"InspiralTools", GENERATEWAVEFORMC);
   ATTATCHSTATUSPTR(status);

  LALInspiralWaveLength(status->statusPtr, &n, params);
  if (params.approximant ==BCV || params.approximant==BCVSpin)
    {
      params.massChoice = psi0Andpsi3;
      n = 32768;             
    }
  
  LALInspiralParameterCalc(status->statusPtr, &params);
  fprintf(stderr, "Testing Inspiral Signal Generation Codes for %d:\n", params.approximant);
  fprintf(stderr, "Signal length=%d, t0=%e, t2=%e, m1=%e, m2=%e, fLower=%e, fUpper=%e\n", 
		  n, params.t0, params.t2, params.mass1, params.mass2, params.fLower, params.fCutoff);


  LALCreateVector(status->statusPtr, &signal1, n);
  LALCreateVector(status->statusPtr, &signal2, n);
  
  
  /*   params.approximant = BCV;*/
  
   if (params.approximant==TaylorF1 || params.approximant==TaylorF2 
		   || params.approximant==BCV || params.approximant==BCVSpin) 
     {
	RealFFTPlan *revp = NULL;
	COMPLEX8Vector *Signal1 = NULL;
	
	LALInspiralWave(status->statusPtr, signal1, &params);
	
	LALCreateReverseRealFFTPlan(status->statusPtr, &revp, n, 0);
	
	LALCCreateVector(status->statusPtr, &Signal1, n/2+1);
	for (i=1; i<n/2; i++) 
	  {
	    Signal1->data[i].re = signal1->data[i];
		Signal1->data[i].im = signal1->data[n-i];
	  }
	Signal1->data[0].re = 0.;
	Signal1->data[0].im = 0.;
	Signal1->data[n/2].re = 0.;
	Signal1->data[n/2].im = 0.;




	
	LALReverseRealFFT(status->statusPtr, signal2, Signal1, revp);
	LALCDestroyVector (status->statusPtr, &Signal1);
	printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
	
	LALREAL4VectorFFT(status->statusPtr, signal2, signal1, revp);
	LALDestroyRealFFTPlan (status->statusPtr, &revp);
	printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
     }
   else
   {
     LALInspiralWave(status->statusPtr, signal2, &params);
     printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
     /*
       REPORTSTATUS(&status);
     */
   }
   if (status->statusCode) 
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

   LALDestroyVector(status->statusPtr, &signal2);
   LALDestroyVector(status->statusPtr, &signal1);
   DETATCHSTATUSPTR(status);
   LALCheckMemoryLeaks();
}






void LastFreq(LALStatus *status, InspiralTemplate params)
{

   static REAL4Vector *signal1;
   static REAL8 dt;
   int n, i;
   dt = 1./params.tSampling;
   
   INITSTATUS(status,"InspiralTools", GENERATEWAVEFORMC);
   ATTATCHSTATUSPTR(status);

  LALInspiralWaveLength(status->statusPtr, &n, params);
  if (params.approximant ==BCV || params.approximant==BCVSpin)
    {
      params.massChoice = psi0Andpsi3;
      n = 32768;             
    }
  
  LALInspiralParameterCalc(status->statusPtr, &params);

  LALCreateVector(status->statusPtr, &signal1, n);
  
   if (params.approximant==TaylorF1 || params.approximant==TaylorF2 
		   || params.approximant==BCV || params.approximant==BCVSpin) 
     {
	
	LALInspiralWave(status->statusPtr, signal1, &params);
	
     }
   else
   {
     LALInspiralWave(status->statusPtr, signal1, &params);
   }
   fprintf(stderr, "M= %lf F= %lf\n", params.totalMass, params.fFinal);

   LALDestroyVector(status->statusPtr, &signal1);
   DETATCHSTATUSPTR(status);
   LALCheckMemoryLeaks();
}






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


