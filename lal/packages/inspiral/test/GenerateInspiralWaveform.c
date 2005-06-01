/*  <lalVerbatim file="LALGenerateInspiralWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Test program \texttt{LALGenerateInspiralWaveform.c}}
Test routine for wave generation codes. 

To get some help just type the name of the executable and the option --h

Basically, you can provide all the arguments from the InspiralTemplate structure such as
--approximant, --order ....

\vfill{\footnotesize\input{LALGenerateInspiralWaveformCV}}

</lalLaTeX> */


#define LALGENERATEINSPIRALWAVEFORMC_ENORM 0
#define LALGENERATEINSPIRALWAVEFORMC_ESUB  1
#define LALGENERATEINSPIRALWAVEFORMC_EARG  2
#define LALGENERATEINSPIRALWAVEFORMC_EVAL  3
#define LALGENERATEINSPIRALWAVEFORMC_EFILE 4
#define LALGENERATEINSPIRALWAVEFORMC_EMEM  5

#define LALGENERATEINSPIRALWAVEFORMC_MSGENORM "Normal exit"
#define LALGENERATEINSPIRALWAVEFORMC_MSGESUB  "Subroutine failed"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEARG  "Error parsing arguments"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEVAL  "Input argument out of valid range"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEFILE "Could not open file"
#define LALGENERATEINSPIRALWAVEFORMC_MSGEMEM  "Out of memory"

#include <stdio.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>

NRCSID( LALGENERATEINSPIRALWAVEFORMC, "$Id$" );

INT4 lalDebugLevel=7;
#define ERROR( code, msg, statement )                                \
do                                                                   \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), program, __FILE__,       \
		 __LINE__, LALGENERATEINSPIRALWAVEFORMC, statement ? statement :  \
                 "", (msg) );                                        \
}                                                                    \
while (0)

#define WARNING( statement )                                         \
do                                                                   \
if ( lalDebugLevel & LALWARNING )                                    \
{                                                                    \
  LALPrintError( "Warning[0]: program %s, file %s, line %d, %s\n"    \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 LALGENERATEINSPIRALWAVEFORMC, (statement) );                         \
}                                                                    \
while (0)

#define INFO( statement )                                            \
do                                                                   \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", program, __FILE__, __LINE__,        \
		 LALGENERATEINSPIRALWAVEFORMC, (statement) );                         \
}                                                                    \
while (0)

#define SUB( func, statusptr )                                       \
do                                                                   \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( LALGENERATEINSPIRALWAVEFORMC_ESUB, LALGENERATEINSPIRALWAVEFORMC_MSGESUB,                      \
         "Function call \"" #func "\" failed:" );                    \
  exit( LALGENERATEINSPIRALWAVEFORMC_ESUB );                                          \
}                                                                    \
while (0)

typedef struct{
  INT4 PrintParameters;
} OtherParamIn;

char *program;

void printf_timeseries (UINT4 n, REAL4 *signal, REAL8 delta, REAL8 t0) ;
void ParseParameters(UINT4 argc, CHAR **argv, OtherParamIn *otherIn);
void LALGenerateInspiralWaveformHelp(void);



/* --- Main part --- */
int main (int argc , char **argv) {
  REAL4Vector *signal1 = NULL; /*storing waveforms */
  REAL4Vector *signal2 = NULL; /*storing waveforms */
  static LALStatus status;
  InspiralTemplate params; /* the parameters */
  REAL8 dt;                /* some sampling */
  UINT4 n,i;
  InspiralInit paramsInit;

  RealFFTPlan *revp =NULL;

  COMPLEX8Vector *Signal1 =NULL;

  OtherParamIn otherIn; /* some extra parameters to parse*/

  program = *argv;

  /* ---  we start real computation here --- */
  otherIn.PrintParameters = 0; /* by default we don't print the parameters */
  ParseParameters(argc, argv, &otherIn);/*let's parse user parameters     */


  SUB( LALInspiralITStructureSetDefault(&status, &params),
       &status); /*  */ 
  
  SUB( LALInspiralITStructureParseParameters(&status, argc, argv, &params), 
       &status);/*parse user inspiral template parameters */
 
  SUB(  LALInspiralInit(&status, &params, &paramsInit), &status);
  /*  params.signalAmplitude =1;*/
  
  
  if (otherIn.PrintParameters){
    SUB( LALInspiralITStructurePrint(&status, params),  &status); 
  }
     
  /* force those parameters */

  dt 	= 1./params.tSampling;
  n 	= paramsInit.nbins;   

  if (n<=10) {  
      LALWarning(&status, "#nothing to compute; length is too short. You might reduce the fLower or masses values.");
    return 0; 
  }
  
  if (otherIn.PrintParameters)
    {
      fprintf(stderr, "#Testing Inspiral Signal Generation Codes:\n");
      fprintf(stderr, "#Signal length=%d, t0=%e, t2=%e, \n", n, params.t0, params.t2);  
      fprintf(stderr,"#size in bins %d\n",n);
      fprintf(stderr,"#size in seconds %f\n",params.tC);
    }
  
  SUB( LALSCreateVector(&status, &(signal1), n), &status);
  SUB( LALSCreateVector(&status, &(signal2), n), &status);
     
  params.ieta = 1;  /*should be zero or 1 ?? */
  /* */
  
  switch (params.approximant){		   
  case TaylorF1:
  case TaylorF2:
  case BCV:
  case BCVSpin:
    SUB( LALInspiralWave(&status, signal1, &params), &status);
    SUB( LALCreateReverseRealFFTPlan(&status, &revp, n, 0), &status);
        
    SUB(   LALCCreateVector(&status, &Signal1, n/2+1), &status);
    for (i=1; i<n/2; i++) 
      {
	Signal1->data[i].re = signal1->data[i];
	Signal1->data[i].im = signal1->data[n-i];
      }
    Signal1->data[0].re = 0.;
    Signal1->data[0].im = 0.;
    Signal1->data[n/2].re = 0.;
    Signal1->data[n/2].im = 0.;

    SUB( LALReverseRealFFT(&status, signal2, Signal1, revp), &status);
    SUB( LALCDestroyVector (&status, &Signal1), &status);
    printf_timeseries(signal2->length, signal2->data, dt, params.startTime);
    SUB( LALDestroyRealFFTPlan (&status, &revp), &status);
    SUB( LALSDestroyVector(&status, &signal2), &status);
    SUB( LALSDestroyVector(&status, &signal1), &status);
    break;
  case SpinTaylorT3:
  case TaylorT1:
  case TaylorT2:
  case TaylorT3:
  case EOB:
  case PadeT1:
    SUB(LALInspiralWave(&status, signal2, &params), &status);
    if (status.statusCode == 0){
      printf_timeseries(signal2->length, signal2->data, dt, params.startTime);   	  
      SUB( LALSDestroyVector(&status, &signal2), &status);
      SUB( LALSDestroyVector(&status, &signal1), &status);
    }
    else 
      {
	SUB( LALSDestroyVector(&status, &signal1), &status);
	SUB( LALSDestroyVector(&status, &signal2), &status);	
      }
    printf("%f %d %f %f %f\n",
	   params.tC*params.tSampling,
	   n ,
	   params.totalMass,params.eta, params.fLower);
    break;
  case PadeF1:
  case SpinTaylor:
  default:
	    fprintf(stderr, " not available\n");
	    break;
  }
     
  
  LALCheckMemoryLeaks();
  return 0;
}





void 
ParseParameters(	UINT4 			argc, 
			CHAR 			**argv,
			OtherParamIn    	*otherIn)
{
  UINT4 		i = 1;
  
  while(i < argc)
    {
      if ( strcmp(argv[i],	"--verbose") 	== 0 ) {
	otherIn->PrintParameters = 1;
      }
      else if( strcmp(argv[i],	"--h") 	== 0 ) {
	LALGenerateInspiralWaveformHelp(); }
      else if( strcmp(argv[i],"-h") 	== 0 ) {
	LALGenerateInspiralWaveformHelp();
      }
      else if( strcmp(argv[i],"-help") 	== 0 ) {
	LALGenerateInspiralWaveformHelp();
      } 
      else if( strcmp(argv[i],"--help")	== 0 ) {
	LALGenerateInspiralWaveformHelp();
      }
      
      i++;
    }
  
}



void LALGenerateInspiralWaveformHelp(void)
{

  fprintf(stderr,"LALGenerateInspiralWaveform Help\n");
  fprintf(stderr, "-----------------------------------------------\n");
  fprintf(stderr, "--h for help\n");
  fprintf(stderr, "--verbose to print Inspiral Template parameters\n");
  fprintf(stderr, "-----------------------------------------------\n");
  LALInspiralITStructureHelp();

}



void printf_timeseries (UINT4 n, REAL4 *signal, REAL8 delta, REAL8 t0) 
{
  UINT4 i=0;
  FILE *outfile1;

  outfile1=fopen("wave1.dat","a");

  do 
     fprintf (outfile1,"%e %e\n", i*delta+t0, *(signal+i)  );
  while (n-++i); 

  fprintf(outfile1,"&\n");
  fclose(outfile1);
}



