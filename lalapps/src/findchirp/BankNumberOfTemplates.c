/* **************************************************
 * Author: Cokelaer, T. Mars 2004
 * Purpose : compute number of templates in a bank 
 * *************************************************/

#include <stdio.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#include <lalapps.h>

NRCSID( BANKNUMBEROFTEMPLATESC, "$Id$");
RCSID("");

/* Some Error messages */
#define BANKNUMBEROFTEMPLATES_ENORM  0
#define BANKNUMBEROFTEMPLATES_ESUB   1
#define BANKNUMBEROFTEMPLATES_EARG   2
#define BANKNUMBEROFTEMPLATES_EVAL   3
#define BANKNUMBEROFTEMPLATES_EFILE  4
#define BANKNUMBEROFTEMPLATES_EINPUT 5
#define BANKNUMBEROFTEMPLATES_EMEM   6

#define BANKNUMBEROFTEMPLATES_MSGENORM  "Normal exit"
#define BANKNUMBEROFTEMPLATES_MSGESUB   "Subroutine failed"
#define BANKNUMBEROFTEMPLATES_MSGEARG   "Error parsing arguments"
#define BANKNUMBEROFTEMPLATES_MSGEVAL   "Input argument out of valid range"
#define BANKNUMBEROFTEMPLATES_MSGEFILE  "Could not open file"
#define BANKNUMBEROFTEMPLATES_MSGEINPUT "Error reading file"
#define BANKNUMBEROFTEMPLATES_MSGEMEM   "Out of memory"

/* Some constantes to fill InspiralTemplate Structure, Bank structure and internal structure */
/* Bank structure first: */

#define BANKNUMBEROFTEMPLATES_FLOWER       		  40.
#define BANKNUMBEROFTEMPLATES_TSAMPLING    		2048.
#define BANKNUMBEROFTEMPLATES_FUPPER       		1000.
#define BANKNUMBEROFTEMPLATES_ORDER_SIGNAL     	twoPN
#define BANKNUMBEROFTEMPLATES_ORDER_TEMPLATE   	twoPN
#define BANKNUMBEROFTEMPLATES_MMCOARSE     		0.8
#define BANKNUMBEROFTEMPLATES_MMFINE       		0.9
#define BANKNUMBEROFTEMPLATES_MMIN            		5.
#define BANKNUMBEROFTEMPLATES_MMAX           		20.
#define BANKNUMBEROFTEMPLATES_ALPHASIGNAL    		0.
#define BANKNUMBEROFTEMPLATES_DALPHA    		0.01
#define BANKNUMBEROFTEMPLATES_ALPHAMIN    		0.
#define BANKNUMBEROFTEMPLATES_ALPHAMAX    		1.
#define BANKNUMBEROFTEMPLATES_ALPHABANK       		0.01
#define BANKNUMBEROFTEMPLATES_SPACE    		Psi0Psi3
#define BANKNUMBEROFTEMPLATES_IFLSO           		0.
#define BANKNUMBEROFTEMPLATES_PSI0MIN        		10.
#define BANKNUMBEROFTEMPLATES_PSI0MAX    		250000.
#define BANKNUMBEROFTEMPLATES_PSI3MIN     		-2200.
#define BANKNUMBEROFTEMPLATES_PSI3MAX       		-10.
#define BANKNUMBEROFTEMPLATES_NFCUT           		5
#define BANKNUMBEROFTEMPLATES_SIGNAL  			TaylorT1
#define BANKNUMBEROFTEMPLATES_BANK          		BCV
#define BANKNUMBEROFTEMPLATES_TEMPLATE        		BCV
#define BANKNUMBEROFTEMPLATES_TYPE            		0
#define BANKNUMBEROFTEMPLATES_SIGNALAMP      		10.
#define BANKNUMBEROFTEMPLATES_IETA            		1.
#define BANKNUMBEROFTEMPLATES_STARTTIME       		0.
#define BANKNUMBEROFTEMPLATES_STARTPHASE    		0.
#define BANKNUMBEROFTEMPLATES_NSTARTPHASE  		1000
#define BANKNUMBEROFTEMPLATES_SIGNALAMPLITUDE 		10.
#define BANKNUMBEROFTEMPLATES_NENDPAD         		0
#define BANKNUMBEROFTEMPLATES_NOISEAMP        		1.
#define BANKNUMBEROFTEMPLATES_NTRIALS         		2
#define BANKNUMBEROFTEMPLATES_USEED        		122888
#define BANKNUMBEROFTEMPLATES_LOWGM                    3
#define BANKNUMBEROFTEMPLATES_HIGHGM                   6
/* Other Parameters  1 = true ; 0 = false	*/
#define BANKNUMBEROFTEMPLATES_QUIETFLAG       	        0 				/* silent 				*/ 
#define BANKNUMBEROFTEMPLATES_AMBIGUITYFUNCTION      	0				/* Print Ambiguity function		*/
#define BANKNUMBEROFTEMPLATES_FMAXIMIZATION      	0				/* Print SNR function of fendBCV	*/
#define BANKNUMBEROFTEMPLATES_PRINTOVERLAP             0				/* Print Best Overlap 			*/
#define BANKNUMBEROFTEMPLATES_PRINTFILTER              0				/* Print corresponding Filter		*/
#define BANKNUMBEROFTEMPLATES_PRINTBANK		0				/* print the bank of template 		*/
#define BANKNUMBEROFTEMPLATES_CHECK                    0				/* Just check that SNR=1 for identical parameters */
#define BANKNUMBEROFTEMPLATES_RANDOMINJECTION		1				/* Type of injection: random  ?		*/		
#define BANKNUMBEROFTEMPLATES_REGULARINJECTION		0				/* Type of injection: regular ?		*/		


#define DeltaT      				256 

typedef enum{
  LIGOI,
    LIGOII,
    GEO,
    TAMA,
    VIRGO, 
    REALPSD
}DetectorName;



/* Choice of the overlap method. InQuadrature maximize over 
 * the phase parameter; AlphaMaximization maximize over both 
 * phase and alpha parameter (BCV templates) */
typedef enum {
  InQuadrature,
  AlphaMaximization
} OverlapMethodIn;


/* Internal parameters for the Bank Efficiency code:
 */
typedef struct	
{
  INT4 signal; /*random signal*/
  INT4 template;/*template and bank*/
  INT4 bank;
  char *inputPSD;
  DetectorName NoiseModel;
  char *filename;
  double alphaMin, alphaMax, dalpha;

} OtherParamIn;



/* Init the CoarseBank structure to dummy values */
void InitInspiralCoarseBankIn	(	InspiralCoarseBankIn 	*icoarseIn);
void InitRandomInspiralSignalIn	(	RandomInspiralSignalIn 	*randIn);
void InitOtherParamIn 		(	OtherParamIn 		*otherIn);
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseIn, 
				RandomInspiralSignalIn	*randIn, 
				OtherParamIn		*otherIn);


/* Parsing function */
void
ParseParameters(int 			*argc, 
		char 			**argv,
		InspiralCoarseBankIn 	*coarseIn,
		RandomInspiralSignalIn 	*randIn,
		OtherParamIn 		*otherIn);		     

/* function to check validity of parameters*/
void 	
CheckParams(InspiralCoarseBankIn 	coarseIn,
	    RandomInspiralSignalIn 	randIn,
	    OtherParamIn 		otherIn);

/* Help Function */
void 	Help(  );


/* ===================================== THE MAIN PROGRAM ===================================================== */
int
main (int argc, char **argv ) 
{
  /* --- Variables ---*/
  INT4    	i,     kk,    nlist;
  REAL4  	temp;
  REAL8    	df, alpha;
  REAL4Vector   signal;
  void   			*noisemodel;
  RandomInspiralSignalIn	randIn;						/* random signal waveform to inject	*/
  OtherParamIn			otherIn;					/* personal structure to parse params	*/
  InspiralTemplateList 		*list=NULL;
  static LALStatus 		status;
  InspiralCoarseBankIn      	coarseIn; 					/* strcture for the bank of templates	*/
  FILE 				*Finput;   



  lalDebugLevel = LALERROR;
  /* --- Initialisation of parameters and variables --- */

  ParametersInitialization(&coarseIn,&randIn,&otherIn);
  ParseParameters(&argc,argv, &coarseIn, &randIn, &otherIn);			/* Read Parameters from user 		*/
  CheckParams(coarseIn, randIn, otherIn);					/* Check validity of some variables. 	*/

  
  randIn.param.massChoice 	= m1Andm2;  					/* Only to compute the length of "signal"*/ 
  signal.length 		= 0.;
  randIn.param.OmegaS  		= 0.;
  randIn.param.Theta   		= 0.;
  randIn.param.approximant 	= EOB;						/* Only to compute the length of "signal"*/
  LALInspiralWaveLength (&status, &signal.length, randIn.param);
   
  randIn.param.approximant = otherIn.signal;					/* Retrieve the actual approximant of the waveform to inject*/
  
  
  /* --- Allocate memory --- */

  randIn.psd.length 	= signal.length/2 + 1;


  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 	= 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), randIn.psd.length );
  /* spectrum */    
  coarseIn.shf.deltaF 	= randIn.param.tSampling / signal.length;
  df = randIn.param.tSampling/(float) signal.length;
  switch (otherIn.NoiseModel)
    {
    case LIGOI:  noisemodel = LALLIGOIPsd ; 
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df), &status);
      break;
    case VIRGO: noisemodel = LALVIRGOPsd;  
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df), &status);
      break;
    case GEO:   noisemodel = LALGEOPsd;    
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df), &status);
      break;
    case TAMA:  noisemodel = LALTAMAPsd;  
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df), &status);
      break;
    case REALPSD:
      /* Read psd dans le fichier InputPSD.dat */

      Finput = fopen(otherIn.filename,"r");

      for (i=0; i<(int)coarseIn.shf.data->length; i++)
	{
	  /* factor DeltaT since the input data have a sampling frequency of 1./DeltaT */
	  fscanf(Finput, "%f %lf\n", &temp,&coarseIn.shf.data->data[i]);
	  for (kk=1; kk < (int)( DeltaT * df); kk++)
	      fscanf(Finput, "%f %f\n", &temp,&temp);
	}
      fclose(Finput);
      break;
        
    default:    noisemodel = LALLIGOIPsd; 
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df), &status);
      break;
    }


  /* if BCV template */
  if (coarseIn.approximant ==BCV){
    for (alpha = otherIn.alphaMin; alpha <= otherIn.alphaMax; alpha+=otherIn.dalpha){

      coarseIn.alpha = alpha;

      LAL_CALL(LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn), &status);
      
      printf( "%lf %d\n",alpha,nlist);
      LALFree(list);
      list=NULL;
    }
  }
  else {
    
    coarseIn.approximant = TaylorF2;
    LAL_CALL(LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn), &status);
      
      printf( " %d\n",nlist);
      LALFree(list);
      list=NULL;
    }

  

  fflush(stdout);
   
   /* --- destroy the plans, correlation and signal --- */
   LALDDestroyVector( &status, &(coarseIn.shf.data) );
   LALCheckMemoryLeaks();    
   return 0;
}




/* ****************************************************************************
 * Initialization of the parameters fromn random signal, bank and others 
 * structures
 * ***************************************************************************/
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseIn, 
				RandomInspiralSignalIn	*randIn, 
				OtherParamIn		*otherIn)
{
  InitInspiralCoarseBankIn(coarseIn);
  InitRandomInspiralSignalIn(randIn);
  InitOtherParamIn(otherIn);
}
/* ****************************************************************************
 * CoarseIn initialization
 ***************************************************************************/
void InitInspiralCoarseBankIn(InspiralCoarseBankIn *coarseIn)
{
  coarseIn->fLower    	= BANKNUMBEROFTEMPLATES_FLOWER;
  coarseIn->fUpper    	= BANKNUMBEROFTEMPLATES_TSAMPLING/2. - 1.;
  coarseIn->tSampling 	= BANKNUMBEROFTEMPLATES_TSAMPLING;
  coarseIn->space     	= BANKNUMBEROFTEMPLATES_SPACE;
  coarseIn->mmCoarse  	= BANKNUMBEROFTEMPLATES_MMCOARSE;
  coarseIn->mmFine    	= BANKNUMBEROFTEMPLATES_MMFINE;
  coarseIn->iflso       = BANKNUMBEROFTEMPLATES_IFLSO ;
  coarseIn->mMin      	= BANKNUMBEROFTEMPLATES_MMIN;
  coarseIn->mMax      	= BANKNUMBEROFTEMPLATES_MMAX;
  coarseIn->MMax      	= BANKNUMBEROFTEMPLATES_MMAX * 2;
  coarseIn->massRange   = MinMaxComponentMass; 
  coarseIn->etamin 	= coarseIn->mMin * 
	  		( coarseIn->MMax - coarseIn->mMin) 
			/ pow(coarseIn->MMax, 2.);
  coarseIn->psi0Min   	= BANKNUMBEROFTEMPLATES_PSI0MIN;
  coarseIn->psi0Max   	= BANKNUMBEROFTEMPLATES_PSI0MAX;
  coarseIn->psi3Min   	= BANKNUMBEROFTEMPLATES_PSI3MIN;
  coarseIn->psi3Max   	= BANKNUMBEROFTEMPLATES_PSI3MAX;
  coarseIn->alpha     	= BANKNUMBEROFTEMPLATES_ALPHABANK;
  /*unsigned variable ... */
  coarseIn->numFcutTemplates = BANKNUMBEROFTEMPLATES_NFCUT;
  coarseIn->approximant = BANKNUMBEROFTEMPLATES_TEMPLATE;  
  coarseIn->order       = BANKNUMBEROFTEMPLATES_ORDER_TEMPLATE;  
  coarseIn->LowGM     	= BANKNUMBEROFTEMPLATES_LOWGM;
  coarseIn->HighGM    	= BANKNUMBEROFTEMPLATES_HIGHGM;
}

/* ****************************************************************************
 * Random initialization 
 * ****************************************************************************/
void InitRandomInspiralSignalIn(RandomInspiralSignalIn *randIn)
{
  randIn->useed         = BANKNUMBEROFTEMPLATES_USEED;					/* seed for random MC simulation 	*/ 
  randIn->type          = BANKNUMBEROFTEMPLATES_TYPE;					/* type of simulation 			*/
  randIn->SignalAmp     = BANKNUMBEROFTEMPLATES_SIGNALAMP;				/* if x = n +h 				*/
  randIn->param.order   = BANKNUMBEROFTEMPLATES_ORDER_SIGNAL;   			/* and its order			*/
  randIn->param.alpha	= BANKNUMBEROFTEMPLATES_ALPHASIGNAL; 				/* alpha value of BCV 			*/
  randIn->param.ieta    = BANKNUMBEROFTEMPLATES_IETA; 					/*					*/
  randIn->param.mass1 	= BANKNUMBEROFTEMPLATES_MMIN;					/* To aalocate memory			*/ 
  randIn->param.mass2  	= BANKNUMBEROFTEMPLATES_MMIN; 					/* idem 				*/
  randIn->param.fLower  = BANKNUMBEROFTEMPLATES_FLOWER;				/* Lower cutoff freq. of the template	*/
  randIn->param.OmegaS	= 0.;							/* EOB parameter 			*/
  randIn->param.Theta	= 0.;							/* EOB parameter 			*/
  randIn->mMin          = BANKNUMBEROFTEMPLATES_MMIN;					/* min mass to inject			*/
  randIn->mMax          = BANKNUMBEROFTEMPLATES_MMAX;					/* max mass to inject 			*/
  randIn->MMax          = BANKNUMBEROFTEMPLATES_MMAX * 2;				/* total mass max 			*/
  randIn->etaMin        = (BANKNUMBEROFTEMPLATES_MMAX - BANKNUMBEROFTEMPLATES_MMIN) 
	  		/BANKNUMBEROFTEMPLATES_MMAX/ BANKNUMBEROFTEMPLATES_MMAX;
  randIn->psi0Min  	= BANKNUMBEROFTEMPLATES_PSI0MIN;				/* psi0 range 				*/
  randIn->psi0Max  	= BANKNUMBEROFTEMPLATES_PSI0MAX;	
  randIn->psi3Min  	= BANKNUMBEROFTEMPLATES_PSI3MIN;	
  randIn->psi3Max  	= BANKNUMBEROFTEMPLATES_PSI3MAX;	
  randIn->param.approximant 	= BANKNUMBEROFTEMPLATES_SIGNAL;			/* approximant of h, the signal		*/
  randIn->param.tSampling 	= BANKNUMBEROFTEMPLATES_TSAMPLING;			/* sampling 				*/
  randIn->param.startTime       = BANKNUMBEROFTEMPLATES_STARTTIME; 
  randIn->param.startPhase      = BANKNUMBEROFTEMPLATES_STARTPHASE; 
  randIn->param.nStartPad       = BANKNUMBEROFTEMPLATES_NSTARTPHASE;
  randIn->param.signalAmplitude = BANKNUMBEROFTEMPLATES_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKNUMBEROFTEMPLATES_NENDPAD;
  randIn->NoiseAmp              = BANKNUMBEROFTEMPLATES_NOISEAMP;
}

/* ****************************************************************************
 * Other Param initialization							
 * ***************************************************************************/
void InitOtherParamIn(OtherParamIn *otherIn)
{
  otherIn->template     = BANKNUMBEROFTEMPLATES_TEMPLATE;
  otherIn->bank         = -1;
  otherIn->signal       = BANKNUMBEROFTEMPLATES_SIGNAL;
  otherIn->inputPSD     = NULL;
  otherIn->dalpha       = BANKNUMBEROFTEMPLATES_DALPHA;
  otherIn->alphaMax     = BANKNUMBEROFTEMPLATES_ALPHAMAX;
  otherIn->alphaMin     = BANKNUMBEROFTEMPLATES_ALPHAMIN;
}



/* --- Function to parse user parameters --- */
void 
ParseParameters(	int 			*argc, 
			char 			**argv,
			InspiralCoarseBankIn 	*coarseIn,
			RandomInspiralSignalIn 	*randIn,
			OtherParamIn    	*otherIn)
{
  INT4 i        = 1;

  while(i < *argc)
    {
      if ( strcmp(argv[i],	"--fend-bcv") 	== 0 ) {
		coarseIn->LowGM      	= atof(argv[++i]);
       	      	coarseIn->HighGM     	= atof(argv[++i]);
	}
      else if ( strcmp(argv[i],	"--fl") 	== 0 ) 
  		coarseIn->fLower = randIn->param.fLower	= atof(argv[++i]);  
      else if ( strcmp(argv[i],	"--sampling") 	== 0 ) {
  		coarseIn->tSampling = randIn->param.tSampling = atof(argv[i]);  
		randIn->param.fCutoff 	= coarseIn->tSampling/2. - 1.;
	}      
      else if ( strcmp(argv[i],	"--mass-range")	== 0 ) 	{	
		coarseIn->mMin = randIn->mMin = atof(argv[++i]); 
		randIn->param.mass1 = randIn->param.mass2 = randIn->mMin;
		coarseIn->mMax = randIn->mMax = atof(argv[++i]); 
	}
      else if (strcmp(argv[i],	"--psi0-range")	==0) {
		coarseIn->psi0Min = randIn->psi0Min = atof(argv[++i]);
		coarseIn->psi0Max = randIn->psi0Max = atof(argv[++i]);
	}
      else if (strcmp(argv[i],	"--psi3-range")	==0){
		coarseIn->psi3Max = randIn->psi0Max = atof(argv[++i]);
		coarseIn->psi3Min = randIn->psi3Min = atof(argv[++i]);
	}
      else if ( strcmp(argv[i],	"--mm")		==0)	
      	  coarseIn->mmCoarse =atof(argv[++i]);	
      else if ( strcmp(argv[i],	"--number-fcut")	==0)  
  	  coarseIn->numFcutTemplates = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--space")	==0) {
	if (strcmp(argv[++i], "Tau0Tau2")== 0)	coarseIn->space = Tau0Tau2;
	else if (strcmp(argv[i], "Tau0Tau3") == 0)  coarseIn->space = Tau0Tau3;
	}
      else if ( strcmp(argv[i],	"--dalpha")	==0)	     
    	      otherIn->dalpha = atof(argv[++i]);
      else if ( strcmp(argv[i],	"--freq-moment-bank")	==0)	     
    	      coarseIn->fUpper = atof(argv[++i]);
      else if ( strcmp(argv[i],	"--alpha-range")==0) {
    	      otherIn->alphaMin = atof(argv[++i]);
	      otherIn->alphaMax = atof(argv[++i]);
      }
      else if ( strcmp(argv[i],	"--bank")==0) {
	if (strcmp(argv[++i], "BCV")== 0)	coarseIn->approximant = BCV;
	else if (strcmp(argv[i], "SPA") == 0) {
	  coarseIn->approximant = TaylorF2;
	}
	else {
	      fprintf(stderr,"#Error with --bank option. Argument must be either BCV or SPA\n");
	      Help();
	      exit(0);
	}
      }
      else if ( strcmp(argv[i],	"--seed")	==0)
	      randIn->useed = atoi(argv[++i])+1;	
      else if ( strcmp(argv[i],	"--noise-model")==0){
	  i++;

	  if (strcmp(argv[i], "LIGOI")		== 0)	otherIn->NoiseModel = LIGOI;
	  else if (strcmp(argv[i], "LIGOII") 	== 0) 	otherIn->NoiseModel = LIGOII;
	  else if (strcmp(argv[i], "VIRGO") 	== 0)  	otherIn->NoiseModel = VIRGO;
	  else if (strcmp(argv[i], "TAMA") 	== 0)   otherIn->NoiseModel = TAMA;
	  else if (strcmp(argv[i], "GEO") 	== 0)   otherIn->NoiseModel = GEO;
	  else {
	      otherIn->NoiseModel = REALPSD;
	      otherIn->filename= argv[i];
	    }
	}
         else 
	{
	  Help();
	  exit(0);
	}
      
      i++;       
    } 
}




/* --- Function to check validity of parsed parameters (TODO)--- */
void CheckParams(InspiralCoarseBankIn coarseIn,
		 RandomInspiralSignalIn randIn,
		 OtherParamIn otherIn)
{

  if (coarseIn.psi0Min <=0 ) {
	  fprintf(stderr,"# Parameter Parse Error: psi0 Min must be > 0\n"); 
      	  Help(coarseIn, randIn, otherIn);
	  exit(0);
  }
  if (coarseIn.psi0Max <=0 ) {
	  fprintf(stderr,"# Parameter Parse Error: psi0 Max must be > 0\n"); 
      	  Help(coarseIn, randIn, otherIn);
	  exit(0);
  }
  if (coarseIn.psi3Max >=0 ) {
	  fprintf(stderr,"# Parameter Parse Error: psi3 Max must be < 0\n"); 
      	  Help(coarseIn, randIn, otherIn);
	  exit(0);
  }
  if (coarseIn.psi3Min >=0 ) {
	  fprintf(stderr,"# Parameter Parse Error: psi3 Min must be < 0\n"); 
      	  Help(coarseIn, randIn, otherIn);
	  exit(0);
  }
  if (coarseIn.psi0Min > coarseIn.psi0Max
      || coarseIn.psi3Min > coarseIn.psi3Max){
	fprintf(stderr, "psi range should be [psiMin psiMax]; take care of the sign. \n");
      	Help(coarseIn, randIn, otherIn);
      	exit(0);
    }

  if (coarseIn.approximant!=BCV && coarseIn.space == Psi0Psi3){
	fprintf(stderr, "If --bank == SPA then --space should be use and fix to Tau0Tau3 or Tau0Tau2\n");
      	Help(coarseIn, randIn, otherIn);
      	exit(0);
  }
  if (coarseIn.approximant==BCV && coarseIn.space != Psi0Psi3){
	fprintf(stderr, "If --bank == BCV then --space should not be Tau0Tau3 or Tau0Tau2\n");
      	Help(coarseIn, randIn, otherIn);
      	exit(0);
  }	  

}



/* --- Documenation (TODO) --- */
void Help()
{
  fprintf(stderr,"\nPURPOSE:\n");
  fprintf(stderr,"Compute the number of Templates in a given bank\n");
  fprintf(stderr,"\n");
  
  fprintf(stderr,"\nUSAGE:  [options]\n");
  fprintf(stderr,"The options are (with default values in brackets)\n");
  fprintf(stderr,"All options should be followed by a number except -quiet\n");
  fprintf(stderr,"         --quiet   : if this flag is present, the output is restricted to the minimum\n");
  fprintf(stderr,"            --fl   : lower frequency cutoff             (%7.2f) Hz\n", BANKNUMBEROFTEMPLATES_FLOWER);
  fprintf(stderr,"    --psi0-range   : 2 arguments requested  (%7.2f) (%7.2f)\n",        BANKNUMBEROFTEMPLATES_PSI0MIN, BANKNUMBEROFTEMPLATES_PSI0MAX);
  fprintf(stderr,"    --psi3-range   : 2 arguments requested  (%7.2f) (%7.2f)\n",        BANKNUMBEROFTEMPLATES_PSI3MIN, BANKNUMBEROFTEMPLATES_PSI3MAX);
  fprintf(stderr,"        --dalpha   : BCV amplitude correction parameter (%7.2f)\n",    BANKNUMBEROFTEMPLATES_DALPHA);
  fprintf(stderr,"   --alpha-range   : 2 arguments requested: range of alpha parameter for BCV bank (%7.2f) (%7.2f)\n",   
		 	 BANKNUMBEROFTEMPLATES_ALPHAMIN, BANKNUMBEROFTEMPLATES_ALPHAMAX);
  fprintf(stderr,"   --number-fcut   : number of layers in Fcut dimension (%7.2d)\n",    BANKNUMBEROFTEMPLATES_NFCUT);
  fprintf(stderr,"            --mm   : minimal match for template bank    (%7.3f)\n",    BANKNUMBEROFTEMPLATES_MMCOARSE);
  fprintf(stderr,"      --fend-bcv   : 2 arguments requested: range of fcut in GM units (%7.3d) (%7.3d)\n", 
	       							  BANKNUMBEROFTEMPLATES_LOWGM, BANKNUMBEROFTEMPLATES_HIGHGM);
  fprintf(stderr,"   --noise-model   : design sensitivity (LIGOI)\n"); 
  fprintf(stderr,"         --bank    : SPA or BCV  (%7.2d)\n",    BANKNUMBEROFTEMPLATES_BANK);



}


		      
