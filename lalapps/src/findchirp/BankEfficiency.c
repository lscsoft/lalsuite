 
/* *****************************************************************************
 * Author: Cokelaer, T. and Sathyaprakash, B. S.
 * ************************************************************************** */

/* to do
 * put the check loop inside the main loop by forcing ntrials to be=1 
 * and force to use only one template equalk to the signal
 * 3 - print in a file 3d correlation*/

/* --- some includes --- */


#include <processtable.h>
#include <stdio.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLHeaders.h>
#include <lalapps.h>
#include <lal/Date.h>
#include <getopt.h>        /* do I still use it ?  */
#include <lal/FindChirp.h> /*get rid of that. I use only one function to read t
			     he xml template bank*/

/* --- version information --- */
NRCSID( BANKEFFICIENCYC, "$Id$");
RCSID(  "$Id$");


#define CVS_ID_STRING      "$Id$"
#define CVS_NAME_STRING    "$Name$"
#define CVS_REVISION       "$Revision$"
#define CVS_SOURCE         "$Source$"
#define CVS_DATE           "$Date$"
#define PROGRAM_NAME       "BankEfficiency"

/* --- Some Error messages --- */
#define BANKEFFICIENCY_ENORM  0
#define BANKEFFICIENCY_ESUB   1
#define BANKEFFICIENCY_EARG   2
#define BANKEFFICIENCY_EVAL   3
#define BANKEFFICIENCY_EFILE  4
#define BANKEFFICIENCY_EINPUT 5
#define BANKEFFICIENCY_EMEM   6

#define BANKEFFICIENCY_MSGENORM  "Normal exit"
#define BANKEFFICIENCY_MSGESUB   "Subroutine failed"
#define BANKEFFICIENCY_MSGEARG   "Error parsing arguments"
#define BANKEFFICIENCY_MSGEVAL   "Input argument out of valid range"
#define BANKEFFICIENCY_MSGEFILE  "Could not open file"
#define BANKEFFICIENCY_MSGEINPUT "Error reading file"
#define BANKEFFICIENCY_MSGEMEM   "Out of memory"
#define BANKEFFICIENCY_MSGPARSER "Missing arguments ??  "

/* --- some constants --- 
 * useful to fill  InspiralTemplate, Bank and otherIn structure */
#define BANKEFFICIENCY_ALPHABANK       		0.01
#define BANKEFFICIENCY_ALPHASIGNAL    		0.
#define BANKEFFICIENCY_BANK          		BCV
#define BANKEFFICIENCY_FLOWER       		40.
#define BANKEFFICIENCY_FUPPER       		1000.
#define BANKEFFICIENCY_HIGHGM                   6 
#define BANKEFFICIENCY_IETA	                1
#define BANKEFFICIENCY_IFLSO           		0.
#define BANKEFFICIENCY_LOWGM                    3
#define BANKEFFICIENCY_MMCOARSE     		0.8
#define BANKEFFICIENCY_MMFINE       		0.9
#define BANKEFFICIENCY_MMIN            		5.
#define BANKEFFICIENCY_MMAX           		22.
#define BANKEFFICIENCY_NOISEMODEL           	"LIGOI"
#define BANKEFFICIENCY_NENDPAD         		0
#define BANKEFFICIENCY_NFCUT           		5
#define BANKEFFICIENCY_NOISEAMPLITUDE  		1.
#define BANKEFFICIENCY_NSTARTPHASE  		1000
#define BANKEFFICIENCY_NTRIALS         		1
#define BANKEFFICIENCY_PSI0MIN        		10.
#define BANKEFFICIENCY_PSI0MAX    		250000.
#define BANKEFFICIENCY_PSI3MIN     		-2200.
#define BANKEFFICIENCY_PSI3MAX       		-10.
#define BANKEFFICIENCY_ORDER_SIGNAL     	twoPN
#define BANKEFFICIENCY_ORDER_TEMPLATE   	twoPN
#define BANKEFFICIENCY_SIGNAL  			TaylorT1
#define BANKEFFICIENCY_SIGNALAMPLITUDE 		10.
#define BANKEFFICIENCY_SPACE    		Psi0Psi3
#define BANKEFFICIENCY_STARTTIME       		0.
#define BANKEFFICIENCY_STARTPHASE    		0.
#define BANKEFFICIENCY_TEMPLATE        		BCV
#define BANKEFFICIENCY_TYPE            		0
#define BANKEFFICIENCY_TSAMPLING    		2048.
#define BANKEFFICIENCY_USEED        		122888

/* Other Parameters  1 = true ; 0 = false	*/
#define BANKEFFICIENCY_QUIETFLAG       	        0 				/* silent 				*/ 
#define BANKEFFICIENCY_AMBIGUITYFUNCTION      	0				/* Print Ambiguity function		*/
#define BANKEFFICIENCY_FASTSIMULATION           0                               /* cheating code (dont use template far from injection; use with care */
#define BANKEFFICIENCY_FMAXIMIZATION      	0				/* Print SNR function of fendBCV	*/
#define BANKEFFICIENCY_PRINTOVERLAP             0				/* Print Best Overlap 			*/
#define BANKEFFICIENCY_PRINTOVERLAP_FILE        "BE_BestOverlap.dat"		/* Print Best Overlap in a file		*/
#define BANKEFFICIENCY_PRINTFILTER              0				/* Print corresponding Filter		*/

#define BANKEFFICIENCY_PRINTBANK		0				/* print the template bank		*/
#define BANKEFFICIENCY_PRINTBANK_FILEASCII	"BE_Bank.dat"			/* print template bank in a file	*/
#define BANKEFFICIENCY_PRINTBANKXML		0				/* print the template bank		*/
#define BANKEFFICIENCY_PRINTBANK_FILEXML	"BE_Bank.xml"			/* print template bank in a file	*/


#define BANKEFFICIENCY_PRINTRESULT		1				/* print the result (ascii)		*/
#define BANKEFFICIENCY_PRINTRESULTXML		0				/* print the template bank		*/
#define BANKEFFICIENCY_PRINTRESULT_FILEXML	"BE_Result.xml"			/* print the result (xml file)  	*/

#define BANKEFFICIENCY_PRINTPROTOTYPE		0				/* print the overlap of the templates	*/
#define BANKEFFICIENCY_PRINTBANKOVERLAP		0				/* print the overlap of the templates	*/

#define BANKEFFICIENCY_PRINTPSD                 0				/* print psd used in <x|x>      	*/
#define BANKEFFICIENCY_PRINTPSD_FILE		"BE_PSD.dat"			/* Print Psd in a file			*/

#define BANKEFFICIENCY_PRINTTEMPLATE    	0				/* print the  BCV final template	*/
#define BANKEFFICIENCY_CHECK                    0				/* Just check that SNR=1 for identical parameters */
#define BANKEFFICIENCY_RANDOMINJECTION		1				/* Type of injection: random  ?		*/		
#define BANKEFFICIENCY_REGULARINJECTION		0				/* Type of injection: regular ?		*/		

#define BANKEFFICIENCY_PRINTPROTO_FILEXML	"BE_Proto.xml"			/* print the result (xml file)  	*/
#define BANKEFFICIENCY_READXMLBANK	        "InputBank.xml"			/* print the result (xml file)  	*/


/* --- temporary flag for the sampling of real psd --- */
#define DeltaT      				256 

/* ==================== local structures ==================== */
/* An enumerate for the desigm sensitivit
 *  */


typedef enum{
   NoUserChoice,
   BBH,
   BNS,
   BHNS
} BinaryInjection;

typedef enum{
    LIGOI,
    LIGOA,
    GEO,
    TAMA,
    VIRGO,
    REALPSD
} DetectorName;

/* Structure to store the value of a 2 by 2 matrix. 
 * This matrix is used in the alpha maximization 
 * process for BCV templates.
 * */
typedef struct{
  double a11, a21, a22;
} BCVMaximizationMatrix;

/* Choice of the overlap method. InQuadrature maximize over 
 * the phase parameter; AlphaMaximization maximize over both 
 * phase and alpha parameter (BCV templates). 
 * */
typedef enum{
  InQuadrature,
  AlphaMaximization
} OverlapMethodIn;

/* Internal parameters for the Bank Efficiency code:
 * PrintOverlap	:
 * PrintFilter	:
 * overlapMethod: InQuadrature(classic overlap)  or AlphaMaximization
 * m1		: mass1 to inject
 * m2		: mass2 to inject
 * psi0		: psi0 to inject
 * psi3		: psi3 to inject
 * inputPSD	: name of an input file for the psd.
 * */
typedef struct{
  INT4 signal;				/* name of the random signal to inject 	*/	
  INT4 template;			/* name of the template in the bank 	*/
  INT4 bank;				/* type of bank to use 			*/
  UINT4 ntrials;				/* number of simulations		*/
  INT4 quietFlag;			/* a flag for verbose mode or not	*/
  INT4 ambiguityFunction;		/* do we want to save the ambiguity function ? */
  INT4 FastSimulation;                  /* target the injection in the bank --> less 
					   computation but Match might be less than the
					   optimal  is noise is injected or template and 
					   signal are not faithful (ie: BCV against physical
					   template */

  INT4 FMaximization;			
  INT4 lengthFactor;			/* multiply estimated length of filters by that factor */
  INT4 PrintOverlap;		
  INT4 PrintFilter;		
		
  INT4 PrintBankOverlap;		/* print match of each templates 	*/
  INT4 PrintBank;			/* print bank of templates 		*/
  INT4 PrintBankXml;			/* print bank of templates 		*/
  INT4 PrintResultXml;			/* print bank of templates 		*/
  INT4 PrintPrototype;
  INT4 PrintPsd;                        /* print the psd used in <x|x>          */
  BinaryInjection binaryInjection; /*injection will be set by the mass-range*/
  INT4 PrintTemplate;  

  OverlapMethodIn overlapMethod;
  INT4 check;				/* compute only one correlation where bothe 
					   template and signal have the same parameters */
  double m1,m2, psi0,psi3;
  char *inputPSD;
  char *inputXMLBank; 
  DetectorName NoiseModel;
  char *filename;

} OtherParamIn;



/* strucutre to output the results */
typedef struct{
  REAL4 psi0_trigger;
  REAL4 psi3_trigger;
  REAL4 psi0_inject;
  REAL4 psi3_inject;
  REAL4 fend_trigger;
  REAL4 fend_inject;
  REAL4 mass1_inject;
  REAL4 mass2_inject;
  REAL4 totalMass_trigger;
  REAL4 eta_trigger;
  INT4 layer;
  REAL4 rho_final;
  REAL4 alpha;
  REAL4 alpha_f;
  INT4 bin;
  REAL4 phase;
  UINT4 ntrial;
} ResultIn;


/* As to be cleaned
 * Function to store the optimal  parameters of the overlap 
 * lmax : number of the layer
 * fMax : frequency cutoff of the template
 * jmax: number of the best template. use if one want to use the FMaximization option 
 * */
void
KeepHighestValues(InspiralWaveOverlapOut in , 
		  INT4 j, INT4 l , double frequency,
		  InspiralWaveOverlapOut *out, 
		  INT4 *jmax, INT4 *lmax, double *fMax);


/* function to create the filters in the BCV overlap.
 * Input: 	- Vectors F^(-5/3) and so on
 * 		- Matrix of the moments for the amplitude
 * 		- indice at which the template starts
 * 		- psi0 and psi3 for the phase
 * Output	- The two orthogonal filters 
 * */
void 
LALCreateFilters(REAL4Vector 	*Filter1,
		 REAL4Vector 	*Filter2,
		 REAL4Vector 	VectorPowerFm5_3,
		 REAL4Vector 	VectorPowerFm2_3,
		 REAL4Vector 	VectorPowerFm7_6,
		 REAL4Vector 	VectorPowerFm1_2,
		 BCVMaximizationMatrix    matrix,
		 UINT4 		kMin,
		 REAL8 		psi0,
		 REAL8 		psi3);


/* Function to compute a orthogonal vector 
 * */
void
LALGetOrthogonal(REAL4Vector *filter);


/* Functon to compute the overlap between an injected signal 
 * and a set of templates based on BCV templates. 
 * It returns the overlap, the phase and alpha parameter. 
 * */
void
LALWaveOverlapBCV(LALStatus 		  *status,
		  REAL4Vector             *correlation,
		  InspiralWaveOverlapOut  *overlapout,
		  InspiralWaveOverlapIn   *overlapin,
		  REAL4Vector             *Filter1, 
		  REAL4Vector             *Filter2,
		  BCVMaximizationMatrix    matrix,
		  OtherParamIn             otherIn);

/* Function to store the moments needed by the BCV overlap process 
 * a11, a22, a21 are three vectors which stored the three components 
 * of the matrix needed in the BCV maximization process. 
 * */
void
LALCreateMomentVector(LALStatus            *status,
		      REAL4Vector          *a11,
		      REAL4Vector          *a21,
		      REAL4Vector          *a22,
		      REAL8FrequencySeries *psd,
		      InspiralTemplate     *params);

/* Function to create Vectors of the form y(f)= f^(-a/b). 
 * where f is the frequency. 
 * */
void 
LALCreateVectorFreqPower( REAL4Vector	 	*vector,
			  InspiralTemplate 	params,
			  int 			a,
			  int 			b);


/* function to generate waveform (has to be tested). 
 * */
void
LALGenerateWaveform(LALStatus 			*status,
		    REAL4Vector 		*signal,
		    RandomInspiralSignalIn 	*params);

/* print a time series 
 * */
void
printf_timeseries(INT4 		n,
		  REAL4 	*signal,
		  REAL8 	delta,
		  REAL8 	t0);


/* Init the CoarseBank structure 
 * */
void InitInspiralCoarseBankIn(InspiralCoarseBankIn 	*coarseIn);
/* Init the random structure 
 * */
void InitRandomInspiralSignalIn(RandomInspiralSignalIn 	*randIn);
/* Init the OtherParamIn structure
 * */
void InitOtherParamIn(OtherParamIn 			*otherIn);
/* Function to initialize the different strucutre */
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseIn, 
				RandomInspiralSignalIn	*randIn, 
				OtherParamIn		*otherIn);
/* Parsing function 
 * */
void
ParseParameters(int 			*argc, 
		char 			**argv,
		InspiralCoarseBankIn 	*coarseIn,
		RandomInspiralSignalIn 	*randIn,
		OtherParamIn 		*otherIn);		     

/* function to check validity of some parameters
 * */
void 	
CheckParams(InspiralCoarseBankIn 	coarseIn,
	    RandomInspiralSignalIn 	randIn,
	    OtherParamIn 		otherIn);


/* Default values
 * */
void 
SetDefault(InspiralCoarseBankIn 	*coarseIn, 
	   RandomInspiralSignalIn 	*randIn,
	   OtherParamIn 		*otherIn);

/* Help Function 
 * */
void 	
Help(	InspiralCoarseBankIn 	coarseIn,
	RandomInspiralSignalIn 	randIn,
	OtherParamIn 		otherIn);

/* Print Function for final results 
 * */
void
PrintResults(InspiralTemplate       	bank,
	    InspiralTemplate       	injected,
	    InspiralWaveOverlapOut 	overlapout,
	    REAL8 			a,
	    UINT4 			b);

/* Print each  template bank overlap 
 * */
void
PrintBankOverlap(InspiralTemplateList 	**list,
		 int 			nlist,
		 float 			*overlap,
		 InspiralCoarseBankIn 	coarseIn);

/* Print the template bank coordinates in ascii format 
 * */
void BEPrintBank(InspiralCoarseBankIn 	coarse, 
		 InspiralTemplateList 	**list,
		 UINT4 			nlist);

/* print the template bank coordinates  in xml format 
 * */
void BEPrintBankXml(InspiralTemplateList *coarseList,
		    UINT4 		 numCoarse,
		    InspiralCoarseBankIn   coarseIn,
		    RandomInspiralSignalIn randIn,
		    OtherParamIn           otherIn
		    );

/* get the matrix involved in alpha maximization 
 * process in the BCV correlation.
 * */
void BEGetMatrixFromVectors(REAL4Vector A11,
			    REAL4Vector A21,
			    REAL4Vector A22,
			    UINT4 	k,
			    BCVMaximizationMatrix *matrix2fill);


/* print an error  message 
 * */
void BEPrintError(char *chaine);

void BEFillProc(ProcessParamsTable     *proc,
		InspiralCoarseBankIn   coarseIn,
		RandomInspiralSignalIn randIn,
		OtherParamIn           otherIn);
void 
BEPrintResultsXml(InspiralCoarseBankIn   coarseIn,
		  RandomInspiralSignalIn randIn,
		  OtherParamIn           otherIn, 
		  ResultIn inject
		  );
void 
BEPrintProtoXml(InspiralCoarseBankIn   coarseIn,
		  RandomInspiralSignalIn randIn,
		  OtherParamIn           otherIn
		);

void
BEReadXmlBank(  LALStatus  *status, 
		CHAR *bankFileName, 
		InspiralTemplateList **list,
		INT4 *nlist, 
		InspiralCoarseBankIn coarseIn);



/* 
 * ==================== THE MAIN PROGRAM ==================== 
 * */
int
main (INT4 argc, CHAR **argv ) 
{
  /* --- Variables ---*/
  /*comments on variables will come later */
  INT4		nn;
  INT4  	nby2;
  INT4  	l=0;
  INT4  	lmax;
  INT4 		nbank	= 0;
  INT4		j	= 0;
  INT4 		n;
  INT4 		k, kk;
  INT4 		jmax;
  INT4  	nlist;
  INT4  	kMin 	= 0;
  INT4  	temporder;
  UINT4         ntrials = 0;
  REAL4 	temp;
  UINT4 	i;
    
  

  INT8 		il;
  
  REAL8		df;
  REAL8		fendBCV	= 0;
  REAL8 	fMax;
  
  REAL4Vector   signal;
  REAL4Vector   correlation;
  REAL4Vector   VectorPhase;
  REAL4Vector   VectorPowerFm5_3;
  REAL4Vector   VectorPowerFm2_3;
  REAL4Vector   VectorPowerFm7_6;
  REAL4Vector   VectorPowerFm1_2;
  REAL4Vector   VectorA11;
  REAL4Vector   VectorA21;
  REAL4Vector   VectorA22;
  REAL4Vector   Filter1; 
  REAL4Vector   Filter2;
  
  ResultIn      Result;

  RandomInspiralSignalIn 	randIn;						/* random signal waveform to inject	*/

  InspiralWaveOverlapIn 	overlapin;					/* structure for Overlap		*/

  InspiralWaveOverlapOut     	overlapout; 
  InspiralWaveOverlapOut	overlapoutmax;

  OtherParamIn 			otherIn;					/**/
  
  RealFFTPlan 			*fwdp = NULL;
  RealFFTPlan 			*revp = NULL;
  
  InspiralTemplateList      	*list = NULL;

  InspiralCoarseBankIn      	coarseBankIn; 					/* strcture for the bank of templates	*/

  static LALStatus 		status;

  BCVMaximizationMatrix 	matrix;

  FILE				*Finput;   					/* to read the input psd 		*/
  FILE                          *Foutput;					/* to print some output in files	*/

  float 			*bankEfficiencyOverlapIn;			
  /* --- main code start here --- */ 
  /* --- Debugging level --- */

  lalDebugLevel = 0;
 
  /* --- Some initialization --- */
  
  ParametersInitialization(	&coarseBankIn, 
		  		&randIn,
			       	&otherIn);					/* Initialization of structure		*/

  ParseParameters(		&argc,
		  		argv,
			       	&coarseBankIn,
			       	&randIn,
			       	&otherIn);					/* Read Parameters from user 		*/

  CheckParams(			coarseBankIn,
		 		randIn,
				otherIn);					/* Check validity of some variables. 	*/



  if (otherIn.PrintPrototype){
    BEPrintProtoXml(coarseBankIn, randIn, otherIn);    
    exit(0);
  }
  
  /* --- Estimate size of the signal --- */
  randIn.param.massChoice 	= m1Andm2;					/* Only to compute the length of "signal"*/ 
  signal.length 		= 0.;
  randIn.param.approximant 	= EOB;						/* Only to compute the length of "signal"*/
  LAL_CALL(LALInspiralWaveLength(&status, 
			  	 &signal.length, randIn.param), &status);
  randIn.param.approximant = otherIn.signal;					/* Retrieve the actual approximant of the 
										   waveform to inject. otherwise it will
										   be EOB all the time. */
  /* --- now that we have the maiximum length of the signal we can go ahead --- */
  signal.length*=otherIn.lengthFactor;
  /* --- Allocate memory --- */
  correlation.length 	= signal.length;					/* affect length of the correlations	*/
  randIn.psd.length 	= signal.length/2 + 1;					/* as well as random data length	*/
  signal.data 		= (REAL4*) LALCalloc(1, sizeof(REAL4)*signal.length);	/* memory for the injected signal	*/
  correlation.data 	= (REAL4*) LALCalloc(1, sizeof(REAL4)*correlation.length); /* and correlation 			*/

  memset( &(coarseBankIn.shf), 0, sizeof(REAL8FrequencySeries) );			/* ?? coarseBankIn. should be after ??	*/

  coarseBankIn.shf.f0 	= 0;
  LAL_CALL(LALDCreateVector( &status, &(coarseBankIn.shf.data), randIn.psd.length ),&status);
  
  coarseBankIn.shf.deltaF 	= randIn.param.tSampling / signal.length;
  randIn.psd.data 	= (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length); /* randin allocation			*/
  
  /* --- Compute Noise Spectral Density --- */
  df = randIn.param.tSampling/(float) signal.length;				/* frequency resolution			*/
  switch (otherIn.NoiseModel)							/* model choice				*/
    {
    case LIGOI:  
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseBankIn.shf.data, &LALLIGOIPsd, df), &status);
      break;
    case LIGOA:  
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseBankIn.shf.data, &LALAdvLIGOPsd, df), &status);
      break;
    case VIRGO: 
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseBankIn.shf.data, &LALVIRGOPsd, df), &status);;
      break;
    case GEO:
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseBankIn.shf.data, &LALGEOPsd, df), &status);
      break;
    case TAMA:
      LAL_CALL(LALNoiseSpectralDensity (&status, coarseBankIn.shf.data, &LALTAMAPsd, df), &status);
      break;
    case REALPSD:								/* has to be done properly. */
      /* Read psd dans le fichier InputPSD.dat */
      Finput = fopen(otherIn.filename,"r");
      for (i=0; i<coarseBankIn.shf.data->length; i++)
	{
	  /* factor DeltaT since the input data have a sampling frequency of 1./DeltaT */
	  fscanf(Finput, "%f %lf\n", &temp,&coarseBankIn.shf.data->data[i]);
	  for (kk=1; kk < (int)( DeltaT * df); kk++)
	      fscanf(Finput, "%f %f\n", &temp,&temp);
	}
      fclose(Finput);
      break;
     }


  /* --- save the psd in the file psd.dat if requested --- */
  if (otherIn.PrintPsd)
   {
     Foutput= fopen(BANKEFFICIENCY_PRINTPSD_FILE,"w");
     for (i=1; i<coarseBankIn.shf.data->length; i++)
       fprintf(Foutput, "%f %e\n",(float)i*df, coarseBankIn.shf.data->data[i]);  
     fclose(Foutput);
   }

  /* --- and affect copy coarseBankIn.psd in randIn.psd --- */  
  for (i=0; i< randIn.psd.length; i++){
    randIn.psd.data[i] = coarseBankIn.shf.data->data[i];
  }
  
  
  
  /* --- And the bank of templates --- */
  
   temporder  = coarseBankIn.order;
   coarseBankIn.order = 4; 			/* force to be 2PN otherwise the bank 
						   is not generated. has to be fixed 
						   in the bank itself*/
   
   LAL_CALL(LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseBankIn),&status);
   if (nlist==0) exit(0);	
   bankEfficiencyOverlapIn = (float*) malloc(sizeof(float*)* nlist);              	/* allocate memory to store the bank 	*/
   for (il=0; il<nlist; il++) 								/* get back the order 			*/
     (list[il]).params.order = temporder;

   /* Do we want to read an xml template bank from a file ? */
   if (otherIn.inputXMLBank){
     BEReadXmlBank(&status, otherIn.inputXMLBank, &list, &nlist, coarseBankIn );
   }

  /* --- do we want to print the bank ? --- */
  if (otherIn.PrintBank) BEPrintBank(coarseBankIn, &list, nlist);				/* ascii format				*/
  if (otherIn.PrintBankXml) BEPrintBankXml(list,  nlist, coarseBankIn, randIn, otherIn);	/* xml format				*/

  /* --- Estimate the fft's plans --- */
  LAL_CALL(LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0), &status);
  LAL_CALL(LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0), &status);
  
  
  /* --- The overlap structure --- */
  overlapin.nBegin 	= 0;
  overlapin.nEnd 	= 0;
  overlapin.psd 	= randIn.psd;
  overlapin.fwdp 	= randIn.fwdp = fwdp;
  overlapin.revp 	= revp;
  
  /* If ambiguity function flag equal to 1 we don't want to compute 
     tons of trials since th number of ouput could be too important*/
  if (otherIn.ambiguityFunction  )      otherIn.ntrials=1;

 


  list[1].params.startPhase = 0;   				/* why ? */
  
  /* --- Compute the frequency vectors and others --- */
  VectorPhase.length      = signal.length / 2;
  VectorA11.length        = signal.length / 2 ;
  VectorA22.length        = signal.length / 2;
  VectorA21.length        = signal.length / 2 ;
  VectorPowerFm2_3.length = signal.length / 2 ;
  VectorPowerFm5_3.length = signal.length / 2;       
  VectorPowerFm7_6.length = signal.length / 2 ;
  VectorPowerFm1_2.length = signal.length / 2;       
   
  Filter1.length          = correlation.length;
  Filter2.length          = correlation.length;
  
  VectorPhase.data        = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPhase.length);   
  VectorA11.data          = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorA11.length);
  VectorA21.data          = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorA21.length);
  VectorA22.data          = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorA22.length);
  VectorPowerFm2_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPowerFm2_3.length);
  VectorPowerFm5_3.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPowerFm5_3.length);
  VectorPowerFm7_6.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPowerFm7_6.length);
  VectorPowerFm1_2.data   = (REAL4*) LALCalloc(1, sizeof(REAL4) * VectorPowerFm1_2.length);
  Filter1.data            = (REAL4*) LALCalloc(1, sizeof(REAL4) * Filter1.length);
  Filter2.data            = (REAL4*) LALCalloc(1, sizeof(REAL4) * Filter2.length);
   
   
  n 	= signal.length / 2;
  df 	= randIn.param.tSampling / (REAL8)n / 2.;      
  
  

  /* --- Create useful vector once for all --- */
  LAL_CALL(LALCreateVectorFreqPower(&VectorPowerFm5_3, randIn.param, -5, 3), &status);
  LAL_CALL(LALCreateVectorFreqPower(&VectorPowerFm2_3, randIn.param, -2, 3), &status);
  LAL_CALL(LALCreateVectorFreqPower(&VectorPowerFm7_6, randIn.param, -7, 6), &status);
  LAL_CALL(LALCreateVectorFreqPower(&VectorPowerFm1_2, randIn.param, -1, 2), &status);

  
  /* --- Create Matrix for BCV Maximization once for all --- */
  LAL_CALL(LALCreateMomentVector(&status, 
				 &VectorA11,
				 &VectorA21,
				 &VectorA22,
				 &coarseBankIn.shf, 
				 &(list[j].params)), &status);
  
  /* If check option is on, we compute the best overlap 
   * (we should get one if template and approximant are the same) 
   * otherwise we jump to the else with the bank efficiency itself
   * */   
  if (otherIn.check) {
    randIn.param.approximant    = otherIn.signal;  				/* The waveform parameter for injection */
    randIn.param.fCutoff 	= coarseBankIn.fUpper; 				/* its cutoff frequency 		*/
    /* What kind of parameters do we use ?*/
    if (otherIn.signal == BCV) 
      randIn.param.massChoice = psi0Andpsi3;
    else
      randIn.param.massChoice = m1Andm2;					/* might be totalMass ans eta		*/
    /* Let's compute the random parameters of the waveform to inject*/
    for (i=0; i<signal.length; i++) signal.data[i] = 0.;       
    
    /* we might force to compute a non random waveform giving the two
       input parameters m1 and m2 */
    if (otherIn.m1==-1 && otherIn.m2 ==-1){
      LAL_CALL(LALRandomInspiralSignal(&status, &signal, &randIn), &status);
    }
    else
      {
	randIn.param.mass1 = otherIn.m1;
	randIn.param.mass2 = otherIn.m2;
	LALGenerateWaveform(&status,&signal, &randIn);
      }
    overlapin.signal 	= signal;
    jmax           	= 0;
    overlapout.max 	= 0.0;
    kMin  		= floor(randIn.param.fLower / df );
    
    overlapin.param 	= randIn.param;
    list[1].params      = randIn.param;
    
    /*
      overlapin.param.fFinal  = list[1].params.fFinal;
      overlapin.param.fCutoff = list[1].params.fFinal;
    */
    k = floor(overlapin.param.fFinal / df);
 
    BEGetMatrixFromVectors(VectorA11, VectorA21, VectorA22, k, &matrix);	/* fill the bcv matrix */
    switch(otherIn.overlapMethod)
      {
      case AlphaMaximization:
	LAL_CALL(LALCreateFilters(&Filter1,
				  &Filter2,
				  VectorPowerFm5_3,
				  VectorPowerFm2_3,
				  VectorPowerFm7_6,
				  VectorPowerFm1_2,
				  matrix,
				  kMin,					       
				  list[1].params.psi0,
				  list[1].params.psi3
				  ), &status);
	/* The overlap given two filters and the input signal*/
	LAL_CALL(LALWaveOverlapBCV(&status, 
				   &correlation,
				   &overlapout,
				   &overlapin,
				   &Filter1,
				   &Filter2,
				   matrix,
				   otherIn), &status);
	/* given overlapout, j, l, fendBCV then we fill overlapoutmax, 
	 * jmax, lmax, fMax.
	 * in that case (checking), j and l are not importamt though. */	
	KeepHighestValues(overlapout, j, l, fendBCV, &overlapoutmax,  &jmax, &lmax, &fMax);
	break;
      case InQuadrature:
	{ 
	  overlapin.param.fCutoff = 1023;
	  for (i=0; i<signal.length; i++) 
		correlation.data[i] = 0.;	   	   	  	 
	  
	  LAL_CALL(LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin), &status);
	  KeepHighestValues(overlapout, j, l, fendBCV, &overlapoutmax,  &jmax, &lmax, &fMax);
	}
      }
    
    printf("#");
    PrintResults(list[1].params, randIn.param,overlapoutmax,  list[1].params.fFinal, 1);

    /* --- if requested we print the correlation in a file--- */
    if (otherIn.PrintOverlap){
     Foutput  = fopen(BANKEFFICIENCY_PRINTOVERLAP_FILE,"w");
     for ( i = 0; i < correlation.length; i++)
	fprintf(Foutput,"%f %f\n", i*df, (correlation.data[i]));
     fclose(Foutput);
    }
  }
  else{										 /* --- the real code is here --- */
    /* --- The main loop --- */
    
    if (otherIn.PrintBankOverlap) /*we just create the empty file here if requested why ?? */
      {
	Foutput = fopen("FF.sr4","w");
	fclose(Foutput); 
      }
    
    while (++ntrials <= otherIn.ntrials) 
      {
        randIn.param.approximant    	= otherIn.signal;  			/* The waveform parameter for injection */
	randIn.param.fCutoff 		= coarseBankIn.fUpper; 			/* its cutoff frequency 		*/
	/* What kind of parameters do we use ?*/
	if (otherIn.signal == BCV) 
	  randIn.param.massChoice = psi0Andpsi3;
	else
	{
          randIn.param.massChoice = m1Andm2;
	   randIn.param.massChoice = totalMassUAndEta;
	}
	/* Let's compute the random parameters of the waveform to inject*/
	 for (i=0; i<signal.length; i++) signal.data[i] = 0.;       
	 
	 /* we might force to compute a non random waveform giving the two
	    input parameters m1 and m2 */


	 if (otherIn.m1==-1 && otherIn.m2 ==-1)
	   {   
	     if (otherIn.binaryInjection == BHNS){	 
	       randIn.param.massChoice = bhns;
	       LAL_CALL(LALRandomInspiralSignal(&status, &signal, &randIn), &status);
	     }
	     else
	       LAL_CALL(LALRandomInspiralSignal(&status, &signal, &randIn), &status);
	     
	   }
	 else
	   {
	     randIn.param.mass1 = otherIn.m1;
	     randIn.param.mass2 = otherIn.m2;
	     LAL_CALL(LALGenerateWaveform(&status,&signal, &randIn), &status);
	   }
	     
	 overlapin.signal 	= signal;
	 jmax           	= 0;
	 overlapoutmax.max 	= 0.0;

	 
	 /* we might force to use only one template giving the input psi0 and psi3. 
	    Then, we don't need to use the bank*/
	 if (otherIn.psi0 != -1 && otherIn.psi3 != -1)
	   {
	     list[0].params.psi0 = otherIn.psi0;
	     list[0].params.psi3 = otherIn.psi3;
	     nlist  = 1;
	   }
	 
	 
	 nby2 	= Filter1.length / 2;
	 nn 	= Filter1.length;
	 kMin  	= floor(randIn.param.fLower /df );
	 df   	= randIn.param.tSampling / (REAL8)(VectorPhase.length) / 2.;
	 

	 /* --- we can process the data through the bank now --- */
	 nbank = 0;
	 for (j = 0; j < nlist; j++) 
	   {	   	     	     	     
	     overlapin.param 	= list[j].params;	   


	     overlapout.max = -11.11; 
	  

	     if (
		 (otherIn.FastSimulation == 1 ) &&
		 ((fabs(randIn.param.t0 - list[j].params.t0)> .1) ||
		 (fabs(randIn.param.t3 - list[j].params.t3)> .1) )
		 )
	       {
		 /*fprintf(stderr, "Fast Simulation: skipped template number %d t0=%lf t3=%lf\n ", j, list[j].params.t0, list[j].params.t3);*/
	       }
	     else
	       {
	     
		 switch(otherIn.overlapMethod)
		   {
		   case AlphaMaximization:	   	   	   
		     nbank++;
		     fendBCV = list[j].params.fFinal;
		     overlapin.param.fFinal   =  fendBCV;
		     overlapin.param.fCutoff  =  fendBCV;
		     
		     /* Extract some values for template creation*/
		     k = floor(fendBCV / df);
		     BEGetMatrixFromVectors(VectorA11, VectorA21, VectorA22, k, &matrix);
		     
		     /* Template creation*/
		     LAL_CALL(LALCreateFilters(&Filter1,
					       &Filter2,
					       VectorPowerFm5_3,
					       VectorPowerFm2_3,
					       VectorPowerFm7_6,
					       VectorPowerFm1_2,
					       matrix,
					       kMin,					       
					       list[j].params.psi0,
					       list[j].params.psi3
					       ), &status);
		     
		     /* The overlap given two filters and the input signal*/
		     LAL_CALL(LALWaveOverlapBCV(&status, 
						&correlation, &overlapout, &overlapin,
						&Filter1, &Filter2, matrix, otherIn), 
			      &status);
		     
		     break;
		   case InQuadrature:
		     if (coarseBankIn.approximant == BCV) /* bcv template but no apha maximization */
		       {
			 fendBCV  			=  list[j].params.fFinal;  
			 overlapin.param.fCutoff 	= fendBCV;
			 overlapin.param.fFinal 	= fendBCV;
			 for (i=0; i<signal.length; i++) correlation.data[i] = 0.;	   	   	  
			 LAL_CALL(LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin), &status);
		       }    
		     else	/* --- classical version of overlap using non bcv templates. 
				   thus flso has to  be fixed manually --- */
		       {
			 /* SHOULD be replace by flso of the template in the bank*/
			 

			 fendBCV   = 1./LAL_PI/pow(6, 1.5)/(list[j].params.mass1+list[j].params.mass2)/LAL_MTSUN_SI;
			 if (coarseBankIn.approximant==EOB && coarseBankIn.order==threePN)
			   fendBCV   = 1./LAL_PI/pow(2.3, 1.5)/(list[j].params.mass1+list[j].params.mass2)/LAL_MTSUN_SI;
			 
			 

			 if (fendBCV > randIn.param.tSampling/2 ) 
			   fendBCV = randIn.param.tSampling/2. - 1;
			 
			 
			 
			 
			 overlapin.param.fFinal  = fendBCV;		    
			 overlapin.param.fCutoff = fendBCV;
		overlapin.param = randIn.param;	 
			 for (i=0; i<signal.length; i++) correlation.data[i] = 0.;	   	   	  
			 LAL_CALL(LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin), &status);

		       }
		     
		     
		     overlapout.alpha = -1;	       /* since alpha is not initalized in that switch we give a  constante value */
		     break;	     
		   }     /* switch end */
		 
		 /* --- if overlap is the largest one, then we keep some 
		  * other information . Later we might just use the previous 
		  * vector. Keep everything and use it outside the bank 
		  * process --- */	     
		 bankEfficiencyOverlapIn[j] = overlapout.max; 
		 /* --- if overlap is the largest one, then we keep some 
		  * other information . Later we might just use the previous 
		  * vector. Keep everything and use it outside the bank 
		  * process --- */	     
		 KeepHighestValues(overlapout, j, l, fendBCV, &overlapoutmax,  &jmax, &lmax, &fMax);
		
	       }/*end of if cheating code*/  

		 
	       }/*end of the bank process*/

	   
	 
	 /* --- Now print some results or other stuffs --- */
	 /* --- first, do we want to print the bank results (ambiguity function)*/
	 if (otherIn.PrintBankOverlap) 
	   PrintBankOverlap(&list, nlist ,  bankEfficiencyOverlapIn, coarseBankIn);
	
	 
	 /* Then print the maximum overlap over the whole bank and the corresponding 
	  * parameter of the templates, injected signal and so on. This is the main results
	  * We'll print one line for each simulation (ntrials parameter)*/
	 LAL_CALL(LALInspiralParameterCalc(&status, &list[jmax].params), &status);
	 if (!otherIn.quietFlag)
	   PrintResults(list[jmax].params, randIn.param, overlapoutmax, fMax, list[jmax].nLayer);


	 

	 if (otherIn.PrintResultXml)
	   {
	     if (otherIn.overlapMethod!=InQuadrature)
	     {
		     Result.psi0_trigger = list[jmax].params.psi0;
		     Result.psi3_trigger = list[jmax].params.psi3;
		     Result.psi0_inject  = randIn.param.psi0;
		     Result.psi3_inject  = randIn.param.psi3;
	     }
	     else
	     {
		     Result.psi0_trigger = list[jmax].params.t0;
		     Result.psi3_trigger = list[jmax].params.t3;
		     Result.psi0_inject  = randIn.param.t0;
		     Result.psi3_inject  = randIn.param.t3;
	     }   
	     Result.mass1_inject = randIn.param.mass1;
	     Result.mass2_inject = randIn.param.mass2;
	     Result.fend_trigger = fendBCV;
	     Result.fend_inject = randIn.param.fFinal;
	     Result.totalMass_trigger = list[jmax].params.totalMass;
	     Result.eta_trigger = list[jmax].params.eta;
	     Result.rho_final   = overlapoutmax.max; 
	     Result.alpha = overlapoutmax.alpha;
	     Result.alpha_f = overlapoutmax.alpha* pow(fendBCV,2./3.);; 
	     Result.bin = overlapoutmax.bin;
	     Result.phase = overlapoutmax.phase;
	     Result.layer = list[jmax].nLayer;
	     Result.ntrial = ntrials;
	     BEPrintResultsXml(coarseBankIn, randIn, otherIn, Result);	 
	   }
	 
	 /* --- here we might do something else which is not needed but 
	  * it is a good check . We get the best templates given by the bank 
	  * process and vary the last frequency --- */
	 fflush(stdout);	


	 /* now some studies related to the probability of detection */



      }

  }


   /* --- destroy the plans, correlation and signal --- */
   
   LAL_CALL(LALFree(VectorPowerFm5_3.data), &status);
   LAL_CALL(LALFree(VectorPowerFm2_3.data), &status);
   LAL_CALL(LALFree(VectorPowerFm1_2.data), &status);
   LAL_CALL(LALFree(VectorPowerFm7_6.data), &status);
   LAL_CALL(LALFree(VectorPhase.data), &status);;  
   LALFree(VectorA11.data);
   LALFree(VectorA21.data);
   LALFree(VectorA22.data);
   
   LALFree(Filter2.data);
   LALFree(Filter1.data);
   
   LALFree(randIn.psd.data);
   LALDDestroyVector( &status, &(coarseBankIn.shf.data) );
   
   LALFree(signal.data);
   LALFree(correlation.data);
   LALFree(list);
   
   LALDestroyRealFFTPlan(&status,&fwdp);
   LALDestroyRealFFTPlan(&status,&revp);
   LALCheckMemoryLeaks();    
   
   
   return 0;
}


/* ****************************************************************************
 * Function to print a real4 vector in a file (wave1.dat)
 * ****************************************************************************/
void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
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

/* ****************************************************************************
 * Initialization of the parameters fromn random signal, bank and others 
 * structures
 * ***************************************************************************/
void ParametersInitialization(	InspiralCoarseBankIn 	*coarseBankIn, 
				RandomInspiralSignalIn	*randIn, 
				OtherParamIn		*otherIn)
{
  InitInspiralCoarseBankIn(coarseBankIn);
  InitRandomInspiralSignalIn(randIn);
  InitOtherParamIn(otherIn);
}

/* ****************************************************************************
 * CoarseIn initialization
 ***************************************************************************/
void InitInspiralCoarseBankIn(InspiralCoarseBankIn *coarseBankIn)
{
  coarseBankIn->fLower    	= BANKEFFICIENCY_FLOWER;
  coarseBankIn->fUpper    	= BANKEFFICIENCY_TSAMPLING/2. - 1.;
  coarseBankIn->tSampling 	= BANKEFFICIENCY_TSAMPLING;
  coarseBankIn->space     	= BANKEFFICIENCY_SPACE;
  coarseBankIn->mmCoarse  	= BANKEFFICIENCY_MMCOARSE;
  coarseBankIn->mmFine    	= BANKEFFICIENCY_MMFINE;
  coarseBankIn->iflso       = BANKEFFICIENCY_IFLSO ;
  coarseBankIn->mMin      	= BANKEFFICIENCY_MMIN;
  coarseBankIn->mMax      	= BANKEFFICIENCY_MMAX;
  coarseBankIn->MMax      	= BANKEFFICIENCY_MMAX * 2;
  coarseBankIn->massRange   = MinMaxComponentMass; 
  coarseBankIn->etamin 	= coarseBankIn->mMin * coarseBankIn->mMax  
			/ pow(coarseBankIn->MMax, 2.);

  coarseBankIn->psi0Min   	= BANKEFFICIENCY_PSI0MIN;
  coarseBankIn->psi0Max   	= BANKEFFICIENCY_PSI0MAX;
  coarseBankIn->psi3Min   	= BANKEFFICIENCY_PSI3MIN;
  coarseBankIn->psi3Max   	= BANKEFFICIENCY_PSI3MAX;
  coarseBankIn->alpha     	= BANKEFFICIENCY_ALPHABANK;
  /*unsigned variable ... */
  coarseBankIn->numFcutTemplates = BANKEFFICIENCY_NFCUT;
  coarseBankIn->approximant = BANKEFFICIENCY_TEMPLATE;  
  coarseBankIn->order       = BANKEFFICIENCY_ORDER_TEMPLATE;  
  coarseBankIn->LowGM     	= BANKEFFICIENCY_LOWGM;
  coarseBankIn->HighGM    	= BANKEFFICIENCY_HIGHGM;
}

/* ****************************************************************************
 * Random initialization 
 * ****************************************************************************/
void InitRandomInspiralSignalIn(RandomInspiralSignalIn *randIn)
{
  randIn->useed         	= BANKEFFICIENCY_USEED;				/* seed for random MC simulation 	*/ 
  randIn->type          	= BANKEFFICIENCY_TYPE;				/* type of simulation 			*/
  randIn->SignalAmp     	= BANKEFFICIENCY_SIGNALAMPLITUDE;		/* if x = n +h 				*/
  randIn->param.order   	= BANKEFFICIENCY_ORDER_SIGNAL;   		/* and its order			*/
  randIn->param.alpha		= BANKEFFICIENCY_ALPHASIGNAL; 			/* alpha value of BCV 			*/
  randIn->param.ieta    	= BANKEFFICIENCY_IETA; 				/*					*/
  randIn->param.mass1 		= BANKEFFICIENCY_MMIN;				/* To aalocate memory			*/ 
  randIn->param.mass2  		= BANKEFFICIENCY_MMIN; 				/* idem 				*/
  randIn->param.fLower  	= BANKEFFICIENCY_FLOWER;			/* Lower cutoff freq. of the template	*/
  randIn->param.OmegaS		= 0.;						/* EOB parameter 			*/
  randIn->param.Theta		= 0.;						/* EOB parameter 			*/
  randIn->mMin          	= BANKEFFICIENCY_MMIN;				/* min mass to inject			*/
  randIn->mMax          	= BANKEFFICIENCY_MMAX;				/* max mass to inject 			*/
  randIn->MMax          	= BANKEFFICIENCY_MMAX * 2;			/* total mass max 			*/
  randIn->etaMin        = (BANKEFFICIENCY_MMIN * BANKEFFICIENCY_MMAX) 
    /BANKEFFICIENCY_MMAX/ BANKEFFICIENCY_MMAX;
  randIn->psi0Min  		= BANKEFFICIENCY_PSI0MIN;			/* psi0 range 				*/
  randIn->psi0Max  		= BANKEFFICIENCY_PSI0MAX;	
  randIn->psi3Min  		= BANKEFFICIENCY_PSI3MIN;	
  randIn->psi3Max  		= BANKEFFICIENCY_PSI3MAX;	
  randIn->param.approximant 	= BANKEFFICIENCY_SIGNAL;			/* approximant of h, the signal		*/
  randIn->param.tSampling 	= BANKEFFICIENCY_TSAMPLING;			/* sampling 				*/
  randIn->param.fCutoff		= BANKEFFICIENCY_TSAMPLING/2.-1;		/* sampling 				*/
  randIn->param.startTime       = BANKEFFICIENCY_STARTTIME; 
  randIn->param.startPhase      = BANKEFFICIENCY_STARTPHASE; 
  randIn->param.nStartPad       = BANKEFFICIENCY_NSTARTPHASE;
  randIn->param.signalAmplitude = BANKEFFICIENCY_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKEFFICIENCY_NENDPAD;
  randIn->NoiseAmp              = BANKEFFICIENCY_NOISEAMPLITUDE;
}

/* ****************************************************************************
 * Other Param initialization							
 * ***************************************************************************/
void InitOtherParamIn(OtherParamIn *otherIn)
{
  otherIn->template		= BANKEFFICIENCY_TEMPLATE;
  otherIn->bank         	= -1;
  otherIn->signal       	= BANKEFFICIENCY_SIGNAL;
  otherIn->m1           	= -1;
  otherIn->m2           	= -1;
  otherIn->lengthFactor	= 1;
  otherIn->psi0         	= -1;
  otherIn->psi3         	= -1;
  otherIn->PrintOverlap 	= BANKEFFICIENCY_PRINTOVERLAP;
  otherIn->PrintFilter  	= BANKEFFICIENCY_PRINTFILTER;
  otherIn->PrintPsd             = BANKEFFICIENCY_PRINTPSD;
  otherIn->overlapMethod	= AlphaMaximization;

  otherIn->PrintBank    	= BANKEFFICIENCY_PRINTBANK;
  otherIn->PrintBankXml    	= BANKEFFICIENCY_PRINTBANKXML;
  otherIn->PrintResultXml    	= BANKEFFICIENCY_PRINTRESULTXML;
  otherIn->PrintPrototype    	= BANKEFFICIENCY_PRINTPROTOTYPE;

  otherIn->PrintBankOverlap	= BANKEFFICIENCY_PRINTBANKOVERLAP;
  otherIn->PrintTemplate	= BANKEFFICIENCY_PRINTTEMPLATE; 
  otherIn->check        	= BANKEFFICIENCY_CHECK;
  otherIn->inputPSD          	= NULL;
  otherIn->quietFlag 	     	= BANKEFFICIENCY_QUIETFLAG;
  otherIn->ambiguityFunction 	= BANKEFFICIENCY_AMBIGUITYFUNCTION;
  otherIn->ntrials 	     	= BANKEFFICIENCY_NTRIALS;
  otherIn->FMaximization     	= BANKEFFICIENCY_FMAXIMIZATION;
  otherIn->FastSimulation       = BANKEFFICIENCY_FASTSIMULATION;
  otherIn->NoiseModel           = LIGOI;
  otherIn->binaryInjection=NoUserChoice; 
  otherIn->inputXMLBank = NULL;
}



/* ***************************************************************************
 *  Function to parse user parameters 
 *  ************************************************************************ */
void 
ParseParameters(	INT4 			*argc, 
			CHAR 			**argv,
			InspiralCoarseBankIn 	*coarseBankIn,
			RandomInspiralSignalIn 	*randIn,
			OtherParamIn    	*otherIn)
{
  INT4 		i = 1;


  while(i < *argc)
    {
      if ( strcmp(argv[i],	"--fend-bcv") 	== 0 ) {
		coarseBankIn->LowGM      	= atof(argv[++i]);
       	      	coarseBankIn->HighGM     	= atof(argv[++i]);
	}
      else if ( strcmp(argv[i],	"--fl") 	== 0 ) 
  		coarseBankIn->fLower = randIn->param.fLower	= atof(argv[++i]);  
      else if ( strcmp(argv[i],	"--sampling") 	== 0 ) {
  		coarseBankIn->tSampling = randIn->param.tSampling = atof(argv[++i]);  
		randIn->param.fCutoff 	= coarseBankIn->tSampling/2. - 1.;
		coarseBankIn->fUpper 	= coarseBankIn->tSampling/2. - 1.;

	}      
      else if ( strcmp(argv[i],	"--mass-range-bank")	== 0 ) 	{	
		coarseBankIn->mMin = atof(argv[++i]);
		coarseBankIn->mMax = atof(argv[++i]); 
		coarseBankIn->MMax = coarseBankIn->mMax * 2.; 
		coarseBankIn->etamin 	= coarseBankIn->mMin *  coarseBankIn->mMax 
		  / pow(coarseBankIn->MMax, 2.);
      }
      else if ( strcmp(argv[i],	"--mass-range")	== 0 ) 	{	
		randIn->mMin = atof(argv[++i]);
		randIn->param.mass1 = randIn->param.mass2 = randIn->mMin;
		randIn->mMax = atof(argv[++i]); 
		randIn->MMax = randIn->mMax * 2.; 
		randIn->etaMin 	= randIn->mMin * randIn->mMax 
			/ pow(randIn->MMax, 2.);
	}
      else if ( strcmp(argv[i], "--m1") 	== 0 )
	otherIn->m1 = atof(argv[++i]);	      
      else if ( strcmp(argv[i], "--m2") 	== 0 )  	
	      otherIn->m2 = atof(argv[++i]); 	       
      else if ( strcmp(argv[i], "--psi0") 	== 0 )	
	      otherIn->psi0 = atof(argv[++i]);  
      else if ( strcmp(argv[i], "--psi3") 	== 0 ) 
	      otherIn->psi3 = atof(argv[++i]); 
      else if (strcmp(argv[i],	"--psi0-range")	==0) {
		coarseBankIn->psi0Min = randIn->psi0Min = atof(argv[++i]);
		coarseBankIn->psi0Max = randIn->psi0Max = atof(argv[++i]);
	}
      else if (strcmp(argv[i],	"--psi3-range")	==0){
		coarseBankIn->psi3Max = randIn->psi0Max = atof(argv[++i]);
		coarseBankIn->psi3Min = randIn->psi3Min = atof(argv[++i]);
	}
      else if ( strcmp(argv[i],	"--mm")		==0){
	     coarseBankIn->mmCoarse = atof(argv[++i]);
      }	  
      else if ( strcmp(argv[i],	"--number-fcut")	==0)  
  	  coarseBankIn->numFcutTemplates = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--simulation-type")	==0)	     	
	      randIn->type = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--signal-amplitude")	==0)	     
    	      randIn->SignalAmp = atof(argv[++i]);	
      else if ( strcmp(argv[i],	"--noise-amplitude")	==0)	     
    	      randIn->NoiseAmp = atof(argv[++i]);	
      else if ( strcmp(argv[i],	"--alpha-bank")	==0)	     
    	      coarseBankIn->alpha = atof(argv[++i]);      
      else if ( strcmp(argv[i],	"--alpha-signal")==0)	     
    	      randIn->param.alpha = atof(argv[++i]);
      else if ( strcmp(argv[i],	"--freq-moment-bank")	==0)	     
	coarseBankIn->fUpper = atof(argv[++i]);      
      else if ( strcmp(argv[i],	"--template-order")==0)     
	      coarseBankIn->order = atoi(argv[++i]); 	
      else if ( strcmp(argv[i],	"--signal-order")==0) 
    	      randIn->param.order = atoi(argv[++i]); 	
      else if ( strcmp(argv[i],	"--ntrial")	==0)
	      otherIn->ntrials = atoi(argv[++i]);       	
      else if ( strcmp(argv[i],	"--n")	==0)
	      otherIn->ntrials = atoi(argv[++i]);       	
      else if ( strcmp(argv[i],	"--seed")	==0)
	      randIn->useed = atoi(argv[++i])+1;	
      else if ( strcmp(argv[i],"--lengthFactor")		==0)   
	      otherIn->lengthFactor = atoi(argv[++i]);	
      else if ( strcmp(argv[i],	"--noise-model")==0){
	  i++;

	  if (strcmp(argv[i], "LIGOI")		== 0)	otherIn->NoiseModel = LIGOI;
	  else if (strcmp(argv[i], "LIGOA") 	== 0)  	otherIn->NoiseModel = LIGOA;
	  else if (strcmp(argv[i], "VIRGO") 	== 0)  	otherIn->NoiseModel = VIRGO;
	  else if (strcmp(argv[i], "TAMA") 	== 0)   otherIn->NoiseModel = TAMA;
	  else if (strcmp(argv[i], "GEO") 	== 0)   otherIn->NoiseModel = GEO;
	  else 
	    {
	      otherIn->NoiseModel = REALPSD;
	      otherIn->filename= argv[i];
	      /*TODO check whether the file exist otherwise quit with error message */
	    }	  
	}
      else if ( strcmp(argv[i],"--quiet")		==0)     otherIn->quietFlag 		= 1;
      else if ( strcmp(argv[i],"--InQuadrature")	==0)     otherIn->overlapMethod 	= InQuadrature;
      else if ( strcmp(argv[i],"--print-ambiguity-function")==0) otherIn->ambiguityFunction 	= 1;
      else if ( strcmp(argv[i],"--FMaximization")	==0) 	 otherIn->FMaximization 	= 1;
      else if ( strcmp(argv[i],"--print-overlap")	==0)  	 otherIn->PrintOverlap 		= 1;
      else if ( strcmp(argv[i],"--print-psd")	        ==0)  	 otherIn->PrintPsd 		= 1;
      else if ( strcmp(argv[i],"--PrintFilter")		==0) 	 otherIn->PrintFilter 		= 1; 
      else if ( strcmp(argv[i],"--print-bank-overlap")	==0)  	 otherIn->PrintBankOverlap 	= 1;
      else if ( strcmp(argv[i],"--print-template")	==0)  	 otherIn->PrintTemplate 	= 1;
      else if ( strcmp(argv[i],"--print-bank")		==0) 	 otherIn->PrintBank		= 1;
      else if ( strcmp(argv[i],"--print-bank-xml")	==0) 	 otherIn->PrintBankXml	        = 1;
      else if ( strcmp(argv[i],"--print-result-xml")	==0) 	 otherIn->PrintResultXml        = 1;
      else if ( strcmp(argv[i],"--print-prototype")	==0) 	 otherIn->PrintPrototype        = 1;
      else if ( strcmp(argv[i],"--read-bank-xml")       ==0){
	otherIn->inputXMLBank = argv[++i];
      }
      else if ( strcmp(argv[i],"--fastSimulation")      ==0)     otherIn->FastSimulation        = 1;
      else if ( strcmp(argv[i],"--binaryInjection")      ==0) {
	  if ( strcmp(argv[++i],	"BHNS")	==0)   {
	    otherIn->binaryInjection = BHNS;
	    /*should force randIn.MMax to be 3+ maxmass and not 2timesmaxmass. But
	     mass-range option might be given after that option. so we should do that later.
	    */
	  }
	  else {
	    fprintf(stderr, "binaryInjection should be follow by [BHNS]\n");
	    Help(*coarseBankIn, *randIn, *otherIn);
	    exit(0);
	  }
	  
      }


      else if ( strcmp(argv[i],"--check")		==0)     otherIn->check 		= 1;
      else if ( strcmp(argv[i],"--InputPSD")		==0)  	 otherIn->inputPSD 		= argv[++i];
      else if ( strcmp(argv[i],"--template")		==0)
	{
	  if ( strcmp(argv[++i],	"TaylorT1")	==0)	otherIn->template = TaylorT1;
	  else if ( strcmp(argv[i],	"TaylorT2")	==0)	otherIn->template = TaylorT2;
	  else if ( strcmp(argv[i],	"TaylorT3")	==0)	otherIn->template = TaylorT3;
	  else if ( strcmp(argv[i],	"TaylorF1")	==0)	otherIn->template = TaylorF1;
	  else if ( strcmp(argv[i],	"TaylorF2")	==0)	otherIn->template = TaylorF2;
	  else if ( strcmp(argv[i],	"PadeT1")	==0)	otherIn->template = PadeT1;
	  else if ( strcmp(argv[i],	"PadeF1")	==0)	otherIn->template = PadeF1;
	  else if ( strcmp(argv[i],	"EOB")		==0)	otherIn->template = EOB;
	  else if ( strcmp(argv[i],	"BCV")		==0)    otherIn->template = BCV;
	  else if ( strcmp(argv[i],	"SpinTaylorT3")	==0)	otherIn->template = SpinTaylorT3;

	  coarseBankIn->approximant = otherIn->template;
	  if ( coarseBankIn->approximant == BCV ) 	
	    coarseBankIn->space = Psi0Psi3; 
	  else 
	    coarseBankIn->space = Tau0Tau3;
	}	
      else if ( (strcmp(argv[i], "--signal")==0))
	{
	  if (strcmp(argv[++i],	  "TaylorT1")	==0) 		otherIn->signal = TaylorT1;
	  else if (strcmp(argv[i],"TaylorT2")	==0)		otherIn->signal = TaylorT2;
	  else if (strcmp(argv[i],"TaylorT3")	==0)	    	otherIn->signal = TaylorT3;
	  else if (strcmp(argv[i],"TaylorF1")	==0)	    	otherIn->signal = TaylorF1;
	  else if (strcmp(argv[i],"TaylorF2")	==0)	    	otherIn->signal = TaylorF2;
	  else if (strcmp(argv[i],"PadeT1")	==0)	    	otherIn->signal = PadeT1;
	  else if (strcmp(argv[i],"PadeF1")	==0)	    	otherIn->signal = PadeF1;
	  else if (strcmp(argv[i],"EOB")	==0)    	otherIn->signal = EOB;
	  else if (strcmp(argv[i],"BCV")	==0)            otherIn->signal = BCV;
	  else if (strcmp(argv[i],"SpinTaylorT3")==0)	    	otherIn->signal = SpinTaylorT3;
	  randIn->param.approximant = otherIn->signal;	
	  
	}	
      else 
	{
	  Help(*coarseBankIn, *randIn, *otherIn);
	  exit(0);
	}
      
      i++;       
    } 
}




/* --- Function to check validity of some parsed parameters (TODO)--- */
void CheckParams(InspiralCoarseBankIn coarseBankIn,
		 RandomInspiralSignalIn randIn,
		 OtherParamIn otherIn)
{
  REAL4 temp;
  temp =randIn.param.mass1; /*just to avoid boring warning*/

  if (coarseBankIn.psi0Min <=0 ) {
  	BEPrintError("Psi0 must be > 0");
	exit(0);
  }
  if (coarseBankIn.psi0Max <=0 ) {
	  BEPrintError("psi0 Max must be > 0"); 
	  exit(0);
  }
  if (coarseBankIn.psi3Max >=0 ) {
	  BEPrintError("psi3 Max must be < 0"); 
	  exit(0);
  }
  if (coarseBankIn.psi3Min >=0 ) {
	BEPrintError("psi3 Min must be < 0"); 
	exit(0);
  }
  if (coarseBankIn.psi0Min > coarseBankIn.psi0Max
      || coarseBankIn.psi3Min > coarseBankIn.psi3Max){
	BEPrintError("psi range should be [psiMin psiMax]; take care of the sign. (i.e. -10 -2000 or 10 2000)\n");
      	exit(0);
    }
  if (coarseBankIn.mMin > coarseBankIn.mMax )
  {
	  BEPrintError("in option --mass-range, first argument shoukd be < to the second one\n");
  }
  if (coarseBankIn.approximant != BCV && otherIn.overlapMethod != InQuadrature)
    {
      BEPrintError("If the template are not BCV, the overlap \nmethod must be InQuadrature (use the \"--InQuadrature\" option)\n");
      exit(0);
    }
  if (coarseBankIn.tSampling <= 2.*coarseBankIn.fUpper)
    {
      BEPrintError("arg given by  --freq-moment-bank should be < to half the value of the argument given by --sampling \n");
      exit(0);
    }
  if (otherIn.binaryInjection==BHNS){
    if (randIn.mMin >3 || randIn.mMax < 3 ){
      BEPrintError("if you want to inject BHNS systems then adjust the mass-range so that the minimum is less than 3 solar mass and the maximum  is greater than 3solar mass !! \n");
      exit(0);
    }

  }
    

}

void
BEPrintError(char *chaine)
{
  fprintf(stderr,"|=============================================================|\n");
  fprintf(stderr,"| BankEfficiency code						\n");
  fprintf(stderr,"| Error while parsing parameters				\n");
  fprintf(stderr,"|=============================================================|\n");
  fprintf(stderr,"| %s \n",chaine);
  fprintf(stderr,"|type BankEfficiency -h to get details\n");
  fprintf(stderr,"|=============================================================|\n");
  
}



/* ****************************************************************************
 *  Documenation  on line
 *  **************************************************************************/
void Help(	InspiralCoarseBankIn   coarseBankIn,
       		RandomInspiralSignalIn randIn,
		OtherParamIn           otherIn)
{
  INT4 		temp;
	
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"BankEfficiency Code\n");
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"\nPURPOSE: Test efficiency of bank of templates using Monte-Carlo simulations\n");

  fprintf(stderr,"\nUSAGE:  [options]\n");
  fprintf(stderr,"There are two categories of options (default values in brackets)\n");
  fprintf(stderr,"First options must be followed by one or two arguments:\n\n");
  fprintf(stderr,"--alpha-bank		: BCV amplitude correction parameter 			(%7.2f)\n", BANKEFFICIENCY_ALPHABANK);
  fprintf(stderr,"--alpha-signal	: BCV amplitude correction parameter 			(%7.2f)\n",BANKEFFICIENCY_ALPHASIGNAL);
  fprintf(stderr,"--binaryInjection[BHNS]: if argument=BHNS then system are Blackhole+neutron star\n");
  fprintf(stderr,"--fl			: lower frequency cutoff             			(%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"--length-factor	: multiply length signal by some factor(should be of unity order) ( %7.2d)\n", 1);
  fprintf(stderr,"--fend-bcv		: Lower and highest values of bcv frequency cutoff 	(%7.2d, %7.2d)\n", BANKEFFICIENCY_LOWGM, BANKEFFICIENCY_HIGHGM);
  fprintf(stderr,"--mass-range		: minimal mass of component stars    			(%7.2f, %7.2f) Mo\n", BANKEFFICIENCY_MMIN, BANKEFFICIENCY_MMAX);
  fprintf(stderr,"--mm			: minimal match for template bank    			(%7.3f)\n", BANKEFFICIENCY_MMCOARSE);
  fprintf(stderr,"--noise-amplitude	: use it with simulation-type <> 0 			(%7.2f)\n", BANKEFFICIENCY_NOISEAMPLITUDE); 
  fprintf(stderr,"--noise-model		: design sensitivity to use(LIGOI ,LIGOA VIRGO, GEO, TAMA	(%s)\n", BANKEFFICIENCY_NOISEMODEL);
  fprintf(stderr,"--number-fcut		: number of layers in Fcut dimension 			(%7.2d)\n", BANKEFFICIENCY_NFCUT);
  fprintf(stderr,"--ntrial		: number of trials                   			(%7.2d)\n", BANKEFFICIENCY_NTRIALS);
  fprintf(stderr,"--psi0-range		: psi0 range in the bank				(%7.2f, %7.2f) \n", BANKEFFICIENCY_PSI0MIN, BANKEFFICIENCY_PSI0MAX);
  fprintf(stderr,"--psi3-range		: psi3 range in the bank 				(%7.2f, %7.2f) \n", BANKEFFICIENCY_PSI3MAX, BANKEFFICIENCY_PSI3MIN);
  fprintf(stderr,"--template		: PN model for template bank, e.g. BCV 			(%7.2d)\n", BANKEFFICIENCY_TEMPLATE);
  fprintf(stderr,"--sampling		: sampling frequency             			(%7.2f)\n", BANKEFFICIENCY_TSAMPLING);
  fprintf(stderr,"--seed			: seed for random generation         		(%7.2d)\n", BANKEFFICIENCY_USEED);
  fprintf(stderr,"--signal		: model to inject in Monte-Carlo, e.g. EOB            	(%7.2d)\n", BANKEFFICIENCY_SIGNAL);
  fprintf(stderr,"--signal-order		: order of PN model                             (%7.2d)\n", BANKEFFICIENCY_ORDER_SIGNAL);
  fprintf(stderr,"--simulation-type	: type of simulation, 0, 1 or 2      			(%7.2d)\n", BANKEFFICIENCY_TYPE);
  fprintf(stderr,"--template-order	: order of signal to injec                            	(%7.2d)\n", BANKEFFICIENCY_ORDER_TEMPLATE);
  fprintf(stderr,"--signal-amplitude	: amplitude of the signal                             	(%7.2f)\n", BANKEFFICIENCY_SIGNALAMPLITUDE);

  fprintf(stderr,"\n --\nThe second category doesn't take any arguments.\n Those options are mainly used to print results in files\n");
  fprintf(stderr,"Verbose Options \n\n");
  fprintf(stderr,"     --check		: For a given random waveform, we compute overlap with same parameter (should get overlap=1)\n");
  /*fprintf(stderr,"     --FMaximization	: Maximize over Ending frequency 			(%d)\n", BANKEFFICIENCY_FMAXIMIZATION);*/
  fprintf(stderr,"     --FastSimulation		: filters only template around the injection (dt0=.5 and dt3=.25)\n");

  fprintf(stderr,"     --PrintOverlap	: Print Overlap given by the best template 		(%d)\n", BANKEFFICIENCY_PRINTOVERLAP);
  /*  fprintf(stderr,"     --PrintFilter	: Print filters giving the best overlap template 	(%d)\n", BANKEFFICIENCY_PRINTFILTER);
  fprintf(stderr,"     --print-ambiguity-function : print ambiguity function (bank overlap)	(%d)\n", BANKEFFICIENCY_AMBIGUITYFUNCTION);*/

  fprintf(stderr,"     --print-bank	: Print the bank in ascii format (%s)		\n",  BANKEFFICIENCY_PRINTBANK_FILEASCII);
  fprintf(stderr,"     --print-bank-xml	: Print the bank in xml format (%s)		\n",  BANKEFFICIENCY_PRINTBANK_FILEXML);
  fprintf(stderr,"    --print-result-xml  : Print result in xml format (%s)		\n",  BANKEFFICIENCY_PRINTRESULT_FILEXML);


  fprintf(stderr,"     --print-overlap	: Print the overlap given by the best template	(%s)	(%d)\n", 
		  BANKEFFICIENCY_PRINTOVERLAP_FILE, 
		  BANKEFFICIENCY_PRINTOVERLAP);
  fprintf(stderr,"     --print-psd	: Print the psd in                              	(%s)	(%d)\n", 
		  BANKEFFICIENCY_PRINTPSD_FILE, 
		  BANKEFFICIENCY_PRINTPSD);

  fprintf(stderr,"     --InQuadrature	: Overlap is compute without alpha maximization\n");
  fprintf(stderr,"     --quiet		: if this flag is present, the output is restricted to the minimum\n");
/*  fprintf(stderr,"     --InputPSD	: to give the name of real PSD\n");*/
  /*to avoid boring warning*/
  temp = coarseBankIn.approximant;
  temp =  randIn.useed;
  temp =  otherIn.signal;
}


		      
/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */
void
KeepHighestValues(InspiralWaveOverlapOut in , 
		  INT4 j, INT4 l, double frequency,
		  InspiralWaveOverlapOut *out, 
		  INT4 *jmax, INT4 *lmax, double *fMax)
{
  if (in.max > out->max)
    { 
      out->max   = in.max;
      out->phase = in.phase;
      out->alpha = in.alpha;
      out->bin   = in.bin;
      *jmax      = j;
      *lmax 	 = l;
      *fMax 	 = frequency;
    }
}


/* ****************************************************************************
 * Some output Results
 * ************************************************************************** */
void
PrintResults(	    	InspiralTemplate       	bank,
	    		InspiralTemplate       	injected,
	    		InspiralWaveOverlapOut 	overlapout,
	    		REAL8 			fendBCV,
	    		UINT4 			layer)
{
 if (bank.approximant==BCV){  
  	fprintf(stdout, "%e %e ", bank.psi0, bank.psi3);  /*triggered */
	fprintf(stdout, "%e %e    ", injected.psi0, injected.psi3);
      	fprintf(stdout, "%e %e %e %e   ", fendBCV,   injected.fFinal, bank.totalMass, bank.eta);
      	fprintf(stdout, "%e %e %e %e   ", injected.mass1, injected.mass2, overlapout.max, overlapout.phase);
      	fprintf(stdout, "%e %e %d %d\n ", overlapout.alpha, overlapout.alpha*pow(fendBCV,2./3.), layer, overlapout.bin);
	}
  else{
  	fprintf(stdout, "%e %e ", bank.t0, bank.t3);  /*triggered */
	fprintf(stdout, "%e %e    ", injected.t0, injected.t3);
	fprintf(stdout, "%e %e %e %e   ", fendBCV,   injected.fFinal, bank.totalMass, bank.eta);
	fprintf(stdout, "%e %e %e %e   ", injected.mass1, injected.mass2, overlapout.max, overlapout.phase);
	fprintf(stdout, "%e %e %d %d\n ", overlapout.alpha, overlapout.alpha*pow(fendBCV,2./3.), layer, overlapout.bin);
	}
}





/* ****************************************************************************
 * Functio to generate and stored the moments in  REAL4vectors. 
 * ************************************************************************* */
void LALCreateMomentVector(
			LALStatus             *status,
			REAL4Vector           *a11,
			REAL4Vector           *a21,
			REAL4Vector           *a22,
			REAL8FrequencySeries  *psd,
			InspiralTemplate      *params)
{
  REAL8 		m7 = 0;							/* the moments */
  REAL8 		m5 = 0;	
  REAL8 		m3 = 0;
  REAL8 		f;
	  
  INT4 			kMin;
  INT4		  	kMax;
  INT4			k;

  InspiralMomentsIn 	in;
    
  INITSTATUS (status, "LALInspiralGetBCVMaximizationMatrix", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  /* inspiral structure */
  in.shf 	= psd;  							/* The spectrum 			*/
  in.xmin 	= params->fLower;						/* Lower frequency of integral?	is it correct or should we put zero? 	*/
  in.xmax 	= params->tSampling/2.;						/* We dont' need to cut here		*/
  in.norm 	= 1./4.;							/* Some normalization			*/


  in.norm*=params->tSampling * params->tSampling;

  /* The minimum and maximum frequency */
  kMin = (UINT8) floor( (in.xmin - psd->f0) / psd->deltaF );
  kMax = (UINT8) floor( (in.xmax - psd->f0) / psd->deltaF );
  
  /* Skip frequency less than fLower.*/
  for (k=0; k<kMin ; k++)
    {
      a11->data[k] = 0;
      a21->data[k] = 0;
      a22->data[k] = 0;
    }
  for ( k = kMin; k < kMax; ++k )
    {
      f = psd->f0 + (REAL8) k * psd->deltaF;

      if (psd->data->data[k])
	{
	  m7 +=  pow( f, -(7./3.) ) / psd->data->data[k] * psd->deltaF / in.norm;
	  m5 +=  pow( f, -(5./3) )  / psd->data->data[k] * psd->deltaF / in.norm;
	  m3 +=  pow( f, -(1.) )    / psd->data->data[k] * psd->deltaF / in.norm;

	  a11->data[k] = 1./sqrt(m7);  
	  a22->data[k] = 1. / sqrt(m3 - m5*m5/m7);
	  a21->data[k] = -m5 / m7 * a22->data[k];
	}
    }

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

  

  

/* ****************************************************************************
 * Create an orthogonal filter. by taking its complex conjuguate
 * ************************************************************************** */
void LALGetOrthogonalFilter2(	
			REAL4Vector *filter
			)
{
  UINT4   		i;
  UINT4			n 	= filter->length;
  UINT4 		nby2    = filter->length / 2;
	  
  REAL4	  temp;
  
  for (i = 1; i < nby2; i++)
    {
      temp 		=  filter->data[i];
      filter->data[i] 	= -filter->data[n-i];
      filter->data[n-i] = temp;
    }
}

  

/* *****************************************************************************
 * Create a vector f^{a/b} 
 * ************************************************************************** */
void 
LALCreateVectorFreqPower(
       			REAL4Vector 		*vector,
			InspiralTemplate 	params,
			INT4 			a,
			INT4 			b)
{
  INT4 			i;
  INT4			n = vector->length;						/* Length of the vector to create 	*/
	  
  REAL8		power = (REAL8)a / (REAL8)b;					/* the frequency power			*/
  REAL8 	f;								/* La frequence				*/
  REAL8		df = params.tSampling/(REAL8)n/2.;				/* sampling frequency			*/

  /* First value equal to zero */
  vector->data[0] = 0.;
  for( i = 1; i < n; i++)
    {
      f = i * df; 
      /* Probably not efficient but this function is 
       * just called once at the beginning */
      vector->data[i] = pow(f, power);
    }
}


/* *****************************************************************************
 * Create the filters for BCV correlation 					
 * ************************************************************************** */
void 
LALCreateFilters(	REAL4Vector 		*Filter1,
		 	REAL4Vector 		*Filter2,
		 	REAL4Vector 		VectorPowerFm5_3,
		 	REAL4Vector 		VectorPowerFm2_3,
		 	REAL4Vector 		VectorPowerFm7_6,
		 	REAL4Vector 		VectorPowerFm1_2,
		 	BCVMaximizationMatrix   matrix,				/* moments 				*/
		 	UINT4 			kMin,				/* position of the integration		*/
		 	double 			psi0,				/* change to real4 or real 8		*/
		 	double 			psi3)				/* idem 				*/
{
  UINT4			i;
  UINT4  		n = Filter1->length;
  UINT4	  		nby2 = Filter1->length / 2;
	  
  REAL8 		amplitude;
  REAL8	  		cos_phase;
  REAL8	  		sin_phase;
  
  /* Create the templates */	       
  for (i = 0; i< kMin-1; i++)
    {
      Filter1->data[i] = 0;
      Filter2->data[i] = 0;
    }
  /* Should add the upper fendBCV here as well */
  for (i=kMin; i< nby2; i++)
    {
      amplitude  =   psi0 * VectorPowerFm5_3.data[i]  				/* The same for the two filters 	*/
	      + psi3 * VectorPowerFm2_3.data[i];
      
      cos_phase  = cos(amplitude);						/* The same for the two filters		*/
      sin_phase  = sin(amplitude);						/* The same for the two filters 	*/							          
      /* Fill the first filter here */
      amplitude  = matrix.a11 * VectorPowerFm7_6.data[i];
      
      Filter1->data[i]   =  amplitude * cos_phase;      
      Filter1->data[n-i] = -amplitude * sin_phase;
      /* Fill the second one here */
      amplitude =  matrix.a21 * VectorPowerFm7_6.data[i] + 
	matrix.a22 * VectorPowerFm1_2.data[i];
      
      Filter2->data[i]   =  amplitude * cos_phase;      
      Filter2->data[n-i] = -amplitude * sin_phase;
    }  
}


/* ****************************************************************************
 * Draft functions to finalize 
 * ************************************************************************* */
void LALGenerateWaveform(LALStatus              *status,
			 REAL4Vector            *signal,
			 RandomInspiralSignalIn *randIn)


{
  REAL8  		norm;
  
  REAL4Vector  		buff;
  
  InspiralWaveNormaliseIn normin;
  
  INITSTATUS (status, "LALGenerateWaveform", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  buff.length = signal->length;
  if (!(buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length))) {
    ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
  }
  
  /* set up the structure for normalising the signal */
  normin.psd = &(randIn->psd);
  normin.df = randIn->param.tSampling / (REAL8) signal->length;
  normin.fCutoff = randIn->param.fCutoff;
  normin.samplingRate = randIn->param.tSampling;


  LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
  
  if (randIn->param.approximant == TaylorF1 ||
      randIn->param.approximant == TaylorF2 ||
      randIn->param.approximant == PadeF1)
    {
      LALInspiralWave(status->statusPtr, signal, &randIn->param);
      CHECKSTATUSPTR(status);
    }
  else
    {
      randIn->param.fFinal=0;
      LALInspiralWave(status->statusPtr, &buff, &randIn->param);
      CHECKSTATUSPTR(status);
      LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
      CHECKSTATUSPTR(status);
    }
  LALInspiralWaveNormaliseLSO(status->statusPtr, signal, &norm, &normin);
  CHECKSTATUSPTR(status);
  
  
  if (buff.data != NULL) LALFree(buff.data);
  DETATCHSTATUSPTR(status);
   RETURN(status);
}






/*  ****************************************************************************
 *  The Overlap function for BCV templates 
 *  ***************************************************************************/
void
LALWaveOverlapBCV(	     LALStatus               *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapOut  *overlapout,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector             *Filter1,
			     REAL4Vector             *Filter2,
			     BCVMaximizationMatrix    matrix,
			     OtherParamIn             otherIn 
		  	     )
     /*  </lalVerbatim>  */
{
  REAL4 
	  rhoMax 	= 0,							/* SNR maximum 				*/ 
	  maxV1 	= 0, maxV2 	= 0, 					/* V1 and V2 value  of triggered bin	*/
	  x1_2, x2_2, x3_2, x4_2, 						/* some temporary data 	      		*/
	  V0, V1, V2, 								/* idem 				*/
	  rho		= 0;							/* the SNR 				*/
	  
  INT4 
	  maxbin 	= 0;							/* position of the maximum SNR 		*/
  
  REAL4Vector 
	  x1, x2, x3, x4;							/* output of the correlations 		*/

  InspiralWaveCorrelateIn 
	  corrin;								/* the correlation input structure 	*/
  UINT4 
	  i, 									
	  nBegin, 								/* beginning of valid correlation 	*/
	  nEnd, 								/* end of valid correlation 		*/
	  n = correlation->length; 						/* length of the vectors		*/
  REAL4 phi, phase, cphase, sphase, cphi, sphi;
	  
  INITSTATUS (status, "LALInspiralWaveOverlapBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  /* size affectation*/
  x1.length = x2.length = x3.length = x4.length = n;
  
  /* memory allocation*/
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);

  /**/
  overlapin->param.nStartPad 	= 0;
  overlapin->param.startPhase 	= 0;
  
  /*  We want to compute <x|h1> and <x|h2>; let's prepare the correlation  
   *  process; let's get the input data namely x and the spectrum namely psd 
   *  */
  corrin.fCutoff        = overlapin->param.fFinal;				/* the integral's cutoff frequency 	*/
  corrin.samplingRate   = overlapin->param.tSampling;				/* the sampling 			*/
  corrin.df             = overlapin->param.tSampling / 
	  		(REAL8) overlapin->signal.length;			/* the resolution in frequency 		*/
  corrin.psd            = overlapin->psd;					/* the spectrum 			*/
  corrin.signal1        = overlapin->signal;					/* the signal x				*/
  corrin.revp           = overlapin->revp;					/* fourier transform flag		*/
  
  /* Filter h1 and its orthonormal h1* */  
  corrin.signal2        = *Filter1;						/* The first template 			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);			/* the correlation			*/
  CHECKSTATUSPTR(status);							
  LALGetOrthogonalFilter2(Filter1);						/* get the orthonormalized template     */
  corrin.signal2        = *Filter1;						/* the second template 		*/
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);
  
  /* Filter h2 and its orthornormal filter */
  corrin.signal2        = *Filter2;						/* The first template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);			/* the correlation 			*/
  CHECKSTATUSPTR(status);		
  LALGetOrthogonalFilter2( Filter2);						/* get the orthonormailzed templates 	*/
  corrin.signal2        = *Filter2;						/* the second template			*/
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);			/* the correlation 			*/
  
  /* nbegin and end for the correlation processus 
   * as well as overlapout affectation */
  nBegin                = overlapin->nBegin;					/* starting bin 			*/
  nEnd              	= Filter1->length - overlapin->nEnd;			/* ending valid bin 			*/
  overlapout->max   	= 0.;							
  overlapout->bin   	= nBegin;
  overlapout->phase 	= 0;

  /* some output for debugging, checking ... */
#if 0
  if (otherIn.PrintOverlap)
    {
      for (i=0; i< x1.length; i++)
	printf("%d %e\n", i,fabs(x1.data[i]));
      printf("&\n");	       
      
      for (i=1; i< Filter1->length; i++)
	printf("%d %e\n", i,fabs(x2.data[i]));
      printf("&\n");	       
      
      for (i=0; i< Filter1->length; i++)
	printf("%d %e\n", i,fabs(x3.data[i]));
      printf("&\n");	       
      
      for (i=0; i< Filter1->length-1; i++)
	printf("%d %e\n", i,fabs(x4.data[i]));
      printf("&\n");	       
    }
#endif


  /* Once we have the 4 correlations, we can search for the
   * maximum of correlation and compute the final SNR */
  for (i = nBegin; i < nEnd; i++) 
    {
      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      
      V0 = x1_2 + x2_2 + x3_2 + x4_2;
      V1 = x1_2 + x3_2 - x2_2 - x4_2;
      V2 = 2*x1.data[i]*x2.data[i] + 2*x3.data[i]*x4.data[i];

      rho = V0 + sqrt(V1*V1+V2*V2);						/* final correlation			*/
      
     if(otherIn.PrintOverlap==1)						/* we don't need to store all the cor- 	*/
	correlation->data[i] = sqrt(rho /2.);					/* relation except if needed		*/
      if ( rho > rhoMax)							/* Keep the position of the max only	*/
	{
	  rhoMax 	= rho;
	  maxbin 	= i;
	  maxV1  	= V1;
	  maxV2 	= V2;
	}          
    }


  /* Let's compute the phase corresponding to the maximum 
   * and store some results.*/ 
  overlapout->bin       = maxbin;						/* position of the maximum 		*/
  overlapout->phase 	= .5 * atan2(maxV2, maxV1);				/* its phase				*/
  overlapout->max 	= sqrt(rhoMax /2.) ; 					/* its maximum				*/

  /* Is it needed ?*/
  while (overlapout->phase > LAL_PI || overlapout->phase > 0 ) 
	  overlapout->phase -= LAL_PI;
  while (overlapout->phase < -LAL_PI || overlapout->phase < 0 ) 
	  overlapout->phase += LAL_PI;

  overlapout->alpha = -(matrix.a22 * tan(overlapout->phase)) 				/* Compute the final alpha parameter 	*/
	  / (matrix.a11 + matrix.a21* tan(overlapout->phase));				/* which optimize the overlap		*/	
  /* The final template */
  #if 0
   for (i=0; i< correlation->length; i++)
	printf("%d %e\n", i,fabs(correlation->data[i]));
   printf("&\n");	 
	phase 		= overlapout->phase;
	  phi   		= overlapin->param.startPhase;
   
    for (i=0; i<(UINT4)correlation->length; i++) 
    {
      cphase = cos(phase);
      cphi   = cos(phi);
      sphase = sqrt(1-cphase*cphase);
      sphi   = sqrt(1-cphi*cphi);
      
      correlation->data[i] = x1.data[i] * cphase * cphi
	+ x2.data[i] * sphase * cphi
	+ x3.data[i] * cphase * sphi
	+ x4.data[i] * sphase * sphi;	
    }
  for (i=0; i< correlation->length; i++)
	printf("%d %e\n", i,fabs(correlation->data[i]));
      printf("&\n");	    
  #endif
  
  /* Free memory */
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}



/* *****************************************************************************
 *  A draft function to finalized later ; purpose: print the bank 
 *  ****************************************************************************/
void
BEPrintBank(	InspiralCoarseBankIn coarseBankIn,
		InspiralTemplateList **list,
      		UINT4 nlist
		) 
{
  UINT4 	i;

  FILE 		*output;


  output = fopen(BANKEFFICIENCY_PRINTBANK_FILEASCII,"w");

  fprintf(output,"#Number of Coarse Bank Templates=%d\n",nlist);
  if (coarseBankIn.approximant == BCV)
    {
      fprintf(output, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n", 
	      coarseBankIn.psi0Min, coarseBankIn.psi0Max, coarseBankIn.psi3Min, coarseBankIn.psi3Max);
      fprintf(output, "#psi0 psi3 nLayer totalMass fFinal\n");
      
    }
  else
    {
      fprintf(output, "#mMin=%e, mMax=%e\n", coarseBankIn.mMin,coarseBankIn.mMax);
      fprintf(output, "#tau0, tau3, mass1, mass2\n");
    }

  for ( i = 0; i < nlist; i++)
    {
      if (coarseBankIn.approximant == BCV)
	{      		   
      		fprintf(output, "%e %e %d %e %e\n", 	(*list)[i].params.psi0, 
				  			(*list)[i].params.psi3,
				  			(*list)[i].nLayer, 
							(*list)[i].params.totalMass,
						       	(*list)[i].params.fFinal);
	} 
      else
	{

      		fprintf(output, "%e %e %e %e\n", 	(*list)[i].params.t0, 
				  			(*list)[i].params.t3, 
							(*list)[i].params.mass1, 
							(*list)[i].params.mass2);
	}
    }	   
  fprintf(output,"&\n");
  fclose(output);
}

/* personal function. Should not be used. output hte bank  overlaps in a binary
 * file to be used with DES3.
 * */

void
PrintBankOverlap(
		 InspiralTemplateList 	**list,
		 int 			nlist,
		 float 			*overlap,
		 InspiralCoarseBankIn 	coarseBankIn
		)
{

}
#if 0
FILE 		*output1;
  FILE 		*output2;
  
  long Npsi0, Npsi3;

  double dx0, dx3;
  double psi0Min = coarseBankIn.psi0Min;
  double psi0Max = coarseBankIn.psi0Max;
  double psi3Min = coarseBankIn.psi3Min;
  double psi3Max = coarseBankIn.psi3Max;
  double psi0, psi3;

  int    numfcut = coarseBankIn.numFcutTemplates;
  int    i,j,l,n;

  float  *a;

  double minimalMatch = coarseBankIn.mmCoarse;
  double theta, myphi, fac;


  output1 = fopen("FF.sr4", "a+");
  output2 = fopen("FF.dim", "w");
   dx0 = sqrt(2.L * (1.L - minimalMatch)/(*list)[0].metric.g00 );
   dx3 = sqrt(2.L * (1.L - minimalMatch)/(*list)[0].metric.g11 );


   fprintf(stderr, "%f %f\n", dx0, dx3);
   
     if ((*list)[0].metric.theta!=0.L)
   {
     myphi = atan2(dx3, dx0);
     theta = fabs((*list)[0].metric.theta);
     if (theta <= myphi) {
       fac = cos(theta);
       dx0  /= fac;
       dx3 *=  fac;
     } 
     else {
       fac = sin(theta);
       dx0 = dx3 / fac;
       dx3 = dx0 * fac;
     }
   }
   
     fprintf(stderr, "%f %f\n", dx0, dx3);
     fprintf(stderr,"Approxi= %d\n", coarseBankIn.approximant); 
     a =( float*)malloc(sizeof(float) ); 

     switch( coarseBankIn.approximant ){
	case BCV: 
       /*some data*/
       
       
       Npsi0 = floor(( psi0Max - psi0Min ) / dx0) + 1 ;
       Npsi3 = floor(-( psi3Min - psi3Max ) / dx3) + 1;

       printf("%f %f %ld %ld\n", dx0, dx3, Npsi0, Npsi3);     
       
/*	Npsi0 = 56; */
	Npsi3 = nlist/Npsi0/5;  
	dx0  =  ( psi0Max - psi0Min)/(Npsi0);
	dx3  =  ( psi3Max - psi3Min)/(Npsi3);
	
       printf("%f %f %ld %ld\n", dx0, dx3, Npsi0, Npsi3);     
       fflush(stdout); 
       /* The dimension file */
       fprintf(output2,"%7ld %14.8f %14.8f  %s\n",Npsi3, -psi3Max+dx3, dx3,"psi3");
       fprintf(output2,"%7ld %14.8f %14.8f  %s\n",Npsi0, psi0Min+dx0, dx0,"psi0");
       fprintf(output2,"%7d %14.8f %14.8f  %s\n",numfcut,1.,1.,"layers");
       
       fclose(output2);
       float tab[Npsi0][Npsi3][numfcut];
       for (l = 0; l < numfcut; l++){
	 for (i = 0; i < Npsi0; i++){
	   for (j = 0; j < Npsi3; j++){
	     tab[i][j][l]=0;
	   }
	 }
       }
       n = 0;
       /*fill to zero */
       while (n< nlist)
	 {
	   psi0 = (*list)[n].params.psi0;
	   psi3 = (*list)[n].params.psi3;
	   
	   i = (int)(( psi0 - psi0Min ) / dx0);
           j = (int)(-( psi3 - psi3Max) / dx3);;
	 if (i > Npsi0 || (long)j > Npsi3 || (*list)[n].nLayer > (unsigned int)numfcut)
		fprintf(stderr," %d %d %d %f %f \n", n, i, j, psi0, psi3); 
	  
  	  tab[i][j][(*list)[n].nLayer-1] = overlap[n];
	  n++;
	}
/*fprintf(stderr, "ok %d %d %d\n",numfcut, Npsi0, Npsi3);*/
      for (l = 1; l <= numfcut; l++){
       for (i = 0; i < Npsi0; i++){
	for (j = 0; j < Npsi3; j++){
/*        fprintf(stderr,"%ld %ld %ld %lf\n", i , j , l , tab[i][j][l]); */
       if (tab[i][j][l-1]!=0)
	a[0] = tab[i][j][l-1];	
	else a[0] = 0;
/*fprintf(stderr,"#%d %d %d %lf\n",i, j, l-1, a[0]);	*/
	fwrite(a,4,1,output1);
	 }
	}	 	
      }
      
    fclose(output1);
      break;
    case BCVSpin: 
    case SpinTaylorT3: 
    case TaylorT1: 
    case TaylorT2: 
    case TaylorT3: 
    case TaylorF1: 
    case TaylorF2: 
    case PadeT1: 
    case PadeF1: 
    case EOB:
    case SpinTaylor: 
      break;      
  }

}
#endif





/* routine to print the bank in a xml file*/
/* should add the status variable since we call sone lal functions here 
 * this routine is not well commented and won't be in close future.
 * */
/* !! tau2 is used to store the Layer number */
void 
BEPrintBankXml(
	       InspiralTemplateList *coarseList, 
	       UINT4 numCoarse,
	       InspiralCoarseBankIn   coarseBankIn,
	       RandomInspiralSignalIn randIn,
	       OtherParamIn           otherIn)
{
#define MAXIFO 2
  LALStatus             status = blank_status;
  MetadataTable         proctable;

  MetadataTable         templateBank;
  SnglInspiralTable     *tmplt  = NULL;
  CHAR  ifo[3];                           /* two character ifo code       */
  CHAR *channelName = NULL;               /* channel string               */
  UINT4 	i;  
  LIGOLwXMLStream       xmlStream;
  CHAR  fname[256];
  LIGOTimeGPS gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  
  strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
  memset( ifo, 0, sizeof(ifo) );
  memcpy( ifo, "MC", sizeof(ifo) - 1 );


  /* --- first we create the filename --- */
  LALSnprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTBANK_FILEXML ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  
  /* --- let's save the bank itself --- */
  templateBank.snglInspiralTable = NULL;

  

  /* --- convert the templates to sngl_inspiral structures 
     in order to  write results in the XML file --- */
  if ( numCoarse )
    {
      /* --- we sae one line --- */
      tmplt = templateBank.snglInspiralTable = (SnglInspiralTable *)
	LALCalloc( 1, sizeof(SnglInspiralTable) );
      LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
      LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		   "BankEfficiency" );
      LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		   channelName );
      tmplt->mass1   = (REAL4) coarseList[0].params.mass1;
      tmplt->mass2   = (REAL4) coarseList[0].params.mass2;
      tmplt->mchirp  = (REAL4) coarseList[0].params.chirpMass;
      tmplt->eta     = (REAL4) coarseList[0].params.eta;
      tmplt->tau0    = (REAL4) coarseList[0].params.t0;
      tmplt->tau2    = (REAL4) coarseList[0].nLayer;
      tmplt->tau3    = (REAL4) coarseList[0].params.t3;
      tmplt->tau4    = (REAL4) coarseList[0].params.t4;
      tmplt->tau5    = (REAL4) coarseList[0].params.t5;
      tmplt->ttotal  = (REAL4) coarseList[0].params.tC;
      tmplt->psi0    = (REAL4) coarseList[0].params.psi0;
      tmplt->psi3    = (REAL4) coarseList[0].params.psi3;
      tmplt->f_final = (REAL4) coarseList[0].params.fFinal;
      
      /* --- then the others --- */
      for ( i = 1; i < numCoarse; ++i )
	{
	  tmplt = tmplt->next = (SnglInspiralTable *)
	    LALCalloc( 1, sizeof(SnglInspiralTable) );
	  LALSnprintf( tmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), ifo );
	  LALSnprintf( tmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR), 
		       "BankEfficiency" );
	  LALSnprintf( tmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
		       channelName );
	  tmplt->mass1   = (REAL4) coarseList[i].params.mass1;
	  tmplt->mass2   = (REAL4) coarseList[i].params.mass2;
	  tmplt->mchirp  = (REAL4) coarseList[i].params.chirpMass;
	  tmplt->eta     = (REAL4) coarseList[i].params.eta;
	  tmplt->tau0    = (REAL4) coarseList[i].params.t0;
	  tmplt->tau2    = (REAL4) coarseList[i].nLayer;
	  tmplt->tau3    = (REAL4) coarseList[i].params.t3;
	  tmplt->tau4    = (REAL4) coarseList[i].params.t4;
	  tmplt->tau5    = (REAL4) coarseList[i].params.t5;
	  tmplt->ttotal  = (REAL4) coarseList[i].params.tC;
	  tmplt->psi0    = (REAL4) coarseList[i].params.psi0;
	  tmplt->psi3    = (REAL4) coarseList[i].params.psi3;
	  tmplt->f_final = (REAL4) coarseList[i].params.fFinal;
	}
    }
  
  /* -- we start to fill the xml file here --- */
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );
   

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
			    &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
				    PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);




  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );



  /* write process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
	       ifoName[0], ifoName[1] );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
			    &accuracy ), &status );
  
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	    &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
				    process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  
  /* write process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				    process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				    process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  
  /* finally write sngl inspiral table */
  if ( templateBank.snglInspiralTable )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sngl_inspiral_table ), 
		&status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
					sngl_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    }


  /* --- desallocate memory here --- */
  while ( templateBank.snglInspiralTable )
    {
      tmplt = templateBank.snglInspiralTable;
      templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
      LALFree( tmplt );
    }
  
  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlStream ), &status );
}



/* function to get matrix components given 
 * the vectors of the different moments   */
void BEGetMatrixFromVectors( 
			    REAL4Vector A11,
			    REAL4Vector A21,
			    REAL4Vector A22,
			    UINT4 k,
			    BCVMaximizationMatrix *matrix2fill
			    )
{
  matrix2fill->a11 = A11.data[k];
  matrix2fill->a21 = A21.data[k];
  matrix2fill->a22 = A22.data[k];   
}






void 
BEFillProc(
	       ProcessParamsTable     *this_proc_param,
	       InspiralCoarseBankIn   coarseBankIn,
	       RandomInspiralSignalIn randIn,
	       OtherParamIn           otherIn)
{


#define ADD_PROCESS_PARAM( pptype, format, ppname, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME ); \
  LALSnprintf( this_proc_param->param,   LIGOMETA_PARAM_MAX,   "%-20s", ppname); \
  LALSnprintf( this_proc_param->type,    LIGOMETA_TYPE_MAX,    "%-10s",   pptype ); \
  LALSnprintf( this_proc_param->value,   LIGOMETA_VALUE_MAX,   format, ppvalue );

  ADD_PROCESS_PARAM("float",	"%12.5f",	 "--alpha-bank",	coarseBankIn.alpha);
  ADD_PROCESS_PARAM("float",	"%12.5f",	 "--alpha-signal",	randIn.param.alpha);
  ADD_PROCESS_PARAM("float",	"%12.5f",	 "--fl",		coarseBankIn.fLower);
  ADD_PROCESS_PARAM("float",	"%12.5f",	 "--fendbcv",		coarseBankIn.LowGM);
  ADD_PROCESS_PARAM("float",	"%12.5f",	 "--fendbcv",		coarseBankIn.HighGM);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--sampling",		coarseBankIn.tSampling);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range",		randIn.mMin);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range",		randIn.mMax);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range-bank",	coarseBankIn.mMin);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mass-range-bank",	coarseBankIn.mMax);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--m1",			otherIn.m1);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--m2",			otherIn.m2);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0",		otherIn.psi0);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3",		otherIn.psi3);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0-range",		coarseBankIn.psi0Min);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi0-range",		coarseBankIn.psi0Max);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3-range",		coarseBankIn.psi3Min);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--psi3-range",		coarseBankIn.psi3Max);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--mm",			coarseBankIn.mmCoarse);
  ADD_PROCESS_PARAM("int",	"%6d",		"--number-fcut",	coarseBankIn.numFcutTemplates);
  ADD_PROCESS_PARAM("int",	"%6d",		"--signal",		otherIn.signal);
  ADD_PROCESS_PARAM("int",	"%6d",		"--simulation-type",	randIn.type);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--signal-amplitude",	randIn.SignalAmp);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--noise-amplitude",	randIn.NoiseAmp);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--alpha-bank",		coarseBankIn.alpha);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--alpha-signal",	randIn.param.alpha);
  ADD_PROCESS_PARAM("float",	"%12.5f",	"--freq-moment-bank",	coarseBankIn.fUpper);
  ADD_PROCESS_PARAM("int",	"%6d",		"--template",		otherIn.template);
  ADD_PROCESS_PARAM("int",	"%6d",		"--template-order",	coarseBankIn.order);
  ADD_PROCESS_PARAM("int",	"%6d",		"--signal-order",	randIn.param.order);
  ADD_PROCESS_PARAM("int",	"%6d",		"--nn",			otherIn.ntrials);
  ADD_PROCESS_PARAM("int",	"%6d",		"--seed",		randIn.useed);
  switch (otherIn.NoiseModel){
  case LIGOI:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"LIGOI");
    break;
  case LIGOA:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"LIGOA");
    break;
  case VIRGO:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"VIRGO");
    break;
  case GEO:
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"GEO");
    break;
  case TAMA:    
    ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"TAMA");    
    break;
  case REALPSD:
     ADD_PROCESS_PARAM("string",	"%6s",		"--noisemodel",		"REALPSD");    
    break;
  }
  ADD_PROCESS_PARAM("int",	"%6d",		"--template",		otherIn.template);
  ADD_PROCESS_PARAM("int",	"%6d",		"--signal",		otherIn.signal);

  switch (otherIn.overlapMethod){
  case InQuadrature:
    ADD_PROCESS_PARAM("string",	"%s", 	"OverlapMethod", 	"InQuadrature");
    break;
  case AlphaMaximization:
    ADD_PROCESS_PARAM("string",	"%s", 	"OverlapMethod", 	"AlphaMaximization");
    break;
  }
  
	  
#undef ADD_PROCESS_PARAM
}


/* xml file for the standalone code */
void 
BEPrintResultsXml(InspiralCoarseBankIn   coarseBankIn,
		  RandomInspiralSignalIn randIn,
		  OtherParamIn           otherIn,
		  ResultIn               trigger	      
		  )
{
#define MAXIFO 2
#define BANKEFFICIENCY_PARAMS_ROW \
"         %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,"



#define LIGOLW_XML_BANKEFFICIENCY \
"   <Table Name=\"bankefficiencygroup:bankefficiency:table\">\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0T\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3T\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fI\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:totalMassT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:etaT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass1I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass2I\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:overlap\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:phase\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alpha\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alpha_f\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:layerT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:bin\" Type=\"real_4\"/>\n" \
"      <Stream Name=\"bankefficiencygroup:bankefficiency:table\" Type=\"Local\" Delimiter=\",\">\n"


  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR  fname[256];
  LIGOTimeGPS gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  LALSnprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTRESULT_FILEXML ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );



  if (trigger.ntrial == 1){
    strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    memset( ifo, 0, sizeof(ifo) );
    memcpy( ifo, "MC", sizeof(ifo) - 1 );
    
    
    
    /* -- we start to fill the xml file here --- */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );
    
    
    /* create the process and process params tables */
    templateBank.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->start_time),
			      &accuracy ), &status );
    LAL_CALL( populate_process_table( &status, templateBank.processTable, 
				      PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    
    
    /* write process table */
    LALSnprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
		 ifoName[0], ifoName[1] );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->end_time),
			      &accuracy ), &status );
    
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
				      process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* write process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				      process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* finally write sngl inspiral table */
    
    fprintf(xmlStream.fp,"%s",LIGOLW_XML_BANKEFFICIENCY);
   
    fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW"\n",
	    trigger.psi0_trigger,
	    trigger.psi3_trigger,
	    trigger.psi0_inject, 
	    trigger.psi3_inject,
	    trigger.fend_trigger, 
	    trigger.fend_inject,
	    trigger.totalMass_trigger,
	    trigger.eta_trigger,
	    trigger.mass1_inject,
	    trigger.mass2_inject,
	    trigger.rho_final,
	    trigger.phase,
	    trigger.alpha,
	    trigger.alpha_f, 
	    trigger.layer,
	    trigger.bin);
     if (trigger.ntrial == otherIn.ntrials){
       fprintf( xmlStream.fp, LIGOLW_XML_TABLE_FOOTER );
       fprintf( xmlStream.fp, LIGOLW_XML_FOOTER );
     }
     fclose( xmlStream.fp );
     xmlStream.fp = NULL;
  }
  else 
    {
      xmlStream.fp = fopen(fname,"a+");
      
      if (trigger.ntrial == otherIn.ntrials){
	
	fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW,
		trigger.psi0_trigger,
		trigger.psi3_trigger,
		trigger.psi0_inject, 
		trigger.psi3_inject,
		trigger.fend_trigger, 
		trigger.fend_inject,
		trigger.totalMass_trigger,
		trigger.eta_trigger,
		trigger.mass1_inject,
		trigger.mass2_inject,
		trigger.rho_final,
		trigger.phase,
		trigger.alpha,
		trigger.alpha_f, 
		trigger.layer,
		trigger.bin);
	
	fprintf( xmlStream.fp, LIGOLW_XML_TABLE_FOOTER );
	fprintf( xmlStream.fp, LIGOLW_XML_FOOTER );
      }
      else
	{
	  fprintf(xmlStream.fp,BANKEFFICIENCY_PARAMS_ROW"\n",
		  trigger.psi0_trigger,
		  trigger.psi3_trigger,
		  trigger.psi0_inject, 
		  trigger.psi3_inject,
		  trigger.fend_trigger, 
		  trigger.fend_inject,
		  trigger.totalMass_trigger,
		  trigger.eta_trigger,
		  trigger.mass1_inject,
		  trigger.mass2_inject,
		  trigger.rho_final,
		  trigger.phase,
		  trigger.alpha,
		  trigger.alpha_f, 
		  trigger.layer,
		  trigger.bin);
	}
     
      fclose( xmlStream.fp );
      xmlStream.fp = NULL;
    }  


  
  /* close the output xml file */

}


/* print prototype for condor code (option --print-result-prototype)*/
void 
BEPrintProtoXml(InspiralCoarseBankIn   coarseBankIn,
		  RandomInspiralSignalIn randIn,
		  OtherParamIn           otherIn
		  )
{
#define MAXIFO 2

  LALStatus             status = blank_status;


  MetadataTable         templateBank;
  CHAR  ifo[3];                           /* two character ifo code       */
  LIGOLwXMLStream       xmlStream;
  CHAR  fname[256];
  LIGOTimeGPS gpsStartTime 	= { 0, 0 };    /* input data GPS start time    */
  LIGOTimeGPS gpsEndTime 	= { 0, 0 };      /* input data GPS end time      */
  LALLeapSecAccuracy    accuracy = 1;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];

  MetadataTable         processParamsTable;
  ProcessParamsTable   *this_proc_param = NULL;

  LALSnprintf( fname, sizeof(fname), BANKEFFICIENCY_PRINTPROTO_FILEXML ,
	       ifo, gpsStartTime.gpsSeconds,
	       gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );




    strncpy( ifoName[0], "no", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    strncpy( ifoName[1], "ne", LIGOMETA_IFO_MAX * sizeof(CHAR) );
    memset( ifo, 0, sizeof(ifo) );
    memcpy( ifo, "MC", sizeof(ifo) - 1 );
    
    
    
    /* -- we start to fill the xml file here --- */
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, fname), &status );
    
    
    /* create the process and process params tables */
    templateBank.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->start_time),
			      &accuracy ), &status );
    LAL_CALL( populate_process_table( &status, templateBank.processTable, 
				      PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
    this_proc_param = processParamsTable.processParamsTable = 
      (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
    
    BEFillProc(this_proc_param, coarseBankIn, randIn, otherIn);
    
    memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
    
    
    
    /* write process table */
    LALSnprintf( templateBank.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
		 ifoName[0], ifoName[1] );
    LAL_CALL( LALGPSTimeNow ( &status, &(templateBank.processTable->end_time),
			      &accuracy ), &status );
    
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
	      &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, templateBank, 
				      process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* write process params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
				      process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
				      process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
    
    /* finally write sngl inspiral table */
    
     
      fclose( xmlStream.fp );
      xmlStream.fp = NULL;


}

void
BEReadXmlBank(  LALStatus  *status, 
		 CHAR *bankFileName, 
		 InspiralTemplateList **list,
		 INT4 *nlist,
		InspiralCoarseBankIn coarseIn)
{
  INT4 i;
  InspiralTemplate *tempPars;
  INT4 startTemplate = -1;
  INT4 stopTemplate  = -1;
  INT4                          numTmplts    = 0;
  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;
  InspiralTemplateNode         *tmpltHead    = NULL;
  InspiralTemplateNode         *tmpltCurrent = NULL;
  InspiralMetric metric;


  INITSTATUS( status, 
	      "BEReadXmlBank", BANKEFFICIENCYC );
  ATTATCHSTATUSPTR( status );
    
  /* how many template in the xml file ? */
  numTmplts = InspiralTmpltBankFromLIGOLw( &bankHead, 
					   bankFileName,
					   startTemplate, 
					   stopTemplate );

  if ( numTmplts < 0 )
    {
      fprintf( stderr, "error: unable to read templates from %s\n", 
	       bankFileName );
      exit( 1 );
    }
  else if ( numTmplts == 0 )
    {
      
      fprintf( stderr, "no templates found in template bank file: %s\n"
	       "exiting without searching for events.\n", bankFileName );
       
    }
   else
     fprintf( stderr, "%d templates found in template bank file: %s\n", numTmplts, bankFileName);
  

  
  if ( !
       (tempPars = (InspiralTemplate *)LALCalloc( 1, sizeof(InspiralTemplate) ))
       )
    {
      ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
    }
  
 
  bankCurrent = bankHead; 

  for (i=0; i<numTmplts; i++){
    *nlist = i;
    *list = (InspiralTemplateList*) 
      LALRealloc( *list, sizeof(InspiralTemplateList) * (*nlist + 1) );
    if ( ! *list )
      {
	LALFree( tempPars );
	ABORT( status, LALINSPIRALBANKH_EMEM, LALINSPIRALBANKH_MSGEMEM );
      }
    memset( *list + *nlist, 0, sizeof(InspiralTemplateList) );
          

    LAL_CALL( LALFindChirpCreateTmpltNode( status->statusPtr, 
					   bankCurrent, 
					   &tmpltCurrent ), status->statusPtr );
    if ( !tmpltHead ) tmpltHead = tmpltCurrent;
  
    
    (*list)[*nlist].ID = *nlist;
    (*list)[*nlist].params = *bankCurrent; 
    /* Extra info which are not in the bank */
    (*list)[*nlist].params.fLower = coarseIn.fLower; 
    (*list)[*nlist].params.tSampling = coarseIn.tSampling; 
    (*list)[*nlist].params.fCutoff = coarseIn.fUpper; 
    (*list)[*nlist].params.fFinal = coarseIn.fUpper; 
    (*list)[*nlist].params.order = coarseIn.order; 
    (*list)[*nlist].params.approximant = coarseIn.approximant; 
    (*list)[*nlist].params.massChoice = t03;
    (*list)[*nlist].params.distance =  1.;
    (*list)[*nlist].params.signalAmplitude= 1.;
    /*  (*list)[*nlist].params.nStartPad = 0;
      (*list)[*nlist].params.nEndPad = 0;
    */
    LALInspiralParameterCalc( status->statusPtr,  &((*list)[*nlist].params ));
    CHECKSTATUSPTR( status );
  
    (*list)[*nlist].metric = metric; 

    CHECKSTATUSPTR( status );
    ++(*nlist); 


    bankCurrent = bankCurrent->next ;

  }

 for ( i = 0; i < *nlist - 1 ; ++i )
  {
    (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
    (*list)[i].params.next     = &((*list)[i+1].params);
    (*list)[i].params.fine     = NULL;
  }
  (*list)[i].params.minMatch = (REAL4) coarseIn.mmCoarse;
  (*list)[i].params.next = NULL;
  (*list)[i].params.fine = NULL;


   DETATCHSTATUSPTR( status );
  RETURN( status );


}
