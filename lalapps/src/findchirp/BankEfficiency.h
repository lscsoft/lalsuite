
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
#include <lal/FindChirp.h>
#include <gsl/gsl_histogram.h>

#define MAXIFO 2
#define BANKEFFICIENCY_PARAMS_ROW \
"         %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%d,%d,%d,"
#define BANKEFFICIENCY_PARAMS_ROW_SPACE \
"         %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %d %d %d"

#define LIGOLW_XML_BANKEFFICIENCY \
"   <Table Name=\"bankefficiencygroup:bankefficiency:table\">\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0T\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3T\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0TC\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3TC\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi0I\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:psi3I\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fFinalT\"    Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fFinalTC\"    Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:fFinalI\"    Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:totalMassT\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:etaT\"       Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:totalMassTC\" Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:etaTC\"       Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass1I\"     Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:mass2I\"     Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:snr\"        Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:phase\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alpha\"      Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alphaF\"     Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:layer\"      Type=\"int_4s\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:bin\"        Type=\"int_4s\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:snrC\"       Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:phaseC\"     Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alphaC\"     Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:alphaFC\"    Type=\"real_4\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:layerC\"     Type=\"int_4s\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:binC\"       Type=\"int_4s\"/>\n" \
"      <Column Name=\"bankefficiencygroup:bankefficiency:coaTime\"    Type=\"int_4s\"/>\n" \
"      <Stream Name=\"bankefficiencygroup:bankefficiency:table\"      Type=\"Local\" Delimiter=\",\">\n"



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
#define BANKEFFICIENCY_ALPHAFCONSTRAINT         1 				/* alphaF between 0 and 1	       */
#define BANKEFFICIENCY_QUIETFLAG       	        0 				/* silent 				*/ 
#define BANKEFFICIENCY_AMBIGUITYFUNCTION      	0				/* Print Ambiguity function		*/
#define BANKEFFICIENCY_FASTSIMULATION           0                               /* cheating code (dont use template far from injection; use with care */
#define BANKEFFICIENCY_FMAXIMIZATION      	0				/* Print SNR function of fendBCV	*/
#define BANKEFFICIENCY_PRINTOVERLAP             0				/* Print Best Overlap 			*/
#define BANKEFFICIENCY_PRINTBESTOVERLAP         0				/* Print Best Overlap 			*/
#define BANKEFFICIENCY_PRINTBESTTEMPLATE        0				/* Print Best Template 			*/
#define BANKEFFICIENCY_PRINTSNRHISTO  0				/* Print histogram of the n template correlation. 
										   Useful with no signal option to see the distribution of false alarm rate */


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
  ALPHAFConstraint,
    ALPHAFUnconstraint,
    BOTHAlphaMethod
} AlphaConstraint;


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
  AlphaConstraint alphaFConstraint;                 /* force alpha_F to be constriant between 0 and 1 */
  INT4 signal;				/* name of the random signal to inject 	*/	
  INT4 template;			/* name of the template in the bank 	*/
  INT4 bank;				/* type of bank to use 			*/
  UINT4 ntrials;				/* number of simulations		*/
  INT4 quietFlag;			/* a flag for verbose mode or not	*/
  INT4 ambiguityFunction;		/* do we want to save the ambiguity function ? */
  INT4 lalDebug; 
  INT4 FastSimulation;                  /* target the injection in the bank --> less 
					   computation but Match might be less than the
					   optimal  is noise is injected or template and 
					   signal are not faithful (ie: BCV against physical
					   template */

  INT4 FMaximization;			
  INT4 lengthFactor;			/* multiply estimated length of filters by that factor */
  INT4 PrintOverlap;		
  INT4 PrintFilter;		
  INT4 PrintSNRHisto;		
  INT4 PrintBankOverlap;		/* print match of each templates 	*/
  INT4 PrintBank;			/* print bank of templates 		*/
  INT4 PrintBankXml;			/* print bank of templates 		*/
  INT4 PrintResultXml;			/* print bank of templates 		*/
  INT4 PrintPrototype;
  INT4 PrintPsd;                        /* print the psd used in <x|x>          */
  BinaryInjection binaryInjection; /*injection will be set by the mass-range*/
  INT4 PrintTemplate;  
  INT4 PrintBestOverlap, PrintBestTemplate, extraFinalPrinting ;
  OverlapMethodIn overlapMethod;
  INT4 check;				/* compute only one correlation where bothe 
					   template and signal have the same parameters */
  double m1,m2, psi0,psi3;
  char *inputPSD;
  char *inputXMLBank; 
  DetectorName NoiseModel;
  char *filename;

} OtherParamIn;


typedef struct{
  REAL4  rhoMaxConstraint;
  INT4   rhoBinConstraint;
  REAL4  alphaConstraint;
  REAL4  phaseConstraint;
  REAL4  freqConstraint;
  INT4   layerConstraint;
  INT4   templateNumberConstraint;
  InspiralTemplate bestConstraintTemplate;

  REAL4  rhoMaxUnconstraint;
  INT4   rhoBinUnconstraint;
  REAL4  alphaUnconstraint;
  REAL4  phaseUnconstraint;
  REAL4  freqUnconstraint;
  INT4   layerUnconstraint;
  INT4   templateNumberUnconstraint;
  InspiralTemplate bestUnconstraintTemplate;



} OverlapOutputIn;

/* strucutre to output the results */
typedef struct{
  REAL4 psi0_trigger;
  REAL4 psi3_trigger;
  REAL4 psi0_triggerC;
  REAL4 psi3_triggerC;
  REAL4 psi0_inject;
  REAL4 psi3_inject;
  REAL4 fend_trigger;
  REAL4 fend_triggerC;
  REAL4 fend_inject;
  REAL4 mass1_inject;
  REAL4 mass2_inject;
  REAL4 totalMass_trigger;
  REAL4 eta_trigger;
  REAL4 totalMass_triggerC;
  REAL4 eta_triggerC;
  INT4 layer;
  REAL4 rho_final;
  REAL4 alpha;
  REAL4 alpha_f;
  INT4 bin;
  REAL4 phase;

  INT4 layerC;
  REAL4 rho_finalC;
  REAL4 alphaC;
  REAL4 alpha_fC;
  INT4 binC;
  REAL4 phaseC;

  INT4 coaTime;

  UINT4 ntrial;
} ResultIn;


/* As to be cleaned
 * Function to store the optimal  parameters of the overlap 
 * lmax : number of the layer
 * fMax : frequency cutoff of the template
 * jmax: number of the best template. use if one want to use the FMaximization option 
 * */
void
KeepHighestValues(OverlapOutputIn in , 
		  OverlapOutputIn *out
		  );


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


void
BEInitOverlapOutputIn(OverlapOutputIn *this);

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
		  InspiralWaveOverlapIn   *overlapin,
		  REAL4Vector             *Filter1, 
		  REAL4Vector             *Filter2,
		  BCVMaximizationMatrix    matrix,
		  OtherParamIn             otherIn, 
		  OverlapOutputIn             *OverlapOutput);

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


void 
GetResult(
	  InspiralTemplateList   	**list,
		InspiralTemplate       	injected,
		OverlapOutputIn 	bestOverlapout, 
		ResultIn                *result,
		OtherParamIn                 otherIn );
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
PrintResults(	 ResultIn result);


/* Print each  template bank overlap 
 * */
void
PrintBankOverlap(InspiralTemplateList 	**list,
		 int 			sizeBank,
		 float 			*overlap,
		 InspiralCoarseBankIn 	coarseIn);

/* Print the template bank coordinates in ascii format 
 * */
void BEPrintBank(InspiralCoarseBankIn 	coarse, 
		 InspiralTemplateList 	**list,
		 UINT4 			sizeBank);

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

/* xml file for the standalone code */
void 
BEPrintResultsXml( InspiralCoarseBankIn   coarseBankIn,
		   RandomInspiralSignalIn randIn,
		   OtherParamIn           otherIn,
		   ResultIn trigger
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
		INT4 *sizeBank, 
		InspiralCoarseBankIn coarseIn);


void
BEFillOverlapOutput(InspiralWaveOverlapOut overlapout, 
		    OverlapOutputIn *this);

