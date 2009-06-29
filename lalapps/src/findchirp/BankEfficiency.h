/*
*  Copyright (C) 2007 Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <processtable.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_matrix.h>
#include <lalapps.h>

#include <lal/LALNoiseModels.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataTables.h>

#include <lal/LIGOLwXMLHeaders.h>
#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralHeaders.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>



/* Here, I defined my own xml table outside the lal strcuture although
   it can be put  into the liXmlHeader files I guess. I dont want to
   use the lal definition for the  time being in order to avoid any 
   overlap with other users. */


#define myfprintf(fp,oldmacro) PRINT_ ## oldmacro(fp)

#define MAXIFO 2

/* --- the macro to read/print the results --- */
#define BANKEFFICIENCY_PARAMS_ROW \
"       %f, %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%u,%u,%f"

/* --- the macro to read/print the results --- */
#define BANKEFFICIENCY_PARAMS_ROW_SPACE \
"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %u %u %f"

/* --- the XML output --- */
#define PRINT_LIGOLW_XML_BANKEFFICIENCY(fp) (\
fputs( "   <Table Name=\"bankefficiencygroup:bankefficiency:table\">\n", fp) == EOF || \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:mass1\"               Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:mass2\"               Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:psi0_sim\"           Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:psi3_sim\"           Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:tau0\"               Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:tau3\"               Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:tau0_sim\"           Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:tau3_sim\"           Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ecc\"                Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ecc_sim\"            Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ecc_sim_fl\"         Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ffinal\"             Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ffinal_sim\"         Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:mass1_sim\"          Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:mass2_sim\"          Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:inclination_sim\"    Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:polarisation_sim\"   Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:phase_sim\"          Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:snr\"                Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:snr_at_ta\"          Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:phase\"              Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:alpha_f\"            Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:time\"               Type=\"int_4s\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:time_sim\"           Type=\"int_4s\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:nfast\"              Type=\"int_4s\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:nfast_max\"          Type=\"int_4s\"/>\n", fp) == EOF ||  \
fputs( "      <Column Name=\"bankefficiencygroup:bankefficiency:ematch\"             Type=\"real_4\"/>\n", fp) == EOF ||  \
fputs( "      <Stream Name=\"bankefficiencygroup:bankefficiency:table\"              Type=\"Local\" Delimiter=\",\">\n", fp) == EOF )



/* --- some CVS macros --- */
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

/* --- Some flags --- */
#define BANKEFFICIENCY_ALPHAFCONSTRAINT     1 /* alphaF between 0 and 1       */
#define BANKEFFICIENCY_ALPHAFUNCONSTRAINT   0 /* alphaF between 0 and 1       */
#define BANKEFFICIENCY_QUIETFLAG            0 /* silent                       */ 
#define BANKEFFICIENCY_FASTSIMULATION       0 /* cheating code (dont use temp
                                                 late far from injection; use 
                                                 with care                    */
#define BANKEFFICIENCY_PRINTOVERLAP         0 /* Print Best Overlap           */
#define BANKEFFICIENCY_PRINTBESTOVERLAP     0 /* Print Best Overlap           */
#define BANKEFFICIENCY_PRINTBESTTEMPLATE    0 /* Print Best Template          */
#define BANKEFFICIENCY_PRINTSNRHISTO        0 /* Print histogram of all SNRs. */
#define BANKEFFICIENCY_PRINTFILTER          0 /* Print corresponding Filter   */
#define BANKEFFICIENCY_PRINTBANK            0 /* print the template bank      */
#define BANKEFFICIENCY_PRINTBANKXML         0 /* print the template bank      */
#define BANKEFFICIENCY_PRINTRESULT          1 /* print the result (ascii)     */
#define BANKEFFICIENCY_PRINTRESULTXML       0 /* print the template bank      */
#define BANKEFFICIENCY_PRINTPROTOTYPE       0 /* print the prototype          */
#define BANKEFFICIENCY_PRINTBANKOVERLAP     0 /* print template's overlap     */
#define BANKEFFICIENCY_PRINTPSD             0 /* print psd used in <x|x>      */
#define BANKEFFICIENCY_PRINTTEMPLATE        0 /* print BCV's final template   */
#define BANKEFFICIENCY_FAITHFULNESS         0               
#define BANKEFFICIENCY_PROTOTYPE            "BankEfficiency-Proto"  
#define BANKEFFICIENCY_ASCIIRESULTS         "Trigger.dat"
#define BANKEFFICIENCY_XMLBANK              "BankEfficiency-Bank.xml"
#define BANKEFFICIENCY_ASCIIBANK            "BankEfficiency-Bank.dat"
#define BANKEFFICIENCY_ASCIISIGNAL          "BankEfficiency-Signal.dat"
#define BANKEFFICIENCY_PRINTOVERLAP_FILE    "BankEfficiency-BestOverlap.dat"
#define BANKEFFICIENCY_PRINTBANK_FILEASCII  "BankEfficiency-Bank.dat"       
#define BANKEFFICIENCY_PRINTBANK_FILEXML    "BankEfficiency-Bank.xml"       
#define BANKEFFICIENCY_XMLRESULTS           "BankEfficiency-Result.xml"     
#define BANKEFFICIENCY_PRINTPSD_FILE        "BankEfficiency-PSD.dat"        
#define BANKEFFICIENCY_SNRHISTOGRAM         "BankEfficiency-SNR_histrogram.dat"


/* ==================== local structures ==================== */
/* 
 *  */
 
/* --- Eccentric bank parameters --- */
typedef struct {
  REAL8 min;
  REAL8 max;
  REAL8 step;
  INT4 bins;
} EccentricBank;
  

/* --- Type of BCV filtering--- */
typedef enum {
  ALPHAFConstraint,
  ALPHAFUnconstraint,
  BOTHAlphaMethod
} AlphaConstraint;


/* --- Type of injections */
typedef enum {
  NoUserChoice,
  BBH,
  BNS,
  BHNS
} BinaryInjection;

/* --- type of noisemodels --- */
typedef enum{
  UNITY,
  LIGOI,
  LIGOA,
  GEO,
  TAMA,
  VIRGO,
  REALPSD,
  READPSD,
  EGO
} DetectorName;

/* --- type of fast option --- */
typedef enum{
  None,
  EMatch,
  Heuristic1    
} FastOption;

typedef struct{
  AlphaConstraint alphaFConstraint; /* force alpha_F to be  between 0 and 1   */
  INT4 signal;                      /* name of the random signal to inject    */  
  INT4 template;                    /* name of the template in the bank       */
  INT4 ntrials;                     /* number of simulations                  */
  FastOption fastSimulation;              /* target the injection in the bank --> less 
                                       computation but Match might be less than 
                                       the optimal  is noise is injected or 
                                       template andsignal are not faithful (ie:
                                       BCV against physical template          */

  INT4 lengthFactor;                /* multiply estimated length by this number*/
  INT4 printSNRHisto;       
  INT4 printBank;                   /* print bank of templates      */
  INT4 printResultXml;
  INT4 printPrototype;
  INT4 printPsd;                    /* print the psd used in <x|x>            */
  BinaryInjection binaryInjection;  /*injection will be set by the mass-range */
  INT4 printBestOverlap, printBestTemplate, extraFinalPrinting ;

  INT4 BankFromFile;
  CHAR *BankFile;

  INT4 faithfulness;
  INT4 snrAtCoaTime;
  REAL8 m1,m2, psi0,psi3, tau0, tau3;
  DetectorName noiseModel;
  REAL4 maxTotalMass;
  REAL4 minTotalMass;
  INT4 startTime;
  INT4 numSeconds;
  INT4 increaseVector;
  REAL4 signalfFinal;
  INT4 startPhase;
  CHAR *inputPSD;
  INT4 useed;
  INT4 ambiguity;
  REAL8 eMatch;
  REAL8 mmFine;
  REAL8 t0FineMin, t0FineMax, t3FineMin, t3FineMax, t0FineBin, t3FineBin;
  EccentricBank eccentricBank;
  EccentricBank eccentricSignal;
  CHAR tag[256];
  INT4 fastParam1;
}
UserParametersIn;

/* struct to save the results */
typedef struct{
  REAL4  rhoMax;
  INT4   rhoBin;
  REAL4  alpha;
  REAL4  phase;
  REAL4  freq;
  INT4   layer;
  INT4   templateNumber;
  InspiralTemplate bestTemplate;
  REAL4 eccentricity;
  REAL4 spin1_x, spin1_y,spin1_z;
  REAL4 spin2_x, spin2_y,spin2_z;
  REAL4 sourceTheta0, sourcePhi0;
  REAL4 polarisationAngle, inclination;
  InspiralTemplate bestUTemplate;
  REAL4 snrAtCoaTime; 
} OverlapOutputIn;

/* structure to output the results */
typedef struct{
  REAL4 tau0_inject;
  REAL4 tau0_trigger;
  REAL4 tau3_inject;
  REAL4 tau3_trigger;
  REAL4 psi0_trigger;
  REAL4 psi3_trigger;
  REAL4 beta_trigger;
  REAL4 beta_inject;
  REAL4 psi0_inject;
  REAL4 psi3_inject;
  REAL4 fend_trigger;
  REAL4 fend_inject;
  REAL4 mass1_trigger;
  REAL4 mass2_trigger;
  REAL4 mass1_inject;
  REAL4 mass2_inject;
  INT4 layer;
  REAL4 rho_final;
  REAL4 alphaF;
  INT4 bin;
  REAL4 phase;
  UINT4 ntrial;
  UINT4 nfast;
  REAL4 snrAtCoaTime; 
  REAL4 eccentricity;
  REAL4 spin1_x, spin1_y,spin1_z;
  REAL4 spin2_x, spin2_y,spin2_z;
  REAL4 sourceTheta0, sourcePhi0;
  REAL4 inclination, polarisationAngle;
} ResultIn;


/* BCV moments */
typedef struct{
  REAL4Vector a11;
  REAL4Vector a21;
  REAL4Vector a22;
} BankEfficiencyMoments;


/* vectors for BCV storage*/
typedef struct{
  REAL4Vector fm5_3;
  REAL4Vector fm2_3;
  REAL4Vector fm7_6;
  REAL4Vector fm1_2;
} BankEfficiencyPowerVector;


typedef struct{
  BankEfficiencyPowerVector powerVector;
  BankEfficiencyMoments     moments;
  REAL4Vector               FilterBCV1;
  REAL4Vector               FilterBCV2;
} BankEfficiencyBCV;
  	
typedef struct{
  INT4  filteringIndex;
  INT4  ntrials;     /* number of simulations done                            */
  INT4  N;           /* number of simulations to be done                      */
  REAL4 lastEMatch;  /* keep the last value of ematch (reset for each simulation)*/
  INT4  filter_processed;/* number of filtering done in each simulation       */
  INT4  stop;            /* flag to carry on the current simulation           */
  REAL4 bestSNR;         /* keep the best SNR                                 */
  INT4  bestSNRIndex;    /* keep the best SNR index                           */
  REAL4 tau0best;
  REAL4 tau3best;
  INT4  SNRMax; /* Index to know how many filtering ago the max SNr was found */
  REAL4 eMatch;        /* keep track of the original requested ematch         */
  REAL4 bestEMatch;    /* keep track of the best ematch for a given simulation*/
  RandomInspiralSignalIn randIn;
  INT4 fastParam1; 
} BankEfficiencySimulation;
	

      
    

/* As to be cleaned
 * Function to store the optimal  parameters of the overlap 
 * lmax : number of the layer
 * fMax : frequency cutoff of the template
 * jmax: number of the best template. use if one want to use the FMaximization 
 * option. 
 * */
void
BankEfficiencyKeepHighestValues(
  OverlapOutputIn  in , 
  OverlapOutputIn *out,
  InspiralTemplate insptmplt);
          
/* function to create the filters in the BCV overlap.
 * Input:   - Vectors F^(-5/3) and so on
 *      - Matrix of the moments for the amplitude
 *      - indice at which the template starts
 *      - psi0 and psi3 for the phase
 * Output   - The two orthogonal filters 
 * */
void 
BankEfficiencyCreateBCVFilters(
  BankEfficiencyBCV *bankEfficiencyBCV,
  UINT4              kMin,
  UINT4              kMax,
  REAL4              psi0,
  REAL4              psi3);
        
void BankEfficiencyBankPrintAscii(
  MetadataTable          templateBank ,
  UINT4                  numCoarse,
  InspiralCoarseBankIn   coarseBankIn );


void BankEfficiencyBankPrintXML(
  MetadataTable          templateBank,
  InspiralCoarseBankIn   coarseBankIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam);


void BankEfficiencyGenerateInputData(
  LALStatus              *status,
  REAL4Vector            *signal,
  RandomInspiralSignalIn *randIn,
  UserParametersIn        userParam);

void BankEfficiencyInitOverlapOutputIn(
  OverlapOutputIn *this);

/* Function to compute a orthogonal vector in BCV
 * */
void BankEfficiencyGetOrthogonalFilter(
  REAL4Vector *filter);


/* Functon to compute the overlap between an injected signal 
 * and a set of templates based on BCV templates. 
 * It returns the overlap, the phase and alpha parameter. 
 * */
void BankEfficiencyWaveOverlapBCV(
  LALStatus               *status,
  REAL4Vector             *correlation,
  InspiralWaveOverlapIn   *overlapin,
  REAL4Vector             *Filter1, 
  REAL4Vector             *Filter2,
  UserParametersIn         userParam, 
  OverlapOutputIn         *OverlapOutput,
  BankEfficiencyMoments   *moments);

/* Function to store the moments needed by the BCV overlap process 
 * a11, a22, a21 are three vectors which stored the three components 
 * of the matrix needed in the BCV maximization process. 
 * */
void
BankEfficiencyCreateBCVMomentVector(
  BankEfficiencyMoments *moments,
  REAL8FrequencySeries  *psd,
  REAL8                  sampling, 
  REAL8                  fLower,
  INT4                   length);

/* Function to create Vectors of the form y(f)= f^(-a/b). 
 * where f is the frequency. 
 * */
void BankEfficiencyCreateVectorFreqPower( 
  REAL4Vector      *vector,
  InspiralTemplate  params,
  INT4              a,
  INT4              b);


/* function to generate waveform (has to be tested). 
 * */
void BankEfficiencyGenerateWaveform(
  LALStatus              *status,
  REAL4Vector            *signal,
  RandomInspiralSignalIn *params);


void BankEfficiencyGetResult(
  LALStatus        *status,
  InspiralTemplate *list,
  InspiralTemplate  injected,
  OverlapOutputIn   bestOverlapout, 
  ResultIn         *result,
  UserParametersIn  userParam );

/* Init the CoarseBank structure 
 * */
void
BankEfficiencyInitInspiralCoarseBankIn(
    InspiralCoarseBankIn    *coarseIn);

/* Init the random structure 
 * */
void BankEfficiencyInitRandomInspiralSignalIn(
  RandomInspiralSignalIn *randIn);

/* Init the UserParametersIn structure
 * */
void BankEfficiencyInitUserParametersIn(
  UserParametersIn *userParam);

/* Function to initialize the different strucutre */
void BankEfficiencyParametersInitialization(
  InspiralCoarseBankIn   *coarseIn, 
  RandomInspiralSignalIn *randIn, 
  UserParametersIn       *userParam);

/* Parsing function 
 * */
void BankEfficiencyParseParameters(
  INT4                    *argc, 
  CHAR                   **argv,
  InspiralCoarseBankIn    *coarseIn,
  RandomInspiralSignalIn  *randIn,
  UserParametersIn        *userParam);             

/* function to check validity of some parameters
 * */
void BankEfficiencyUpdateParams(
  InspiralCoarseBankIn   *coarseIn,
  RandomInspiralSignalIn *randIn,
  UserParametersIn       *userParam);


/* Default values
 * */
void BankEfficiencySetDefault(
  InspiralCoarseBankIn   *coarseIn, 
  RandomInspiralSignalIn *randIn,
  UserParametersIn       *userParam);

/* Help Function 
 * */
void BankEfficiencyHelp(void);

/* Print Function for final results 
 * */
void BankEfficiencyPrintResults(
  ResultIn                    result,
  RandomInspiralSignalIn      randIn,
  BankEfficiencySimulation    simulation);


/* Print each  template bank overlap 
 * */
void BankEfficiencyPrintBankOverlap(
  InspiralTemplateList **list,
  INT4                   sizeBank,
  REAL4                 *overlap,
  InspiralCoarseBankIn   coarseIn);

/* Print the template bank coordinates in ascii format 
 * */
void BankEfficiencyPrintBank(
  InspiralCoarseBankIn   coarse, 
  InspiralTemplateList **list,
  UINT4                  sizeBank);

/* print the template bank coordinates  in xml format 
 * */
void BankEfficiencyPrintBankXml(
  InspiralTemplateList  *coarseList,
  UINT4                  numCoarse,
  InspiralCoarseBankIn   coarseIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam);

void BankEfficiencyGetMaximumSize(
  LALStatus             *status,               
  RandomInspiralSignalIn randIn,
  InspiralCoarseBankIn   coarseBankIn, 
  UserParametersIn       userParam, 
  UINT4                 *length);


void BankEfficiencyCreatePsd(
  LALStatus                   *status, 
  InspiralCoarseBankIn        *coarseBankIn, 
  RandomInspiralSignalIn      *randIn,
  UserParametersIn             userParam);

/* print an error  message 
 * */
void BankEfficiencyPrintError(char *error);

void BankEfficiencyFillProc(
  ProcessParamsTable          *proc,
  InspiralCoarseBankIn        coarseIn,
  RandomInspiralSignalIn      randIn,
  UserParametersIn            userParam);

/* xml file for the standalone code */
void BankEfficiencyPrintResultsXml( 
  InspiralCoarseBankIn     coarseBankIn,
  RandomInspiralSignalIn   randIn,
  UserParametersIn         userParam,
  ResultIn                 trigger,
  BankEfficiencySimulation simulation);

void BankEfficiencyPrintProtoXml(
  InspiralCoarseBankIn   coarseIn,
  RandomInspiralSignalIn randIn,
  UserParametersIn       userParam);

void BankEfficiencyReadXmlBank(  
  LALStatus             *status, 
  CHAR                  *bankFileName, 
  InspiralTemplateList **list,
  INT4                  *sizeBank, 
  InspiralCoarseBankIn   coarseIn);


void BankEfficiencyCreateBank(
  LALStatus             *status, 
  InspiralCoarseBankIn  *coarseBankIn,  
  InspiralTemplateList **list,
  INT4                  *sizeBank);


void BankEfficiencyCreatePowerVector(
  LALStatus                 *status, 
  BankEfficiencyPowerVector *powerVector,
  RandomInspiralSignalIn     randIn,
  INT4                       length);


void 
BankEfficiencyInspiralOverlapBCV(
  LALStatus                   *status,
  InspiralTemplate            *list,
  UserParametersIn             userParam, 
  RandomInspiralSignalIn      *randIn,
  InspiralWaveOverlapIn       *overlapin,
  OverlapOutputIn             *output,
  REAL4Vector                 *correlation,
  BankEfficiencyBCV           *bankefficiencyBCV);

void BankEfficiencyParseGetInt(CHAR **argv, INT4  *index, INT4 *data);
void BankEfficiencyParseGetDouble(CHAR **argv, INT4  *index, REAL8 *data);
void BankEfficiencyParseGetDouble2(CHAR **argv, INT4  *index, REAL8 *data1, 
  REAL8 *data2);

void BankEfficiencyParseGetString(CHAR **argv, INT4  *index );
CHAR* BankEfficiencyGetStringFromGridType(INT4 input);
CHAR* BankEfficiencyGetStringFromSimulationType(INT4 input);
CHAR* BankEfficiencyGetStringFromDetector(INT4 input);
CHAR* BankEfficiencyGetStringFromTemplate(INT4 input);
CHAR* BankEfficiencyGetStringFromNoiseModel(INT4 input);
CHAR* BankEfficiencyGetStringFromScientificRun(INT4 input);
CHAR* BankEfficiencyGetStringFromFastSimulation(INT4 input);


void BankEfficiencyAscii2Xml(void);

void BankEfficiencyInspiralBankGeneration(
  LALStatus            *status,
  InspiralCoarseBankIn *input,
  SnglInspiralTable   **first,
  INT4                 *ntiles,
  UserParametersIn      userParam);

void BankEfficiencyInspiralCreateFineBank(
  LALStatus             *status,
  InspiralTemplateList **outlist,
  INT4                  *nlist,
  InspiralFineBankIn     fineIn, 
  UserParametersIn       param);


void BankEfficiencyCreateTemplateBank(
  LALStatus              *status,
  InspiralCoarseBankIn   *coarseBankIn,
  MetadataTable          *templateBank,
  SnglInspiralTable      **tmpltHead,
  UserParametersIn        userParam,
  RandomInspiralSignalIn *randIn,
  INT4                   *sizeBank);

void BankEfficiencyUpdateSNRHistogram(
  REAL4Vector   *correlation,
  gsl_histogram *histogramNoise);

static int vrbflg = 0;

REAL4 GetMaxVector(REAL4 *vect, INT4 n);
REAL4 GetMinVectorNonZero(REAL4 *vect, INT4 n);

void BankEfficiencyWaveOverlap(
  LALStatus             *status,
  REAL4Vector            *correlation, 
  InspiralWaveOverlapIn  *overlapin,
  OverlapOutputIn        *test,
  INT4                    startPad);
  
void BankEfficiencySaveVector(
  const char  *filename, 
  REAL4Vector  correlation, 
  REAL4        tSampling);

void BankEfficiencySaveVectorAppend(
  const char  *filename, 
  REAL4Vector  correlation, 
  REAL4        tSampling);

void BankEfficiencyPrintMessage(const char *str);


  
  
void BankEfficiencyPopulateAmbiguityFunction(
  gsl_matrix      *amb1,
  REAL4Vector      correlation,
  INT4             tmpltIndex,
  OverlapOutputIn  outputTemplate,
  INT4             startPad, 
  InspiralTemplate insptmplt);


typedef struct {
REAL4 *mass1;
REAL4 *mass2;
REAL4 *fFinal;
REAL4 *tau0;
REAL4 *tau3;
REAL4 *psi0;
REAL4 *psi3;
REAL4 *alpha;
UINT4 *index;
UINT4 *used;
UINT4 size;
REAL4 *eccentricity;
REAL4 *snr;
REAL4 *gamma0;
REAL4 *gamma1;
REAL4 *gamma2;
REAL4 *gamma3;
REAL4 *gamma4;
REAL4 *gamma5;
Approximant approximant;
} Mybank;    

  
/* --- initialise and create a template bank structure for the fast option ---*/   
void BankEfficiencyInitMyBank(
  Mybank            *mybank, 
  INT4              *sizeBank,
  SnglInspiralTable *tmpltHead,
  UserParametersIn   userParam);
  
  
/* --- initialise the eccentric template bank parameters --- */  
void BankEfficiencyEccentricBankInit(
  UserParametersIn *userParam);
  
/* --- print the ambiguity function in a file --- */  
void BankEfficiencyPrintAmbiguity(
  UserParametersIn userParam,
  INT4             sizebank,
  gsl_matrix       *amb1  
);

void BankEfficiencyError(const char * str);
void BankEfficiencyCompare(REAL4 a, REAL4 b, const char *str);
void BankEfficiencyValidity(REAL4 a,  REAL4 min,  REAL4 max, const char * str);
REAL4 BankEfficiencyRandom(REAL4 min, REAL4 max);


void GetClosestValidTemplate(Mybank bank, RandomInspiralSignalIn randIn, 
    UINT4 *fast_index);
  
void GetTau03FromMasses(RandomInspiralSignalIn randIn,REAL4 *tau0, REAL4 *tau3);
  
REAL4 BankEfficiencyComputeEMatch(
  RandomInspiralSignalIn *randIn, Mybank mybank, INT4 index);

void BankEfficiencyCreateListFromTmplt(
  LALStatus         *status,
  InspiralTemplate  *insptmplt, 
  Mybank             mybank,
  INT4               index);

void  BankEfficiencyFinalise(
  LALStatus              *status,
  Mybank                  mybank,
  OverlapOutputIn         overlapOutputBestTemplate,
  RandomInspiralSignalIn  randIn,
  UserParametersIn        userParam,
  BankEfficiencySimulation simulation,
  InspiralCoarseBankIn    coarseBankIn);

void BankEfficiencyInitAmbiguity(gsl_matrix *amb, INT4 sizeBank);

 
void BankEfficiencyReadBankFromFile (
  LALStatus            *status,
  SnglInspiralTable    **first,
  INT4                 *nlist, 
  InspiralCoarseBankIn *coarseIn,
  UserParametersIn     userParam);
