/*-----------------------------------------------------------------------
 *
 * File Name: SFTCleanTest.c
 * Authors:  Krishnan, B. 
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan August 2005
 *            
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="NormalizeSFTRngMedTestCV">
Author: Krishnan, B.  
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*********************************************** </lalLaTeX> */


#include <lal/NormalizeSFTRngMed.h>
#include <glob.h> 
#include <lal/GeneratePulsarSignal.h>

NRCSID (NORMALIZESFTRNGMEDC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="NORMALIZESFTRNGMEDCErrorTable"> */
#define NORMALIZESFTRNGMEDC_ENORM 0
#define NORMALIZESFTRNGMEDC_ESUB  1
#define NORMALIZESFTRNGMEDC_EARG  2
#define NORMALIZESFTRNGMEDC_EBAD  3
#define NORMALIZESFTRNGMEDC_EFILE 4

#define NORMALIZESFTRNGMEDC_MSGENORM "Normal exit"
#define NORMALIZESFTRNGMEDC_MSGESUB  "Subroutine failed"
#define NORMALIZESFTRNGMEDC_MSGEARG  "Error parsing arguments"
#define NORMALIZESFTRNGMEDC_MSGEBAD  "Bad argument values"
#define NORMALIZESFTRNGMEDC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;


#define MAXFILENAMELENGTH 256
#define INPUTSFTDIR "/home/badkri/fakesfts/"
#define OUTPUTSFTDIR "./test/"
#define FMIN 255.0
#define FMAX 255.1
#define BLKSIZE 25

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, NORMALIZESFTRNGMEDC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              NORMALIZESFTRNGMEDC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( NORMALIZESFTRNGMEDC_ESUB, NORMALIZESFTRNGMEDC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return NORMALIZESFTRNGMEDC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

#define TRUE (1==1)
#define FALSE (1==0)

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){ 
  static LALStatus    status;  /* LALStatus pointer */ 

  static SFTVector  *sft=NULL;

  CHAR *fname=NULL;

  /* user input variables */
  BOOLEAN uvar_help;
  CHAR *uvar_inputSFTDir;    /* directory for unclean sfts */
  CHAR *uvar_outputSFTDir;   /* directory for cleaned sfts */
  REAL8 uvar_fmin, uvar_fmax;
  INT4 uvar_blockSize;
  /* set defaults */

  lalDebugLevel = 0;
  /* LALDebugLevel must be called before anything else */
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  uvar_help = FALSE;

  uvar_fmin = FMIN;
  uvar_fmax = FMAX;

  uvar_blockSize = BLKSIZE;
  
  uvar_inputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_inputSFTDir, INPUTSFTDIR);

  uvar_outputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_outputSFTDir, OUTPUTSFTDIR);

  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",         'h', UVAR_HELP,     "Print this message",    &uvar_help),         &status);  
  SUB( LALRegisterSTRINGUserVar( &status, "inputSFTDir",  'i', UVAR_OPTIONAL, "Input SFT Directory",   &uvar_inputSFTDir),  &status);
  SUB( LALRegisterSTRINGUserVar( &status, "outputSFTDir", 'o', UVAR_OPTIONAL, "Output SFT Directory",  &uvar_outputSFTDir), &status);
  SUB( LALRegisterREALUserVar(   &status, "fmin",         'm', UVAR_OPTIONAL, "Minimum Frequency",     &uvar_fmin),         &status);
  SUB( LALRegisterREALUserVar(   &status, "fmax",         'M', UVAR_OPTIONAL, "Maximum Frequency",     &uvar_fmax),         &status);
  SUB( LALRegisterINTUserVar(    &status, "blockSize",    'b', UVAR_OPTIONAL, "Rng Med block size",    &uvar_blockSize),    &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  fname = (CHAR *)LALMalloc(256*sizeof(CHAR));
  strcpy(fname, uvar_inputSFTDir);
  strcat(fname, "/*SFT*");
  SUB ( LALReadSFTfiles( &status, &sft, uvar_fmin, uvar_fmax,0, fname), &status);

  SUB ( LALNormalizeSFTVect(&status, sft, uvar_blockSize), &status);

  SUB ( LALDestroySFTVector (&status, &sft), &status);

  SUB (LALDestroyUserVars(&status), &status);

  LALFree(fname);

  LALCheckMemoryLeaks(); 

  INFO( NORMALIZESFTRNGMEDC_MSGENORM );
  return NORMALIZESFTRNGMEDC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













