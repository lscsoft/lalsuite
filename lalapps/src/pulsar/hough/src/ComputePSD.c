/*-----------------------------------------------------------------------
 *
 * File Name: ComputePSD.c
 * Authors:  Krishnan, B.  ; Sintes, A. M.
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan 
 *            Modified by Sintes 
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="ComputePSDCV">
Author: Krishnan, B. , Sintes, A.M. 
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*********************************************** </lalLaTeX> */


#include <glob.h> 
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SFTfileIO.h>
#include <lal/Random.h>
#include <lal/PulsarDataTypes.h>
#include <lal/UserInput.h>

NRCSID (COMPUTEPSDC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="ComputePSDCErrorTable"> */
#define COMPUTEPSDC_ENORM 0
#define COMPUTEPSDC_ESUB  1
#define COMPUTEPSDC_EARG  2
#define COMPUTEPSDC_EBAD  3
#define COMPUTEPSDC_EFILE 4

#define COMPUTEPSDC_MSGENORM "Normal exit"
#define COMPUTEPSDC_MSGESUB  "Subroutine failed"
#define COMPUTEPSDC_MSGEARG  "Error parsing arguments"
#define COMPUTEPSDC_MSGEBAD  "Bad argument values"
#define COMPUTEPSDC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;


#define MAXFILENAMELENGTH 256
/* defaults chosen for L1 */
/*#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs"*/
#define INPUTSFTDIR "/home/badkri/L1sfts/"
/*#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean"*/
#define OUTPUTPSDFILE "./psd"
#define STARTFREQ 200.0
#define BANDFREQ 20.0
#define MAXFILES 3000 /* maximum number of files to read in a directory */


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, COMPUTEPSDC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              COMPUTEPSDC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( COMPUTEPSDC_ESUB, COMPUTEPSDC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return COMPUTEPSDC_ESUB;                                  \
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
  static LALStatus       status;  /* LALStatus pointer */ 
  static SFTtype         *sft;
  
  INT4 j, nBins; 
  INT4 count;
  INT4 mObsCoh;  
  REAL4 *periodo=NULL;
  REAL8 f0, deltaF;
  FILE *fpOut=NULL;
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  COMPLEX8 *inputData;

  /* log file and strings */
  FILE   *fpLog=NULL;
  CHAR   *fnameLog=NULL; 
  CHAR   *logstr=NULL; 

  /* user input variables */
  BOOLEAN uvar_help;
  CHAR *uvar_inputSFTDir;    /* directory for unclean sfts */
  CHAR *uvar_outputPSDFILE;   /* directory for cleaned sfts */
  REAL8 uvar_fStart, uvar_fBand;

  /* set defaults */

  lalDebugLevel = 0;
  /* LALDebugLevel must be called before anything else */
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  uvar_help = FALSE;
  
  uvar_inputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_inputSFTDir, INPUTSFTDIR);

  uvar_outputPSDFILE = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_outputPSDFILE, OUTPUTPSDFILE);

  uvar_fStart = STARTFREQ;
  uvar_fBand = BANDFREQ;  


  /* register user input variables */
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                          &uvar_help),            &status);  
  SUB( LALRegisterSTRINGUserVar( &status, "inputSFTDir",     'i', UVAR_OPTIONAL, "Input SFT Directory",                         &uvar_inputSFTDir),     &status);
  SUB( LALRegisterSTRINGUserVar( &status, "outputFILE",      'o', UVAR_OPTIONAL, "Output ASCII PSD file",                       &uvar_outputPSDFILE),   &status);
  SUB( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Frequency to start from",                     &uvar_fStart),          &status);
  SUB( LALRegisterREALUserVar(   &status, "fBand",           'b', UVAR_OPTIONAL, "Frequency Band",                              &uvar_fBand),           &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /********logging the user input variables*********************/
  /* open log file for writing */
  fnameLog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
  strcpy(fnameLog,uvar_outputPSDFILE);
  strcat(fnameLog, "_log");
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }

  /* get the log string */
  SUB( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  

  fprintf( fpLog, "## LOG FILE FOR SFT Cleaning\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);



  /* append an ident-string defining the exact CVS-version of the code used */
  {
    CHAR command[1024] = "";
    fprintf (fpLog, "\n\n# CVS-versions of executable:\n");
    fprintf (fpLog, "# -----------------------------------------\n");
    fclose (fpLog);
    
    sprintf (command, "ident %s | sort -u >> %s", argv[0], fnameLog);
    system (command);	/* we don't check this. If it fails, we assume that */
    			/* one of the system-commands was not available, and */
    			/* therefore the CVS-versions will not be logged */

    LALFree(fnameLog); 
  }
  /*end of logging*********************************************************/

 /*********************************************************************/
  { 
    CHAR     command[256];
    glob_t   globbuf;
    INT4    jj;
     
    strcpy(command, uvar_inputSFTDir);
    strcat(command, "/*SFT*.*");
    
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    if(globbuf.gl_pathc==0)
      {
	fprintf(stderr,"No SFTs in directory %s ... Exiting.\n", uvar_inputSFTDir);
	return 1;  /* stop the program */
      }
    
    /* we will read up to a certain number of SFT files, but not all 
       if there are too many ! */ 
    mObsCoh = (MAXFILES < globbuf.gl_pathc) ? MAXFILES : globbuf.gl_pathc ;
    
    /* Remember to do the following: 
       globfree(&globbuf); after reading the file names. The file names are 
       globbuf.gl_pathv[fileno]   that one can copy into whatever as:
       strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);  */
    
    for (jj=0; jj < mObsCoh; jj++){
      strcpy(filelist[jj],globbuf.gl_pathv[jj]);
    }
    globfree(&globbuf);	
  }

  /*********************************************************************/

  sft=NULL;
  SUB (LALReadSFTfile (&status, &sft, uvar_fStart, uvar_fStart + uvar_fBand, filelist[0]), &status);

  /* read the first sft and calculate its periodogram */
  /* this is because we don't yet know the length of the sfts */
  /* find the length and allocate memory for the periodogram */
  nBins = sft->data->length;
  periodo = (REAL4 *)LALMalloc(nBins * sizeof(REAL4));


  /* copy values of f0 and deltaF */
  f0 = sft->f0;
  deltaF = sft->deltaF;

  for (count= 0; count < nBins; count++)
    {     
      REAL4 re, im;
      inputData = sft->data->data + count;
      re = inputData->re;
      im = inputData->re;
      periodo[count] = re*re + im*im;
    }
  LALDestroySFTtype (&status, &sft);

  /* now read the other sfts and add to periodo */
  for (j=1; j < mObsCoh; j++)
    {
      sft=NULL;
      SUB (LALReadSFTfile (&status, &sft, uvar_fStart, uvar_fStart + uvar_fBand, filelist[j]), &status);

      /* for debugging */
      fprintf(stdout, "Working on SFT No. %d\n", j);
      for (count= 0; count < nBins; count++)
	{
	  REAL4 re, im;
	  inputData = sft->data->data + count;
	  re = inputData->re;
	  im = inputData->re;
	  periodo[count] += re*re + im*im;
	}

      LALDestroySFTtype (&status, &sft);
    }

   /*********************************************************************/
   /* renormalization for one-sided PSD */
   
   {
     REAL8 factor, Tsft;
     INT8  Nsamples;
     
      Tsft= 1.0/deltaF;
      Nsamples = 2 * nBins;
      factor= 2.0*deltaF*Tsft*Tsft/Nsamples/Nsamples/mObsCoh;
      
      for (count= 0; count < nBins; count++){ periodo[count] *= factor; }
    }
   /*********************************************************************/

  /* write periodo to the output file */
  fprintf(stdout, "Now writing output...\n");
  if (  (fpOut = fopen(uvar_outputPSDFILE, "w")) == NULL)
    {
      fprintf(stderr, "Unable to open output file %s for writing...exiting \n", uvar_outputPSDFILE);
      exit(1);
    }
  for (count = 0; count < nBins; count++)
    fprintf(fpOut, "%f   %e\n", f0 + count*deltaF, periodo[count]);
  fclose(fpOut);

  SUB (LALDestroyUserVars(&status), &status);

  LALFree(periodo);

  LALCheckMemoryLeaks(); 

  INFO( COMPUTEPSDC_MSGENORM );
  return COMPUTEPSDC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













