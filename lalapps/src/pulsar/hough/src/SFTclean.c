/*-----------------------------------------------------------------------
 *
 * File Name: SFTclean.c
 * Authors:  Krishnan, B. 
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan July 12, 2004
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="SFTbinTestCV">
Author: Krishnan, B.  
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*********************************************** </lalLaTeX> */


#include <lal/SFTClean.h>
#include <glob.h> 
#include <lalapps.h>

RCSID ( "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="SFTCLEANCErrorTable"> */
#define SFTCLEANC_ENORM 0
#define SFTCLEANC_ESUB  1
#define SFTCLEANC_EARG  2
#define SFTCLEANC_EBAD  3
#define SFTCLEANC_EFILE 4

#define SFTCLEANC_MSGENORM "Normal exit"
#define SFTCLEANC_MSGESUB  "Subroutine failed"
#define SFTCLEANC_MSGEARG  "Error parsing arguments"
#define SFTCLEANC_MSGEBAD  "Bad argument values"
#define SFTCLEANC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

extern int lalDebugLevel;


#define MAXFILENAMELENGTH 256
/* defaults chosen for L1 */
#define HARMONICFILE "./S3lines_H1_xavi.txt"
/*#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs"*/
#define INPUTSFTDIR "./S3H1data/"
/*#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean"*/
#define OUTPUTSFTDIR "./test/"
#define STARTFREQ 950.0
#define BANDFREQ 20.0
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define WINDOWSIZE 100
#define MAXBINS 20



/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, SFTCLEANC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              SFTCLEANC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( SFTCLEANC_ESUB, SFTCLEANC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return SFTCLEANC_ESUB;                                  \
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

  static SFTVector       *inputSFTs;
  static SFTCatalog *catalog = NULL;
  static SFTConstraints constraints;

  static LineNoiseInfo   lines, lines2;
  static LineHarmonicsInfo harmonics; 
  
  INT4 nLines=0, count1, nHarmonicSets;

  /* log file and strings */
  FILE   *fpLog=NULL;
  CHAR   *fnameLog=NULL; 
  CHAR   *logstr=NULL; 

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  FILE *fp=NULL;   
  INT4 seed, ranCount;  
  RandomParams *randPar=NULL; 

  /* user input variables */
  BOOLEAN uvar_help;
  CHAR *uvar_harmonicfname;        /* file with harmonics info */
  CHAR *uvar_inputSFTDir;    /* directory for unclean sfts */
  CHAR *uvar_outputSFTDir;   /* directory for cleaned sfts */
  REAL8 uvar_fStart, uvar_fBand;
  INT4  uvar_window, uvar_maxBins;

  /* set defaults */

  lalDebugLevel = 0;
  /* LALDebugLevel must be called before anything else */
  LAL_CALL( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

  uvar_help = FALSE;

  uvar_harmonicfname = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_harmonicfname,HARMONICFILE); 
  
  uvar_inputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_inputSFTDir, INPUTSFTDIR);

  uvar_outputSFTDir = (CHAR *)LALMalloc(256 * sizeof(CHAR));
  strcpy(uvar_outputSFTDir, OUTPUTSFTDIR);

  uvar_fStart = STARTFREQ;
  uvar_fBand = BANDFREQ;  
  uvar_window = WINDOWSIZE;
  uvar_maxBins = MAXBINS;

  /* register user input variables */
  LAL_CALL( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                          &uvar_help),            &status);  
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "inputSFTDir",     'i', UVAR_OPTIONAL, "Input SFT Directory",                         &uvar_inputSFTDir),     &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "outputSFTDir",    'o', UVAR_OPTIONAL, "Output SFT Directory",                        &uvar_outputSFTDir),    &status);
  LAL_CALL( LALRegisterSTRINGUserVar( &status, "harmonicfname",   'H', UVAR_OPTIONAL, "File with list of lines",                     &uvar_harmonicfname),   &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Frequency to start cleaning",                 &uvar_fStart),          &status);
  LAL_CALL( LALRegisterREALUserVar(   &status, "fBand",           'b', UVAR_OPTIONAL, "Frequency Band",                              &uvar_fBand),           &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "window",          'w', UVAR_OPTIONAL, "No. of bins for generating random numbers",   &uvar_window),          &status);
  LAL_CALL( LALRegisterINTUserVar(    &status, "maxBins",         'm', UVAR_OPTIONAL, "Max. No. of bins to clean",                   &uvar_maxBins),         &status);

  /* read all command line variables */
  LAL_CALL( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 

  /********logging the user input variables*********************/
  /* open log file for writing */
  fnameLog = (CHAR *)LALMalloc( 512*sizeof(CHAR));
  strcpy(fnameLog,uvar_outputSFTDir);
  strcat(fnameLog, "/CLEAN_log");
  if ((fpLog = fopen(fnameLog, "w")) == NULL) {
    fprintf(stderr, "Unable to open file %s for writing\n", fnameLog);
    LALFree(fnameLog);
    exit(1);
  }

  /* get the log string */
  LAL_CALL( LALUserVarGetLog(&status, &logstr, UVAR_LOGFMT_CFGFILE), &status);  

  fprintf( fpLog, "## LOG FILE FOR SFT Cleaning\n\n");
  fprintf( fpLog, "# User Input:\n");
  fprintf( fpLog, "#-------------------------------------------\n");
  fprintf( fpLog, logstr);
  LALFree(logstr);



  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  fp=fopen("/dev/urandom", "r");
  /*   if (!fp) { */
  /*      fprintf(stderr,"Error in SFTCleanTest. Error code %i; desc: %s \n",SFTCLEANTESTC_ERANDFILE,SFTCLEANTESTC_MSGERANDFILE); */
  /*      exit(1); */
  /*   } */
  ranCount = fread(&seed, sizeof(seed), 1, fp);
  /*   if (!(ranCount==1)) { */
  /*      fprintf(stderr,"Error in SFTCleanTest. Error code %i; desc: %s \n",SFTCLEANTESTC_ERANDSEED,SFTCLEANTESTC_MSGERANDSEED); */
  /*      exit(1); */
  /*   } */
  fclose(fp);
  LAL_CALL ( LALCreateRandomParams (&status, &randPar, seed), &status );

  /* copy contents of harmonics file into logfile */
  fprintf(fpLog, "\n\n# Contents of harmonics file:\n");
  fclose(fpLog);
  {
    CHAR command[1024] = "";
    sprintf(command, "cat %s >> %s", uvar_harmonicfname, fnameLog);
    system(command);
  }

  /* append an ident-string defining the exact CVS-version of the code used */
  fpLog = fopen(fnameLog, "a");
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



 
  LAL_CALL( LALFindNumberHarmonics (&status, &harmonics, uvar_harmonicfname), &status); 
  nHarmonicSets = harmonics.nHarmonicSets; 

  if (nHarmonicSets > 0)
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
    

      LAL_CALL( LALReadHarmonicsInfo( &status, &harmonics, uvar_harmonicfname ), &status);
      
      nLines = 0;
      for (count1=0; count1 < nHarmonicSets; count1++)
	{
	  nLines += *(harmonics.numHarmonics + count1);
	}
      
      lines.nLines = nLines;
      lines.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      LAL_CALL( LALHarmonics2Lines( &status, &lines, &harmonics), &status);


      lines2.nLines = nLines;
      lines2.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      LAL_CALL( LALChooseLines( &status, &lines2, &lines, uvar_fStart, uvar_fStart + uvar_fBand), &status);
      nLines = lines2.nLines;
    }


  {
    CHAR *tempDir = NULL;
    UINT4 k;
    /* catalog containing single sft -- since here we loop over sfts */
    SFTCatalog catalog1;
    
    catalog1.length = 1;
    
    constraints.detector = NULL;
    
    tempDir = (CHAR *)LALCalloc( 512, sizeof(CHAR));
    strcpy( tempDir, uvar_inputSFTDir);
    strcat( tempDir, "/*SFT*.*");
    
    LAL_CALL( LALSFTdataFind( &status, &catalog, tempDir, &constraints), &status);
    
    /* loop over sfts and clean and write them */
    for ( k = 1; k < catalog->length; k++) {
      
      catalog1.data = catalog->data + k;
      
      LAL_CALL( LALLoadSFTs ( &status, &inputSFTs, &catalog1, uvar_fStart, uvar_fStart + uvar_fBand), &status);
      
      /* clean the sft vector -- in this case vector contains just a single sft */
      LAL_CALL ( LALCleanSFTVector( &status, inputSFTs, uvar_maxBins, uvar_window, &lines2, randPar), &status);
      
      /* write the sft vector */  
      LAL_CALL (LALWriteSFTVector2Dir ( &status, inputSFTs, uvar_outputSFTDir, "Cleaned SFT", "clean"), &status);
      
      /* destroy input sft */
      LAL_CALL ( LALDestroySFTVector( &status, &inputSFTs), &status);
      
    } /* loop over sfts */
    
    LALFree(tempDir);
  }
 
 
  

  /* Free memory */

  LAL_CALL( LALDestroySFTCatalog( &status, &catalog ), &status);

  if (nLines > 0)
    {
      LALFree(lines.lineFreq);
      LALFree(lines.leftWing);
      LALFree(lines.rightWing);

      LALFree(lines2.lineFreq);
      LALFree(lines2.leftWing);
      LALFree(lines2.rightWing);

    }
 
  if (nHarmonicSets > 0)
    {
      LALFree(harmonics.startFreq);
      LALFree(harmonics.gapFreq);
      LALFree(harmonics.numHarmonics);
      LALFree(harmonics.leftWing);
      LALFree(harmonics.rightWing);
    }

  LAL_CALL (LALDestroyUserVars(&status), &status);

  /* 09/09/05 gam; randPar now a parameter for LALCleanCOMPLEX8SFT */
  LAL_CALL ( LALDestroyRandomParams (&status, &randPar), &status);

  LALCheckMemoryLeaks(); 

  return SFTCLEANC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













