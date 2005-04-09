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


#include "./SFTbin.h"
#include <glob.h> 


NRCSID (SFTCLEANC, "$Id$");

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

INT4 lalDebugLevel=0;


#define MAXFILENAMELENGTH 256
/* defaults chosen for L1 */

#define HARMONICFILE "./harmonicsS2LLO4KC.txt" 
/*#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs"*/
#define INPUTSFTDIR "/home/badkri/L1sfts/"
/*#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean"*/
#define OUTPUTSFTDIR "/home/badkri/lscsoft/lalapps/src/pulsar/hough/src/temp/"
#define STARTFREQ 150.0
#define BANDFREQ 300.0
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
  static SFTtype         *sft;
  static LineNoiseInfo   lines, lines2;
  static LineHarmonicsInfo harmonics; 
  
  INT4 j; 
  INT4 nLines=0, count1, nHarmonicSets;
  INT4 mObsCoh;  
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  CHAR tempstr1[256], tempstr2[256]; 

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
  SUB( LALGetDebugLevel( &status, argc, argv, 'd'), &status);

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
  SUB( LALRegisterBOOLUserVar(   &status, "help",            'h', UVAR_HELP,     "Print this message",                          &uvar_help),            &status);  
  SUB( LALRegisterSTRINGUserVar( &status, "inputSFTDir",     'i', UVAR_OPTIONAL, "Input SFT Directory",                         &uvar_inputSFTDir),     &status);
  SUB( LALRegisterSTRINGUserVar( &status, "outputSFTDir",    'o', UVAR_OPTIONAL, "Output SFT Directory",                        &uvar_outputSFTDir),    &status);
  SUB( LALRegisterSTRINGUserVar( &status, "harmonicfname",   'H', UVAR_OPTIONAL, "File with list of lines",                     &uvar_harmonicfname),   &status);
  SUB( LALRegisterREALUserVar(   &status, "fStart",          'f', UVAR_OPTIONAL, "Frequency to start cleaning",                 &uvar_fStart),          &status);
  SUB( LALRegisterREALUserVar(   &status, "fBand",           'b', UVAR_OPTIONAL, "Frequency Band",                              &uvar_fBand),           &status);
  SUB( LALRegisterINTUserVar(    &status, "window",          'w', UVAR_OPTIONAL, "No. of bins for generating random numbers",   &uvar_window),          &status);
  SUB( LALRegisterINTUserVar(    &status, "maxBins",         'm', UVAR_OPTIONAL, "Max. No. of bins to clean",                   &uvar_maxBins),         &status);

  /* read all command line variables */
  SUB( LALUserVarReadAllInput(&status, argc, argv), &status);

  /* exit if help was required */
  if (uvar_help)
    exit(0); 
 
  SUB( FindNumberHarmonics (&status, &harmonics, uvar_harmonicfname), &status); 
  nHarmonicSets = harmonics.nHarmonicSets; 

  if (nHarmonicSets > 0)
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
    

      SUB( ReadHarmonicsInfo( &status, &harmonics, uvar_harmonicfname ), &status);
      
      nLines = 0;
      for (count1=0; count1 < nHarmonicSets; count1++)
	{
	  nLines += *(harmonics.numHarmonics + count1);
	}
      
      lines.nLines = nLines;
      lines.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      SUB( Harmonics2Lines( &status, &lines, &harmonics), &status);


      lines2.nLines = nLines;
      lines2.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines2.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      SUB( ChooseLines( &status, &lines2, &lines, uvar_fStart, uvar_fStart + uvar_fBand), &status);
      nLines = lines2.nLines;
    }


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


  for (j=0; j < mObsCoh; j++)
    { 
      sft=NULL;
      SUB (LALReadSFTfile (&status, &sft, uvar_fStart, uvar_fStart + uvar_fBand, filelist[j]), &status);

      /* clean the sft */
      if (nLines > 0)
	SUB( CleanCOMPLEX8SFT( &status, sft, uvar_maxBins, uvar_window, &lines2), &status);
      
      /* make the output sft filename */
      sprintf(tempstr1, "%d", sft->epoch.gpsSeconds);
      strcpy(tempstr2,uvar_outputSFTDir);
      strcat(tempstr2, "/CLEAN_SFT.");
      strcat(tempstr2, tempstr1);

      /* write the sft */
      SUB( LALWriteSFTfile( &status, sft, tempstr2),  &status );

      LALDestroySFTtype (&status, &sft);
    }


  /* Free memory */
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

  SUB (LALDestroyUserVars(&status), &status);

  LALCheckMemoryLeaks(); 

  INFO( SFTCLEANC_MSGENORM );
  return SFTCLEANC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













