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
#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs"
#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean"
#define STARTFREQ 150.0
#define BANDFREQ 300.0
#define MAXFILES 3000 /* maximum number of files to read in a directory */
#define WINDOWSIZE 100

#define USAGE "Usage: %s [-d debuglevel] (0)\n [-w window size] (100)\n [-H harmonics file name] (./harmonicsS2LLO4KC.txt)\n [-i input sft dir] (/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_L1_Funky-v3Calv5DQ30MinSFTs)\n [-o output sft dir] (/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_L1_Funky-v3Calv5DQ30MinSFTs-clean) \n [-f start frequency] (150.0 Hz) \n [-b bandwidth] (300Hz) \n [-h print usage] \n"


/* Usage format string. */


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


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* vvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------------------------------ */
int main(int argc, char *argv[]){ 
  static LALStatus       status;  /* LALStatus pointer */ 
  static SFTtype         *sft;
  static LineNoiseInfo   lines;
  static LineHarmonicsInfo harmonics; 
  
  CHAR *harmonicfname;        /* file with harmonics info */
  CHAR *inputSFTDir;    /* directory for unclean sfts */
  CHAR *outputSFTDir;   /* directory for cleaned sfts */
  INT4 j, window, arg; 
  INT4 nLines=0, count1, nHarmonicSets;
  INT4 mObsCoh;  
  REAL8 fStart, fBand;
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  CHAR tempstr1[256], tempstr2[256]; 

  /* set defaults */
  harmonicfname = HARMONICFILE; 
  fStart = STARTFREQ;
  fBand = BANDFREQ;  
  inputSFTDir = INPUTSFTDIR;
  outputSFTDir = OUTPUTSFTDIR;
  window = WINDOWSIZE;

  /********************************************************/  
  /* Parse argument list.  i stores the current position. */
  /********************************************************/
  arg = 1;
  while ( arg < argc ) {
    /* Parse debuglevel option. */
    if ( !strcmp( argv[arg], "-d" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        lalDebugLevel = atoi( argv[arg++] );
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse harmonicsfile */
    else if ( !strcmp( argv[arg], "-H" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        harmonicfname = argv[arg++];
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse input SFT directory */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        inputSFTDir = argv[arg++];
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse output SFT dir */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        outputSFTDir = argv[arg++];
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse start freq dir */
    else if ( !strcmp( argv[arg], "-f" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fStart = atof(argv[arg++]);
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse start freq dir */
    else if ( !strcmp( argv[arg], "-w" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        window = atof(argv[arg++]);
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse bandwidth  */
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fBand = atof(argv[arg++]);
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* Unrecognized option. */
    else {
      ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return SFTCLEANC_EARG;
    }
  } 

  /* End of argument parsing loop. */
  /******************************************************************/   

  SUB( FindNumberHarmonics (&status, &harmonics, harmonicfname), &status); 
  nHarmonicSets = harmonics.nHarmonicSets; 

  if (nHarmonicSets > 0)
    {
      harmonics.startFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.gapFreq = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.numHarmonics = (INT4 *)LALMalloc(harmonics.nHarmonicSets * sizeof(INT4));
      harmonics.leftWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
      harmonics.rightWing = (REAL8 *)LALMalloc(harmonics.nHarmonicSets * sizeof(REAL8));
    

      SUB( ReadHarmonicsInfo( &status, &harmonics, harmonicfname ), &status);
      
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
    }


{ 
    CHAR     command[256];
    glob_t   globbuf;
    INT4    j;
     
    strcpy(command, inputSFTDir);
    strcat(command, "/*SFT*.*");
    
    globbuf.gl_offs = 1;
    glob(command, GLOB_ERR, NULL, &globbuf);
    
    if(globbuf.gl_pathc==0)
      {
	fprintf(stderr,"No SFTs in directory %s ... Exiting.\n", inputSFTDir);
	return 1;  /* stop the program */
      }
    
    /* we will read up to a certain number of SFT files, but not all 
       if there are too many ! */ 
    mObsCoh = (MAXFILES < globbuf.gl_pathc) ? MAXFILES : globbuf.gl_pathc ;
    
    /* Remember to do the following: 
       globfree(&globbuf); after reading the file names. The file names are 
       globbuf.gl_pathv[fileno]   that one can copy into whatever as:
       strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);  */
    
    for (j=0; j < mObsCoh; j++){
      strcpy(filelist[j],globbuf.gl_pathv[j]);
    }
    globfree(&globbuf);	
  }


  for (j=0; j < mObsCoh; j++)
    { 
      sft=NULL;
      SUB (LALReadSFTfile (&status, &sft, fStart, fStart + fBand, filelist[j]), &status);

      /* clean the sft */
      if (nLines > 0)
	SUB( CleanCOMPLEX8SFT( &status, sft, 2, window, &lines), &status);
      
      /* make the output sft filename */
      sprintf(tempstr1, "%d", sft->epoch.gpsSeconds);
      strcpy(tempstr2,outputSFTDir);
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
    }
 
  if (nHarmonicSets > 0)
    {
      LALFree(harmonics.startFreq);
      LALFree(harmonics.gapFreq);
      LALFree(harmonics.numHarmonics);
      LALFree(harmonics.leftWing);
      LALFree(harmonics.rightWing);
    }

  LALCheckMemoryLeaks(); 

  INFO( SFTCLEANC_MSGENORM );
  return SFTCLEANC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */













