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
#define LINEFILE "./linenoiseS2LHO4KC.txt"
#define HARMONICFILE "./harmonicsS2LHO4KC.txt" 
#define INPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_H1_Funky-v3Cal30MinSFTs"
#define OUTPUTSFTDIR "/nfs/morbo/geo600/hannover/sft/S2-LIGO-clean/S2_H1_FunkyClean-v3Cal30MinSFTs"
#define STARTFREQ 150.0
#define BANDFREQ 300.0
#define MAXFILES 3000 /* maximum number of files to read in a directory */

/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-l linefile] [-h harmonics file name] [-o output file] [-f start frequency] [-b bandwidth] \n"


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
  static SFTHeader1      header;
  static COMPLEX8SFTData1 sft;
  static LineNoiseInfo  lines1, lines2;
  static LineHarmonicsInfo harmonics; 
  
  CHAR *fname;
  CHAR *linefname;           /* file with line noise info */
  CHAR *harmonicfname;        /* file with harmonics info */
  CHAR *inputSFTDir;    /* directory for unclean sfts */
  CHAR *outputSFTDir;   /* directory for cleaned sfts */
  INT4 j, arg;                         /* Argument counter */
  INT4 nLines1, nLines2, count1, nHarmonicSets;
  INT4 mObsCoh;  
  REAL8 fStart, fBand, timeBase, deltaF;
  INT8 fStartBin, fBandBin;  
  CHAR filelist[MAXFILES][MAXFILENAMELENGTH];
  CHAR tempstr1[256], tempstr2[256];
  INT4 temp; 

  /* set defaults */
  linefname = LINEFILE;
  harmonicfname = HARMONICFILE; 
  fStart = STARTFREQ;
  fBand = BANDFREQ;  
  inputSFTDir = INPUTSFTDIR;
  outputSFTDir = OUTPUTSFTDIR;

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
    /* parse linefile */
    else if ( !strcmp( argv[arg], "-l" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        linefname = argv[arg++];
      } else {
        ERROR( SFTCLEANC_EARG, SFTCLEANC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTCLEANC_EARG;
      }
    }  
    /* parse harmonicsfile */
    else if ( !strcmp( argv[arg], "-h" ) ) {
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
  
  {  
    CHAR     command[256];
    glob_t   globbuf;
    UINT4    k; 

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
    mObsCoh = (MAXFILES < globbuf.gl_pathc)? MAXFILES : globbuf.gl_pathc;
    
    /* Remember to do the following: 
       globfree(&globbuf); after reading the file names. The file names are 
       globbuf.gl_pathv[fileno]   that one can copy into whatever as:
       strcpy(filelist[fileno],globbuf.gl_pathv[fileno]);  */
  
    for (k=0; k < mObsCoh; k++)
      {
	strcpy(filelist[k],globbuf.gl_pathv[k]);
      }
    globfree(&globbuf);	
  }
   
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
      
      nLines1 = 0;
      for (count1=0; count1 < nHarmonicSets; count1++)
	{
	  nLines1 += *(harmonics.numHarmonics + count1);
	}
      
      lines1.nLines = nLines1;
      lines1.lineFreq = (REAL8 *)LALMalloc(nLines1 * sizeof(REAL8));
      lines1.leftWing = (REAL8 *)LALMalloc(nLines1 * sizeof(REAL8));
      lines1.rightWing = (REAL8 *)LALMalloc(nLines1 * sizeof(REAL8));

      SUB( Harmonics2Lines( &status, &lines1, &harmonics), &status);
    }

  SUB( FindNumberLines( &status, &lines2, linefname ), &status);
  nLines2 = lines2.nLines;
 
  if (nLines2 > 0)
    {
      lines2.lineFreq = (REAL8 *)LALMalloc(nLines2 * sizeof(REAL8));
      lines2.leftWing = (REAL8 *)LALMalloc(nLines2 * sizeof(REAL8));
      lines2.rightWing = (REAL8 *)LALMalloc(nLines2 * sizeof(REAL8));
      
      SUB( ReadLineInfo( &status, &lines2, linefname ), &status);
    }


  SUB( ReadSFTbinHeader1( &status, &header, filelist[0] ), &status );
  timeBase = header.timeBase; /* Coherent integration time */
  deltaF = 1./timeBase;  /* The frequency resolution */
  fStartBin = (INT8)(fStart * timeBase);
  fBandBin = (INT8)(fBand * timeBase);

  sft.length = fBandBin; 
  sft.fminBinIndex = fStartBin; 
  sft.data = NULL;
  sft.data = (COMPLEX8 *)LALMalloc(fBandBin * sizeof(COMPLEX8));

  for (j=0; j < mObsCoh; j++)
    { 
      fname = filelist[j];
      SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, fname ), &status );
      SUB( CleanCOMPLEX8SFT( &status, &sft, 2, &lines1), &status);
      SUB( CleanCOMPLEX8SFT( &status, &sft, 2, &lines2), &status);
      
      /* make the output sft filename */
      sprintf(tempstr1, "%d",j);
      strcpy(tempstr2,outputSFTDir);
      strcat(tempstr2, "/CLEAN_SFT.");  
      for (temp = 0; temp < 5-strlen(tempstr1); temp++)
	strcat(tempstr2, "0");
      strcat(tempstr2, tempstr1);
      /* write the sft to the file */
      SUB( WriteCOMPLEX8SFT( &status, &sft, tempstr2),  &status );
    }


  /* Free memory */
  LALFree(lines1.lineFreq);
  LALFree(lines1.leftWing);
  LALFree(lines1.rightWing);
 
  LALFree(lines2.lineFreq);
  LALFree(lines2.leftWing);
  LALFree(lines2.rightWing);

  LALFree(harmonics.startFreq);
  LALFree(harmonics.gapFreq);
  LALFree(harmonics.numHarmonics);
  LALFree(harmonics.leftWing);
  LALFree(harmonics.rightWing);


  LALFree(sft.data);
  LALCheckMemoryLeaks(); 

  INFO( SFTCLEANC_MSGENORM );
  return SFTCLEANC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



