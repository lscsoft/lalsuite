/*-----------------------------------------------------------------------
 *
 * File Name: SFTbinTest.c
 * Authors: Sintes, A.M., 
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="SFTbinTestCV">
Author: Sintes, A.M.,
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*********************************************** </lalLaTeX> */


#include "./SFTbinIO.h"

NRCSID (SFTBINTESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="SFTBINTESTCErrorTable"> */
#define SFTBINTESTC_ENORM 0
#define SFTBINTESTC_ESUB  1
#define SFTBINTESTC_EARG  2
#define SFTBINTESTC_EBAD  3
#define SFTBINTESTC_EFILE 4

#define SFTBINTESTC_MSGENORM "Normal exit"
#define SFTBINTESTC_MSGESUB  "Subroutine failed"
#define SFTBINTESTC_MSGEARG  "Error parsing arguments"
#define SFTBINTESTC_MSGEBAD  "Bad argument values"
#define SFTBINTESTC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;

#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define SFTBASEFILENAME "./data1/SFT."
#define OUTFILE "./outsft.000000"
#define FILEIN "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_H1_FunkyCal30MinSFTs/CAL_SFT.734365561"
/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-i infile] \n"


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, SFTBINTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              SFTBINTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( SFTBINTESTC_ESUB, SFTBINTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return SFTBINTESTC_ESUB;                                  \
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
  static SFTtype          sft2;
  
  CHAR *fname;               /* The input filename */
  CHAR *outfname;           /* output sft file name */
  INT4 arg;                         /* Argument counter */
  INT4  length;
  INT4 k;  
 
  fname = FILEIN;
  outfname = OUTFILE;
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
        ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTBINTESTC_EARG;
      }
    }
    /* Parse input sft file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTBINTESTC_EARG;
      }
    }   
    /* Unrecognized option. */
    else {
      ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return SFTBINTESTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/
  
  SUB( ReadSFTbinHeader1( &status, &header, fname),  &status );
  printf("I am able to read the header \n");
  printf("length = %d \n", header.length);
  printf("time = %lf \n", header.timeBase);
 
  length = 10;
  sft.length = length;
/*
 *   sft.fminBinIndex = floor(151.298 * header.timeBase + 0.5);
 */
  sft.fminBinIndex = header.fminBinIndex + 5;
  sft.data = NULL;
  sft.data = (COMPLEX8 *)LALMalloc(length* sizeof(COMPLEX8));
    
  SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, fname),  &status );
  printf("I am able to read the data \n");
  printf(" sft fminBinIndex = %d \n", sft.fminBinIndex);
  printf(" sft timeBase = %lf \n", sft.timeBase);

  length = sft.length;  
  k=0;
  while(length-- >0){
    printf("frequency = %9.6g\n", (1.0*sft.fminBinIndex + 1.0*k)/ header.timeBase);
    printf("re = %g im = %g \n"  , sft.data[k].re,sft.data[k].im );
    k++;
  }

  /* write the sft to a file */
  SUB( WriteCOMPLEX8SFTbinData1( &status, &sft, outfname),  &status );
  printf("finished writing the sft\n");

  printf("now trying to reread the sft...\n");
  SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, outfname),  &status );  
  printf("..successful\n");
  printf(" sft fminBinIndex = %d \n", sft.fminBinIndex);
  printf(" sft timeBase = %lf \n", sft.timeBase);
  
  length = sft.length;    
  k=0;
  while(length-- >0){
    printf("frequency = %9.6g\n", (1.0*sft.fminBinIndex + 1.0*k)/ header.timeBase);
    printf("re = %g im = %g \n"  , sft.data[k].re,sft.data[k].im );
    k++;
  }

  length = sft.length;    
  sft2.data = NULL;
  printf("now trying to reread the original sft with another function.\n");
  
  SUB( LALCCreateVector (&status, &(sft2.data), length),  &status ); 
  SUB( ReadSFTtype( &status, &sft2, fname, sft.fminBinIndex),  &status );
  printf("..successful\n");
  printf(" sft2 f0 = %lf \n", sft2.f0);
  printf(" sft2 deltaF = %lf \n", sft2.deltaF);

  length = sft.length;    
  k=0;
  while(length-- >0){
    printf("frequency = %9.6g\n", sft2.f0 + sft2.deltaF*k);
    printf("re = %g im = %g \n"  , sft2.data->data[k].re,sft2.data->data[k].im );
    k++;
  }

 
  SUB( LALCDestroyVector (&status, &(sft2.data) ),  &status );

  LALFree(sft.data);
  LALCheckMemoryLeaks(); 

  INFO( SFTBINTESTC_MSGENORM );
  return SFTBINTESTC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



