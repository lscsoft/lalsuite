/*
*  Copyright (C) 2007 Badri Krishnan, Alicia Sintes Olives
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

/**
 * \file
 * \ingroup pulsarApps
 * \author Alicia Sintes Olives, Badri Krishnan
 */

#include "./SFTbin.h"

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
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
/*@}*/


/* Default parameters. */


#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define SFTBASEFILENAME "./data1/SFT."
#define OUTFILE "./outsft.000000"
#define FILEIN "/nfs/morbo/geo600/hannover/sft/S2-LIGO/S2_H1_FunkyCal30MinSFTs/CAL_SFT.734365561"
#define LINEFILE "./linenoiseS2LHO4KC.txt"
/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-o infile] \n"


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, "$Id$", statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
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
  static   COMPLEX8SFTData1 sft;
  static   REAL8Periodogram1 peri;
  static LineNoiseInfo  lines;
  
  CHAR *fname;               /* The input filename */
  CHAR *linefname;           /* file with line noise info */
  CHAR *outfname;           /* output sft file name */
  INT4 arg;                         /* Argument counter */
  INT4  length, nLines;
  INT4 k;  
  /*   FILE *fp=NULL;                     input file */
  /*  UINT4 i,j;                          Index counter, etc */
  /* INT4 k; */

  fname = FILEIN;
  linefname = LINEFILE;
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
      } else {
        ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SFTBINTESTC_EARG;
      }
    }
    /* Parse input sft file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SFTBINTESTC_EARG;
      }
    }  
    else if ( !strcmp( argv[arg], "-l" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        linefname = argv[arg++];
      } else {
        ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return SFTBINTESTC_EARG;
      }
    }  
    /* Unrecognized option. */
    else {
      ERROR( SFTBINTESTC_EARG, SFTBINTESTC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return SFTBINTESTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/


  
  SUB( ReadSFTbinHeader1( &status, &header, fname),  &status );
  printf("I am able to read the header \n");
  printf("length = %d \n", header.length);
  printf("time = %lf \n", header.timeBase);
 
  length = 1810;
  sft.length = length;
  sft.fminBinIndex = floor(301.0 * header.timeBase + 0.5);
  sft.data = NULL;
  sft.data = (COMPLEX8 *)LALMalloc(length* sizeof(COMPLEX8));
    
  SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, fname),  &status );
  printf("I am able to read the data \n");
  printf(" sft fminBinIndex = %d \n", sft.fminBinIndex);
  printf(" sft timeBase = %lf \n", sft.timeBase);
  
  /*  k=0;
  while(length-- >0){
    printf("frequency = %9.6g\n", (1.0*sft.fminBinIndex + 1.0*k)/ header.timeBase);
    printf("re = %g im = %g \n"  , sft.data[k].re,sft.data[k].im );
    k++;
    }*/
  // length = sft.length;

  peri.length = length;
  peri.data = NULL;
  peri.data = (REAL8 *)LALMalloc(length* sizeof(REAL8));

  SUB( COMPLEX8SFT2Periodogram1( &status, &peri, &sft),  &status );
  printf("periodogram computed \n");

  k=0;
  while(length-- >0){
   printf("periodogram[x] =  %g \n"  , peri.data[k] );
   k++;
  }


  /*printf("now start the cleaning routines...\n");

  SUB( FindNumberLines( &status, &lines, linefname ), &status);

  nLines = lines.nLines;
  if (nLines > 0)
    {
      lines.lineFreq = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.leftWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));
      lines.rightWing = (REAL8 *)LALMalloc(nLines * sizeof(REAL8));

      SUB( ReadLineInfo( &status, &lines, linefname ), &status);

      SUB( CleanCOMPLEX8SFT( &status, &sft, 2, &lines), &status);

      LALFree(lines.lineFreq);
      LALFree(lines.leftWing);
      LALFree(lines.rightWing);
    }

  printf("...sft has been cleaned\n");

  length = sft.length;
  SUB( COMPLEX8SFT2Periodogram1( &status, &peri, &sft),  &status );
  printf("cleaned periodogram computed \n");

  k=0;
  while(length-- >0){
   printf("clean periodogram[x] =  %g \n"  , peri.data[k] );
   k++;
   } */

  /* write the sft to a file */
  /*SUB( WriteCOMPLEX8SFT( &status, &sft, outfname),  &status );
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
  } */


  LALFree(sft.data);
  LALFree(peri.data);
  LALCheckMemoryLeaks(); 

  INFO( SFTBINTESTC_MSGENORM );
  return SFTBINTESTC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



