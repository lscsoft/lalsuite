/*
*  Copyright (C) 2007 Alicia Sintes Olives
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
 * \author Alicia Sintes Olives
 */

#include "./PeakSelect.h"

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define PEAKSELECTTESTC_ENORM 0
#define PEAKSELECTTESTC_ESUB  1
#define PEAKSELECTTESTC_EARG  2
#define PEAKSELECTTESTC_EBAD  3
#define PEAKSELECTTESTC_EFILE 4

#define PEAKSELECTTESTC_MSGENORM "Normal exit"
#define PEAKSELECTTESTC_MSGESUB  "Subroutine failed"
#define PEAKSELECTTESTC_MSGEARG  "Error parsing arguments"
#define PEAKSELECTTESTC_MSGEBAD  "Bad argument values"
#define PEAKSELECTTESTC_MSGEFILE "Could not create output file"
/*@}*/


/* Default parameters. */


#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define SFTBASEFILENAME "./data1/SFT."
#define OUTFILENAME "./OutHough.asc"
#define FILEIN "./data1/SFT.00001"
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
    ERROR( PEAKSELECTTESTC_ESUB, PEAKSELECTTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return PEAKSELECTTESTC_ESUB;                                  \
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
  static   LALStatus         status;  /* LALStatus pointer */
  static   SFTHeader1        header;
  static   COMPLEX8SFTData1  sft;
  static   REAL8Periodogram1 peri;
  static   REAL8Periodogram1 psd;
  static   UCHARPeakGram     pg1;
  static   HOUGHPeakGram     pgOut;

  CHAR   *fname = NULL;               /* The output filename */
  INT4   arg;                         /* Argument counter */
  INT4   length;
  REAL8  mean;
  REAL8  threshold;
  INT4   nPeaks;
  UINT2  block;

  /* FILE *fp=NULL;  input file */
  /* UINT4 i,j;   Index counter, etc */
  /* INT4 k; */

  fname = FILEIN;

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
        ERROR( PEAKSELECTTESTC_EARG, PEAKSELECTTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return PEAKSELECTTESTC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( PEAKSELECTTESTC_EARG, PEAKSELECTTESTC_MSGEARG, 0 );
        XLALPrintError( USAGE, *argv );
        return PEAKSELECTTESTC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( PEAKSELECTTESTC_EARG, PEAKSELECTTESTC_MSGEARG, 0 );
      XLALPrintError( USAGE, *argv );
      return PEAKSELECTTESTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/


  SUB( ReadSFTbinHeader1( &status, &header, fname),  &status );
  printf("I am able to read the header \n");

  /*  length = 5; */
  length = 20;
  block = 5;
  
  sft.length = length;
  sft.fminBinIndex= header.fminBinIndex +3;
  sft.data = NULL;
  sft.data = (COMPLEX8 *)LALMalloc(length* sizeof(COMPLEX8));

  SUB( ReadCOMPLEX8SFTbinData1( &status, &sft, fname),  &status );
  printf("I am able to read the data \n");

  peri.length = length;
  peri.data = NULL;
  peri.data = (REAL8 *)LALMalloc(length* sizeof(REAL8));
  
  psd.length = length;
  psd.data = NULL;
  psd.data = (REAL8 *)LALMalloc(length* sizeof(REAL8));

  SUB( COMPLEX8SFT2Periodogram1( &status, &peri, &sft),  &status );
  printf("periodogram computed \n");

  SUB( LALComputeMeanPower ( &status, &mean, &peri),  &status );
  printf("Mean of periodogram computed = %g\n" , mean);

  SUB( LALPeriodo2PSDrng ( &status, &psd, &peri, &block),  &status );
  printf("rng median computed \n");

  /* threshold = 1.5; */
  threshold = 1.5*mean;
  
  pg1.length = length;
  pg1.data = NULL;
  pg1.data = (UCHAR *)LALMalloc(length* sizeof(UCHAR));

  SUB( LALSelectPeakWhiteNoise( &status, &pg1, &threshold, &peri), &status );

  nPeaks = pg1.nPeaks;
  printf("Number of Peaks selected = %d\n" , nPeaks);

  pgOut.length = nPeaks;
  pgOut.peak = NULL;
  pgOut.peak = (INT4 *)LALMalloc(nPeaks* sizeof(INT4));

  SUB( LALUCHAR2HOUGHPeak( &status, &pgOut, &pg1), &status );
  printf("Peaks converted \n" );
  printf("First peak is  = %d\n" , pgOut.peak[0] );

  LALFree(sft.data);
  LALFree(peri.data);
  LALFree(psd.data);
  LALFree(pg1.data);
  LALFree(pgOut.peak);
  LALCheckMemoryLeaks();

  INFO( PEAKSELECTTESTC_MSGENORM );
  return PEAKSELECTTESTC_ENORM;
}

