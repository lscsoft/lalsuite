/*-----------------------------------------------------------------------
 *
 * File Name: SFTfileIOTest.c
 * Authors: Sintes, A.M., 
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="SFTfileIOTestCV">
Author: Sintes, A.M., Prix, R.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************

\subsection{Program \texttt{SFTfileIOTest.c}}
\label{s:SFTfileIOTest.c}

Tests the routines in \verb@SFTfileIO.h@.

\subsubsection*{Usage}
\begin{verbatim}
SFTfileIOTest -d debugLevel -i inputSFT
\end{verbatim}

\subsubsection*{Description}

This program has ugly code for testing the SFT file reading and writing.

\subsubsection*{Exit codes}

\input{SFTFILEIOTESTCErrorTable}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SFTfileIOTestCV}}

</lalLaTeX> */


#include <lal/SFTfileIO.h>

NRCSID (SFTFILEIOTESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="SFTFILEIOTESTCErrorTable"> */
#define SFTFILEIOTESTC_ENORM 0
#define SFTFILEIOTESTC_ESUB  1
#define SFTFILEIOTESTC_EARG  2
#define SFTFILEIOTESTC_EBAD  3
#define SFTFILEIOTESTC_EFILE 4

#define SFTFILEIOTESTC_MSGENORM "Normal exit"
#define SFTFILEIOTESTC_MSGESUB  "Subroutine failed"
#define SFTFILEIOTESTC_MSGEARG  "Error parsing arguments"
#define SFTFILEIOTESTC_MSGEBAD  "Bad argument values"
#define SFTFILEIOTESTC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=3;

#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define SFTBASEFILENAME "./data1/SFT."
#define OUTFILE "./outsft.0"
#define FILEIN "./inputsft.0"

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
              __LINE__, SFTFILEIOTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              SFTFILEIOTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return SFTFILEIOTESTC_ESUB;                                  \
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
  static SFTHeader      header;
  static SFTtype          sft2;
  SFTtype *sft1 = NULL;
  REAL8 fmin, fmax;
  
  CHAR *fname;               /* The input filename */
  CHAR *outfname;           /* output sft file name */
  INT4 arg;                         /* Argument counter */
  INT4  length;
  INT4 k;  
  UINT4 fminindex;
 
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
        ERROR( SFTFILEIOTESTC_EARG, SFTFILEIOTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTFILEIOTESTC_EARG;
      }
    }
    /* Parse input sft file option. */
    else if ( !strcmp( argv[arg], "-i" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( SFTFILEIOTESTC_EARG, SFTFILEIOTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return SFTFILEIOTESTC_EARG;
      }
    }   
    /* Unrecognized option. */
    else {
      ERROR( SFTFILEIOTESTC_EARG, SFTFILEIOTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return SFTFILEIOTESTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/
  
  SUB( LALReadSFTheader( &status, &header, fname),  &status );
  printf("I am able to read the header \n");
  printf("length = %d \n", header.length);
  printf("time = %lf \n", header.timeBase);
 
  length = 10;
/*
 *   sft.fminBinIndex = floor(151.298 * header.timeBase + 0.5);
 */
  fminindex = header.fminBinIndex + 5;
  length = 10;
  sft2.data = NULL;
  printf("now trying to reread the original sft with another function.\n");
  
  SUB( LALCCreateVector (&status, &(sft2.data), length),  &status ); 
  SUB( LALReadSFTtype( &status, &sft2, fname, fminindex),  &status );
  printf("..successful\n");
  printf(" sft2 f0 = %lf \n", sft2.f0);
  printf(" sft2 deltaF = %lf \n", sft2.deltaF);


  for (k=0; k < length; k++) {
    printf("k=%d: frequency = %9.6g: ", k, sft2.f0 + sft2.deltaF*k);
    printf("re = %g im = %g \n"  , sft2.data->data[k].re,sft2.data->data[k].im );
  }

  /* write SFT to disk */
  printf ("\nWriting this SFT to disk\n");
  SUB (LALWriteSFTtoFile (&status, &sft2, OUTFILE), &status);


  /* try to do the same thing with ReadSFTfile() */
  printf ("\nNow testing LALReadSFTfile():\n");
  fmin = fminindex / header.timeBase;
  fmax = (fminindex + length) / header.timeBase;

  printf ("\nfmin = %f, fmax = %f\n", fmin, fmax);
  SUB ( LALReadSFTfile (&status, &sft1, fmin, fmax, OUTFILE), &status);

  printf(" sft1 f0 = %lf \n", sft1->f0);
  printf(" sft1 deltaF = %lf \n", sft1->deltaF);
  for (k=0; k < sft1->data->length; k++)
    {
      printf("k=%d: frequency = %9.6g: ", k, sft1->f0 + k * sft1->deltaF);
      printf("re = %g im = %g \n"  , sft1->data->data[k].re, sft1->data->data[k].im );
    }
  

 
  SUB( LALCDestroyVector (&status, &(sft2.data) ),  &status );
  LALCDestroyVector (&status, &(sft1->data) );
  LALFree ( sft1 );

  LALCheckMemoryLeaks(); 

  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
