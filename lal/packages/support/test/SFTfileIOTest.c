/*-----------------------------------------------------------------------
 *
 * File Name: SFTfileIOTest.c
 * Authors: Sintes, A.M., Machenschalk, B.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="SFTfileIOTestCV">
Author: Sintes, A.M., Prix, R., Machenschalk, B.
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


#include <config.h>
#include <lal/SFTfileIO.h>

NRCSID (SFTFILEIOTESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="SFTFILEIOTESTCErrorTable"> */
#define SFTFILEIOTESTC_ENORM 	0
#define SFTFILEIOTESTC_ESUB  	1
#define SFTFILEIOTESTC_EARG  	2
#define SFTFILEIOTESTC_EBAD  	3
#define SFTFILEIOTESTC_EFILE 	4
#define SFTFILEIOTESTC_ESFTDIFF 5

#define SFTFILEIOTESTC_MSGENORM "Normal exit"
#define SFTFILEIOTESTC_MSGESUB  "Subroutine failed"
#define SFTFILEIOTESTC_MSGEARG  "Error parsing arguments"
#define SFTFILEIOTESTC_MSGEBAD  "Bad argument values"
#define SFTFILEIOTESTC_MSGEFILE "Could not create output file"
#define SFTFILEIOTESTC_MSGESFTDIFF "initial and final SFTs differ"

/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=3;

#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define INFILE  "inputsft.0" /* little endian sft */
#define INFILE2 "inputsft.1" /* big endian sft */
#ifndef _MSC_VER
#define OUTFILE1 "./TestOutputSFT.0"
#define OUTFILE2 "./TestOutputSFT.1"
#define FPATTERN "./Test*SFT?[0-9]"
#else
#define OUTFILE1 ".\\TestOutputSFT.0"
#define OUTFILE2 ".\\TestOutputSFT.1"
#define FPATTERN ".\\Test*SFT?[0-9]"
#endif

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


/* ----------------------------------------------------------------------*/
int main(int argc, char *argv[]){ 
  static LALStatus       status;  /* LALStatus pointer */ 
  SFTtype *sft1 = NULL;
  SFTtype *sft2 = NULL;
  SFTtype *sft3 = NULL;
  SFTVector *sftvect = NULL;
  REAL8 fmin, fmax;
    
  const CHAR *fname;              /* The input filename */
  const CHAR *outfname;           /* output sft file name */
  INT4 arg;                       /* Argument counter */
 
  fname = INFILE;
  outfname = OUTFILE1;

  /**********************************************************/  
  /* Parse argument list.  arg stores the current position. */
  /**********************************************************/
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

  fmin = 1008.5;
  fmax = 1009.1;

  printf ("Testing LALWriteSFTfile() and LALReadSFTfile()\n");
  /* read SFT from disk into sft1 */
  /* try to do the same thing with ReadSFTfile() */
  SUB (LALReadSFTfile (&status, &sft1, fmin, fmax, fname), &status);

  /* write SFT to disk */
  SUB (LALWriteSFTfile (&status, sft1, OUTFILE1), &status);

  /* read the written SFT into sft2 */
  SUB (LALReadSFTfile (&status, &sft2, fmin, fmax, OUTFILE1), &status);

  /* write it out again */
  SUB (LALWriteSFTfile (&status, sft2, OUTFILE2), &status);

  /* compare sft1 and sft2 */
  if ( (sft1->epoch.gpsSeconds != sft2->epoch.gpsSeconds)
       || (sft1->epoch.gpsNanoSeconds != sft2->epoch.gpsNanoSeconds)  )
    {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }
  if ( (sft1->f0 != sft2->f0) || (sft1->deltaF != sft2->deltaF) ) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }
  if ( sft1->data->length != sft2->data->length) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }
  if ( memcmp (sft1->data->data, sft2->data->data, sft1->data->length) ) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }

  /* read the big-endian sft if default input file was used */
  if (0 == strcmp(INFILE,fname)) {
    printf ("Testing LALReadSFTfile() with big-endian SFT\n");

    SUB (LALReadSFTfile (&status, &sft3, fmin, fmax, INFILE2), &status);

    /* compare sft3 and sft2 */
    if ( (sft3->epoch.gpsSeconds != sft2->epoch.gpsSeconds)
	 || (sft3->epoch.gpsNanoSeconds != sft2->epoch.gpsNanoSeconds)  )
    {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }
    if ( (sft3->f0 != sft2->f0) || (sft3->deltaF != sft2->deltaF) ) {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }
    if ( sft3->data->length != sft2->data->length) {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }
    if ( memcmp (sft3->data->data, sft2->data->data, sft3->data->length) ) {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }

    LALDestroySFTtype (&status, &sft3);
  }

  LALDestroySFTtype (&status, &sft1);
  LALDestroySFTtype (&status, &sft2);

  /*----------------------------------------------------------------------*/
  printf ("Testing LALReadSFTfiles()\n");

  /* now re-read OUTFILE1 and OUTFILE2 using a pattern */
  SUB (LALReadSFTfiles (&status, &sftvect, fmin, fmax, FPATTERN), &status);
  
  if(sftvect->length != 2)
    {
      LALPrintError ("Exactly TWO SFTs were expected with pattern '%s', but got only %d\n", FPATTERN, sftvect->length);
      return SFTFILEIOTESTC_ESUB;
    }

  /* compare sft1 and sft2 */
  if ( (sftvect->data[0].epoch.gpsSeconds != sftvect->data[1].epoch.gpsSeconds)
       || (sftvect->data[0].epoch.gpsNanoSeconds != sftvect->data[1].epoch.gpsNanoSeconds)  )
    {
      ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
      return SFTFILEIOTESTC_ESFTDIFF;
    }
  if ( (sftvect->data[0].f0 != sftvect->data[1].f0) || (sftvect->data[0].deltaF != sftvect->data[1].deltaF) ) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }
  if ( sftvect->data[0].data->length != sftvect->data[1].data->length) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }
  if ( memcmp (sftvect->data[0].data->data, sftvect->data[1].data->data, sftvect->data[0].data->length) ) {
    ERROR (SFTFILEIOTESTC_ESFTDIFF, SFTFILEIOTESTC_MSGESFTDIFF, 0);
    return SFTFILEIOTESTC_ESFTDIFF;
  }

  LALDestroySFTVector (&status, &sftvect);

  LALCheckMemoryLeaks(); 

  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}
