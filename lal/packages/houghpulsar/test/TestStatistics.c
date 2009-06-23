/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton, Alicia Sintes Olives
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

/*-----------------------------------------------------------------------
 *
 * File Name: TestStatistics.c
 *
 * Authors: Krishnan, B., Sintes, A.M.,
 *
 * Revision: $Id$
 *
 * History:   Created by Badri Krishnan May 24, 2003
 *
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="TestStatisticsCV">
Author: Krishnan, B., Sintes, A.M.
$Id$
********************************************************* </lalVerbatim> */

/* ************************************************************<lalLaTeX>
\subsection{Program \ \texttt{TestStatistics.c}}
\label{s:TestStatistics.c}
Tests the statistics and the histogram number count of
a given total Hough map.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Usage}
\begin{verbatim}
TestStatistics [-d debuglevel] [-o outfile]
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}
This program creates a Hough map and ...

The \verb@-d@ option sets the debug level to the specified value
\verb@debuglevel@.  The \verb@-o@ flag tells the program to print the histogram
of the Hough number counts to the specified data file \verb@outfile@.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Exit codes}
\vspace{0.1in}
\input{TESTSTATISTICSCErrorTable}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHoughStatistics()
LALHoughHistogram()
LALPrintError()
LALMalloc()
LALFree()
LALCheckMemoryLeaks()
\end{verbatim}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{TestStatisticsCV}}

********************************************</lalLaTeX> */



#include <lal/Statistics.h>
/*#include "./Statistics.h"*/


NRCSID (TESTSTATISTICSC, "$Id$");


/* Error codes and messages */

/************** <lalErrTable file="TESTSTATISTICSCErrorTable"> */
#define TESTSTATISTICSC_ENORM 0
#define TESTSTATISTICSC_ESUB  1
#define TESTSTATISTICSC_EARG  2
#define TESTSTATISTICSC_EBAD  3
#define TESTSTATISTICSC_EFILE 4

#define TESTSTATISTICSC_MSGENORM "Normal exit"
#define TESTSTATISTICSC_MSGESUB  "Subroutine failed"
#define TESTSTATISTICSC_MSGEARG  "Error parsing arguments"
#define TESTSTATISTICSC_MSGEBAD  "Bad argument values"
#define TESTSTATISTICSC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;
#define FILEOUT "OutHistogram.asc"      /* file output */

/* Usage format string. */
#define USAGE "Usage: %s [-d debuglevel] [-o outfile]\n"

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, TESTSTATISTICSC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              TESTSTATISTICSC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( TESTSTATISTICSC_ESUB, TESTSTATISTICSC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return TESTSTATISTICSC_ESUB;                                  \
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

  static LALStatus       status;
  static HOUGHMapTotal   ht;
  static HoughStats      stats;
  static UINT8Vector     hist;

  const CHAR   *fname = NULL;               /* The output filename */
  INT4   arg;                         /* Argument counter */
  UINT4   i;
  FILE   *fp=NULL;                    /* Output file */

  fname = FILEOUT;
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
        ERROR( TESTSTATISTICSC_EARG, TESTSTATISTICSC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTSTATISTICSC_EARG;
      }
    }
    /* Parse output file option. */
    else if ( !strcmp( argv[arg], "-o" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        fname = argv[arg++];
      } else {
        ERROR( TESTSTATISTICSC_EARG, TESTSTATISTICSC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return TESTSTATISTICSC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( TESTSTATISTICSC_EARG, TESTSTATISTICSC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return TESTSTATISTICSC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/

  /* ------------------------------------------------------- */

  /* create a hough map */
  ht.mObsCoh = 1000;
  ht.xSide = 100;
  ht.ySide = 100;
  /* allocate memory for hough map */
  ht.map = (HoughTT *)LALMalloc(ht.xSide * ht.ySide * sizeof(HoughTT));
  /* make up a hough map */
  for (i = 0; i < floor(ht.xSide * ht.ySide / 2.0); i++) ht.map[i] = 500;
  for (i = ceil(ht.xSide * ht.ySide / 2.0); i < ht.xSide * ht.ySide; i++) ht.map[i] = 900;
  SUB( LALHoughStatistics ( &status, &stats, &ht), &status );

  printf(" Maximum number count: %d\n", stats.maxCount);
  printf(" Location: %d  %d\n", stats.maxIndex[0], stats.maxIndex[1]);
  printf(" Minimum number count: %d\n", stats.minCount);
  printf(" Location: %d  %d\n", stats.minIndex[0], stats.minIndex[1]);
  printf(" Average number count: %f\n", stats.avgCount);
  printf(" Standard deviation of number count: %f\n", stats.stdDev);

  /* allocate memory for histogram */

  hist.length = ht.mObsCoh +1;
  hist.data= NULL;
  hist.data = (UINT8 *)LALMalloc((hist.length)*sizeof(UINT8));

  SUB( LALHoughHistogram ( &status, &hist, &ht), &status);

  /* write histogram to a file */
  fp = fopen(fname, "w");

  if ( !fp ){
    ERROR( TESTSTATISTICSC_EFILE, TESTSTATISTICSC_MSGEFILE, 0 );
    return TESTSTATISTICSC_EFILE;
  }

  for (i=0; i < hist.length; i++){
    fprintf(fp,"%d  %llu\n", i, hist.data[i]);
  }


  fclose(fp);
  LALFree(ht.map);
  LALFree(hist.data);
  LALCheckMemoryLeaks();

  INFO(TESTSTATISTICSC_MSGENORM);
  return TESTSTATISTICSC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */






