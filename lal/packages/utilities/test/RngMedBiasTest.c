/*
*  Copyright (C) 2007 Badri Krishnan, Jolien Creighton
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
 * File Name: RngMedBiasTest.c
 * Authors: Krishnan, B., Itoh, Y.,
 *
 * Revision: $Id$
 *
 * History:   Created by Krishnan Mar 2, 2004
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */
/************************************ <lalVerbatim file="RngMedBiasTestCV">
Author: Krishnan, B. and Itoh, Y.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*********************************************** </lalLaTeX> */


#include <lal/RngMedBias.h>

NRCSID (RNGMEDBIASTESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="RNGMEDBIASTESTCErrorTable"> */
#define RNGMEDBIASTESTC_ENORM 0
#define RNGMEDBIASTESTC_ESUB  1
#define RNGMEDBIASTESTC_EARG  2
#define RNGMEDBIASTESTC_EBAD  3
#define RNGMEDBIASTESTC_EFILE 4

#define RNGMEDBIASTESTC_MSGENORM "Normal exit"
#define RNGMEDBIASTESTC_MSGESUB  "Subroutine failed"
#define RNGMEDBIASTESTC_MSGEARG  "Error parsing arguments"
#define RNGMEDBIASTESTC_MSGEBAD  "Bad argument values"
#define RNGMEDBIASTESTC_MSGEFILE "Could not create output file"
/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=0;
/* Usage format string. */

#define USAGE "Usage: %s [-d debuglevel] [-b block size] \n"


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, RNGMEDBIASTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              RNGMEDBIASTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( RNGMEDBIASTESTC_ESUB, RNGMEDBIASTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return RNGMEDBIASTESTC_ESUB;                                  \
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
  INT4 arg, blkSize;                         /* Argument counter */
  REAL8 bias;

  /* default values */
  blkSize=7;
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
        ERROR( RNGMEDBIASTESTC_EARG, RNGMEDBIASTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RNGMEDBIASTESTC_EARG;
      }
    }
    else if ( !strcmp( argv[arg], "-b" ) ) {
      if ( argc > arg + 1 ) {
        arg++;
        blkSize = atof(argv[arg++]);
      } else {
        ERROR( RNGMEDBIASTESTC_EARG, RNGMEDBIASTESTC_MSGEARG, 0 );
        LALPrintError( USAGE, *argv );
        return RNGMEDBIASTESTC_EARG;
      }
    }
    /* Unrecognized option. */
    else {
      ERROR( RNGMEDBIASTESTC_EARG, RNGMEDBIASTESTC_MSGEARG, 0 );
      LALPrintError( USAGE, *argv );
      return RNGMEDBIASTESTC_EARG;
    }
  } /* End of argument parsing loop. */
  /******************************************************************/



  SUB( LALRngMedBias( &status, &bias, blkSize), &status);

  printf("The correction factor for block size %d is:  %1.15lf\n", blkSize, bias);

  INFO( RNGMEDBIASTESTC_MSGENORM );
  return RNGMEDBIASTESTC_ENORM;
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



