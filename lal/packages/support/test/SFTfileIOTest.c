/*
 * Copyright (C) 2004, 2005 R. Prix, B. Machenschalk, A.M. Sintes
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


/** \file
 * \ingroup SFTfileIO
 * \author R. Prix, B. Machenschalk, A.M. Sintes
 * 
 * \brief Test-code for SFT-fileIO library 
 *
 * $Id$ 
 *
 */

/*---------- INCLUDES ----------*/
#include <config.h>
#include <lal/SFTfileIO.h>

NRCSID (SFTFILEIOTESTC, "$Id$");

/*---------- DEFINES ----------*/

/** \name Error codes */
/*@{*/
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
/*@}*/


#ifndef SRCDIR
#define SRCDIR "."
#endif

#define TESTDIR SRCDIR "/"


/* Default parameters. */
#define MAXFILENAMELENGTH 64
#define NFSIZE 5
#define THRESHOLD 2.0
#define INFILE   	"inputsft.0" /* little endian sft */
#define INFILE2 	"inputsft.1" /* big endian sft */
#ifndef _MSC_VER
#define OUTFILE1 "./TestOutputSFT.0"
#define OUTFILE2 "./TestOutputSFT.1"
#define FPATTERN "./Test*SFT?[0-9]"
#else
#define OUTFILE1 ".\\TestOutputSFT.0"
#define OUTFILE2 ".\\TestOutputSFT.1"
#define FPATTERN ".\\Test*SFT?[0-9]"
#endif

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


#define SHOULD_FAIL( func, statusptr )							\
do { 											\
  if ( func, ! (statusptr)->statusCode ) {						\
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' should have failed for this SFT but didn't!");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define SHOULD_WORK( func, statusptr )							\
do { 											\
  if ( func, (statusptr)->statusCode ) {						\
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' failed but should have worked for this SFT!");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)




/*---------- empty initializers ---------- */
LALStatus empty_status;

/*---------- Global variables ----------*/

INT4 lalDebugLevel = 3;

/* ----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{ 
  LALStatus status = empty_status;

  SFTCatalog *catalog = NULL;
  SFTVector *sft_vect = NULL;

  /* band to read from infile.* SFTs */
  REAL8 fMin = 1008.5;
  REAL8 fMax = 1009.1;
    
  /* check that mal-formated SFTs are properly detected */
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad1", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad2", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad3", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad4", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad5", NULL ), &status );

  /* the following (SFT-bad6) has a wrong CRC64 checksum. However, this is 
   * not checked in LALSFTdataFind, so it should succeed! */
  SHOULD_WORK( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad6", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad7", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad8", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad9", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad10", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad11", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad12", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad13", NULL ), &status );
  SHOULD_FAIL( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad14", NULL ), &status );

  /* check that proper v2-SFTs are read-in properly */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test1", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test2", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test3", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test4", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test5", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test6", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test7", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* now completely read-in a v2 merged-SFT */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-good", NULL ), &status );
  SHOULD_WORK ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* 3 SFTs with 4 frequency-bins should have been read */
  if ( (sft_vect->length != 3) || ( sft_vect->data[0].data->length != 4 ) )
    {
      LALPrintError ( "\nFailed to read in 3 SFTs from merged-SFTfile 'SFT-good'!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* write v1-SFT to disk */
  SUB (LALWriteSFTfile (&status, &(sft_vect->data[2]), "outputsft_v1.sft"), &status);
  /* write v2-SFT to disk */
  SHOULD_WORK (LALWriteSFT2file( &status, &(sft_vect->data[2]), "outputsft_v2.sft", "A v2-SFT file for testing!"), &status );

  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  /* read the previous two SFTs back */
  SUB ( LALSFTdataFind ( &status, &catalog, "outputsft_*.sft", NULL ), &status );
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );

  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to read back in 'outputsft_*.sft'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }
  /* try to write one of those v1-SFTs as v2: should fail without detector-info ! */
  SHOULD_FAIL (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );
  /* put detector there */
  strcpy ( sft_vect->data[0].name, "H1" );
  SHOULD_WORK (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );


  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  /* read v1-SFTs: 'inputsft.0' and 'inputsft.1' (one is big-endian, the other little-endian!) */
  SUB ( LALSFTdataFind (&status, &catalog, TESTDIR "inputsft.?", NULL ), &status );
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, fMin, fMax ), &status );
  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to read in v1-SFTs 'inputsft.0' and 'inputsft.1'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }
  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  LALCheckMemoryLeaks(); 

  LALPrintError ("\n\n--------------------------------------------------------------------------------\n");
  LALPrintError ("\n    OK. All tests passed correctly ! (error-messages above are OK!)\n");
  LALPrintError ("\n--------------------------------------------------------------------------------\n");


  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}
