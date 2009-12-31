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
          "Function call '" #func "' should have failed for this SFT but didn't!\n");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define SHOULD_FAIL_WITH_CODE( func, statusptr, code )					\
do { 											\
  if ( func, (statusptr)->statusCode != code) {						\
    LALPrintError( "Function call '" #func "' should have failed with code " #code ", but returned %d instead.\n",	\
		   (statusptr)->statusCode );						\
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
SFTConstraints empty_constraints;
/*---------- Global variables ----------*/

INT4 lalDebugLevel = 3;

/* ----------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  LALStatus status = empty_status;

  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_constraints;
  SFTVector *sft_vect = NULL;
  MultiSFTVector *multsft_vect = NULL;
  CHAR detector[2] = "H1";
  INT4 crc_check;

  /* band to read from infile.* SFTs */
  REAL8 fMin = 1008.5;
  REAL8 fMax = 1009.1;

  if ( argc == 1)	/* avoid warning */
    argc = 1;

  /* check that mal-formated SFTs are properly detected */
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad1", NULL ), &status, SFTFILEIO_EMERGEDSFT );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad2", NULL ), &status, SFTFILEIO_EMERGEDSFT );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad3", NULL ), &status, SFTFILEIO_EMERGEDSFT );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad4", NULL ), &status, SFTFILEIO_EMERGEDSFT );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad5", NULL ), &status, SFTFILEIO_EMERGEDSFT );

  /* the following (SFT-bad6) has a wrong CRC64 checksum. However, this is
   * not checked in LALSFTdataFind, so it should succeed! */
  SHOULD_WORK( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad6", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad7", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad8", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad9", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad10", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad11", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad12", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad13", NULL ), &status, SFTFILEIO_EHEADER );
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-bad14", NULL ), &status, SFTFILEIO_EHEADER );

  /* now check some crc-checksums */
  SHOULD_WORK( LALCheckSFTs ( &status, &crc_check, TESTDIR "SFT-test1", NULL ), &status );
  if ( crc_check != 0 )
    {
      LALPrintError ("\nLALCheckSFTs(): SFT-test1 has correct checksum but LALCheckSFTs claimed it hasn't.\n\n");
      return crc_check;
    }
  SHOULD_WORK( LALCheckSFTs ( &status, &crc_check, TESTDIR "SFT-bad6", NULL ), &status );
  if ( crc_check != SFTFILEIO_ECRC64 )
    {
      LALPrintError ( "\nLALCheckSFTs() failed to catch invalid CRC checksum in SFT-bad6 \n\n");
      return SFTFILEIOTESTC_ESUB;
    }

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
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test*", NULL ), &status,  SFTFILEIO_EDIFFTSFT );
  /* skip sft nr 4 with has Tsft=50 instead of Tsft=60 */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test[123567]*", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  /* try the same with a ";" separated list of files and of patterns */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog,
				 TESTDIR "SFT-test1;"
				 TESTDIR "SFT-test2;"
				 TESTDIR "SFT-test3;"
				 TESTDIR "SFT-test5;"
				 TESTDIR "SFT-test6;"
				 TESTDIR "SFT-test7", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TESTDIR "SFT-test[123]*;" TESTDIR "SFT-test[567]*", NULL ), &status );
  /* load once as a single SFT-vector (mix of detectors) */
  SHOULD_WORK ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );

/*   SHOULD_WORK ( _LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status ); */
/*   { */
/*     UINT4 alpha; */
/*     CHAR name[100], comment[200]; */

/*     for ( alpha=0; alpha < sft_vect->length; alpha ++ ) */
/*       { */
/* 	REAL8 Tsft = 1.0 / sft_vect->data[alpha].deltaF; */
/* 	sprintf ( name, "SFT-test%d.new", alpha + 1 ); */
/* 	sprintf ( comment, "This is the SFT '%s'!", name ); */
/* 	sft_vect->data[alpha].epoch.gpsSeconds += alpha * Tsft ; */
/* 	SHOULD_WORK ( LALWriteSFT2file( &status, &(sft_vect->data[alpha]), name, comment), &status ); */
/*       } */
/*   } */

  /* load once as a multi-SFT vector */
  SHOULD_WORK ( LALLoadMultiSFTs ( &status, &multsft_vect, catalog, -1, -1 ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* 6 SFTs from 2 IFOs should have been read */
  if ( (sft_vect->length != 6) 	/* either as a single SFTVector */
       || (multsft_vect->length != 2) 	/* or separated by detector */
       || (multsft_vect->data[0]->length != 5) || ( multsft_vect->data[1]->length != 1 ) )
    {
      LALPrintError ( "\nFailed to read in multi-SFT from 2 IFOs 'SFT-test*'!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* ----- v2 SFT writing ----- */
  /* write v2-SFT to disk */
  SHOULD_WORK ( LALWriteSFT2file( &status, &(multsft_vect->data[0]->data[0]), "outputsftv2_v2.sft", "A v2-SFT file for testing!"), &status );

  SHOULD_WORK ( LALWriteSFTVector2Dir( &status, multsft_vect->data[0], ".", "A v2-SFT file for testing!", "test"), &status);

  /* write v2-SFt as a v1-SFT to disk (correct normalization) */
  multsft_vect->data[0]->data[0].epoch.gpsSeconds += 60;	/* shift start-time so they don't look like segmented SFTs! */
  SHOULD_WORK ( LALWrite_v2SFT_to_v1file( &status, &(multsft_vect->data[0]->data[0]), "outputsftv2_v1.sft"), &status );

  SUB ( LALDestroySFTVector ( &status, &sft_vect ), &status );
  SUB ( LALDestroyMultiSFTVector (&status, &multsft_vect ), &status );
  /* ----- read the previous two SFTs back */
  SHOULD_FAIL_WITH_CODE ( LALSFTdataFind ( &status, &catalog, "outputsftv2_*.sft", NULL ), &status, SFTFILEIO_EDETECTOR );
  /* need to set proper detector! */
  constraints.detector = detector;
  SUB ( LALSFTdataFind ( &status, &catalog, "outputsftv2_*.sft", &constraints ), &status);
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );

  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to read back in 'outputsftv2_*.sft'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }
  /* the data of 'outputsftv2_v2.sft' and 'outputsftv2_v1.sft' should agree, as the normalization
   * should be corrected again when reading-in
   */
  {
    UINT4 i;
    UINT4 numBins = sft_vect->data[0].data->length;
    for ( i=0; i < numBins; i++)
      {
	COMPLEX8 *data1 = &(sft_vect->data[0].data->data[i]);
	COMPLEX8 *data2 = &(sft_vect->data[1].data->data[i]);

	if ( (data1->re != data2->re) || (data1->im != data2->im) )
	  {
	    LALPrintError ("\nv1- and v2- SFT differ after writing/reading\n\n");
	    return SFTFILEIOTESTC_ESFTDIFF;
	  }
      } /* for i < numBins */
  }
  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* `----- v1 SFT writing */

  /* read v1-SFTs: 'inputsft.0' and 'inputsft.1' (one is big-endian, the other little-endian!) */
  SUB ( LALSFTdataFind (&status, &catalog, TESTDIR "inputsft.?", &constraints ), &status );
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, fMin, fMax ), &status );
  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) LALPrintError ("\nFailed to read in v1-SFTs 'inputsft.0' and 'inputsft.1'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* write v1-SFT to disk */
  SUB ( LALWriteSFTfile (&status, &(sft_vect->data[0]), "outputsft_v1.sft"), &status);

  /* try to write this v1-SFTs as v2: should fail without detector-info ! */
  strncpy( sft_vect->data[0].name, "??", 2 );
  SHOULD_FAIL (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );

  /* put detector there */
  strcpy ( sft_vect->data[0].name, "H1" );
  SHOULD_WORK (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );

  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  LALCheckMemoryLeaks();

  LALPrintError ("\n\n--------------------------------------------------------------------------------\n");
  LALPrintError ("\n    OK. All tests passed correctly ! (error-messages above are OK!)\n");
  LALPrintError ("\n--------------------------------------------------------------------------------\n");


  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}
