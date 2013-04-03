/*
 * Copyright (C) 2010 Karl Wette
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

/*---------- INCLUDES ----------*/
#include <config.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTutils.h>
#include <lal/Units.h>

/*---------- DEFINES ----------*/
/** \file
 * \ingroup SFTfileIO_h
 * \author R. Prix, B. Machenschalk, A.M. Sintes
 *
 * \brief Test-code for SFT-fileIO library
 *
 *
 */

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


/** \cond DONT_DOXYGEN */
/* Default parameters. */

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )					\
  do {									\
  if ( lalDebugLevel & LALERROR )					\
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"	\
		    "        %s %s\n", (code), *argv, __FILE__,		\
		    __LINE__, "$Id$", statement ? statement :	\
		    "", (msg) );					\
} while (0)

#define INFO( statement )					      \
  do {								      \
    if ( lalDebugLevel & LALINFO )				      \
      XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"   \
		      "        %s\n", *argv, __FILE__, __LINE__,      \
		      "$Id$", (statement) );		      \
  } while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,		     \
           "Function call \"" #func "\" failed:" );                  \
    return SFTFILEIOTESTC_ESUB;					     \
  }								     \
} while (0)


#define SHOULD_FAIL( func, statusptr )							\
do { 											\
  if ( func, ! (statusptr)->statusCode ) {						\
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' should have failed for this SFT but didn't!\n");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
  else   xlalErrno = 0;							                \
} while(0)

#define SHOULD_FAIL_WITH_CODE( func, statusptr, code )					\
do { 											\
  if ( func, (statusptr)->statusCode != code) {						\
    XLALPrintError( "Function call '" #func "' should have failed with code " #code ", but returned %d instead.\n",	\
		   (statusptr)->statusCode );						\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
  else   xlalErrno = 0;							                \
} while(0)


#define SHOULD_WORK( func, statusptr )							\
do { 											\
  if ( func, (statusptr)->statusCode ) {						\
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' failed but should have worked for this SFT!");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

/*---------- empty initializers ---------- */
LALStatus empty_status;
SFTConstraints empty_constraints;
/*---------- Global variables ----------*/

/* ----------------------------------------------------------------------*/

static int CompareSFTVectors(SFTVector *sft_vect, SFTVector *sft_vect2);
static int CompareSFTVectors(SFTVector *sft_vect, SFTVector *sft_vect2)
{
  UINT4 sft,bin;
  if (sft_vect->length != sft_vect2->length) {
    XLALPrintError ( "CompareSFTVectors(): vector lengths differ!\n");
    return(-1);
  }
  for(sft=0; sft < sft_vect->length; sft++) {
    SFTtype sft1 = sft_vect->data[sft];
    SFTtype sft2 = sft_vect2->data[sft];
    if ((sft1.epoch.gpsSeconds != sft1.epoch.gpsSeconds) ||
	(sft1.epoch.gpsNanoSeconds != sft1.epoch.gpsNanoSeconds)) {
      XLALPrintError ( "CompareSFTVectors(): SFT#%u epochs differ (%f/%f)!\n",
		       sft, GPS2REAL8(sft1.epoch), GPS2REAL8(sft2.epoch) );
      return(-1);
    }
    if (!sft1.name || !sft2.name || strcmp(sft1.name,sft2.name)) {
      XLALPrintError ( "CompareSFTVectors(): SFT#%u names differ!\n", sft);
      return(-1);
    }
    if (sft1.f0 != sft2.f0) {
      XLALPrintError ( "CompareSFTVectors(): f0 of SFT#%u differ (%f/%f)!\n",
		       sft, sft1.f0, sft2.f0 );
      return(-1);
    }
    if (sft1.deltaF != sft2.deltaF) {
      XLALPrintError ( "CompareSFTVectors(): deltaF of SFT#%u differ (%f/%f)!\n",
		       sft, sft1.deltaF, sft2.deltaF );
      return(-1);
    }
    if (XLALUnitCompare(&sft1.sampleUnits,&sft2.sampleUnits)) {
      CHAR buf1[256], buf2[256];
      if(!XLALUnitAsString(buf1,256,&sft1.sampleUnits))
	*buf1 = '\0';
      if(!XLALUnitAsString(buf2,256,&sft2.sampleUnits))
	*buf2 = '\0';
      XLALPrintError ( "CompareSFTVectors(): Units of SFT#%u differ (%s/%s)!\n",
		       sft,buf1,buf2 );
      return(-1);
    }
    if (sft1.data->length != sft2.data->length) {
      XLALPrintError ( "CompareSFTVectors(): lengths of SFT#%u differ!\n", sft);
      return(-1);
    }
    for(bin=0; bin < sft1.data->length; bin++) {
      if((crealf(sft1.data->data[bin]) != crealf(sft2.data->data[bin])) ||
	 (cimagf(sft1.data->data[bin]) != cimagf(sft2.data->data[bin]))) {
	XLALPrintError ( "CompareSFTVectors(): bins %u of SFT#%u differ!\n", sft, bin);
	return(-1);
      }
    }
  }
  return(0);
}

int main(int argc, char *argv[])
{
  const char *fn = __func__;
  LALStatus status = empty_status;

  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_constraints;
  SFTVector *sft_vect = NULL;
  SFTVector *sft_vect2 = NULL;
  MultiSFTVector *multsft_vect = NULL;
  MultiSFTVector *multsft_vect2 = NULL;
  CHAR detector[2] = "H1";
  INT4 crc_check;

  /* band to read from infile.* SFTs */
  REAL8 fMin = 1008.5;
  REAL8 fMax = 1009.1;


  if ( argc == 1)	/* avoid warning */
    argc = 1;

  /* check that mal-formated SFTs are properly detected */
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad1", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad2", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad3", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad4", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad5", NULL ), &status);

  /* the following (SFT-bad6) has a wrong CRC64 checksum. However, this is
   * not checked in LALSFTdataFind, so it should succeed! */
  SHOULD_WORK( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad6", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad7", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad8", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad9", NULL ), &status);
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad10", NULL ), &status );
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad11", NULL ), &status );
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad12", NULL ), &status );
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad13", NULL ), &status );
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-bad14", NULL ), &status );

  /* now check some crc-checksums */
  SHOULD_WORK( LALCheckSFTs ( &status, &crc_check, TEST_DATA_DIR "SFT-test1", NULL ), &status );
  if ( crc_check != 0 )
    {
      XLALPrintError ("\nLALCheckSFTs(): SFT-test1 has correct checksum but LALCheckSFTs claimed it hasn't.\n\n");
      return crc_check;
    }
  SHOULD_WORK( LALCheckSFTs ( &status, &crc_check, TEST_DATA_DIR "SFT-bad6", NULL ), &status );
  if ( crc_check != SFTFILEIO_ECRC64 )
    {
      XLALPrintError ( "\nLALCheckSFTs() failed to catch invalid CRC checksum in SFT-bad6 \n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* check that proper v2-SFTs are read-in properly */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test1", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test2", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test3", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test4", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test5", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test6", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test7", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* now completely read-in a v2 merged-SFT */
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test*", NULL ), &status );
  /* skip sft nr 4 with has Tsft=50 instead of Tsft=60 */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test[123567]*", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  /* try the same with a ";" separated list of files and of patterns */
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog,
				 TEST_DATA_DIR "SFT-test1;"
				 TEST_DATA_DIR "SFT-test2;"
				 TEST_DATA_DIR "SFT-test3;"
				 TEST_DATA_DIR "SFT-test5;"
				 TEST_DATA_DIR "SFT-test6;"
				 TEST_DATA_DIR "SFT-test7", NULL ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );
  SHOULD_WORK ( LALSFTdataFind ( &status, &catalog, TEST_DATA_DIR "SFT-test[123]*;" TEST_DATA_DIR "SFT-test[5]*", NULL ), &status );

  /* load once as a single SFT-vector (mix of detectors) */
  SHOULD_WORK ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );

  /* load once as a multi-SFT vector */
  SHOULD_WORK ( LALLoadMultiSFTs ( &status, &multsft_vect, catalog, -1, -1 ), &status );
  /* load again, using XLAL API */
  if ( ( multsft_vect2 = XLALLoadMultiSFTs ( catalog, -1, -1 )) == NULL ) {
    XLALPrintError ("%s: XLALLoadMultiSFTs (cat, -1, -1) failed with xlalErrno = %d\n", fn, xlalErrno );
    return SFTFILEIOTESTC_ESUB;
  }
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* 6 SFTs from 2 IFOs should have been read */
  if ( (sft_vect->length != 4) 	/* either as a single SFTVector */
       || (multsft_vect->length != 2) 	/* or separated by detector */
       || (multsft_vect->data[0]->length != 3) || ( multsft_vect->data[1]->length != 1 ) )
    {
      XLALPrintError ( "\nFailed to read in multi-SFT from 2 IFOs 'SFT-test*'!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* compare results from LALLoadMultiSFTs() and XLALLoadMultiSFTs() */
  {
    UINT4 numIFOs = multsft_vect->length;
    UINT4 X;
    for ( X=0; X < numIFOs; X ++ )
      {
        if( CompareSFTVectors ( multsft_vect->data[X], multsft_vect2->data[X] ) ) {
          XLALPrintError ("%s: comparing (X)LALLoadMultiSFTs(): sft-vectors differ for X=%d\n", fn, X );
          return SFTFILEIOTESTC_ESUB;
        }
      } /* for X < numIFOs */
  } /* ------ */

  /* ----- v2 SFT writing ----- */
  /* write v2-SFT to disk */
  SHOULD_WORK ( LALWriteSFT2file( &status, &(multsft_vect->data[0]->data[0]), "outputsftv2_v2.sft", "A v2-SFT file for testing!"), &status );

  SHOULD_WORK ( LALWriteSFTVector2Dir( &status, multsft_vect->data[0], ".", "A v2-SFT file for testing!", "test"), &status);

  /* write v2-SFT to single file */
  {
    const CHAR *currSingleSFT = NULL;
    UINT4 i = 0;
    FILE *fpConcat = NULL, *fpSingle = NULL;
    int concat = 0, single = 0;

    xlalErrno = 0;
    if (XLAL_SUCCESS != XLALWriteSFTVector2File ( multsft_vect->data[0], ".", "A v2-SFT file for testing!", "test_concat" )) {
      LALPrintError ( "\n XLALWriteSFTVector2File failed to write multi-SFT vector to file!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }
    /* check that the single file SFT is the same as the single SFTs */
    const UINT4 numSingleSFTs = 3;
    const CHAR *singleSFTs[] = {
      "H-1_H1_60SFT_test-000012345-61.sft",
      "H-1_H1_60SFT_test-000012465-61.sft",
      "H-1_H1_60SFT_test-000012585-61.sft"
    };
    printf("*** Comparing single and concatenated SFTs ***\n");
    /* try to open concatenated SFT */
    const CHAR *concatSFT = "H-3_H1_60SFT_test_concat-000012345-302.sft";
    if ( ( fpConcat = fopen(concatSFT, "rb" ) ) == NULL ) {
      LALPrintError ( "\n Cound not open SFT '%s'!\n\n", concatSFT);
      return SFTFILEIOTESTC_ESUB;
    }
    /* do loop while concat. SFT has data */
    while (!feof(fpConcat)) {
      /* get character from concat. SFT */
      concat = fgetc(fpConcat);
      if ( ferror(fpConcat) ) {
	LALPrintError ( "\n IO error reading '%s'!\n\n", concatSFT);
	return SFTFILEIOTESTC_ESUB;
      }
      /* get character from single SFT */
      while (1) {
	/* need to open next single SFT file */
	if (fpSingle == NULL) {
	  /* break if we've run out of single SFTs */
	  if (i == numSingleSFTs)
	    break;
	  /* try to open single SFT */
	  if ( ( fpSingle = fopen(singleSFTs[i], "rb" ) ) == NULL ) {
	    LALPrintError ( "\n Cound not open SFT '%s'!\n\n", singleSFTs[i]);
	    return SFTFILEIOTESTC_ESUB;
	  }
	  currSingleSFT = singleSFTs[i];
	}
	/* get character from single SFT */
	single = fgetc(fpSingle);
	if ( ferror(fpSingle) ) {
	  LALPrintError ( "\n IO error reading '%s'!\n\n", singleSFTs[i]);
	  return SFTFILEIOTESTC_ESUB;
	}
	/* if single SFT is out of data, close it (open next one at beginning of loop) */
	if (feof(fpSingle)) {
	  fclose(fpSingle);
	  fpSingle = NULL;
	  ++i;
	}
	/* otherwise we have a valid character */
	else
	  break;
      }
      /* do character-by-character comparison */
      if ( concat != single ) {
	LALPrintError ( "\n Comparison failed between '%s'(last char = %i) and '%s'(last char = %i)!!\n\n",
			concatSFT, concat, currSingleSFT, single );
	return SFTFILEIOTESTC_ESFTDIFF;
      }
    }
    fclose(fpConcat);
    printf( "*** Comparing was successful!!! ***\n");
  }

  /* write v2-SFt as a v1-SFT to disk (correct normalization) */
  multsft_vect->data[0]->data[0].epoch.gpsSeconds += 60;	/* shift start-time so they don't look like segmented SFTs! */
  SHOULD_WORK ( LALWrite_v2SFT_to_v1file( &status, &(multsft_vect->data[0]->data[0]), "outputsftv2_v1.sft"), &status );

  SUB ( LALDestroySFTVector ( &status, &sft_vect ), &status );
  SUB ( LALDestroyMultiSFTVector (&status, &multsft_vect ), &status );
  SUB ( LALDestroyMultiSFTVector (&status, &multsft_vect2 ), &status );

  /* ----- read the previous two SFTs back */
  SHOULD_FAIL ( LALSFTdataFind ( &status, &catalog, "outputsftv2_*.sft", NULL ), &status );
  /* need to set proper detector! */
  constraints.detector = detector;
  SUB ( LALSFTdataFind ( &status, &catalog, "outputsftv2_*.sft", &constraints ), &status);
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, -1, -1 ), &status );

  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) XLALPrintError ("\nFailed to read back in 'outputsftv2_*.sft'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  sft_vect2 = XLALLoadSFTs ( catalog, -1, -1 );
  if (!sft_vect2)
    {
      XLALPrintError ( "\nXLALLoadSFTs() call failed (where it should have succeeded)!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* compare the SFT vectors just read */
  if(CompareSFTVectors(sft_vect, sft_vect2))
    return SFTFILEIOTESTC_ESUB;

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

	if ( (crealf(*data1) != crealf(*data2)) || (cimagf(*data1) != cimagf(*data2)) )
	  {
	    XLALPrintError ("\nv1- and v2- SFT differ after writing/reading\n\n");
	    return SFTFILEIOTESTC_ESFTDIFF;
	  }
      } /* for i < numBins */
  }
  SUB ( LALDestroySFTVector (&status, &sft_vect2 ), &status );
  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* `----- v1 SFT writing */

  /* read v1-SFTs: 'inputsft.0' and 'inputsft.1' (one is big-endian, the other little-endian!) */
  SUB ( LALSFTdataFind (&status, &catalog, TEST_DATA_DIR "inputsft.?", &constraints ), &status );
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, fMin, fMax ), &status );
  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) XLALPrintError ("\nFailed to read in v1-SFTs 'inputsft.0' and 'inputsft.1'\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* read with XLALLoadSFTs() */
  sft_vect2 = XLALLoadSFTs ( catalog, fMin, fMax );
  if (!sft_vect2)
    {
      XLALPrintError ( "\nXLALLoadSFTs() call failed (where it should have succeeded)!\n\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* compare the SFT vectors just read */
  if(CompareSFTVectors(sft_vect, sft_vect2))
    return SFTFILEIOTESTC_ESUB;

  /* write v1-SFT to disk */
  SUB ( LALWriteSFTfile (&status, &(sft_vect->data[0]), "outputsft_v1.sft"), &status);

  /* try to write this v1-SFTs as v2: should fail without detector-info ! */
  strncpy( sft_vect->data[0].name, "??", 2 );
  SHOULD_FAIL (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );

  /* put detector there */
  strcpy ( sft_vect->data[0].name, "H1" );
  SHOULD_WORK (LALWriteSFT2file( &status, &(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!"), &status );

  SUB ( LALDestroySFTVector (&status, &sft_vect2 ), &status );
  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  /* ---------- test timestamps-reading functions by comparing LAL- and XLAL-versions against each other ---------- */
  {
#define TS_FNAME "testTimestamps.dat"
    LIGOTimeGPSVector *ts1 = NULL, *ts2 = NULL;

    /* ----- load timestamps with deprecated LAL function  */
    SUB ( LALReadTimestampsFile ( &status, &ts1, TEST_DATA_DIR TS_FNAME ), &status );
    /* ----- load timestamps w new XLAL function */
    if ( (ts2 = XLALReadTimestampsFile ( TEST_DATA_DIR TS_FNAME )) == NULL ) {
      XLALPrintError ("XLALReadTimestampsFile() failed to read timestamps from file '%s'. xlalErrno = %d\n", TS_FNAME );
      return SFTFILEIOTESTC_ESUB;
    }
    /* ----- compare the two */
    if ( ts1->length != ts2->length ) {
      XLALPrintError ("Read timestamps-lists differ in length %d != %d\n", ts1->length, ts2->length );
      return 1;
    }
    if ( ts1->deltaT != ts2->deltaT ) {
      XLALPrintError ("Read timestamps-lists differ in deltaT %g != %g\n", ts1->deltaT, ts2->deltaT );
      return 1;
    }
    UINT4 i, numTS = ts1->length;
    for ( i = 0; i < numTS; i ++ )
      {
        if ( XLALGPSDiff( &ts1->data[i], &ts2->data[i]) != 0 ) {
          XLALPrintError ("Read timestamps-lists differ in entry %d: { %d, %d } != { %d, %d }\n",
                          i + 1,
                          ts1->data[i].gpsSeconds, ts1->data[i].gpsNanoSeconds,
                          ts2->data[i].gpsSeconds, ts2->data[i].gpsNanoSeconds );
          return 1;
        }
      } /* for i < numTS */

    /* free mem */
    XLALDestroyTimestampVector ( ts1 );
    XLALDestroyTimestampVector ( ts2 );
  }

  /* ------------------------------ */
  LALCheckMemoryLeaks();

  XLALPrintError ("\n\n--------------------------------------------------------------------------------\n");
  XLALPrintError ("\n    OK. All tests passed correctly ! (error-messages above are OK!)\n");
  XLALPrintError ("\n--------------------------------------------------------------------------------\n");


  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}
/** \endcond */
