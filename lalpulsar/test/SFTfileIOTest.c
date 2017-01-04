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

#include <lal/LALStdio.h>
#include <lal/SFTfileIO.h>
#include <lal/SFTutils.h>
#include <lal/Units.h>

/*---------- DEFINES ----------*/
/**
 * \file
 * \ingroup SFTfileIO_h
 * \author R. Prix, B. Machenschalk, A.M. Sintes
 *
 * \brief Test-code for SFT-fileIO library
 *
 */

/** \cond DONT_DOXYGEN */

#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

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
    if ((sft1.epoch.gpsSeconds != sft2.epoch.gpsSeconds) ||
	(sft1.epoch.gpsNanoSeconds != sft2.epoch.gpsNanoSeconds)) {
      XLALPrintError ( "CompareSFTVectors(): SFT#%u epochs differ (%f/%f)!\n",
		       sft, GPS2REAL8(sft1.epoch), GPS2REAL8(sft2.epoch) );
      return(-1);
    }
    if ( strncmp(sft1.name,sft2.name, sizeof(sft1.name)) ) {
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

int main( void )
{
  const char *fn = __func__;

  SFTCatalog *catalog = NULL;
  SFTConstraints XLAL_INIT_DECL(constraints);
  SFTVector *sft_vect = NULL;
  SFTVector *sft_vect2 = NULL;
  MultiSFTVector *multsft_vect = NULL;
  MultiSFTVector *multsft_vect2 = NULL;
  CHAR detector[2] = "H1";
  BOOLEAN crc_check;

  /* band to read from infile.* SFTs */
  REAL8 fMin = 1008.5;
  REAL8 fMax = 1009.1;


  /* check that mal-formated SFTs are properly detected */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad1", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad2", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad3", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad4", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad5", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();

  /* the following (SFT-bad6) has a wrong CRC64 checksum. However, this is
   * not checked in XLALSFTdataFind, so it should succeed! */
  XLAL_CHECK_MAIN( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad6", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);

  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad7", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad8", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad9", NULL ) ) == NULL, XLAL_EFUNC); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad10", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad11", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad12", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad13", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad14", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();

  /* now check some crc-checksums */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test1", NULL ) ) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALCheckCRCSFTCatalog (&crc_check, catalog ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);
  if ( !crc_check )
    {
      XLALPrintError ("\nLALCheckSFTs(): SFT-test1 has correct checksum but LALCheckSFTs claimed it hasn't.\n\n");
      return EXIT_FAILURE;
    }
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-bad6", NULL ) ) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALCheckCRCSFTCatalog (&crc_check, catalog ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);
  if ( crc_check )
    {
      XLALPrintError ( "\nLALCheckSFTs() failed to catch invalid CRC checksum in SFT-bad6 \n\n");
      return EXIT_FAILURE;
    }

  /* check that proper v2-SFTs are read-in properly */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test1", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test2", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test3", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test4", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test5", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test6", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test7", NULL ) ) != NULL, XLAL_EFUNC ); XLALClearErrno();
  XLALDestroySFTCatalog(catalog);

  /* now completely read-in a v2 merged-SFT */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test*", NULL ) ) == NULL, XLAL_EFUNC ); XLALClearErrno();
  /* skip sft nr 4 with has Tsft=50 instead of Tsft=60 */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test[123567]*", NULL ) ) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);
  /* try the same with a ";" separated list of files and of patterns */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind (
				 TEST_DATA_DIR "SFT-test1;"
				 TEST_DATA_DIR "SFT-test2;"
				 TEST_DATA_DIR "SFT-test3;"
				 TEST_DATA_DIR "SFT-test5;"
				 TEST_DATA_DIR "SFT-test6;"
				 TEST_DATA_DIR "SFT-test7", NULL ) ) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "SFT-test[123]*;" TEST_DATA_DIR "SFT-test[5]*", NULL ) ) != NULL, XLAL_EFUNC );

  /* load once as a single SFT-vector (mix of detectors) */
  XLAL_CHECK_MAIN ( ( sft_vect = XLALLoadSFTs ( catalog, -1, -1 ) ) != NULL, XLAL_EFUNC );

  /* load once as a multi-SFT vector */
  XLAL_CHECK_MAIN ( ( multsft_vect = XLALLoadMultiSFTs ( catalog, -1, -1 ) ) != NULL, XLAL_EFUNC );
  /* load again, using XLAL API */
  if ( ( multsft_vect2 = XLALLoadMultiSFTs ( catalog, -1, -1 )) == NULL ) {
    XLALPrintError ("%s: XLALLoadMultiSFTs (cat, -1, -1) failed with xlalErrno = %d\n", fn, xlalErrno );
    return EXIT_FAILURE;
  }
  XLALDestroySFTCatalog(catalog);

  /* 6 SFTs from 2 IFOs should have been read */
  if ( (sft_vect->length != 4) 	/* either as a single SFTVector */
       || (multsft_vect->length != 2) 	/* or separated by detector */
       || (multsft_vect->data[0]->length != 3) || ( multsft_vect->data[1]->length != 1 ) )
    {
      XLALPrintError ( "\nFailed to read in multi-SFT from 2 IFOs 'SFT-test*'!\n\n");
      return EXIT_FAILURE;
    }

  /* compare results from XLALLoadMultiSFTs() and XLALLoadMultiSFTs() */
  {
    UINT4 numIFOs = multsft_vect->length;
    UINT4 X;
    for ( X=0; X < numIFOs; X ++ )
      {
        if( CompareSFTVectors ( multsft_vect->data[X], multsft_vect2->data[X] ) ) {
          XLALPrintError ("%s: comparing (X)XLALLoadMultiSFTs(): sft-vectors differ for X=%d\n", fn, X );
          return EXIT_FAILURE;
        }
      } /* for X < numIFOs */
  } /* ------ */

  /* ----- v2 SFT writing ----- */
  /* write v2-SFT to disk */
  XLAL_CHECK_MAIN ( XLALWriteSFT2file(&(multsft_vect->data[0]->data[0]), "outputsftv2_r1.sft", "A v2-SFT file for testing!") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( XLALWriteSFTVector2Dir(multsft_vect->data[0], ".", "A v2-SFT file for testing!", "test") == XLAL_SUCCESS, XLAL_EFUNC);

  /* write v2-SFT to single file */
  {
    const CHAR *currSingleSFT = NULL;
    UINT4 i = 0;
    FILE *fpConcat = NULL, *fpSingle = NULL;
    int concat = 0, single = 0;

    xlalErrno = 0;
    if (XLAL_SUCCESS != XLALWriteSFTVector2File ( multsft_vect->data[0], ".", "A v2-SFT file for testing!", "test_concat" )) {
      LALPrintError ( "\n XLALWriteSFTVector2File failed to write multi-SFT vector to file!\n\n");
      return EXIT_FAILURE;
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
      return EXIT_FAILURE;
    }
    /* do loop while concat. SFT has data */
    while (!feof(fpConcat)) {
      /* get character from concat. SFT */
      concat = fgetc(fpConcat);
      if ( ferror(fpConcat) ) {
	LALPrintError ( "\n IO error reading '%s'!\n\n", concatSFT);
	return EXIT_FAILURE;
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
	    return EXIT_FAILURE;
	  }
	  currSingleSFT = singleSFTs[i];
	}
	/* get character from single SFT */
	single = fgetc(fpSingle);
	if ( ferror(fpSingle) ) {
	  LALPrintError ( "\n IO error reading '%s'!\n\n", singleSFTs[i]);
	  return EXIT_FAILURE;
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
	return EXIT_FAILURE;
      }
    }
    fclose(fpConcat);
    printf( "*** Comparing was successful!!! ***\n");
  }

  /* write v2-SFT again */
  multsft_vect->data[0]->data[0].epoch.gpsSeconds += 60;       /* shift start-time so they don't look like segmented SFTs! */
  XLAL_CHECK_MAIN ( XLALWriteSFT2file(&(multsft_vect->data[0]->data[0]), "outputsftv2_r2.sft", "A v2-SFT file for testing!") == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroySFTVector ( sft_vect );
  sft_vect = NULL;
  XLALDestroyMultiSFTVector ( multsft_vect );
  multsft_vect = NULL;
  XLALDestroyMultiSFTVector ( multsft_vect2 );
  multsft_vect2 = NULL;

  /* ----- read the previous SFTs back */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( "outputsftv2_r*.sft", NULL ) ) != NULL, XLAL_EFUNC );
  XLALDestroySFTCatalog(catalog);
  constraints.detector = detector;
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( "outputsftv2_r*.sft", &constraints ) ) != NULL, XLAL_EFUNC);
  XLAL_CHECK_MAIN ( ( sft_vect = XLALLoadSFTs ( catalog, -1, -1 ) ) != NULL, XLAL_EFUNC );

  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) XLALPrintError ("\nFailed to read back in 'outputsftv2_r*.sft'\n\n");
      return EXIT_FAILURE;
    }

  sft_vect2 = XLALLoadSFTs ( catalog, -1, -1 );
  if (!sft_vect2)
    {
      XLALPrintError ( "\nXLALLoadSFTs() call failed (where it should have succeeded)!\n\n");
      return EXIT_FAILURE;
    }

  /* compare the SFT vectors just read */
  if(CompareSFTVectors(sft_vect, sft_vect2))
    return EXIT_FAILURE;

  XLALDestroySFTVector ( sft_vect2 );
  sft_vect2 = NULL;
  XLALDestroySFTVector ( sft_vect );
  sft_vect = NULL;
  XLALDestroySFTCatalog(catalog);

  /* read v1-SFTs: 'inputsft.0' and 'inputsft.1' (one is big-endian, the other little-endian!) */
  XLAL_CHECK_MAIN ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "inputsft.?", &constraints ) ) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( ( sft_vect = XLALLoadSFTs ( catalog, fMin, fMax ) ) != NULL, XLAL_EFUNC );
  if ( sft_vect->length != 2 )
    {
      if ( lalDebugLevel ) XLALPrintError ("\nFailed to read in v1-SFTs 'inputsft.0' and 'inputsft.1'\n\n");
      return EXIT_FAILURE;
    }

  /* read with XLALLoadSFTs() */
  sft_vect2 = XLALLoadSFTs ( catalog, fMin, fMax );
  if (!sft_vect2)
    {
      XLALPrintError ( "\nXLALLoadSFTs() call failed (where it should have succeeded)!\n\n");
      return EXIT_FAILURE;
    }

  /* compare the SFT vectors just read */
  if(CompareSFTVectors(sft_vect, sft_vect2))
    return EXIT_FAILURE;

  /* try to write this v1-SFTs as v2: should fail without detector-info ! */
  strncpy( sft_vect->data[0].name, "??", 2 );
  XLAL_CHECK_MAIN ( XLALWriteSFT2file(&(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!") != XLAL_SUCCESS, XLAL_EFUNC ); XLALClearErrno();

  /* put detector there */
  strcpy ( sft_vect->data[0].name, "H1" );
  XLAL_CHECK_MAIN ( XLALWriteSFT2file(&(sft_vect->data[0]), "outputsft_v2.sft", "Another v2-SFT file for testing!") == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroySFTVector ( sft_vect2 );
  sft_vect2 = NULL;
  XLALDestroySFTVector ( sft_vect );
  sft_vect = NULL;
  XLALDestroySFTCatalog(catalog);

  /* ---------- test timestamps-reading functions by comparing LAL- and XLAL-versions against each other ---------- */
  {
#define TS_FNAME "testTimestamps.dat"
#define TS_FNAME_NEW "testTimestampsNew.dat"
    LIGOTimeGPSVector *ts2 = NULL, *ts3 = NULL;

    /* ----- load timestamps w new XLAL function */
    XLAL_CHECK_MAIN ( (ts2 = XLALReadTimestampsFile ( TEST_DATA_DIR TS_FNAME )) != NULL, XLAL_EFUNC );
    XLAL_CHECK_MAIN ( (ts3 = XLALReadTimestampsFile ( TEST_DATA_DIR TS_FNAME_NEW )) != NULL, XLAL_EFUNC );

    /* ----- compare the 3 */
    XLAL_CHECK_MAIN ( ts2->length == ts3->length, XLAL_EFAILED, "Read timestamps-lists differ in length %d != %d\n", ts2->length, ts3->length );

    XLAL_CHECK_MAIN ( ts2->deltaT == ts3->deltaT, XLAL_EFAILED, "Read timestamps-lists differ in deltaT %g != %g\n", ts2->deltaT, ts3->deltaT );

    UINT4 numTS = ts2->length;
    char buf1[256], buf2[256];
    for ( UINT4 i = 0; i < numTS; i ++ )
      {
        XLAL_CHECK_MAIN ( XLALGPSDiff( &ts2->data[i], &ts3->data[i]) == 0, XLAL_EFAILED,
                          "Timestamps-lists differ in entry %" LAL_UINT4_FORMAT ": %s != %s\n", i + 1, XLALGPSToStr ( buf1, &ts2->data[i] ), XLALGPSToStr ( buf2, &ts3->data[i] ) );
      } /* for i < numTS */

    /* free mem */
    XLALDestroyTimestampVector ( ts2 );
    XLALDestroyTimestampVector ( ts3 );
  }

  /* ------------------------------ */
  LALCheckMemoryLeaks();

  XLALPrintError ("\n\n--------------------------------------------------------------------------------\n");
  XLALPrintError ("\n    OK. All tests passed correctly ! (error-messages above are OK!)\n");
  XLALPrintError ("\n--------------------------------------------------------------------------------\n");


  return EXIT_SUCCESS;
}
/** \endcond */
