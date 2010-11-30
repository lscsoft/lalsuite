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
#include <stdlib.h>

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


/* Default parameters. */

/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )					\
  do {									\
  if ( lalDebugLevel & LALERROR )					\
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"	\
		    "        %s %s\n", (code), *argv, __FILE__,		\
		    __LINE__, SFTFILEIOTESTC, statement ? statement :	\
		    "", (msg) );					\
} while (0)

#define INFO( statement )					      \
  do {								      \
    if ( lalDebugLevel & LALINFO )				      \
      XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"   \
		      "        %s\n", *argv, __FILE__, __LINE__,      \
		      SFTFILEIOTESTC, (statement) );		      \
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
  xlalErrno = 0;							                \
  if ( func, ! (statusptr)->statusCode ) {						\
    ERROR( SFTFILEIOTESTC_ESUB, SFTFILEIOTESTC_MSGESUB,      				\
          "Function call '" #func "' should have failed for this SFT but didn't!\n");	\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)

#define SHOULD_FAIL_WITH_CODE( func, statusptr, code )					\
do { 											\
  xlalErrno = 0;							                \
  if ( func, (statusptr)->statusCode != code) {						\
    XLALPrintError( "Function call '" #func "' should have failed with code " #code ", but returned %d instead.\n",	\
		   (statusptr)->statusCode );						\
    return SFTFILEIOTESTC_ESUB;   			                               	\
   }											\
} while(0)


#define SHOULD_WORK( func, statusptr )							\
do { 											\
  xlalErrno = 0;							                \
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

/* ----------------------------------------------------------------------*/

static int CompareSFTVectors(SFTVector *sft_vect, SFTVector *sft_vect2);
static int CompareSFTVectors(SFTVector *sft_vect, SFTVector *sft_vect2)
{
  UINT4 sft,bin;
  if (sft_vect->length != sft_vect2->length) {
    XLALPrintError ( "\nCompareSFTVectors(): vector lengths differ!\n\n");
    return(-1);
  }
  for(sft=0; sft < sft_vect->length; sft++) {
    SFTtype sft1 = sft_vect->data[sft];
    SFTtype sft2 = sft_vect2->data[sft];
    if (sft1.data->length != sft2.data->length) {
      XLALPrintError ( "\nCompareSFTVectors(): lengths of SFT#%u differ!\n\n", sft);
      return(-1);
    }
    for(bin=0; bin < sft1.data->length; bin++) {
      if((sft1.data->data[bin].re != sft2.data->data[bin].re) ||
	 (sft1.data->data[bin].im != sft2.data->data[bin].im)) {
	XLALPrintError ( "\nCompareSFTVectors(): bins %u of SFT#%u differ!\n\n", sft, bin);
	return(-1);
      }
    }
  }
  return(0);
}

int main(int argc, char *argv[])
{
  LALStatus status = empty_status;
  SFTCatalog *catalog = NULL;
  SFTConstraints constraints = empty_constraints;
  SFTVector *sft_vect = NULL;
  SFTVector *sft_vect2 = NULL;
  CHAR detector[2] = "H1";
  INT4 crc_check;
  /* band to read from infile.* SFTs */
  REAL8 fMin = -1;
  REAL8 fMax = -1;

  lalDebugLevel = 3;

  if(argc != 4) {
    XLALPrintError ( "Usage: %s <files> <fmin> <fmax>\n", argv[0]);
      return SFTFILEIOTESTC_EARG;
  } else {
    fMin = atof(argv[2]);
    fMax = atof(argv[3]);
  }

  SUB ( LALSFTdataFind ( &status, &catalog, argv[1], &constraints ), &status);
  SUB ( LALLoadSFTs ( &status, &sft_vect, catalog, fMin, fMax ), &status );

  sft_vect2 = XLALLoadSFTs ( catalog, fMin, fMax );
  if (!sft_vect2)
    {
      XLALPrintError ( "ERROR: XLALLoadSFTs() call failed!\n");
      return SFTFILEIOTESTC_ESUB;
    }

  /* compare the SFT vectors just read */
  if(CompareSFTVectors(sft_vect, sft_vect2))
    return SFTFILEIOTESTC_ESUB;

  SUB ( LALDestroySFTVector (&status, &sft_vect2 ), &status );
  SUB ( LALDestroySFTVector (&status, &sft_vect ), &status );
  SUB ( LALDestroySFTCatalog( &status, &catalog), &status );

  LALCheckMemoryLeaks();

  INFO( SFTFILEIOTESTC_MSGENORM );
  return SFTFILEIOTESTC_ENORM;
}
