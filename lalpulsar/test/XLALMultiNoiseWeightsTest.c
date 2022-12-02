/*
 * Copyright (C) 2013 John T. Whelan
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <config.h>
#include <math.h>
#include <sys/times.h>

#include <lal/NormalizeSFTRngMed.h>
#include <lal/SFTfileIO.h>

/**
 * \author John T. Whelan
 * \file
 * \ingroup SFTfileIO_h
 * \brief Tests for XLALComputeMultiNoiseWeights()
 *
 * PSDs are calculated using the test SFTs created for
 * SFTfileIOTest.c
 *
 */

/*---------- macros ---------- */
#define FRACERR(x,y) (fabs((x)-(y))/(0.5*((x)+(y))))

/* ----- internal prototypes ---------- */
int XLALCompareMultiNoiseWeights ( MultiNoiseWeights *multiWeights1, MultiNoiseWeights *multiWeights2, REAL8 tolerance );

/* ----- function definitions ---------- */
int
main (  void )
{
  LALStatus XLAL_INIT_DECL(status);
  SFTCatalog *catalog = NULL;
  SFTConstraints XLAL_INIT_DECL(constraints);
  MultiSFTVector *multiSFTs = NULL;
  MultiPSDVector *multiPSDs = NULL;
  MultiNoiseWeights *multiWeightsXLAL = NULL;
  MultiNoiseWeights *multiWeightsCorrect = NULL;
  UINT4 rngmedBins = 11;
  REAL8 tolerance = 2e-6;	/* same algorithm, should be basically identical results */

  /* Construct the "correct" weights, calculated using the old LAL routines */
  UINT4 numIFOsCorrect = 2;
  XLAL_CHECK ( ( multiWeightsCorrect = XLALCalloc ( 1, sizeof(*multiWeightsCorrect ) ) ) != NULL, XLAL_ENOMEM );
  multiWeightsCorrect->length = numIFOsCorrect;
  multiWeightsCorrect->Sinv_Tsft = 1.980867126449e+52;
  XLAL_CHECK ( ( multiWeightsCorrect->data = XLALCalloc ( numIFOsCorrect, sizeof(*multiWeightsCorrect->data ) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK ( ( multiWeightsCorrect->data[0] = XLALCreateREAL8Vector(4) ) != NULL, XLAL_ENOMEM );
  multiWeightsCorrect->data[0]->data[0] = 6.425160659487e-05;
  multiWeightsCorrect->data[0]->data[1] = 7.259453662367e-06;
  multiWeightsCorrect->data[0]->data[2] = 9.838893684664e-04;
  multiWeightsCorrect->data[0]->data[3] = 5.043766789923e-05;
  XLAL_CHECK ( ( multiWeightsCorrect->data[1] = XLALCreateREAL8Vector(3) ) != NULL, XLAL_ENOMEM );
  multiWeightsCorrect->data[1]->data[0] = 1.582309910283e-04;
  multiWeightsCorrect->data[1]->data[1] = 5.345673753744e-04;
  multiWeightsCorrect->data[1]->data[2] = 6.998201363537e+00;

  /* Construct the catalog */
  XLAL_CHECK ( ( catalog = XLALSFTdataFind ( TEST_DATA_DIR "MultiNoiseWeightsTest*.sft", &constraints ) ) != NULL, XLAL_EFUNC, " XLALSFTdataFind failed\n" );

  /* Load the SFTs */
  XLAL_CHECK ( ( multiSFTs = XLALLoadMultiSFTs ( catalog, -1, -1 ) ) != NULL, XLAL_EFUNC, " XLALLoadMultiSFTs failed\n" );

  /* calculate the psd and normalize the SFTs */
  XLAL_CHECK ( ( multiPSDs = XLALNormalizeMultiSFTVect ( multiSFTs, rngmedBins, NULL ) ) != NULL, XLAL_EFUNC, " XLALNormalizeMultiSFTVect failed\n" );

  /* Get weights using XLAL function */
  XLAL_CHECK ( ( multiWeightsXLAL = XLALComputeMultiNoiseWeights ( multiPSDs, rngmedBins, 0 ) ) != NULL, XLAL_EFUNC, " XLALComputeMultiNoiseWeights failed\n" );

  /* Compare XLAL weights to reference */
  XLAL_CHECK ( XLALCompareMultiNoiseWeights ( multiWeightsXLAL, multiWeightsCorrect, tolerance ) == XLAL_SUCCESS, XLAL_EFAILED, "Comparison between XLAL and reference MultiNoiseWeights failed\n" );

  /* Also check copying some weights */
  MultiNoiseWeights *copiedMultiWeights = NULL;
  XLAL_CHECK ( ( copiedMultiWeights = XLALCopyMultiNoiseWeights ( multiWeightsCorrect ) ) != NULL, XLAL_EFUNC, "XLALCopyMultiNoiseWeights failed\n" );
  XLAL_CHECK ( XLALCompareMultiNoiseWeights ( copiedMultiWeights, multiWeightsCorrect, tolerance ) == XLAL_SUCCESS, XLAL_EFAILED, "Comparison between reference MultiNoiseWeights and a copy failed\n" );

  /* Clean up memory */
  XLALDestroyMultiNoiseWeights ( multiWeightsCorrect );
  XLALDestroyMultiNoiseWeights ( multiWeightsXLAL );
  XLALDestroyMultiNoiseWeights ( copiedMultiWeights );
  XLALDestroyMultiPSDVector ( multiPSDs );
  XLALDestroyMultiSFTVector ( multiSFTs );
  XLALDestroySFTCatalog ( catalog );
  /* check for memory-leaks */
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} /* main() */

/**
 * Comparison function for two multiNoiseWeights vectors, return success or failure for given tolerance.
 * The fractional error is checked for the weights and normalization factors.
 *
 */
int
XLALCompareMultiNoiseWeights ( MultiNoiseWeights *multiWeights1, MultiNoiseWeights *multiWeights2, REAL8 tolerance )
{

  XLALPrintInfo("Sinv_Tsft1 %.12e; Sinv_Tsft2 %.12e\n", multiWeights1->Sinv_Tsft, multiWeights2->Sinv_Tsft);
  XLAL_CHECK ( FRACERR( multiWeights1->Sinv_Tsft, multiWeights2->Sinv_Tsft ) <= tolerance, XLAL_EFAILED, "%s: Sinv_Tsft differs by more than tolerance %g multiWeights1 = %g, multiWeights2 = %g\n", __func__, tolerance, multiWeights1->Sinv_Tsft, multiWeights2->Sinv_Tsft );

  XLAL_CHECK( multiWeights1->length == multiWeights2->length, XLAL_EFAILED, "%s: numbers of detectors differ multiWeights1 = %d, multiWeights2 = %d\n", __func__, multiWeights1->length, multiWeights2->length );
  UINT4 numIFOs = multiWeights1->length;
  for ( UINT4 X = 0; X < numIFOs; X++)
    {
      REAL8Vector *weights1 = multiWeights1->data[X];
      REAL8Vector *weights2 = multiWeights2->data[X];

      XLAL_CHECK( weights1->length == weights2->length, XLAL_EFAILED, "%s: numbers of SFTs for detector %d differ multiWeights1 = %d, multiWeights2 = %d\n", __func__, X, weights1->length, weights2->length );
      UINT4 numSFTs = weights1->length;

      for ( UINT4 alpha = 0; alpha < numSFTs; alpha++)
	{
	  XLALPrintInfo("IFO %d; SFT %d; weight1 %.12e; weight2 %.12e\n", X, alpha, weights1->data[alpha], weights2->data[alpha]);
	  XLAL_CHECK ( FRACERR( weights1->data[alpha], weights2->data[alpha] ) <= tolerance, XLAL_EFAILED, "%s: weights for IFO %d, SFT %d differ by more than tolerance %g multiWeights1 = %g, multiWeights2 = %g\n", __func__, X, alpha, tolerance, weights1->data[alpha], weights2->data[alpha] );
	}
    }

  return XLAL_SUCCESS;

} /* XLALCompareMultiNoiseWeights() */
