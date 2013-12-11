/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

/**
 * \author UTB Relativity Group; contact whelan\@phys.utb.edu
 * \file
 * \ingroup StochasticInverseNoise_c
 *
 * \brief Test suite for <tt>LALStochasticInverseNoise()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./StochasticInverseNoiseTest [options]
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -n length      frequency series contain length points
 * -w filename    read uncalibrated noise PSD from file filename
 * -f filename    read response function from file filename
 * -u filename    print calibrated inverse noise PSD to file filename
 * -m filename    print half-calibrated inverse noise PSD to file filename
 * \endcode
 *
 * ### Description ###
 *
 * This program tests the function <tt>LALStochasticInverseNoise()</tt>,
 * which outputs an uncalibrated and "half-calibrated" inverse noise spectra
 * from a uncalibrated data stream and a response function.
 *
 * First, it tests that the correct error codes
 * (cf. \ref StochasticCrossCorrelation_h)
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to output structure</em></li>
 * <li> <em>null pointer to input structure</em></li>
 * <li> <em>null pointer to uncalibrated noise</em></li>
 * <li> <em>null pointer to response function</em></li>
 * <li> <em>null pointer to calibrated inverse noise</em></li>
 * <li> <em>null pointer to half-calibrated inverse noise</em></li>
 * <li> <em>null pointer to data member of uncalibrated noise</em></li>
 * <li> <em>null pointer to data member of response function</em></li>
 * <li> <em>null pointer to data member of calibrated inverse noise</em></li>
 * <li> <em>null pointer to data member of half-calibrated inverse noise</em></li>
 * <li> <em>null pointer to data member of data member of uncalibrated noise</em></li>
 * <li> <em>null pointer to data member of data member of response function</em></li>
 * <li> <em>null pointer to data member of data member of calibrated inverse noise</em></li>
 * <li> <em>null pointer to data member of data member of half-calibrated inverse noise</em></li>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * <li> negative start frequency</li>
 * <li> length mismatch between uncalibrated noise and response function</li>
 * <li> length mismatch between uncalibrated noise and calibrated inverse noise</li>
 * <li> length mismatch between uncalibrated noise and half-calibrated inverse noise</li>
 * <li> frequency spacing mismatch between uncalibrated noise and response function</li>
 * <li> start frequency mismatch between uncalibrated noise and response function</li>
 * </ul>
 *
 * It then verifies that the correct uncalibrated and half-calibrated inverse
 * noise are generated for a simple test case:
 * <ol>
 * <li> \f$\tilde{R}(f)=(1+i)f^2\f$, \f$P(f)=f^3\f$.  The
 * expected results are \f$1/P^{\mathrm{C}}(f)=2f\f$,
 * \f$1/P^{\mathrm{HC}}(f)=(1-i)f^{-1}\f$.</li>
 * </ol>
 *
 * For each successful test (both of these valid data and the invalid ones
 * described above), it prints "\c PASS" to standard output; if a
 * test fails, it prints "\c FAIL".
 *
 * If the four \c filename arguments are present, it also
 * calculates a spectrum based on user-specified data and it prints the
 * noise spectra to the files specified by the user.
 *
 * ### Uses ###
 *
 * \code
 * getopt()
 * LALStochasticInverseNoise()
 * LALSCreateVector()
 * LALCCreateVector()
 * LALSDestroyVector()
 * LALCDestroyVector()
 * LALSReadFrequencySeries()
 * LALCReadFrequencySeries()
 * LALSPrintFrequencySeries()
 * LALCPrintFrequencySeries()
 * LALUnitAsString()
 * LALUnitCompare()
 * LALCheckMemoryLeaks()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 * <li> No specific error checking is done on user-specified data.  If
 * the \c length argument missing, the resulting defaults
 * will cause a bad data error.</li>
 * <li> The length of the user-provided series must be specified, even
 * though it could in principle be deduced from the input file, because
 * the data sequences must be allocated before the
 * <tt>LALCReadFrequencySeries()</tt> function is called.</li>
 * <li> If one \c filename argument, but not both, is present,
 * the user-specified data will be silently ignored.</li>
 * </ul>
 *
 */

/**\name Error Codes */ /*@{*/
#define STOCHASTICINVERSENOISETESTC_ENOM 0	/**< Nominal exit */
#define STOCHASTICINVERSENOISETESTC_EARG 1	/**< Error parsing command-line arguments */
#define STOCHASTICINVERSENOISETESTC_ECHK 2	/**< Error checking failed to catch bad data */
#define STOCHASTICINVERSENOISETESTC_EFLS 3	/**< Incorrect answer for valid data */
#define STOCHASTICINVERSENOISETESTC_EUSE 4	/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define STOCHASTICINVERSENOISETESTC_MSGENOM "Nominal exit"
#define STOCHASTICINVERSENOISETESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICINVERSENOISETESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICINVERSENOISETESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICINVERSENOISETESTC_MSGEUSE "Bad user-entered data"


#include <lal/LALStdlib.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>

#include "CheckStatus.h"

#define STOCHASTICINVERSENOISETESTC_TRUE     1
#define STOCHASTICINVERSENOISETESTC_FALSE    0
#define STOCHASTICINVERSENOISETESTC_DELTAF   1.0
#define STOCHASTICINVERSENOISETESTC_F0       0.0
#define STOCHASTICINVERSENOISETESTC_LENGTH   8
#define STOCHASTICINVERSENOISETESTC_TOL      1e-6

extern char *optarg;
extern int   optind;


BOOLEAN optVerbose = STOCHASTICINVERSENOISETESTC_FALSE;
UINT4 optLength    = 0;
CHAR  optWNoiseFile[LALNameLength]   = "";
CHAR  optWFilterFile[LALNameLength]     = "";
CHAR  optInvNoiseFile[LALNameLength] = "";
CHAR  optHWInvNoiseFile[LALNameLength] = "";

INT4 code;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

int main(int argc, char *argv[])
{

  static LALStatus                status;

  StochasticInverseNoiseInput        input;
  StochasticInverseNoiseOutput       output;

  LIGOTimeGPS              epoch = {1234,56789};

  REAL4FrequencySeries     wNoise;
  COMPLEX8FrequencySeries  wFilter;
  REAL4FrequencySeries     invNoise;
  COMPLEX8FrequencySeries  hwInvNoise;

  REAL8      f;
  INT4       i;
  REAL4      expectedReal;
  REAL4      expectedImag;

  LALUnitPair              unitPair;
  LALUnit                  expectedUnit;
  BOOLEAN                  result;

  CHARVector               *unitString = NULL;


  ParseOptions (argc, argv);

  /* define valid parameters */
  wNoise.f0     = STOCHASTICINVERSENOISETESTC_F0;
  wNoise.deltaF = STOCHASTICINVERSENOISETESTC_DELTAF;
  wNoise.epoch  = epoch;
  wNoise.data   = NULL;
  invNoise.data = NULL;
#ifndef LAL_NDEBUG
  REAL4FrequencySeries     realBadData = wNoise;
#endif

  wFilter.f0     = wNoise.f0;
  wFilter.deltaF  = wNoise.deltaF;
  wFilter.epoch   = wNoise.epoch;
  wFilter.data    = NULL;
  hwInvNoise.data = NULL;
#ifndef LAL_NDEBUG
  COMPLEX8FrequencySeries  complexBadData  = wFilter;
#endif

  /******** Set Testing  Units ********/
  /* response function */
  wFilter.sampleUnits = lalDimensionlessUnit;
  wFilter.sampleUnits.unitNumerator[LALUnitIndexADCCount] = 1;
  wFilter.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;

  /* uncalibrated noise */
  wNoise.sampleUnits = lalDimensionlessUnit;
  wNoise.sampleUnits.unitNumerator[LALUnitIndexADCCount] = 2;
  wNoise.sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
  /**************************************************/

  /* allocate memory */
  /* printf("About to create\n"); */
  LALSCreateVector(&status, &(wNoise.data),
		   STOCHASTICINVERSENOISETESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  /* printf("Just created\n"); */

  LALSCreateVector(&status, &(invNoise.data),
		   STOCHASTICINVERSENOISETESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(wFilter.data),
		   STOCHASTICINVERSENOISETESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(hwInvNoise.data),
		   STOCHASTICINVERSENOISETESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.unCalibratedNoisePSD = &wNoise;
  input.responseFunction = &wFilter;
  output.calibratedInverseNoisePSD = &invNoise;
  output.halfCalibratedInverseNoisePSD = &hwInvNoise;

 /* TEST INVALID DATA HERE -------------------------------------- */

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to input structure */
    LALStochasticInverseNoise(&status, &output, NULL);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to output structure */
    LALStochasticInverseNoise(&status, NULL, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to output structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to wNoise member of input structure */
    input.unCalibratedNoisePSD = NULL;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to uncalibrated noise results in error:\n       \"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unCalibratedNoisePSD = &wNoise;
    /* test behavior for null pointer to wFitler member of input structure */
    input.responseFunction = NULL;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to response function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.responseFunction = &wFilter;

    /* test behavior for null pointer to invNoise member of output structure */
    output.calibratedInverseNoisePSD = NULL;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to calibrated inverse noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.calibratedInverseNoisePSD = &invNoise;

    /* test behavior for null pointer to half-calibrated inverse noise member */
    /* of output structure */
    output.halfCalibratedInverseNoisePSD = NULL;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to half-calibrated inverse noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.halfCalibratedInverseNoisePSD = &hwInvNoise;

    /* test behavior for null pointer to data member of wnoise */
    input.unCalibratedNoisePSD = &realBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of uncalibrated noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unCalibratedNoisePSD = &wNoise;

    /* test behavior for null pointer to data member of wFilter */
    input.responseFunction = &complexBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of response function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.responseFunction = &wFilter;

    /* test behavior for null pointer to data member of invNoise  */
    output.calibratedInverseNoisePSD = &realBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of calibrated inverse noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.calibratedInverseNoisePSD = &invNoise;

    /* test behavior for null pointer to data member of half-calibrated */
    /* inverse noise */
    output.halfCalibratedInverseNoisePSD = &complexBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of half-calibrated noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.halfCalibratedInverseNoisePSD = &hwInvNoise;

    /* Create a vector for testing REAL4 null data-data pointers */
    LALSCreateVector(&status, &(realBadData.data), STOCHASTICINVERSENOISETESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    REAL4                   *sPtr;
    sPtr = realBadData.data->data;
    realBadData.data->data = NULL;

    /* test behavior for null pointer to data-data member of wNoise */
    input.unCalibratedNoisePSD = &realBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data-data member of uncalibrated noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.unCalibratedNoisePSD = &wNoise;

    /* test behavior for null pointer to data-data member of invNoise */
    output.calibratedInverseNoisePSD = &realBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data-data member of calibrated inverse noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.calibratedInverseNoisePSD = &invNoise;

    /* Create a vector for testing COMPLEX8 null data-data pointers */
    LALCCreateVector(&status, &(complexBadData.data), STOCHASTICINVERSENOISETESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    COMPLEX8                *cPtr;
    cPtr = complexBadData.data->data;
    complexBadData.data->data = NULL;

    /* test behavior for null pointer to data-data member of wFilter */
    input.responseFunction = &complexBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data-data member of response function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.responseFunction = &wFilter;

    /* test behavior for null pointer to data-data member of hwInvNoise */
    output.halfCalibratedInverseNoisePSD = &complexBadData;
    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data-data member of half-calibrated inverse noise results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    output.halfCalibratedInverseNoisePSD = &hwInvNoise;

    /** clean up **/
    realBadData.data->data = sPtr;
    LALSDestroyVector(&status, &(realBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }

    complexBadData.data->data = cPtr;
    LALCDestroyVector(&status, &(complexBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* test behavior for zero length */
    wNoise.data->length =
      wFilter.data->length =
      invNoise.data->length =
      hwInvNoise.data->length = 0;

    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
			      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    wNoise.data->length =
      wFilter.data->length =
      invNoise.data->length =
      hwInvNoise.data->length =STOCHASTICINVERSENOISETESTC_LENGTH;

    /* test behavior for negative frequency spacing */
    wNoise.deltaF =
      wFilter.deltaF = - STOCHASTICINVERSENOISETESTC_DELTAF;

    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    wNoise.deltaF =
      wFilter.deltaF = STOCHASTICINVERSENOISETESTC_DELTAF;

    /* test behavior for zero frequency spacing */
    wNoise.deltaF =
    wFilter.deltaF = 0.0;

    LALStochasticInverseNoise(&status, &output, &input);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICINVERSENOISETESTC_ECHK,
			      STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
	   STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    wNoise.deltaF =
      wFilter.deltaF = STOCHASTICINVERSENOISETESTC_DELTAF;
  } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  wNoise.f0 = wFilter.f0 = -3.0;

  LALStochasticInverseNoise(&status, &output, &input);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
			    STOCHASTICINVERSENOISETESTC_ECHK,
			    STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: negative start  frequency results in error:\n     \"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  /* reasign valid f0 */
  wNoise.f0 =
    wFilter.f0 = STOCHASTICINVERSENOISETESTC_F0;

  /* test behavior for length mismatch between wNoise and wFilter */
  wFilter.data->length
    = STOCHASTICINVERSENOISETESTC_LENGTH - 1;
  LALStochasticInverseNoise(&status, &output, &input);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICINVERSENOISETESTC_ECHK,
			    STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between uncalibrated noise and response function results in error:\n       \"%s\"\n",
	 STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  wFilter.data->length = STOCHASTICINVERSENOISETESTC_LENGTH;

  /* test behavior length mismatch between wNoise and invNoise */
  invNoise.data->length
    = STOCHASTICINVERSENOISETESTC_LENGTH - 1;
  LALStochasticInverseNoise(&status, &output, &input);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICINVERSENOISETESTC_ECHK,
			    STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
  {
    return code;
  }
   printf("  PASS: length mismatch between calibrated inverse noise and calibrated inverse noise results in error:\n       \"%s\"\n",
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
   invNoise.data->length = STOCHASTICINVERSENOISETESTC_LENGTH;

   /* test behavior length mismatch between wNoise and hwInvNoise */
   hwInvNoise.data->length
    = STOCHASTICINVERSENOISETESTC_LENGTH - 1;
   LALStochasticInverseNoise(&status, &output, &input);
   if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			     STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			     STOCHASTICINVERSENOISETESTC_ECHK,
			     STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
   {
       return code;
   }
   printf("  PASS: length mismatch between whtiened inverse noise and half-calibrated inverse noise results in error:\n       \"%s\"\n",
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
   hwInvNoise.data->length = STOCHASTICINVERSENOISETESTC_LENGTH;

   /* test behavior for initial frequency mismatch between wNoise and wFilter*/
   wFilter.f0 = 30;
   LALStochasticInverseNoise(&status, &output, &input);
   if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			     STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			     STOCHASTICINVERSENOISETESTC_ECHK,
			     STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
   {
       return code;
   }
   printf("  PASS: initial frequency mismatch between uncalibrated noise and response function results in error:\n       \"%s\"\n",
          STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
   wFilter.f0 = STOCHASTICINVERSENOISETESTC_F0;

   /* test behavior for frequency spacing mismatch between wNoise and wFilter*/
   wFilter.deltaF =
     2.0 * STOCHASTICINVERSENOISETESTC_DELTAF;
   LALStochasticInverseNoise(&status, &output, &input);
   if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			     STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			     STOCHASTICINVERSENOISETESTC_ECHK,
			     STOCHASTICINVERSENOISETESTC_MSGECHK) ) )
   {
       return code;
   }
   printf("  PASS: frequency spacing mismatch between uncalibrated noise and response function results in error:\n       \"%s\"\n",
          STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
   wFilter.deltaF = STOCHASTICINVERSENOISETESTC_DELTAF;


 /* VALID TEST DATA HERE ----------------------------------------- */

   /* create input to test */
   for (i=0; i < STOCHASTICINVERSENOISETESTC_LENGTH; i++)
   {
     f = i*STOCHASTICINVERSENOISETESTC_DELTAF;

     wNoise.data->data[i]     = f*f*f;
     wFilter.data->data[i] = crectf( f*f, f*f );
   }

   /* fill inverse noise input and output */
   input.unCalibratedNoisePSD      = &(wNoise);
   input.responseFunction                 = &(wFilter);
   output.calibratedInverseNoisePSD         = &(invNoise);
   output.halfCalibratedInverseNoisePSD       = &(hwInvNoise);

   /*calculate the calibrated inverse noise and half-calibrated inverse noise*/
   LALStochasticInverseNoise(&status, &output, &input );
   if ( ( code = CheckStatus (&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     printf("  Valid Data Test:\n");
     printf("  Checking half-calibrated inverse noise...\n");
   }
   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", hwInvNoise.f0,
	    STOCHASTICINVERSENOISETESTC_F0);
   }
   if ( fabs(hwInvNoise.f0-STOCHASTICINVERSENOISETESTC_F0)
	> STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", hwInvNoise.deltaF,
            STOCHASTICINVERSENOISETESTC_DELTAF);
   }
   if ( fabs(hwInvNoise.deltaF-STOCHASTICINVERSENOISETESTC_DELTAF)
        / STOCHASTICINVERSENOISETESTC_DELTAF
	> STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

   /* check output units */

   expectedUnit = lalDimensionlessUnit;
   expectedUnit.unitNumerator[LALUnitIndexADCCount] = -1;
   expectedUnit.unitNumerator[LALUnitIndexStrain] = -1;
   expectedUnit.unitNumerator[LALUnitIndexSecond] = -1;
   unitPair.unitOne = &expectedUnit;
   unitPair.unitTwo = &(hwInvNoise.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     STOCHASTICINVERSENOISETESTC_EFLS,
			     STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units of 1/PHW(f) are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
     {
       return code;
     }
   }

   if (!result)
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
	      STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("1/PHW(0)=%g + %g i, should be 0\n",
            crealf(hwInvNoise.data->data[0]), cimagf(hwInvNoise.data->data[0]));
   }
   if ( fabs(crealf(hwInvNoise.data->data[0])) > STOCHASTICINVERSENOISETESTC_TOL
        || fabs(cimagf(hwInvNoise.data->data[0]))
	> STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

  for (i=1; i < STOCHASTICINVERSENOISETESTC_LENGTH; i++)
  {
    f = i*STOCHASTICINVERSENOISETESTC_DELTAF;
    expectedReal = 1/f;
    expectedImag = - expectedReal;
    if (optVerbose)
    {
      printf("1/PHW(%f Hz)=%g + %g i, should be %g + %g i\n",
	     f, crealf(hwInvNoise.data->data[i]), cimagf(hwInvNoise.data->data[i]),
	     expectedReal, expectedImag);
    }
    if (fabs(crealf(hwInvNoise.data->data[i]) - expectedReal)/expectedReal
	> STOCHASTICINVERSENOISETESTC_TOL
	|| fabs(cimagf(hwInvNoise.data->data[i]) - expectedImag)/expectedImag
	> STOCHASTICINVERSENOISETESTC_TOL)
    {
      printf("  FAIL: Valid data test\n");
      if (optVerbose)
      {
	printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
      }
      return STOCHASTICINVERSENOISETESTC_EFLS;
    }
  }
  /****** check valid calibrated inverse noise *******/
   if (optVerbose)
   {
     printf("  Checking calibrated inverse noise...\n");
   }
   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", invNoise.f0,
	    STOCHASTICINVERSENOISETESTC_F0);
   }
   if ( fabs(invNoise.f0-STOCHASTICINVERSENOISETESTC_F0)
	> STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", invNoise.deltaF,
            STOCHASTICINVERSENOISETESTC_DELTAF);
   }
   if ( fabs(invNoise.deltaF-STOCHASTICINVERSENOISETESTC_DELTAF)
        / STOCHASTICINVERSENOISETESTC_DELTAF
	> STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

   /* check output units */
  expectedUnit = lalDimensionlessUnit;
  expectedUnit.unitNumerator[LALUnitIndexStrain] = -2;
  expectedUnit.unitNumerator[LALUnitIndexSecond] = -1;
  unitPair.unitOne = &expectedUnit;
  unitPair.unitTwo = &(invNoise.sampleUnits);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EFLS,
			      STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units of 1/P(f) are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICINVERSENOISETESTC_EFLS,
			    STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
  }

  if (!result)
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICINVERSENOISETESTC_MSGEFLS);
    }
    return STOCHASTICINVERSENOISETESTC_EFLS;
  }

   /* check output values */
   if (optVerbose)
   {
     printf("1/P(0)=%g, should be 0\n",
            invNoise.data->data[0]);
   }
   if ( fabs(invNoise.data->data[0]) > STOCHASTICINVERSENOISETESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
     }
     return STOCHASTICINVERSENOISETESTC_EFLS;
   }

  for (i=1; i < STOCHASTICINVERSENOISETESTC_LENGTH; i++)
  {
    f = i*STOCHASTICINVERSENOISETESTC_DELTAF;
    expectedReal = 2*f;
    if (optVerbose)
    {
      printf("1/P(%f Hz)=%g, should be %g\n",
	     f, invNoise.data->data[i], expectedReal);
    }
    if ( fabs(invNoise.data->data[i] - expectedReal)/expectedReal
	 > STOCHASTICINVERSENOISETESTC_TOL )
    {
      printf("  FAIL: Valid data test\n");
      if (optVerbose)
      {
	printf("Exiting with error: %s\n", STOCHASTICINVERSENOISETESTC_MSGEFLS);
      }
      return STOCHASTICINVERSENOISETESTC_EFLS;
    }
  }

  printf("  PASS: Valid data test #1\n");

  /******* clean up valid data ******/
  LALSDestroyVector(&status, &(wNoise .data));
  if ( ( code = CheckStatus (&status, 0 , "",
			     STOCHASTICINVERSENOISETESTC_EFLS,
			     STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSDestroyVector(&status, &(invNoise.data));
  if ( ( code = CheckStatus (&status, 0 , "",
			     STOCHASTICINVERSENOISETESTC_EFLS,
			     STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCDestroyVector(&status, &(wFilter.data));
  if ( ( code = CheckStatus (&status, 0 , "",
			     STOCHASTICINVERSENOISETESTC_EFLS,
			     STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCDestroyVector(&status, &(hwInvNoise.data));
  if ( ( code = CheckStatus (&status, 0 , "",
			     STOCHASTICINVERSENOISETESTC_EFLS,
			     STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
  {
    return code;
  }

   LALCheckMemoryLeaks();

   printf("PASS: all tests\n");


 /* VALID USER TEST DATA HERE ----------------------------------------- */

  if (optWNoiseFile[0] && optWFilterFile[0] &&
      optInvNoiseFile[0] && optHWInvNoiseFile[0])
  {
    /* allocate memory */
    wNoise.data      = NULL;
    wFilter.data     = NULL;
    invNoise.data    = NULL;
    hwInvNoise.data  = NULL;

    LALSCreateVector(&status, &(wNoise.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EUSE,
			      STOCHASTICINVERSENOISETESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EUSE,
			      STOCHASTICINVERSENOISETESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(wFilter.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EUSE,
			      STOCHASTICINVERSENOISETESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(hwInvNoise.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICINVERSENOISETESTC_EUSE,
			      STOCHASTICINVERSENOISETESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* Read input files */
    LALSReadFrequencySeries(&status, &wNoise, optWNoiseFile);
    LALCReadFrequencySeries(&status, &wFilter, optWFilterFile);

    /* fill inverse noise input and output */
    input.unCalibratedNoisePSD = &(wNoise);
    input.responseFunction            = &(wFilter);
    output.calibratedInverseNoisePSD    = &(invNoise);
    output.halfCalibratedInverseNoisePSD  = &(hwInvNoise);

    /*calculate the calibrated inverse noise and half-calibrated inverse noise*/
    LALStochasticInverseNoise(&status, &output, &input );
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EUSE,
			       STOCHASTICINVERSENOISETESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* print output files */
    LALSPrintFrequencySeries(output.calibratedInverseNoisePSD,
                             optInvNoiseFile);
    printf("====== Calibrated Inverse Noise PSD Written to File %s ======\n",
	   optInvNoiseFile);
    LALCPrintFrequencySeries(output.halfCalibratedInverseNoisePSD,
                             optHWInvNoiseFile);
    printf("===== Half-Calibrated Inverse Noise PSD Written to File %s =====\n",
	   optHWInvNoiseFile);


    /* clean up valid data */
    LALSDestroyVector(&status, &(wNoise .data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(wFilter.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(hwInvNoise.data));
    if ( ( code = CheckStatus (&status, 0 , "",
			       STOCHASTICINVERSENOISETESTC_EFLS,
			       STOCHASTICINVERSENOISETESTC_MSGEFLS) ) )
    {
      return code;
    }

  }
  LALCheckMemoryLeaks();
  return 0;

}


/*----------------------------------------------------------------------*/


/* Usage () Message */
/* Prints a usage message for program program and exits with code exitcode.*/

static void Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  fprintf (stderr, "  -n length      frequency series contain length points\n");
  fprintf (stderr, "  -w filename    read uncalibrated noise PSD from file filename\n");
  fprintf (stderr, "  -f filename    read response function from file filename \n");
  fprintf (stderr, "  -u filename    print calibrated inverse noise PSD to file filename\n");
  fprintf (stderr, "  -m filename    print half-calibrated inverse noise PSD to file filename\n");
  exit (exitcode);
}

/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  FILE *fp;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:n:w:f:u:m:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {

      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;

      case 'w': /* specify uncalibrated noise file */
        strncpy (optWNoiseFile, optarg, LALNameLength);
        break;

      case 'f': /* specify response function file */
        strncpy (optWFilterFile, optarg, LALNameLength);
        break;

      case 'u': /* specify calibrated inverse noise file */
        strncpy (optInvNoiseFile, optarg, LALNameLength);
        break;

      case 'm': /* specify hwInvNoise file */
        strncpy (optHWInvNoiseFile, optarg, LALNameLength);
        break;

      case 'd': /* set debug level */
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICINVERSENOISETESTC_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen ("/dev/null", "w", stderr);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen ("/dev/null", "w", stdout);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

/** \endcond */
