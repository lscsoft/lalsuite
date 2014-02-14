/*
*  Copyright (C) 2007 John Whelan
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
 * \author UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
 * \file
 * \ingroup StochasticCrossCorrelation_c
 *
 * \brief A program to test <tt>LALStochasticHeterodynedCrossCorrelationStatistic()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./StochasticHeterodynedCrossCorrelationStatisticTest [options]
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -i filename    read first data stream from file filename
 * -j filename    read second data stream from file filename
 * -k filename    read optimal filter from file filename
 * -n length      frequency series contain length points
 * -t             epochs need not match
 * \endcode
 *
 * This program tests the function
 * <tt>LALStochasticHeterodynedCrossCorrelationStatistic()</tt>, which calculates
 * the cross-correlation statistic given two zero-padded and
 * Fourier-transformed data streams and a (frequency domain) optimal
 * filter.
 *
 * First, it tests that the correct error codes
 * (cf. \ref StochasticCrossCorrelation_h)
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to output structure</em></li>
 * <li> <em>null pointer to input structure</em></li>
 * <li> <em>null pointer to first data stream</em></li>
 * <li> <em>null pointer to second data stream</em></li>
 * <li> <em>null pointer to optimal filter</em></li>
 * <li> <em>null pointer to data member of first data stream</em></li>
 * <li> <em>null pointer to data member of second data stream</em></li>
 * <li> <em>null pointer to data member of optimal filter</em></li>
 * <li> <em>null pointer to data member of data member of first data stream</em></li>
 * <li> <em>null pointer to data member of data member of second data stream</em></li>
 * <li> <em>null pointer to data member of data member of optimal filter</em></li>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * <li> negative start frequency</li>
 * <li> length mismatch between optimal filter and first data stream</li>
 * <li> length mismatch between optimal filter and second data stream</li>
 * <li> frequency spacing mismatch between optimal filter and first data stream</li>
 * <li> frequency spacing mismatch between optimal filter and second data stream</li>
 * <li> start frequency mismatch between optimal filter and first data stream</li>
 * <li> start frequency mismatch between optimal filter and second data stream</li>
 * <li> mismatch between epochs of data streams</li>
 * </ul>
 *
 * It then verifies that the correct cross-correlation statistic (value
 * and units) is generated for each of the following simple test cases:
 * <ol>
 * <li> \f$\widetilde{Q}(f) = \frac{f(N\,\delta f - f)}{(N\,\delta
 * f/2)^2}\f$; \f$\widetilde{\bar{h}}_1(f)=f^2+if\f$,
 * \f$\widetilde{\bar{h}}_2(f)=f^{-2}-if^{-1}\f$.  With \f$f_0=\delta
 * f=80\,\textrm{Hz}\f$ and \f$N=9\f$, the expected value is
 * \f$-1248i\f$.</li>
 * <li> \f$\widetilde{Q}(f) = 1\f$ for
 * \f$300\,\textrm{Hz}<f<500\,\textrm{Hz}\f$, 0 otherwise;
 * \f$\widetilde{\bar{h}}_1(f)=1-\widetilde{\bar{h}}_2(f)=f/800\,\textrm{Hz}\f$.
 * With \f$f_0=\delta f=80\,\textrm{Hz}\f$ and \f$N=9\f$, the expected value is
 * \f$58.4\f$.</li>
 * </ol>
 * For each successful test
 * (both of these valid data and the invalid ones described above), it
 * prints "\c PASS" to standard output; if a test fails, it
 * prints "\c FAIL".
 *
 * If the \c filename arguments are present, it also reads in the
 * optimal filter and the two data streams from the specified files and
 * use the specified parameters to calculate the cross-correlation
 * statistic.  The result is printed to standard output along with the
 * resulting units in terms of the basic SI units.
 *
 * ### Uses ###
 *
 * \code
 * LALStochasticHeterodynedCrossCorrelationStatistic()
 * LALCheckMemoryLeaks()
 * LALCReadFrequencySeries()
 * LALCCreateVector()
 * LALCDestroyVector()
 * LALCHARCreateVector()
 * LALCHARDestroyVector()
 * LALUnitAsString()
 * LALUnitCompare()
 * getopt()
 * printf()
 * fprintf()
 * freopen()
 * fabs()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 * <li> No specific error checking is done on user-specified data.  If
 * \c length is missing, the resulting default will cause a bad
 * data error.
 * </li><li> The length of the user-provided series must be specified, even
 * though it could in principle be deduced from the input file,
 * because the data sequences must be allocated before the
 * <tt>LALCReadFrequencySeries()</tt> function is called.
 * </li><li> If some, but not all, of the \c filename arguments are
 * present, the user-specified data will be silently ignored.</li>
 * </ul>
 *
 */

/**\name Error Codes */ /*@{*/
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ENOM 0	/**< Nominal exit */
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EARG 1	/**< Error parsing command-line arguments */
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK 2	/**< Error checking failed to catch bad data */
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS 3	/**< Incorrect answer for valid data */
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE 4	/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGENOM "Nominal exit"
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE "Bad user-entered data"



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
#include <lal/Units.h>

#include "CheckStatus.h"
#include "CheckStatus.c"

#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH    9
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0        80.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF    80.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TOL       1e-6

#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_WINMIN   300.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_WINMAX   500.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FLIM     800.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP1   -1248.0
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP2      58.4

#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE     1
#define STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FALSE    0

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FALSE;
BOOLEAN optMatch   = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE;
UINT4 optLength     = 0;
CHAR optData1File[LALNameLength] = "";
CHAR optData2File[LALNameLength] = "";
CHAR optFilterFile[LALNameLength] = "";

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

int main( int argc, char *argv[] )
{
  static LALStatus                status;

  StochasticCrossCorrelationInput           input;
  COMPLEX8WithUnits           output;

  COMPLEX8FrequencySeries  goodData1;
  COMPLEX8FrequencySeries  goodData2;
  COMPLEX8FrequencySeries  goodFilter;

  LIGOTimeGPS              epoch0 = {0,0};
  LIGOTimeGPS              epoch1 = {630720000,123456789};
  LIGOTimeGPS              epoch2 = {630720000,987654321};
  LIGOTimeGPS              epoch3 = {630722222,123456789};

  LALUnitPair              unitPair;
  BOOLEAN                  result;

  CHARVector               *unitString = NULL;

  UINT4 i;
  REAL4 f, x;
  INT4 code;


  ParseOptions( argc, argv );

  /* define valid parameters */

  goodFilter.f0     = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;
  goodFilter.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;
  goodFilter.epoch  = epoch0;
  goodFilter.data   = NULL;

#ifndef LAL_NDEBUG
  COMPLEX8FrequencySeries  badFilter = goodFilter;
#endif

  goodData1 = goodFilter;

  goodData1.epoch = epoch1;
  goodData2 = goodData1;
#ifndef LAL_NDEBUG
  COMPLEX8FrequencySeries  badData1 = goodData2;
  COMPLEX8FrequencySeries  badData2 = badData1;
#endif

  LALCCreateVector(&status, &(goodData1.data),
                          STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(goodData2.data),
                          STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(goodFilter.data),
                          STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.hBarTildeOne  = &goodData1;
  input.hBarTildeTwo  = &goodData2;
  input.optimalFilter = &goodFilter;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to output structure */
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, NULL, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
      {
        return code;
      }
    printf("  PASS: null pointer to output structure results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to input structure */
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, NULL, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to first data stream */
    input.hBarTildeOne = NULL;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to second data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to second data stream */
    input.hBarTildeTwo = NULL;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to optimal filter */
    input.optimalFilter = NULL;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to optimal filter */
    input.optimalFilter = &goodFilter;

    /* test behavior for null pointer to data member of first data stream */
    input.hBarTildeOne = &badData1;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of first data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to data member of second data stream */
    input.hBarTildeTwo = &badData2;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to data member of optimal filter */
    input.optimalFilter = &badFilter;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of optimal filter */
    input.optimalFilter = &goodFilter;

    /* Create a vector for testing null data-data pointers */
    LALCCreateVector(&status, &(badFilter.data),
                          STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    COMPLEX8               *tempPtr;
    tempPtr = badFilter.data->data;
    badFilter.data->data = NULL;
    badData1.data = badData2.data = badFilter.data;

    /* test behavior for null pointer to data member of data member of first data stream */
    input.hBarTildeOne = &badData1;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of first data stream */
    input.hBarTildeOne = &goodData1;

    /* test behavior for null pointer to data member of data member of second data stream */
    input.hBarTildeTwo = &badData2;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second data stream results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of second data stream */
    input.hBarTildeTwo = &goodData2;

    /* test behavior for null pointer to data member of data member of optimal filter */
    input.optimalFilter = &badFilter;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
			      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of optimal filter results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* assign valid pointer to data member of data member of optimal filter */
    input.optimalFilter = &goodFilter;

    /* clean up */

    badFilter.data->data = tempPtr;
    LALCDestroyVector(&status, &(badFilter.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    badData1.data = badData2.data = badFilter.data;

    /* test behavior for zero length */
    goodData1.data->length = goodData2.data->length
      = goodFilter.data->length = 0;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
			      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    goodData1.data->length = goodData2.data->length
      = goodFilter.data->length = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH;

    /* test behavior for negative frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = -STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

    /* test behavior for zero frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = 0;
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
			      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    goodData1.deltaF = goodData2.deltaF
      = goodFilter.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;
   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  goodData1.f0 = goodData2.f0
    = goodFilter.f0 = -20.0;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);

  /* reassign valid start frequency */
  goodData1.f0 = goodData2.f0
    = goodFilter.f0 = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;

  /* test behavior for length mismatch
     between optimal filter and first data stream */
  goodData1.data->length = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH - 1;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between optimal filter and first data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);

  /* reassign correct length */
  goodData1.data->length = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH;

  /* test behavior for length mismatch
     between optimal filter and first data stream */
  goodData2.data->length = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH - 1;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between optimal filter and second data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);

  /* reassign correct length */
  goodData2.data->length = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH;

  /* test behavior for frequency spacing mismatch
     between optimal filter and first data stream */
  goodData1.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF * 2.0;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: frequency spacing mismatch between optimal filter and first data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);

  /* reassign correct frequency spacing */
  goodData1.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;

  /* test behavior for frequency spacing mismatch
     between optimal filter and second data stream */
  goodData2.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF * 2.0;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: frequency spacing mismatch between optimal filter and second data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);

  /* reassign correct frequency spacing */
  goodData2.deltaF = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;

  /* test behavior for start frequency mismatch
     between optimal filter and first data stream */
  goodData1.f0 = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0 + 2.0;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between optimal filter and first data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);

  /* reassign correct start frequency */
  goodData1.f0 = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;

  /* test behavior for start frequency mismatch
     between optimal filter and second data stream */
  goodData2.f0 = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0 + 2.0;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between optimal filter and second data stream results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);

  /* reassign correct start frequency */
  goodData2.f0 = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;

  /* test behavior for mismatch between epochs of data streams */
  goodData2.epoch = epoch2;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMTIME,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMTIME,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  goodData2.epoch = epoch3;
  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMTIME,
			    STOCHASTICCROSSCORRELATIONH_MSGEMMTIME,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ECHK,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: mismatch between epochs of data streams results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMTIME);

  /* reassign correct epoch */
  goodData2.epoch = epoch1;

  /******************** Test Valid Data Case #1 ***********************/
  goodData1.sampleUnits = lalDimensionlessUnit;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexStrain] = 1;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
  goodData2.sampleUnits = goodData1.sampleUnits;
  goodFilter.sampleUnits = lalDimensionlessUnit;
  goodFilter.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;

  goodData1.f0 = goodData2.f0 = goodFilter.f0
    = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;

  for (i=0; i<STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH; ++i)
  {
    f = goodData1.f0
      + i * STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;
    x = f
      / (STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FLIM / 2.0);
    /* printf ("%f\n",x); */
    goodData1.data->data[i] = crectf( x*x, x );
    goodData2.data->data[i] = crectf( 1.0/crealf(goodData1.data->data[i]), -1.0/cimagf(goodData1.data->data[i]) );
    goodFilter.data->data[i] = crectf( x * (2-x), 0.0 );
    /*    printf ("%f + %f i    %f + %f i    %f + %f i\n",
	    goodData1.data->data[i].re, goodData1.data->data[i].im,
	    goodData2.data->data[i].re, goodData2.data->data[i].im,
	    goodFilter.data->data[i].re, goodFilter.data->data[i].im
	    ); */
  }

  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    printf("Y=%g + %g i, should be %g i\n", crealf(output.value),
	   cimagf(output.value),
	   STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP1);
  }
  if ( ( fabs(crealf(output.value)/STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP1)
	 > STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TOL )
       || ( fabs((cimagf(output.value) - STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP1)
		 / STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP1)
	    > STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TOL )
       )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS);
    }
    return STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS;
  }

  unitPair.unitOne = &(goodData1.sampleUnits);
  unitPair.unitTwo = &(output.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
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
             STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS);
    }
    return STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS;
  }

  printf("  PASS: Valid data test #1\n");


  /********************** Test Valid Data Case #2 *****************/
  goodData1.sampleUnits = lalDimensionlessUnit;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexStrain] = 1;
  goodData1.sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
  goodData2.sampleUnits = lalDimensionlessUnit;
  goodData2.sampleUnits.unitNumerator[LALUnitIndexADCCount] = 1;
  goodData2.sampleUnits.unitNumerator[LALUnitIndexSecond] = 1;
  goodFilter.sampleUnits = lalDimensionlessUnit;
  goodFilter.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;
  goodFilter.sampleUnits.unitNumerator[LALUnitIndexADCCount] = -1;

  goodData1.f0 = goodData2.f0 = goodFilter.f0
    = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0;

  for (i=0; i<STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_LENGTH; ++i)
  {
    f = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_F0
      + i * STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_DELTAF;
    goodData1.data->data[i] = crectf( f/STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FLIM, 0.0 );
    goodData2.data->data[i] = crectf( 1 - crealf(goodData1.data->data[i]), 0.0 );
    if ( f > STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_WINMIN
         && f < STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_WINMAX )
    {
      goodFilter.data->data[i] = 1.0;
    }
    else
    {
      goodFilter.data->data[i] = 0.0;
    }
  }

  LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    printf("Y = %g + %g i, should be %g\n", crealf(output.value), cimagf(output.value),
           STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP2);
  }
  if ( ( fabs(crealf(output.value)-STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP2)
         / STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP2
         > STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TOL )
       || ( fabs(cimagf(output.value))
            / STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EXP2
            > STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TOL )
       )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS);
    }
    return STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS;
  }


  unitPair.unitOne = &lalSecondUnit;
  unitPair.unitTwo = &(output.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }

  if (!result)
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS);
    }
    return STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS;
  }

  printf("  PASS: Valid data test #2\n");


  /* clean up */
  LALCDestroyVector(&status, &(goodFilter.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCDestroyVector(&status, &(goodData1.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCDestroyVector(&status, &(goodData2.data));
  if ( ( code = CheckStatus(&status, 0 , "",
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			    STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
  {
    return code;
  }

  printf("PASS: all tests\n");
  LALCheckMemoryLeaks();

  if (optData1File[0] && optData2File[0] && optFilterFile[0])
  {

    /* Allocate Memory */
    LALCCreateVector(&status, &(goodFilter.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(goodData1.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(goodData2.data), optLength);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    /* Read Data From Files */
    LALCReadFrequencySeries(&status, &(goodFilter), optFilterFile);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCReadFrequencySeries(&status, &(goodData1), optData1File);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    LALCReadFrequencySeries(&status, &(goodData2), optData2File);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }
    /* Calculate CC Statistic */
    LALStochasticHeterodynedCrossCorrelationStatistic(&status, &output, &input, STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EUSE,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* Convert Unit Structure to String */
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, &(output.units) );
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }

    printf("=========== Cross-Correlation Statistic for User-Specified Data Is =======\n");
    printf("     %g + %g i %s\n", crealf(output.value), cimagf(output.value),
           unitString->data);

    /* Deallocate Memory */
    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodFilter.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodData1.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(goodData2.data));
    if ( ( code = CheckStatus(&status, 0 , "",
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_EFLS,
			      STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }


  /* normal exit */
  LALCheckMemoryLeaks();
  return STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_ENOM;
}

/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  fprintf (stderr, "  -i filename    read first data stream from file filename\n");
  fprintf (stderr, "  -j filename    read second data stream from file filename\n");
  fprintf (stderr, "  -k filename    read optimal filter from file filename\n");
  fprintf (stderr, "  -n length      frequency series contain length points\n");
  fprintf (stderr, "  -t             epochs need not match\n");
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

    c = getopt (argc, argv, "hqvd:i:j:k:n:t");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 't': /* epochs need not match */
        optMatch = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_FALSE;
        break;

      case 'i': /* specify file containing first data stream */
        strncpy (optData1File, optarg, LALNameLength);
        break;

      case 'j': /* specify file containing second data stream */
        strncpy (optData2File, optarg, LALNameLength);
        break;

      case 'k': /* specify file containing optimal filter */
        strncpy (optFilterFile, optarg, LALNameLength);
        break;

      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;

      case 'd': /* set debug level */
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICHETERODYNEDCROSSCORRELATIONSTATISTICTESTC_TRUE;
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
