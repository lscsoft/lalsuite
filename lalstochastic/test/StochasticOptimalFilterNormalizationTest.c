/*
*  Copyright (C) 2007 Kaice T. Reilly, John Whelan
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
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 * \file
 * \ingroup StochasticOptimalFilterNormalization_c
 *
 * \brief Test suite for <tt>LALStochasticOptimalFilterNormalization()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./StochasticOptimalFilterNormalizationTest [options]
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -n length      frequency series contain length points
 * -f fRef        set normalization reference frequency to fRef
 * -w filename    read gravitational-wave spectrum from file filename
 * -g filename    read overlap reduction function from file filename
 * -i filename    read first inverse noise PSD from file filename
 * -j filename    read second inverse noise PSD from file filename
 * -y             use normalization appropriate to heterodyned data
 * \endcode
 *
 * ### Description ###
 *
 * This program tests the function
 * <tt>LALStochasticOptimalFilterNormalization()</tt>, which calculates
 * the normalization factor for an optimal filter and the expected
 * variance per unit time of the cross-correlation statistic, given a
 * stochastic gravitational-wave background spectrum
 * \f$h_{100}^2\Omega_{\mathrm{GW}}(f)\f$, an overlap reduction
 * function \f$\gamma(f)\f$, and unwhitened noise power spectral densities
 * \f$\{P_i(f)\f$ for a pair of detectors.
 *
 * First, it tests that the correct error codes
 * (cf. \ref StochasticCrossCorrelation_h)
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to input structure</em></li>
 * <li> <em>null pointer to output structure</em></li>
 * <li> <em>null pointer to overlap reduction function</em></li>
 * <li> <em>null pointer to gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to first inverse noise PSD</em></li>
 * <li> <em>null pointer to second inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of overlap reduction function</em></li>
 * <li> <em>null pointer to data member of gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to data member of first inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of second inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of data member of overlap reduction function</em></li>
 * <li> <em>null pointer to data member of data member of gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to data member of data member of first inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of data member of second inverse noise PSD</em></li>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * <li> negative start frequency</li>
 * <li> length mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> length mismatch between overlap reduction function and first inverse noise PSD</li>
 * <li> length mismatch between overlap reduction function and second inverse noise PSD</li>
 * <li> frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> frequency spacing mismatch between overlap reduction function and first inverse noise PSD</li>
 * <li> frequency spacing mismatch between overlap reduction function and second inverse noise PSD</li>
 * <li> start frequency mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> start frequency mismatch between overlap reduction function and first inverse noise PSD</li>
 * <li> start frequency mismatch between overlap reduction function and second inverse noise PSD</li>
 * <li> reference frequency less than frequency spacing</li>
 * <li> reference frequency greater than maximum frequency</li>
 * </ul>
 *
 * It then verifies that the correct optimal filter is generated
 * (checking the normalization by verifying that \eqref{stochastic_e_mu}
 * is satisfied) for each of the following simple test cases (with
 * \f$f_0=0\f$, \f$f_{\mathrm{R}}=1\f$, \f$\delta f=1\f$ and \f$N=8\f$):
 * <ol>
 * <li> \f$\gamma(f) = h_{100}^2\Omega_{\mathrm{GW}}(f)
 * = P^{\mathrm{C}}_1(f)
 * = P^{\mathrm{C}}_2(f) = 1\f$;
 * The expected results in this case are \f$\lambda=3.079042427975\f$,
 * \f$\sigma^2\,T^{-1}=4.822422518205\f$</li>
 * <li> \f$\gamma(f) = P^{\mathrm{C}}_1(f)
 * = P^{\mathrm{C}}_2(f) = 1\f$;
 * \f$h_{100}^2\Omega_{\mathrm{GW}}(f)=f^3\f$.
 * The expected results in this case are \f$\lambda=.4474881517327\f$,
 * \f$\sigma^2\,T^{-1}=.700859760794\f$</li>
 * </ol>
 *
 * ### Exit codes ###
 *
 *
 * ### Uses ###
 *
 * \code
 * LALStochasticOptimalFilterNormalization()
 * LALCheckMemoryLeaks()
 * LALCReadFrequencySeries()
 * LALSCreateVector()
 * LALSDestroyVector()
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
 * data error.  If \c fRef is unspecified, a default value of
 * 1\,Hz is used.</li>
 * <li> The length of the user-provided series must be specified, even
 * though it could in principle be deduced from the input file, because
 * the data sequences must be allocated before the
 * <tt>LALCReadFrequencySeries()</tt> function is called.</li>
 * <li> If some, but not all, of the \c filename arguments are
 * present, the user-specified data will be silently ignored.</li>
 * </ul>
 *
 */

/**\name Error Codes */ /*@{*/
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ENOM 0	/**< Nominal exit */
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EARG 1	/**< Error parsing command-line arguments */
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK 2	/**< Error checking failed to catch bad data */
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS 3	/**< Incorrect answer for valid data */
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EUSE 4	/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGENOM "Nominal exit"
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEUSE "Bad user-entered data"



#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

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

#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TRUE     1
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FALSE    0
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF   1.0
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0       0.0
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF     1.0
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH   8
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TOL      1e-6

#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA1     3.079042427975
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT1    4.822422518205
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA2     .4474881517327
#define STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT2    .700859760794

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FALSE;
BOOLEAN optHetero = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FALSE;
REAL8 optFRef      = 1.0;
UINT4 optLength    = 0.0;
CHAR  optOverlapFile[LALNameLength]   = "";
CHAR  optOmegaFile[LALNameLength]     = "";
CHAR  optInvNoise1File[LALNameLength] = "";
CHAR  optInvNoise2File[LALNameLength] = "";
CHAR  optOptimalFile[LALNameLength]     = "";

INT4 code;

static void Usage (const char *program, int exitflag);
static void ParseOptions (int argc, char *argv[]);

int main(int argc, char *argv[])
{

  static LALStatus         status;

  StochasticOptimalFilterNormalizationInput             input;
  StochasticOptimalFilterNormalizationOutput            output;
  StochasticOptimalFilterNormalizationParameters        params;

  REAL4WithUnits           normOut, varOut;

  LIGOTimeGPS              epoch = {1,0};

  REAL4FrequencySeries     overlap;
  REAL4FrequencySeries     omegaGW;
  REAL4FrequencySeries     invNoise1;
  REAL4FrequencySeries     invNoise2;

  INT4       i;
  REAL8      f;

  LALUnitPair              unitPair;
  LALUnit                  expectedUnit;
  BOOLEAN                  result;

  CHARVector               *unitString = NULL;


  ParseOptions (argc, argv);

  /* define valid params */

  overlap.name[0] = '\0';
  overlap.f0     = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0;
  overlap.deltaF = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;
  overlap.epoch  = epoch;
  overlap.data   = NULL;
  overlap.sampleUnits = lalDimensionlessUnit;
  omegaGW = invNoise1 = invNoise2 = overlap;
#ifndef LAL_NDEBUG
  REAL4FrequencySeries     realBadData = omegaGW;
#endif

  invNoise1.sampleUnits.unitNumerator[LALUnitIndexStrain] = -2;
  invNoise1.sampleUnits.unitNumerator[LALUnitIndexSecond] = -1;
  invNoise1.sampleUnits.powerOfTen = 36;
  invNoise2.sampleUnits = invNoise1.sampleUnits;

  overlap.sampleUnits.unitNumerator[LALUnitIndexStrain] = 2;

  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF;
  params.heterodyned = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FALSE;
  params.window1 = params.window2 = NULL;

  strncpy(overlap.name, "", LALNameLength);

  output.normalization = &normOut;
  output.variance = &varOut;

  /* allocate memory */
  LALSCreateVector(&status, &(overlap.data),
                   STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALSCreateVector(&status, &(omegaGW.data),
                   STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise1.data),
                   STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise2.data),
                   STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.overlapReductionFunction = &overlap;
  input.omegaGW = &omegaGW;
  input.inverseNoisePSD1 = &invNoise1;
  input.inverseNoisePSD2 = &invNoise2;

  /* TEST INVALID DATA HERE -------------------------------------- */
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to input structure */
    LALStochasticOptimalFilterNormalization(&status, &output, NULL, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to output structure */
    LALStochasticOptimalFilterNormalization(&status, NULL, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to output structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to overlap reduction function */
    input.overlapReductionFunction = NULL;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to overlap reduction function results in error:\n      \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to gravitational-wave spectrum */
    input.omegaGW = NULL;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to gravitational-wave spectrum results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* test behavior for null pointer to inverse noise 1 member */
    /* of input structure */
    input.inverseNoisePSD1 = NULL;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD1 = &invNoise1;

    /* test behavior for null pointer to inverse noise 2 member */
    /* of input structure*/
    input.inverseNoisePSD2 = NULL;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD2 = &invNoise2;

    /* test behavior for null pointer to data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* test behavior for null pointer to data member of first inverse noise PSD */
    input.inverseNoisePSD1 = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD1 = &invNoise1;

    /* test behavior for null pointer to data member of second inverse noise PSD */
    input.inverseNoisePSD2 = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD2 = &invNoise2;

    /* Create a vector for testing REAL4 null data-data pointers */
    LALSCreateVector(&status, &(realBadData.data), STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    REAL4*                   realTempPtr;
    realTempPtr = realBadData.data->data;
    realBadData.data->data = NULL;

    /* test behavior for null pointer to data member of data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to data member of data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* test behavior for null pointer to data member of data member of */
    /* first inverse noise PSD */
    input.inverseNoisePSD1 = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD1 = &invNoise1;

    /* test behavior for null pointer to data member of data member of */
    /* second inverse noise PSD */
    input.inverseNoisePSD2 = &realBadData;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.inverseNoisePSD2 = &invNoise2;

    /** clean up **/
    realBadData.data->data = realTempPtr;
    LALSDestroyVector(&status, &(realBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* test behavior for zero length */
    overlap.data->length =
      omegaGW.data->length = invNoise1.data->length =
      invNoise2.data->length = 0;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
                              STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    overlap.data->length =
      omegaGW.data->length = invNoise1.data->length =
      invNoise2.data->length = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH;

    /* test behavior for negative frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF = -3.5;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
                              STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF
      = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

    /* test behavior for zero frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF = 0;
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
                              STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF =
      omegaGW.deltaF = invNoise1.deltaF =
      invNoise2.deltaF = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;
  } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  overlap.f0 = omegaGW.f0
    = invNoise1.f0
    = invNoise2.f0 = -3.0;

  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);

  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n      \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  overlap.f0 = omegaGW.f0 =
    invNoise1.f0 = invNoise2.f0 = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0;

  /* test behavior for length mismatch between overlap reduction function and omega */
  omegaGW.data->length = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH - 1;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  omegaGW.data->length = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH;

  /* test behavior length mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.data->length
    = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH - 1;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
   {
     return code;
   }
  printf("  PASS: length mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  invNoise1.data->length
    = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH;

  /* test behavior length mismatch between overlap reduction function and second inverse noise PSD */
  invNoise2.data->length = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH - 1;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  invNoise2.data->length = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH;

  /* test behavior for frequency spacing mismatch between overlap reduction function and omega */
  omegaGW.deltaF = 2.0 * STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  omegaGW.deltaF = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

  /* test behavior frequency spacing mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.deltaF = 2.0 * STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  invNoise1.deltaF = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

  /* test behavior frequency spacing mismatch between overlap reduction function and inverse */
  /* noise 2 */
  invNoise2.deltaF = 2.0 * STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  invNoise2.deltaF = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

  /* test behavior for start frequency mismatch between overlap reduction function and omega */
  omegaGW.f0 = 30;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n         \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  omegaGW.f0 = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0;

  /* test behavior start frequency mismatch between overlap reduction function and first inverse noise PSD */
  invNoise1.f0 = 30;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and first inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  invNoise1.f0 = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0;

  /* test behavior start frequency mismatch between overlap reduction function and second inverse noise PSD */
  invNoise2.f0 = 30;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and second inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  invNoise2.f0 = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_F0;

  /* test behavior for reference frequency to be less than frequency spacing */
  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF/2;
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EOORFREF,
                            STOCHASTICCROSSCORRELATIONH_MSGEOORFREF,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: reference frequency less than frequency spacing results in error:\n       \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF;

  /* test behavior for reference frequency to be greater than its maximum */
  params.fRef = (STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH*STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF);
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EOORFREF,
                            STOCHASTICCROSSCORRELATIONH_MSGEOORFREF,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: reference frequency greater than maximum frewquency results in error:\n       \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF;


 /* VALID TEST DATA HERE ----------------------------------------- */

  /* create input to test */
  overlap.data->data[0] = 1;
  omegaGW.data->data[0] = 0;
  invNoise1.data->data[0] = 0;
  invNoise2.data->data[0] = 0;

  /** Test 1 **/
  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF;
  for (i=1; i < STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

    overlap.data->data[i] = 1;
    omegaGW.data->data[i] = 1;
    invNoise1.data->data[i] = 1;
    invNoise2.data->data[i] = 1;
  }

  /* fill optimal input */
  input.overlapReductionFunction    = &(overlap);
  input.omegaGW                     = &(omegaGW);
  input.inverseNoisePSD1             = &(invNoise1);
  input.inverseNoisePSD2             = &(invNoise2);

  /* calculate optimal filter normalization */
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output normalization */
  if (optVerbose)
  {
    printf("lambda=%g, should be %g\n", normOut.value,
           STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA1);
  }
  if ( fabs(normOut.value
            - STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA1)
       / STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA1
       > STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check units of normalization*/

  expectedUnit = lalDimensionlessUnit;
  expectedUnit.unitNumerator[LALUnitIndexSecond] = -1;
  expectedUnit.powerOfTen = -36;
  unitPair.unitOne = &expectedUnit;
  unitPair.unitTwo = &(normOut.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check output variance per unit time */
  if (optVerbose)
  {
    printf("sigma^2/T=%g, should be %g\n", varOut.value,
           STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT1);
  }
  if ( fabs(varOut.value
            - STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT1)
       / STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT1
       > STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check units of variance per unit time */

  unitPair.unitOne = &lalSecondUnit;
  unitPair.unitTwo = &(varOut.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  printf("  PASS: Valid data test #1\n");

  /** Test 2 **/
  params.fRef = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_FREF;
  for (i=0; i < STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LENGTH; i++)
     {
       f = i*STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_DELTAF;

       overlap.data->data[i] = 1;
       omegaGW.data->data[i] = pow(f,3);
       invNoise1.data->data[i] = 1;
       invNoise2.data->data[i] = 1;
    }

  /* fill optimal input */
  input.overlapReductionFunction    = &(overlap);
  input.omegaGW                     = &(omegaGW);
  input.inverseNoisePSD1             = &(invNoise1);
  input.inverseNoisePSD2             = &(invNoise2);

  /* calculate optimal filter */
  LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output normalization */
  if (optVerbose)
  {
    printf("lambda=%g, should be %g\n", normOut.value,
           STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA2);
  }
  if ( fabs(normOut.value
            - STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA2)
       / STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_LAMBDA2
       > STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check units of normalization*/

  expectedUnit = lalDimensionlessUnit;
  expectedUnit.unitNumerator[LALUnitIndexSecond] = -1;
  expectedUnit.powerOfTen = -36;
  unitPair.unitOne = &expectedUnit;
  unitPair.unitTwo = &(normOut.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check output variance per unit time */
  if (optVerbose)
  {
    printf("sigma^2/T=%g, should be %g\n", varOut.value,
           STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT2);
  }
  if ( fabs(varOut.value
            - STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT2)
       / STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_VARPERT2
       > STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  /* check units of variance per unit time */

  unitPair.unitOne = &lalSecondUnit;
  unitPair.unitTwo = &(varOut.units);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS;
  }

  printf("  PASS: Valid data test #2\n");

  /* clean up valid data */
  LALSDestroyVector(&status, &(overlap.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(omegaGW.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise1.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise2.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
  {
       return code;
  }

  LALCheckMemoryLeaks();

  printf("PASS: all tests\n");


 /* VALID USER TEST DATA HERE ----------------------------------------- */

  if (optOverlapFile[0] &&  optOmegaFile[0] && optInvNoise1File[0] && optInvNoise2File[0])
  {

    /* allocate memory */
    overlap.data     = NULL;
    omegaGW.data     = NULL;
    invNoise1.data   = NULL;
    invNoise2.data   = NULL;

    LALSCreateVector(&status, &(overlap.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALSCreateVector(&status, &(omegaGW.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise1.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise2.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* Read input files */
    LALSReadFrequencySeries(&status, &overlap,   optOverlapFile);
    LALSReadFrequencySeries(&status, &omegaGW,   optOmegaFile);
    LALSReadFrequencySeries(&status, &invNoise1, optInvNoise1File);
    LALSReadFrequencySeries(&status, &invNoise2, optInvNoise2File);

    /* fill optimal input */
    input.overlapReductionFunction    = &(overlap);
    input.omegaGW                     = &(omegaGW);
    input.inverseNoisePSD1             = &(invNoise1);
    input.inverseNoisePSD2             = &(invNoise2);

    params.fRef = optFRef;
    params.heterodyned = optHetero;

    /* calculate optimal filter normalization */
    LALStochasticOptimalFilterNormalization(&status, &output, &input, &params);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EUSE,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* Convert Unit Structure to String */
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, &(normOut.units) );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    printf("=========== Optimal Filter Normalization for User-Specified Data Is =======\n");
    printf("     %g %s\n", normOut.value, unitString->data);

    /* Convert Unit Structure to String */
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, &(varOut.units) );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }

    printf("=========== Variance per Unit Time for User-Specified Data Is =======\n");
    printf("     %g %s\n", varOut.value, unitString->data);

    /* clean up */
    LALSDestroyVector(&status, &(overlap.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(omegaGW.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise1.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise2.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_MSGEFLS) ) )
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
  fprintf (stderr, "  -f fRef        set normalization reference frequency to fRef\n");
  fprintf (stderr, "  -w filename    read gravitational-wave spectrum from file filename\n");
  fprintf (stderr, "  -g filename    read overlap reduction function from file filename\n");
  fprintf (stderr, "  -i filename    read first inverse noise PSD from file filename\n");
  fprintf (stderr, "  -j filename    read second inverse noise PSD from file filename\n");
  fprintf (stderr, "  -o filename    print optimal filter to file filename\n");
  fprintf (stderr, "  -y             use normalization appropriate to heterodyned data\n");
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

    c = getopt (argc, argv, "hqvd:n:f:w:g:i:j:s:t:o:y");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'o': /* specify output file */
        strncpy (optOptimalFile, optarg, LALNameLength);
        break;

      case 'f': /* specify refernce frequency */
        optFRef = atoi (optarg);
        break;

      case 'n': /* specify number of points in frequency series */
        optLength = atoi (optarg);
        break;

      case 'w': /* specify omegaGW file */
        strncpy (optOmegaFile, optarg, LALNameLength);
        break;

      case 'g': /* specify overlap file */
        strncpy (optOverlapFile, optarg, LALNameLength);
        break;

      case 'i': /* specify InvNoise1 file */
        strncpy (optInvNoise1File, optarg, LALNameLength);
        break;

      case 'j': /* specify InvNoise2 file */
        strncpy (optInvNoise2File, optarg, LALNameLength);
        break;

      case 'd': /* set debug level */
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TRUE;
        break;

      case 'y': /* optHetero */
        optHetero = STOCHASTICOPTIMALFILTERNORMALIZATIONTESTC_TRUE;
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
