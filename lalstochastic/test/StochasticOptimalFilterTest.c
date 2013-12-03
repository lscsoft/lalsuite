/*
*  Copyright (C) 2007 Jolien Creighton, Kaice T. Reilly, John Whelan
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
 * \ingroup StochasticOptimalFilter_c
 *
 * \brief Test suite for <tt>LALStochasticOptimalFilter()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./StochasticOptimalFilterTest [options]
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -n length      frequency series contain length points
 * -f fRef        set normalization reference frequency to fRef
 * -w filename    read gravitational-wave spectrum from file filename
 * -g filename    read overlap reduction function from file filename
 * -i filename    read first calibrated inverse noise PSD from file filename
 * -j filename    read second calibrated inverse noise PSD from file filename
 * -s filename    read first half-calibrated inverse noise PSD from file filename
 * -t filename    read second half-calibrated inverse noise PSD from file filename
 * -o filename    print optimal filter to file filename
 * -y             use normalization appropriate to heterodyned data
 * \endcode
 *
 * ### Description ###
 *
 * This program tests the function <tt>LALStochasticOptimalFilter()</tt>,
 * which generates a normalized optimal filter from a stochastic
 * gravitational-wave background spectrum
 * \f$h_{100}^2\Omega_{\mathrm{GW}}(f)\f$, an overlap reduction
 * function \f$\gamma(f)\f$, and calibrated and half-calibrated noise power
 * spectral densities \f$\{P^{\mathrm{C}}_i(f),P^{\mathrm{HC}}_i(f)\}\f$ for a
 * pair of detectors.
 *
 * First, it tests that the correct error codes
 * (cf. \ref StochasticCrossCorrelation_h)
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to input structure</em></li>
 * <li> <em>null pointer to output series</em></li>
 * <li> <em>null pointer to overlap reduction function</em></li>
 * <li> <em>null pointer to gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to first half-calibrated inverse noise PSD</em></li>
 * <li> <em>null pointer to second half-calibrated inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of output series</em></li>
 * <li> <em>null pointer to data member of overlap reduction function</em></li>
 * <li> <em>null pointer to data member of gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to data member of first half-calibrated inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of second half-calibrated inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of data member of output series</em></li>
 * <li> <em>null pointer to data member of data member of overlap reduction function</em></li>
 * <li> <em>null pointer to data member of data member of gravitational-wave spectrum</em></li>
 * <li> <em>null pointer to data member of data member of first half-calibrated inverse noise PSD</em></li>
 * <li> <em>null pointer to data member of data member of second half-calibrated inverse noise PSD</em></li>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * <li> negative start frequency</li>
 * <li> length mismatch between overlap reduction function and output series</li>
 * <li> length mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> length mismatch between overlap reduction function and first half-calibrated inverse noise PSD</li>
 * <li> length mismatch between overlap reduction function and second half-calibrated inverse noise PSD</li>
 * <li> frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> frequency spacing mismatch between overlap reduction function and first half-calibrated inverse noise PSD</li>
 * <li> frequency spacing mismatch between overlap reduction function and second half-calibrated inverse noise PSD</li>
 * <li> start frequency mismatch between overlap reduction function and gravitational-wave spectrum</li>
 * <li> start frequency mismatch between overlap reduction function and first half-calibrated inverse noise PSD</li>
 * <li> start frequency mismatch between overlap reduction function and second half-calibrated inverse noise PSD</li>
 * <li> reference frequency less than frequency spacing</li>
 * <li> reference frequency greater than maximum frequency</li>
 * </ul>
 *
 * It then verifies that the correct optimal filter is generated
 * [calculating the normalization with
 * <tt>LALStochasticOptimalFilterNormalization()</tt> as described in
 * \ref StochasticOptimalFilterNormalization_c, and
 * checking the normalization by verifying that \eqref{stochastic_e_mu}
 * is satisfied] for each of the following simple test cases:
 * <ol>
 * <li> \f$\gamma(f) = h_{100}^2\Omega_{\mathrm{GW}}(f) = P^{\mathrm{C}}_1(f)
 * =P^{\mathrm{C}}_2(f)=P^{\mathrm{HC}}_1(f)=P^{\mathrm{HC}}_2(f)=1\f$;
 * The expected optimal filter in this case is
 * \f$\widetilde{Q}(f)\propto f^{-3}\f$.</li>
 * <li> \f$\gamma(f) = P^{\mathrm{C}}_1(f) = P^{\mathrm{C}}_2(f) = P^{\mathrm{HC}}_1(f)
 * = P^{\mathrm{HC}}_2(f)=1\f$;
 * \f$h_{100}^2\Omega_{\mathrm{GW}}(f)=f^3\f$.
 * The expected optimal filter in this case is
 * \f$\widetilde{Q}(f)=\textrm{constant}\f$.</li>
 * </ol>
 *
 * ### Uses ###
 *
 * \code
 * LALStochasticOptimalFilter()
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
#define STOCHASTICOPTIMALFILTERTESTC_ENOM 0	/**< Nominal exit */
#define STOCHASTICOPTIMALFILTERTESTC_EARG 1	/**< Error parsing command-line arguments */
#define STOCHASTICOPTIMALFILTERTESTC_ECHK 2	/**< Error checking failed to catch bad data */
#define STOCHASTICOPTIMALFILTERTESTC_EFLS 3	/**< Incorrect answer for valid data */
#define STOCHASTICOPTIMALFILTERTESTC_EUSE 4	/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define STOCHASTICOPTIMALFILTERTESTC_MSGENOM "Nominal exit"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEARG "Error parsing command-line arguments"
#define STOCHASTICOPTIMALFILTERTESTC_MSGECHK "Error checking failed to catch bad data"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEFLS "Incorrect answer for valid data"
#define STOCHASTICOPTIMALFILTERTESTC_MSGEUSE "Bad user-entered data"

#define LAL_USE_OLD_COMPLEX_STRUCTS
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

#define STOCHASTICOPTIMALFILTERTESTC_TRUE     1
#define STOCHASTICOPTIMALFILTERTESTC_FALSE    0
#define STOCHASTICOPTIMALFILTERTESTC_DELTAF   1.0
#define STOCHASTICOPTIMALFILTERTESTC_F0       0.0
#define STOCHASTICOPTIMALFILTERTESTC_FREF     1.0
#define STOCHASTICOPTIMALFILTERTESTC_LENGTH   8
#define STOCHASTICOPTIMALFILTERTESTC_TOL      1e-6

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = STOCHASTICOPTIMALFILTERTESTC_FALSE;
BOOLEAN optHetero = STOCHASTICOPTIMALFILTERTESTC_FALSE;
REAL8 optFRef      = 1.0;
UINT4 optLength    = 0.0;
CHAR  optOverlapFile[LALNameLength]   = "";
CHAR  optOmegaFile[LALNameLength]     = "";
CHAR  optInvNoise1File[LALNameLength] = "";
CHAR  optInvNoise2File[LALNameLength] = "";
CHAR  optHwInvNoise1File[LALNameLength] = "";
CHAR  optHwInvNoise2File[LALNameLength] = "";
CHAR  optOptimalFile[LALNameLength]     = "";

INT4 code;

static void Usage (const char *program, int exitflag);
static void ParseOptions (int argc, char *argv[]);
static REAL8 mu(const REAL4FrequencySeries*, const REAL4FrequencySeries*,
                const COMPLEX8FrequencySeries*);

int main(int argc, char *argv[])
{

  static LALStatus         status;

  StochasticOptimalFilterInput             input;

  StochasticOptimalFilterNormalizationInput  normIn;
  StochasticOptimalFilterNormalizationOutput normOut;
  StochasticOptimalFilterNormalizationParameters normParams;

  LIGOTimeGPS              epoch = {1,0};

  REAL4FrequencySeries     overlap;
  REAL4FrequencySeries     omegaGW;
  REAL4FrequencySeries     invNoise1;
  REAL4FrequencySeries     invNoise2;
  COMPLEX8FrequencySeries  hcInvNoise1;
  COMPLEX8FrequencySeries  hcInvNoise2;
  COMPLEX8FrequencySeries  optimal;

  REAL8                    omegaRef;

  INT4       i;
  REAL8      f;
  REAL8      muTest;  /*refers to the calculated mu to test normalization*/
  REAL8      testNum;     /*temporary value used to check optimal output */

  LALUnitPair              unitPair;
  LALUnit                  expectedUnit;
  BOOLEAN                  result;

  CHARVector               *unitString = NULL;

  REAL4WithUnits           lambda;


  ParseOptions (argc, argv);

  /* define valid params */
  overlap.name[0] = '\0';
  overlap.f0     = STOCHASTICOPTIMALFILTERTESTC_F0;
  overlap.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;
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

  lambda.value = 0;
  lambda.units = lalDimensionlessUnit;

  strncpy(overlap.name, "", LALNameLength);
  hcInvNoise1.name[0] = '\0';
  hcInvNoise1.f0     = overlap.f0;
  hcInvNoise1.deltaF = overlap.deltaF;
  hcInvNoise1.epoch  = overlap.epoch;
  hcInvNoise1.data   = NULL;
  hcInvNoise1.sampleUnits = lalDimensionlessUnit;
#ifndef LAL_NDEBUG
  COMPLEX8FrequencySeries  complexBadData = optimal = hcInvNoise2 = hcInvNoise1;
#endif

  hcInvNoise1.sampleUnits.unitNumerator[LALUnitIndexStrain] = -1;
  hcInvNoise1.sampleUnits.unitNumerator[LALUnitIndexADCCount] = -1;
  hcInvNoise1.sampleUnits.unitNumerator[LALUnitIndexSecond] = -1;
  hcInvNoise1.sampleUnits.powerOfTen = 18;
  hcInvNoise2.sampleUnits = hcInvNoise1.sampleUnits;

  /* allocate memory */
  LALSCreateVector(&status, &(overlap.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALSCreateVector(&status, &(omegaGW.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise1.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALSCreateVector(&status, &(invNoise2.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  LALCCreateVector(&status, &(hcInvNoise1.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(hcInvNoise2.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }
  LALCCreateVector(&status, &(optimal.data),
                   STOCHASTICOPTIMALFILTERTESTC_LENGTH);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  input.overlapReductionFunction = &overlap;
  input.omegaGW = &omegaGW;
  input.halfCalibratedInverseNoisePSD1 = &hcInvNoise1;
  input.halfCalibratedInverseNoisePSD2 = &hcInvNoise2;

  /* TEST INVALID DATA HERE -------------------------------------- */
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    /* test behavior for null pointer to input structure */
    LALStochasticOptimalFilter(&status, &optimal, NULL, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to input structure results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to output series */
    LALStochasticOptimalFilter(&status, NULL, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to output series results in error:\n \"%s\"\n",STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to overlap reduction function */
    input.overlapReductionFunction = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to overlap reduction function results in error:\n      \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to gravitational-wave spectrum */
    input.omegaGW = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to gravitational-wave spectrum results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* test behavior for null pointer to half-calibrated inverse noise 1 */
    /* of input structure */

    input.halfCalibratedInverseNoisePSD1 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to first half-calibrated inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD1 = &hcInvNoise1;

    /* test behavior for null pointer to second half-calibrated inverse */
    /* noise PSD member of input structure */
    input.halfCalibratedInverseNoisePSD2 = NULL;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to second half-calibrated inverse noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD2 = &hcInvNoise2;

    /* test behavior for null pointer to data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* test behavior for null pointer to data member of half-calibrated */
    /* first calibrated inverse noise PSD */
    input.halfCalibratedInverseNoisePSD1 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of first half-calibrated noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD1 = &hcInvNoise1;

    /* test behavior for null pointer to data member of half-calibrated */
    /* second calibrated inverse noise PSD */
    input.halfCalibratedInverseNoisePSD2 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of second half-calibrated noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD2 = &hcInvNoise2;

    /* test behavior for null pointer to data member of output */
    /* frequency series */
    LALStochasticOptimalFilter(&status, &complexBadData, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);


    /* Create a vector for testing REAL4 null data-data pointers */
    LALSCreateVector(&status, &(realBadData.data), STOCHASTICOPTIMALFILTERTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    REAL4*                   realTempPtr;
    realTempPtr = realBadData.data->data;
    realBadData.data->data = NULL;

    /* test behavior for null pointer to data member of data member of */
    /* output frequency series */

    LALStochasticOptimalFilter(&status, &complexBadData, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* test behavior for null pointer to data member of data member of overlap */
    input.overlapReductionFunction = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of overlap reduction function results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.overlapReductionFunction = &overlap;

    /* test behavior for null pointer to data member of data member of omega */
    input.omegaGW = &realBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of gravitational-wave spectrum results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.omegaGW = &omegaGW;

    /* Create a vector for testing COMPLEX8 null data-data pointers */
    LALCCreateVector(&status, &(complexBadData.data), STOCHASTICOPTIMALFILTERTESTC_LENGTH);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    COMPLEX8*                complexTempPtr;
    complexTempPtr = complexBadData.data->data;
    complexBadData.data->data = NULL;

    /* test behavior for null pointer to data member of data member of */
    /* first half-calibrated inverse noise PSD */
    input.halfCalibratedInverseNoisePSD1 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of first half-calibrated noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD1 = &hcInvNoise1;

    /* test behavior for null pointer to data member of data member of */
    /* second half-calibrated inverse noise PSD */
    input.halfCalibratedInverseNoisePSD2 = &complexBadData;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENULLPTR,
                              STOCHASTICCROSSCORRELATIONH_MSGENULLPTR,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: null pointer to data member of data member of second half-calibrated noise PSD results in error:\n       \"%s\"\n", STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    input.halfCalibratedInverseNoisePSD2 = &hcInvNoise2;

    /** clean up **/
    realBadData.data->data = realTempPtr;
    LALSDestroyVector(&status, &(realBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    complexBadData.data->data = complexTempPtr;
    LALCDestroyVector(&status, &(complexBadData.data));
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* test behavior for zero length */
    overlap.data->length =
      omegaGW.data->length = invNoise1.data->length =
      invNoise2.data->length =
      hcInvNoise1.data->length =
      hcInvNoise2.data->length =
      optimal.data->length = 0;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EZEROLEN,
                              STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero length results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    /* reassign valid length */
    overlap.data->length =
      omegaGW.data->length = invNoise1.data->length =
      invNoise2.data->length =
      hcInvNoise1.data->length =
      hcInvNoise2.data->length =
      optimal.data->length =STOCHASTICOPTIMALFILTERTESTC_LENGTH;

    /* test behavior for negative frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF =
      hcInvNoise1.deltaF =
      hcInvNoise2.deltaF = -3.5;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
                              STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF =
      hcInvNoise1.deltaF =
      hcInvNoise2.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;

    /* test behavior for zero frequency spacing */
    overlap.deltaF = omegaGW.deltaF =
      invNoise1.deltaF = invNoise2.deltaF =
      hcInvNoise1.deltaF =
      hcInvNoise2.deltaF = 0;
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
                              STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF,
                              STOCHASTICOPTIMALFILTERTESTC_ECHK,
                              STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
    printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
           STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
    /* reassign valid frequency spacing */
    overlap.deltaF =
      omegaGW.deltaF = invNoise1.deltaF =
      invNoise2.deltaF =
      hcInvNoise1.deltaF =
      hcInvNoise2.deltaF =  STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

  /* test behavior for negative start frequency */
  overlap.f0 = omegaGW.f0
    = invNoise1.f0
    = invNoise2.f0
    = hcInvNoise1.f0
    = hcInvNoise2.f0 = -3.0;

  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);

  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: negative start frequency results in error:\n      \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  overlap.f0 = omegaGW.f0 =
    invNoise1.f0 = invNoise2.f0 =
    hcInvNoise1.f0 =
    hcInvNoise2.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;

  /* test behavior length mismatch between overlap reduction function and output */
  optimal.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and output series results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  optimal.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;

  /* test behavior for length mismatch between overlap reduction function and omega */
  omegaGW.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: length mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  omegaGW.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;

  /* test behavior length mismatch between overlap reduction function and half-calibrated */
  /* first calibrated inverse noise PSD */
  hcInvNoise1.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and first half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  hcInvNoise1.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;

  /* test behavior length mismatch between overlap reduction function and half-calibrated */
  /* first calibrated inverse noise PSD */
  hcInvNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH - 1;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMLEN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMLEN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: length mismatch between overlap reduction function and second half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  hcInvNoise2.data->length = STOCHASTICOPTIMALFILTERTESTC_LENGTH;

  /* test behavior for frequency spacing mismatch between overlap reduction function and omega */
  omegaGW.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  omegaGW.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;

  /* test behavior frequency spacing mismatch between overlap reduction function and */
  /* first half-calibrated inverse noise PSD */
  hcInvNoise1.deltaF = 2.0 * STOCHASTICOPTIMALFILTERTESTC_DELTAF;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and first half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  hcInvNoise1.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;

  /* test behavior frequency spacing mismatch between overlap reduction function and */
  /* second half-calibrated inverse noise PSD */
  hcInvNoise2.deltaF = 305;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
    {
      return code;
    }
  printf("  PASS: frequency spacing mismatch between overlap reduction function and second half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  hcInvNoise2.deltaF = STOCHASTICOPTIMALFILTERTESTC_DELTAF;

  /* test behavior for start frequency mismatch between overlap reduction function and omega */
  omegaGW.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and gravitational-wave spectrum results in error:\n         \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  omegaGW.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;

  /* test behavior start frequency mismatch between overlap reduction function and half-calibrated */
  /* first calibrated inverse noise PSD */
  hcInvNoise1.f0 = 30;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and first half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  hcInvNoise1.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;

  /* test behavior start frequency mismatch between overlap reduction function and half-calibrated */
  /* second calibrated inverse noise PSD */
  hcInvNoise2.f0 = 305;
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status, STOCHASTICCROSSCORRELATIONH_EMMFMIN,
                            STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN,
                            STOCHASTICOPTIMALFILTERTESTC_ECHK,
                            STOCHASTICOPTIMALFILTERTESTC_MSGECHK) ) )
  {
    return code;
  }
  printf("  PASS: start frequency mismatch between overlap reduction function and second half-calibrated inverse noise PSD results in error:\n       \"%s\"\n",
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  hcInvNoise2.f0 = STOCHASTICOPTIMALFILTERTESTC_F0;

 /* VALID TEST DATA HERE ----------------------------------------- */

  /* create input to test */
  overlap.data->data[0] = 1;
  omegaGW.data->data[0] = 0;
  invNoise1.data->data[0] = 0;
  invNoise2.data->data[0] = 0;
  hcInvNoise1.data->data[0].realf_FIXME = 0;
  hcInvNoise1.data->data[0].imagf_FIXME = 0;
  hcInvNoise2.data->data[0].realf_FIXME = 0;
  hcInvNoise2.data->data[0].imagf_FIXME = 0;

  /** Test 1 **/
  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;

    overlap.data->data[i] = 1;
    omegaGW.data->data[i] = 1;
    invNoise1.data->data[i] = 1;
    invNoise2.data->data[i] = 1;
    hcInvNoise1.data->data[i].realf_FIXME = 1;
    hcInvNoise1.data->data[i].imagf_FIXME = 0;
    hcInvNoise2.data->data[i].realf_FIXME = 1;
    hcInvNoise2.data->data[i].imagf_FIXME = 0;
  }

  /* fill normalization output */
  normOut.variance = NULL;
  normOut.normalization = &lambda;

  /* fill normalization input */
  normIn.overlapReductionFunction     = &(overlap);
  normIn.omegaGW                      = &(omegaGW);
  normIn.inverseNoisePSD1             = &(invNoise1);
  normIn.inverseNoisePSD2             = &(invNoise2);

  /* fill normalziation parameters */
  normParams.fRef        = STOCHASTICOPTIMALFILTERTESTC_FREF;
  normParams.heterodyned = STOCHASTICOPTIMALFILTERTESTC_FALSE;
  normParams.window1     = NULL;
  normParams.window2     = NULL;

  /* calculate normalization */
  LALStochasticOptimalFilterNormalization(&status, &normOut,
                                          &normIn, &normParams);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /* fill optimal input */
  input.overlapReductionFunction     = &(overlap);
  input.omegaGW                      = &(omegaGW);
  input.halfCalibratedInverseNoisePSD1 = &(hcInvNoise1);
  input.halfCalibratedInverseNoisePSD2 = &(hcInvNoise2);

  /* calculate optimal filter */
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output f0 */
  if (optVerbose)
  {
    printf("f0=%g, should be %g\n", optimal.f0,
           STOCHASTICOPTIMALFILTERTESTC_F0);
  }
  if ( fabs(optimal.f0-STOCHASTICOPTIMALFILTERTESTC_F0)
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output deltaF */
  if (optVerbose)
  {
    printf("deltaF=%g, should be %g\n", optimal.deltaF,
           STOCHASTICOPTIMALFILTERTESTC_DELTAF);
  }
  if ( fabs(optimal.deltaF-STOCHASTICOPTIMALFILTERTESTC_DELTAF)
       / STOCHASTICOPTIMALFILTERTESTC_DELTAF
        > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output units */

  expectedUnit = lalDimensionlessUnit;
  expectedUnit.unitNumerator[LALUnitIndexADCCount] = -2;
  unitPair.unitOne = &expectedUnit;
  unitPair.unitTwo = &(optimal.sampleUnits);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output values */
  if (optVerbose)
  {
    printf("Q(0)=%g + %g i, should be 0\n",
           crealf(optimal.data->data[0]), cimagf(optimal.data->data[0]));
  }
  if ( fabs(crealf(optimal.data->data[0])) > STOCHASTICOPTIMALFILTERTESTC_TOL
       || fabs(cimagf(optimal.data->data[0]))
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #1\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    testNum = (STOCHASTICOPTIMALFILTERTESTC_DELTAF
               *STOCHASTICOPTIMALFILTERTESTC_DELTAF
               *STOCHASTICOPTIMALFILTERTESTC_DELTAF)
      /(f*f*f);
    if (optVerbose)
    {
      printf("Q(%g Hz)/Re(Q(%g Hz))=%g + %g i, should be %g\n",
             f, STOCHASTICOPTIMALFILTERTESTC_DELTAF,
             crealf(optimal.data->data[i])/crealf(optimal.data->data[1]),
             cimagf(optimal.data->data[i])/crealf(optimal.data->data[1]),
             testNum);
    }
    if (fabs(crealf(optimal.data->data[i])/crealf(optimal.data->data[1])
             - testNum)/testNum
        > STOCHASTICOPTIMALFILTERTESTC_TOL
        || fabs(cimagf(optimal.data->data[i])/crealf(optimal.data->data[1]))
        > STOCHASTICOPTIMALFILTERTESTC_TOL)
    {
      printf("  FAIL: Valid data test #1\n");
      if (optVerbose)
      {
        printf("Exiting with error: %s\n",
               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
      }
      return STOCHASTICOPTIMALFILTERTESTC_EFLS;
    }
  }

  /* normalization costant */
  omegaRef = 1.0;
  muTest = mu(&omegaGW, &overlap, &optimal);
  if (optVerbose)
  {
      printf("mu=%g, should be %g\n", muTest, omegaRef);
  }
  if ( fabs(muTest - omegaRef)/omegaRef > STOCHASTICOPTIMALFILTERTESTC_TOL)
  {
     printf("  FAIL: Valid data test #1\n" );
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  printf("  PASS: Valid data test #1\n");

  /** Test 2 **/
  for (i=0; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;

    overlap.data->data[i] = 1;
    omegaGW.data->data[i] = pow(f,3);
    invNoise1.data->data[i] = 1;
    invNoise2.data->data[i] = 1;
    hcInvNoise1.data->data[i].realf_FIXME = 1;
    hcInvNoise1.data->data[i].imagf_FIXME = 0;
    hcInvNoise2.data->data[i].realf_FIXME = 1;
    hcInvNoise2.data->data[i].imagf_FIXME = 0;
  }

  /* fill normalization output */
  normOut.variance = NULL;
  normOut.normalization = &lambda;

  /* fill normalization input */
  normIn.overlapReductionFunction     = &(overlap);
  normIn.omegaGW                      = &(omegaGW);
  normIn.inverseNoisePSD1             = &(invNoise1);
  normIn.inverseNoisePSD2             = &(invNoise2);

  /* fill normalziation parameters */
  normParams.fRef = STOCHASTICOPTIMALFILTERTESTC_FREF;
  normParams.heterodyned = STOCHASTICOPTIMALFILTERTESTC_FALSE;

  /* calculate normalization */
  LALStochasticOptimalFilterNormalization(&status, &normOut,
                                          &normIn, &normParams);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /* fill optimal input */
  input.overlapReductionFunction     = &(overlap);
  input.omegaGW                      = &(omegaGW);
  input.halfCalibratedInverseNoisePSD1 = &(hcInvNoise1);
  input.halfCalibratedInverseNoisePSD2 = &(hcInvNoise2);

  /* calculate optimal filter */
  LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
  if ( ( code = CheckStatus(&status,0, "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  /******   check output   ******/

  /* check output f0 */
  if (optVerbose)
  {
    printf("f0=%g, should be %g\n", optimal.f0,
           STOCHASTICOPTIMALFILTERTESTC_F0);
  }
  if ( fabs(optimal.f0-STOCHASTICOPTIMALFILTERTESTC_F0)
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n",
             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output deltaF */
  if (optVerbose)
  {
    printf("deltaF=%g, should be %g\n", optimal.deltaF,
           STOCHASTICOPTIMALFILTERTESTC_DELTAF);
  }
  if ( fabs(optimal.deltaF-STOCHASTICOPTIMALFILTERTESTC_DELTAF)
       / STOCHASTICOPTIMALFILTERTESTC_DELTAF
        > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output units */
  expectedUnit = lalDimensionlessUnit;
  expectedUnit.unitNumerator[LALUnitIndexADCCount] = -2;
  unitPair.unitOne = &expectedUnit;
  unitPair.unitTwo = &(optimal.sampleUnits);
  LALUnitCompare(&status, &result, &unitPair);
  if ( ( code = CheckStatus(&status, 0 , "",
                            STOCHASTICOPTIMALFILTERTESTC_EFLS,
                            STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
    return code;
  }

  if (optVerbose)
  {
    LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALUnitAsString( &status, unitString, unitPair.unitTwo );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "Units are \"%s\", ", unitString->data );

    LALUnitAsString( &status, unitString, unitPair.unitOne );
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    printf( "should be \"%s\"\n", unitString->data );

    LALCHARDestroyVector(&status, &unitString);
    if ( ( code = CheckStatus(&status, 0 , "",
                              STOCHASTICOPTIMALFILTERTESTC_EFLS,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
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
             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  /* check output values */
  if (optVerbose)
  {
    printf("Q(0)=%g + %g i, should be 0\n",
           crealf(optimal.data->data[0]), cimagf(optimal.data->data[0]));
  }
  if ( fabs(crealf(optimal.data->data[0])) > STOCHASTICOPTIMALFILTERTESTC_TOL
       || fabs(cimagf(optimal.data->data[0]))
       > STOCHASTICOPTIMALFILTERTESTC_TOL )
  {
    printf("  FAIL: Valid data test #2\n");
    if (optVerbose)
    {
      printf("Exiting with error: %s\n", STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
    }
    return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  testNum = 1.0;
  for (i=1; i < STOCHASTICOPTIMALFILTERTESTC_LENGTH; i++)
  {
    f = i*STOCHASTICOPTIMALFILTERTESTC_DELTAF;
    if (optVerbose)
    {
      printf("Q(%g Hz)/Re(Q(%g Hz))=%g + %g i, should be %g\n",
             f, STOCHASTICOPTIMALFILTERTESTC_DELTAF,
             crealf(optimal.data->data[i])/crealf(optimal.data->data[1]),
             cimagf(optimal.data->data[i])/crealf(optimal.data->data[1]),
             testNum);
    }
    if (fabs(crealf(optimal.data->data[i])/crealf(optimal.data->data[1])
             - testNum)/testNum
        > STOCHASTICOPTIMALFILTERTESTC_TOL
        || fabs(cimagf(optimal.data->data[i])/crealf(optimal.data->data[1]))
        > STOCHASTICOPTIMALFILTERTESTC_TOL)
    {
      printf("  FAIL: Valid data test #2\n");
      if (optVerbose)
      {
        printf("Exiting with error: %s\n",
               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS);
      }
      return STOCHASTICOPTIMALFILTERTESTC_EFLS;
    }
  }

  /* normalization costant */
  omegaRef = normParams.fRef * normParams.fRef * normParams.fRef;
  muTest = mu(&omegaGW, &overlap, &optimal);
  if (optVerbose)
  {
      printf("mu=%g, should be %g\n", muTest, omegaRef);
  }
  if ( fabs(muTest - omegaRef)/omegaRef > STOCHASTICOPTIMALFILTERTESTC_TOL)
  {
     printf("  FAIL: Valid data test #2\n" );
     return STOCHASTICOPTIMALFILTERTESTC_EFLS;
  }

  printf("  PASS: Valid data test #2\n");

  /* clean up valid data */
  LALSDestroyVector(&status, &(overlap.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(omegaGW.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise1.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALSDestroyVector(&status, &(invNoise2.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(hcInvNoise1.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(hcInvNoise2.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
       }
  LALCDestroyVector(&status, &(optimal.data));
  if ( ( code = CheckStatus (&status, 0 , "",
                             STOCHASTICOPTIMALFILTERTESTC_EFLS,
                             STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
  {
       return code;
  }
  LALCheckMemoryLeaks();

 printf("PASS: all tests\n");


 /* VALID USER TEST DATA HERE ----------------------------------------- */

  if (optOverlapFile[0] &&  optOmegaFile[0] && optInvNoise1File[0] && optInvNoise2File[0]
      && optHwInvNoise1File[0] && optHwInvNoise2File[0] && optOptimalFile[0])
  {

    /* allocate memory */
    overlap.data     = NULL;
    omegaGW.data     = NULL;
    invNoise1.data   = NULL;
    invNoise2.data   = NULL;
    hcInvNoise1.data = NULL;
    hcInvNoise2.data = NULL;
    optimal.data     = NULL;

    LALSCreateVector(&status, &(overlap.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    LALSCreateVector(&status, &(omegaGW.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise1.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSCreateVector(&status, &(invNoise2.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(hcInvNoise1.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(hcInvNoise2.data), optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCCreateVector(&status, &(optimal.data),optLength);
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }

    /* Read input files */
    LALSReadFrequencySeries(&status, &overlap,   optOverlapFile);
    LALSReadFrequencySeries(&status, &omegaGW,   optOmegaFile);
    LALSReadFrequencySeries(&status, &invNoise1, optInvNoise1File);
    LALSReadFrequencySeries(&status, &invNoise2, optInvNoise2File);
    LALCReadFrequencySeries(&status, &hcInvNoise1, optHwInvNoise1File);
    LALCReadFrequencySeries(&status, &hcInvNoise2, optHwInvNoise2File);


    /* fill normalization output */
    normOut.variance = NULL;
    normOut.normalization = &lambda;

    /* fill normalization input */
    normIn.overlapReductionFunction     = &(overlap);
    normIn.omegaGW                      = &(omegaGW);
    normIn.inverseNoisePSD1             = &(invNoise1);
    normIn.inverseNoisePSD2             = &(invNoise2);

    /* fill normalziation parameters */
    normParams.fRef = optFRef;
    normParams.heterodyned = optHetero;

    /* calculate normalization */
    LALStochasticOptimalFilterNormalization(&status, &normOut,
                                            &normIn, &normParams);
    if ( ( code = CheckStatus(&status,0, "",
                              STOCHASTICOPTIMALFILTERTESTC_EUSE,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEUSE) ) )
    {
      return code;
    }

    /* fill optimal input */
    input.overlapReductionFunction     = &(overlap);
    input.omegaGW                      = &(omegaGW);
    input.halfCalibratedInverseNoisePSD1 = &(hcInvNoise1);
    input.halfCalibratedInverseNoisePSD2 = &(hcInvNoise2);

    /* calculate optimal filter */
    LALStochasticOptimalFilter(&status, &optimal, &input, &lambda);
    if ( ( code = CheckStatus(&status,0, "",
                              STOCHASTICOPTIMALFILTERTESTC_EUSE,
                              STOCHASTICOPTIMALFILTERTESTC_MSGEUSE) ) )
    {
      return code;
    }


    LALCPrintFrequencySeries(&optimal, optOptimalFile);
    printf("=========== Optimal Filter Written to File %s ===========\n",
           optOptimalFile);

    /* clean up */
    LALSDestroyVector(&status, &(overlap.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(omegaGW.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise1.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALSDestroyVector(&status, &(invNoise2.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(hcInvNoise1.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(hcInvNoise2.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
    LALCDestroyVector(&status, &(optimal.data));
    if ( ( code = CheckStatus (&status, 0 , "",
                               STOCHASTICOPTIMALFILTERTESTC_EFLS,
                               STOCHASTICOPTIMALFILTERTESTC_MSGEFLS) ) )
    {
      return code;
    }
  }
  LALCheckMemoryLeaks();
  return 0;

}

/*----------------------------------------------------------------------*/
static REAL8 mu(const REAL4FrequencySeries* omegaGW,
                const REAL4FrequencySeries* overlap,
                const COMPLEX8FrequencySeries* optimal)
{
  INT4       i;
  REAL8      f;
  REAL8      constant;
  REAL8      f3;
  REAL8      deltaF;
  INT4       length;
  REAL8      muval = 0.0;

  deltaF = omegaGW->deltaF;
  length = omegaGW->data->length;

  constant = deltaF * 3.0L * ( (LAL_H0FAC_SI*1.0e+18) * (LAL_H0FAC_SI*1.0e+18) )
             / ( 20.0L * (LAL_PI*LAL_PI) );

  /* calculate mu */
  for(i=1; i < (length);i++)
    {
      f = i*deltaF;
      f3 = f*f*f;
      muval += (2.0  * constant * (omegaGW->data->data[i]) *
            (overlap->data->data[i]) * (crealf(optimal->data->data[i])))/f3;
    }
    return muval;
}

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
  fprintf (stderr, "  -i filename    read first calibrated inverse noise PSD from file filename\n");
  fprintf (stderr, "  -j filename    read second calibrated inverse noise PSD from file filename\n");
  fprintf (stderr, "  -s filename    read first half-calibrated inverse noise PSD from file filename\n");
  fprintf (stderr, "  -t filename    read second half-calibrated inverse noise PSD from file filename\n");
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

      case 's': /* specify hcInvNoise1 file */
        strncpy (optHwInvNoise1File, optarg, LALNameLength);
        break;

      case 't': /* specify hcInvNoise2 file */
        strncpy (optHwInvNoise2File, optarg, LALNameLength);
        break;

      case 'd': /* set debug level */
        break;

      case 'v': /* optVerbose */
        optVerbose = STOCHASTICOPTIMALFILTERTESTC_TRUE;
        break;

      case 'y': /* optHetero */
        optHetero = STOCHASTICOPTIMALFILTERTESTC_TRUE;
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
