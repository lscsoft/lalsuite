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

#include <complex.h>
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

#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/AVFactories.h>
#include <lal/ReadFTSeries.h>
#include <lal/PrintFTSeries.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>

#include "CheckStatus.h"

/**
 * \file
 * \ingroup CoarseGrainFrequencySeries_h
 *
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 *
 * \brief Test suite for <tt>LALCCoarseGrainFrequencySeries()</tt>.
 *
 * ### Usage ###
 *
 * \code
 * ./CCoarseGrainFrequencySeriesTest
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * -i filename    read fine grained series from file filename
 * -o filename    print coarse grained  series to file filename
 * -n length      input series contains length points
 * -m length      output series contains length points
 * -e deltaF      set coarse grained frequency spacing to deltaF
 * -f f0          set start frequency of output to f0
 * \endcode
 *
 * ### Description ###
 *
 * This program tests the routine
 * <tt>LALCCoarseGrainFrequencySeries()</tt>, which coarse-grains a
 * frequency series.
 *
 * First, it tests that the correct error codes
 * (cf \ref CoarseGrainFrequencySeries_h
 * are generated for the following error conditions (tests in
 * \e italics are not performed if \c LAL_NDEBUG is set, as
 * the corresponding checks in the code are made using the ASSERT macro):
 * <ul>
 * <li> <em>null pointer to output series</em></li>
 * <li> <em>null pointer to input series</em></li>
 * <li> <em>null pointer to data member of output series</em></li>
 * <li> <em>null pointer to data member of input series</em></li>
 * <li> <em>null pointer to data member of data member of input series</em></li>
 * <li> <em>null pointer to data member of data member of output series</em>
 * <li> <em>zero length</em></li>
 * <li> <em>negative frequency spacing</em></li>
 * <li> <em>zero frequency spacing</em></li>
 * </ul>
 *
 * It then verifies that the correct values are obtained for some simple
 * test cases
 * <ul>
 * <li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
 * f=2\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{1/2,2,4,6\}\f$.</li>
 * <li> \f$\{h_\ell'\}=\{0,1,2,3,4,5,6,7\}\f$, \f$f'_0=f_0\f$, \f$\delta f'=1\f$, \f$\delta
 * f=3\f$, \f$N=3\f$; the expected output is \f$\{h_k\}=\{2/3,3,6\}\f$.</li>
 * <li> \f$f_0'=40\f$, \f$\delta f'= 1\f$,
 * \f$\{h_k\}=\{f_k+i\,f_k^{-1}|k=0,\ldots,4\}\f$, \f$f_0=41\f$,
 * \f$f_0=f'_0\f$, \f$\delta f=2\f$.
 * \f$\delta f'=3\f$, \f$N=\f$ ; the expected output is
 * \f[
 * \{h'_\ell\}=\left\{41+i\left(\frac{1}{40}+\frac{2}{41}+\frac{1}{42}\right),
 * 43+i\left(\frac{1}{42}+\frac{2}{43}+\frac{1}{44}\right)
 * \right\}
 * \f]
 * % '</li>
 * </ul>
 * For each successful test (both of these valid data and the invalid
 * ones described above), it prints \c PASS to standard output;
 * if a test fails, it prints \c FAIL.
 *
 * If the \c filename arguments are present, it also reads a
 * frequency series from a file, calls
 * <tt>LALCCoarseGrainFrequencySeries()</tt>, and writes the results to
 * the specified output file.
 *
 * ### Notes ###
 *
 * <ul>
 * <li> In addition to the error checks tested in this routine, the
 * function checks for errors related to inconsistency of coarse
 * graining parameters, as well as duplicate input and output pointers.
 * Tests of these error checks are still to be
 * added to this test program.</li>
 * <li> No specific error checking is done on user-specified data.  If
 * \c length is missing, the resulting default will cause a bad
 * data error.</li>
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
#define CCOARSEGRAINFREQUENCYSERIESTESTC_ENOM 0		/**< Nominal exit */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EARG 1		/**< Error parsing command-line arguments */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK 2		/**< Error checking failed to catch bad data */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS 3		/**< Incorrect answer for valid data */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE 4		/**< Bad user-entered data */
/*@}*/

/** \cond DONT_DOXYGEN */
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGENOM "Nominal exit"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEARG "Error parsing command-line arguments"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK "Error checking failed to catch bad data"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS "Incorrect answer for valid data"
#define CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE "Bad user-entered data"




#define CCOARSEGRAINFREQUENCYSERIESTESTC_TOL           1e-6

#define CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL           1e-5

#define CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC      1234
#define CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS       56789

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0    1.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F00        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0    8

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1    2.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F01        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1    4

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2    3.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F02        0.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2    3

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3    1.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F03        40.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3    5

#define CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4    2.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_F04        41.0
#define CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4    2

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TSDT     (1.0 / 1024.0)
#define CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN     16384

#define CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN     50
#define CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX     500

#define CCOARSEGRAINFREQUENCYSERIESTESTC_RESRATIO 9.14

#define CCOARSEGRAINFREQUENCYSERIESTESTC_TRUE     1
#define CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE    0

extern char *optarg;
extern int   optind;

BOOLEAN optVerbose = CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE;
UINT4 optInLength    = 0;
UINT4 optOutLength   = 0;
REAL8 optDeltaF     = -1.0;
REAL8 optF0       = 0.0;

INT4 optSeed = 1682;

CHAR optInputFile[LALNameLength] = "";
CHAR optOutputFile[LALNameLength] = "";
INT4 code;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

int
main( int argc, char *argv[] )
{

   static LALStatus         status;

   UINT4      i;
   REAL8      f;

   const COMPLEX8  testInputDataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0]
     = {0.0, 1.0, 2.0, 3.0,
        4.0, 5.0, 6.0, 7.0};

   const COMPLEX8
     expectedOutput1DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1]
     = {0.5, 2.0, 4.0, 6.0};

   const COMPLEX8
     expectedOutput2DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2]
     = {(2.0/3.0), 3.0, 6.0};

   const COMPLEX8
     testInput3DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3]
     = {40.0 + I * 1.0/40.0, 41.0 + I * 1.0/41.0, 42.0 + I * 1.0/42.0,
        43.0 + I * 1.0/43.0, 44.0 + I * 1.0/44.0};

   const COMPLEX8
     expectedOutput4DataData[CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4]
     = {41.0 + I * (1.0/40.0+2.0/41.0+1.0/42.0) / 4.0,
        43.0 + I * (1.0/42.0+2.0/43.0+1.0/44.0) / 4.0};

   COMPLEX8FrequencySeries             goodInput;
   COMPLEX8FrequencySeries     goodOutput;

   BOOLEAN                result;
   LALUnitPair            unitPair;

   CHARVector             *unitString;

   FrequencySamplingParams     params;

   COMPLEX8               coarseTotal, fineTotal;
   COMPLEX8               cError;
   REAL4TimeSeries        timeSeries;

   RealFFTPlan            *fftPlan;

   RandomParams           *randomParams;


   ParseOptions( argc, argv );

   /* TEST INVALID DATA HERE ------------------------------------------- */

   /* define valid parameters */
   goodInput.f0                   = CCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   goodInput.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   goodInput.epoch.gpsSeconds     = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC;
   goodInput.epoch.gpsNanoSeconds = CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS;
   goodInput.data                 = NULL;
   goodOutput.data                = NULL;

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F00;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

#ifndef LAL_NDEBUG
   COMPLEX8FrequencySeries badInput = goodInput;
   COMPLEX8FrequencySeries badOutput = goodOutput;
#endif

   /* allocate input and output vectors */
   LALCCreateVector(&status, &(goodInput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }
   LALCCreateVector(&status, &(goodOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

#ifndef LAL_NDEBUG
   if ( ! lalNoDebug )
   {
     /* test behavior for null pointer to output series */
     LALCCoarseGrainFrequencySeries(&status, NULL, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, NULL, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to parameter structure */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, NULL);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to parameter structure results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of output series */
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of input series */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     COMPLEX8                   *cPtr = badOutput.data->data;
     badOutput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = cPtr;
     LALCDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for null pointer to data member of data member of output series */
     LALCCreateVector(&status, &(badInput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     cPtr = badInput.data->data;
     badInput.data->data = NULL;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &badInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: null pointer to data member of data member of input series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badInput.data->data = cPtr;
     LALCDestroyVector(&status, &(badInput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for duplicate pointers */

     /* input and output series */
     LALCCoarseGrainFrequencySeries(&status, &goodInput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput = goodInput;
     badOutput.data = goodInput.data;

     /* data members of input and output series */
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

     badOutput.data = NULL;

     /* data members of data members of input and output series */
     LALCCreateVector(&status, &(badOutput.data),
		      CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     cPtr = badOutput.data->data;
     badOutput.data->data = goodInput.data->data;
     LALCCoarseGrainFrequencySeries(&status, &badOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
			       COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: duplicate pointers to data members of data members of input and output series results in error:\n       \"%s\"\n", COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);
     badOutput.data->data = cPtr;
     LALCDestroyVector(&status, &(badOutput.data));
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     /* test behavior for zero length */
     /* input */

     goodInput.data->length = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in input results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodInput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;
     goodOutput.data->length = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

     /* output */

     goodOutput.data->length = params.length = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status, COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
			       COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero length in output results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

     goodOutput.data->length = params.length
       = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0;

     /* test behavior for negative frequency spacing */
     goodInput.deltaF = params.deltaF
       = -CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: negative frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* test behavior for zero frequency spacing */
     goodInput.deltaF = params.deltaF = 0;
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus(&status,
			       COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
			       COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_ECHK,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGECHK) ) )
     {
       return code;
     }
     printf("  PASS: zero frequency spacing results in error:\n       \"%s\"\n",
	    COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

     /* reassign valid frequency spacing */
     goodInput.deltaF = params.deltaF
       = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF0;

   } /* if ( ! lalNoDebug ) */
#endif /* LAL_NDEBUG */

   LALCDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* TEST VALID DATA HERE --------------------------------------------- */

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F01;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1;

   /* allocate input and output vectors */
   LALCCreateVector(&status, &(goodOutput.data), CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input time-series parameters */
   strncpy(goodInput.name,"Dummy test data",LALNameLength);
   goodInput.sampleUnits  = lalDimensionlessUnit;

     /* fill input data */
   for (i=0; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH0; ++i)
   {
     goodInput.data->data[i] = testInputDataData[i];
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F01);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1
        > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            crealf(goodOutput.data->data[0]), cimagf(goodOutput.data->data[0]),
            crealf(expectedOutput1DataData[0]), cimagf(expectedOutput1DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput1DataData[0]))
        /creal(expectedOutput1DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput1DataData[0]))
        /creal(expectedOutput1DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #1\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH1; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF1;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              crealf(goodOutput.data->data[i]), cimagf(goodOutput.data->data[i]),
              crealf(expectedOutput1DataData[i]), cimagf(expectedOutput1DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput1DataData[i]))
          /creal(expectedOutput1DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput1DataData[i]))
          /creal(expectedOutput1DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         )
     {
       printf("  FAIL: Valid data test #1\n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }


   /*-------Test #2-------*/

   params.f0                      = CCOARSEGRAINFREQUENCYSERIESTESTC_F02;
   params.deltaF               = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
   params.length               = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2;

   LALCCreateVector(&status, &(goodOutput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F02);
   }
   if (goodOutput.f0)
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
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
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            crealf(goodOutput.data->data[0]), cimagf(goodOutput.data->data[0]),
            crealf(expectedOutput2DataData[0]), cimagf(expectedOutput2DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput2DataData[0]))
        /creal(expectedOutput2DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput2DataData[0]))
        /creal(expectedOutput2DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #2\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH2; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF2;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              crealf(goodOutput.data->data[i]), cimagf(goodOutput.data->data[i]),
              crealf(expectedOutput2DataData[i]), cimagf(expectedOutput2DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput2DataData[i]))
          /creal(expectedOutput2DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput2DataData[i]))
          /creal(expectedOutput2DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #2 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALCDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   /*-------Test #3-------*/

   goodInput.f0               = CCOARSEGRAINFREQUENCYSERIESTESTC_F03;
   goodInput.deltaF           = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF3;

   LALCCreateVector(&status, &(goodInput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* fill input data */
   for (i=0; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH3; ++i)
   {
     goodInput.data->data[i] = testInput3DataData[i];
   }

   params.f0                  = CCOARSEGRAINFREQUENCYSERIESTESTC_F04;
   params.deltaF              = CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
   params.length              = CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4;

   LALCCreateVector(&status, &(goodOutput.data),
                    CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* check output f0 */
   if (optVerbose)
   {
     printf("f0=%g, should be %g\n", goodOutput.f0,
            CCOARSEGRAINFREQUENCYSERIESTESTC_F04);
   }
   if ( fabs(goodOutput.f0-CCOARSEGRAINFREQUENCYSERIESTESTC_F04)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_F04 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output deltaF */
   if (optVerbose)
   {
     printf("deltaF=%g, should be %g\n", goodOutput.deltaF,
            CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4);
   }
   if ( fabs(goodOutput.deltaF-CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4)
        / CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4 > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n", CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output epoch */
   if (optVerbose)
   {
     printf("epoch=%d seconds, %d nanoseconds; should be %d seconds, %d nanoseconds\n",
            goodOutput.epoch.gpsSeconds, goodOutput.epoch.gpsNanoSeconds,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC,
            CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS);
   }
   if ( goodOutput.epoch.gpsSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHSEC
        || goodOutput.epoch.gpsNanoSeconds
        != CCOARSEGRAINFREQUENCYSERIESTESTC_EPOCHNS )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output units */
   unitPair.unitOne = &(goodInput.sampleUnits);
   unitPair.unitTwo = &(goodOutput.sampleUnits);
   LALUnitCompare(&status, &result, &unitPair);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   if (optVerbose)
   {
     unitString = NULL;
     LALCHARCreateVector(&status, &unitString, LALUnitTextSize);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }

     LALUnitAsString( &status, unitString, unitPair.unitTwo );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "Units are \"%s\", ", unitString->data );

     LALUnitAsString( &status, unitString, unitPair.unitOne );
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
     printf( "should be \"%s\"\n", unitString->data );

     LALCHARDestroyVector(&status, &unitString);
     if ( ( code = CheckStatus(&status, 0 , "",
			       CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
     {
       return code;
     }
   }

   if (!result)
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
              CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   /* check output values */
   if (optVerbose)
   {
     printf("hBarTilde(0)=%g + %g i, should be %g + %g i\n",
            crealf(goodOutput.data->data[0]), cimagf(goodOutput.data->data[0]),
            crealf(expectedOutput4DataData[0]), cimagf(expectedOutput4DataData[0]));
   }
   if ((fabs(creal(goodOutput.data->data[0] - expectedOutput4DataData[0]))
        /creal(expectedOutput4DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
        ||
       (fabs(cimag(goodOutput.data->data[0] - expectedOutput4DataData[0]))
        /creal(expectedOutput4DataData[0]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
   {
     printf("  FAIL: Valid data test #3\n");
     if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   for (i=1; i<CCOARSEGRAINFREQUENCYSERIESTESTC_LENGTH4; ++i)
   {
     f = i * CCOARSEGRAINFREQUENCYSERIESTESTC_DELTAF4;
     if (optVerbose)
     {
       printf("hBarTilde(%f Hz)=%g + %g i, should be %g + %g i\n", f,
              crealf(goodOutput.data->data[i]), cimagf(goodOutput.data->data[i]),
              crealf(expectedOutput4DataData[i]), cimagf(expectedOutput4DataData[i]));
     }
     if ((fabs(creal(goodOutput.data->data[i] - expectedOutput4DataData[i]))
          /creal(expectedOutput4DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
         ||
         (fabs(cimag(goodOutput.data->data[i] - expectedOutput4DataData[i]))
          /creal(expectedOutput4DataData[i]) > CCOARSEGRAINFREQUENCYSERIESTESTC_TOL )
       )
     {
       printf("  FAIL: Valid data test #3 \n");
       if (optVerbose)
       {
         printf("Exiting with error: %s\n",
                CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
       }
       return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
     }
   }

   /* clean up valid data */
   LALCDestroyVector(&status, &goodInput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &goodOutput.data);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   /*-------Test #4-------*/

   strncpy(timeSeries.name,"Random test data",LALNameLength);
   timeSeries.deltaT = CCOARSEGRAINFREQUENCYSERIESTESTC_TSDT;
   timeSeries.f0 = 0.0;
   timeSeries.sampleUnits = lalDimensionlessUnit;
   timeSeries.epoch = goodOutput.epoch;

   timeSeries.data = NULL;

   LALSCreateVector(&status, &(timeSeries.data),
		    CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   goodInput.data = NULL;

   LALCCreateVector(&status, &(goodInput.data),
		    ( CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN / 2 + 1 ) );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* construct plan */
   fftPlan = NULL;

   LALCreateForwardRealFFTPlan(&status, &fftPlan,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_TSLEN,
			       CCOARSEGRAINFREQUENCYSERIESTESTC_FALSE);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS ) ) )
   {
     return code;
   }

   randomParams = NULL;

   LALCreateRandomParams( &status, &randomParams, optSeed );

   /* fill time series with normal deviates */

   LALNormalDeviates( &status, timeSeries.data, randomParams );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALTimeFreqRealFFT( &status, &goodInput, &timeSeries, fftPlan );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   for ( i = 0 ;
	 i * goodInput.deltaF <= CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN ;
	 ++i )
   {
     goodInput.data->data[i] = 0.0;
   }

   for ( i = CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX / goodInput.deltaF ;
	 i < goodInput.data->length ;
	 ++i )
   {
     goodInput.data->data[i] = 0.0;
   }

   fineTotal = 0.0;
   for ( i = 0 ; i < goodInput.data->length ; ++i )
   {
     fineTotal += goodInput.data->data[i];
   }
   fineTotal *= goodInput.deltaF;

   params.deltaF
     = goodInput.deltaF * CCOARSEGRAINFREQUENCYSERIESTESTC_RESRATIO;
   params.f0  = CCOARSEGRAINFREQUENCYSERIESTESTC_FMIN - params.deltaF;
   params.length = (CCOARSEGRAINFREQUENCYSERIESTESTC_FMAX - params.f0)
                   / params.deltaF + 3;

   goodOutput.data = NULL;
   LALCCreateVector(&status, &(goodOutput.data),
		    params.length );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   /* coarse grain */
   LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
   if ( ( code = CheckStatus( &status, 0 , "",
			      CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   coarseTotal = 0.0;
   for ( i = 0 ; i < goodOutput.data->length ; ++i )
   {
     coarseTotal += goodOutput.data->data[i];
   }
   coarseTotal *= goodOutput.deltaF;

   cError = ( creal(coarseTotal - fineTotal) * creal(fineTotal)
		 + cimag(coarseTotal - fineTotal) * cimag(fineTotal) )
               / ( creal(fineTotal) * creal(fineTotal) + cimag(fineTotal) * cimag(fineTotal) );
   cError += I * ( cimag(coarseTotal - fineTotal) * creal(fineTotal)
		 - creal(coarseTotal - fineTotal) * cimag(fineTotal) )
               / ( creal(fineTotal) * creal(fineTotal) + cimag(fineTotal) * cimag(fineTotal) );

   if (optVerbose)
   {
     printf("Integral of fine-grained frequency series = %g + %g i\n",
	    crealf(fineTotal), cimagf(fineTotal));
     printf("Integral of coarse-grained frequency series = %g + %g i\n",
	    crealf(coarseTotal), cimagf(coarseTotal));
     printf( "Fractional error is %e + % e i\n",
	     crealf(cError), cimagf(cError) );

     LALSPrintTimeSeries( &timeSeries, "ccg_TimeSeries.dat" );
     LALCPrintFrequencySeries( &goodInput, "ccg_fgFreqSeries.dat" );
     LALCPrintFrequencySeries( &goodOutput, "ccg_cgFreqSeries.dat" );
   }

   if ( creal(cError) * creal(cError) + cimag(cError) * cimag(cError)
	> ( CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL
	    * CCOARSEGRAINFREQUENCYSERIESTESTC_RANTOL )
      )
   {
     printf("  FAIL: Valid data test #4 \n");
     if (optVerbose)
     {
       printf("Exiting with error: %s\n",
	      CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS);
     }
     return CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS;
   }

   LALDestroyRandomParams( &status, &randomParams );
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALDestroyRealFFTPlan(&status, &fftPlan);
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALSDestroyVector(&status, &(timeSeries.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &(goodInput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCDestroyVector(&status, &(goodOutput.data));
   if ( ( code = CheckStatus(&status, 0 , "",
			     CCOARSEGRAINFREQUENCYSERIESTESTC_EFLS,
			     CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEFLS) ) )
   {
     return code;
   }

   LALCheckMemoryLeaks();

   printf("PASS: all tests\n");

   /**************** Process User-Entered Data, If Any **************/

   if (optInputFile[0] && optOutputFile[0]) {

     params.f0 = optF0;
     params.length = optOutLength;
     params.deltaF = optDeltaF;

     goodInput.data  = NULL;
     goodOutput.data = NULL;

     LALCCreateVector(&status, &goodInput.data, optInLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCCreateVector(&status, &goodOutput.data, optOutLength);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* Read input file */
     LALCReadFrequencySeries(&status, &goodInput, optInputFile);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     /* coarse grain */
     LALCCoarseGrainFrequencySeries(&status, &goodOutput, &goodInput, &params);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }

     LALCPrintFrequencySeries(&goodOutput, optOutputFile);

     printf("===== Coarse-Graining of User-Specified Series Written to File %s =====\n", optOutputFile);

     /* clean up valid data */
     LALCDestroyVector(&status, &goodInput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCDestroyVector(&status, &goodOutput.data);
     if ( ( code = CheckStatus( &status, 0 , "",
				CCOARSEGRAINFREQUENCYSERIESTESTC_EUSE,
				CCOARSEGRAINFREQUENCYSERIESTESTC_MSGEUSE) ) )
     {
       return code;
     }
     LALCheckMemoryLeaks();
   }
   return CCOARSEGRAINFREQUENCYSERIESTESTC_ENOM;
}

/*------------------------------------------------------------------------*/

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
  fprintf (stderr, "  -s ranseed        use random number seed ranseed\n");
  fprintf (stderr, "  -i filename    read fine grained series from file filename\n");
  fprintf (stderr, "  -o filename    print coarse grained series to file filename\n");
  fprintf (stderr, "  -n length      input series contains length points\n");
  fprintf (stderr, "  -m length      output series contains length points\n");
  fprintf (stderr, "  -e deltaF      set coarse grained frequency spacing to deltaF\n");
  fprintf (stderr, "  -f f0          set start frequency of output to f0\n");
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

    c = getopt (argc, argv, "hqvd:s:i:o:n:m:e:f:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'i': /* specify input file */
        strncpy (optInputFile, optarg, LALNameLength);
        break;

      case 'o': /* specify output file */
        strncpy (optOutputFile, optarg, LALNameLength);
        break;

      case 'n': /* specify number of points in input series */
        optInLength = atoi (optarg);
        break;

      case 'm': /* specify number of points in output series */
        optOutLength = atoi (optarg);
        break;

      case 'e': /* specify frequency resolution */
        optDeltaF = atof (optarg);
        break;

      case 'f': /* specify start frequency */
        optF0 = atof (optarg);
        break;

      case 's': /* set random number seed */
        optSeed = atoi (optarg);
        break;

      case 'd': /* set debug level */
        break;

      case 'v': /* optVerbose */
        optVerbose = CCOARSEGRAINFREQUENCYSERIESTESTC_TRUE;
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
